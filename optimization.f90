module optimization

   implicit none

   contains

      subroutine maxdLqdt
         USE global_variables
         use fftwfunction
         USE data_ops
         USE function_ops
         USE mpi
         IMPLICIT NONE

         REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: uvec0, gradJ0, gradJ1, gradJproj, gradJused0, gradJused1, diff_gradJ, unit_normal, d0, d1, vecTransported_GradJused0, vecTransported_d0, normalL2, normalHs
         REAL(pr), DIMENSION(:,:,:),   ALLOCATABLE :: f_scalar
         logical :: optimizationSuccessful

         REAL(pr) :: J0, J1, deltaJ, tau0, tau1, tau0Init, beta, inner, norm, test, alpha, HildertOrderS, gradJ0norm
         REAL(pr), DIMENSION(1:2) :: tau_brack
         real(pr), dimension(3) :: testVec, testVec2

         INTEGER :: iter, gradType, mnbrak_flag, FixConstr_flag, while_flag

         ALLOCATE( uvec0(1:n(1),1:n(2),1:local_N,1:3) )

         ALLOCATE( gradJ0(1:n(1),1:n(2),1:local_N,1:3) )
         ALLOCATE( gradJ1(1:n(1),1:n(2),1:local_N,1:3) )

         ALLOCATE( gradJproj(1:n(1),1:n(2),1:local_N,1:3) )


         ALLOCATE( normalL2(1:n(1),1:n(2),1:local_N,1:3) )
         ALLOCATE( normalHs(1:n(1),1:n(2),1:local_N,1:3) )

         ALLOCATE( gradJused0(1:n(1),1:n(2),1:local_N,1:3) )
         ALLOCATE( gradJused1(1:n(1),1:n(2),1:local_N,1:3) )

         ALLOCATE( d0(1:n(1),1:n(2),1:local_N,1:3) )
         ALLOCATE( d1(1:n(1),1:n(2),1:local_N,1:3) )


         ALLOCATE( vecTransported_GradJused0(1:n(1),1:n(2),1:local_N,1:3) )
         ALLOCATE( vecTransported_d0(1:n(1),1:n(2),1:local_N,1:3) )

         ALLOCATE( diff_gradJ(1:n(1),1:n(2),1:local_N,1:3) )
         ALLOCATE( f_scalar(1:n(1),1:n(2),1:local_N) )

         if(rank==0 .and. verboseOptimization) print*, "starting maxdLqdt"
         
         call initializeConstraintDirectory
         

         !====================================
         ! rescale to match the potentially changed constraint value
         !====================================
         if(rank==0 .and. verboseOptimization) print*, "rescaling"
         call rescaleLqNorm(Uvec,constraintB)
         call divAvg_free(uvec)


         
         !====================================
         ! potential kappa test
         !====================================
         if(kappaTest) then
            gradJ1 = GradL2ForLq(Uvec)
            call kappa_test(uvec, gradJ1, .true., "maxdLqdt", "initial grad L2", "L2")  
         end if

         !======================================================
         ! initialize loop variables
         !======================================================
         if(rank==0 .and. verboseOptimization) print*, "initializing loop variables"
         J0 = eval_J(Uvec, "maxdLqdt")
         J1 = 0.0_pr
         deltaJ = ABS( (J1-J0)/J0 )                            ! just to have something > OPTIM_TOL
         if(B_list_iterator==bIterOffset+1) then
            iter = optimizationIterOffset
         else
            iter = 0
         end if
         write(optimizationIterationTxt, '(i4.4)') iter
         d0 = 0.0_pr
         gradJused0 = 0.0_pr
         vecTransported_GradJused0 = 0.0_pr
         vecTransported_d0 = 0.0_pr
         tau0Init = 10.0_pr**(2.0_pr)
         tau0 = tau0Init
         tau1 = tau0Init
         !===============================
         ! PRINT INITIAL VALUES
         !===============================
         if(rank==0) print*, iter, "iteration"
         if(rank==0) print '(7A16, A12, A7, A12, 2A16)', achar(9), "tau_brack(a)", "tau used", "tau_brack(b)", "J0", "deltaJ", "J1", "lambda1", "bIter", "||u||_q", "Jv/(|Jv|+|Jn|)", "Jn/(|Jv|+|Jn|)"
         test = calc_global_Lq_norm(Uvec)
         testVec = calc_global_dLqdt_inclParts(uvec, lebesgueQ)
         testVec(1) = abs(testVec(2)) + abs(testVec(3))
         if(rank==0) print '(A16, 6ES16.7, ES12.3, I7, F12.7, 2F16.7)', achar(9), tau_brack(1), tau1, tau_brack(2), J0, deltaJ, J1, lambda1, B_list_iterator, test, testVec(2)/testVec(1), testVec(3)/testVec(1)


         IF (save_diag_Optim) THEN
            IF (rank==0) THEN
               CALL save_diagnostics_optim("maxdLqdt", iter, 0.0_pr, 0.0_pr, J0, (/0.0_pr, 0.0_pr, 0.0_pr/), (/0.0_pr, 0.0_pr, 0.0_pr/), 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr)
            END IF
            if(save_spectraEveryXiteration>0) then
               if(modulo(iter,save_spectraEveryXiteration) == 0) then
                  CALL calculateSaveSpectrum(uvec,"uvec")
               end if
            end if
         END IF

         while_flag = 0

         !lambda1 = 10.0_pr

         DO WHILE ( (ABS(deltaJ) > OPTIM_TOL) .AND. (iter<MAX_ITER) .AND. (while_flag<1) )
            iter = iter + 1
            write(optimizationIterationTxt, '(i4.4)') iter
            if(rank==0) print*, iter, "iteration"


            !======================================================
            ! CALCULATE BASIC GRADIENT IN THE H^s TOPOLOGY
            !======================================================
            if(rank==0 .and. verboseOptimization) print*, "calculating gradient"
            gradJ1 = GradL2ForLq(Uvec)
            if(useBanachGradient) then
               if(rank==0) print*, "todo change back to L2 gradient here"
               gradJ1 = BanachGradient(uvec)
            else
               hildertOrderS = (3.0_pr/2.0_pr)-(1.0_pr/lebesgueQ)
               call HilbertGradient(gradJ1, hildertOrderS)
            end if
            
            call divAvg_free(gradJ1)

            !======================================================
            ! ORTHOGONAL GRADIENT
            !======================================================
            !======================================================
            ! CALCULATE Normal Of Tangent Space
            !======================================================
            if(rank==0 .and. verboseOptimization) print*, "calculating normal of tangent"
            normalL2 = calcConstraintDerivativeL2(Uvec)     ! normal = q |u|^{q-2} u = nabla( ||u||_q^q )
            normalHs = normalL2
            if(useBanachGradient) then
               normalHs = BanachGradient(normalHs)
            else
               call HilbertGradient(normalHs, hildertOrderS)
            end if            
            !testVec = calc0thFourierModes(normalHs)
            !call divAvg_free(normalHs)
            !testVec = calc0thFourierModes(normalHs)
            !if(rank==0) print*, "0th order fourier mode of divfree normal HS", testVec(1), testVec(2), testVec(3)
            call divAvg_free(normalHs)
            alpha = global_summed_field_inner_product(normalHs,gradJ1,"H_l^(3/2-1/q)")
            alpha = alpha/global_summed_field_inner_product(normalHs,normalHs,"H_l^(3/2-1/q)")
            
            !======================================================
            ! PROJECT to tangent space, gradJ1 = Projection_Tang(u) (nabla^{H^\dots} J)
            !======================================================
            if(rank==0 .and. verboseOptimization) print*, "projecting"
            gradJproj(:,:,:,:) = gradJ1(:,:,:,:) - alpha*normalHs(:,:,:,:)
            !call divAvg_free(gradJproj)

            test = abs(global_summed_field_inner_product(gradJproj, normalHs, "H_l^(3/2-1/q)")/global_summed_field_inner_product(gradJproj+normalHs,gradJproj+normalHs, "H_l^(3/2-1/q)"))
            if(test>checkNormalTolerance) then
               if(rank==0) print*, "WARNING", achar(9), achar(9), "abs(inner(gradJproj,normalHs)/||gradJproj+normalHs||^2)", test, ">", checkNormalTolerance, "should be 0"
               call optim_msg_handle(42)
               !call kappa_test(uvec, normalHs, .false., "||u||_q", "projection normal", "H_l^(3/2-1/q)")
            !else
               !if(iter == 1 .and. kappaTest) call kappa_test(uvec, normalHs, .true., "||u||_q", "projection normal", "H_l^(3/2-1/q)")       !kappa test for projection onto tangent space normal = derivative
            end if
            call divergence(gradJproj, f_scalar)
            test = global_inner_product(f_scalar,f_scalar,"L2")/global_summed_field_inner_product(gradJproj,gradJproj,"H_l^(3/2-1/q)")
            if (test>checkDivergenceTolerance) then
               if(rank==0) print*, "WARNING", achar(9), achar(9), "||div(gradJproj)||_2^2/||gradJproj||_H^s", test, ">", checkDivergenceTolerance, "should be 0"
               call optim_msg_handle(43)
            end if

            
            
            
            !======================================================
            ! Base Descend Direction
            !======================================================
            if(useOrthogonalGradient) then
               gradJused1(:,:,:,:) = gradJproj(:,:,:,:)
            else
               gradJused1(:,:,:,:) = gradJ1(:,:,:,:)
            end if

            !======================================================
            ! CONJUGATE GRADIENT
            !======================================================
            if(useConjugateGradient) then
               if(rank==0 .and. verboseOptimization) print*, "conjugate gradient"
               !======================================================
               ! Calculate Momentum Term (Polak-Ribière)
               !======================================================
               beta = global_summed_field_inner_product(gradJused1,gradJused1-vecTransported_d0,"H_l^(3/2-1/q)")       !!! d0 or usedgraj0??? (8.29) in absil et al
               gradJ0norm = global_summed_field_inner_product(gradJused0,gradJused0,"H_l^(3/2-1/q)")
               !if(rank==0 .and. verboseOptimization) print*, "gradJ0norm", gradJ0norm
               if(gradJ0norm<MACH_EPSILON) then
                  beta = 0.0_pr
               else
                  beta = beta/gradJ0norm
               end if

               if(modulo(iter,resetMomentumTermEveryXiterations) == 0) then
                  if(rank==0) print*, achar(9), "resetting momentum term"
                  beta = 0.0_pr
               end if
               
               d1(:,:,:,:) = gradJused1(:,:,:,:) + beta*vecTransported_d0(:,:,:,:)
            else
               d1(:,:,:,:) = gradJused1(:,:,:,:)
            end if


            !======================================================
            ! normalize d = || u ||_H^s/|| d ||_H^s d
            ! so that tau measures ratio between old and new
            !======================================================
            
            if(normalizeDirection) then
               if(rank==0 .and. verboseOptimization) print*, "normalize direction"
               d1(:,:,:,:) =  sqrt(global_summed_field_inner_product(uvec,uvec,"H_l^(3/2-1/q)")/global_summed_field_inner_product(d1,d1,"H_l^(3/2-1/q)")) * d1(:,:,:,:)
            end if

            if(rank==0 .and. verboseOptimization) print*, "check average zero d1"
            testVec = calc0thFourierModes(d1)
            if(sum(abs(testVec(:)))>checkAverageTolerance*constraintB) then
               if(rank==0) print*, "WARNING", achar(9), achar(9), checkAverageTolerance*constraintB,"< abs(average(gradJproj)) =", testVec(1), testVec(2), testVec(3)
               call optim_msg_handle(44)
            end if
            if(rank==0 .and. verboseOptimization) print*, "check average zero uvec"
            testVec = calc0thFourierModes(uvec)
            if(sum(abs(testVec(:)))>checkAverageTolerance*constraintB) then
               if(rank==0) print*, "WARNING", achar(9), achar(9), checkAverageTolerance*constraintB,"< abs(average(gradJproj)) =", testVec(1), testVec(2), testVec(3)
               call optim_msg_handle(44)
            end if
            
            !======================================
            ! FIND OPTIMAL tau BY ARC OPTIMIZATION
            !====================================== 
            if(rank==0 .and. verboseOptimization) print*, "find tau"

            CALL optim_msg_handle(20) 
            tau_brack = mnbrak(iter, "maxdLqdt", Uvec, d1, 0.0_pr, tau0, mnbrak_flag)

            IF (mnbrak_flag /= 0) THEN
               IF (rank==0) print*, "mnbrak error"
               select case (mnbrak_flag)
                  CASE (1)
                     IF (rank==0) print*, achar(9), "Going uphill... Verify gradient!"   !achar(9) = ascii code 9 = tab
                     call optim_msg_handle(22)
                  CASE (2)
                     IF (rank==0) print*, achar(9), "Could not bracket minimum."
                  CASE (3)
                     IF (rank==0) print*, achar(9), "Decreasing tau..."
               end select
               CALL optim_error_handle(mnbrak_flag)            
               IF (save_diag_Optim) THEN
                  IF (rank==0) THEN
                     CALL save_diagnostics_optim("maxdLqdt", iter, tau1, beta, J0, (/0.0_pr, 0.0_pr, 0.0_pr/), (/0.0_pr, 0.0_pr, 0.0_pr/), 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr)
                  END IF
               END IF
               while_flag = 1
               exit  !exits loop 
            ELSE
               CALL optim_msg_handle(21)
            END If
            CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

            CALL optim_msg_handle(30)

            tau1 = brent(iter, "maxdLqdt", Uvec, d1, tau_brack)    ! I add the new variable iter
            if(tau1>tau_brack(2)) then
               if(rank==0 .and. tauDebugToConsole) print*, "warning! brent tau", tau1, "> tau brak(2)", tau_brack(2)
               if(tau1>2.0*tau_brack(2)) then
                  CALL optim_msg_handle(32)     ! Optimal tau is too large
                  if(rank==0 .and. tauDebugToConsole) print*, "WARNING! brent tau", tau1, "is too large!"
               end if
            else
               CALL optim_msg_handle(31)
            end if


            !IF (rank==0) print*, "tau1", tau1
            if(tau1>tau_max) then
               if (rank==0) print*, "WARNING: tau1", tau1, "> TAU_MAX", tau_max, " terminating optimzation, not smooth enough?"
               call optim_msg_handle(32)
               while_flag = 1
               tau1 = 0.0_pr
            elseif(tau1>sqrt(tau_max)) then
               if (rank==0) print*, achar(9), achar(9), "warning: tau1", tau1, ">", sqrt(tau_max)
               call optim_msg_handle(41)
            end if


            !======================================
            ! update vector transport FOR NEXT STEP with "old" tangentspace at "old" u
            !======================================
            if(useRiemannianGeometry) then
               if(rank==0 .and. verboseOptimization) print*, "update vector transport"
               vecTransported_GradJused0 = vectorTransport(Uvec, lebesgueQ, constraintB, tau1*d1, gradJused1)
               vecTransported_d0 = vectorTransport(Uvec, lebesgueQ, constraintB, tau1*d1, d1)
            else
               vecTransported_GradJused0 = gradJused1
               vecTransported_d0 = d1
            end if
            if(iter == 1 .and. kappaTest) call kappa_test(uvec, d1, .true., "maxdLqdt", "start_d1", "H_l^(3/2-1/q)")  

            !======================================
            ! update u
            !======================================
            if(rank==0 .and. verboseOptimization) print*, "update u"
            uvec0 = uvec
            Uvec = Uvec + tau1*d1
            call rescaleLqNorm(uvec, constraintB)
                     

            J1 = eval_J(Uvec, "maxdLqdt")
            deltaJ = (J1-J0)/ABS(J0)
            IF (deltaJ < -MACH_EPSILON) THEN      ! Change om March 30, 2017
                IF (save_diag_Optim) THEN
                    IF (rank==0) THEN
                        CALL save_diagnostics_optim("maxdLqdt", iter+1, tau1, beta, J0, (/0.0_pr, 0.0_pr, 0.0_pr/), (/0.0_pr, 0.0_pr, 0.0_pr/), 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr)   
                    END IF
                END IF
                CALL optim_msg_handle(0)
                if(rank==0) print*, "WARNING! Cost functional not increasing, exiting"
                if(rank==0 .and. tauDebugToConsole) then
                  print '(4A20)', achar(9), "tau_brack(a)", "tau used", "tau_brack(b)"
                  print '(A20, 3ES20.12)', achar(9), tau_brack(1), tau1, tau_brack(2)
                end if
                !save_data_optim = .FALSE.
                while_flag = 1
                EXIT
            ELSE
               !save_data_optim = .TRUE.
            END IF
            CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)     ! for everyone to know save_data_optim = .FALSE., otherwise next if barrier might wait indefinitely
                

            IF (save_diag_Optim) THEN
               if(rank==0 .and. verboseOptimization) print*, "saving diag optim"
               IF (rank==0) THEN
                  CALL save_diagnostics_optim("maxdLqdt", iter, tau1, beta, J1, (/0.0_pr, 0.0_pr, 0.0_pr/), (/0.0_pr, 0.0_pr, 0.0_pr/), 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr)
               END IF
               if(save_spectraEveryXiteration>0) then
                  if(modulo(iter,save_spectraEveryXiteration) == 0) then
                     CALL calculateSaveSpectrum(uvec,"uvec")
                     !CALL calculateSaveSpectrum(gradJ1,"gradJ")
                     !CALL calculateSaveSpectrum(gradJproj,"gradJproj")
                     CALL calculateSaveSpectrum(d1,"d")
                  end if
               end if
               if(save_scalarFieldsEveryXiteration>0) then
                  if(modulo(iter, save_scalarFieldsEveryXiteration) == 0) call diagnosticScalars(Uvec)
               end if
               if(save_uvecEveryXiteration>0) then
                  if(modulo(iter, save_uvecEveryXiteration) == 0) call save_Ctrl(Uvec, iter, "maxdLqdt_result")
               end if
            END IF
            CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)


            !===============================
            ! PRINT RESULTS
            !===============================
            if(rank==0) print '(7A16, A12, A7, A12, 2A16)', achar(9), "tau_brack(a)", "tau used", "tau_brack(b)", "J0", "deltaJ", "J1", "lambda1", "bIter", "||u||_q", "Jv/(|Jv|+|Jn|)", "Jn/(|Jv|+|Jn|)"
            test = calc_global_Lq_norm(Uvec)
            testVec = calc_global_dLqdt_inclParts(uvec, lebesgueQ)
            testVec(1) = abs(testVec(2)) + abs(testVec(3))
            if(rank==0) print '(A16, 6ES16.7, ES12.3, I7, F12.7, 2F16.7)', achar(9), tau_brack(1), tau1, tau_brack(2), J0, deltaJ, J1, lambda1, B_list_iterator, test, testVec(2)/testVec(1), testVec(3)/testVec(1)


            !===============================
            ! UPDATE OLD VARIABLES
            !===============================   
            if(rank==0 .and. verboseOptimization) print*, "update old variables"
            J0 = J1
            gradJ0 = gradJ1
            gradJused0 = gradJused1
            tau0 = tau1
            d0 = d1



            !if(lambda1>0.1) lambda1=lambda1/1.1

         END DO


         if(iter<MAX_ITER .and. while_flag<1) then
            optimizationSuccessful = .true.
            if(rank==0) then
               print*, "optimization terminated successful after", iter, "iterations"
               print*, " "
            end if
            if (kappaTest) then
               !CALL kappa_test(uvec0, gradJ1, "end_gradJ1", "H_l^(3/2-1/q)")
               CALL kappa_test(uvec0, d1, .true.,"maxdLqdt", "end_d1", "H_l^(3/2-1/q)")
            end if
         else
            optimizationSuccessful = .false.
            if(rank==0) then
               print*, "optimization terminated unsuccessful", iter, MAX_ITER
            end if
            CALL kappa_test(uvec0, d1, .true., "maxdLqdt", "end_d1", "H_l^(3/2-1/q)")
            CALL kappa_test(uvec0, gradJ1, .true., "maxdLqdt", "end_gradJ", "H_l^(3/2-1/q)")
         end if


         IF (save_data_Optim) THEN
            optimizationIterationTxt = "end"
            if(rank==0 .and. verboseOptimization) print*, "saving data optim"
            CALL diagnosticScalars(Uvec)
            CALL save_Ctrl(Uvec, iter, "maxdLqdt_result")
            CALL calculateSaveSpectrum(uvec,"uvec")
            CALL calculateSaveSpectrum(gradJproj,"gradJproj")
            CALL calculateSaveSpectrum(gradJ0,"gradJ")
            CALL calculateSaveSpectrum(d0,"d")
         END IF

         CALL optim_msg_handle(1)

         !===============================
         ! save results
         !===============================
         call save_to_optimizationResultList(constraintB, J1, iter, optimizationSuccessful)
         
         deallocate(uvec0)
         deallocate(gradJ0)
         deallocate(gradJ1)
         deallocate(normalL2)
         deallocate(normalHs)
         deallocate(gradJproj)
         deallocate(gradJused0)
         deallocate(gradJused1)
         deallocate(d0)
         deallocate(d1)
         deallocate(vecTransported_GradJused0)
         deallocate(vecTransported_d0)
         DEALLOCATE(diff_gradJ)
         DEALLOCATE(f_scalar)


         CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
         if(rank==0 .and. verboseOptimization) print*, "ending maxdLqdt"

      end subroutine


   !=================================================
   ! FUNCTION THAT CALCULATES COST FUNCTIONAL
   !=================================================
   RECURSIVE FUNCTION eval_J(myfield, mysystem) RESULT (J)
      USE global_variables
      USE function_ops
      USE mpi
      IMPLICIT NONE

      REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(IN) :: myfield
      REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: aux1
      CHARACTER(len=*), INTENT(IN) :: mysystem
      REAL(pr) :: local_J, J

      allocate( aux1(1:n(1),1:n(2),1:local_N,1:3) )
      aux1 = myfield
      SELECT CASE (mysystem)
         case ("maxdLqdt")
            local_J = calc_local_dLqdt(aux1, lebesgueQ)
            call mpi_allreduce(local_J, J, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
         CASE ("LineMin")
            CALL rescaleLqNorm(aux1, constraintB)
            J = -1.0_pr*eval_J(aux1, "maxdLqdt")
         case ("||u||_q")
            J = calc_global_Lq_norm(aux1)
         case DEFAULT
            IF (rank==0) print*, "WARNING, case ", mysystem, " for mysystem not found in eval_grad_J"
            J = 0.0_pr
      END SELECT
      deallocate(aux1)
   END FUNCTION eval_J




   !=================================================
   ! iteratively calculate the banach gradient by solving
   ! |v|^{s-2} ((s-2)e_{v_l} e_{v_m} + delta_lm) v_l^{k+1}
   !        - partial_i(|nabla v|^{s-2} ((s-2) e_{partial_i v_m} e_{partial_j v_l} + delta_ij delta_lm) partial_j v_l^{k+1})
   !        = s |v|^{s-2} v_m - s partial_i(|nabla v|^{s-2} partial_i v_m) - 1/lambda nabla rho + 1/lambda L2Grad
   ! and
   ! Delta rho = lambda nabla cdot (|v^{k+1}|^{s-2} v^{k+1} - partial_i(|nabla v^{k+1}|^{s-2} partial_i v^{k+1})
   ! in an alternating fashion
   !=================================================
   function banachGradient(l2Grad) RESULT (v)
      USE global_variables
      USE function_ops
      USE mpi
      IMPLICIT NONE

      REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(IN) :: l2Grad
      real(pr) :: lambda
      real(pr) :: s                    ! s = 3q/(q+1)
      real(pr), dimension(1:n(1),1:n(2),1:local_N,1:3) :: v, v_old
      real(pr), dimension(1:n(1),1:n(2),1:local_N) :: rho
      real(pr) :: tau         ! step size
      integer :: bgIter
      real(pr) :: residual = 99999.9_pr
      real(pr) :: glob_norm_v_sq


      s = 3.0_pr*lebesgueQ/(lebesgueQ+1.0_pr)

      rho = 0.0_pr
      lambda = 1.0_pr
      v_old = l2Grad
      tau = 1.0_pr

      do while (bgIter <= banachGradIterMax .and. residual>banachIterTol)
         bgIter = bgIter + 1

         !if(tau>0.01_pr) then
         !   tau = 0.9_pr*tau
         !else
         !   tau = 0.01_pr
         !end if
         v = BanachGradientIterationOld(l2Grad, v_old, lambda, rho, tau)
         rho = BanachGradientCalcRho(v, lambda)


         residual = global_summed_field_inner_product(v-v_old,v-v_old,"L2")/global_summed_field_inner_product(v,v,"L2")
         glob_norm_v_sq = global_summed_field_inner_product(v,v,"L2")
         if(rank==0) print*, "banach gradient", " iter", bgIter, "||v-v_old||_2^2/||v||_2^2", residual, "||v||_2^2", glob_norm_v_sq

         v_old = v

      end do
      
   END FUNCTION banachGradient

   !==================================================
   ! BRACKET THE LOCATION OF OPTIMAL TAU
   ! Press Teukolsky Vetterling Flannery - Numerical Recipes - 10.1 Initially Bracketing a Minimum
   !==================================================

   FUNCTION mnbrak(optimizationIter, mysystem, phi, grad, tA0, tB0, myflag) RESULT (tau_brack)
      USE global_variables
      USE data_ops
      USE function_ops
      IMPLICIT NONE
  
      CHARACTER(len=*), INTENT(IN) :: mysystem
      REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(IN) :: phi
      REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(IN) :: grad
      REAL(pr), INTENT(IN) :: tA0, tB0
      INTEGER, INTENT(INOUT) :: myflag
      REAL(pr), DIMENSION(1:2) :: tau_brack

      REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: phi_bar
      REAL(pr) :: aux, tP, FP, Pmax, R, Q
      REAL(pr) :: FA, FB, FC, tA, tB, tC
 
      REAL(pr), PARAMETER :: GOLD = (1.0_pr + SQRT(5.0_pr))/2.0_pr
      REAL(pr), PARAMETER :: CGOLD = 1.0_pr/GOLD  
      REAL(pr), PARAMETER :: GLIMIT = 10.0_pr
      REAL(pr), PARAMETER :: tMAX = 10.0_pr
      INTEGER, PARAMETER :: ITMAX = 100       ! maximal iterations in mnbrak method
      INTEGER :: FuncEval, iter, optimizationIter
      LOGICAL :: saveLineMin

      saveLineMin = .TRUE.

      ALLOCATE( phi_bar(1:n(1),1:n(2),1:local_N,1:3) )



      if(mnbra_calcSaveAllJvalues) then
         tB = MAX(100.0*tB0, MACH_EPSILON)
         phi_bar = phi + tB*grad
         FB = eval_J(phi_bar, "LineMin")
         IF (saveLineMin) CALL save_linemin_data(0.0_pr, tB, 0.0_pr, 0.0_pr, FB, 0.0_pr, -1, optimizationIter, mysystem, "replace")
         DO WHILE (tB > 1.0e-5_pr) 
            tB = CGOLD*tB
            phi_bar = phi + tB*grad
            FB = eval_J(phi_bar, "LineMin")
            FuncEval = FuncEval+1
            IF (saveLineMin) CALL save_linemin_data(0.0_pr, tB, 0.0_pr, 0.0_pr, FB, 0.0_pr, -1, optimizationIter, mysystem, "append")
         END DO
         tB = MAX(tB0, MACH_EPSILON)
      end if


      FuncEval = 0
      iter = 0      
 
      tA = tA0
      tB = MAX(tB0, MACH_EPSILON)

      phi_bar = phi + tA*grad
      FA = eval_J(phi_bar, "LineMin")
      FuncEval = FuncEval+1

      phi_bar = phi + tB*grad
      FB = eval_J(phi_bar, "LineMin")
      FuncEval = FuncEval+1


      IF (saveLineMin) CALL save_linemin_data(tA, tB, tC, FA, FB, FC, iter, optimizationIter, mysystem, "replace")

      DO WHILE (FB > FA .AND. tB > MACH_EPSILON) 
         tB = CGOLD*tB
         phi_bar = phi + tB*grad
         FB = eval_J(phi_bar, "LineMin")
         FuncEval = FuncEval+1
         IF (saveLineMin) CALL save_linemin_data(tA, tB, tC, FA, FB, FC, iter, optimizationIter, mysystem, "append")
      END DO

      IF (tB .LE. MACH_EPSILON) THEN
         myflag = 1
         RETURN
      END IF

      tC = GOLD*tB
      phi_bar = phi + tC*grad
      FC = eval_J(phi_bar, "LineMin")
      FuncEval = FuncEval+1

      IF (saveLineMin) CALL save_linemin_data(tA, tB, tC, FA, FB, FC, iter, optimizationIter, mysystem, "append")

      DO WHILE (FB>=FC .AND. iter<ITMAX)
         iter = iter+1
         tC = GOLD*tC
         phi_bar = phi + tC*grad
         FC = eval_J(phi_bar, "LineMin")
         FuncEval = FuncEval+1
                  
         R = (tB-tA)*(FB-FC)
         Q = (tB-tC)*(FB-FA)
         tP = tB - 0.5_pr*((tB-tC)*Q - (tB-tA)*R)/( SIGN( MAX(ABS(Q-R),MACH_EPSILON), Q-R) )
    
         Pmax = tB + GLIMIT*(tC-tB)
    
         IF ( (tB-tP)*(tP-tC)>0 ) THEN
            phi_bar = phi + tP*grad
            FP = eval_J(phi_bar, "LineMin")
        
            IF (FP<FC) THEN
               tA = tB
               FA = FB
               tB = tP
               FB = FP
               EXIT 
            ELSEIF (FP>FB) THEN
               tC = tP
               FC = FP
               EXIT
            END IF
        
            tP = tC + GOLD*(tC-tB)
            phi_bar = phi + tP*grad
            FP = eval_J(phi_bar, "LineMin")
        
         ELSEIF ( (tC-tP)*(tP-Pmax)>0 ) THEN
            phi_bar = phi + tP*grad
            FP = eval_J(phi_bar, "LineMin")
       
            IF (FP<FC) THEN
               tB = tC
               tC = tP
               FB = FC
               FC = FP
               tP = tC+GOLD*(tC-tB)
               phi_bar = phi + tP*grad
               FP = eval_J(phi_bar, "LineMin")
            END IF
        
        ELSEIF ( (tP-Pmax)*(Pmax-tC)>=0 ) THEN
            tP = Pmax
            phi_bar = phi + tP*grad
            FP = eval_J(phi_bar, "LineMin")
        
         ELSE
            tP = tC + GOLD*(tC-tB)
            phi_bar = phi + tP*grad
            FP = eval_J(phi_bar, "LineMin")
         END IF
    
         tA = tB
         tB = tC
         tC = tP
         FA = FB
         FB = FC
         FC = FP
        
         IF (saveLineMin) CALL save_linemin_data(tA, tB, tC, FA, FB, FC, iter, optimizationIter, mysystem, "append")
 
      END DO

      tau_brack(1) = tA
      tau_brack(2) = tC

      IF (iter .GE. ITMAX) THEN
         myflag = 2
      ELSE
         myflag = 0
      END IF

      DEALLOCATE( phi_bar )
    

   END FUNCTION mnbrak


   !============================================
   ! BRENT ALGORITHM FOR LINE OPTIMIZATION
   !============================================
   FUNCTION brent(iteration, mysystem, phi, grad, tau_brack) RESULT (X)   ! I add new variable iteration on Feb. 2017
      USE global_variables
      USE data_ops
      USE function_ops
      use mpi
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: iteration
      CHARACTER(len=*), INTENT(IN) :: mysystem
      REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(IN) :: phi, grad
      REAL(pr), DIMENSION(1:2), INTENT(IN) :: tau_brack
      REAL(pr) :: X
   
      REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: phi_bar
      REAL(pr) :: D, A, B, V, W, E, ETEMP, P, Q, R, U, XM
      REAL(pr) :: FV, FW, FU, FX
      REAL(pr) :: TOL1, TOL2
      INTEGER :: FLAG, j


      INTEGER, PARAMETER :: ITMAX = 300       ! Maximal iterations in brent method
      REAL(pr), PARAMETER :: TOL = 1E-6
      REAL(pr), PARAMETER :: ZEPS = 1E-8
      REAL(pr), PARAMETER :: CGOLD = .381966

      CHARACTER(2) :: K0txt, E0txt, IGtxt
      CHARACTER(4) :: itertxt
      CHARACTER(100) :: filename
      WRITE(K0txt, '(i2.2)') K0_index
      WRITE(E0txt, '(i2.2)') E0_index
      WRITE(IGtxt, '(i2.2)') iguess
      WRITE(itertxt, '(i4.4)') iteration


      ALLOCATE( phi_bar(1:n(1),1:n(2),1:local_N,1:3) )
 
      D = 0.0_pr
      A = MIN(tau_brack(1),tau_brack(2))
      B = MAX(tau_brack(1),tau_brack(2))
      V = tau_brack(2)*CGOLD 
      W = V
      X = V 
      E = 0.0_pr


      call createDirectoryIfNonExistent(ConstraintDir//"tau-data")

      !filename = HomeDir//"/brent_info"//".dat"
      filename = ConstraintDir//"tau-data/"//"brent-info-"//itertxt//".dat"
      !filename = HomeDir//"/brent_info_maxdEdtHeli_Nu_E"//E0txt//"_IG"//IGtxt//"_brent_info.dat"
 
      phi_bar = phi + D*grad
      FX = eval_J(phi_bar, "LineMin")

      if(rank==0) then
         OPEN(10, FILE = filename, FORM = 'FORMATTED', STATUS = 'REPLACE')
         WRITE(10, "(3 G20.12)") "#", "Tau", "J"
         WRITE(10, "(I20, 2 ES20.12)") 0, D, FX
         CLOSE(10)
      end if
      call mpi_barrier(mpi_comm_world, statinfo)
      
      phi_bar = phi + X*grad
      FX = eval_J(phi_bar, "LineMin")
 
      FV = FX 
      FW = FX


      DO j=1,ITMAX


         if(rank==0) then
            OPEN(10, FILE = filename, FORM = 'FORMATTED', STATUS = 'OLD', POSITION = 'APPEND')
            WRITE(10, "(I20, 2 ES20.12)") j, X, FX
            CLOSE(10)
         end if
         call mpi_barrier(mpi_comm_world, statinfo)

         XM = 0.5_pr*(A+B)
         TOL1 = TOL*ABS(X)+ZEPS
         TOL2 = 2.0_pr*TOL1
    
         IF ( ABS(X-XM) <= (TOL2-0.5*(B-A)) ) EXIT
    
         FLAG = 1
         IF ( ABS(E) > TOL1 ) THEN
            R = (X-W)*(FX-FV)
            Q = (X-V)*(FX-FW)
            P = (X-V)*Q - (X-W)*R
            Q = 2.0_pr*(Q-R)
            IF ( Q > 0.0_pr ) P = -P 
      
            Q = ABS(Q)
            ETEMP = E
            E = D
        
            IF ( (ABS(P) >= ABS(0.5_pr*Q*ETEMP)) .OR. (P <= Q*(A-X)) .OR. (P >= Q*(B-X)) ) THEN
               FLAG = 1
            ELSE
               FLAG = 2
            END IF

         END IF
    
         SELECT CASE (FLAG)
            CASE (1)
               IF (X >= XM) THEN
                  E = A-X
               ELSE
                  E=B-X
               END IF
               D = CGOLD*E
            CASE (2)
               D = P/Q
               U = X+D
               IF ( (U-A < TOL2) .OR. (B-U < TOL2) ) D = SIGN(TOL1, XM-X)
         END SELECT
    
         IF ( ABS(D) >= TOL1 ) THEN
            U = X+D
         ELSE
            U = X + SIGN(TOL1,D)
         END IF
    
         phi_bar = phi + U*grad
         FU = eval_J(phi_bar, "LineMin")
 
         IF ( FU <= FX ) THEN
            IF ( U >= X ) THEN
               A = X
            ELSE
               B = X 
            END IF
            V = W
            FV = FW
            W = X
            FW = FX
            X = U
            FX = FU
         ELSE
            IF ( U < X ) THEN
               A = U
            ELSE
               B = U
            END IF
        
            IF ( (FU <= FW) .OR. (W == X) ) THEN
               V = W
               FV = FW
               W = U
               FW = FU
            ELSEIF ( (FU <= FV) .OR. (V==X) .OR. (V==W)) THEN 
               V = U
               FV = FU
            END IF
         END IF
      END DO


      if(rank==0) then
         OPEN(10, FILE = filename, FORM = 'FORMATTED', STATUS = 'OLD', POSITION = 'APPEND')
         WRITE(10, "(I20, 2 ES20.12)") 999, X, FX
         CLOSE(10)
      end if
      call mpi_barrier(mpi_comm_world, statinfo)

      DEALLOCATE( phi_bar )

   END FUNCTION brent


   !============================================
   ! Perform kappa_test
   !============================================
   SUBROUTINE kappa_test(phi, gradJ, useAdjustedKappaTest, mysystem, nameOfKappaTest, inner_product_space)
      USE global_variables
      USE function_ops
      USE data_ops
      USE mpi
      use fftwfunction
      IMPLICIT NONE

      !REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(IN) :: phi
      CHARACTER(len=*), INTENT(IN) :: mysystem
      logical, intent(in) :: useAdjustedKappaTest
      REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(IN) :: phi
      REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(IN) :: gradJ
      !REAL(pr), INTENT(IN) :: J0
      REAL(pr) :: J0
      CHARACTER(len=*), INTENT(IN) :: nameOfKappaTest, inner_product_space
      REAL(pr) :: myepsilon, kappa, kappa_org, deltaJ
      real(pr), dimension(:,:), allocatable :: kappaArray
      REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: phi_pert, phi_bar
      REAL(pr) :: myexp, J1, innerProd
      REAL(pr), DIMENSION(1:3) :: dx
      integer :: kappaTestSize = 16
      INTEGER :: ii, jj
      
      ! factor adjusted kappa
      integer :: iAdjMin, adjAvgSize
      real(pr) :: fracAdjMin, frac

      CHARACTER(2) :: K0txt, E0txt, IGtxt
      CHARACTER(4) :: strLebesgueQ
      CHARACTER(3) :: tempStr
      character(50) :: phiPertText		! load random b/divfree sine/load te0080/...
      


      complex(pr), dimension(1:n(1),1:n(2),1:local_n,1:3) :: testRemoveAfterwards,testRemoveAfterwards2
      real(pr), dimension(1:n(1),1:n(2),1:local_n) :: testRemoveAfterwardsCompX, testRemoveAfterwardsCompY, testRemoveAfterwardsCompZ
         


      IF (rank==0) THEN
         print*, "kappaTest ", nameOfKappaTest
         !print*, "q", lebesgueQ
      END If
      WRITE(tempStr,'(F3.1)') lebesgueQ-int(lebesgueQ)              ! might result in rounding errors for 1.999999999999
      WRITE(strLebesgueQ,'(I2.2,a2)') int(lebesgueQ), tempStr(2:)
      
   
      IF (rank==0) THEN
      !   print*, "warning, fixing lq norm of uvec"
      END If
      !call Fix_Lq(uvec, 1.0_pr)


      WRITE(K0txt,'(i2.2)') K0_index
      WRITE(E0txt,'(i2.2)') E0_index
      WRITE(IGtxt,'(i2.2)') iguess
      ALLOCATE( phi_pert(1:n(1),1:n(2),1:local_N,1:3) )
      ALLOCATE( phi_bar(1:n(1),1:n(2),1:local_N,1:3) )




      !phiPertText = "divfree-sine"
      !phiPertText = "load-random-b"
      !phiPertText = "load-random-exp-b"
      !phiPertText = "load-k-random-a"
      phiPertText = "load-random-poly-b"
      !phiPertText = "load-random-smooth-b"
      !phiPertText = "save-random"	! generate new random field, stop afterwards and copy to it input folder
      !phiPertText = "load-te0080"
      phiPertText = trim(phiPertText)
      call kappa_test_pert(phi_pert, phiPertText, -5.0_pr, 0.0_pr, 0.0_pr)
      

      

      !call calculateSaveSpectrum(phi, "phi")
      !call calculateSaveSpectrum(phi_pert, trim("phiPert-"//phiPertText))
      !call calculateSaveSpectrum(gradJ, "gradJ"//"-q-"//strLebesgueQ)


      innerProd = global_summed_field_inner_product(gradJ, phi_pert, inner_product_space)


      J0 = eval_J(phi, mysystem)
      allocate( kappaArray(1:kappaTestSize,1:5) )

      DO ii=1,kappaTestSize
         
         myexp = -7.0_pr - real(kappaTestSize, pr)/2.0 + real(ii,pr)
         myepsilon = 10.0_pr**myexp         
         
         phi_bar = phi + myepsilon*phi_pert
         CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

         J1 = eval_J(phi_bar, mysystem)

         deltaJ = J1-J0

         kappa = (J1-J0)/(myepsilon*innerProd)

         CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

         kappaArray(ii,:) = (/myepsilon, innerProd, deltaJ, kappa, kappa/)

         if (rank==0) print*, achar(9), achar(9), "kappa test", ii, "/", kappaTestSize, kappa, LOG10(ABS(kappa - 1.0_pr))

      END DO


      !!!!! CALCULATE ADJUSTED KAPPA TEST !!!!!
      !!! tries to find a factor c such that c*gradJ is the acutal gradient
      !!! can be used for d = projected( grad J ) + beta vector(transport)
      !!! since this is probably not normalized in the correct way for the kappa test
      iAdjMin = 1
      adjAvgSize = 4
      fracAdjMin = 9999
      if(rank==0) print*, ""
      do ii=1,kappaTestSize-(adjAvgSize-1)-1
         frac = 0.0_pr
         do jj=0,adjAvgSize-1
            frac = frac + abs(1-kappaArray(ii+jj,4)/kappaArray(ii+jj+1,4))
         end do
         if(frac<fracAdjMin) then
            fracAdjMin = frac
            iAdjMin = ii
         end if
      end do

      ! best guess for adjusting
      frac = 1.0_pr
      do jj=0,adjAvgSize-1
         frac = frac*abs(kappaArray(iAdjMin+jj,4))
      end do
      frac = sign( (frac)**(1.0_pr/(adjAvgSize)) , kappaArray(iAdjMin,4) )
      kappaArray(:,4) = kappaArray(:,4)/frac


      CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
      
      

      if(rank==0) then
         DO ii=1,kappaTestSize
            !(/myepsilon, kappa/) = kappaArray(ii,:) !WHY FORTRAN CAN'T I DO SOMETHING LIKE THIS!!!
            myepsilon = kappaArray(ii,1)
            innerProd = kappaArray(ii,2)
            deltaJ = kappaArray(ii,3)
            kappa = kappaArray(ii,4)
            kappa_org = kappaArray(ii,5)
            if(useAdjustedKappaTest) then
               print*, achar(9), achar(9), "kappa test adjusted", ii, "/", kappaTestSize, kappa, LOG10(ABS(kappa - 1.0_pr))
               CALL save_kappa_test(myepsilon, innerProd, deltaJ, kappa, kappa_org, frac, ii, nameOfKappaTest)
            end if
         END DO
      end if

      deallocate( kappaArray )

      DEALLOCATE(phi_pert)
      DEALLOCATE(phi_bar)

   END SUBROUTINE kappa_test

END MODULE
