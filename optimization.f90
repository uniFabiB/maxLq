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

         REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: gradJ0, gradJ1, diff_gradJ, unit_normal, d0, d1, vecTransported_GradJ0, vecTransported_d0
         REAL(pr), DIMENSION(:,:,:),   ALLOCATABLE :: f_scalar

         REAL(pr) :: J0, J1, deltaJ, ell1, ell2, tau0, tau1, beta, local_scalar_L2norm, divU_L2, divGradJ_L2, inner, norm, test
         REAL(pr), DIMENSION(1:3) :: local_field_L2norm, local_gradJK0, local_gradJK1, K, E, gradJ_K0, gradJ_K1
         REAL(pr), DIMENSION(1:2) :: tau_brack
         REAL(pr), DIMENSION(1:3) :: dLqdt

         INTEGER :: iter, gradType, mnbrak_flag, FixConstr_flag

         ALLOCATE( gradJ0(1:n(1),1:n(2),1:local_N,1:3) )
         ALLOCATE( gradJ1(1:n(1),1:n(2),1:local_N,1:3) )

         ALLOCATE( d0(1:n(1),1:n(2),1:local_N,1:3) )
         ALLOCATE( d1(1:n(1),1:n(2),1:local_N,1:3) )


         ALLOCATE( vecTransported_GradJ0(1:n(1),1:n(2),1:local_N,1:3) )
         ALLOCATE( vecTransported_d0(1:n(1),1:n(2),1:local_N,1:3) )

         ALLOCATE( diff_gradJ(1:n(1),1:n(2),1:local_N,1:3) )
         ALLOCATE( f_scalar(1:n(1),1:n(2),1:local_N) )

         if(rank==0 .and. verboseOptimization) print*, "starting maxdLqdt"
         
         call initializeConstraintDirectory
         

         !====================================
         ! rescale to match the potentially changed constraint value
         !====================================
         if(rank==0 .and. verboseOptimization) print*, "rescaling"
         call rescaleLqNorm(Uvec,lebesgueQ,constraintB)




         !====================================
         ! CALCULATE DIAGNOSTICS OF CONTROL
         !====================================
         if(rank==0 .and. verboseOptimization) print*, "calculating diagnostics"
         local_field_L2norm = Energy(Uvec)
         CALL MPI_ALLREDUCE(local_field_L2norm, K, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)



         CALL divergence(Uvec, f_scalar)
         local_scalar_L2norm = inner_product(f_scalar, f_scalar, "L2")
         CALL MPI_ALLREDUCE(local_scalar_L2norm, divU_L2, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo) 
   
         local_field_L2norm = Enstrophy(Uvec)
         CALL MPI_ALLREDUCE(local_field_L2norm, E, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)


         dLqdt = calc_dLqdt(Uvec, lebesgueQ)
         

         !======================================================
         ! initialize loop variables
         !======================================================
         if(rank==0 .and. verboseOptimization) print*, "initializing loop variables"
         J0 = eval_J(Uvec, "maxdLqdt")
         J1 = 1.5_pr*J0
         deltaJ = ABS( (J1-J0)/J0 )                            ! just to have something > OPTIM_TOL
         iter = 0
         d0 = 0.0_pr
         vecTransported_GradJ0 = 0.0_pr
         vecTransported_d0 = 0.0_pr
         tau0 = 10.0_pr**(-4.0_pr)
         tau1 = tau0
         
         if (kappaTest) then
            J0 = eval_J(Uvec, "maxdLqdt")
            gradJ1 = GradL2ForLq(Uvec, lebesgueQ)
            CALL kappa_test(Uvec, gradJ1, J0, "maxdLqdt")   ! The function kappa_test_pert is not defined properly, on May 4, 2017
         end if

         DO WHILE ( (ABS(deltaJ) > OPTIM_TOL) .AND. (iter<MAX_ITER) )
            iter = iter + 1

            !======================================================
            ! CALCULATE BASIC GRADIENT IN THE H^s TOPOLOGY
            !======================================================
            if(rank==0 .and. verboseOptimization) print*, "calculating gradient"
            gradJ1 = GradL2ForLq(Uvec, lebesgueQ)
            CALL SobolevGradient(gradJ1, (3.0_pr*lebesgueQ-1.0_pr)/(2.0_pr*lebesgueQ))


            !======================================================
            ! CALCULATE Normal Of Tangent Space
            !======================================================
            if(rank==0 .and. verboseOptimization) print*, "calculating normal of tangent"
            unit_normal = calcConstraintDerivativeL2(Uvec, lebesgueQ)     ! unit_normal = q |u|^{q-2} u = nabla( ||u||_q^q )
            call sobolevGradient(unit_normal, (3.0_pr*lebesgueQ-1.0_pr)/(2.0_pr*lebesgueQ))  ! nabla^(H^...) (||u||_q^q)
            inner = global_summed_field_inner_product(unit_normal,unit_normal,"H_l^((3q-1)/(2q))")
            norm = sqrt(inner)
            unit_normal(:,:,:,:) = unit_normal(:,:,:,:)/norm               ! n = nabla( ||u||_q^q )/||nabla( ||u||_q^q )||_{H_l^...}


            !======================================================
            ! PROJECT to tangent space, gradJ1 = Projection_Tang(u) (nabla^{H^\dots} J)
            !======================================================
            if(rank==0 .and. verboseOptimization) print*, "projecting"
            inner = global_summed_field_inner_product(gradJ1,unit_normal,"H_l^((3q-1)/(2q))")                  ! < nabla J, n >_{H_l^...}
            gradJ1(:,:,:,:) = gradJ1(:,:,:,:) - inner*unit_normal(:,:,:,:)                                     ! nabla J = nabla J - < nabla J, n >_{H_l^...} n
            call div_free(gradJ1)                                                                              ! project average free, div free
            

            !======================================================
            ! Calculate Momentum Term (Polak-Ribière)
            !======================================================
            if(rank==0 .and. verboseOptimization) print*, "momentum term"
            beta = global_summed_field_inner_product(gradJ1,gradJ1-vecTransported_GradJ0,"H_l^((3q-1)/(2q))")
            beta = beta/global_summed_field_inner_product(gradJ1,gradJ1,"H_l^((3q-1)/(2q))")
            
            
            !======================================================
            ! Descend Direction
            !======================================================
            if(rank==0 .and. verboseOptimization) print*, "direction"
            d1(:,:,:,:) = gradJ1(:,:,:,:) + beta*vecTransported_d0

            !===========================================
            ! GET DIAGNOSTICS OF ASCENT DIRECTION AND SAVE
            !===========================================
            if(rank==0 .and. verboseOptimization) print*, "get diagnostics"
            local_field_L2norm = Energy(gradJ1)
            CALL MPI_ALLREDUCE(local_field_L2norm, gradJ_K1, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
            CALL divergence(gradJ1, f_scalar)
            local_scalar_L2norm = inner_product(f_scalar, f_scalar, "L2")
            CALL MPI_ALLREDUCE(local_scalar_L2norm, divGradJ_L2, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
            
            !======================================
            ! FIND OPTIMAL tau BY ARC OPTIMIZATION
            !====================================== 
            if(rank==0 .and. verboseOptimization) print*, "find tau"

            CALL optim_msg_handle(20) 
            tau_brack = mnbrak(iter, "maxdLqdt", Uvec, d1, 0.0_pr, tau0, mnbrak_flag)
            if(rank==0 .and. tauDebugToConsole) print*, achar(9), "iter", iter, "tau_brack",tau_brack

            IF (mnbrak_flag /= 0) THEN
               IF (rank==0) THEN
                  print*, "mnbrak error"
                  select case (mnbrak_flag)
                     CASE (1)
                        print*, achar(9), "Going uphill... Verify gradient!"   !achar(9) = ascii code 9 = tab
                     CASE (2)
                        print*, achar(9), "Could not bracket minimum."
                     CASE (3)
                        print*, achar(9), "Decreasing tau..."
                  end select
               END IF
               CALL optim_error_handle(mnbrak_flag)            
               IF (save_diag_Optim) THEN
                  IF (rank==0) THEN
                     CALL save_diagnostics_optim("maxdLqdt", iter, tau1, beta, J0, K, E, divU_L2, dLqdt(1), dLqdt(2), dLqdt(3))   
                  END IF
               END IF
               RETURN
            ELSE
               CALL optim_msg_handle(21)
            END If
            CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

            CALL optim_msg_handle(30)

            tau1 = brent(iter, "maxdLqdt", Uvec, d1, tau_brack)    ! I add the new variable iter
            if(rank==0 .and. tauDebugToConsole) print*, achar(9), "iter", iter, "used tau (tau1)",tau1
            CALL optim_msg_handle(31)


            !IF (rank==0) print*, "tau1", tau1
            !tau1 = MIN(tau1, tau_max)                                  ! TAU_MAX = 10.0_pr
            !IF (tau1 == TAU_MAX) CALL optim_msg_handle(32)


            !======================================
            ! update vector transport for next step with "old" tangentspace at "old" u
            !======================================
            if(rank==0 .and. verboseOptimization) print*, "update vector transport"
            vecTransported_GradJ0 = vectorTransport(Uvec, lebesgueQ, constraintB, tau1*d1, gradJ1)
            vecTransported_d0 = vectorTransport(Uvec, lebesgueQ, constraintB, tau1*d1, d1)

            !======================================
            ! update u
            !======================================
            if(rank==0 .and. verboseOptimization) print*, "update u"
            Uvec = Uvec + tau1*d1
            call rescaleLqNorm(uvec, lebesgueQ, constraintB)
            !test = calc_global_Lq_norm(Uvec, lebesgueQ)
            !if(rank==0) print*, "uvec Lq norm", test, "B", constraintB

            !====================================
            ! CALCULATE DIAGNOSTICS OF CONTROL; Compute related quantities' values of new velocity
            !====================================
            if(rank==0 .and. verboseOptimization) print*, "calc diagnostics"
            local_field_L2norm = Energy(Uvec)
            CALL MPI_ALLREDUCE(local_field_L2norm, K, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)

            CALL divergence(Uvec, f_scalar)
            local_scalar_L2norm = inner_product(f_scalar, f_scalar, "L2")
            CALL MPI_ALLREDUCE(local_scalar_L2norm, divU_L2, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo) 

            local_field_L2norm = Enstrophy(Uvec)
            CALL MPI_ALLREDUCE(local_field_L2norm, E, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)

            dLqdt = calc_dLqdt(Uvec, lebesgueQ)

            J1 = eval_J(Uvec, "maxdLqdt")
            deltaJ = (J1-J0)/ABS(J0)
            IF (deltaJ < -MACH_EPSILON) THEN      ! Change om March 30, 2017
                IF (save_diag_Optim) THEN
                    IF (rank==0) THEN
                        CALL save_diagnostics_optim("maxdLqdt", iter+1, tau1, beta, J0, K, E, divU_L2, dLqdt(1), dLqdt(2), dLqdt(3))   
                    END IF
                END IF
                CALL optim_msg_handle(0)
                if(rank==0) print*, "WARNING! Cost functional not increasing, exiting"
                !save_data_optim = .FALSE.
                EXIT
            ELSE
               !save_data_optim = .TRUE.
            END IF
            CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)     ! for everyone to know save_data_optim = .FALSE., otherwise next if barrier might wait indefinitely
                

            IF (save_diag_Optim) THEN
               if(rank==0 .and. verboseOptimization) print*, "saving diag optim"
               IF (rank==0) THEN
                  CALL save_diagnostics_optim("maxdLqdt", iter, tau1, beta, J1, K, E, divU_L2, dLqdt(1), dLqdt(2), dLqdt(3))   
               END IF
            END IF
            CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

            !===============================
            ! UPDATE OLD VARIABLES
            !===============================   
            if(rank==0 .and. verboseOptimization) print*, "update old variables"
            if(rank==0) print*, achar(9), "iter =", iter, "J1", J1, "tau", tau1, "deltaJ", deltaJ
            J0 = J1
            gradJ0 = gradJ1
            tau0 = tau1
            d0 = d1

         END DO


         if(rank==0) then
            if(iter<MAX_ITER) then
               print*, "optimization terminated successful after", iter, "iterations"
            if(rank==0) print*, achar(9), "iter final =", iter, "J1", J1, "tau", tau1, "deltaJ", deltaJ
            else
               print*, "optimization terminated by max iterations", iter, MAX_ITER
            end if
         end if


         IF (save_data_Optim) THEN
            if(rank==0 .and. verboseOptimization) print*, "saving data optim"
            CALL vel2vort(Uvec, Wvec)
            CALL diagnosticScalars(Uvec, Wvec, iter)

            IF (MOD(iter,1)==0) THEN
               CALL save_Ctrl(Uvec, Wvec, iter, "maxdLqdt")
            END IF

         END IF

         CALL optim_msg_handle(1)

         !===============================
         ! save results
         !===============================    
         call save_to_optimizationResultList(constraintB, J1, iter)
         



         DEALLOCATE(gradJ0)
         DEALLOCATE(gradJ1)
         deallocate(d0)
         deallocate(d1)
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
            local_J = calc_dLqdt(aux1, lebesgueQ)
            call mpi_allreduce(local_J, J, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
         CASE ("LineMin")
            CALL rescaleLqNorm(aux1, lebesgueQ, constraintB)
            J = -1.0_pr*eval_J(aux1, "maxdLqdt")
         case DEFAULT
            IF (rank==0)  print*, "WARNING, case ", mysystem, " for mysystem not found in eval_grad_J"
            J = 0.0_pr
      END SELECT
      deallocate(aux1)
   END FUNCTION eval_J

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
 

      if(mnbra_calcSaveAllJvalues) then
         DO WHILE (tB > MACH_EPSILON) 
            tB = CGOLD*tB
            phi_bar = phi + tB*grad
            FB = eval_J(phi_bar, "LineMin")
            FuncEval = FuncEval+1
            IF (saveLineMin) CALL save_linemin_data(tA, tB, tC, FA, FB, FC, iter, optimizationIter, mysystem, "append")
         END DO
         tB = MAX(tB0, MACH_EPSILON)
      end if


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


      call createDirectoryIfNonExistent(ConstraintDir//"tau_data")

      !filename = HomeDir//"/brent_info"//".dat"
      filename = ConstraintDir//"tau_data/"//"brent_info_"//itertxt//".dat"
      !filename = HomeDir//"/brent_info_maxdEdtHeli_Nu_E"//E0txt//"_IG"//IGtxt//"_brent_info.dat"
      OPEN(10, FILE = filename, FORM = 'FORMATTED', STATUS = 'REPLACE')
      WRITE(10, "(3 G20.12)") "#", "Tau", "J" 
      phi_bar = phi + D*grad
      FX = eval_J(phi_bar, "LineMin")
      WRITE(10, "(3 G20.12)") 0.0, D, FX
      CLOSE(10)
      
      phi_bar = phi + X*grad
      FX = eval_J(phi_bar, "LineMin")
 
      FV = FX 
      FW = FX


      DO j=1,ITMAX

         OPEN(10, FILE = filename, FORM = 'FORMATTED', STATUS = 'OLD', POSITION = 'APPEND')
         WRITE(10, "(3 G20.12)") j, X, FX
         CLOSE(10)

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

      DEALLOCATE( phi_bar )

   END FUNCTION brent


   !============================================
   ! Perform kappa_test
   !============================================
   SUBROUTINE kappa_test(phi, gradJ, J0, mysystem)
      USE global_variables
      USE function_ops
      USE data_ops
      USE mpi
      use fftwfunction
      IMPLICIT NONE

      !REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(IN) :: phi
      REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(IN) :: phi
      REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(IN) :: gradJ
      !REAL(pr), INTENT(IN) :: J0
      REAL(pr), INTENT(INOUT) :: J0
      CHARACTER(len=*), INTENT(IN) :: mysystem
      REAL(pr) :: myepsilon, kappa, deltaJ
      REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: phi_pert, phi_bar
      REAL(pr) :: myexp, J1
      REAL(pr), DIMENSION(1:3) :: dx, local_inner_prod, global_inner_prod
      INTEGER :: ii      
      CHARACTER(2) :: K0txt, E0txt, IGtxt
      CHARACTER(4) :: strLebesgueQ
      CHARACTER(3) :: tempStr
      character(50) :: phiPertText		! load random b/divfree sine/load te0080/...
      real(pr), dimension(1:n(1),1:n(2),1:local_N) :: testScalarField
      real(pr) :: testScalar


      complex(pr), dimension(1:n(1),1:n(2),1:local_n,1:3) :: testRemoveAfterwards,testRemoveAfterwards2
      real(pr), dimension(1:n(1),1:n(2),1:local_n) :: testRemoveAfterwardsCompX, testRemoveAfterwardsCompY, testRemoveAfterwardsCompZ
         


      IF (rank==0) THEN
         print*, "kappaTest"
         print*, "q", lebesgueQ
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
      call kappa_test_pert(phi_pert, phiPertText, -3.0_pr, 0.0_pr, 0.0_pr)
      
      call divergence(phi_pert,testScalarField)
      testScalar = inner_product(testScalarField,testScalarField,"L2")
      if(rank==0) then
      	print*, "|| nabla cdot phi_pert ||_2^2 = ", testScalar
      end if

      local_inner_prod = field_inner_product(gradJ, phi_pert, "L2")
      
      
      

      call calculateSaveSpectrum(phi, "phi")
      call calculateSaveSpectrum(phi_pert, trim("phiPert-"//phiPertText))
      call calculateSaveSpectrum(gradJ, "gradJ"//"_q-"//strLebesgueQ)


      CALL MPI_ALLREDUCE(local_inner_prod, global_inner_prod, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)

      !DO ii=1,21
      DO ii=1,30
         if (rank==0) then
               print*, "kappa test index", ii, "/30"
         end if
         
         myexp = -14.0_pr + 0.5_pr*REAL(ii,pr)
         myepsilon = 10.0_pr**myexp         
         
         phi_bar = phi + myepsilon*phi_pert

         J1 = eval_J(phi_bar, "maxdLqdt")
         
         kappa = (J1-J0)/(myepsilon*SUM(global_inner_prod))
         deltaJ = J1-J0  

         IF (rank==0) THEN
            CALL save_kappa_test(myepsilon, SUM(global_inner_prod), deltaJ, kappa, ii, "kappa_q-"//strLebesgueQ//"_"//mysystem//".dat")
         END IF 
         

         CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
      END DO

      if( rank == 0 ) then
         print*, "kmax", kmax
      end if

      DEALLOCATE(phi_pert)
      DEALLOCATE(phi_bar)

   END SUBROUTINE kappa_test

END MODULE
