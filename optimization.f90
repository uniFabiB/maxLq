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

         REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: gradJ0, gradJ1, diff_gradJ
         REAL(pr), DIMENSION(:,:,:),   ALLOCATABLE :: f_scalar

         REAL(pr) :: J0, J1, deltaJ, ell1, ell2, tau, beta, local_scalar_L2norm, divU_L2, divGradJ_L2
         REAL(pr), DIMENSION(1:3) :: local_field_L2norm, local_gradJK0, local_gradJK1, K, E, gradJ_K0, gradJ_K1
         REAL(pr), DIMENSION(1:2) :: tau_brack
         REAL(pr), DIMENSION(1:3) :: dLqdt   ! Modified on April 24, 2017

         INTEGER :: iter, gradType, mnbrak_flag, FixConstr_flag, kappaTestLebesgueQIter
         LOGICAL :: iterMAX = .FALSE.

         ALLOCATE( gradJ0(1:n(1),1:n(2),1:local_N,1:3) )
         ALLOCATE( gradJ1(1:n(1),1:n(2),1:local_N,1:3) )
         ALLOCATE( diff_gradJ(1:n(1),1:n(2),1:local_N,1:3) )
         ALLOCATE( f_scalar(1:n(1),1:n(2),1:local_N) )

         ell2 = lambda2
         ell1 = lambda1

         CALL Fix_Constr(Uvec, FixConstr_flag)
         IF (FixConstr_flag /= 0) THEN
            CALL optim_msg_handle(13)
            RETURN
         END IF
         !====================================
         ! CALCULATE DIAGNOSTICS OF CONTROL
         !====================================
         local_field_L2norm = Energy(Uvec)
         CALL MPI_ALLREDUCE(local_field_L2norm, K, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)

         CALL divergence(Uvec, f_scalar)
         local_scalar_L2norm = inner_product(f_scalar, f_scalar, "L2")
         CALL MPI_ALLREDUCE(local_scalar_L2norm, divU_L2, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo) 
   
         local_field_L2norm = Enstrophy(Uvec)
         CALL MPI_ALLREDUCE(local_field_L2norm, E, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)

         dLqdt = calc_dLqdt(Uvec, lebesgueQ)

         J0 = eval_J(Uvec, "maxdLqdt", "K0E0")
         J1 = 1.5_pr*J0
         deltaJ = ABS( (J1-J0)/J0 )    
         iter = 0

         IF (save_diag_Optim) THEN
            IF (rank==0) THEN
               CALL save_diagnostics_optim("maxdLqdt", iter, 0.0_pr, 0.0_pr, J0, K, E, divU_L2, dLqdt(1), dLqdt(2), dLqdt(3))   ! Modify April 24, 2017
            END IF
            CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
         END IF
         IF (save_data_Optim) THEN
            CALL vel2vort(Uvec, Wvec)
            CALL diagnosticScalars(Uvec, Wvec, iter)
            CALL save_Ctrl(Uvec, Wvec, iter, "maxdLqdt")      ! Save velocity and vorticity
         END IF

         !!! kappa test !!!
         if (.true.) then
            do kappaTestLebesgueQIter=1, size(lebesgueQlist)
               lebesgueQ = lebesgueQlist(kappaTestLebesgueQIter)
               J0 = eval_J(Uvec, "maxdLqdt", "K0E0")
               CALL eval_grad_J(Uvec, gradJ1, iter, "maxdLqdt")   ! Calculate L2 gradient of dEdtHeli, modified on April 24, 2017
               CALL kappa_test(Uvec, gradJ1, J0, "maxdLqdt")   ! The function kappa_test_pert is not defined properly, on May 4, 2017
            end do
         end if

         if (.false.) then
            IF (kappaTest) THen
               CALL eval_grad_J(Uvec, gradJ1, iter, "maxdLqdt")   ! Calculate L2 gradient of dEdtHeli, modified on April 24, 2017
               CALL kappa_test(Uvec, gradJ1, J0, "maxdLqdt")   ! The function kappa_test_pert is not defined properly, on May 4, 2017
               kappaTest = .FALSE.
            end if
         end if


         DO WHILE ( (ABS(deltaJ) > OPTIM_TOL) .AND. (iter<MAX_ITER) )
            IF (rank==0) THEN
               print *, "iter =", iter
            END IF


            !call test(iter)
            !======================================================
            ! CALCULATE BASIC GRADIENT IN THE H^s TOPOLOGY
            !======================================================
            CALL eval_grad_J(Uvec, gradJ1, iter, "maxdLqdt")   ! Calculate L2 gradient of dEdtHeli, modified on April 24, 2017

            !CALL SobolevGradient(gradJ1, 2, ell1, ell2)   ! I use order=2; it is "1" by Diago's code.

            !======================================================
            ! DECIDE DIRECTION OF INCREASE
            !====================================================== 
            !            IF (MOD(iter+1,2)/=0) THEN                ! change on March 16, 2017
            !               gradType = 0
            !            ELSE
            !               gradType = 2
            !            END IF

            !            IF (MOD(iter,10)==0) THEN                  ! change on May 24, 2017
            !               gradType = 0
            !            ELSE
            !               gradType = 2
            !            END IF

            gradType = 0

            SELECT CASE (gradType)
               CASE (0) !REGULAR STEEPEST ASCENT METHOD
                  local_field_L2norm = 2.0_pr*Energy(gradJ1)
                  CALL MPI_ALLREDUCE(local_field_L2norm, gradJ_K1, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
                  beta = 0.0_pr
                  gradJ_K0 = gradJ_K1
                  diff_gradJ = gradJ1
   
               CASE (1) !FLETCHER-REEVES METHOD
                  local_field_L2norm = 2.0_pr*Energy(gradJ1)
                  CALL MPI_ALLREDUCE(local_field_L2norm, gradJ_K1, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
                  beta = SUM(gradJ_K1)/SUM(gradJ_K0)
                  gradJ_K0 = gradJ_K1
                  diff_gradJ = gradJ1
                  gradJ1 = gradJ1 + beta*gradJ0
    
               CASE (2) !POLAK-RIBIERE METHOD
                  diff_gradJ = gradJ1 - diff_gradJ
                  local_field_L2norm = field_inner_product(gradJ1, diff_gradJ, "L2")
                  CALL MPI_ALLREDUCE(local_field_L2norm, gradJ_K1, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
                  beta = SUM(gradJ_K1)/SUM(gradJ_K0)
                  local_field_L2norm = 2.0_pr*Energy(gradJ1)
                  CALL MPI_ALLREDUCE(local_field_L2norm, gradJ_K0, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
                  diff_gradJ = gradJ1
                  gradJ1 = gradJ1 + beta*gradJ0
   
            END SELECT

            !===========================================
            ! GET DIAGNOSTICS OF ASCENT DIRECTION AND SAVE
            !===========================================
            local_field_L2norm = Energy(gradJ1)
            CALL MPI_ALLREDUCE(local_field_L2norm, gradJ_K1, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
            CALL divergence(gradJ1, f_scalar)
            local_scalar_L2norm = inner_product(f_scalar, f_scalar, "L2")
            CALL MPI_ALLREDUCE(local_scalar_L2norm, divGradJ_L2, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
            
            !======================================
            ! FIND OPTIMAL tau BY ARC OPTIMIZATION
            !====================================== 
            IF (iter==0) THEN 
               tau = alpha0
            END IF
            tau_brack(1) = 0.0_pr

            CALL optim_msg_handle(20) 
            tau_brack = mnbrak("maxdLqdt", Uvec, gradJ1, tau_brack(1), tau, mnbrak_flag)

            IF (mnbrak_flag /= 0) THEN
               IF (rank==0) THEN
                  print *, "Brent iteration beyond maximum, the maxdLqdt stops iterating ... "
               END IF
               CALL optim_error_handle(mnbrak_flag)            
               IF (save_diag_Optim) THEN
                  IF (rank==0) THEN
                     CALL save_diagnostics_optim("maxdLqdt", iter, tau, beta, J0, K, E, divU_L2, dLqdt(1), dLqdt(2), dLqdt(3))   
                  END IF
                  CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
               END IF
               RETURN
            ELSE
               CALL optim_msg_handle(21)
            END IF

            CALL optim_msg_handle(30)

            tau = brent(iter, "maxdLqdt", Uvec, gradJ1, tau_brack)    ! I add the new variable iter
            CALL optim_msg_handle(31)
            tau = MIN(tau, TAU_MAX)                                  ! TAU_MAX = 10.0_pr

            !======================================
            ! UPDATE CONTROL VARIABLE
            !======================================
            IF (tau == TAU_MAX) THEN
               CALL optim_msg_handle(32)
            END IF

            Uvec = Uvec + tau*gradJ1 
            CALL Fix_Constr(Uvec, FixConstr_flag)
            IF (FixConstr_flag /= 0) THEN
               CALL optim_msg_handle(13)
               RETURN
            END IF

            !====================================
            ! CALCULATE DIAGNOSTICS OF CONTROL; Compute related quantities' values of new velocity
            !====================================
            local_field_L2norm = Energy(Uvec)
            CALL MPI_ALLREDUCE(local_field_L2norm, K, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)

            CALL divergence(Uvec, f_scalar)
            local_scalar_L2norm = inner_product(f_scalar, f_scalar, "L2")
            CALL MPI_ALLREDUCE(local_scalar_L2norm, divU_L2, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo) 

            local_field_L2norm = Enstrophy(Uvec)
            CALL MPI_ALLREDUCE(local_field_L2norm, E, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)

            dLqdt = calc_dLqdt(Uvec, lebesgueQ)

            J1 = eval_J(Uvec, "maxdLqdt", "K0E0")
            deltaJ = (J1-J0)/ABS(J0)
            IF (deltaJ < -1E-18) THEN      ! Change om March 30, 2017
                IF (save_diag_Optim) THEN
                    IF (rank==0) THEN
                        CALL save_diagnostics_optim("maxdLqdt", iter+1, tau, beta, J0, K, E, divU_L2, dLqdt(1), dLqdt(2), dLqdt(3))   
                    END IF
                    CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
                END IF
                CALL optim_msg_handle(0)
                save_data_optim = .FALSE.
                EXIT
            ELSE
               save_data_optim = .TRUE.
            END IF

            !===============================
            ! UPDATE VARIABLES
            !===============================     
            J0 = J1
            gradJ0 = gradJ1
            iter = iter + 1

            !call testFFT(gradJ0(:,:,:,1))

            IF (save_diag_Optim) THEN
               IF (rank==0) THEN
                  CALL save_diagnostics_optim("maxdLqdt", iter, tau, beta, J1, K, E, divU_L2, dLqdt(1), dLqdt(2), dLqdt(3))   
               END IF
               CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
            END IF
  
            IF (save_data_Optim) THEN
               CALL vel2vort(Uvec, Wvec)
               CALL diagnosticScalars(Uvec, Wvec, iter)

               IF (MOD(iter,1)==0) THEN
                  CALL save_Ctrl(Uvec, Wvec, iter, "maxdLqdt")
               END IF

            END IF
 
         END DO

         CALL optim_msg_handle(1)

         DEALLOCATE(gradJ0)
         DEALLOCATE(gradJ1)
         DEALLOCATE(diff_gradJ)
         DEALLOCATE(f_scalar)

      end subroutine
      
   !================================================= 
   ! SUBROUTINE THAT PROJECTS A GIVEN FIELD ONTO THE 
   ! CONSTRAINT MANIFOLD (CO-DIMENSION 1 or 2)
   !=================================================
   SUBROUTINE Fix_Constr(myfield, myflag)
       USE global_variables
       USE function_ops
       USE data_ops
       USE mpi
       IMPLICIT NONE

       REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(INOUT) :: myfield
       INTEGER, INTENT(OUT) :: myflag
            !       INTEGER :: optim_msg_flag

       SELECT CASE (ConsType)
          CASE (1)
             CALL div_free(myfield)
             CALL Fix_E0(myfield)
             myflag = 0 
          CASE (2)
            !             CALL optim_msg_handle(10)
            !             CALL Fix_K0E0(myfield, optim_msg_flag)
            !             CALL optim_msg_handle(optim_msg_flag)
            !             IF (optim_msg_flag == 11) THEN
            !                myflag = 0
            !             ELSE
            !                myflag = 1
            !             END IF
       END SELECT
   END SUBROUTINE Fix_Constr

   !=================================================
   ! FUNCTION THAT CALCULATES COST FUNCTIONAL
   !=================================================
   RECURSIVE FUNCTION eval_J(myfield, mysystem, myJ) RESULT (J)
      USE global_variables
      USE data_ops
      USE fftwfunction
      USE function_ops
      USE mpi
      IMPLICIT NONE

      REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(IN) :: myfield
      CHARACTER(len=*), INTENT(IN) :: mysystem
      CHARACTER(len=*), INTENT(IN) :: myJ
      REAL(pr) :: J

      REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: phi, aux1, aux2, aux3, aux4, aux5   ! Newly added aux5, April 24, 2017
      REAL(pr), DIMENSION(1:3) :: local_VecField_norm, global_VecField_norm
      REAL(pr) :: local_Helicity, global_Helicity                                      ! Newly added on May 2, 2017

      real(pr) :: local_J              ! added Nov 5, 2024
      
      INTEGER :: constr_flag

      J = 0.0_pr

      SELECT CASE (mysystem)
         case ("maxdLqdt")
            allocate( aux1(1:n(1),1:n(2),1:local_N,1:3) )
            aux1 = myfield
            !call div_free(aux1)
            local_J = calc_dLqdt(aux1, lebesgueQ)

            !print*, "eval_j local_j", local_j
            call mpi_allreduce(local_J, J, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
            deallocate(aux1)

               !         CASE ("FixK0E0")
               !            SELECT CASE (myJ)
               !               CASE ("Enstrophy")
               !                  local_vecField_norm = Enstrophy(myfield)
               !                  CALL MPI_ALLREDUCE(local_VecField_norm, global_VecField_norm, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)                 
               !
               !                  J = 0.5_pr*( ( SUM(global_VecField_norm) - E0 )/E0 )**2
               !
               !               CASE ("Energy")
               !                  local_VecField_norm = Energy(myfield) 
               !                  CALL MPI_ALLREDUCE(local_VecField_norm, global_VecField_norm, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)                 
               !
               !                  J = 0.5_pr*( ( SUM(global_VecField_norm) - K0 )/K0 )**2  

               !              END SELECT 
             
         CASE ("LineMin")
            ALLOCATE( aux1(1:n(1),1:n(2),1:local_N,1:3) )
            aux1 = myfield
            SELECT CASE (myJ)            ! Need further modification
                  !               CASE ("FixK0E0")

                  !                 CALL div_free(aux1)
                  !                 CALL Fix_K0(aux1)
                  !                 J = eval_J(aux1, "FixK0E0", "Enstrophy")
                  
               CASE ("maxdLqdt")

                  CALL div_free(aux1)
                  CALL Fix_Constr(aux1, constr_flag)
                  IF (constr_flag == 0) THEN
                     J = -1.0_pr*eval_J(aux1, "maxdLqdt", "K0E0")
                  ELSE
                     CALL optim_msg_handle(14)
                     J = -1.0_pr*eval_J(aux1, "maxdLqdt", "K0E0")
                  END IF

            END SELECT
         DEALLOCATE(aux1)
      END SELECT
   END FUNCTION eval_J
 
   !================================================
   ! Obtain gradient of cost functional
   !================================================
   SUBROUTINE eval_grad_J(phi, gradJ, iter_index, myJ)
      USE global_variables
      USE data_ops
      USE fftwfunction
      USE function_ops
      USE mpi
      IMPLICIT NONE

      REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(IN) :: phi
      REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(OUT) :: gradJ
      INTEGER, INTENT(IN) :: iter_index
      CHARACTER(len=*), INTENT(IN) :: myJ

      !REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: aux1, aux2, aux3, aux4, aux5   ! Newly added aux5, April 24, 2017
      !REAL(pr), DIMENSION(3) :: local_VecField_norm, global_VecField_norm

      !INTEGER :: ii

      !print* , myJ
      SELECT CASE (myJ)
         case ("maxdLqdt")
            gradJ = GradL2ForLq(phi, lebesgueQ)


!         CASE ("Energy")

!         CASE ("Enstrophy")
!            DO ii=1,3
!               gradJ(:,:,:,ii) = phi(:,:,:,ii)
!            END DO
!            local_VecField_norm = Enstrophy(phi)
!            CALL MPI_ALLREDUCE(local_VecField_norm, global_VecField_norm, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)                 
!
!            CALL div_free(gradJ)
!            CALL laplacian(gradJ)
!
!            gradJ = (SUM(global_VecField_norm) - E0)*gradJ/E0**2

      END SELECT      

   END SUBROUTINE eval_grad_J

   !==================================================
   ! BRACKET THE LOCATION OF OPTIMAL TAU
   !==================================================

   FUNCTION mnbrak(mysystem, phi, grad, tA0, tB0, myflag) RESULT (tau_brack)
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
      INTEGER :: FuncEval, iter 
      LOGICAL :: saveLineMin

      saveLineMin = .TRUE.

      ALLOCATE( phi_bar(1:n(1),1:n(2),1:local_N,1:3) )

      FuncEval = 0
      iter = 0
 
      tA = tA0
      tB = MAX(tB0, MACH_EPSILON)

      phi_bar = phi + tA*grad
      FA = eval_J(phi_bar, "LineMin", mysystem)
      FuncEval = FuncEval+1

      phi_bar = phi + tB*grad
      FB = eval_J(phi_bar, "LineMin", mysystem)
      FuncEval = FuncEval+1

      IF (saveLineMin) CALL save_linemin_data(tA, tB, tC, FA, FB, FC, iter, mysystem, "replace")
 
      DO WHILE (FB > FA .AND. tB > MACH_EPSILON) 
         tB = CGOLD*tB
         phi_bar = phi + tB*grad
         FB = eval_J(phi_bar, "LineMin", mysystem)
         FuncEval = FuncEval+1
         IF (saveLineMin) CALL save_linemin_data(tA, tB, tC, FA, FB, FC, iter, mysystem, "append")
      END DO

      IF (tB .LE. MACH_EPSILON) THEN
         myflag = 1
         RETURN
      END IF

      tC = GOLD*tB
      phi_bar = phi + tC*grad
      FC = eval_J(phi_bar, "LineMin", mysystem)
      FuncEval = FuncEval+1

      IF (saveLineMin) CALL save_linemin_data(tA, tB, tC, FA, FB, FC, iter, mysystem, "append")

      DO WHILE (FB>=FC .AND. iter<ITMAX)
         iter = iter+1
         tC = GOLD*tC
         phi_bar = phi + tC*grad
         FC = eval_J(phi_bar, "LineMin", mysystem)
         FuncEval = FuncEval+1
                  
         R = (tB-tA)*(FB-FC)
         Q = (tB-tC)*(FB-FA)
         tP = tB - 0.5_pr*((tB-tC)*Q - (tB-tA)*R)/( SIGN( MAX(ABS(Q-R),MACH_EPSILON), Q-R) )
    
         Pmax = tB + GLIMIT*(tC-tB)
    
         IF ( (tB-tP)*(tP-tC)>0 ) THEN
            phi_bar = phi + tP*grad
            FP = eval_J(phi_bar, "LineMin", mysystem)
        
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
            FP = eval_J(phi_bar, "LineMin", mysystem)
        
         ELSEIF ( (tC-tP)*(tP-Pmax)>0 ) THEN
            phi_bar = phi + tP*grad
            FP = eval_J(phi_bar, "LineMin", mysystem)
       
            IF (FP<FC) THEN
               tB = tC
               tC = tP
               FB = FC
               FC = FP
               tP = tC+GOLD*(tC-tB)
               phi_bar = phi + tP*grad
               FP = eval_J(phi_bar, "LineMin", mysystem)
            END IF
        
        ELSEIF ( (tP-Pmax)*(Pmax-tC)>=0 ) THEN
            tP = Pmax
            phi_bar = phi + tP*grad
            FP = eval_J(phi_bar, "LineMin", mysystem)
        
         ELSE
            tP = tC + GOLD*(tC-tB)
            phi_bar = phi + tP*grad
            FP = eval_J(phi_bar, "LineMin", mysystem)
         END IF
    
         tA = tB
         tB = tC
         tC = tP
         FA = FB
         FB = FC
         FC = FP
        
         IF (saveLineMin) CALL save_linemin_data(tA, tB, tC, FA, FB, FC, iter, mysystem, "append")
 
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
      CHARACTER(5) :: itertxt
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

      filename = HomeDir//"/brent_info_maxdEdtHeli_Nu_E"//E0txt//"_IG"//IGtxt//"_brent_info.dat"
      OPEN(10, FILE = filename, FORM = 'FORMATTED', STATUS = 'REPLACE')
      WRITE(10,*) "# Iter  Tau" 
      phi_bar = phi + D*grad
      FX = eval_J(phi_bar, "LineMin", mysystem)
      WRITE(10, "(G20.12, G20.12)") D, FX
      CLOSE(10)
      
      phi_bar = phi + X*grad
      FX = eval_J(phi_bar, "LineMin", mysystem)
 
      FV = FX 
      FW = FX


      DO j=1,ITMAX

         OPEN(10, FILE = filename, FORM = 'FORMATTED', STATUS = 'OLD', POSITION = 'APPEND')
         WRITE(10, "(G20.12, G20.12)") X, FX
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
         FU = eval_J(phi_bar, "LineMin", mysystem)
 
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
      REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(INOUT) :: phi
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


      complex(pr), dimension(1:n(1),1:n(2),1:local_n,1:3) :: testRemoveAfterwards,testRemoveAfterwards2
      real(pr), dimension(1:n(1),1:n(2),1:local_n) :: testRemoveAfterwardsCompX, testRemoveAfterwardsCompY, testRemoveAfterwardsCompZ
         


      IF (rank==0) THEN
         print*, "kappaTest"
         print*, "q", lebesgueQ
      END If
      
      WRITE(K0txt,'(i2.2)') K0_index
      WRITE(E0txt,'(i2.2)') E0_index
      WRITE(IGtxt,'(i2.2)') iguess
      ALLOCATE( phi_pert(1:n(1),1:n(2),1:local_N,1:3) )
      ALLOCATE( phi_bar(1:n(1),1:n(2),1:local_N,1:3) )


      !CALL kappa_test_pert(phi_pert, "divfree sine", 1.0_pr, 2.0_pr, 4.0_pr)
      !CALL kappa_test_pert(phi_pert, "random smooth", 1.0_pr, 2.0_pr, 4.0_pr)
      !CALL kappa_test_pert(phi_pert, "random", 0.0_pr, 0.0_pr, 0.0_pr)
      CALL kappa_test_pert(phi_pert, "load te0080", 0.0_pr, 0.0_pr, 0.0_pr)
          
      local_inner_prod = field_inner_product(gradJ, phi_pert, "L2")

      CALL MPI_ALLREDUCE(local_inner_prod, global_inner_prod, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)

      DO ii=1,21
         if (rank==0) then
               print*, "kappa test index", ii, "/21"
         end if
         
         myexp = -14.0_pr + 0.5_pr*REAL(ii,pr)
         myepsilon = 10.0_pr**myexp         
         
         phi_bar = phi + myepsilon*phi_pert


         !!! TEST REMOVE AFTERWARDS !!!
         ! output the fourier transforms of phi and phi_pert
         testRemoveAfterwards = dcmplx(phi)
         call fftfwdv(testRemoveAfterwards,testRemoveAfterwards2)
         testRemoveAfterwardsCompX = real(testRemoveAfterwards2(:,:,:,1),pr)
         testRemoveAfterwardsCompY = real(testRemoveAfterwards2(:,:,:,2),pr)
         testRemoveAfterwardsCompZ = real(testRemoveAfterwards2(:,:,:,3),pr)
         CALL save_field_R3toR3_ncdf(testRemoveAfterwardsCompX, testRemoveAfterwardsCompY, testRemoveAfterwardsCompZ, "Ux", "Uy", "Uz", HomeDir//"phi.nc", "netCDF")
         testRemoveAfterwards = dcmplx(phi_pert)
         call fftfwdv(testRemoveAfterwards,testRemoveAfterwards2)
         testRemoveAfterwardsCompX = real(testRemoveAfterwards2(:,:,:,1),pr)
         testRemoveAfterwardsCompY = real(testRemoveAfterwards2(:,:,:,2),pr)
         testRemoveAfterwardsCompZ = real(testRemoveAfterwards2(:,:,:,3),pr)
         CALL save_field_R3toR3_ncdf(testRemoveAfterwardsCompX, testRemoveAfterwardsCompY, testRemoveAfterwardsCompZ, "Ux", "Uy", "Uz", HomeDir//"phi_pert.nc", "netCDF")
         !!! TEST REMOVE AFTERWARDS !!!




         J1 = eval_J(phi_bar, mysystem, "maxdLqdt")
         
         kappa = (J1-J0)/(myepsilon*SUM(global_inner_prod))
         deltaJ = J1-J0  

         IF (rank==0) THEN
            WRITE(tempStr,'(F3.1)') lebesgueQ-int(lebesgueQ)              ! might result in rounding errors for 1.999999999999
            WRITE(strLebesgueQ,'(I2.2,a2)') int(lebesgueQ), tempStr(2:)
            CALL save_kappa_test(myepsilon, SUM(global_inner_prod), deltaJ, kappa, ii, "KappaTest/"//"q-"//strLebesgueQ//"_"//mysystem//".dat")
         END IF 

         CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
      END DO

      DEALLOCATE(phi_pert)
      DEALLOCATE(phi_bar)

   END SUBROUTINE kappa_test

END MODULE
