module function_ops
  
   implicit none
   CONTAINS

      !===================================
      !--Initial guess
      !===================================
      SUBROUTINE initial_guess
                        !         use, intrinsic :: iso_c_binding                           ! Added on March 19, 2017
         USE global_variables
         USE data_ops
         use fftwfunction
         USE mpi
         IMPLICIT NONE

         REAL(pr), DIMENSION(1:3) :: dx, local_K, K
         INTEGER :: ii, jj, kk, nn
         INTEGER, DIMENSION(:), ALLOCATABLE :: seed
         REAL(pr) :: X, Y, Z
         LOGICAL :: read_from_file
         CHARACTER(200) :: filename
         CHARACTER(2) :: K0txt, E0txt, IGtxt, Fx_txt, Fy_txt, Fz_txt
                        !         CHARACTER(2) :: WEIGHTtxt                           ! WEIGHTtxt is newly added on April 24, 2017
         CHARACTER(4) :: Ntxt
         REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: Ux, Uy, Uz
         COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: Ux_hat, Uy_hat, Uz_hat
         COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: aux, faux
         real(pr), dimension(:,:,:,:), allocatable :: auxVec
         logical :: loadSuccessful
         
         complex(pr), dimension(:,:,:,:), allocatable :: auxVec2, fauxVec2
         real(pr) :: norm_k

                        !         ALLOCATE( aux(1:n(1),1:n(2),1:local_N) ) 
                        !         ALLOCATE( faux(1:n(1),1:n(2),1:local_N) )
         ALLOCATE( Ux(1:n(1),1:n(2),1:local_N) )
         ALLOCATE( Uy(1:n(1),1:n(2),1:local_N) )
         ALLOCATE( Uz(1:n(1),1:n(2),1:local_N) )


         dx = 1.0_pr/REAL(n, pr)

         SELECT CASE (iguess)
            case (0)
               filename = inputDir//"FRT_N256E500T017_Uvec_fwdTE0220.nc"
               CALL read_field_R3toR3_ncdf(Uvec, filename, "Ux", "Uy", "Uz", loadSuccessful)
               !x = 2.0_pr*PI*REAL(n(1),pr)/16.0_pr     ! wave number from which on to cut fourier modes
               !call divAvg_free(uvec)

            case (1)
               allocate( auxVec(1:n(1),1:n(2),1:local_N,1:3) )

               call kappa_test_pert(uvec,"sine",4.0_pr,2.0_pr,5.0_pr, loadSuccessful)
               call kappa_test_pert(auxvec,"sine",2.0_pr,3.0_pr,1.0_pr, loadSuccessful)
               uvec(:,:,:,:) = uvec(:,:,:,:) + auxvec(:,:,:,:)
               call kappa_test_pert(auxvec,"sine",7.0_pr,2.0_pr,6.0_pr, loadSuccessful)
               uvec(:,:,:,:) = uvec(:,:,:,:) + auxvec(:,:,:,:)
               call kappa_test_pert(auxvec,"sine",8.0_pr,1.0_pr,2.0_pr, loadSuccessful)
               uvec(:,:,:,:) = uvec(:,:,:,:) + auxvec(:,:,:,:)
               call kappa_test_pert(auxvec,"sine",1.0_pr,9.0_pr,0.0_pr, loadSuccessful)
               uvec(:,:,:,:) = uvec(:,:,:,:) + auxvec(:,:,:,:)
               call kappa_test_pert(auxvec,"sine",5.0_pr,7.0_pr,2.0_pr, loadSuccessful)
               uvec(:,:,:,:) = uvec(:,:,:,:) + auxvec(:,:,:,:)
               call kappa_test_pert(auxvec,"sine",1.0_pr,2.0_pr,2.0_pr, loadSuccessful)
               uvec(:,:,:,:) = uvec(:,:,:,:) + auxvec(:,:,:,:)
               call kappa_test_pert(auxvec,"sine",2.0_pr,4.0_pr,3.0_pr, loadSuccessful)
               uvec(:,:,:,:) = uvec(:,:,:,:) + auxvec(:,:,:,:)
               call kappa_test_pert(auxvec,"sine",7.0_pr,2.0_pr,3.0_pr, loadSuccessful)
               uvec(:,:,:,:) = uvec(:,:,:,:) + auxvec(:,:,:,:)
               call divAvg_free(uvec)
               
               deallocate( auxVec )


            case (2)
               call kappa_test_pert(uvec, "load-random-a", 0.0_pr, 0.0_pr, 0.0_pr, loadSuccessful)

            case (3)
               call kappa_test_pert(uvec, "load-random-smooth-a", 0.0_pr, 0.0_pr, 0.0_pr, loadSuccessful)
               
            case (4)
               call kappa_test_pert(uvec, "load-random-exp-a", 0.0_pr, 0.0_pr, 0.0_pr, loadSuccessful)
            
            case (5)
               call kappa_test_pert(uvec, "load-random-poly-a", -3.0_pr, 0.0_pr, 0.0_pr, loadSuccessful)

            case (6)
               call kappa_test_pert(uvec, "load-k-random-a", 1000.0_pr, 0.0_pr, 0.0_pr, loadSuccessful)

            case (7)
               call kappa_test_pert(uvec, "save-random", 0.0_pr, 0.0_pr, 0.0_pr, loadSuccessful)

            CASE (9)                                 ! Can be used when recover from the terminated code
               !filename = "/work/yund0050/maxdEdtHeli_100_06/3_005_WEIGHT100_N0256_E37_IG10_DoubleResolution_u0.nc"                           ! Added on March 24, 2017, only work once
               !filename = inputDir//"u_result_B32_0512.nc"
               filename = inputDir//loadTempFunctionName
               CALL read_field_R3toR3_ncdf(Uvec, filename, "Ux", "Uy", "Uz", loadSuccessful)
               
            CASE (50)                                   ! Arnold-Beltrami-Childress (ABC) flow
               DO kk=1, local_N
                  DO jj=1, n(2)
                     DO ii=1, n(1)
                        X = REAL(ii-1,pr)*dx(1)
                        Y = REAL(jj-1,pr)*dx(2)
                        Z = REAL(kk+local_k_offset-1,pr)*dx(3)
                        !                        Z = REAL(kk-1,pr)*dx(3)
                        Ux(ii,jj,kk) = 1.0_pr/1.0_pr*(cos(2.0_pr*PI*Z) + cos(2.0_pr*PI*Y))
                        Uy(ii,jj,kk) = 1.0_pr/1.0_pr*(cos(2.0_pr*PI*X) + cos(2.0_pr*PI*Z))
                        Uz(ii,jj,kk) = 1.0_pr/1.0_pr*(cos(2.0_pr*PI*Y) + cos(2.0_pr*PI*X))
                     END DO
                  END DO
               END DO
               uvec(:,:,:,1) = ux(:,:,:)
               uvec(:,:,:,2) = uy(:,:,:)
               uvec(:,:,:,3) = uz(:,:,:)

            CASE (51)                                   ! Arnold-Beltrami-Childress (ABC) flow, plus perturbation
               DO kk=1, local_N
                  DO jj=1, n(2)
                     DO ii=1, n(1)
                        X = REAL(ii-1,pr)*dx(1)
                        Y = REAL(jj-1,pr)*dx(2)
                        Z = REAL(kk+local_k_offset-1,pr)*dx(3)
                        !                        Z = REAL(kk-1,pr)*dx(3)
                        Ux(ii,jj,kk) = 1.0_pr/1.0_pr*(cos(2.0_pr*PI*Z) + cos(2.0_pr*PI*Y))
                        Uy(ii,jj,kk) = 1.0_pr/1.0_pr*(cos(2.0_pr*PI*X) + cos(2.0_pr*PI*Z))
                        Uz(ii,jj,kk) = 1.0_pr/1.0_pr*(cos(2.0_pr*PI*Y) + cos(2.0_pr*PI*X))
                     END DO
                  END DO
               END DO

                                       ! ADD PERTURBATION TO (50)
               DO nn = 5,10
                  DO kk=1, local_N
                     DO jj=1, n(2)
                        DO ii=1, n(1)
                           X = REAL(ii-1,pr)*dx(1)
                           Y = REAL(jj-1,pr)*dx(2)
                           Z = REAL(kk+local_k_offset-1,pr)*dx(3)
                        !                           Z = REAL(kk-1,pr)*dx(3)
                           Ux(ii,jj,kk) = Ux(ii,jj,kk) + (1.0_pr/REAL(nn,pr)**2)*SIN(2.0_pr*PI*REAL(nn,pr)*X)*COS(2.0_pr*PI*REAL(nn,pr)*Y)*SIN(2.0_pr*PI*REAL(nn,pr)*Z)
                           Uy(ii,jj,kk) = Uy(ii,jj,kk) + (1.0_pr/REAL(nn,pr)**2)*COS(2.0_pr*PI*REAL(nn,pr)*X)*COS(2.0_pr*PI*REAL(nn,pr)*Y)*COS(2.0_pr*PI*REAL(nn,pr)*Z)
                           Uz(ii,jj,kk) = Uz(ii,jj,kk) + (1.0_pr/REAL(nn,pr)**2)*COS(2.0_pr*PI*REAL(nn,pr)*X)*SIN(2.0_pr*PI*REAL(nn,pr)*Y)*SIN(2.0_pr*PI*REAL(nn,pr)*Z)
                        END DO
                     END DO
                  END DO
               END DO

               Uvec(:,:,:,1) = Ux
               Uvec(:,:,:,2) = Uy
               Uvec(:,:,:,3) = Uz
               CALL divAvg_free(Uvec)   

            CASE (60)                              ! I added on Feb 23, 2017; Arnold-Beltrami-Childress (ABC) flow. Modified on March 16, 2017, add read file data; Discuss with Diego
               if (E0_index .eq. 12) then
                  DO kk=1, local_N 
                     DO jj=1, n(2)
                        DO ii=1, n(1)
                           X = REAL(ii-1,pr)*dx(1)
                           Y = REAL(jj-1,pr)*dx(2)
                           Z = REAL(kk+local_k_offset-1,pr)*dx(3)
                        !                           Z = REAL(kk-1,pr)*dx(3)
                           Ux(ii,jj,kk) = 1.0_pr/1.0_pr*(cos(2.0_pr*PI*Z) + cos(2.0_pr*PI*Y))
                           Uy(ii,jj,kk) = 1.0_pr/1.0_pr*(cos(2.0_pr*PI*X) + cos(2.0_pr*PI*Z))
                           Uz(ii,jj,kk) = 1.0_pr/1.0_pr*(cos(2.0_pr*PI*Y) + cos(2.0_pr*PI*X))
                        END DO
                     END DO
                  END DO
           
                  Uvec(:,:,:,1) = Ux
                  Uvec(:,:,:,2) = Uy
                  Uvec(:,:,:,3) = Uz
               else
                  WRITE(K0txt,'(i2.2)') K0_index
                  WRITE(E0txt,'(i2.2)') E0_index-1
                  WRITE(IGtxt,'(i2.2)') iguess
                        !                  WRITE(WEIGHTtxt, '(i2.2)') int_WEIGHT
                        !                  filename = "/work/yund0050/MultiObjective_095_01/WEIGHT"//WEIGHTtxt//"_E"//E0txt//"_maxdEdtHeli_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_u0.nc"
                  filename = HomeDir//"_E"//E0txt//"_IG"//IGtxt//"_u0.nc"                           ! Newly added on May 8, 2017
                  CALL read_field_R3toR3_ncdf(Uvec, filename, "Ux", "Uy", "Uz", loadSuccessful)

                  Ux = Uvec(:,:,:,1)
                  !CALL dealiasing(Ux) ! commented out by fb
                  Uvec(:,:,:,1) = Ux

                  Uy = Uvec(:,:,:,2)
                  !CALL dealiasing(Uy) ! commented out by fb
                  Uvec(:,:,:,2) = Uy

                  Uz = Uvec(:,:,:,3)
                  !CALL dealiasing(Uz) ! commented out by fb
                  Uvec(:,:,:,3) = Uz

                  CALL divAvg_free(Uvec)
               end if

            CASE (61)                              !Initial guess is a perturbation of Case(60), with start from obtained velocity
               if (E0_index .eq. 12) then                           ! Changed from .eq. 1
                  DO kk=1, local_N
                     DO jj=1, n(2)
                        DO ii=1, n(1)
                           X = REAL(ii-1,pr)*dx(1)
                           Y = REAL(jj-1,pr)*dx(2)
                           Z = REAL(kk+local_k_offset-1,pr)*dx(3)
                        !                           Z = REAL(kk-1,pr)*dx(3)
                           Ux(ii,jj,kk) = 1.0_pr/1.0_pr*(cos(2.0_pr*PI*Z) + cos(2.0_pr*PI*Y))
                           Uy(ii,jj,kk) = 1.0_pr/1.0_pr*(cos(2.0_pr*PI*X) + cos(2.0_pr*PI*Z))
                           Uz(ii,jj,kk) = 1.0_pr/1.0_pr*(cos(2.0_pr*PI*Y) + cos(2.0_pr*PI*Z))
                        END DO
                     END DO
                  END DO
             
                                          ! ADD PERTURBATION TO (60)
                  DO nn = 3,10
                     DO kk=1, local_N
                        DO jj=1, n(2)
                           DO ii=1, n(1)
                              X = REAL(ii-1,pr)*dx(1)
                              Y = REAL(jj-1,pr)*dx(2)
                              Z = REAL(kk+local_k_offset-1,pr)*dx(3)
                        !                              Z = REAL(kk-1,pr)*dx(3)
                              Ux(ii,jj,kk) = Ux(ii,jj,kk) + (1.0_pr/REAL(nn,pr)**2)*SIN(2.0_pr*PI*REAL(nn,pr)*X)*COS(2.0_pr*PI*REAL(nn,pr)*Y)*SIN(2.0_pr*PI*REAL(nn,pr)*Z)
                              Uy(ii,jj,kk) = Uy(ii,jj,kk) + (1.0_pr/REAL(nn,pr)**2)*COS(2.0_pr*PI*REAL(nn,pr)*X)*COS(2.0_pr*PI*REAL(nn,pr)*Y)*COS(2.0_pr*PI*REAL(nn,pr)*Z)
                              Uz(ii,jj,kk) = Uz(ii,jj,kk) + (1.0_pr/REAL(nn,pr)**2)*COS(2.0_pr*PI*REAL(nn,pr)*X)*SIN(2.0_pr*PI*REAL(nn,pr)*Y)*SIN(2.0_pr*PI*REAL(nn,pr)*Z)
                           END DO
                        END DO
                     END DO
                  END DO
               
                  Uvec(:,:,:,1) = Ux
                  Uvec(:,:,:,2) = Uy
                  Uvec(:,:,:,3) = Uz
                  CALL divAvg_free(Uvec)   
               else
                  WRITE(K0txt,'(i2.2)') K0_index
                  WRITE(E0txt,'(i2.2)') E0_index-1
                  WRITE(IGtxt,'(i2.2)') iguess
                        !                  WRITE(WEIGHTtxt, '(i2.2)') int_WEIGHT
                        !                  filename = "/work/yund0050/MultiObjective_095_01/WEIGHT"//WEIGHTtxt//"_E"//E0txt//"_maxdEdtHeli_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_u0.nc"
                  filename = HomeDir//"_E"//E0txt//"_IG"//IGtxt//"_u0.nc"                           ! Newly added on May 8, 2017
                  CALL read_field_R3toR3_ncdf(Uvec, filename, "Ux", "Uy", "Uz", loadSuccessful)

                  Ux = Uvec(:,:,:,1)
                  !CALL dealiasing(Ux) ! commented out by fb
                  Uvec(:,:,:,1) = Ux

                  Uy = Uvec(:,:,:,2)
                  !CALL dealiasing(Uy) ! commented out by fb
                  Uvec(:,:,:,2) = Uy

                  Uz = Uvec(:,:,:,3)
                  !CALL dealiasing(Uz) ! commented out by fb
                  Uvec(:,:,:,3) = Uz

                  CALL divAvg_free(Uvec)   
               end if

            CASE (62)                           !Initial guess is a perturbation of Case(60), with start from obtained velocity, plus perturbation
               if (E0_index .eq. 1) then
                  DO kk=1, local_N
                     DO jj=1, n(2)
                        DO ii=1, n(1)
                           X = REAL(ii-1,pr)*dx(1)
                           Y = REAL(jj-1,pr)*dx(2)
                           Z = REAL(kk+local_k_offset-1,pr)*dx(3)
                        !                           Z = REAL(kk-1,pr)*dx(3)
                           Ux(ii,jj,kk) = 1.0_pr/1.0_pr*(cos(2.0_pr*PI*Z) + cos(2.0_pr*PI*Y))
                           Uy(ii,jj,kk) = 1.0_pr/1.0_pr*(cos(2.0_pr*PI*X) + cos(2.0_pr*PI*Z))
                           Uz(ii,jj,kk) = 1.0_pr/1.0_pr*(cos(2.0_pr*PI*Y) + cos(2.0_pr*PI*Z))
                        END DO
                     END DO
                  END DO
             
                                          ! ADD PERTURBATION TO (60)
                  DO nn = 2,10
                     DO kk=1, local_N
                        DO jj=1, n(2)
                           DO ii=1, n(1)
                              X = REAL(ii-1,pr)*dx(1)
                              Y = REAL(jj-1,pr)*dx(2)
                              Z = REAL(kk+local_k_offset-1,pr)*dx(3)
                        !                              Z = REAL(kk-1,pr)*dx(3)
                              Ux(ii,jj,kk) = Ux(ii,jj,kk) + (1.0_pr/REAL(nn,pr)**2)*SIN(2.0_pr*PI*REAL(nn,pr)*X)*COS(2.0_pr*PI*REAL(nn,pr)*Y)*SIN(2.0_pr*PI*REAL(nn,pr)*Z)
                              Uy(ii,jj,kk) = Uy(ii,jj,kk) + (1.0_pr/REAL(nn,pr)**2)*COS(2.0_pr*PI*REAL(nn,pr)*X)*COS(2.0_pr*PI*REAL(nn,pr)*Y)*COS(2.0_pr*PI*REAL(nn,pr)*Z)
                              Uz(ii,jj,kk) = Uz(ii,jj,kk) + (1.0_pr/REAL(nn,pr)**2)*COS(2.0_pr*PI*REAL(nn,pr)*X)*SIN(2.0_pr*PI*REAL(nn,pr)*Y)*SIN(2.0_pr*PI*REAL(nn,pr)*Z)
                           END DO
                        END DO
                     END DO
                  END DO
               
                  Uvec(:,:,:,1) = Ux
                  Uvec(:,:,:,2) = Uy
                  Uvec(:,:,:,3) = Uz
                  CALL divAvg_free(Uvec)   
               else
                  WRITE(K0txt,'(i2.2)') K0_index
                  WRITE(E0txt,'(i2.2)') E0_index-1
                  WRITE(IGtxt,'(i2.2)') iguess
                        !                  WRITE(WEIGHTtxt, '(i2.2)') int_WEIGHT
                        !                  filename = "/work/yund0050/MultiObjective_095_01/WEIGHT"//WEIGHTtxt//"_E"//E0txt//"_maxdEdtHeli_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_u0.nc"
                  filename = HomeDir//"_E"//E0txt//"_IG"//IGtxt//"_u0.nc"                           ! Newly added on May 8, 2017
                  CALL read_field_R3toR3_ncdf(Uvec, filename, "Ux", "Uy", "Uz", loadSuccessful)

                  Ux = Uvec(:,:,:,1)
                  !CALL dealiasing(Ux) ! commented out by fb
                  Uvec(:,:,:,1) = Ux

                  Uy = Uvec(:,:,:,2)
                  !CALL dealiasing(Uy) ! commented out by fb
                  Uvec(:,:,:,2) = Uy

                  Uz = Uvec(:,:,:,3)
                  !CALL dealiasing(Uz) ! commented out by fb
                  Uvec(:,:,:,3) = Uz

                  DO nn = 2,10                                     ! Newly added on March 23, 2017, add perturbation even for reading data
                     DO kk=1, local_N
                        DO jj=1, n(2)
                           DO ii=1, n(1)
                              X = REAL(ii-1,pr)*dx(1)
                              Y = REAL(jj-1,pr)*dx(2)
                              Z = REAL(kk+local_k_offset-1,pr)*dx(3)
                        !                              Z = REAL(kk-1,pr)*dx(3)
                              Ux(ii,jj,kk) = Ux(ii,jj,kk) + 0.01_pr*(1.0_pr/REAL(nn,pr)**2)*SIN(2.0_pr*PI*REAL(nn,pr)*X)*COS(2.0_pr*PI*REAL(nn,pr)*Y)*SIN(2.0_pr*PI*REAL(nn,pr)*Z)
                              Uy(ii,jj,kk) = Uy(ii,jj,kk) + 0.01_pr*(1.0_pr/REAL(nn,pr)**2)*COS(2.0_pr*PI*REAL(nn,pr)*X)*COS(2.0_pr*PI*REAL(nn,pr)*Y)*COS(2.0_pr*PI*REAL(nn,pr)*Z)
                              Uz(ii,jj,kk) = Uz(ii,jj,kk) + 0.01_pr*(1.0_pr/REAL(nn,pr)**2)*COS(2.0_pr*PI*REAL(nn,pr)*X)*SIN(2.0_pr*PI*REAL(nn,pr)*Y)*SIN(2.0_pr*PI*REAL(nn,pr)*Z)
                           END DO
                        END DO
                     END DO
                  END DO

                  CALL divAvg_free(Uvec)   
               end if

           END SELECT

                        ! Deallocate variables
           DEALLOCATE( Ux )
           DEALLOCATE( Uy )
           DEALLOCATE( Uz )

      END SUBROUTINE initial_guess

      subroutine save_fourier_transform_of_uvec()
         USE global_variables
         USE fftwfunction
         use data_ops
         use mpi
         IMPLICIT NONE

         character(len=200) :: filename_real, filename_imag

         complex(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3) :: aux, faux
         real(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3) :: uvec_fourier_real, uvec_fourier_imag
         real(pr), dimension(1:n(1),1:n(2),1:local_N) :: fx, fy, fz



         filename_real = inputDir//"uvec_fourier_realPart_new.nc"
         filename_imag = inputDir//"uvec_fourier_imagPart_new.nc"

         if(rank==0) print*, "starting save fourier"

         aux = dcmplx(uvec, 0.0_pr)		
         CALL fftfwdv(aux, faux)

         uvec_fourier_real = real(faux)
         uvec_fourier_imag = aimag(faux)

         if(rank==0) print*, "rank0", uvec_fourier_real(1,1,1,1), uvec_fourier_real(2,2,2,2), uvec_fourier_imag(1,1,1,1), uvec_fourier_imag(2,2,2,2)
         if(rank==1) print*, "rank1", uvec_fourier_real(1,1,1,1), uvec_fourier_real(2,2,2,2), uvec_fourier_imag(1,1,1,1), uvec_fourier_imag(2,2,2,2)

         fx = uvec_fourier_real(:,:,:,1)
         fy = uvec_fourier_real(:,:,:,2)
         fz = uvec_fourier_real(:,:,:,3)
         CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
         CALL save_field_R3toR3_ncdf(fx,fy,fz,"FU_r_x", "FU_r_y", "FU_r_z", trim(adjustl(filename_real)), "netCDF")
         CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
         
         fx = uvec_fourier_imag(:,:,:,1)
         fy = uvec_fourier_imag(:,:,:,2)
         fz = uvec_fourier_imag(:,:,:,3)
         CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
         CALL save_field_R3toR3_ncdf(fx,fy,fz,"FU_i_x", "FU_i_y", "FU_i_z", trim(adjustl(filename_imag)), "netCDF")
         CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

         if(rank==0) print*, "exiting program"
         if(rank==0) print*, "to refine restart with according function call"
         call exit
         if(rank==0) print*, "ending save fourier"
         
      end subroutine save_fourier_transform_of_uvec


      subroutine load_fourier_transform_of_uvec_refining(resol_pre)
         USE global_variables
         use netCDF
         USE fftwfunction
         use data_ops
         use mpi
         IMPLICIT NONE


         integer, intent(in) :: resol_pre

         character(len=200) :: filename_real, filename_imag

         complex(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3) :: aux, faux
         real(pr), dimension(:,:,:,:), allocatable :: uvec_fourier_real, uvec_fourier_imag
         real(pr), dimension(:,:,:), allocatable :: local_f
         real(pr), dimension(1:n(1),1:n(2),1:local_N) :: fx, fy, fz

         real(pr) :: adjustment_factor

         integer, dimension(1:3) :: n_pre
         integer(c_intptr_t), dimension(1:3) :: C_n_pre

         integer :: ii, jj, kk, ll, mm, nn
         INTEGER :: ncout, ncid, fid, dimids(3)
         INTEGER :: x_dimid, y_dimid, z_dimid
         INTEGER :: fname_len, nx_ncdf, ny_ncdf, nz_ncdf
         INTEGER, DIMENSION(1:3) :: starts, counts

         integer(c_intptr_t) :: C_local_alloc_pre, C_local_N_pre, C_local_k_offset_pre       ! Newly Added Jan 17; Modified on March 18, 2017
         integer :: local_alloc_pre, local_N_pre, local_k_offset_pre                         ! Newly Added March 18, 2017

         filename_real = inputDir//"uvec_fourier_realPart_n64.nc"
         filename_imag = inputDir//"uvec_fourier_imagPart_n64.nc"

         if(rank==0) print*, "starting load fourier"

         ! !if(resol_pre <= 0 .or. resol_pre==resol) then
         ! !   allocate( uvec_fourier_real(1:n(1),1:n(2),1:local_N,1:3) )
         ! !   allocate( uvec_fourier_imag(1:n(1),1:n(2),1:local_N,1:3) )
         ! !   CALL read_field_R3toR3_ncdf(uvec_fourier_real, filename_real, "FU_r_x", "FU_r_y", "FU_r_z")
         ! !   CALL read_field_R3toR3_ncdf(uvec_fourier_imag, filename_imag, "FU_i_x", "FU_i_y", "FU_i_z")
         ! !else
         ! if(.true.) then
         !    ! refine from resol_pre
         !    n_pre(1) = resol_pre
         !    n_pre(2) = resol_pre
         !    n_pre(3) = resol_pre
         !    C_n_pre(1) = resol_pre
         !    C_n_pre(2) = resol_pre
         !    C_n_pre(3) = resol_pre

            
         !    C_local_alloc_pre = fftw_mpi_local_size_3d(C_n_pre(3), C_n_pre(2), C_n_pre(1), MPI_COMM_WORLD, C_local_N_pre, C_local_k_offset_pre)   ! Newly Added March 18, 2017, use "C_..."
         !    local_N_pre = int( C_local_N_pre )
         !    local_k_offset_pre = int( C_local_k_offset_pre )
            
         !    print*, "rank", rank, "n_pre(:)", n_pre(:), "local_n_pre", local_n_pre

         !    starts = (/ 1, 1, rank*local_N_pre+1 /)
         !    counts = (/ n_pre(1), n_pre(2), local_N_pre /)

            
         !    allocate( uvec_fourier_real(1:n(1),1:n(2),1:local_N,1:3) )
         !    allocate( uvec_fourier_imag(1:n(1),1:n(2),1:local_N,1:3) )
         !    allocate( local_f(1:n_pre(1),1:n_pre(2),1:local_N_pre) )

         !    uvec_fourier_real = 0.0_pr
         !    uvec_fourier_imag = 0.0_pr

         !    !--------------------------
         !    ! START netCDF ROUTINES
         !    !--------------------------
         !    DO ii=0,np-1
         !       IF (rank == ii) THEN 
         !          ncout = nf90_open(filename_real, NF90_NOWRITE, ncid)
         !          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

         !          ncout = nf90_inq_dimid(ncid, "x", x_dimid)
         !          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
         !          ncout = nf90_inq_dimid(ncid, "y", y_dimid)
         !          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
         !          ncout = nf90_inq_dimid(ncid, "z", z_dimid)
         !          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

         !          ncout = nf90_inquire_dimension(ncid, x_dimid, len = nx_ncdf)
         !          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
         !          ncout = nf90_inquire_dimension(ncid, y_dimid, len = ny_ncdf)
         !          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
         !          ncout = nf90_inquire_dimension(ncid, z_dimid, len = nz_ncdf)
         !          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

         !          ncout = nf90_inq_varid(ncid, "FU_r_x", fid)
         !          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
         !          ncout = nf90_get_var(ncid, fid, local_f, start = starts, count = counts)
         !          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

         !          do mm=0,local_N_pre/2
         !             do kk=0,n_pre(2)/2

         !                do kk = 0, n_pre(2)-1
         !                   if (kk<=n_pre(2)/2) then
         !                   else
         !                   end if
         !                end do 
         !                do jj = 0, n_pre(1)-1
         !                   if (jj<=n_pre(1)/2) then
         !                   else
         !                   end if
         !                end do 
         !                DO i = 0, n(1)-1
         !                   IF (i<=n(1)/2) THEN
         !                      K1(i+1) = 2.0_pr*PI*REAL(i,pr)
         !                   ELSE
         !                      K1(i+1) = 2.0_pr*PI*REAL(i-n(1),pr)
         !                   END IF
         !                end do
         !                do jj=0,n_pre(1)/2
         !                   uvec_fourier_real(jj,kk,mm,1) = local_f(jj,kk,mm,1)
         !                   uvec_fourier_real(n(1)-jj,kk,mm,1) = local_f(n_pre(1)-jj,kk,mm,1)
         !                   uvec_fourier_real(jj,n(2)-kk,mm,1) = local_f(jj,n_pre(2)-kk,mm,1)
         !                   uvec_fourier_real(n(1)-jj,n(2)-kk,mm,1) = local_f(n_pre(1)-jj,n_pre(2)-kk,mm,1)
         !                end do
         !             end do
         !          end do



         !          ncout = nf90_inq_varid(ncid, "FU_r_y", fid)
         !          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
         !          ncout = nf90_get_var(ncid, fid, local_f, start = starts, count = counts)
         !          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
         !          uvec_fourier_real(1:n_pre(1),1:n_pre(2),1:local_N_pre,2) = local_f

         !          ncout = nf90_inq_varid(ncid, "FU_r_z", fid)
         !          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
         !          ncout = nf90_get_var(ncid, fid, local_f, start = starts, count = counts)
         !          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
         !          uvec_fourier_real(1:n_pre(1),1:n_pre(2),1:local_N_pre,3) = local_f

         !          ncout = nf90_close(ncid)

         !       END IF
         !       CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
         !    END DO     
         !    CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)



         !    !--------------------------
         !    ! START netCDF ROUTINES
         !    !--------------------------
         !    DO ii=0,np-1
         !       IF (rank == ii) THEN 
         !          ncout = nf90_open(filename_imag, NF90_NOWRITE, ncid)
         !          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

         !          ncout = nf90_inq_dimid(ncid, "x", x_dimid)
         !          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
         !          ncout = nf90_inq_dimid(ncid, "y", y_dimid)
         !          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
         !          ncout = nf90_inq_dimid(ncid, "z", z_dimid)
         !          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

         !          ncout = nf90_inquire_dimension(ncid, x_dimid, len = nx_ncdf)
         !          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
         !          ncout = nf90_inquire_dimension(ncid, y_dimid, len = ny_ncdf)
         !          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
         !          ncout = nf90_inquire_dimension(ncid, z_dimid, len = nz_ncdf)
         !          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

         !          ncout = nf90_inq_varid(ncid, "FU_i_x", fid)
         !          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
         !          ncout = nf90_get_var(ncid, fid, local_f, start = starts, count = counts)
         !          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
         !          uvec_fourier_imag(1:n_pre(1),1:n_pre(2),1:local_N_pre,1) = local_f

         !          ncout = nf90_inq_varid(ncid, "FU_i_y", fid)
         !          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
         !          ncout = nf90_get_var(ncid, fid, local_f, start = starts, count = counts)
         !          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
         !          uvec_fourier_imag(1:n_pre(1),1:n_pre(2),1:local_N_pre,2) = local_f

         !          ncout = nf90_inq_varid(ncid, "FU_i_z", fid)
         !          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
         !          ncout = nf90_get_var(ncid, fid, local_f, start = starts, count = counts)
         !          IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
         !          uvec_fourier_imag(1:n_pre(1),1:n_pre(2),1:local_N_pre,3) = local_f

         !          ncout = nf90_close(ncid)

         !       END IF
         !       CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
         !    END DO     
         !    CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

         !    DEALLOCATE(local_f)

         ! end if



         ! faux = dcmplx(uvec_fourier_real, uvec_fourier_imag)
         ! call fftbwdv(faux, aux)
         ! uvec = real(aux)

         ! adjustment_factor = real(resol,pr)/real(resol_pre,pr)
         ! uvec(:,:,:,:) = adjustment_factor*uvec(:,:,:,:)

         if(rank==0) print*, "ending load fourier"
         
      end subroutine load_fourier_transform_of_uvec_refining




      !===================================
      !  PERTURBATION FOR KAPPA TEST
      !===================================
      SUBROUTINE kappa_test_pert(phi_pert, mytype, m1, m2, m3, successful)
         USE global_variables
         USE fftwfunction
         use data_ops
         use mpi
         IMPLICIT NONE

         REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(OUT) :: phi_pert
         complex(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3) :: aux, faux
         complex(pr), DIMENSION(1:n(1),1:n(2),1:local_N) :: auxScalar, fauxScalar
         CHARACTER(len=*), INTENT(IN) :: mytype
         CHARACTER(100) :: filename
         CHARACTER(4) :: auxStr
         REAL(pr), INTENT(IN) :: m1, m2, m3
         REAL(pr), DIMENSION(1:3) :: dx
         INTEGER :: ii, jj, kk, ll
         real(pr), dimension(1:3) :: l2NormTestLoc, l2NormTestGlob
         REAL(pr) :: X, Y, Z, ampl, norm_K, normalization_const
         logical :: loadSuccessful
         logical, intent(out) :: successful

         successful = .true.



         
         if (n(1)<100) then
            WRITE(auxStr, '(i2)') n(1)
            auxStr = '00'//auxStr(1:2)
         elseif (n(1)<1000) then
            WRITE(auxStr, '(i3)') n(1)
            auxStr = '0'//auxStr(1:3)
         elseif (n(1)<10000) then
            WRITE(auxStr, '(i4)') n(1)
         end if

         dx = 1.0/REAL(n, pr)
         SELECT CASE (mytype)

            case ("load-te0080")
               filename = inputDir//"FRT_N256E500T017_Uvec_fwdTE0080.nc"
               CALL read_field_R3toR3_ncdf(phi_pert, filename, "Ux", "Uy", "Uz", loadSuccessful)
               successful = loadSuccessful


            CASE ("sine")
               DO kk=1, local_N
                  DO jj=1, n(2)
                     DO ii=1, n(1)
                        X = REAL(ii-1,pr)*dx(1)
                        Y = REAL(jj-1,pr)*dx(2)
                        Z = REAL(kk+local_k_offset-1,pr)*dx(3)
                        phi_pert(ii,jj,kk,1) = SIN(2.0_pr*PI*m1*X)*COS(2.0_pr*PI*m2*Y)*COS(2.0_pr*PI*m3*Z)
                        phi_pert(ii,jj,kk,2) = -COS(2.0_pr*PI*m1*X)*SIN(2.0_pr*PI*m2*Y)*COS(2.0_pr*PI*m3*Z)
                        phi_pert(ii,jj,kk,3) = COS(2.0_pr*PI*m1*X)*COS(2.0_pr*PI*m2*Y)*SIN(2.0_pr*PI*m3*Z)
                     END DO
                  END DO
               END DO               
               call divAvg_free(phi_pert)
            
            case ("divfree-sine")
               DO kk=1, local_N
                  DO jj=1, n(2)
                     DO ii=1, n(1)
                        X = REAL(ii-1,pr)*dx(1)
                        Y = REAL(jj-1,pr)*dx(2)
                        Z = REAL(kk+local_k_offset-1,pr)*dx(3)
                        !m2 = m1 + m3      ! such that phi_pert is divfree
                        phi_pert(ii,jj,kk,1) = SIN(2.0_pr*PI*m1*X)*COS(2.0_pr*PI*(m1 + m3)*Y)*COS(2.0_pr*PI*m3*Z)
                        phi_pert(ii,jj,kk,2) = -COS(2.0_pr*PI*m1*X)*SIN(2.0_pr*PI*(m1 + m3)*Y)*COS(2.0_pr*PI*m3*Z)
                        phi_pert(ii,jj,kk,3) = COS(2.0_pr*PI*m1*X)*COS(2.0_pr*PI*(m1 + m3)*Y)*SIN(2.0_pr*PI*m3*Z)
                     END DO
                  END DO
               END DO


            case ("load-random-a")
               filename = inputDir//"n"//auxStr//"_random_a.nc"
               CALL read_field_R3toR3_ncdf(phi_pert, filename, "Ux", "Uy", "Uz", loadSuccessful)
               successful = loadSuccessful
               call divAvg_free(phi_pert)

            case ("load-random-b")
               filename = inputDir//"n"//auxStr//"_random_b.nc"
               CALL read_field_R3toR3_ncdf(phi_pert, filename, "Ux", "Uy", "Uz", loadSuccessful)
               successful = loadSuccessful
               call divAvg_free(phi_pert)

            CASE ("load-random-smooth-a")
               filename = inputDir//"n"//auxStr//"_random_a.nc"
               CALL read_field_R3toR3_ncdf(phi_pert, filename, "Ux", "Uy", "Uz", loadSuccessful)
               successful = loadSuccessful
               if(.not.loadSuccessful) phi_pert = 0.0_pr
               do kk=1,3
                  call dealias_scalar(phi_pert(:,:,:,kk), 20.0_pr)
               end do
               call divAvg_free(phi_pert)
               
            CASE ("load-k-random-a")
               filename = inputDir//"n"//auxStr//"_random_a.nc"
               CALL read_field_R3toR3_ncdf(phi_pert, filename, "Ux", "Uy", "Uz", loadSuccessful)
               successful = loadSuccessful
               if(.not.loadSuccessful) phi_pert = 0.0_pr
               aux = dcmplx(phi_pert, 0.0_pr)		
               CALL fftfwdv(aux, faux)
               DO kk=1, local_N
                  DO jj=1, n(2)
                     DO ii=1, n(1)
                        norm_K = SQRT( K1(ii)**2 + K2(jj)**2 + K3(kk+local_k_offset)**2 )
                        If (norm_k < MACH_EPSILON) then
                           do ll=1,3
                              faux(ii,jj,kk,ll) = 0.0_pr
                           end do
                        else
                           do ll=1,3
                              faux(ii,jj,kk,ll) = faux(ii,jj,kk,ll)*m1
                           end do
                        END IF
                     END DO
                  END DO
               END DO
               CALL fftbwdv(faux, aux)
               phi_pert = real(aux,pr)
               call divAvg_free(phi_pert)
               
               
            CASE ("load-k-random-b")
               filename = inputDir//"n"//auxStr//"_random_b.nc"
               CALL read_field_R3toR3_ncdf(phi_pert, filename, "Ux", "Uy", "Uz", loadSuccessful)
               successful = loadSuccessful
               if(.not.loadSuccessful) phi_pert = 0.0_pr
               aux = dcmplx(phi_pert, 0.0_pr)		
               CALL fftfwdv(aux, faux)
               DO kk=1, local_N
                  DO jj=1, n(2)
                     DO ii=1, n(1)
                        norm_K = SQRT( K1(ii)**2 + K2(jj)**2 + K3(kk+local_k_offset)**2 )
                        If (norm_k < MACH_EPSILON) then
                           do ll=1,3
                              faux(ii,jj,kk,ll) = 0.0_pr
                           end do
                        else
                           do ll=1,3                        	
                              faux(ii,jj,kk,ll) = faux(ii,jj,kk,ll)*m1
                           end do
                        END IF
                     END DO
                  END DO
               END DO
               CALL fftbwdv(faux, aux)
               phi_pert = real(aux,pr)
               call divAvg_free(phi_pert)
               
            CASE ("load-random-poly-a")
               filename = inputDir//"n"//auxStr//"_random_a.nc"
               CALL read_field_R3toR3_ncdf(phi_pert, filename, "Ux", "Uy", "Uz", loadSuccessful)
               successful = loadSuccessful
               if(.not.loadSuccessful) phi_pert = 0.0_pr
               aux = dcmplx(phi_pert, 0.0_pr)		
               CALL fftfwdv(aux, faux)
               DO kk=1, local_N
                  DO jj=1, n(2)
                     DO ii=1, n(1)
                        norm_K = SQRT( K1(ii)**2 + K2(jj)**2 + K3(kk+local_k_offset)**2 )
                        If (norm_k < MACH_EPSILON) then
                           do ll=1,3
                              faux(ii,jj,kk,ll) = 0.0_pr
                           end do
                        else
                           do ll=1,3                        	
                              faux(ii,jj,kk,ll) = faux(ii,jj,kk,ll)*(10.0_pr**m2)*(norm_K/(real(n(1),pr)/4.0_pr))**(m1)
                           end do
                        END IF
                     END DO
                  END DO
               END DO
               CALL fftbwdv(faux, aux)
               phi_pert = real(aux,pr)
               call divAvg_free(phi_pert)
               
            CASE ("create-random-poly-b")
               call random_number(phi_pert)
               phi_pert = 2.0_pr*phi_pert-1.0_pr
               aux = dcmplx(phi_pert, 0.0_pr)

               CALL fftfwdv(aux, faux)
               DO kk=1, local_N
                  DO jj=1, n(2)
                     DO ii=1, n(1)
                        norm_K = SQRT( K1(ii)**2 + K2(jj)**2 + K3(kk+local_k_offset)**2 )
                        If (norm_k < MACH_EPSILON) then
                           do ll=1,3
                              faux(ii,jj,kk,ll) = 0.0_pr
                           end do
                        END IF
                     END DO
                  END DO
               END DO

               CALL fftbwdv(faux, aux)

               phi_pert = real(aux,pr)
               call divAvg_free(phi_pert)

               ! normalization
               normalization_const = 100.0_pr*global_summed_field_inner_product(phi_pert,phi_pert,"L2")

               phi_pert(:,:,:,:) = phi_pert(:,:,:,:)/normalization_const

               aux = dcmplx(phi_pert, 0.0_pr)		
               CALL fftfwdv(aux, faux)
               DO kk=1, local_N
                  DO jj=1, n(2)
                     DO ii=1, n(1)
                        norm_K = SQRT( K1(ii)**2 + K2(jj)**2 + K3(kk+local_k_offset)**2 )
                        If (norm_k < MACH_EPSILON) then
                           do ll=1,3
                              faux(ii,jj,kk,ll) = 0.0_pr
                           end do
                        else
                           do ll=1,3                        	
                              faux(ii,jj,kk,ll) = faux(ii,jj,kk,ll)*(10.0_pr**m2)*(norm_K/(real(n(1),pr)/4.0_pr))**(m1)
                           end do
                        END IF
                     END DO
                  END DO
               END DO
               CALL fftbwdv(faux, aux)
               phi_pert = real(aux,pr)
               call divAvg_free(phi_pert)

            CASE ("load-random-poly-b")
               filename = inputDir//"n"//auxStr//"_random_b.nc"
               CALL read_field_R3toR3_ncdf(phi_pert, filename, "Ux", "Uy", "Uz", loadSuccessful)
               successful = loadSuccessful
               if(.not.loadSuccessful) phi_pert = 0.0_pr
               aux = dcmplx(phi_pert, 0.0_pr)		
               CALL fftfwdv(aux, faux)
               DO kk=1, local_N
                  DO jj=1, n(2)
                     DO ii=1, n(1)
                        norm_K = SQRT( K1(ii)**2 + K2(jj)**2 + K3(kk+local_k_offset)**2 )
                        If (norm_k < MACH_EPSILON) then
                           do ll=1,3
                              faux(ii,jj,kk,ll) = 0.0_pr
                           end do
                        else
                           do ll=1,3                        	
                              faux(ii,jj,kk,ll) = faux(ii,jj,kk,ll)*(10.0_pr**m2)*(norm_K/(real(n(1),pr)/4.0_pr))**(m1)
                           end do
                        END IF
                     END DO
                  END DO
               END DO
               CALL fftbwdv(faux, aux)
               phi_pert = real(aux,pr)
               call divAvg_free(phi_pert)
               
            CASE ("load-random-exp-a")
               filename = inputDir//"n"//auxStr//"_random_a.nc"
               CALL read_field_R3toR3_ncdf(phi_pert, filename, "Ux", "Uy", "Uz", loadSuccessful)
               successful = loadSuccessful
               if(.not.loadSuccessful) phi_pert = 0.0_pr
               aux = dcmplx(phi_pert, 0.0_pr)		
               CALL fftfwdv(aux, faux)
               DO kk=1, local_N
                  DO jj=1, n(2)
                     DO ii=1, n(1)
                        norm_K = SQRT( K1(ii)**2 + K2(jj)**2 + K3(kk+local_k_offset)**2 )
                        If (norm_k < MACH_EPSILON) then
                           do ll=1,3
                              faux(ii,jj,kk,ll) = 0.0_pr
                           end do
                        else
                           do ll=1,3                        	
                              faux(ii,jj,kk,ll) = faux(ii,jj,kk,ll)*10.0_pr**(-norm_K/5.0_pr)
                           end do
                        END IF
                     END DO
                  END DO
               END DO
               CALL fftbwdv(faux, aux)
               phi_pert = real(aux,pr)
               call divAvg_free(phi_pert)
               
            CASE ("load-random-exp-b")
               filename = inputDir//"n"//auxStr//"_random_b.nc"
               CALL read_field_R3toR3_ncdf(phi_pert, filename, "Ux", "Uy", "Uz", loadSuccessful)
               successful = loadSuccessful
               if(.not.loadSuccessful) phi_pert = 0.0_pr
               aux = dcmplx(phi_pert, 0.0_pr)		
               CALL fftfwdv(aux, faux)
               DO kk=1, local_N
                  DO jj=1, n(2)
                     DO ii=1, n(1)
                        norm_K = SQRT( K1(ii)**2 + K2(jj)**2 + K3(kk+local_k_offset)**2 )
                        If (norm_k < MACH_EPSILON) then
                           do ll=1,3
                              faux(ii,jj,kk,ll) = 0.0_pr
                           end do
                        else
                           do ll=1,3                        	
                              faux(ii,jj,kk,ll) = faux(ii,jj,kk,ll)*10.0_pr**(-norm_K/5.0_pr)
                           end do
                        END IF
                     END DO
                  END DO
               END DO
               CALL fftbwdv(faux, aux)
               phi_pert = real(aux,pr)
               call divAvg_free(phi_pert)

            CASE ("load-random-smooth-b")
               filename = inputDir//"n"//auxStr//"_random_b.nc"
               CALL read_field_R3toR3_ncdf(phi_pert, filename, "Ux", "Uy", "Uz", loadSuccessful)
               successful = loadSuccessful
               if(.not.loadSuccessful) phi_pert = 0.0_pr
               do kk=1,3
                  call dealias_scalar(phi_pert(:,:,:,kk), 20.0_pr)
               end do

            CASE ("save-random")
               call random_number(phi_pert)
               phi_pert = 2.0_pr*phi_pert-1.0_pr
               aux = dcmplx(phi_pert, 0.0_pr)

               CALL fftfwdv(aux, faux)
               DO kk=1, local_N
                  DO jj=1, n(2)
                     DO ii=1, n(1)
                        norm_K = SQRT( K1(ii)**2 + K2(jj)**2 + K3(kk+local_k_offset)**2 )
                        If (norm_k < MACH_EPSILON) then
                           do ll=1,3
                              faux(ii,jj,kk,ll) = 0.0_pr
                           end do
                        END IF
                     END DO
                  END DO
               END DO

               CALL fftbwdv(faux, aux)

               phi_pert = real(aux,pr)
               call divAvg_free(phi_pert)

               ! normalization
               normalization_const = 100.0_pr*global_summed_field_inner_product(phi_pert,phi_pert,"L2")

               phi_pert(:,:,:,:) = phi_pert(:,:,:,:)/normalization_const

               filename = HomeDir//"n"//auxStr//"_random.nc"
               CALL save_field_R3toR3_ncdf(phi_pert(:,:,:,1), phi_pert(:,:,:,2), phi_pert(:,:,:,3), "Ux", "Uy", "Uz", filename, "netCDF")


               CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

               if(rank==0) print*, "finished saving random field, exiting"
               call exit


         END SELECT
   
      END SUBROUTINE kappa_test_pert
      
      !=========================================================
      ! Calculate the L^2 derivative of ||u||_q = B
      ! q ||u||_q^{q-1} nabla ||u||_q = nabla ||u||_q^q = q |u|^{q-2} u
      ! nabla ||u||_q = ||u||_q^{1-q} |u|^{q-2} u
      !=========================================================
      function calcConstraintDerivativeL2(u) result(resultVec)
         use global_variables
         implicit none
         real(pr), dimension(1:n(1),1:n(2),1:local_N,1:3), intent(in) :: u
         real(pr), dimension(1:n(1),1:n(2),1:local_N,1:3) :: resultVec
         real(pr), dimension(1:n(1),1:n(2),1:local_N) :: aux_u_q2
         real(pr) :: norm
         integer :: ii

         norm = calc_global_Lq_norm(u)

         call calc_uk(u,lebesgueQ-2.0_pr,aux_u_q2)                           ! aux_u_q2 = |u|^{q-2}

         resultVec = 0.0_pr
         do ii = 1,3
            resultVec(:,:,:,ii) = aux_u_q2(:,:,:)*u(:,:,:,ii)         ! aux_uq2u = |u|^{q-2}u
            if (toDealias .and. (lebesgueQ-2.0_pr > mach_epsilon)) call dealias_scalar(resultVec(:,:,:,ii), 2.0_pr)
         end do

         resultVec(:,:,:,:) = (norm**(1.0_pr - lebesgueQ)) * resultVec(:,:,:,:)

         
      end function calcConstraintDerivativeL2

      !=========================================================
      ! Calculate and save spectrum of a vectorfield
      !=========================================================
      subroutine calculateSaveSpectrum(vec, fileName)
         use global_variables
         use fftwfunction
         use data_ops
         use mpi
         implicit none
         real(pr), dimension(1:n(1),1:n(2),1:local_n,1:3), intent(in) :: vec
         CHARACTER(len=*), INTENT(IN) :: fileName
         real(pr), dimension(1:n(1),1:2) :: spectrum

         CALL calculate_spectral_data(vec, Spectrum)
         CALL save_spectral_data(Spectrum, fileName)
         
      end subroutine calculateSaveSpectrum

      !=========================================================
      ! Calculate and save spectrum of a matrixfield
      !=========================================================
      subroutine calculateSaveSpectrumMatrix(matrix, fileName)
         use global_variables
         use fftwfunction
         use data_ops
         use mpi
         implicit none
         real(pr), dimension(1:n(1),1:n(2),1:local_n,1:3,1:3), intent(in) :: matrix
         CHARACTER(len=*), INTENT(IN) :: fileName
         real(pr), dimension(1:n(1),1:2,1:3) :: spectrumPart
         real(pr), dimension(1:n(1),1:2) :: spectrum
         integer :: ii
         
         spectrum = 0.0_pr
         do ii = 1,3
            CALL calculate_spectral_data(matrix(:,:,:,:,ii), spectrumPart(:,:,ii))
            spectrum(:,:) = spectrum(:,:) + spectrumPart(:,:,ii)
         end do
         CALL save_spectral_data(Spectrum, fileName)
      end subroutine calculateSaveSpectrumMatrix


      !=========================================================
      ! Calculate L2-gradient of Lp right hand side in physical space
      !=========================================================
      ! nabla^L2 R = (q-2)|u|^{q-4} (u cdot (nu Delta u - nabla p)) u
      !              - |u|^{q-2} nabla p
      !              - 2 nabla (u cdot nabla) Delta^{-1} nabla cdot (|u|^{q-2}u)
      !              + nu |u|^{q-2} Delta u + nu Delta (|u|^{q-2}u)
      ! where
      !     Delta^{-1} f is the solution to Delta v = f, int v = 0
      !
      !     nabla (u cdot nabla) (...) = partial_i u_j partial_j (...)
      !
      !     p = Delta^{-1} (nabla u column nabla u^T)
      !       = Delta^{-1} (partial_i u_j partial_j u_i)
      !=========================================================

      function GradL2ForLq(U) result (V)
         use global_variables
         use fftwfunction
         implicit none
         integer :: ii
         real(pr), dimension(1:n(1),1:n(2),1:local_N,1:3), intent(in) :: U
         real(pr), dimension(1:n(1),1:n(2),1:local_N,1:3) :: V
         !real(pr), dimension(:,:,:), allocatable :: aux, aux2, aux3, aux4, aux_p, aux_u_q2, aux_u_q4
         !real(pr), dimension(:,:,:,:), allocatable :: aux_gradP, aux_DeltaU, aux_uq2u, aux3_vec, aux4_vec, aux5_vec, aux6_vec, aux7_vec, gradU1, gradU2, gradU3, e_u
         
         real(pr), dimension(1:n(1),1:n(2),1:local_N) :: aux, aux2, aux3, aux4, aux_p, aux_u_q2, aux_u_q4
         real(pr), dimension(1:n(1),1:n(2),1:local_N,1:3) :: aux_gradP, aux_DeltaU, aux_uq2u, aux3_vec, aux4_vec, aux5_vec, aux6_vec, aux7_vec, gradU1, gradU2, gradU3, e_u
         

         !allocate( aux(1:n(1),1:n(2),1:local_N) )
         !allocate( aux2(1:n(1),1:n(2),1:local_N) )
         !allocate( aux3(1:n(1),1:n(2),1:local_N) )
         !allocate( aux4(1:n(1),1:n(2),1:local_N) )
         !allocate( aux_p(1:n(1),1:n(2),1:local_N) )
         !allocate( aux_u_q2(1:n(1),1:n(2),1:local_N) )
         !allocate( aux_u_q4(1:n(1),1:n(2),1:local_N) )

         !allocate( aux_gradP(1:n(1),1:n(2),1:local_N,1:3) )
         !allocate( aux_DeltaU(1:n(1),1:n(2),1:local_N,1:3) )
         !allocate( aux_uq2u(1:n(1),1:n(2),1:local_N,1:3) )
         !allocate( aux3_vec(1:n(1),1:n(2),1:local_N,1:3) )
         !allocate( aux4_vec(1:n(1),1:n(2),1:local_N,1:3) )
         !allocate( aux5_vec(1:n(1),1:n(2),1:local_N,1:3) )
         !allocate( aux6_vec(1:n(1),1:n(2),1:local_N,1:3) )
         !allocate( aux7_vec(1:n(1),1:n(2),1:local_N,1:3) )
         !allocate( gradU1(1:n(1),1:n(2),1:local_N,1:3) )
         !allocate( gradU2(1:n(1),1:n(2),1:local_N,1:3) )
         !allocate( gradU3(1:n(1),1:n(2),1:local_N,1:3) )
         !allocate( e_u(1:n(1),1:n(2),1:local_N,1:3) )


         call calc_uk(u,lebesgueQ-2.0_pr,aux_u_q2)                           ! aux_u_q2 = |u|^{q-2}
         

         call calc_nablaUnablaUt(u, aux)                             ! aux = nabla u : nabla u^T

         call solve_poisson(aux, 1.0_pr, aux_p)                      ! aux_p = p = Delta^{-1} (nabla u : nabla u^T) = Delta^{-1} (aux)
         
         call gradient(aux_p, aux_gradP)                             ! aux_gradP = nabla p

         aux_DeltaU = u
         call laplacian(aux_DeltaU)                                  ! aux_DeltaU = Delta u
         

         aux2 = 0.0_pr
         if(use_e_u_instead_of_uqMinus4) then
            e_u = calc_unitVectorInUdirection(u)
            do ii = 1,3
               aux2(:,:,:) = aux2(:,:,:) + e_u(:,:,:,ii)*(viscCoefficient*visc*aux_DeltaU(:,:,:,ii)-pressureCoefficient*aux_gradP(:,:,:,ii)) ! aux2 = e_u cdot (nu Delta u - nabla p)
            end do
            if (toDealias .and. dealiase_if_mult_by_e_u) call dealias_scalar(aux2, 2.0_pr)

            aux2(:,:,:) = aux_u_q2(:,:,:) * aux2(:,:,:)                 ! aux2 = |u|^{q-2} (e_u cdot (nu Delta u - nabla p))
            if (toDealias) call dealias_scalar(aux2, 2.0_pr)

            do ii = 1,3
               aux3_vec(:,:,:,ii) = aux2(:,:,:) * e_u(:,:,:,ii)         ! aux3_vec = |u|^{q-2} (e_u cdot (nu Delta u - nabla p)) e_u
               if (toDealias .and. dealiase_if_mult_by_e_u) call dealias_scalar(aux3_vec(:,:,:,ii), 2.0_pr)
            end do
         else
            do ii = 1,3
               aux2(:,:,:) = aux2(:,:,:) + u(:,:,:,ii)*(viscCoefficient*visc*aux_DeltaU(:,:,:,ii)-pressureCoefficient*aux_gradP(:,:,:,ii)) ! aux2 = u cdot (nu Delta u - nabla p)
            end do
            if (toDealias) call dealias_scalar(aux2, 2.0_pr)

            call calc_uk(u,lebesgueQ-4.0_pr,aux_u_q4)                   ! aux_u_q4 = |u|^{q-4}
            aux2(:,:,:) = aux_u_q4(:,:,:) * aux2(:,:,:)                 ! aux2 = |u|^{q-4} (u cdot (nu Delta u - nabla p))
            if (toDealias) call dealias_scalar(aux2, 2.0_pr)

            do ii = 1,3
               aux3_vec(:,:,:,ii) = aux2(:,:,:) * u(:,:,:,ii)           ! aux3_vec = |u|^{q-4} (u cdot (nu Delta u - nabla p)) u
               if (toDealias) call dealias_scalar(aux3_vec(:,:,:,ii), 2.0_pr)
            end do
         end if
         ! aux3_vec = |u|^{q-4} (u cdot (nu Delta u - nabla p)) u
         ! aux3_vec = |u|^{q-2} (e_u cdot (nu Delta u - nabla p)) e_u


         
         do ii = 1,3
            aux_uq2u(:,:,:,ii) = aux_u_q2(:,:,:)*u(:,:,:,ii)         ! aux_uq2u = |u|^{q-2}u
            if (toDealias) call dealias_scalar(aux_uq2u(:,:,:,ii), 2.0_pr)
         end do

         call divergence(aux_uq2u, aux3)                             ! aux3 = nabla cdot (|u|^{q-2}u)
         call solve_poisson(aux3, 1.0_pr, aux4)                      ! aux4 = Delta^{-1} nabla cdot (|u|^{q-2}u)

         call gradient(aux4, aux4_vec)                               ! aux4_vec = nabla (Delta^{-1} nabla cdot (|u|^{q-2}u))
         
         call gradient(u(:,:,:,1),gradU1)                            ! gradU1 = nabla u_1
         call gradient(u(:,:,:,2),gradU2)                            ! gradU2 = nabla u_2
         call gradient(u(:,:,:,3),gradU3)                            ! gradU3 = nabla u_3


         do ii= 1,3
            !aux5_vec = nabla (u cdot nabla) Delta^{-1} nabla cdot (|u|^{q-2}u)
            aux5_vec(:,:,:,ii) = gradU1(:,:,:,ii)*aux4_vec(:,:,:,1)+gradU2(:,:,:,ii)*aux4_vec(:,:,:,2)+gradU3(:,:,:,ii)*aux4_vec(:,:,:,3)
            if (toDealias) call dealias_scalar(aux5_vec(:,:,:,ii), 2.0_pr)
         end do


         aux7_vec= aux_uq2u                                          ! aux7_vec = |u|^{q-2} u
         call laplacian(aux7_vec)                                    ! aux7_vec = Delta (|u|^{q-2} u)




         do ii = 1,3
            
            aux4_vec(:,:,:,ii) = aux_u_q2(:,:,:)*aux_gradP(:,:,:,ii) ! aux4_vec = |u|^{q-2} nabla p
            if (toDealias) call dealias_scalar(aux4_vec(:,:,:,ii), 2.0_pr)

            aux6_vec(:,:,:,ii) = aux_u_q2(:,:,:)*aux_DeltaU(:,:,:,ii)! aux6_vec = |u|^{q-2} Delta u
            if (toDealias) call dealias_scalar(aux6_vec(:,:,:,ii), 2.0_pr)

            v(:,:,:,ii) = (lebesgueQ-2.0_pr) * aux3_vec(:,:,:,ii) &
               - pressureCoefficient * aux4_vec(:,:,:,ii) &
               - pressureCoefficient * 2.0_pr * aux5_vec(:,:,:,ii) &
               + viscCoefficient * visc * aux6_vec(:,:,:,ii) &
               + viscCoefficient * visc * aux7_vec(:,:,:,ii) 

            ! nabla^L2 R = (q-2)|u|^{q-2} (e_u cdot (nu Delta u - nabla p)) e_u
            !              - |u|^{q-2} nabla p
            !              - 2 nabla (u cdot nabla) Delta^{-1} nabla cdot (|u|^{q-2}u)
            !              + nu |u|^{q-2} Delta u
            !              + nu Delta (|u|^{q-2}u)

            ! nabla^L2 R = (q-2)|u|^{q-4} (u cdot (nu Delta u - nabla p)) u
            !              - |u|^{q-2} nabla p
            !              - 2 nabla (u cdot nabla) Delta^{-1} nabla cdot (|u|^{q-2}u)
            !              + nu |u|^{q-2} Delta u
            !              + nu Delta (|u|^{q-2}u)

         end do



         !deallocate( aux, aux2, aux3, aux4, aux_p, aux_u_q2, aux_u_q4 )
         !deallocate( aux_gradP, aux_DeltaU, aux_uq2u, aux3_vec, aux4_vec, aux5_vec, aux6_vec, aux7_vec, gradU1, gradU2, gradU3, e_u)

      end function GradL2ForLq


      !=========================================================
      ! Calculate d/dt Lq right hand side in physical space
      !     THIS ONLY CALCULATES THE LOCAL PART
      !=========================================================
      ! d/dt Lq = - nu || |u|^((q-2)/2) |nabla u| ||_2^2 - 4(q-2)/q^2 nu ||nabla |u|^(q/2)||_2^2                        + (q-2) int p |u|^(q-4) u cdot (u cdot nabla) u
      !         = - nu int |u|^(q-2) |nabla u|^2         - (q-2) nu int |u|^(q-4) u_i partial_j u_i u_k partial_j u_k   + (q-2) int p |u|^(q-4) u cdot (u cdot nabla) u
      !         =                 R_1                    +                              R_2                             +                       R_3
      !         = R
      ! p = Delta^{-1} (nabla u cdot nabla u^T)
      ! returns (/ R, R1+R2, R3 /)  = (/ localResult, viscosityContribution, nonlinearPressureContribution  /)
      !=========================================================
      function calc_local_dLqdt_inclParts(u, q) result (local_result)
         use global_variables
         implicit none
         real(pr) :: q, R_1, R_2, R_3
         real(pr), dimension(3) :: local_result
         real(pr), dimension(1:n(1),1:n(2),1:local_n,1:3), intent(in) :: u
         real(pr), dimension(1:n(1),1:n(2),1:local_n,1:3) :: e_u
         real(pr), dimension(1:n(1),1:n(2),1:local_n) :: aux_u_q2, aux_p, aux_eegradu, aux_uugradu, aux, aux_u_q4
         

         call calc_uk(u, q-2.0_pr, aux_u_q2)             ! aux_u_q2 = |u|^(q-2)
         
         call calc_nablaModSquared(u, aux)               ! aux = |nabla u|^2
         
         R_1 = - viscCoefficient * visc * inner_product(aux_u_q2, aux, "L2")
         !R_1  = - nu int |u|^(q-2) |nabla u|^2



         e_u = calc_unitVectorInUdirection(u)
         if(use_e_u_instead_of_uqMinus4) then
            ! this is the one that somehow makes using e_u much worse !
            ! all other occurances are not that impactful !
            aux = calc_vgraduvgradu(u, e_u)                  ! aux = (e_u)_i partial_j u_i (e_u)_k partial_j u_k
            R_2 = - viscCoefficient * (q-2.0_pr) * visc * inner_product(aux_u_q2, aux, "L2")
            !R_2 = - (q-2) nu int |u|^(q-2) (e_u)_i partial_j u_i (e_u)_k partial_j u_k
         else
            call calc_uk(u, q-4.0_pr, aux_u_q4)             ! aux_u_q4 = |u|^(q-4)
            aux = calc_vgraduvgradu(u, u)                     ! aux = u_i partial_j u_i u_k partial_j u_k
            R_2 = - viscCoefficient * (q-2.0_pr) * visc * inner_product(aux_u_q4, aux, "L2")
            !R_2 = - (q-2) nu int |u|^(q-4) u_i partial_j u_i u_k partial_j u_k
         end if
         
         

         call calc_nablaUnablaUt(u, aux)                 ! aux = nabla u : nabla u^T
         call solve_poisson(aux, 1.0_pr, aux_p)          ! aux_p = p = Delta^{-1} (nabla u : nabla u^T) = Delta^{-1} (aux2)

         if(use_e_u_instead_of_uqMinus4) then
            aux(:,:,:) = aux_p(:,:,:)*aux_u_q2(:,:,:)       ! aux = p |u|^(q-2)
            if (toDealias) call dealias_scalar(aux, 2.0_pr)
            aux_eegradu = calc_vvgradu(u, e_u)              ! aux_eegradu = e_u (e_u cdot nabla) u
            R_3 = pressureCoefficient * (q-2.0_pr)*inner_product(aux, aux_eegradu, "L2")
            ! R_3  = (q-2) int p |u|^(q-2) e_u cdot (e_u cdot nabla) u
         else
            aux(:,:,:) = aux_p(:,:,:)*aux_u_q4(:,:,:)       ! aux = p |u|^(q-4)
            if (toDealias) call dealias_scalar(aux, 2.0_pr)
            aux_uugradu = calc_vvgradu(u,u)                 ! aux_uugradu = u cdot (u cdot nabla) u
            R_3 = pressureCoefficient * (q-2.0_pr)*inner_product(aux, aux_uugradu, "L2")
            !R_3  = (q-2) int p |u|^(q-4) u cdot (u cdot nabla) u
         end if

         local_result = (/R_1 + R_2 + R_3, R_1+R_2, R_3/)

      end function calc_local_dLqdt_inclParts

      !=========================================================
      ! Calculate d/dt Lq right hand side in physical space
      !                    local
      !=========================================================
      ! returns localResult (scalar)
      !=========================================================
      function calc_local_dLqdt(u, q) result (local_result)
         use global_variables

         real(pr), dimension(3) :: local_result_inclParts
         real(pr) :: local_result
         real(pr), intent(in) :: q
         real(pr), dimension(1:n(1),1:n(2),1:local_n,1:3), intent(in) :: u
         
         local_result_inclParts = calc_local_dLqdt_inclParts(u,q)
         local_result = local_result_inclParts(1)

      end function calc_local_dLqdt


      !=========================================================
      ! Calculate d/dt Lq right hand side in physical space
      !                    global
      !=========================================================
      ! d/dt Lq = - nu || |u|^((q-2)/2) |nabla u| ||_2^2 - 4(q-2)/q^2 nu ||nabla |u|^(q/2)||_2^2                        + (q-2) int p |u|^(q-4) u cdot (u cdot nabla) u
      !         = - nu int |u|^(q-2) |nabla u|^2         - (q-2) nu int |u|^(q-4) u_i partial_j u_i u_k partial_j u_k   + (q-2) int p |u|^(q-4) u cdot (u cdot nabla) u
      !         =                 R_1                    +                              R_2                             +                       R_3
      !         = R
      ! p = Delta^{-1} (nabla u cdot nabla u^T)
      ! returns (/ R, R1+R2, R3 /)  = (/ globalResult, globalViscosityContribution, globalNonlinearPressureContribution  /)
      !=========================================================
      function calc_global_dLqdt_inclParts(u, q) result (global_result_R_Rvisc_Rnonlin)
         use global_variables
         use mpi
         implicit none
         real(pr), intent(in) :: q
         real(pr) :: R_ges, R_visc, R_nonLin
         real(pr), dimension(3) :: local_result, global_result_R_Rvisc_Rnonlin
         real(pr), dimension(1:n(1),1:n(2),1:local_n,1:3), intent(in) :: u
         real(pr), dimension(1:n(1),1:n(2),1:local_n) :: aux_u_q2, aux_u_q4, aux_p, aux_uugradu, aux
         
         local_result = calc_local_dLqdt_inclParts(u,q)
         call mpi_allreduce(local_result, global_result_R_Rvisc_Rnonlin, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)

      end function calc_global_dLqdt_inclParts




      !==========================================
      ! TEST FFT
      !==========================================
      SUBROUTINE testFFT(f)
         USE global_variables
         USE fftwfunction
         IMPLICIT NONE

         REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N), INTENT(IN) :: f

         COMPLEX(pr), DIMENSION(1:n(1),1:n(2),1:local_N) :: aux, faux, test, ftest
         real(pr), DIMENSION(1:n(1),1:n(2),1:local_N) :: aux_r, test_r
         complex(pr), dimension(1:n(1)/2+1,1:n(2),1:local_N) :: faux_r, ftest_r
         INTEGER :: i1
         real(pr) :: result_local, result_global

         aux = dcmplx(f,0.0_pr)
         !aux = cmplx(f,0.0_pr)

         

         CALL fftfwd(aux, faux)
         CALL fftbwd(faux,aux)

         test = aux

         do i1 = 1,19
            CALL fftfwd(test,ftest)
            CALL fftbwd(ftest,test)
         end do

         if (rank==0) then
            print*
            print*, "fft test"
         end if
         result_local = inner_product(f,f,"L2")
         CALL MPI_ALLREDUCE(result_local, result_global, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo) 
         if (rank==0) then
            print*, "L2 norm u            ", result_global
         end if


         if (rank==0) then
            print*
            print*, "complex fft"
         end if
         
         result_local = inner_product(real(aux,pr),real(aux,pr),"L2")
         CALL MPI_ALLREDUCE(result_local, result_global, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo) 
         if (rank==0) then
            print*, "L2 norm f-1(f(u))    ", result_global
         end if

         result_local = inner_product(real(test,pr),real(test,pr),"L2")
         CALL MPI_ALLREDUCE(result_local, result_global, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo) 
         if (rank==0) then
            print*, "L2 norm (f-1f)^20(u) ", result_global
         end if


         aux = f-aux
         result_local = maxval(abs(real(aux,pr)))
         CALL MPI_ALLREDUCE(result_local, result_global, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo) 
         if (rank==0) then
            print*, "max of u-f f-1 u           ", result_global
         end if
         
         test = f-test
         result_local = maxval(abs(real(test,pr)))
         CALL MPI_ALLREDUCE(result_local, result_global, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo) 
         if (rank==0) then
            print*, "max of u-(f f-1)^20 u      ", result_global
         end if

         test = faux-ftest
         result_local = maxval(abs(real(test,pr)))
         CALL MPI_ALLREDUCE(result_local, result_global, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo) 
         if (rank==0) then
            print*, "max of f u - (f f-1)^20 f u", result_global
         end if
   



         if (rank==0) then
            print*
            print*, "real fft"
         end if

         aux_r = f

         CALL fftfwdr(aux_r, faux_r)
         CALL fftbwdr(faux_r, aux_r)

         test_r = aux_r

         do i1 = 1,19
            CALL fftfwdr(test_r, ftest_r)
            CALL fftbwdr(ftest_r, test_r)
         end do



         result_local = inner_product(aux_r,aux_r,"L2")
         CALL MPI_ALLREDUCE(result_local, result_global, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo) 
         if (rank==0) then
            print*, "L2 norm f-1(f(u))    ", result_global
         end if
         
         result_local = inner_product(test_r,test_r,"L2")
         CALL MPI_ALLREDUCE(result_local, result_global, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo) 
         if (rank==0) then
            print*, "L2 norm (f-1f)^20(u) ", result_global
         end if
         
         aux_r = f-aux_r
         result_local = maxval(abs(aux_r))
         CALL MPI_ALLREDUCE(result_local, result_global, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo) 
         if (rank==0) then
            print*, "max of u-f f-1 u           ", result_global
         end if
         
         test_r = f-test_r
         result_local = maxval(abs(test_r))
         CALL MPI_ALLREDUCE(result_local, result_global, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo) 
         if (rank==0) then
            print*, "max of u-(f f-1)^20 u      ", result_global
         end if

         test_r = faux-ftest
         result_local = maxval(abs(test_r))
         CALL MPI_ALLREDUCE(result_local, result_global, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo) 
         if (rank==0) then
            print*, "max of f u - (f f-1)^20 f u", result_global
         end if

      END SUBROUTINE testFFT



      !=========================================================
      ! Calculate W^{1,s} norm
      ! result = (int (|u|^s+|nabla u|^s)^(1/s)
      !                 where |u|^s and |nabla u|^s is calculated using dealiasing
      !                 |nabla u| = (sum_ij |partial_i u_j|^2)^(1/2)
      !=========================================================
      function calc_global_W1s_norm(u,s) result (norm)
         use global_variables
         implicit none

         real(pr), dimension(1:n(1),1:n(2),1:local_n,1:3), intent(in) :: u
         real(pr), dimension(1:n(1),1:n(2),1:local_n,1:3,1:3) :: gradu                 ! gradu(:,:,:,ii,jj) = partial_j u_i
         real(pr), dimension(1:n(1),1:n(2),1:local_n) :: usHalf, nablausHalf
         real(pr), dimension(1:n(1),1:n(2),1:local_n,1:3) :: aux
         real(pr), intent(in) :: s
         real(pr) :: norm1, norm
         integer :: ii


         do ii=1,3
            call gradient(u(:,:,:,ii), aux)
            gradu(:,:,:,ii,:) = aux(:,:,:,:)
         end do

         ! usHalf = |u|^(s/2)
         call calc_uk(u,s/2.0_pr,usHalf)

         ! nablausHalf = |nabla u|^(s/2)
         call calc_nablauk(gradu, s/2.0_pr, nablausHalf)


         ! norm1 = int |u|^s
         norm1 = global_inner_product(usHalf,usHalf,"L2")

         ! norm = int |nabla u|^s
         norm = global_inner_product(nablausHalf, nablausHalf,"L2")

         ! norm = int (|u|^s + |nabla u|^s)
         norm = norm1+norm

         ! norm = (int (|u|^s + |nabla u|^s))^(1/s)
         norm = norm**(1.0_pr/s)
         
      end function

      !=========================================================
      ! Calculate L^q norm
      ! result = (int |u|^q)^(1/q) where |u|^q is calculated using dealiasing
      !=========================================================
      function calc_global_Lq_norm(u) result (norm)
         use global_variables
         implicit none

         real(pr), dimension(1:n(1),1:n(2),1:local_n,1:3), intent(in) :: u
         real(pr), dimension(:,:,:), allocatable :: aux
         real(pr) :: norm
         real(pr) :: inner_prod

         allocate ( aux(1:n(1),1:n(2),1:local_n) )

         call calc_uk(u,lebesgueQ/(2.0_pr),aux)
         inner_prod = global_inner_product(aux,aux,"L2")
         norm = inner_prod**(1.0_pr/lebesgueQ)

         deallocate(aux)
         
      end function

      
      !=========================================================
      ! rescale to ||u||_q = B
      ! u = B u/||u||_q
      !=========================================================
      subroutine rescaleLqNorm(u, B)
         use global_variables
         implicit none
         real(pr), intent(in) :: B
         real(pr), dimension(1:n(1),1:n(2),1:local_N,1:3), intent(inout) :: u
         real(pr) :: LqNorm
         
         LqNorm = calc_global_Lq_norm(u)
         u(:,:,:,:) = B/LqNorm*u(:,:,:,:)

      end subroutine rescaleLqNorm
      
      
      !=========================================================
      ! vector transport
      ! Gamma_{q,B}(u, eta, xi) = B/||u+eta||_q [ xi-(u+eta)||u+eta||_q^{-q} \int |u+eta|^{q-2}(u+eta)cdot xi ]
      !=========================================================
      function vectorTransport(u,q,B,eta,xi) result (resultVec)
         use global_variables
         use mpi
         implicit none
         real(pr), dimension(1:n(1),1:n(2),1:local_n,1:3), intent(in) :: u, eta, xi
         real(pr), dimension(1:n(1),1:n(2),1:local_n,1:3) :: resultVec
         real(pr), dimension(:,:,:,:), allocatable :: uPlusEta, aux
         real(pr), dimension(:,:,:), allocatable :: uPlusEta_Qm2
         real(pr), intent(in) :: B, q
         real(pr) :: Lq_norm, intResult
         integer :: ii

         allocate( uPlusEta(1:n(1),1:n(2),1:local_N,1:3) )
         allocate( aux(1:n(1),1:n(2),1:local_N,1:3) )
         allocate( uPlusEta_Qm2(1:n(1),1:n(2),1:local_N) )

         uPlusEta(:,:,:,:) = u(:,:,:,:) + eta(:,:,:,:)

         call calc_uk(uPlusEta,q-2.0_pr,uPlusEta_Qm2)

         do ii = 1,3
            aux(:,:,:,ii) = uPlusEta_Qm2(:,:,:)*uPlusEta(:,:,:,ii)         ! aux = |u+eta|^{q-2}(u+eta)
            if (toDealias .and. (q-2.0_pr > mach_epsilon)) call dealias_scalar(aux(:,:,:,ii), 2.0_pr)
         end do

         intResult = global_summed_field_inner_product(aux,xi,"L2")              ! \int |u+eta|^{q-2}(u+eta)cdot xi
         
         Lq_norm = calc_global_Lq_norm(uPlusEta)

         resultVec(:,:,:,:) = xi(:,:,:,:) - Lq_norm**(-q)*intResult*uPlusEta(:,:,:,:)  !                 xi-(u+eta)||u+eta||_q^{-q} \int |u+eta|^{q-2}(u+eta)cdot xi
         resultVec(:,:,:,:) = B/Lq_norm*resultVec(:,:,:,:)                                ! B/||u+eta||_q [ xi-(u+eta)||u+eta||_q^{-q} \int |u+eta|^{q-2}(u+eta)cdot xi ]

         deallocate(uPlusEta)
         deallocate(aux)
         deallocate(uPlusEta_Qm2)

      end function vectorTransport



      !=========================================================
      ! Calculate |u|(x) = sqrt(sum_i u_i^2)(x)
      !=========================================================
      function calc_pointwiseVectorNorm(u) result (norm)
         use global_variables
         implicit none

         integer :: ii
         real(pr), dimension(1:n(1),1:n(2),1:local_n,1:3), intent(in) :: u
         real(pr), dimension(1:n(1),1:n(2),1:local_n) :: norm

         norm(:,:,:) = 0.0_pr
         do ii=1,3
            norm(:,:,:) = norm(:,:,:) + u(:,:,:,ii)*u(:,:,:,ii)
         end do
         
         if(toDealias) call dealias_scalar(norm,2.0_pr)

         norm(:,:,:) = sqrt(abs(norm(:,:,:)))
      end function



      !=========================================================
      ! Calculate |M|(x) = sqrt(sum_ij M_ij^2)(x)
      !=========================================================
      function calc_pointwiseFrobeniusNorm(M) result (norm)
         use global_variables
         implicit none

         integer :: ii, jj
         real(pr), dimension(1:n(1),1:n(2),1:local_n,1:3,1:3), intent(in) :: M
         real(pr), dimension(1:n(1),1:n(2),1:local_n) :: norm

         norm(:,:,:) = 0.0_pr
         do jj=1,3
            do ii=1,3
               norm(:,:,:) = norm(:,:,:) + M(:,:,:,ii,jj)*M(:,:,:,ii,jj)
            end do
         end do
         
         if(toDealias) call dealias_scalar(norm,2.0_pr)

         norm(:,:,:) = sqrt(abs(norm(:,:,:)))

      end function

      !=========================================================
      ! Calculate e_u = u/|u|, where |u| is the euklidean norm -> resultVec
      !=========================================================
      function calc_unitVectorInUdirection(u) result (resultVec)
         use global_variables
         implicit none

         integer :: ii, a1, a2, a3
         real(pr), dimension(1:n(1),1:n(2),1:local_n,1:3), intent(in) :: u
         real(pr), dimension(1:n(1),1:n(2),1:local_n,1:3) :: resultVec
         real(pr), dimension(:,:,:), allocatable :: u_norm

         allocate( u_norm(1:n(1),1:n(2),1:local_N) )

         u_norm = calc_pointwiseVectorNorm(u)      

         do ii=1,3
            where (u_norm(:,:,:) < mach_epsilon)      ! uk can be negative because of rounding errors or 1/0 -> results in NaN values
               resultVec(:,:,:,ii) = 0.0_pr
            elsewhere
               resultVec(:,:,:,ii) = u(:,:,:,ii)/u_norm(:,:,:)
            end where
         end do

         do ii=1,3
            do a3=1,local_n
               do a2=1,n(2)
                  do a1=1,n(1)
                     if(isnan(resultVec(a1,a2,a3,ii))) then
                        print*, "unitvector is nan", ii, u(a1,a2,a3,:), u_norm(a1,a2,a3)
                        stop 1
                     end if
                  end do
               end do
            end do
         end do

         deallocate(u_norm)

      end function calc_unitVectorInUdirection



      !=========================================================
      ! Calculate e_M = M/|M|, i.e. {e_M}_ij = M_ij/|M|, where |M| is the Frobenius norm (=L^2 matrix norm) -> resultMat
      !=========================================================
      function calc_unitMatrixInMdirection(M) result (resultMatrix)
         use global_variables
         implicit none

         integer :: ii, jj, a1, a2, a3
         real(pr), dimension(1:n(1),1:n(2),1:local_n,1:3,1:3), intent(in) :: M
         real(pr), dimension(1:n(1),1:n(2),1:local_n,1:3,1:3) :: resultMatrix
         real(pr), dimension(:,:,:), allocatable :: matrix_norm

         allocate( matrix_norm(1:n(1),1:n(2),1:local_N) )

         matrix_norm = calc_pointwiseFrobeniusNorm(M)

         do jj=1,3
            do ii=1,3
               where (matrix_norm(:,:,:) < mach_epsilon)
                  resultMatrix(:,:,:,ii,jj) = 0.0_pr
               elsewhere
                  resultMatrix(:,:,:,ii,jj) = M(:,:,:,ii,jj)/matrix_norm(:,:,:)
               end where
            end do
         end do
            
         do jj=1,3
            do ii=1,3
               do a3=1,local_n
                  do a2=1,n(2)
                     do a1=1,n(1)
                        if(isnan(resultMatrix(a1,a2,a3,ii,jj))) then
                           print*, "resultMatrix is nan", ii, jj, M(a1,a2,a3,:,:), ", matrix_norm", matrix_norm(a1,a2,a3)
                           stop 1
                        end if
                     end do
                  end do
               end do
            end do
         end do

         deallocate(matrix_norm)

      end function calc_unitMatrixInMdirection



      !=========================================================
      ! Calculate nabla u : nabla u^T = partial_i u_j partial_j u_i -> nablaUnablaUt
      !=========================================================
      subroutine calc_nablaUnablaUt(u, nablaUnablaUt)
         use global_variables
         implicit none

         integer :: ii, jj
         real(pr), dimension(1:n(1),1:n(2),1:local_n,1:3), intent(in) :: u
         real(pr), dimension(1:n(1),1:n(2),1:local_n), intent(out) :: nablaUnablaUt
         real(pr), dimension(:,:,:), allocatable :: uComp
         real(pr), dimension(:,:,:,:), allocatable :: gradComp
         real(pr), dimension(:,:,:,:,:), allocatable :: gradU     ! this should be faster than calculating the gradient in a double loop

         allocate( uComp(1:n(1),1:n(2),1:local_N) )
         allocate( gradComp(1:n(1),1:n(2),1:local_N,1:3) )
         allocate( gradU(1:n(1),1:n(2),1:local_N,1:3,1:3) )

         do ii=1,3
            uComp(:,:,:) = u(:,:,:,ii)
            call gradient(uComp, gradComp)
            gradU(:,:,:,ii,:) = gradComp(:,:,:,:)
         end do

         nablaUnablaUt = 0.0_pr
         do ii=1,3
            do jj=1,3
               nablaUnablaUt(:,:,:) = nablaUnablaUt(:,:,:) + gradU(:,:,:,jj,ii)*gradU(:,:,:,ii,jj)
            end do
         end do
         if(toDealias) call dealias_scalar(nablaUnablaUt,2.0_pr)

         deallocate(uComp)
         deallocate(gradComp)
         deallocate(gradU)

      end subroutine calc_nablaUnablaUt

      !=========================================================
      ! Solve (scalar) Poisson Equation: l^2 Delta phi = f -> phi
      !     (l = cut off parameter for nummerical purposes)
      !=========================================================
      subroutine solve_poisson(f, l, phi)
         use global_variables
         use fftwfunction
         implicit none

         integer :: ii, jj, kk
         integer, dimension(3) :: zeroWaveNumberIndex
         real(pr), dimension(1:n(1),1:n(2),1:local_n), intent(in) :: f
         real(pr), dimension(1:n(1),1:n(2),1:local_n), intent(out) :: phi
         real(pr), intent(in) :: l
         real(pr) :: ksq
         complex(pr), dimension(:,:,:), allocatable :: fHat, phiHat                 ! fourier transform of f and phi
         complex(pr), dimension(:,:,:), allocatable :: fComplex, phiComplex         ! f as a complex function fComplex = f + 0*i

         allocate( fHat(1:n(1),1:n(2),1:local_N) )
         allocate( phiHat(1:n(1),1:n(2),1:local_N) )
         allocate( fComplex(1:n(1),1:n(2),1:local_N) )
         allocate( phiComplex(1:n(1),1:n(2),1:local_N) )

         fComplex = dcmplx(f(:,:,:),0.0_pr)
         call fftfwd(fComplex,fHat)

         do kk=1,local_N
            do jj=1,n(2)
               do ii=1,n(1)
                  ksq = K1(ii)**2 + K2(jj)**2 + K3(kk+local_k_offset)**2
                  if (ksq == 0) then
                     phiHat(ii,jj,kk) = 0                                           ! enforces average free condition of phi (and avoids dividing by 0)
                  else
                     phiHat(ii,jj,kk) = - fHat(ii,jj,kk)/((l)**2.0_pr*ksq)
                  end if
               end do
            end do
         end do



         CALL fftbwd(phiHat,phiComplex)
         phi = REAL(phiComplex,pr)
         
         deallocate(fHat)
         deallocate(phiHat)
         deallocate(fComplex)
         deallocate(phiComplex)

      end subroutine solve_poisson
         
      !=========================================================
      ! Calculate v cdot (v cdot nabla) u -> res
      !=========================================================
      function calc_vvgradu(u, v) result (res)
         use global_variables
         implicit none

         real(pr) :: k
         real(pr), dimension(1:n(1),1:n(2),1:local_n,1:3), intent(in) :: u, v
         real(pr), dimension(1:n(1),1:n(2),1:local_n) :: res
         real(pr), dimension(:,:,:), allocatable :: uComponent
         real(pr), dimension(:,:,:,:), allocatable :: gradUComponent
         real(pr), dimension(:,:,:,:), allocatable :: aux
         real(pr), dimension(:,:,:,:,:), allocatable :: gradU         ! gradU(:,:,:,ii,jj) = partial_j u_i
     
         integer :: ii, jj

         allocate( uComponent(1:n(1),1:n(2),1:local_N) )
         allocate( gradUComponent(1:n(1),1:n(2),1:local_N,1:3) )
         allocate( gradU(1:n(1),1:n(2),1:local_N,1:3,1:3) )

         do ii=1,3
            uComponent = u(:,:,:,ii)
            call gradient(uComponent, gradUComponent)
            do jj=1,3
               gradU(:,:,:,ii,jj) = gradUComponent(:,:,:,jj)         ! gradU(:,:,:,ii,jj) = partial_j u_i
               !uugradu(:,:,:) = uugradu(:,:,:) + u(:,:,:,ii)*u(:,:,:,jj)*gradUComponent(:,:,:,jj)
            end do
         end do

         deallocate(uComponent)
         deallocate(gradUComponent)
         allocate( aux(1:n(1),1:n(2),1:local_N,1:3) )

         aux = 0.0_pr
         do jj=1,3
            do ii=1,3
               aux(:,:,:,jj) = aux(:,:,:,jj) + v(:,:,:,ii)*gradU(:,:,:,ii,jj)
            end do
         end do

         do ii=1,3
            if(toDealias) then
               if(.not. use_e_u_instead_of_uqMinus4) then
                  call dealias_scalar(aux(:,:,:,ii),2.0_pr)
               else
                  if(dealiase_if_mult_by_e_u) call dealias_scalar(aux(:,:,:,ii),2.0_pr)
               end if
            end if
         end do
         

         res = 0.0_pr
         do jj=1,3
            res(:,:,:) = res(:,:,:) + v(:,:,:,jj)*aux(:,:,:,jj)
         end do
         if(toDealias) then
            if(.not. use_e_u_instead_of_uqMinus4) then
               call dealias_scalar(res(:,:,:),2.0_pr)
            else
               if(dealiase_if_mult_by_e_u) call dealias_scalar(res(:,:,:),2.0_pr)
            end if
         end if

         deallocate( aux )


      end function calc_vvgradu


      !=========================================================
      ! Calculate u cdot (u cdot nabla) u -> uugradu
      !=========================================================
      subroutine calc_uugradu(U, uugradu)
         use global_variables
         implicit none

         real(pr) :: k
         real(pr), dimension(1:n(1),1:n(2),1:local_n,1:3), intent(in) :: u
         real(pr), dimension(1:n(1),1:n(2),1:local_n), INTENT(OUT) :: uugradu
         real(pr), dimension(:,:,:), allocatable :: uComponent
         real(pr), dimension(:,:,:,:), allocatable :: gradUComponent
         real(pr), dimension(:,:,:,:), allocatable :: aux
         real(pr), dimension(:,:,:,:,:), allocatable :: gradU         ! gradU(:,:,:,ii,jj) = partial_j u_i
     
         integer :: ii, jj

         allocate( uComponent(1:n(1),1:n(2),1:local_N) )
         allocate( gradUComponent(1:n(1),1:n(2),1:local_N,1:3) )
         allocate( gradU(1:n(1),1:n(2),1:local_N,1:3,1:3) )

         do ii=1,3
            uComponent = u(:,:,:,ii)
            call gradient(uComponent, gradUComponent)
            do jj=1,3
               gradU(:,:,:,ii,jj) = gradUComponent(:,:,:,jj)         ! gradU(:,:,:,ii,jj) = partial_j u_i
               !uugradu(:,:,:) = uugradu(:,:,:) + u(:,:,:,ii)*u(:,:,:,jj)*gradUComponent(:,:,:,jj)
            end do
         end do

         deallocate(uComponent)
         deallocate(gradUComponent)
         allocate( aux(1:n(1),1:n(2),1:local_N,1:3) )

         aux = 0.0_pr
         do ii=1,3
            do jj=1,3
               aux(:,:,:,jj) = aux(:,:,:,jj) + u(:,:,:,ii)*gradU(:,:,:,ii,jj)
            end do
         end do

         do ii=1,3
            if(toDealias) call dealias_scalar(aux(:,:,:,ii),2.0_pr)
         end do
         

         uugradu = 0.0_pr
         do jj=1,3
            uugradu(:,:,:) = uugradu(:,:,:) + u(:,:,:,jj)*aux(:,:,:,jj)
         end do
         if(toDealias) call dealias_scalar(uugradu,2.0_pr)


      end subroutine calc_uugradu

      !=========================================================
      ! Calculate v_i partial_j u_i v_k partial_j u_k -> res
      !=========================================================
      function calc_vgraduvgradu(u, v) result (res)
         use global_variables
         implicit none

         real(pr), dimension(1:n(1),1:n(2),1:local_n,1:3), intent(in) :: u, v
         real(pr), dimension(1:n(1),1:n(2),1:local_n) :: res
         real(pr), dimension(:,:,:,:,:), allocatable :: gradU         ! gradU(:,:,:,ii,jj) = partial_j u_i
         real(pr), dimension(:,:,:), allocatable :: uComponent
         real(pr), dimension(:,:,:,:), allocatable :: gradUComponent
         real(pr), dimension(:,:,:,:), allocatable :: vGradu
     
         integer :: ii, jj, kk

         allocate( uComponent(1:n(1),1:n(2),1:local_N) )
         allocate( gradUComponent(1:n(1),1:n(2),1:local_N,1:3) )
         allocate( gradU(1:n(1),1:n(2),1:local_n,1:3,1:3) )
         allocate( vGradu(1:n(1),1:n(2),1:local_n,1:3) )

         do ii=1,3
            uComponent(:,:,:) = u(:,:,:,ii)
            call gradient(uComponent,gradUComponent)
            gradU(:,:,:,ii,:) = gradUComponent(:,:,:,:)         ! gradU(:,:,:,ii,jj) = partial_j u_i
         end do

         vGradu = 0.0_pr
         do jj=1,3
            do ii=1,3
               vGradu(:,:,:,jj) = vGradu(:,:,:,jj) + v(:,:,:,ii)*gradU(:,:,:,ii,jj)
            end do
         end do

         do jj=1,3
            if(toDealias) then
               if(.not. use_e_u_instead_of_uqMinus4) then
                  call dealias_scalar(vGradu(:,:,:,jj),2.0_pr)
               else
                  if(dealiase_if_mult_by_e_u) call dealias_scalar(vGradu(:,:,:,jj),2.0_pr)
               end if
            end if
         end do

         res = 0.0_pr
         do jj=1,3
            res(:,:,:) = res(:,:,:) + vGradu(:,:,:,jj)*vGradu(:,:,:,jj)
         end do
         if(toDealias) call dealias_scalar(res,2.0_pr)
         

         deallocate( uComponent )
         deallocate( gradUComponent )
         deallocate( gradU )
         deallocate( vGradu )

      end function calc_vgraduvgradu

      !=========================================================
      ! Calculate u_i partial_j u_i u_k partial_j u_k -> ugradu2
      !=========================================================
      subroutine calc_ugraduugradu(U, ugradu2)
         use global_variables
         implicit none

         real(pr), dimension(1:n(1),1:n(2),1:local_n,1:3), intent(in) :: u
         real(pr), dimension(1:n(1),1:n(2),1:local_n), INTENT(OUT) :: ugradu2
         real(pr), dimension(:,:,:,:,:), allocatable :: gradU         ! gradU(:,:,:,ii,jj) = partial_j u_i
         real(pr), dimension(:,:,:), allocatable :: uComponent
         real(pr), dimension(:,:,:,:), allocatable :: gradUComponent
         real(pr), dimension(:,:,:,:), allocatable :: uGradu
     
         integer :: ii, jj, kk

         allocate( uComponent(1:n(1),1:n(2),1:local_N) )
         allocate( gradUComponent(1:n(1),1:n(2),1:local_N,1:3) )
         allocate( gradU(1:n(1),1:n(2),1:local_n,1:3,1:3) )
         allocate( uGradu(1:n(1),1:n(2),1:local_n,1:3) )

         do ii=1,3
            uComponent(:,:,:) = u(:,:,:,ii)
            call gradient(uComponent,gradUComponent)
            do jj=1,3
               gradU(:,:,:,ii,jj) = gradUComponent(:,:,:,jj)         ! gradU(:,:,:,ii,jj) = partial_j u_i
            end do
         end do

         uGradu = 0.0_pr
         do ii=1,3
            do jj=1,3
               uGradu(:,:,:,jj) = uGradu(:,:,:,jj) + u(:,:,:,ii)*gradU(:,:,:,ii,jj)
            end do
         end do

         do jj=1,3
            if(toDealias) call dealias_scalar(uGradu(:,:,:,jj),2.0_pr)
         end do



         uGradu2 = 0.0_pr
         do jj=1,3
            uGradu2(:,:,:) = uGradu2(:,:,:) + uGradu(:,:,:,jj)*uGradu(:,:,:,jj)
         end do
         if(toDealias) call dealias_scalar(uGradu2,2.0_pr)
         

         deallocate( uComponent )
         deallocate( gradUComponent )
         deallocate( gradU )
         deallocate( uGradu )

      end subroutine calc_ugraduugradu



      !=========================================================
      ! Calculate f^k -> fk using maximal products of order 2
      !=========================================================
      recursive function calc_gk_order2(g,k) result (gk)
         use global_variables
         implicit none

         real(pr), intent(in) :: k
         real(pr), dimension(1:n(1),1:n(2),1:local_n), intent(in) :: g
         real(pr), dimension(1:n(1),1:n(2),1:local_n) :: gk
         real(pr), dimension(:,:,:), allocatable :: g2, g_km2
         integer :: ii,jj,kk

         if (k+mach_epsilon<0) then
            if(rank == 0 .and. abs(lebesgueQ-2.0_pr)>mach_epsilon .and. dividingByZeroWarnings>0) then
               dividingByZeroWarnings = dividingByZeroWarnings - 1
               if(dividingByZeroWarnings==0) print*, "WARNING: STOPPING THE DIVIDING BY 0 WARNING"
               print*, "warning (in calc_gk_order2) calculating g^{-|x|}, setting y/0 terms to 0"
            end if
            where (abs(g) < mach_epsilon)      ! uk can be negative because of rounding errors or 1/0 -> results in NaN values
               gk = 0.0_pr
            elsewhere
               gk = g**(k)
            end where
         elseif (k+mach_epsilon < 2.0_pr) then
            gk(:,:,:) = g(:,:,:)**(k)
            if(toDealias) call dealias_scalar(gk,k)
         else
            allocate( g2(1:n(1),1:n(2),1:local_N) )
            allocate( g_km2(1:n(1),1:n(2),1:local_N) )

            g2(:,:,:) = g(:,:,:)*g(:,:,:)
            if(toDealias) call dealias_scalar(g2,2.0_pr)
            g_km2 = calc_gk_order2(g,k-2.0_pr)
            gk(:,:,:) = g2(:,:,:) * g_km2(:,:,:)
            if(toDealias) call dealias_scalar(gk,2.0_pr)
            
            deallocate( g2 )
            deallocate( g_km2 )
         end if
         
      end function calc_gk_order2

      !=========================================================
      ! Calculate |u|^k -> uk
      !=========================================================
      subroutine calc_uk(U, k, uk)
         use global_variables
         implicit none

         real(pr), intent(in) :: k
         real(pr), dimension(1:n(1),1:n(2),1:local_n,1:3), intent(in) :: u
         real(pr), dimension(1:n(1),1:n(2),1:local_n), INTENT(OUT) :: uk         
         
         uk(:,:,:) = u(:,:,:,1)**2.0_pr+u(:,:,:,2)**2.0_pr+u(:,:,:,3)**2.0_pr
         if(toDealias) call dealias_scalar(uk, 2.0_pr)

         where (uk <= 0.0_pr)      ! uk can be negative because of rounding errors -> results in NaN values
            uk = 0.0_pr
         elsewhere
            uk = uk**(0.5_pr)
         end where

         uk = calc_gk_order2(uk, k)

      end subroutine calc_uk

      
      !=========================================================
      ! Calculate |nabla u|^k -> nablauk
      !=========================================================
      subroutine calc_nablauk(nablau, k, nablauk)
         use global_variables
         implicit none

         real(pr), intent(in) :: k
         real(pr), dimension(1:n(1),1:n(2),1:local_n,1:3,1:3), intent(in) :: nablau
         real(pr), dimension(1:n(1),1:n(2),1:local_n), INTENT(OUT) :: nablauk         
         
         nablauk = calc_pointwiseFrobeniusNorm(nablau)
         nablauk = calc_gk_order2(nablauk, k)

      end subroutine calc_nablauk


      !=========================================================
      ! Calculate |nabla u|^2 -> nablaU2
      !=========================================================
      subroutine calc_nablaModSquared(u, nablaU2)
         use global_variables
         implicit none

         real(pr), dimension(1:n(1),1:n(2),1:local_n,1:3), intent(in) :: u
         real(pr), dimension(1:n(1),1:n(2),1:local_n), intent(out) :: nablaU2

         real(pr), dimension(1:n(1),1:n(2),1:local_N) :: uTemp
         real(pr), dimension(1:n(1),1:n(2),1:local_N,1:3) :: gradUTemp
     
         integer :: ii, jj

         nablaU2 = 0.0_pr

         do ii=1,3
            uTemp(:,:,:) = u(:,:,:,ii)
            CALL gradient(uTemp, gradUTemp)
            do jj=1,3
               nablaU2(:,:,:) = nablaU2(:,:,:) + gradUTemp(:,:,:,jj)*gradUTemp(:,:,:,jj)
            end do
         end do
         if(toDealias) call dealias_scalar(nablaU2,2.0_pr)

      end subroutine calc_nablaModSquared

      !=========================================================
      ! Calculate kinetic energy from function in physical space
      !=========================================================
      FUNCTION Energy(U) RESULT (kin_ener)
         USE global_variables
         IMPLICIT NONE
                     ! 
         REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(IN) :: U
         INTEGER ::  i1, i2, i3
         REAL(pr), DIMENSION(1:3) :: kin_ener
                     ! 
         kin_ener = 0.0_pr
         DO i3=1,local_N
            DO i2=1,n(2) 
               DO i1=1,n(1)  
                  kin_ener(1) = kin_ener(1) + 0.5_pr*U(i1,i2,i3,1)**2*dV
                  kin_ener(2) = kin_ener(2) + 0.5_pr*U(i1,i2,i3,2)**2*dV
                  kin_ener(3) = kin_ener(3) + 0.5_pr*U(i1,i2,i3,3)**2*dV
               END DO
            END DO
         END DO
      END FUNCTION Energy



      !=========================================================
      ! Calculate enstrophy from function in physical space
      !=========================================================
      FUNCTION Enstrophy(U) RESULT (ens)
         USE global_variables
         IMPLICIT NONE

         REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(IN) :: U
         REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: W
         REAL(pr), DIMENSION(1:3) :: ens
         INTEGER ::  i1, i2, i3
         REAL(pr) :: mode

         ALLOCATE( W(1:n(1),1:n(2),1:local_N,1:3) )
         CALL vel2vort(U,W)
         ens = 0.0_pr
         DO i3 = 1,local_N
            DO i2=1,n(2)
               DO i1 = 1,n(1)
                  ens(1) = ens(1) + 0.5_pr*W(i1,i2,i3,1)**2*dV
                  ens(2) = ens(2) + 0.5_pr*W(i1,i2,i3,2)**2*dV
                  ens(3) = ens(3) + 0.5_pr*W(i1,i2,i3,3)**2*dV
               END DO
            END DO
         END DO
         DEALLOCATE(W)
      END FUNCTION Enstrophy





      !====================================================
      ! CALCULATE HELICITY FROM U AND W
      ! H = /int( U /cdot W)
      !====================================================
      FUNCTION Helicity(U,W) RESULT (H)
         USE global_variables
         IMPLICIT NONE
         REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(IN) :: U,W    
         REAL(pr) :: H

         INTEGER :: ii, jj, kk

         H = 0.0_pr

         DO kk=1,local_N
            DO jj=1,n(2)
               DO ii=1,n(1)
                  H = H + ( U(ii,jj,kk,1)*W(ii,jj,kk,1) + U(ii,jj,kk,2)*W(ii,jj,kk,2) + U(ii,jj,kk,3)*W(ii,jj,kk,3) )*dV
               END DO
            END DO
         END DO

      END FUNCTION Helicity 





      !============================================================
      !--Diagnostics: |U|, |W|, stretching factor
      !============================================================
      SUBROUTINE diagnosticScalars(U)
         USE global_variables
         USE data_ops
         USE fftwfunction
         USE mpi
         IMPLICIT NONE

         REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(INOUT) :: U
         REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3) :: W
         REAL(pr), DIMENSION(:,:), ALLOCATABLE :: Spectrum
         REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: aux
         REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: auxVec, allDiagFields
         REAL(pr), DIMENSION(1:3) :: local_q3, K, E, Umax, Wmax, vorCoreData
         REAL(pr) :: local_q1, dEdt, H, magUmax, magWmax, maxHel, minHel, minStretch, maxStretch, divU_L2norm
         REAL(pr), DIMENSION(1:3) :: myPoint, myNormal
         INTEGER :: nn, i1, i2, i3 
         INTEGER, PARAMETER :: numDiagFields = 4 

         ALLOCATE( Spectrum(1:n(1),1:2) )
         CALL calculate_spectral_data(U, Spectrum)

         call vel2vort(u,w)

                                 !ALLOCATE( auxVec(1:n(1),1:n(2),1:local_nlast,1:3) )
         ALLOCATE( allDiagFields(1:n(1),1:n(2),1:local_N,1:numDiagFields) )

                                 !local_q3 = (/ 1.0_pr, -1.0_pr, 0.0_pr /) 
                                 !myPoint = FindMinPoint(magU,magW, local_q3 )
                                 !IF (myIter==0) THEN           
                                 !   myPoint = (/ REAL(120,pr)/REAL(n(1),pr), REAL(116,pr)/REAL(n(2),pr), REAL(120,pr)/REAL(n(3),pr) /)
                                 !   CALL shift_space(U, myPoint)
                                 !   CALL shift_space(W, myPoint)

                                    !myPoint = (/ PI/2.0_pr, PI/2.0_pr, PI/2.0_pr /)
                                    !CALL rotate_space(U, myPoint)
                                    !CALL vel2vort(U,W) 
                                 !END IF

                                 !CALL advection(W, U, auxVec)
         allDiagFields(:,:,:,1) = SQRT( U(:,:,:,1)**2 + U(:,:,:,2)**2 + U(:,:,:,3)**2 )
         allDiagFields(:,:,:,2) = SQRT( W(:,:,:,1)**2 + W(:,:,:,2)**2 + W(:,:,:,3)**2 )
                     !allDiagFields(:,:,:,3) = ( W(:,:,:,1)*auxVec(:,:,:,1) + W(:,:,:,2)*auxVec(:,:,:,2) + W(:,:,:,3)*auxVec(:,:,:,3) )/( W(:,:,:,1)**2 + W(:,:,:,2)**2 + W(:,:,:,3)**2 ) 
         allDiagFields(:,:,:,3) = U(:,:,:,1)*W(:,:,:,1) + U(:,:,:,2)*W(:,:,:,2) + U(:,:,:,3)*W(:,:,:,3)
                                 !DEALLOCATE( auxVec )

         ALLOCATE( aux(1:n(1),1:n(2),1:local_N) ) 
         CALL divergence(U, aux)          
         local_q1 = inner_product(aux, aux, "L2")
         CALL MPI_ALLREDUCE(local_q1, divU_L2norm, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo) 

         local_q3 = Energy(U)
         CALL MPI_ALLREDUCE(local_q3, K, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
         local_q3 = Energy(W)
         CALL MPI_ALLREDUCE(local_q3, E, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)

         local_q3(1) = MAXVAL(U(:,:,:,1))
         local_q3(2) = MAXVAL(U(:,:,:,2))
         local_q3(3) = MAXVAL(U(:,:,:,3))
         CALL MPI_ALLREDUCE(local_q3, Umax, 3, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)
      
         local_q3(1) = MAXVAL(W(:,:,:,1))        
         local_q3(2) = MAXVAL(W(:,:,:,2))
         local_q3(3) = MAXVAL(W(:,:,:,3))
         CALL MPI_ALLREDUCE(local_q3, Wmax, 3, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)

         local_q1 = Helicity(U, W)
         CALL MPI_ALLREDUCE(local_q1, H, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo) 

         aux = allDiagFields(:,:,:,1)
         local_q1 = MAXVAL(aux)
         CALL MPI_ALLREDUCE(local_q1, magUmax, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)
         aux = allDiagFields(:,:,:,2)
         local_q1 = MAXVAL(aux)
         CALL MPI_ALLREDUCE(local_q1, magWmax, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)
         aux = allDiagFields(:,:,:,3)
         local_q1 = MAXVAL(aux)
         CALL MPI_ALLREDUCE(local_q1, maxHel, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)
         local_q1 = MINVAL(aux)
         CALL MPI_ALLREDUCE(local_q1, minHel, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, Statinfo)

         CALL calculate_ring_data(allDiagFields(:,:,:,2))

         CALL calculate_geometric_data(U, W, aux, "Q", vorCoreData)

         allDiagFields(:,:,:,4) = aux


         CALL save_diagnosticScalars(allDiagFields, numDiagFields, "magU,magW,Helicity,VortexCore")

         IF (rank==0 .and. save_diag_fields_values) THEN
            CALL save_diagnosticFields_global("maxdEdt", K, E, Umax, Wmax, magUmax, magWmax, H, maxHel, minHel, vorCoreData)
         END IF
         CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

         DEALLOCATE( Spectrum )
         DEALLOCATE( aux )
         DEALLOCATE( allDiagFields )

      END SUBROUTINE diagnosticScalars 







      !====================================================
      ! Calculate the spectrum of the velocity field
      ! using spherical shells
      !====================================================
      SUBROUTINE calculate_spectral_data(u, spectral_data)
         USE global_variables
         USE fftwfunction
         USE mpi
         IMPLICIT NONE

         REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(IN) :: u
         REAL(pr), DIMENSION(1:n(1), 1:2), INTENT(OUT) :: spectral_data
         COMPLEX(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: fu
         COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: aux, faux
         INTEGER :: i, i1, i2, i3, mode_count
         REAL(pr) :: kk_min, kk_max, norm_K
         REAL(pr) :: local_spectral_data, global_spectral_data, spectral_Ener, L2norm
         REAL(pr), DIMENSION(1:3) :: local_q3, Ener

         ALLOCATE( fu(1:n(1),1:n(2),1:local_N,1:3) )
         ALLOCATE( aux(1:n(1),1:n(2),1:local_N) )
         ALLOCATE( faux(1:n(1),1:n(2),1:local_N) )


         DO i=1,3
            aux = dcmplx(u(:,:,:,i),0.0_pr)
            CALL fftfwd(aux, faux)
            fu(:,:,:,i) = faux
         END DO

         DO i=0,n(1)-1     ! careful with index out of bounds errors, mpif90 compiler does not care about them!!!
            kk_min = 2.0_pr*PI*REAL(i,pr)
            kk_max = 2.0_pr*PI*REAL(i+1,pr)
            spectral_data(i+1,1) = kk_max
            local_spectral_data = 0.0_pr
            global_spectral_data = 0.0_pr
            mode_count = 0
            DO i3=1,local_N
               DO i2=1,n(2)
                  DO i1=1,n(1)
                     norm_K = SQRT( K1(i1)**2 + K2(i2)**2 + K3(i3+local_k_offset)**2 )
                     IF (  (kk_min < norm_K) .AND. (norm_K <= kk_max) ) THEN
                        local_spectral_data = local_spectral_data + ABS(fu(i1,i2,i3,1))**2 + ABS(fu(i1,i2,i3,2))**2 + ABS(fu(i1,i2,i3,3))**2
                     END IF
                  END DO
               END DO
            END DO
            CALL MPI_ALLREDUCE(local_spectral_data, global_spectral_data, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
            spectral_data(i+1,2) = global_spectral_data
         END DO

         !local_q3 = Energy(u)
         !CALL MPI_ALLREDUCE(local_q3, Ener, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
         !spectral_Ener = 2.0_pr*PI*SUM(spectral_data(:,2))
         !if(rank==0) print*, "spectral factor, should it be 1?", SUM(Ener)/spectral_Ener
         !spectral_data(:,2) = SUM(Ener)/spectral_Ener*spectral_data(:,2)

         if(normalizeSpectrumByL2Norm) then
            L2norm = global_summed_field_inner_product(u,u,"L2")
            spectral_Ener = 2.0_pr*PI*SUM(spectral_data(:,2))
            spectral_data(:,2) = L2norm/spectral_Ener*spectral_data(:,2)
         end if

         DEALLOCATE( fu )
         DEALLOCATE( aux )
         DEALLOCATE( faux )

      END SUBROUTINE calculate_spectral_data

      !====================================================================
      ! OBTAIN GEOMETRIC INFORMATION FROM VORTICITY FIELD
      !====================================================================
      SUBROUTINE calculate_geometric_data(U, W, vortexCore, vortexCriterion, diagScalars)
         USE global_variables
         USE data_ops
         USE mpi
         IMPLICIT NONE

         REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(IN) :: U, W
         REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N), INTENT(OUT) :: vortexCore
         CHARACTER(len=*), INTENT(IN) :: vortexCriterion
         REAL(pr), DIMENSION(1:3), INTENT(OUT) :: diagScalars          

         REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: grad_Ux, grad_Uy, grad_Uz
         REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: aux

         INTEGER, PARAMETER :: lwmax = 1000

         REAL(pr), DIMENSION(1:3,1:3) :: A, B
         REAL(pr), DIMENSION(1:3,1:3) :: dP, vR, vL
         REAL(pr), DIMENSION(1:3) :: wr, wi, lambda
         REAL(pr), DIMENSION(1:lwmax) :: work
         INTEGER :: info, lwork, ii, jj, kk, nn, mm

         REAL(pr) :: local_real, global_real
      
         ALLOCATE( grad_Ux(1:n(1),1:n(2),1:local_N,1:3) )  
         ALLOCATE( grad_Uy(1:n(1),1:n(2),1:local_N,1:3) )  
         ALLOCATE( grad_Uz(1:n(1),1:n(2),1:local_N,1:3) )  
         ALLOCATE( aux(1:n(1),1:n(2),1:local_N) )

         aux = U(:,:,:,1)
         CALL gradient(aux, grad_Ux)
         aux = U(:,:,:,2)
         CALL gradient(aux, grad_Uy) 
         aux = U(:,:,:,3)
         CALL gradient(aux, grad_Uz)
         DEALLOCATE( aux ) 

         local_real = 0.0_pr
         global_real = 0.0_pr

         SELECT CASE (vortexCriterion)

            CASE ("Q")  
               DO kk=1,local_N
                  DO jj=1,n(2)
                     DO ii=1,n(1)
                        vortexCore(ii,jj,kk) = MAX(0.25_pr*( W(ii,jj,kk,1)**2 + W(ii,jj,kk,2)**2 + W(ii,jj,kk,3)**2 ) - &
                                                   0.25_pr*( grad_Ux(ii,jj,kk,2)+grad_Uy(ii,jj,kk,1) )**2 - &
                                                   0.25_pr*( grad_Ux(ii,jj,kk,3)+grad_Uz(ii,jj,kk,1) )**2 - &
                                                   0.25_pr*( grad_Uy(ii,jj,kk,3)+grad_Uz(ii,jj,kk,2) )**2 - &
                                                   0.50_pr*( grad_Ux(ii,jj,kk,1)**2 + grad_Uy(ii,jj,kk,2)**2 + grad_Uz(ii,jj,kk,3)**2 ), 0.0_pr)                         

                        IF ( vortexCore(ii,jj,kk) > 0.0_pr ) THEN
                           local_real = local_real + dV
                        END IF 

                     END DO
                  END DO
               END DO 

         END SELECT

         CALL MPI_ALLREDUCE(local_real, global_real, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo) 
         diagScalars(1) = global_real

         local_real = MAXVAL(vortexCore)
         CALL MPI_ALLREDUCE(local_real, global_real, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)
         diagScalars(2) = global_real

         vortexCore = vortexCore/diagScalars(2)

         DEALLOCATE( grad_Ux )
         DEALLOCATE( grad_Uy )
         DEALLOCATE( grad_Uz )

      END SUBROUTINE calculate_geometric_data 

      !==========================================================
      ! Calculate location of the ring-like vorticity structure
      !==========================================================
      SUBROUTINE calculate_ring_data(magW)
         USE global_variables
         USE data_ops
         USE mpi
         IMPLICIT NONE

         REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N), INTENT(IN) :: magW

         INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: ringMask
         REAL(pr), DIMENSION(:,:), ALLOCATABLE :: local_ringLoc, global_ringLoc
         REAL(pr), DIMENSION(1:3) :: dx 
         REAL(pr) :: local_maxWmag, global_maxWmag
         INTEGER :: mm, nn, ii, jj, kk, myflag
         INTEGER :: local_numPointsRing, global_numPointsRing

         local_maxWmag = MAXVAL(magW)
         CALL MPI_ALLREDUCE(local_maxWmag, global_maxWmag, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, Statinfo)

         dx = 1.0_pr/REAL(n,pr)

         ALLOCATE( ringMask(1:n(1),1:n(2),1:local_N) )

         DO kk=1,local_N
            DO jj=1,n(2)
               DO ii=1,n(1)
                  IF ( magW(ii,jj,kk)/global_maxWmag > 0.95_pr ) THEN
                     ringMask(ii,jj,kk) = 1
                  ELSE
                     ringMask(ii,jj,kk) = 0
                  END IF
               END DO
            END DO
         END DO
                     
         local_numPointsRing = SUM(ringMask)

         ALLOCATE( local_ringLoc(1:local_numPointsRing,3) )
         nn = 1

         DO kk=1,local_N
            DO jj=1,n(2)
               DO ii=1,n(1)
                  IF ( ringMask(ii,jj,kk) == 1 ) THEN
                     local_ringLoc(nn,1) = REAL(ii-1,pr)*dx(1)
                     local_ringLoc(nn,2) = REAL(jj-1,pr)*dx(2)
                     local_ringLoc(nn,3) = REAL(local_k_offset+kk-1,pr)*dx(3)
            
                     nn = nn+1

                  END IF
               END DO
            END DO
         END DO
         CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
                                 !CALL MPI_ALLREDUCE(local_numPointsRing, global_numPointsRing, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, Statinfo)

         myflag = 0
         DO mm=0,np-1 
            IF (rank==mm) THEN
               IF (local_numPointsRing > 10) THEN
                  CALL save_ring_data(local_ringLoc, local_numPointsRing, myflag)
                  myflag = 1
               END IF
            END IF
            CALL MPI_BCAST(myflag, 1, MPI_INTEGER, mm, MPI_COMM_WORLD, Statinfo) 
         END DO

         DEALLOCATE( ringMask )
         DEALLOCATE( local_ringLoc ) 

      END SUBROUTINE calculate_ring_data







      !=====================================================================
      ! OBTAINS THE divAvg_free AND AVG_FREE PORJECTION OF A GIVEN VECTOR FIELD
      !=====================================================================
      SUBROUTINE divAvg_free(myfield)
         USE global_variables
         use fftwfunction
         IMPLICIT NONE

         REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(INOUT) :: myfield
         REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: grad_phi
         REAL(pr), DIMENSION(:,:,:),   ALLOCATABLE :: divU 
         COMPLEX(pr), DIMENSION(:,:,:),   ALLOCATABLE :: aux, faux, divU_hat
         COMPLEX(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: grad_phi_hat, U_hat

         INTEGER  :: i1, i2, i3, ii
         REAL(pr) :: ksq
         REAL(pr), DIMENSION (1:3) :: k
         COMPLEX(pr), DIMENSION(1:3) :: tmp

         ALLOCATE( grad_phi(1:n(1),1:n(2),1:local_N,1:3) )
         ALLOCATE( divU(1:n(1),1:n(2),1:local_N) )
         ALLOCATE( aux(1:n(1),1:n(2),1:local_N) )
         ALLOCATE( faux(1:n(1),1:n(2),1:local_N) )
         ALLOCATE( divU_hat(1:n(1),1:n(2),1:local_N) )
         ALLOCATE( grad_phi_hat(1:n(1),1:n(2),1:local_N,1:3) )
         allocate( U_hat(1:n(1),1:n(2),1:local_N,1:3) )

         CALL divergence(myfield, divU)
         aux = dcmplx(divU,0.0_pr)
         CALL fftfwd(aux,faux)
         divU_hat = faux


         DO ii=1,3
            aux(:,:,:) = dcmplx(myfield(:,:,:,ii), 0.0_pr)
            CALL fftfwd(aux, faux)
            U_hat(:,:,:,ii) = faux
         END DO

         DO i3 = 1, local_N
            DO i2 = 1, n(2)
               DO i1 = 1, n(1)        
                  ksq = K1(i1)**2 + K2(i2)**2 + K3(i3+local_k_offset)**2
                  IF (ksq > MACH_EPSILON) THEN
                     grad_phi_hat(i1,i2,i3,1) = -dcmplx(0.0_pr,K1(i1))*divU_hat(i1,i2,i3)/ksq
                     grad_phi_hat(i1,i2,i3,2) = -dcmplx(0.0_pr,K2(i2))*divU_hat(i1,i2,i3)/ksq
                     grad_phi_hat(i1,i2,i3,3) = -dcmplx(0.0_pr,K3(i3+local_k_offset))*divU_hat(i1,i2,i3)/ksq
                  ELSE
                     grad_phi_hat(i1,i2,i3,1) = 0.0_pr
                     grad_phi_hat(i1,i2,i3,2) = 0.0_pr
                     grad_phi_hat(i1,i2,i3,3) = 0.0_pr 
                     U_hat(i1,i2,i3,1) = 0.0_pr
                     U_hat(i1,i2,i3,2) = 0.0_pr
                     U_hat(i1,i2,i3,3) = 0.0_pr
                  END IF
               END DO 
            END DO
         END DO

         DO ii=1,3
            faux = grad_phi_hat(:,:,:,ii)
            CALL fftbwd(faux,aux)
            grad_phi(:,:,:,ii) = REAL(aux, pr)

            faux = U_hat(:,:,:,ii)
            CALL fftbwd(faux,aux)
            myfield(:,:,:,ii) = REAL(aux, pr)
         END DO

         myfield = myfield - grad_phi

         DEALLOCATE( grad_phi )
         DEALLOCATE( divU )
         DEALLOCATE( aux )
         DEALLOCATE( faux )
         DEALLOCATE( divU_hat )
         DEALLOCATE( grad_phi_hat )
         deallocate( U_hat )
      END SUBROUTINE divAvg_free

      !======================================================
      ! FIX THE ENSTROPHY OF A GIVEN VELOCITY FIELD
      !======================================================
      SUBROUTINE Fix_E0(myfield)
                     !           use, intrinsic :: iso_c_binding                                          ! Commented on May 1st, not sure necessary?
         USE global_variables
         use fftwfunction
         USE mpi
         IMPLICIT NONE
                     !           include 'fftw3-mpi.f03'

         REAL(pr), DIMENSION(1:n(1), 1:n(2), 1:local_N,1:3), INTENT(INOUT) :: myfield
         COMPLEX(pr), DIMENSION(1:n(1), 1:n(2), 1:local_N) :: aux, faux
         REAL(pr), DIMENSION(1:3) :: local_E, global_E
         REAL(pr), DIMENSION(1:3) :: local_K, global_K
         REAL(pr) :: alpha
         INTEGER :: ii

         local_E = Enstrophy(myfield)
         CALL MPI_ALLREDUCE(local_E, global_E, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)
         alpha = SQRT(E0/SUM(global_E))
         myfield = alpha*myfield 

      END SUBROUTINE Fix_E0


      !======================================================
      ! FIX L^q NORM (or equiv the (L^q)^q value) OF A GIVEN VELOCITY FIELD
      !======================================================
      SUBROUTINE Fix_Lq(myfield, targetLq)
         !           use, intrinsic :: iso_c_binding                                          ! Commented on May 1st, not sure necessary?
      USE global_variables
      use fftwfunction
      USE mpi
      IMPLICIT NONE
               !           include 'fftw3-mpi.f03'

      REAL(pr), DIMENSION(1:n(1), 1:n(2), 1:local_N,1:3), INTENT(INOUT) :: myfield
      real(pr), dimension(1:n(1), 1:n(2), 1:local_n) :: u_qHalfs ! calc u^(q/2), potentially using dealiasing
      real(pr), intent(in) :: targetLq

      
      REAL(pr) :: local_Lqq, global_Lqq, global_Lq
      REAL(pr) :: factor
      INTEGER :: ii

      call calc_uk(myfield, lebesgueQ/2.0_pr, u_qHalfs)        ! u_qHalfs = |u|^(q/2)

      local_Lqq = inner_product(u_qHalfs, u_qHalfs, "L2")                ! local int |u|^q
      CALL MPI_ALLREDUCE(local_Lqq, global_Lqq, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)

      global_Lq = global_Lqq**(1.0_pr/lebesgueQ)
      factor = targetLq/global_Lq

      myfield(:,:,:,:) = factor*myfield(:,:,:,:)

      END SUBROUTINE Fix_Lq

      !==============================================================
      ! Calculates the vorticity w in physical space given
      ! velocity u in physical space
      !==============================================================
      SUBROUTINE vel2vort(vel, vort)
                     !           use, intrinsic :: iso_c_binding                                          ! Commentted on May 1st, necessary?
         USE global_variables
         use fftwfunction
         USE mpi
         IMPLICIT NONE

         REAL(pr), DIMENSION(1:n(1), 1:n(2), 1:local_N, 1:3), INTENT(IN)  :: vel
         REAL(pr), DIMENSION(1:n(1), 1:n(2), 1:local_N, 1:3), INTENT(OUT) :: vort
         INTEGER  :: i1, i2, i3
         COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: Ux_hat, Uy_hat, Uz_hat, Wx_hat, Wy_hat, Wz_hat
         COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: aux, faux

         ALLOCATE( Ux_hat(1:n(1),1:n(2),1:local_N) ) 
         ALLOCATE( Uy_hat(1:n(1),1:n(2),1:local_N) ) 
         ALLOCATE( Uz_hat(1:n(1),1:n(2),1:local_N) ) 
         ALLOCATE( Wx_hat(1:n(1),1:n(2),1:local_N) ) 
         ALLOCATE( Wy_hat(1:n(1),1:n(2),1:local_N) ) 
         ALLOCATE( Wz_hat(1:n(1),1:n(2),1:local_N) ) 
         ALLOCATE( aux(1:n(1),1:n(2),1:local_N) ) 
         ALLOCATE( faux(1:n(1),1:n(2),1:local_N) )

         aux(:,:,:) = dcmplx(vel(:,:,:,1), 0.0_pr)
         CALL fftfwd(aux, faux)
         Ux_hat = faux
         aux(:,:,:) = dcmplx(vel(:,:,:,2), 0.0_pr)
         CALL fftfwd(aux, faux)
         Uy_hat = faux
         aux(:,:,:) = dcmplx(vel(:,:,:,3), 0.0_pr)
         CALL fftfwd(aux, faux)
         Uz_hat = faux

         DO i3 = 1,local_N
            DO i2 = 1,n(2)
               DO i1 = 1,n(1)
                  Wx_hat(i1,i2,i3) = dcmplx(0.0_pr, K2(i2))*Uz_hat(i1,i2,i3) - dcmplx(0.0_pr, K3(i3+local_k_offset))*Uy_hat(i1,i2,i3)
                  Wy_hat(i1,i2,i3) = dcmplx(0.0_pr, K3(i3+local_k_offset))*Ux_hat(i1,i2,i3) - dcmplx(0.0_pr, K1(i1))*Uz_hat(i1,i2,i3)
                  Wz_hat(i1,i2,i3) = dcmplx(0.0_pr, K1(i1))*Uy_hat(i1,i2,i3) - dcmplx(0.0_pr, K2(i2))*Ux_hat(i1,i2,i3)
               END DO
            END DO
         END DO
         !==============================================================
         !--Transform back
         !==============================================================
         CALL fftbwd (Wx_hat, aux)
         vort(:,:,:,1) = REAL(aux,pr)
         CALL fftbwd (Wy_hat, aux)
         vort(:,:,:,2) = REAL(aux,pr)
         CALL fftbwd (Wz_hat, aux)
         vort(:,:,:,3) = REAL(aux,pr)

         DEALLOCATE( Ux_hat )
         DEALLOCATE( Uy_hat )
         DEALLOCATE( Uz_hat )
         DEALLOCATE( Wx_hat )
         DEALLOCATE( Wy_hat )
         DEALLOCATE( Wz_hat )
         DEALLOCATE( aux )
         DEALLOCATE( faux )

      END SUBROUTINE vel2vort


      !========================================
      ! CALCULATE 0th fourier modes
      !========================================
      function calc0thFourierModes(fullVec) result (globalResVec)
         USE global_variables
         use mpi
         use fftwfunction
         IMPLICIT NONE
                     ! 
         REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(IN) :: fullVec
         real(pr), dimension(3) :: localResVec, globalResVec
         COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: Ux_hat, Uy_hat, Uz_hat
         COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: aux, faux
         real(pr) :: ksq

         integer :: i1, i2, i3

         ALLOCATE( Ux_hat(1:n(1),1:n(2),1:local_N) ) 
         ALLOCATE( Uy_hat(1:n(1),1:n(2),1:local_N) ) 
         ALLOCATE( Uz_hat(1:n(1),1:n(2),1:local_N) ) 
         ALLOCATE( aux(1:n(1),1:n(2),1:local_N) ) 
         ALLOCATE( faux(1:n(1),1:n(2),1:local_N) )

         aux(:,:,:) = dcmplx(fullVec(:,:,:,1), 0.0_pr)
         CALL fftfwd(aux, faux)
         Ux_hat = faux
         aux(:,:,:) = dcmplx(fullVec(:,:,:,2), 0.0_pr)
         CALL fftfwd(aux, faux)
         Uy_hat = faux
         aux(:,:,:) = dcmplx(fullVec(:,:,:,3), 0.0_pr)
         CALL fftfwd(aux, faux)
         Uz_hat = faux

         localResVec = 0.0_pr
         DO i3 = 1, local_N
            DO i2 = 1, n(2)
               DO i1 = 1, n(1)        
                  ksq = K1(i1)**2 + K2(i2)**2 + K3(i3+local_k_offset)**2
                  if (ksq < MACH_EPSILON) THEN
                     localResVec(1) = Ux_hat(i1,i2,i3)
                     localResVec(2) = Uy_hat(i1,i2,i3)
                     localResVec(3) = Uz_hat(i1,i2,i3)
                  end if
               END DO 
            END DO
         END DO

         call mpi_allreduce(localResVec, globalResVec, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)

      end function


      !========================================
      ! CALCULATE DERIVATIVE WRT ii-th VARIABLE
      !========================================
      SUBROUTINE derivative(u, ii)
         USE global_variables
         USE fftwfunction
         USE mpi
         IMPLICIT NONE

         REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N), INTENT(INOUT) :: u
         INTEGER, INTENT(IN) :: ii

         COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: fu
         COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: aux, faux
         INTEGER :: i1, i2, i3
         REAL(pr) :: mode           

         ALLOCATE( fu(1:n(1),1:n(2),1:local_N) )
         ALLOCATE( aux(1:n(1),1:n(2),1:local_N) )
         ALLOCATE( faux(1:n(1),1:n(2),1:local_N) )

         aux = dcmplx(u,0.0_pr)
         CALL fftfwd(aux, fu)

         DO i3=1, local_N
            DO i2=1,n(2)
               DO i1=1, n(1)
                  IF (ii .EQ. 1) THEN
                     faux(i1,i2,i3) = dcmplx(0.0_pr, K1(i1))*fu(i1,i2,i3)
                  ELSEIF (ii .EQ. 2) THEN
                     faux(i1,i2,i3) = dcmplx(0.0_pr, K2(i2))*fu(i1,i2,i3)
                  ELSEIF (ii .EQ. 3) THEN
                     faux(i1,i2,i3) = dcmplx(0.0_pr, K3(i3+local_k_offset))*fu(i1,i2,i3)
                  END IF
               END DO
            END DO
         END DO
         
         CALL fftbwd(faux, aux)
         u = REAL(aux,pr)

         DEALLOCATE(fu)
         DEALLOCATE(aux)
         DEALLOCATE(faux)

      END SUBROUTINE derivative

      !====================================
      ! CALCULATE GRADIENT OF FUNCTION
      !====================================
      SUBROUTINE gradient(u, grad_u)
         USE global_variables
         use fftwfunction
         USE mpi
         IMPLICIT NONE

         REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N), INTENT(IN) :: u
         REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(OUT) :: grad_u

         COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: fu
         COMPLEX(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: fgrad_u
         COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: aux, faux
         INTEGER :: i1, i2, i3, ii
         REAL(pr) :: mode           

         ALLOCATE( fu(1:n(1),1:n(2),1:local_N) )
         ALLOCATE( fgrad_u(1:n(1),1:n(2),1:local_N,1:3) )
         ALLOCATE( aux(1:n(1),1:n(2),1:local_N) )
         ALLOCATE( faux(1:n(1),1:n(2),1:local_N) )

         aux = dcmplx(u,0.0_pr)


         CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)   
         CALL fftfwd(aux, fu)

         DO i3 = 1, local_N
            DO i2 = 1, n(2)
               DO i1 = 1, n(1)
                  fgrad_u(i1,i2,i3,1) = dcmplx(0.0_pr, K1(i1))*fu(i1,i2,i3)
                  fgrad_u(i1,i2,i3,2) = dcmplx(0.0_pr, K2(i2))*fu(i1,i2,i3) 
                  fgrad_u(i1,i2,i3,3) = dcmplx(0.0_pr, K3(i3+local_k_offset))*fu(i1,i2,i3)
               END DO
            END DO
         END DO

         DO ii=1,3
            faux = fgrad_u(:,:,:,ii)
            CALL fftbwd(faux, aux)
            grad_u(:,:,:,ii) = REAL(aux,pr)
         END DO


         DEALLOCATE(fu)
         DEALLOCATE(fgrad_u)
         DEALLOCATE(aux)
         DEALLOCATE(faux)
         
      END SUBROUTINE gradient
      
      !===========================================
      ! CALCULATE DIVERGENCE OF A FIELD
      !===========================================
      SUBROUTINE divergence(u, divU)
                     !          use, intrinsic :: iso_c_binding                                          ! Commentted on May 1st, necessary?
         USE global_variables
         use fftwfunction
         USE mpi
         IMPLICIT NONE
                     !           include 'fftw3-mpi.f03'

         REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(IN) :: u
         REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N), INTENT(OUT) :: divU

         COMPLEX(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: u_hat
         COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: divU_hat
         COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: aux, faux
         INTEGER :: i1, i2, i3, ii
         REAL(pr) :: local_max_mode, max_mode, mode

         ALLOCATE( u_hat(1:n(1),1:n(2),1:local_N,1:3) )
         ALLOCATE( divU_hat(1:n(1),1:n(2),1:local_N) )
         ALLOCATE( aux(1:n(1),1:n(2),1:local_N) )
         ALLOCATE( faux(1:n(1),1:n(2),1:local_N) )
         
         DO ii=1,3
            aux = dcmplx(u(:,:,:,ii), 0.0_pr)
            CALL fftfwd(aux, faux)
            u_hat(:,:,:,ii) = faux
         END DO

         DO i3 = 1,local_N
            DO i2 = 1,n(2)
               DO i1 = 1,n(1)
                  divU_hat(i1,i2,i3) = dcmplx(0.0_pr, K1(i1))*u_hat(i1,i2,i3,1) + dcmplx(0.0_pr, K2(i2))*u_hat(i1,i2,i3,2) + dcmplx(0.0_pr, K3(i3+local_k_offset))*u_hat(i1,i2,i3,3)
               END DO
            END DO
         END DO

         CALL fftbwd(divU_hat, aux)
         divU = REAL(aux,pr)

         DEALLOCATE( u_hat )
         DEALLOCATE( divU_hat )
         DEALLOCATE( aux )
         DEALLOCATE( faux )
      END SUBROUTINE divergence

      !=======================================
      ! CALCULATE LAPLACIAN OF A VECTOR FIELD
      !=======================================
      SUBROUTINE laplacian(f)
         USE global_variables
         USE fftwfunction
         USE mpi
         IMPLICIT NONE

         REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(INOUT) :: f
         COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: f_hat, aux
         INTEGER :: i1, i2, i3, ii
         REAL(pr) :: ksq, mode, max_mode, local_max_mode
      
         ALLOCATE( f_hat(1:n(1),1:n(2),1:local_N) )
         ALLOCATE( aux(1:n(1),1:n(2),1:local_N) )

         DO ii=1,3
            aux = dcmplx(f(:,:,:,ii), 0.0_pr)
            CALL fftfwd(aux, f_hat)

            DO i3 = 1,local_N
               DO i2 = 1,n(2)
                  DO i1 = 1,n(1)
                     mode = SQRT(K1(i1)**2 + K2(i2)**2 + K3(i3+local_k_offset)**2)
                     ksq = K1(i1)**2 + K2(i2)**2 + K3(i3+local_k_offset)**2
                     !f_hat(i1,i2,i3) = -ksq*f_hat(i1,i2,i3)*EXP(-36.0_pr*(mode/Kcut)**36)
                     f_hat(i1,i2,i3) = -ksq*f_hat(i1,i2,i3)
                  END DO
               END DO
            END DO
            CALL fftbwd(f_hat, aux)
            f(:,:,:,ii) = REAL(aux,pr)
         END DO

         DEALLOCATE( f_hat )
         DEALLOCATE( aux )
      END SUBROUTINE laplacian  

      !==========================================
      ! Function that calculates bilaplacian
      !========================================== 
      SUBROUTINE bilaplacian(u)
         USE global_variables
         USE fftwfunction
         USE mpi
         IMPLICIT NONE

         REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(INOUT) :: u
         COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: fu, aux, faux
         INTEGER :: i1, i2, i3, ii
         REAL(pr) :: k4, mode, max_mode, local_max_mode 

         ALLOCATE( fu(1:n(1),1:n(2),1:local_N) )
         ALLOCATE( aux(1:n(1),1:n(2),1:local_N) )
         ALLOCATE( faux(1:n(1),1:n(2),1:local_N) )

         DO ii=1,3
            aux = dcmplx(u(:,:,:,ii), 0.0_pr)
            CALL fftfwd(aux, faux)
            DO i3 = 1,local_N
               DO i2 = 1,n(2)
                  DO i1 = 1,n(1)
                     mode = SQRT(K1(i1)**2 + K2(i2)**2 + K3(i3+local_k_offset)**2) 
                     k4 = ( K1(i1)**2 + K2(i2)**2 + K3(i3+local_k_offset)**2 )**2
                     !fu(i1,i2,i3) = k4*faux(i1,i2,i3)*EXP(-36.0_pr*(mode/Kcut)**36)
                     fu(i1,i2,i3) = k4*faux(i1,i2,i3)
                  END DO
               END DO
            END DO

            CALL fftbwd(fu, aux)
            u(:,:,:,ii) = REAL(aux,pr)
         END DO 

         DEALLOCATE(fu)
         DEALLOCATE(aux)
         DEALLOCATE(faux) 

      END SUBROUTINE bilaplacian



      !==========================================
      ! PERFORM DEALIASING USING CUT OFF
      ! n_cut = 2/(p+1) n
      !==========================================
      subroutine dealias_scalar(g, nonlinOrder)
         USE global_variables
         USE fftwfunction
         IMPLICIT NONE
         
         real(pr), intent(in) :: nonlinOrder           ! quadratic (advection term) 2, cubic 3, ...
         real(pr), DIMENSION(1:n(1),1:n(2),1:local_N), INTENT(INOUT) :: g
         complex(pr), dimension(1:n(1),1:n(2),1:local_N) :: aux, faux
         logical :: err

         integer :: ii,jj,kk

         aux = dcmplx(g, 0.0_pr)
         call fftfwd(aux,faux)
         call dealiasing_cutoff_scalar_complex(faux,nonlinOrder)
         call fftbwd(faux,aux)
         g = real(aux, pr)

      end subroutine dealias_scalar



      !==========================================
      ! PERFORM DEALIASING USING CUT OFF
      ! n_cut = 2/(p+1) n
      !==========================================
      SUBROUTINE dealiasing_cutoff_scalar_complex(f, nonlinOrder)
         USE global_variables
         USE fftwfunction
         IMPLICIT NONE
         
         complex(pr), DIMENSION(1:n(1),1:n(2),1:local_N), INTENT(INOUT) :: f
         
         real(pr), intent(in) :: nonlinOrder
         INTEGER :: i1, i2, i3
         REAL(pr), DIMENSION(1:3) :: k, k_cut_deal_temp
         real(pr) :: norm_k, k_cut_deal

         !if( rank == 0 .and. nonlinOrder < 0) then
         !   print*, "WARNING DEALIASING ORDER ", nonlinOrder, " < 0"
         !end if
 
         k_cut_deal = 2.0_pr*PI*real(n(1),pr)/(nonlinOrder+1.0_pr)         ! n_cut = 2n/(order+1) >> wavenumber_cut = 4 pi n / (order+1) >> going from -k_max to k_max ( |(-pi,pi)| = 2 pi ) >> |k| < k_cut = 2 pi n / (order + 1)

         if(nonlinOrder > 1.0_pr) then
            if (kmax < 0.0_pr) then
               kmax = k_cut_deal/2.0_pr
            else
               kmax = max(kmax, k_cut_deal/2.0_pr)
            end if
         end if

         DO i3 = 1,local_N
            DO i2 = 1, n(2)
               DO i1 = 1, n(1)
                  norm_k = SQRT(K1(i1)**2 + K2(i2)**2 + K3(i3+local_k_offset)**2)
                  !if( norm_k > k_cut_deal ) then
                  !   f(i1,i2,i3) = dcmplx(0.0_pr,0.0_pr)
                  !end if
                  if( max(abs(K1(i1)),abs(K2(i2)),abs(K3(i3+local_k_offset))) > k_cut_deal ) then
                     f(i1,i2,i3) = dcmplx(0.0_pr,0.0_pr)
                  end if
               END DO
            END DO
         END DO

      END SUBROUTINE dealiasing_cutoff_scalar_complex
   

      !=======================================================
      ! CALCULATE THE INNER PRODUCT BETWEEN TWO (scalar) FUNCTIONS
      !=======================================================
      RECURSIVE FUNCTION inner_product(f, g, mytype) RESULT (inn_prod)
         USE global_variables
         USE fftwfunction
         IMPLICIT NONE
                                 !INCLUDE "mpif.h"
   
         REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N), INTENT(IN) :: f, g
         CHARACTER(len=*), INTENT(IN) :: mytype

         REAL(pr) :: real_loc_prod, inn_prod, ksq
         REAL(pr) :: local_inn_prod , order, factor

         COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: fhat, ghat, aux
         REAL(pr),    DIMENSION(:,:,:), ALLOCATABLE :: f_aux, g_aux
         INTEGER :: i1, i2, i3

         ALLOCATE( fhat(1:n(1),1:n(2),1:local_N) )
         ALLOCATE( ghat(1:n(1),1:n(2),1:local_N) )
         ALLOCATE( aux(1:n(1),1:n(2),1:local_N) )
         ALLOCATE( f_aux(1:n(1),1:n(2),1:local_N) )
         ALLOCATE( g_aux(1:n(1),1:n(2),1:local_N) )

         SELECT CASE (mytype)
         CASE ("L2")
            local_inn_prod = 0.0_pr
            DO i3=1,local_N
               DO i2=1,n(2)
                  DO i1=1,n(1)
                     local_inn_prod = local_inn_prod + f(i1,i2,i3)*g(i1,i2,i3)
                  END DO
               END DO
            END DO

            inn_prod = local_inn_prod*dV

         case ("H^(3/2-1/q)")
            aux = dcmplx(f,0.0_pr)
            CALL fftfwd(aux, fhat)
            aux = dcmplx(g,0.0_pr)
            CALL fftfwd(aux, ghat)
            order = (3.0_pr/2.0_pr)-(1.0_pr/lebesgueQ)
            DO i3=1,local_N
               DO i2=1,n(2)
                  DO i1=1,n(1)
                     ksq = sqrt( K1(i1)**2 + K2(i2)**2 + K3(i3+local_k_offset)**2 )
                     !factor = 1.0_pr + ksq**order
                     factor = (1.0_pr + ksq)**order
                     fhat(i1,i2,i3) = factor * fhat(i1,i2,i3)
                     ghat(i1,i2,i3) = factor * ghat(i1,i2,i3)
                  END DO
               END DO
            END DO
            CALL fftbwd(fhat, aux)
            f_aux = REAL(aux,pr)
            CALL fftbwd(ghat, aux)
            g_aux = REAL(aux,pr)
            inn_prod = inner_product(f_aux,g_aux,"L2")

         case ("H_l^(3/2-1/q)")
            aux = dcmplx(f,0.0_pr)
            CALL fftfwd(aux, fhat)
            aux = dcmplx(g,0.0_pr)
            CALL fftfwd(aux, ghat)
            order = (3.0_pr/2.0_pr)-(1.0_pr/lebesgueQ)
            DO i3=1,local_N
               DO i2=1,n(2)
                  DO i1=1,n(1)
                     ksq = sqrt( K1(i1)**2 + K2(i2)**2 + K3(i3+local_k_offset)**2 )
                     factor = 1.0_pr + (lambda1*ksq)**order
                     !factor = (1.0_pr + lambda1*ksq)**order
                     fhat(i1,i2,i3) = factor * fhat(i1,i2,i3)
                     ghat(i1,i2,i3) = factor * ghat(i1,i2,i3)
                  END DO
               END DO
            END DO
            CALL fftbwd(fhat, aux)
            f_aux = REAL(aux,pr)
            CALL fftbwd(ghat, aux)
            g_aux = REAL(aux,pr)
            inn_prod = inner_product(f_aux,g_aux,"L2")

         case ("H^((3q-1)/(2q))")
            if(rank==0) print*, "ERROR IN function_ops.f90 inner_product case H^((3q-1)/(2q)): you are probably using an old (wrong) exponent"
            call exit
            inn_prod = 0

         case ("H_l^((3q-1)/(2q))")
            if(rank==0) print*, "ERROR IN function_ops.f90 inner_product case H_l^((3q-1)/(2q)): you are probably using an old (wrong) exponent"
            call exit
            inn_prod = 0
         case default
            print*, "WARNING", "TYPE ", mytype, " NOT DEFINED (as inner product)"
         END SELECT

         DEALLOCATE( fhat )
         DEALLOCATE( ghat )
         DEALLOCATE( aux )
         DEALLOCATE( f_aux )
         DEALLOCATE( g_aux )
      END FUNCTION inner_product 

      !=======================================================
      ! CALCULATE THE GLOBAL INNER PRODUCT BETWEEN TWO (scalar) FIELDS
      !=======================================================
      FUNCTION global_inner_product(f,g,mytype) RESULT (inn_prod)
         USE global_variables
         use mpi
         IMPLICIT NONE
   
         REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N), INTENT(IN) :: f,g
         real(pr) :: local_inner, global_inner
         CHARACTER(len=*), INTENT(IN) :: mytype
         REAL(pr) :: inn_prod

         local_inner = inner_product(f,g,mytype)
         call mpi_allreduce(local_inner, inn_prod, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, Statinfo)

      END FUNCTION global_inner_product

      !=======================================================
      ! DEPRECATED   (just use global_summed_field_inner_product)
      ! CALCULATE THE INNER PRODUCT BETWEEN TWO (vector) FIELDS
      !=======================================================
      !FUNCTION field_inner_product(f,g,mytype) RESULT (inn_prod)
      !   USE global_variables
      !   IMPLICIT NONE   
      !   REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(IN) :: f,g
      !   CHARACTER(len=*), INTENT(IN) :: mytype
      !
      !   REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: faux, gaux
      !   REAL(pr), DIMENSION(1:3) :: inn_prod
      !
      !   INTEGER :: ii
      !
      !   ALLOCATE( faux(1:n(1),1:n(2),1:local_N) )
      !   ALLOCATE( gaux(1:n(1),1:n(2),1:local_N) )
      !
      !
      !   inn_prod = 0.0_pr
      !   DO ii=1,3
      !      faux(:,:,:) = f(:,:,:,ii)
      !      gaux(:,:,:) = g(:,:,:,ii)
      !      inn_prod(ii) =  inner_product(faux, gaux, mytype)
      !   END DO
      !
      !   DEALLOCATE(faux)
      !   DEALLOCATE(gaux)       
      !
      !END FUNCTION field_inner_product 

      !=======================================================
      ! CALCULATE THE GLOBAL SUMMED INNER PRODUCT BETWEEN TWO (vector) FIELDS
      !=======================================================
      FUNCTION global_summed_field_inner_product(f,g,mytype) RESULT (inn_prod)
         USE global_variables
         use mpi
         IMPLICIT NONE
   
         REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(IN) :: f,g
         REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N) :: faux,gaux
         integer :: ii
         real(pr), dimension(1:3) :: local_inner, global_inner
         CHARACTER(len=*), INTENT(IN) :: mytype
         REAL(pr) :: inn_prod


         inn_prod = 0.0_pr
         DO ii=1,3
            faux(:,:,:) = f(:,:,:,ii)
            gaux(:,:,:) = g(:,:,:,ii)
            inn_prod = inn_prod + global_inner_product(faux, gaux, mytype)
         END DO

      END FUNCTION global_summed_field_inner_product 

      !======================================================
      ! CALCULATE THE HILBERT SPACE GRADIENT OF ORDER order, GIVEN THE L2 GRADIENT
      ! i.e. Fourier (gradOut) = Fourier (gradIn) / (1+l|xi|)^(2*order)
      !======================================================
      SUBROUTINE HilbertGradient(grad, order)
         USE global_variables
         USE fftwfunction
         IMPLICIT NONE

         REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(INOUT) :: grad
         real(pr), INTENT(IN) :: order

         COMPLEX(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: grad_hat
         COMPLEX(pr), DIMENSION(:,:,:), ALLOCATABLE :: aux, faux

         REAL(pr) :: ksq
         INTEGER :: ii,jj,kk

         ALLOCATE( grad_hat(1:n(1),1:n(2),1:local_N,1:3) )
         ALLOCATE( aux(1:n(1),1:n(2),1:local_N) )
         ALLOCATE( faux(1:n(1),1:n(2),1:local_N) )

         DO ii=1,3
            aux = dcmplx(grad(:,:,:,ii),0.0_pr)
            CALL fftfwd(aux,faux)
            grad_hat(:,:,:,ii) = faux
         END DO

         DO kk=1,local_N
            DO jj=1,n(2)
               DO ii=1,n(1)
                  ksq = SQRT( K1(ii)**2 + K2(jj)**2 + K3(kk+local_k_offset)**2 )
                  grad_hat(ii,jj,kk,:) = grad_hat(ii,jj,kk,:)/( 1.0_pr + (lambda1*ksq)**(2.0_pr*order) )
                  !grad_hat(ii,jj,kk,:) = grad_hat(ii,jj,kk,:)/( (1.0_pr + lambda1*ksq)**(2.0_pr*order) )
               END DO
            END DO
         END DO

         DO ii=1,3
            faux = grad_hat(:,:,:,ii)
            CALL fftbwd(faux,aux)
            grad(:,:,:,ii) = REAL(aux,pr)
         END DO 

         DEALLOCATE( grad_hat )
         DEALLOCATE( aux )
         DEALLOCATE( faux )

      END SUBROUTINE HilbertGradient


      !======================================================
      ! NOT IMPLEMENTED BECAUSE NOT FEASIBLE FOR REASONABLE n values
      ! NEEDS TO SOLVE Ax = b where A is a (3*n^3) times (3*n^3) matrix
      ! CALCULATE THE NEXT ITERATION OF THE BANACH GRADIENT
      ! given the L2 gradient, the previous iteration (v_old), lambda (the lagrange multiplier fixing the norm), rho (the lagrange multiplier fixing the divergence free condition)
      !======================================================
      function BanachGradientIteration(l2Grad, v_old, lambda, rho, tau) result (v_new)
         USE global_variables
         USE fftwfunction
         IMPLICIT NONE

         real(pr), dimension(1:n(1),1:n(2),1:local_N,1:3), intent(in) :: l2Grad, v_old
         real(pr), dimension(1:n(1),1:n(2),1:local_N), intent(in) :: rho ! lagrange multiplier for divergence free
         real(pr), intent(in) :: lambda   ! lagrange multiplier for norm
         real(pr), intent(in) :: tau      ! step size

         real(pr) :: s                    ! s = 3q/(q+1)
         real(pr) :: X, Y, Z
         real(pr), dimension(1:3) :: dx
         real(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3) :: v_new


         real(pr), dimension(1:n(1),1:n(2),1:local_N,1:3) :: rhs
         real(pr), dimension(1:n(1),1:n(2),1:local_N,1:3,1:3) :: matrix_f, test_matrix
         real(pr), dimension(1:n(1),1:n(2),1:local_N,1:3,1:3,1:3) :: partial_sG_sjml
         real(pr), dimension(1:n(1),1:n(2),1:local_N,1:3,1:3,1:3,1:3) :: tensor_g
         real(pr) :: norm ! to normalize

         complex(pr), dimension(:,:,:,:), allocatable :: rhs_hat

         real(pr), dimension(:,:,:), allocatable :: tempSca
         real(pr), dimension(:,:,:,:), allocatable :: tempVec
         real(pr), dimension(:,:,:,:,:,:), allocatable :: temp_tensor3
         
         complex(pr), dimension(:,:,:,:,:), allocatable :: temp_matrix, matrix_1, matrix_1_inverted
         complex(pr), dimension(:,:,:,:), allocatable :: temp_vector, ik_vector
         complex(pr), dimension(:,:), allocatable :: basis_function_x, basis_function_y, basis_function_z ! each 1/n_j e^{ik_j x_j) to have 3*res^2 instead of res^6 arrays

         integer :: ii, jj, kk, ll, mm

         complex(pr), dimension(:,:,:), allocatable :: aux,faux

         


         !REAL(pr) :: ksq
         !INTEGER :: ii,jj,kk

         s = 3.0_pr*lebesgueQ/(lebesgueQ+1.0_pr)

         !if(rank==0) print*, "tau", tau
         call BanachGradient_calcFmatrixGtensorRhs(l2Grad, v_old, lambda, rho, matrix_f, tensor_g, rhs, tau)

         !!! define basis functions !!!

         dx = 1.0_pr/REAL(n, pr)
      
         allocate(basis_function_x(1:n(1),1:n(1)))
         allocate(basis_function_y(1:n(2),1:n(2)))
         allocate(basis_function_z(1:local_n,1:local_n))
         ! basis_func_x = 1/n_x e^{ik_x x}           b_f_x(ii,jj) = ... x(ii) k(jj) ...
         do ii=1,n(1)
            X = REAL(ii-1,pr)*dx(1)
            do jj=1,n(1)
               basis_function_x(ii,jj) = 1.0_pr/real(n(1),pr)*exp(dcmplx(0.0_pr, 1.0_pr)*(K1(jj)*X))
            end do
         end do
         ! basis_func_y = 1/n_y e^{ik_y y}           b_f_y(ii,jj) = ... y(ii) k(jj) ...
         do ii=1,n(2)
            Y = REAL(ii-1,pr)*dx(2)
            do jj=1,n(2)
               basis_function_y(ii,jj) = 1.0_pr/real(n(2),pr)*exp(dcmplx(0.0_pr, 1.0_pr)*(K2(jj)*Y))
            end do
         end do
         ! basis_func_z = 1/n_y e^{ik_z z}           b_f_z(ii,jj) = ... z(ii) k(jj) ...
         do ii=1,local_N
            Z = REAL(ii+local_k_offset-1,pr)*dx(3)
            do jj=1,local_N
               basis_function_z(ii,jj) = 1.0_pr/real(n(3),pr)*exp(dcmplx(0.0_pr, 1.0_pr)*(K3(jj+local_k_offset)*Z))
            end do
         end do






         allocate(temp_tensor3(1:n(1),1:n(2),1:local_N,1:3,1:3,1:3))
         allocate(tempVec(1:n(1),1:n(2),1:local_N,1:3))
         allocate(tempSca(1:n(1),1:n(2),1:local_N))

         do ll=1,3
            do kk=1,3
               do jj=1,3
                     tempVec(:,:,:,:) = tensor_g(:,:,:,:,jj,kk,ll)
                     call divergence(tempVec,tempSca)
                     temp_tensor3(:,:,:,jj,kk,ll) = tempSca       ! temp_tensor3(:,:,:,jj,kk,ll) = partial_i g_ijkl
                     partial_sG_sjml(:,:,:,jj,kk,ll) = tempSca
               end do
            end do
         end do


         deallocate(tempSca)
         deallocate(tempVec)

         ! temp_vector_j = i k_j
         allocate(temp_vector(1:n(1),1:n(2),1:local_N,1:3))
         DO kk = 1, local_N
            DO jj = 1, n(2)
               DO ii = 1, n(1)
                  temp_vector(ii,jj,kk,1) = dcmplx(0.0_pr, K1(ii))
                  temp_vector(ii,jj,kk,2) = dcmplx(0.0_pr, K2(jj))
                  temp_vector(ii,jj,kk,3) = dcmplx(0.0_pr, K3(kk+local_k_offset))
               END DO
            END DO
         END DO

         ik_vector = temp_vector

         ! tempmatrix_kl = i partial_s g_sjkl k_j
         allocate(temp_matrix(1:n(1),1:n(2),1:local_N,1:3,1:3))
         temp_matrix = 0.0_pr
         do ll=1,3
            do kk=1,3
               do jj=1,3
                  temp_matrix(:,:,:,kk,ll) = temp_matrix(:,:,:,kk,ll) + temp_tensor3(:,:,:,jj,kk,ll)*temp_vector(:,:,:,jj)
                  ! no dealiasing since in gradient f -> i k f_hat there is also no dealiasing
               end do
            end do
         end do


         allocate(matrix_1(1:n(1),1:n(2),1:local_N,1:3,1:3))
         ! matrix_1 = f_lm - i partial_s g_sjml k_j + ...(later)...
         matrix_1(:,:,:,:,:) = matrix_f(:,:,:,:,:) - temp_matrix(:,:,:,:,:)


         ! temp_tensor3_jml = g_sjml i k_s
         temp_tensor3 = 0.0_pr
         do ll=1,3
            do kk=1,3
               do jj=1,3
                  do ii=1,3
                     temp_tensor3(:,:,:,jj,kk,ll) = temp_tensor3(:,:,:,jj,kk,ll) + tensor_g(:,:,:,ii,jj,kk,ll)*temp_vector(:,:,:,ii)
                     ! no dealiasing since in gradient f -> i k f_hat there is also no dealiasing
                  end do
               end do
            end do
         end do

         ! temp_matrix = g_sjml (i k_s) (i k_j) = temp_tensor3_jml (i k_j)
         temp_matrix = 0.0_pr
         do ll=1,3
            do kk=1,3
               do jj=1,3
                  do ii=1,3
                     temp_matrix(:,:,:,kk,ll) = temp_matrix(:,:,:,kk,ll) + temp_tensor3(:,:,:,jj,kk,ll)*temp_vector(:,:,:,jj)
                     ! no dealiasing since in gradient f -> i k f_hat there is also no dealiasing
                  end do
               end do
            end do
         end do


         deallocate(temp_tensor3)


         ! matrix_1 = f_lm - i partial_s g_sjml k_j + g_sjml k_s k_j
         !             = f_lm - i partial_s g_sjml k_j - g_sjml (i k_s) (i k_j)
         matrix_1(:,:,:,:,:) = matrix_1(:,:,:,:,:) - temp_matrix(:,:,:,:,:)


         temp_vector = solve_hugh_system(matrix_f, partial_sG_sjml, tensor_g, ik_vector, basis_function_x, basis_function_y, basis_function_z, dcmplx(rhs,0.0_pr))
      
         
         call mpi_barrier(mpi_comm_world, statinfo)
         
         allocate(aux(1:n(1),1:n(2),1:local_N))
         allocate(faux(1:n(1),1:n(2),1:local_N))
         do ii=1,3
            faux(:,:,:) = temp_vector(:,:,:,ii)
            call fftbwd(faux,aux)
            v_new(:,:,:,ii) = real(aux,pr)
         end do

         norm = calc_global_W1s_norm(v_new, s)
         v_new(:,:,:,:) = v_new(:,:,:,:)/norm


         !!! testing !!!

         allocate(tempvec(1:n(1),1:n(2),1:local_n,1:3))

         temp_vector = 0.0_pr

         do ii=1,3
            call gradient(v_new(:,:,:,ii), tempVec)
            matrix_1(:,:,:,ii,:) = tempVec(:,:,:,:)
         end do

         
         temp_matrix = 0.0_pr
         do ll=1,3
            do kk=1,3
               do jj=1,3
                  do ii=1,3
                     ! temp matrix(l,k) = g_(l,j,i,k)*partial_j u_i
                     temp_matrix(:,:,:,ll,kk) = temp_matrix(:,:,:,ll,kk) + tensor_g(:,:,:,ll,jj,ii,kk)*matrix_1(:,:,:,ii,jj)
                  end do
               end do
            end do
         end do
         

         allocate(tempSca(1:n(1),1:n(2),1:local_n))

         temp_vector = 0.0_pr
         do ii=1,3
            do jj=1,3
               tempVec(:,:,:,jj) = temp_matrix(:,:,:,jj,ii)
               call divergence(tempVec, tempSca)
            end do
            temp_vector(:,:,:,ii) = tempSca(:,:,:)
         end do
         tempVec = temp_vector

         temp_vector = 0.0_pr
         do ii=1,3
            do jj=1,3
               temp_vector(:,:,:,ii) = temp_vector(:,:,:,ii) + matrix_f(:,:,:,jj,ii)*v_new(:,:,:,jj)
            end do
         end do
         
         temp_vector = temp_vector+ tempVec
         tempVec = temp_vector-rhs

         if(rank==0) print*, "testing if iterationStep is correct: the following should be 0: ", achar(9), global_summed_field_inner_product(tempVec,tempVec,"L2") 

         !!! end testing !!!
         
         deallocate(aux)
         deallocate(faux)

      end function BanachGradientIteration


      function solve_hugh_system(matrix_A, partial_sG_sjml, g_sjml, ik_vector, phi_x, phi_y, phi_z, rhs) result(resultVec)
         use global_variables
         use fftwfunction
         implicit none

         complex(pr), dimension(1:n(1),1:n(2),1:local_N,1:3,1:3):: matrix_1
         
         real(pr), dimension(1:n(1),1:n(2),1:local_N,1:3,1:3,1:3,1:3), intent(in) :: g_sjml
         real(pr), dimension(1:n(1),1:n(2),1:local_N,1:3,1:3,1:3), intent(in) :: partial_sG_sjml
         real(pr), dimension(1:n(1),1:n(2),1:local_N,1:3,1:3), intent(in) :: matrix_A
         complex(pr), dimension(1:n(1),1:n(1)), intent(in) :: phi_x
         complex(pr), dimension(1:n(2),1:n(2)), intent(in) :: phi_y
         complex(pr), dimension(1:local_n,1:local_n), intent(in) :: phi_z
         complex(pr), dimension(1:n(1),1:n(2),1:local_N,1:3), intent(in) :: rhs, ik_vector
         complex(pr), dimension(1:n(1),1:n(2),1:local_N,1:3) :: resultVec

         complex(pr), dimension(:,:,:,:,:,:,:,:), allocatable :: hugh_tensor   ! tensor(x,y,z,k_x,k_y,k_z,l,m)
         complex(pr), dimension(:,:), allocatable :: hugh_matrix   ! tensor(xyzl,k_xk_yk_zlm)
         complex(pr), dimension(:), allocatable :: hugh_rhs, hugh_x, useless_s, rhs_real, rhs_imag
         real(pr), dimension(:,:), allocatable :: hugh_matrix_real_re, hugh_matrix_real_im
         integer :: rank_out = 0
         integer, dimension(:), allocatable :: useless_ipiv


         !!! delete this afterwards !!!
         complex(pr), dimension(:,:,:,:,:,:,:,:,:), allocatable :: k_sG_sjml
         complex(pr), dimension(:,:,:,:,:,:,:,:), allocatable :: myMatrixMixedDerivatives, k_sk_jg_sjml
         complex(pr), dimension(:,:,:,:), allocatable :: myTemp3p1d
         complex(pr), dimension(:,:), allocatable :: myTemp2d, myTempHughMatrixBackup
         complex(pr), dimension(:), allocatable :: myTemp1d, myTemp1d2, myTempHughRHSbackup
         complex(pr) :: myTempCheckVar
         real(pr) :: temptemp, temptemp2

         real(pr), dimension(:), allocatable :: uselessHoch5
         !!! end delete his afterwards !!!

         integer :: ii, jj, kk, ll, mm, nn, oo, pp, qq, rr, indL, indM
         integer :: hughDim
         integer :: info

         hughDim = n(1)*n(2)*local_n*3
         if(rank==0) print*, "hughDim = ", hughDim
         allocate(hugh_matrix(1:hughDim,1:hughDim))
         allocate(hugh_tensor(1:n(1),1:n(2),1:local_N,1:n(1),1:n(2),1:local_N,1:3,1:3))
         

         call mpi_barrier(mpi_comm_world, statinfo)
         if(rank==0) print*, "starting hugh matrix assembly ..."
         call mpi_barrier(mpi_comm_world, statinfo)

         allocate(k_sG_sjml(1:n(1),1:n(2),1:local_n,1:n(1),1:n(2),1:local_n,1:3,1:3,1:3))
         k_sG_sjml = 0.0_pr
         do rr=1,3
            do qq=1,3
               do pp=1,3
                  do oo=1,3
                     do nn=1,local_n
                        do mm=1,n(2)
                           do ll=1,n(1)
                              do kk=1,local_n
                                 do jj=1,n(2)
                                    do ii=1,n(1)
                                       k_sG_sjml(ii,jj,kk,ll,mm,nn,oo,pp,qq) = k_sG_sjml(ii,jj,kk,ll,mm,nn,oo,pp,qq) + ik_vector(ll,mm,nn,oo)*g_sjml(ii,jj,kk,oo,pp,qq,rr)
                                    end do
                                 end do
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do 
         
         allocate(myMatrixMixedDerivatives(1:n(1),1:n(2),1:local_n,1:n(1),1:n(2),1:local_n,1:3,1:3))
         myMatrixMixedDerivatives = 0.0_pr
         allocate(k_sk_jg_sjml(1:n(1),1:n(2),1:local_n,1:n(1),1:n(2),1:local_n,1:3,1:3))
         k_sk_jg_sjml = 0.0_pr

         do qq=1,3
            do pp=1,3
               do oo=1,3
                  do nn=1,local_n
                     do mm=1,n(2)
                        do ll=1,n(1)
                           do kk=1,local_n
                              do jj=1,n(2)
                                 do ii=1,n(1)
                                    myMatrixMixedDerivatives(ii,jj,kk,ll,mm,nn,oo,pp) = myMatrixMixedDerivatives(ii,jj,kk,ll,mm,nn,oo,pp) + ik_vector(ll,mm,nn,qq)*partial_sG_sjml(ii,jj,kk,qq,oo,pp)
                                    k_sk_jg_sjml(ii,jj,kk,ll,mm,nn,oo,pp) = k_sk_jg_sjml(ii,jj,kk,ll,mm,nn,oo,pp) + ik_vector(ll,mm,nn,qq)*k_sG_sjml(ii,jj,kk,ll,mm,nn,qq,oo,pp)
                                 end do
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do


         do pp=1,3
            do oo=1,3
               do nn=1,local_n
                  do mm=1,n(2)
                     do ll=1,n(1)
                        do kk=1,local_n
                           do jj=1,n(2)
                              do ii=1,n(1)
                                 indL = ii + (jj-1)*n(1) + (kk-1)*n(2)*n(1) + (oo-1)*local_n*n(2)*n(1)
                                 indM = ll + (mm-1)*n(1) + (nn-1)*n(2)*n(1) + (pp-1)*local_n*n(2)*n(1)
                                 hugh_tensor(ii,jj,kk,ll,mm,nn,oo,pp) = matrix_A(ii,jj,kk,oo,pp) + myMatrixMixedDerivatives(ii,jj,kk,ll,mm,nn,oo,pp) + k_sk_jg_sjml(ii,jj,kk,ll,mm,nn,oo,pp)
                                 call random_number(temptemp)
                                 call random_number(temptemp2)
                                 !hugh_tensor(ii,jj,kk,ll,mm,nn,oo,pp) = dcmplx(temptemp,temptemp2)
                                 hugh_matrix(indL, indM) = hugh_tensor(ii,jj,kk,ll,mm,nn,oo,pp)
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do

         allocate(hugh_matrix_real_re(1:hughDim,1:hughDim))
         allocate(hugh_matrix_real_im(1:hughDim,1:hughDim))
         hugh_matrix_real_re(:, :) = real(hugh_matrix(:,:))
         hugh_matrix_real_im(:, :) = aimag(hugh_matrix(:,:))



         allocate(myTempHughMatrixBackup(1:hughDim,1:hughDim))

         myTempHughMatrixBackup = hugh_matrix


         call mpi_barrier(mpi_comm_world, statinfo)
         if(rank==0) print*, achar(9), "done"
         call mpi_barrier(mpi_comm_world, statinfo)
         

         allocate(hugh_rhs(1:hughDim))
         allocate(hugh_x(1:hughDim))


         do oo=1,3
            do kk=1,local_n
               do jj=1,n(2)
                  do ii=1,n(1)
                     indL = ii + (jj-1)*n(1) + (kk-1)*n(2)*n(1) + (oo-1)*local_n*n(2)*n(1)
                     hugh_rhs(indL) = rhs(ii,jj,kk,oo)
                  end do
               end do
            end do
         end do

         myTempHughRHSbackup = hugh_rhs
         rhs_real = hugh_rhs
         rhs_imag = hugh_rhs



         allocate(useless_ipiv(1:hughDim))


         call mpi_barrier(mpi_comm_world, statinfo)
         if(rank==0) print*, "starting zgesv ..."
         call mpi_barrier(mpi_comm_world, statinfo)
         
         
         call zgesv(hughDim, 1, hugh_matrix, hughDim, useless_ipiv, hugh_rhs, hughDim, info)

         call mpi_barrier(mpi_comm_world, statinfo)
         if(rank==0) then
            if (info==0) then
               print*, achar(9), "complex done"
            else
               print*, achar(9), "complex done with errors! zgesv error info ", info
            end if
         end if
         call mpi_barrier(mpi_comm_world, statinfo)


         
         call dgesv(hughDim, 1, hugh_matrix_real_re, hughDim, useless_ipiv, rhs_real, hughDim, info)


         call mpi_barrier(mpi_comm_world, statinfo)
         if(rank==0) then
            if (info==0) then
               print*, achar(9), "real part done"
            else
               print*, achar(9), "real part done with errors! dgesv error info ", info
            end if
         end if
         call mpi_barrier(mpi_comm_world, statinfo)

         
         call dgesv(hughDim, 1, hugh_matrix_real_im, hughDim, useless_ipiv, rhs_imag, hughDim, info)

         call mpi_barrier(mpi_comm_world, statinfo)
         if(rank==0) then
            if (info==0) then
               print*, achar(9), "imag part done"
            else
               print*, achar(9), "imag part done with errors! dgesv error info ", info
            end if
         end if
         call mpi_barrier(mpi_comm_world, statinfo)

         ! condition number
         call mpi_barrier(mpi_comm_world, statinfo)
         if(rank==0) print*, "starting zgecon ..."
         call mpi_barrier(mpi_comm_world, statinfo)
         
         allocate(useless_s(1:2*hughDim))
         
         allocate(uselessHoch5(1:2*hughDim))
         temptemp = 0.0_pr
         do jj = 1,hughDim
            do ii = 1,hughDim
               temptemp = temptemp + abs(hugh_matrix(ii,jj))
               !if(abs(hugh_matrix(ii,jj))>temptemp) temptemp = (abs(hugh_matrix(ii,jj)))
            end do         
         end do
         call zgecon('1',hughDim, hugh_matrix, hughDim, temptemp, temptemp2, useless_s, uselessHoch5, info)
         call mpi_barrier(mpi_comm_world, statinfo)
         if(rank==0) then
            if (info==0) then
               print*, achar(9), "done, condition number ", 1.0_pr/temptemp2, "      norm(A) ", temptemp
            else
               print*, achar(9), "done with errors! zgecon error info ", info
            end if
         end if
         call mpi_barrier(mpi_comm_world, statinfo)

         temptemp2 = 0.0_pr
         do ii = 1,20
            print*, hugh_rhs(ii), rhs_real(ii), rhs_imag(ii)
            !temptemp2 = temptemp2 + abs(hugh_rhs(ii)-dcmplx())
         end do


         if(rank==0) print*, "starting vector reconstruction from hugh vector"

         do oo=1,3
            do kk=1,local_n
               do jj=1,n(2)
                  do ii=1,n(1)
                     indL = ii + (jj-1)*n(1) + (kk-1)*n(2)*n(1) + (oo-1)*local_n*n(2)*n(1)
                     resultVec(ii,jj,kk,oo) = hugh_rhs(indL)
                  end do
               end do
            end do
         end do

         if(rank==0) print*, achar(9), "done"


         !if(rank==0) print*, "starting hugh result testing ..."
         allocate(myTemp1d(1:hughDim))
         allocate(myTemp1d2(1:hughDim))
         
         !myTemp1d = 0.0_pr
         !myTemp1d2 = 0.0_pr
         !do jj=1,hughDim
         !   do ii=1,hughDim
         !      myTemp1d(ii) = myTemp1d(ii) + myTempHughMatrixBackup(ii,jj)*hugh_rhs(jj)
         !      !myTemp1d2(ii) = myTemp1d2(ii) + hugh_matrix(ii,jj)+hugh_rhs(jj)
         !   end do
         !end do
         !print*, achar(9), achar(9), "myTemp1d(ii) ", "myTempHughRHSbackup(ii)"
         !do ii=1,3
         !   print*, achar(9), achar(9), myTemp1d(ii), myTempHughRHSbackup(ii)
         !end do
         !print*, achar(9), achar(9), "maxval(abs(myTempHughMatrixBackup(:,:))) ", maxval(abs(myTempHughMatrixBackup(:,:))), ",   minval(abs(myTempHughMatrixBackup(:,:))) ", minval(abs(myTempHughMatrixBackup(:,:)))
         !myTempCheckVar = maxval(abs(myTemp1d(:)-myTempHughRHSbackup(:)))
         !print*, achar(9), achar(9), "Check L^infty error HughMatrix ", myTempCheckVar, " =0?"

         
         deallocate(myTemp1d)
         deallocate(myTemp1d2)

         allocate(myTemp3p1d(1:n(1),1:n(2),1:local_n,1:3))
         
         myTemp3p1d = 0.0_pr

         do pp=1,3
            do oo=1,3
               do nn=1,local_n
                  do mm=1,n(2)
                     do ll=1,n(1)
                        do kk=1,local_n
                           do jj=1,n(2)
                              do ii=1,n(1)
                                 myTemp3p1d(ii,jj,kk,oo) = myTemp3p1d(ii,jj,kk,oo) + hugh_tensor(ii,jj,kk,ll,mm,nn,oo,pp)*resultVec(ll,mm,nn,pp)
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do


         print*, achar(9), achar(9), "myTemp3p1d(ii,ii,ii,ii) ", "rhs(ii,ii,ii,ii) "
         do ii=1,3
            print*, achar(9), achar(9), myTemp3p1d(ii,ii,ii,ii), rhs(ii,ii,ii,ii)
         end do
         myTempCheckVar = maxval(abs(myTemp3p1d-rhs))
         print*, achar(9), achar(9), "Check L^infty error matrix_ijkl ", myTempCheckVar, " =0?"

         !if(rank==0) print*, achar(9), "hugh result testing done"
         
      end function solve_hugh_system


      !======================================================
      ! CALCULATE BANACH GRADIENT TERMS
      !  f = |v|^{s-2} ((s-2) e_v e_v + unitMatrix)
      !  g_ijkl = |nabla v|^{s-2} ((s-2) e_{partial_i v_k} e_{partial_j v_l} + delta_ij delta_lk )
      ! or 
      !  f_ij = |v|^{s-4} ((s-2) v_i v_j + |v|^2 delta_ij)
      !  g_ijkl = |nabla v|^{s-4} ((s-2) partial_i v_k partial_j v_l + |nabla v|^2 delta_ij delta_lk )
      ! rhs = s|v|^{s-2} v - s partial_j (|nabla v|^{s-2} partial_j v) - 1/lambda nabla rho + 1/lambda L2Grad
      !======================================================
      subroutine BanachGradient_calcFmatrixGtensorRhs(l2Grad, v_old, lambda, rho, matrix_f, tensor_g, rhs, tau)
         USE global_variables
         USE fftwfunction
         IMPLICIT NONE

         real(pr), dimension(1:n(1),1:n(2),1:local_N,1:3), intent(in) :: l2Grad, v_old
         real(pr), dimension(1:n(1),1:n(2),1:local_N), intent(in) :: rho ! lagrange multiplier for divergence free
         real(pr), intent(in) :: lambda   ! lagrange multiplier for norm
         real(pr), intent(in) :: tau      ! step size
         real(pr) :: s                    ! s = 3q/(q+1)

         real(pr), dimension(1:n(1),1:n(2),1:local_N,1:3), intent(out) :: rhs
         real(pr), dimension(1:n(1),1:n(2),1:local_N,1:3,1:3), intent(out) :: matrix_f
         real(pr), dimension(1:n(1),1:n(2),1:local_N,1:3,1:3,1:3,1:3), intent(out) :: tensor_g

         real(pr), dimension(:,:,:,:,:), allocatable :: gradv_old          !gradv_old(:,:,:,ii,jj) = partial_j v_i

         real(pr), dimension(:,:,:), allocatable :: normv_old, normv_oldPow, normgradv_old, normgradv_oldPow

         real(pr), dimension(:,:,:), allocatable :: tempSca
         real(pr), dimension(:,:,:,:), allocatable :: tempVec
         real(pr), dimension(:,:,:,:,:), allocatable :: tempTensor
         
         real(pr), dimension(:,:,:,:), allocatable :: e_v_old  ! unit vector in direction banachGrad_old
         real(pr), dimension(:,:,:,:,:), allocatable :: e_gradv_old  ! unit (norm) matrix in direction nabla (banachGrad_old) with convention (:,:,:,i,j) = partial_j v_i


         integer :: ii, jj, kk, ll


         s = 3.0_pr*lebesgueQ/(lebesgueQ+1.0_pr)

         allocate(gradv_old(1:n(1),1:n(2),1:local_N,1:3,1:3))

         
         allocate(tempVec(1:n(1),1:n(2),1:local_N,1:3))
         do ii=1,3
            call gradient(v_old(:,:,:,ii),tempVec)
            gradv_old(:,:,:,ii,:) = tempVec(:,:,:,:)
         end do
         deallocate(tempVec)

         allocate(normv_old(1:n(1),1:n(2),1:local_N))
         allocate(normv_oldPow(1:n(1),1:n(2),1:local_N))
         allocate(normgradv_old(1:n(1),1:n(2),1:local_N))
         allocate(normgradv_oldPow(1:n(1),1:n(2),1:local_N))

         
         normv_old(:,:,:) = calc_pointwiseVectorNorm(v_old)
         normgradv_old(:,:,:) = calc_pointwiseFrobeniusNorm(gradv_old)

         

         !!! calc tensor_g and matrix_f !!!

         if(use_e_u_instead_of_uqMinus4) then
            allocate( e_v_old(1:n(1),1:n(2),1:local_N,1:3) )
            e_v_old = calc_unitVectorInUdirection(v_old)
            normv_oldPow = calc_gk_order2(normv_old, s-2.0_pr)


            ! f = |v|^{s-2} (s-2) e_v e_v + ...(later)...
            if(toDealias .and. dealiase_if_mult_by_e_u) then
               do jj=1,3
                  do ii=1,3
                     matrix_f(:,:,:,ii,jj) = (s-2.0_pr) * e_v_old(:,:,:,ii) * e_v_old(:,:,:,jj)
                     call dealias_scalar(matrix_f(:,:,:,ii,jj),2.0_pr)
                     matrix_f(:,:,:,ii,jj) = matrix_f(:,:,:,ii,jj)*normv_oldPow(:,:,:)
                     call dealias_scalar(matrix_f(:,:,:,ii,jj),2.0_pr)
                  end do
               end do
            else
               do jj=1,3
                  do ii=1,3
                     matrix_f(:,:,:,ii,jj) = normv_oldPow(:,:,:)*(s-2.0_pr) * e_v_old(:,:,:,ii) * e_v_old(:,:,:,jj)
                  end do
               end do
            end if

            ! f = |v|^{s-2} ((s-2) e_v e_v + unitMatrix)
            do ii=1,3
               matrix_f(:,:,:,ii,ii) = matrix_f(:,:,:,ii,ii) + normv_oldPow(:,:,:)
            end do
            deallocate(e_v_old)





            ! g_ijkl = |nabla v|^{s-2} ((s-2) e_{partial_i v_k} e_{partial_j v_l} + delta_ij delta_lk )


            allocate( e_gradv_old(1:n(1),1:n(2),1:local_N,1:3,1:3) )
            e_gradv_old = calc_unitMatrixInMdirection(gradv_old)
            normgradv_oldPow = calc_gk_order2(normgradv_old, s-2.0_pr)

            if(toDealias .and. dealiase_if_mult_by_e_u) then
               do ll=1,3
                  do kk=1,3
                     do jj=1,3
                        do ii=1,3
                           ! g_ijkl = |nabla v|^{s-2} (s-2) e_{partial_i v_k} e_{partial_j v_l} + ...(later)...
                           tensor_g(:,:,:,ii,jj,kk,ll) = (s-2.0_pr) * e_gradv_old(:,:,:,kk,ii) * e_gradv_old(:,:,:,ll,jj)
                           call dealias_scalar(tensor_g(:,:,:,ii,jj,kk,ll),2.0_pr)
                           tensor_g(:,:,:,ii,jj,kk,ll) = normgradv_oldPow(:,:,:) * tensor_g(:,:,:,ii,jj,kk,ll)
                           call dealias_scalar(tensor_g(:,:,:,ii,jj,kk,ll),2.0_pr)
                        end do
                     end do
                  end do
               end do
            else
               do ll=1,3
                  do kk=1,3
                     do jj=1,3
                        do ii=1,3
                           ! g_ijkl = |nabla v|^{s-2} (s-2) e_{partial_i v_k} e_{partial_j v_l} + ...(later)...
                           tensor_g(:,:,:,ii,jj,kk,ll) = normgradv_oldPow(:,:,:) * (s-2.0_pr) * e_gradv_old(:,:,:,kk,ii) * e_gradv_old(:,:,:,ll,jj)
                        end do
                     end do
                  end do
               end do
            end if

            do kk=1,3
               do ii=1,3
                  ! g_ijkl = |nabla v|^{s-2} ((s-2) e_{partial_i v_k} e_{partial_j v_l} + delta_ij delta_lk )
                  tensor_g(:,:,:,ii,ii,kk,kk) = tensor_g(:,:,:,ii,ii,kk,kk) + normgradv_oldPow(:,:,:)
               end do
            end do

            deallocate(e_gradv_old)

         else

            normv_oldPow = calc_gk_order2(normv_old, s-4.0_pr)

            ! f = |v|^{s-4} (s-2) v v + ...(later)...
            if(toDealias) then
               do jj=1,3
                  do ii=1,3
                     matrix_f(:,:,:,ii,jj) = (s-2.0_pr) * v_old(:,:,:,ii) * v_old(:,:,:,jj)
                     call dealias_scalar(matrix_f(:,:,:,ii,jj),2.0_pr)
                     matrix_f(:,:,:,ii,jj) = normv_oldPow(:,:,:) * matrix_f(:,:,:,ii,jj)
                     call dealias_scalar(matrix_f(:,:,:,ii,jj),2.0_pr)
                  end do
               end do
            else
               do jj=1,3
                  do ii=1,3
                     matrix_f(:,:,:,ii,jj) = normv_oldPow(:,:,:) * (s-2.0_pr) * v_old(:,:,:,ii) * v_old(:,:,:,jj)
                  end do
               end do
            end if

            ! f_ij = |v|^{s-4} ((s-2) v_i v_j + |v|^2 delta_ij)
            !      = |v|^{s-4} (s-2) v_i v_j + |v|^{s-2} delta_ij
            normv_oldPow = calc_gk_order2(normv_old, s-2.0_pr)
            
            do ii=1,3
               matrix_f(:,:,:,ii,ii) = matrix_f(:,:,:,ii,ii) + normv_oldPow(:,:,:)
            end do

            ! g_ijkl = |nabla v|^{s-4} ((s-2) partial_i v_k partial_j v_l + |nabla v|^2 delta_ij delta_lk )
            normgradv_oldPow = calc_gk_order2(normgradv_old, s-4.0_pr)

            if(toDealias) then
               do ll=1,3
                  do kk=1,3
                     do jj=1,3
                        do ii=1,3
                           ! g_ijkl = |nabla v|^{s-4} (s-2) partial_i v_k partial_j v_l + ...(later)...
                           tensor_g(:,:,:,ii,jj,kk,ll) = (s-2.0_pr) * gradv_old(:,:,:,kk,ii) * gradv_old(:,:,:,ll,jj)
                           call dealias_scalar(tensor_g(:,:,:,ii,jj,kk,ll), 2.0_pr)
                           tensor_g(:,:,:,ii,jj,kk,ll) = normgradv_oldPow(:,:,:) * tensor_g(:,:,:,ii,jj,kk,ll)
                           call dealias_scalar(tensor_g(:,:,:,ii,jj,kk,ll), 2.0_pr)                           
                        end do
                     end do
                  end do
               end do
            else
               do ll=1,3
                  do kk=1,3
                     do jj=1,3
                        do ii=1,3
                           ! g_ijkl = |nabla v|^{s-4} (s-2) partial_i v_k partial_j v_l + ...(later)...
                           tensor_g(:,:,:,ii,jj,kk,ll) = normgradv_oldPow(:,:,:) * (s-2.0_pr) * gradv_old(:,:,:,kk,ii) * gradv_old(:,:,:,ll,jj)
                        end do
                     end do
                  end do
               end do
            end if

            normgradv_oldPow = calc_gk_order2(normgradv_old, s-2.0_pr)
            do kk=1,3
               do ii=1,3
                  ! g_ijkl = |nabla v|^{s-4} ((s-2) partial_i v_k partial_j v_l + |nabla v|^2 delta_ij delta_lk )
                  tensor_g(:,:,:,ii,ii,kk,kk) = tensor_g(:,:,:,ii,ii,kk,kk) + normgradv_oldPow(:,:,:)
               end do
            end do
         end if

         !!! tensor_g and matrix_f done !!!

         
         !!! calculating rhs !!!
         ! rhs = (s-1-tau)|v|^{s-2} v - (s-1-tau) partial_j (|nabla v|^{s-2} partial_j v) - tau/lambda nabla rho + tau/lambda L2Grad
         normv_oldPow = calc_gk_order2(normv_old, s-2.0_pr)
         normgradv_oldPow = calc_gk_order2(normgradv_old, s-2.0_pr)
         allocate(tempTensor(1:n(1),1:n(2),1:local_N,1:3,1:3))
         
         do ii=1,3
            rhs(:,:,:,ii) = (s-1.0_pr-tau) * normv_oldPow(:,:,:) * v_old(:,:,:,ii)
            if (toDealias) call dealias_scalar(rhs(:,:,:,ii), 2.0_pr)
         end do

         do jj=1,3
            do ii=1,3
               ! tempTensor(i,j) = |nabla v|^{s-2} partial_j v_i
               tempTensor(:,:,:,ii,jj) = normgradv_oldPow(:,:,:)*gradv_old(:,:,:,ii,jj)    ! gradv(:,:,:,i,j) = partial_j v_i
               if (toDealias) call dealias_scalar(tempTensor(:,:,:,ii,jj), 2.0_pr)
            end do
         end do

         allocate(tempSca(1:n(1),1:n(2),1:local_N))
         allocate(tempVec(1:n(1),1:n(2),1:local_N,1:3))
         do kk=1,3
            tempVec(:,:,:,:) = tempTensor(:,:,:,kk,:)   ! tempTensor(:,:,:,i,j) = ... * partial_j v_i
            call divergence(tempVec, tempSca)
            rhs(:,:,:,kk) = rhs(:,:,:,kk) - (s-1.0_pr-tau) * tempSca(:,:,:)
         end do

         deallocate(tempSca)
         deallocate(tempTensor)

         call gradient(rho,tempVec)
         rhs(:,:,:,:) = rhs(:,:,:,:) - tau/lambda * tempVec(:,:,:,:) + tau/lambda*l2Grad(:,:,:,:)

         deallocate(tempVec)
         !!! calculating f, g, rhs done !!!

      end subroutine BanachGradient_calcFmatrixGtensorRhs


      !======================================================
      ! CALCULATE THE RHO IN THE BANACH GRADIENT
      ! given the things solve rho = lambda nabla cdot (|v|^{s-2} v - partial_j |nabla v|^{s-2} partial_j v)
      !======================================================
      function BanachGradientCalcRho(v, lambda) result (rho)
         USE global_variables
         USE fftwfunction
         IMPLICIT NONE

         real(pr), dimension(1:n(1),1:n(2),1:local_N,1:3), intent(in) :: v
         real(pr), intent(in) :: lambda   ! lagrange multiplier for norm
         real(pr) :: s                    ! s = 3q/(q+1)
         real(pr), dimension(1:n(1),1:n(2),1:local_N) :: rho

         real(pr), dimension(1:n(1),1:n(2),1:local_N,1:3) :: poisson_rhs_vec
         real(pr), dimension(1:n(1),1:n(2),1:local_N) :: poisson_rhs

         real(pr), dimension(1:n(1),1:n(2),1:local_N) :: temp_sca
         real(pr), dimension(1:n(1),1:n(2),1:local_N,1:3) :: temp_vec
         real(pr), dimension(1:n(1),1:n(2),1:local_N,1:3,1:3) :: temp_matrix
         real(pr), dimension(1:n(1),1:n(2),1:local_N,1:3,1:3) :: gradv

         real(pr), dimension(1:n(1),1:n(2),1:local_N) :: v_norm, gradv_norm, v_normPow, gradv_normPow

         integer :: ii, jj


         s = 3.0_pr*lebesgueQ/(lebesgueQ+1.0_pr)

         do ii=1,3
            call gradient(v(:,:,:,ii),temp_vec)
            do jj=1,3
               gradv(:,:,:,ii,jj) = temp_vec(:,:,:,jj)         ! grad v(:,:,:,ii,jj) = partial_j v_i
            end do
         end do

         v_norm(:,:,:) = calc_pointwiseVectorNorm(v)
         gradv_norm(:,:,:) = calc_pointwiseFrobeniusNorm(gradv)

         v_normPow = calc_gk_order2(v_norm, s-2.0_pr)
         gradv_normPow = calc_gk_order2(gradv_norm, s-2.0_pr)

         ! poisson_rhs_vec = |v|^{s-2} v - ...(later)...
         do ii=1,3
            poisson_rhs_vec(:,:,:,ii) = v_normPow(:,:,:)*v(:,:,:,ii)
            if(toDealias) call dealias_scalar(poisson_rhs_vec(:,:,:,ii), 2.0_pr)
         end do

         !temp_matrix(:,:,:,i,j) = |nabla v|^{s-2} partial_j v_i
         do jj=1,3
            do ii=1,3
               temp_matrix(:,:,:,ii,jj) = gradv_normPow(:,:,:)*gradv(:,:,:,ii,jj)
               if(toDealias) call dealias_scalar(temp_matrix(:,:,:,ii,jj), 2.0_pr)
            end do
         end do


         ! poisson_rhs_vec = lambda ( |v|^{s-2} v - partial_k (|nabla v|^{s-2} partial_k v))
         do ii=1,3
            ! temp_matrix(:,:,:,ii,jj) = ...*partial_j v_i
            ! temp_sca(:,:,:,jj) = partial_j ( ...* partial_j v_i )
            temp_vec(:,:,:,:) = temp_matrix(:,:,:,ii,:)
            call divergence(temp_vec, temp_sca)
            poisson_rhs_vec(:,:,:,ii) = lambda*poisson_rhs_vec(:,:,:,ii) - lambda*temp_sca(:,:,:)            
         end do

         call divergence(poisson_rhs_vec, poisson_rhs)

         call solve_poisson(poisson_rhs, 1.0_pr, rho)

      end function BanachGradientCalcRho



      !=======================================================
      ! CALCULATE THE GLOBAL SUMMED INNER PRODUCT BETWEEN TWO (vector) FIELDS
      !=======================================================
      FUNCTION matrixInversion3by3(matrix) RESULT (result)
         USE global_variables
         use mpi
         IMPLICIT NONE
   
         complex(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3,1:3), intent(in) :: matrix
         complex(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3,1:3) :: result, test
         complex(pr), DIMENSION(1:n(1),1:n(2),1:local_N) :: a,b,c,d,e,f,g,h,i, fac
         integer :: ii, jj, kk

         a(:,:,:) = matrix(:,:,:,1,1)
         b(:,:,:) = matrix(:,:,:,1,2)
         c(:,:,:) = matrix(:,:,:,1,3)
         d(:,:,:) = matrix(:,:,:,2,1)
         e(:,:,:) = matrix(:,:,:,2,2)
         f(:,:,:) = matrix(:,:,:,2,3)
         g(:,:,:) = matrix(:,:,:,3,1)
         h(:,:,:) = matrix(:,:,:,3,2)
         i(:,:,:) = matrix(:,:,:,3,3)

         !fac = (a(ei-fh)+b(fg-di)+c(dh-eg))^{-1}
         fac(:,:,:) = 1.0_pr/(a(:,:,:)*(e(:,:,:)*i(:,:,:)-f(:,:,:)*h(:,:,:))+b(:,:,:)*(f(:,:,:)*g(:,:,:)-d(:,:,:)*i(:,:,:))+c(:,:,:)*(d(:,:,:)*h(:,:,:)-e(:,:,:)*g(:,:,:)))

         ! ei-fh  ch-bi  bf-ce
         ! fg-di  ai-cg  cd-af
         ! dh-eg  bg-ah  ae-bd

         result(:,:,:,1,1) = e(:,:,:)*i(:,:,:)-f(:,:,:)*h(:,:,:)
         result(:,:,:,1,2) = c(:,:,:)*h(:,:,:)-b(:,:,:)*i(:,:,:)
         result(:,:,:,1,3) = b(:,:,:)*f(:,:,:)-c(:,:,:)*e(:,:,:)
         result(:,:,:,2,1) = f(:,:,:)*g(:,:,:)-d(:,:,:)*i(:,:,:)
         result(:,:,:,2,2) = a(:,:,:)*i(:,:,:)-c(:,:,:)*g(:,:,:)
         result(:,:,:,2,3) = c(:,:,:)*d(:,:,:)-a(:,:,:)*f(:,:,:)
         result(:,:,:,3,1) = d(:,:,:)*h(:,:,:)-e(:,:,:)*g(:,:,:)
         result(:,:,:,3,2) = b(:,:,:)*g(:,:,:)-a(:,:,:)*h(:,:,:)
         result(:,:,:,3,3) = a(:,:,:)*e(:,:,:)-b(:,:,:)*d(:,:,:)

         do jj=1,3
            do ii=1,3
               result(:,:,:,ii,jj) = fac(:,:,:)*result(:,:,:,ii,jj)
            end do
         end do

         test = 0.0_pr
         do kk=1,3
            do jj=1,3
               do ii=1,3
                  test(:,:,:,ii,kk) = test(:,:,:,ii,kk) + result(:,:,:,ii,jj)*matrix(:,:,:,jj,kk)
               end do
            end do
         end do

         if(.false. .and. rank == 0) then
            print*, "matrix inversion test start"
   
            print*, "a"
            print*, test(1,1,1,1,:)
            print*, test(1,1,1,2,:)
            print*, test(1,1,1,3,:)
   
            print*, "b"
            print*, test(2,2,1,1,:)
            print*, test(2,2,1,2,:)
            print*, test(2,2,1,3,:)
   
            print*, "c"
            print*, test(3,2,4,1,:)
            print*, test(3,2,4,2,:)
            print*, test(3,2,4,3,:)
            print*, "matrix inversion test end"
         end if
      END FUNCTION matrixInversion3by3 


      subroutine set_q_resol_bIterOffset_optimIterOffsets(fileName)
         use global_variables
         implicit none
         character(len=*), intent(in) :: fileName
         character(len=:), allocatable :: tempString1, tempString2
         integer :: startIndex, endIndex

         ! extract lebesgueQ value
         startIndex = index(fileName,'_q')   ! first occurance
         tempString1 = fileName(startIndex+2:)
         endIndex = index(tempString1,"_")
         tempString2 = tempString1(:endIndex-1)
         read(tempString2,*) lebesgueQ

         ! extract resol value
         startIndex = index(fileName,'_n')   ! first occurance
         tempString1 = fileName(startIndex+2:)
         endIndex = index(tempString1,"_")
         tempString2 = tempString1(:endIndex-1)
         read(tempString2,*) resol

         ! extract bIterOffset value
         startIndex = index(fileName,'_B')   ! first occurance
         tempString1 = fileName(startIndex+2:)
         endIndex = index(tempString1,"_")
         tempString2 = tempString1(:endIndex-1)
         read(tempString2,*) bIterOffset
         !!! might be reduced by 1 when optimizationIterOffset if not iterend !!!


         ! extract optimizationIterOffset value
         startIndex = index(fileName,'_iter')   ! first occurance
         tempString1 = fileName(startIndex+5:)
         endIndex = index(tempString1,".nc")
         tempString2 = tempString1(:endIndex-1)
         if(index(tempString2,"end")>0) then
            optimizationIterOffset = 0
         else
            bIterOffset = bIterOffset - 1
            read(tempString2,*) optimizationIterOffset
         end if

         if(rank==0) print*, "assigned lebesgueQ = ", lebesgueQ
         if(rank==0) print*, "assigned resol = ", resol
         if(rank==0) print*, "assigned bIterOffset = ", bIterOffset
         if(rank==0) print*, "assigned optimizationIterOffset = ", optimizationIterOffset
         
      end subroutine set_q_resol_bIterOffset_optimIterOffsets


END MODULE
