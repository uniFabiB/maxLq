MODULE fftwfunction   ! Newly added on March 20, 2017

      use, intrinsic :: iso_c_binding    
      use global_variables
      use mpi
      implicit none
      include 'fftw3-mpi.f03'

      type(c_ptr) :: fwdplan, bwdplan, tmpdata                    ! complex fft
      complex(c_double_complex), pointer :: tmppointer(:,:,:)     ! complex fft

      type(c_ptr) :: fwdplan_r, bwdplan_r, tmpdata_r, tmpdata_rcx ! real fft
      complex(c_double_complex), pointer :: tmppointer_rcx(:,:,:) ! real fft
      real(c_double), pointer :: tmppointer_r(:,:,:)              ! real fft
      



      CONTAINS

      !=======================================
      ! Forward Fourier transform (complex, vector)
      !=======================================      
      SUBROUTINE fftfwdV(vec,fvec)
            implicit none

            complex(pr), dimension (1:n(1),1:n(2),1:local_N,1:3), intent(in) :: vec
            complex(pr), dimension (1:n(1),1:n(2),1:local_N,1:3), intent(out) :: fvec
            complex(pr), dimension (1:n(1),1:n(2),1:local_N) :: aux, faux
            integer :: ii
            
            do ii=1,3
                  aux(:,:,:) = vec(:,:,:,ii)
                  call fftfwd(aux, faux)
                  fvec(:,:,:,ii) = faux(:,:,:)
            end do

      END SUBROUTINE fftfwdV

      !=======================================
      ! Backward Fourier transform (complex, vector)
      !=======================================      
      SUBROUTINE fftbwdV(fvec, vec)
            implicit none

            complex(pr), dimension (1:n(1),1:n(2),1:local_N,1:3), intent(in) :: fvec
            complex(pr), dimension (1:n(1),1:n(2),1:local_N,1:3), intent(out) :: vec
            complex(pr), dimension (1:n(1),1:n(2),1:local_N) :: aux, faux
            integer :: ii
            
            do ii=1,3
                  faux(:,:,:) = fvec(:,:,:,ii)
                  call fftbwd(faux, aux)
                  vec(:,:,:,ii) = aux(:,:,:)
            end do

      END SUBROUTINE fftbwdV

      !=======================================
      ! Forward Fourier transform (complex)
      !=======================================      
      SUBROUTINE fftfwd(u,fu)
            implicit none 

            complex(pr), dimension (1:n(1),1:n(2),1:local_N), intent(in) :: u
            complex(pr), dimension (1:n(1),1:n(2),1:local_N), intent(out) :: fu

            tmppointer(:,:,:) = u(:,:,:)

            call fftw_mpi_execute_dft(fwdplan, tmppointer, tmppointer)

            fu(:,:,:) = tmppointer(:,:,:)



            if (isnan(real(fu(1,1,1)))) then
                  print*, "nan in fftwfunction"
            end if

            call mpi_barrier(mpi_comm_world, statinfo)

      END SUBROUTINE fftfwd

      !=======================================
      ! Inverse Fourier transform (complex)
      !=======================================
      SUBROUTINE fftbwd(fu,u)
            implicit none

            complex(pr), dimension (1:n(1),1:n(2),1:local_N), intent(in) :: fu
            complex(pr), dimension (1:n(1),1:n(2),1:local_N), intent(out) :: u
            
            tmppointer(:,:,:) = fu(:,:,:)

            call fftw_mpi_execute_dft(bwdplan, tmppointer, tmppointer)

            u(:,:,:) = tmppointer(:,:,:)
            u(:,:,:) = u(:,:,:)/(product(real(n,pr)))

            call mpi_barrier(mpi_comm_world, statinfo)

      END SUBROUTINE fftbwd

      !=======================================
      ! Forward Fourier transform (real)
      !=======================================      
      SUBROUTINE fftfwdr(u,fu)
            implicit none 

            real(pr), dimension (1:n(1),1:n(2),1:local_N), intent(in) :: u
            complex(pr), dimension (1:n(1)/2+1,1:n(2),1:local_N), intent(out) :: fu

            tmppointer_r(1:n(1),:,:) = u(:,:,:)
            tmppointer_r(n(1)+1:n(1)+2,:,:) = 0.0_pr

            call fftw_mpi_execute_dft_r2c(fwdplan_r, tmppointer_r, tmppointer_rcx)

            fu(:,:,:) = tmppointer_rcx(:,:,:)

            call mpi_barrier(mpi_comm_world, statinfo)

      END SUBROUTINE fftfwdr

      !=======================================
      ! Inverse Fourier transform (real)
      !=======================================
      SUBROUTINE fftbwdr(fu,u)
            implicit none

            complex(pr), dimension (1:n(1)/2+1,1:n(2),1:local_N), intent(in) :: fu
            real(pr), dimension (1:n(1),1:n(2),1:local_N), intent(out) :: u
            !integer :: i,j,k
            
            tmppointer_rcx(:,:,:) = fu(:,:,:)

            call fftw_mpi_execute_dft_c2r(bwdplan_r, tmppointer_rcx, tmppointer_r)




            !do k = 1, local_N
            !      do j = 1, n(2)
            !            do i = 1, n(1)
            !                  u(i,j,k) = tmppointer_r(i,j,k)
            !            end do
            !      end do
            !end do
            u(1:n(1),:,:) = tmppointer_r(1:n(1),:,:)
            u = u/(product(real(n,pr)))

            call mpi_barrier(mpi_comm_world, statinfo)

      END SUBROUTINE fftbwdr

      !=======================================
      ! innitialize
      !=======================================
      SUBROUTINE init_fft()
            implicit none            

            ! complex fft
            tmpdata = fftw_alloc_complex(C_local_alloc)
            call c_f_pointer(tmpdata, tmppointer, [C_n(1),C_n(2),C_local_N])
            fwdplan = fftw_mpi_plan_dft_3d(C_n(3), C_n(2), C_n(1), tmppointer, tmppointer, mpi_comm_world, fftw_forward, fftw_measure)
            bwdplan = fftw_mpi_plan_dft_3d(C_n(3), C_n(2), C_n(1), tmppointer, tmppointer, mpi_comm_world, fftw_backward, fftw_measure)

            ! real fft
            tmpdata_r = fftw_alloc_real(2*C_local_alloc)
            tmpdata_rcx = fftw_alloc_complex(C_local_alloc)
            call c_f_pointer(tmpdata_r, tmppointer_r, [C_n(1)+2, C_n(2),C_local_N])
            call c_f_pointer(tmpdata_rcx, tmppointer_rcx, [C_n(1)/2+1, C_n(2),C_local_N])
            fwdplan_r = fftw_mpi_plan_dft_r2c_3d(C_n(3), C_n(2), C_n(1), tmppointer_r, tmppointer_rcx, mpi_comm_world, fftw_measure)
            bwdplan_r = fftw_mpi_plan_dft_c2r_3d(C_n(3), C_n(2), C_n(1), tmppointer_rcx, tmppointer_r, mpi_comm_world, fftw_measure)

      END SUBROUTINE init_fft


      !=======================================
      ! deallocate
      !=======================================
      SUBROUTINE fft_deallocate()
            implicit none            

            ! complex fft
            call fftw_destroy_plan(fwdplan)
            call fftw_destroy_plan(bwdplan)
            call fftw_free(tmpdata)

            ! real fft
            call fftw_destroy_plan(fwdplan_r)
            call fftw_destroy_plan(bwdplan_r)
            call fftw_free(tmpdata_r)
            call fftw_free(tmpdata_rcx)

            call fftw_mpi_cleanup()
            
      END SUBROUTINE fft_deallocate



END MODULE