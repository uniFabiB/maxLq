MODULE data_ops

  IMPLICIT NONE
 
  CONTAINS
        !==============================
        ! SAVE CONTROL VARIABLE
        !==============================
        SUBROUTINE save_field(field, fieldName)
         USE global_variables
         IMPLICIT NONE

         REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(IN) :: field
         CHARACTER(len=*), INTENT(IN) :: fieldName

         REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: fx, fy, fz

         character(10) :: tempResolTxt
         character(len=:), allocatable :: resolTxt
         character(len=:), allocatable :: filename
   
         ALLOCATE( fx(1:n(1),1:n(2),1:local_N) )
         ALLOCATE( fy(1:n(1),1:n(2),1:local_N) )
         ALLOCATE( fz(1:n(1),1:n(2),1:local_N) )

         WRITE(tempResolTxt, '(i10)') resol

         resolTxt = trim(adjustl(tempResolTxt))
         

         filename = ncDir//fieldName//"_q"//lebesgueQTxt//"_n"//resolTxt//"_B"//bIterTxt//"_iter"//trim(optimizationIterationTxt)//".nc"   ! Newly added on May 8, 2017

         fx = field(:,:,:,1)
         fy = field(:,:,:,2)
         fz = field(:,:,:,3)
         CALL save_field_R3toR3_ncdf(fx,fy,fz,"Ux", "Uy", "Uz", filename, "netCDF")

         DEALLOCATE( fx )
         DEALLOCATE( fy )
         DEALLOCATE( fz )


        END SUBROUTINE save_field

        !============================
        !    SAVE CONTROL VELOCITY
        !============================
        SUBROUTINE save_gradient(myfield, myindex, mysystem)
          USE global_variables  
          IMPLICIT NONE

          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(IN) :: myfield
          INTEGER, INTENT(IN) :: myindex
          CHARACTER(len=*), INTENT(IN) :: mysystem

          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N) :: fx, fy, fz
  
          CHARACTER(2) :: K0txt
          CHARACTER(2) :: E0txt
          CHARACTER(2) :: IGtxt
          CHARACTER(4) :: itertxt
!          CHARACTER(2) :: WEIGHTtxt
          CHARACTER(200) :: filename
      
          WRITE(K0txt, '(i2.2)') K0_index
          WRITE(E0txt, '(i2.2)') E0_index
          WRITE(IGtxt, '(i2.2)') iguess
          WRITE(itertxt, '(i4.4)') myindex
!          WRITE(WEIGHTtxt, '(i2.2)') int_WEIGHT

          SELECT CASE (mysystem)
            CASE ("maxdLqdt") 
   !              filename = "/work/yund0050/MultiObjective_095_01/WEIGHT"//WEIGHTtxt//"_E"//E0txt//"_maxdEdt_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_gradJ.nc"
               filename = ConstraintDir//"gradJ"//"_B"//bIterTxt//".nc"   ! Newly added on May 8, 2017
               fx = myfield(:,:,:,1)
               fy = myfield(:,:,:,2)
               fz = myfield(:,:,:,3)

               CALL save_field_R3toR3_ncdf(fx,fy,fz,"gradJ_x", "gradJ_y", "gradJ_z", filename, "netCDF")

            CASE ("maxdEdt") 
!              filename = "/work/yund0050/MultiObjective_095_01/WEIGHT"//WEIGHTtxt//"_E"//E0txt//"_maxdEdt_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_gradJ.nc"
              filename = HomeDir//"_E"//E0txt//"_IG"//IGtxt//"_gradJ.nc"   ! Newly added on May 8, 2017
              fx = myfield(:,:,:,1)
              fy = myfield(:,:,:,2)
              fz = myfield(:,:,:,3)

              CALL save_field_R3toR3_ncdf(fx,fy,fz,"gradJ_x", "gradJ_y", "gradJ_z", filename, "netCDF")

            CASE ("FixK0E0")
!              filename = "/work/yund0050/MultiObjective_095_01/WEIGHT"//WEIGHTtxt//"_E"//E0txt//"_FixK0E0_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_gradJ.nc"
              filename = HomeDir//"_E"//E0txt//"_IG"//IGtxt//"_gradJ.nc"   ! Newly added on May 8, 2017
              fx = myfield(:,:,:,1)
              fy = myfield(:,:,:,2)
              fz = myfield(:,:,:,3)

              CALL save_field_R3toR3_ncdf(fx,fy,fz,"gradJ_x", "gradJ_y", "gradJ_z", filename, "netCDF")

          END SELECT 

        END SUBROUTINE save_gradient

        !============================
        ! SAVE DIAGNOSTIC SCALARS 
        !============================
        SUBROUTINE save_diagnosticScalars(myField, numFields, myFieldNames)
          USE global_variables  
          IMPLICIT NONE

          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:numFields), INTENT(IN) :: myField
          INTEGER, INTENT(IN) :: numFields
          CHARACTER(len=*), INTENT(IN) :: myFieldNames

          CHARACTER(2) :: K0txt
          CHARACTER(2) :: E0txt
          CHARACTER(2) :: IGtxt
!          CHARACTER(2) :: WEIGHTtxt
          CHARACTER(200) :: filename

          WRITE(K0txt, '(i2.2)') K0_index
          WRITE(E0txt, '(i2.2)') E0_index
          WRITE(IGtxt, '(i2.2)') iguess
!          WRITE(WEIGHTtxt, '(i2.2)') int_WEIGHT
              
!          filename = "/work/yund0050/MultiObjective_095_01/WEIGHT"//WEIGHTtxt//"_E"//E0txt//"_maxdEdt_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_diagScalar.nc"
          filename = ncDir//"diagScalar"//"-B"//bIterTxt//"-iter-"//trim(optimizationIterationTxt)//".nc"   ! Newly added on May 8, 2017

          CALL save_field_R3toRn_ncdf(myfield, numFields, myFieldNames, filename)


        END SUBROUTINE save_diagnosticScalars

        !====================================
        ! SAVE GLOBAL DIAGNOSTICS OF FIELDS 
        !====================================
        SUBROUTINE save_diagnosticFields_global(mysystem, K, E, Umax, Wmax, magUmax, magWmax, &
                                                H, maxHel, minHel, vorCoreData)
          USE global_variables  
          IMPLICIT NONE

          CHARACTER(len=*), INTENT(IN) :: mysystem
          REAL(pr), DIMENSION(1:3), INTENT(IN) :: K, E, Umax, Wmax, vorCoreData
          REAL(pr), INTENT(IN) :: magUmax, magWmax, H, maxHel, minHel

          CHARACTER(2) :: K0txt
          CHARACTER(2) :: E0txt
          CHARACTER(2) :: IGtxt
!          CHARACTER(2) :: WEIGHTtxt
          CHARACTER(200) :: filename
      
          WRITE(K0txt, '(i2.2)') K0_index
          WRITE(E0txt, '(i2.2)') E0_index
          WRITE(IGtxt, '(i2.2)') iguess
!          WRITE(WEIGHTtxt, '(i2.2)') int_WEIGHT
          
!          filename = "/work/yund0050/MultiObjective_095_01/WEIGHT"//WEIGHTtxt//"_E"//E0txt//"_maxdEdt_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_diagFields.dat"
          filename = ConstraintDir//"diagFields"//"-q"//lebesgueQTxt//"-B"//bIterTxt//"-iter"//trim(optimizationIterationTxt)//".dat"   ! Newly added on May 8, 2017

          OPEN(10, FILE = filename, FORM = 'FORMATTED', STATUS = 'REPLACE')
          WRITE(10,*) "# K, E, Umax, Wmax, magUmax, magWmax, H, MaxminH, MaxminS "
          WRITE(10,*) "# x   y   z " 
          WRITE(10, "(3 ES20.12)") K(1), K(2), K(3) 
          WRITE(10, "(3 ES20.12)") E(1), E(2), E(3) 
          WRITE(10, "(3 ES20.12)") Umax(1), Umax(2), Umax(3) 
          WRITE(10, "(3 ES20.12)") Wmax(1), Wmax(2), Wmax(3)
          WRITE(10, "(3 ES20.12)") magUmax, magWmax, H 
          WRITE(10, "(2 ES20.12)") maxHel, minHel
          WRITE(10, "(2 ES20.12)") vorCoreData(1), vorCoreData(2)
          CLOSE(10)
 
        END SUBROUTINE save_diagnosticFields_global


         !============================================================
         !          SAVE USED PARAMETERS
         !============================================================
         subroutine saveUsedParams()
            use global_variables
            character(200) :: filename
            
            filename = HomeDir//"params"//".dat"

            if(rank==0) then
               open(10, file = filename, form = 'FORMATTED', status = 'REPLACE')


               !ES = nice exponential
               !A = char
               !I = integer
               !G = general
               write(10, "(A20, I40)") "resol", n(1)
               write(10, "(A20, ES40.12)") "lebesgueQ", lebesgueQ
               write(10, "(A20, G40.12)") "iguess", iguess
               write(10, "(A20, A40)") "loadTempFunctionName", loadTempFunctionName
               write(10, "(A20, I40)") "bIterOffset", bIterOffset
               write(10, "(A20, I40)") "optimizationIterOffset", optimizationIterOffset
               
               
               write(10, "(A20)") " "
               write(10, "(A20, G40.12)") "useOrthogonalGradient", useOrthogonalGradient
               write(10, "(A20, G40.12)") "useConjugateGradient", useConjugateGradient
               write(10, "(A20, G40.12)") "useRiemannianGeometry", useRiemannianGeometry
               write(10, "(A20)") " "
               write(10, "(A20, G40.12)") "normalizeGradient", normalizeGradient
               write(10, "(A20, G40.12)") "use_e_u_auto_for_q_less_4", use_e_u_auto_for_q_less_4
               write(10, "(A20, G40.12)") "use_e_u_instead_of_uqMinus4", use_e_u_instead_of_uqMinus4
               write(10, "(A20, G40.12)") "dealiase_if_mult_by_e_u", dealiase_if_mult_by_e_u
               write(10, "(A20)") " "
               write(10, "(A20, G40.12)") "MAX_ITER", MAX_ITER
               write(10, "(A20, ES40.12)") "OPTIM_TOL", OPTIM_TOL
               write(10, "(A20, I40)") "toleranceRollingAverageSize", toleranceRollingAverageSize
               write(10, "(A20, I40)") "resetMomentumTermEveryXiterations", resetMomentumTermEveryXiterations
               write(10, "(A20, ES40.12)") "MACH_EPSILON", MACH_EPSILON
               write(10, "(A20, ES40.12)") "TAU_MAX", TAU_MAX
               write(10, "(A20, ES40.12)") "lambda1", lambda1
               write(10, "(A20)") " "
               write(10, "(A20, ES40.12)") "visc", visc
               write(10, "(A20, G40.12)") "viscCoefficient", viscCoefficient
               write(10, "(A20, G40.12)") "pressureCoefficient", pressureCoefficient
               write(10, "(A20)") " "
               write(10, "(A20, I40)") "local_N", local_N
               write(10, "(A20, I40)") "total_local_size", total_local_size

               close(10)
            end if
         end subroutine
        !============================================================
        !          SAVE SPECTRAL DATA
        !============================================================
        SUBROUTINE save_spectral_data(mydata, name)
          USE global_variables
          IMPLICIT NONE
 
          REAL(pr), DIMENSION(1:n(1),1:2), INTENT(IN) :: mydata
          character(len=:), allocatable :: spectraDir
          CHARACTER(200) :: filename
          !character(10) :: dealiasing_str
          CHARACTER(*) :: name
          INTEGER :: i


          !if(toDealias) then
          !  dealiasing_str = "deal_"
          !else
          !  dealiasing_str = "noDeal_"
          !end if

          !filename = HomeDir//trim(dealiasing_str)//name//"_spectrum.dat"



          spectraDir = ConstraintDir//"spectra/"
          call createDirectoryIfNonExistent(spectraDir)

          filename = spectraDir//name//"-spectrum"//"-B"//bIterTxt//"-iter"//trim(optimizationIterationTxt)//".dat"

          if(rank==0) then
            OPEN(10, FILE = filename, FORM = 'FORMATTED', STATUS = 'REPLACE')
            DO i=1,n(1)
               WRITE(10, "(2 ES20.12)") mydata(i,1), mydata(i,2)
            END DO
            CLOSE(10)
          end if
        END SUBROUTINE save_spectral_data

        !==========================================
        ! SAVE 3D SCALAR IN netCDF FORMAT
        !==========================================
        SUBROUTINE save_field_R3toR1_ncdf(myfield, field_name, file_name)
        
          USE global_variables
          USE netcdf
          IMPLICIT NONE
          INCLUDE "mpif.h"
        

          REAL(pr), DIMENSION(:,:,:), INTENT(IN) :: myfield
          CHARACTER(len=*) :: field_name
          CHARACTER(len=*) :: file_name
         
          REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: global_field
          INTEGER, DIMENSION(1:3) :: starts, counts
          INTEGER :: ncout, ncid, varid, dimids(3)
          INTEGER :: x_dimid, y_dimid, z_dimid

          INTEGER :: ii
          
!          print *,  "starting save_field_R3toR1_ncdf to file ", file_name

          IF (parallel_data) THEN
             IF (rank==0) THEN
                ncout = nf90_create(file_name, NF90_CLOBBER, ncid=ncid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_def_dim(ncid, "x", n(1), x_dimid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_def_dim(ncid, "y", n(2), y_dimid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_def_dim(ncid, "z", NF90_UNLIMITED, z_dimid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                dimids = (/ x_dimid, y_dimid, z_dimid /)

                ncout = nf90_def_var(ncid, TRIM(field_name), NF90_DOUBLE, dimids, varid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

                ncout = nf90_enddef(ncid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_close(ncid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
             END IF
             CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

             starts = (/ 1, 1, rank*local_N+1 /)
             counts = (/ n(1), n(2), local_N /)            
 
             !!--------------------------
             !! START netCDF ROUTINES
             !!--------------------------
             DO ii=0,np-1
                IF (rank==ii) THEN 
                   ncout = nf90_open(file_name, NF90_WRITE, ncid)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

                   ncout = nf90_inq_varid(ncid, TRIM(field_name), varid)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                   ncout = nf90_put_var(ncid, varid, myfield, start = starts, count = counts)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                   
                   ncout = nf90_close(ncid)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                END IF
                CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo) 
             END DO
             CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
          ELSE
             IF (rank == 0) THEN 
                ALLOCATE( global_field(1:n(1),1:n(2),1:n(3)) )
             END IF
             CALL MPI_GATHER(myfield, total_local_size, MPI_REAL8, global_field, total_local_size, MPI_REAL8, 0, MPI_COMM_WORLD, Statinfo)

             IF (rank==0) THEN
                ncout = nf90_create(file_name, NF90_CLOBBER, ncid)
                ncout = nf90_def_dim(ncid, "x", n(1), x_dimid)
                ncout = nf90_def_dim(ncid, "y", n(2), y_dimid)
                ncout = nf90_def_dim(ncid, "z", n(3), z_dimid)
          
                dimids =  (/ x_dimid, y_dimid, z_dimid /)

                ncout = nf90_def_var(ncid, field_name, NF90_DOUBLE, dimids, varid)
                ncout = nf90_enddef(ncid)
          
                ncout = nf90_put_var(ncid, varid, global_field)
                ncout = nf90_close(ncid)
             
                DEALLOCATE( global_field ) 
             END IF  
             CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

          END IF

        END SUBROUTINE save_field_R3toR1_ncdf

        !==========================================
        ! SAVE FIELD IN R3
        !==========================================
        SUBROUTINE save_field_R3toR3_ncdf(f1, f2, f3, f1_name, f2_name, f3_name, file_name, myformat)
          
          use, intrinsic :: iso_c_binding
          USE global_variables
          USE netcdf
          IMPLICIT NONE
          INCLUDE "mpif.h"

          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N), INTENT(IN) :: f1, f2, f3
          CHARACTER(len=*) :: file_name
          CHARACTER(len=*) :: f1_name, f2_name, f3_name
          CHARACTER(len=*) :: myformat
         
          REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: u1, u2, u3
          REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: myfield

          REAL(pr) :: local_maxf1, local_maxf2, local_maxf3, maxf1, maxf2, maxf3 
          
          INTEGER :: ncout, ncid, varids(3), dimids(3), f_id
          INTEGER :: x_dimid, y_dimid, z_dimid, ux_id, uy_id, uz_id, uvec_id, maxUx_id, maxUy_id, maxUz_id       

          CHARACTER(100) :: parallel_file
          CHARACTER(2) :: RANKtxt
          INTEGER :: fname_len, ii
          INTEGER, DIMENSION(1:3) :: starts, counts
          
!          print *,  "starting save_field_R3toR3_ncdf to file ", file_name

          IF (parallel_data) THEN
             IF (rank==0) THEN
                ncout = nf90_create(file_name, NF90_CLOBBER, ncid=ncid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_def_dim(ncid, "x", n(1), x_dimid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_def_dim(ncid, "y", n(2), y_dimid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_def_dim(ncid, "z", NF90_UNLIMITED, z_dimid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                dimids = (/ x_dimid, y_dimid, z_dimid /)

                ncout = nf90_def_var(ncid, TRIM(f1_name), NF90_DOUBLE, dimids, ux_id)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_def_var(ncid, TRIM(f2_name), NF90_DOUBLE, dimids, uy_id)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_def_var(ncid, TRIM(f3_name), NF90_DOUBLE, dimids, uz_id)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                     
                ncout = nf90_enddef(ncid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_close(ncid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
             END IF
             CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

             starts = (/ 1, 1, rank*local_N+1 /)
             counts = (/ n(1), n(2), local_N /)            
 
             !!--------------------------
             !! START netCDF ROUTINES
             !!--------------------------
             DO ii=0,np-1
                IF (rank==ii) THEN 
                   ncout = nf90_open(file_name, NF90_WRITE, ncid)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

                   ncout = nf90_inq_varid(ncid, TRIM(f1_name), f_id)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                   ncout = nf90_put_var(ncid, f_id, f1, start = starts, count = counts)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                   
                   ncout = nf90_inq_varid(ncid, TRIM(f2_name), f_id)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                   ncout = nf90_put_var(ncid, f_id, f2, start = starts, count = counts)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                   
                   ncout = nf90_inq_varid(ncid, TRIM(f3_name), f_id)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                   ncout = nf90_put_var(ncid, f_id, f3, start = starts, count = counts)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

                   ncout = nf90_close(ncid)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                END IF
                CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo) 
             END DO
             CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
 
          ELSE 
             IF (rank == 0) THEN 
                ALLOCATE( u1(1:n(1),1:n(2),1:n(3)) )
                ALLOCATE( u2(1:n(1),1:n(2),1:n(3)) )
                ALLOCATE( u3(1:n(1),1:n(2),1:n(3)) )
             END IF
             CALL MPI_GATHER(f1, total_local_size, MPI_REAL8, u1, total_local_size, MPI_REAL8, 0, MPI_COMM_WORLD, Statinfo)
             CALL MPI_GATHER(f2, total_local_size, MPI_REAL8, u2, total_local_size, MPI_REAL8, 0, MPI_COMM_WORLD, Statinfo)
             CALL MPI_GATHER(f3, total_local_size, MPI_REAL8, u3, total_local_size, MPI_REAL8, 0, MPI_COMM_WORLD, Statinfo)

             IF (rank==0) THEN
                ncout = nf90_create(file_name, NF90_CLOBBER, ncid)
                ncout = nf90_def_dim(ncid, "x", n(1), x_dimid)
                ncout = nf90_def_dim(ncid, "y", n(2), y_dimid)
                ncout = nf90_def_dim(ncid, "z", n(3), z_dimid)
                dimids =  (/ x_dimid, y_dimid, z_dimid /)
                
                ncout = nf90_def_var(ncid, TRIM(f1_name), NF90_DOUBLE, dimids, ux_id)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

                IF (ncout /= NF90_NOERR) PRINT *, 'data_ops; save_field_R3toR3_ncdf; NOT RIGHT'

                ncout = nf90_def_var(ncid, TRIM(f2_name), NF90_DOUBLE, dimids, uy_id)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_def_var(ncid, TRIM(f3_name), NF90_DOUBLE, dimids, uz_id)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_enddef(ncid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
          
                ncout = nf90_put_var(ncid, ux_id, u1)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_put_var(ncid, uy_id, u2)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_put_var(ncid, uz_id, u3)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

                ncout = nf90_close(ncid)
             
                DEALLOCATE( u1 ) 
                DEALLOCATE( u2 )
                DEALLOCATE( u3 )

             END IF  
             CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

          END IF 
 
        END SUBROUTINE save_field_R3toR3_ncdf

        !==========================================
        ! SAVE FIELD IN netCDF FORMAT
        !==========================================
        SUBROUTINE save_field_R3toRn_ncdf(myfield, dimRange, mynames, file_name)
          
         USE global_variables
         USE netcdf
         IMPLICIT NONE
         INCLUDE "mpif.h"

         REAL(pr), DIMENSION(:,:,:,:), INTENT(IN) :: myfield
         INTEGER, INTENT(IN) :: dimRange
         CHARACTER(len=*), INTENT(IN) :: mynames
         CHARACTER(len=*), INTENT(IN) :: file_name
        
         REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: f
         REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: local_f

         REAL(pr) :: local_maxf, local_minf, maxf, minf
         INTEGER :: ncout, ncid, ndims, nvars, include_parents, dimids(3)
         INTEGER :: x_dimid, y_dimid, z_dimid, f_id       
         INTEGER :: ii, kk, coma_index, start_index, mynames_length
         CHARACTER(10) :: varName
         CHARACTER(50) :: auxName
         CHARACTER(100) :: parallel_file
         CHARACTER(2) :: RANKtxt
         INTEGER :: fname_len, local_nlast_LR
         INTEGER, DIMENSION(1:3) :: starts, counts

!          print *,  "starting save_field_R3toRn_ncdf to file ", file_name
         
         mynames_length = LEN(mynames)
         coma_index = 0
         start_index = 1

         ALLOCATE(local_f(1:n(1),1:n(2),1:local_N))

         if (.not. parallel_data) then
            if(rank==0) print*, "slight warning: save field R3 -> Rn not implemented for non parallel dataops, using parallel dataops (might be inefficient)"
         end if
         IF (rank==0) THEN
            ncout = nf90_create(file_name, NF90_CLOBBER, ncid=ncid)
            IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
            ncout = nf90_def_dim(ncid, "x", n(1), x_dimid)
            IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
            ncout = nf90_def_dim(ncid, "y", n(2), y_dimid)
            IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
            ncout = nf90_def_dim(ncid, "z", NF90_UNLIMITED, z_dimid)
            IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
            dimids = (/ x_dimid, y_dimid, z_dimid /)

            auxName = mynames
            DO kk=1,dimRange
               IF ( kk < dimRange ) THEN
                  coma_index = SCAN(auxName, ",", .FALSE.)
                  varName = auxName(1:coma_index-1)
                  start_index = coma_index+1
               ELSE
                  varName = auxName
               END IF
               ncout = nf90_def_var(ncid, TRIM(varName), NF90_DOUBLE, dimids, f_id)
               IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
               auxName = auxName(start_index:mynames_length)
            END DO
                  
            ncout = nf90_enddef(ncid)
            IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
            ncout = nf90_close(ncid)
            IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
         END IF
         CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

         starts = (/ 1, 1, rank*local_N+1 /)
         counts = (/ n(1), n(2), local_N /)            

         !!--------------------------
         !! START netCDF ROUTINES
         !!--------------------------
         DO ii=0,np-1
            IF (rank==ii) THEN 
               ncout = nf90_open(file_name, NF90_WRITE, ncid)
               IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

               coma_index = 0
               start_index = 1
               auxName = mynames
               DO kk=1,dimRange
                  IF ( kk < dimRange ) THEN
                     coma_index = SCAN(auxName, ",", .FALSE.)
                     varName = auxName(1:coma_index-1)
                     start_index = coma_index+1
                  ELSE
                     varName = auxName
                  END IF
                  auxName = auxName(start_index:mynames_length)

                  local_f = myfield(:,:,:,kk)

                  ncout = nf90_inq_varid(ncid, TRIM(varName), f_id)
                  IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                  ncout = nf90_put_var(ncid, f_id, local_f, start = starts, count = counts)
                  IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
               END DO  

               ncout = nf90_close(ncid)
               IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
            END IF
            CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo) 
         END DO

         DEALLOCATE( local_f )
         CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
         !else
         !end if
       END SUBROUTINE save_field_R3toRn_ncdf
       !!! OLD START !!!
        SUBROUTINE save_field_R3toRn_ncdfOLD(myfield, dimRange, mynames, file_name)
          
          USE global_variables
          USE netcdf
          IMPLICIT NONE
          INCLUDE "mpif.h"

          REAL(pr), DIMENSION(:,:,:,:), INTENT(IN) :: myfield
          INTEGER, INTENT(IN) :: dimRange
          CHARACTER(len=*), INTENT(IN) :: mynames
          CHARACTER(len=*), INTENT(IN) :: file_name
         
          REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: f
          REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: local_f

          REAL(pr) :: local_maxf, local_minf, maxf, minf
          INTEGER :: ncout, ncid, ndims, nvars, include_parents, dimids(3)
          INTEGER :: x_dimid, y_dimid, z_dimid, f_id       
          INTEGER :: ii, kk, coma_index, start_index, mynames_length
          CHARACTER(10) :: varName
          CHARACTER(50) :: auxName
          CHARACTER(100) :: parallel_file
          CHARACTER(2) :: RANKtxt
          INTEGER :: fname_len, local_nlast_LR
          INTEGER, DIMENSION(1:3) :: starts, counts

!          print *,  "starting save_field_R3toRn_ncdf to file ", file_name
          
          mynames_length = LEN(mynames)
          coma_index = 0
          start_index = 1

          ALLOCATE(local_f(1:n(1),1:n(2),1:local_N))

          IF (rank==0) THEN
             ncout = nf90_create(file_name, NF90_CLOBBER, ncid=ncid)
             IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
             ncout = nf90_def_dim(ncid, "x", n(1), x_dimid)
             IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
             ncout = nf90_def_dim(ncid, "y", n(2), y_dimid)
             IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
             ncout = nf90_def_dim(ncid, "z", NF90_UNLIMITED, z_dimid)
             IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
             dimids = (/ x_dimid, y_dimid, z_dimid /)

             auxName = mynames
             DO kk=1,dimRange
                IF ( kk < dimRange ) THEN
                   coma_index = SCAN(auxName, ",", .FALSE.)
                   varName = auxName(1:coma_index-1)
                   start_index = coma_index+1
                ELSE
                   varName = auxName
                END IF
                ncout = nf90_def_var(ncid, TRIM(varName), NF90_DOUBLE, dimids, f_id)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                auxName = auxName(start_index:mynames_length)
             END DO
                     
             ncout = nf90_enddef(ncid)
             IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
             ncout = nf90_close(ncid)
             IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
          END IF
          CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

          starts = (/ 1, 1, rank*local_N+1 /)
          counts = (/ n(1), n(2), local_N /)            
 
          !!--------------------------
          !! START netCDF ROUTINES
          !!--------------------------
          DO ii=0,np-1
             IF (rank==ii) THEN 
                ncout = nf90_open(file_name, NF90_WRITE, ncid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

                coma_index = 0
                start_index = 1
                auxName = mynames
                DO kk=1,dimRange
                   IF ( kk < dimRange ) THEN
                      coma_index = SCAN(auxName, ",", .FALSE.)
                      varName = auxName(1:coma_index-1)
                      start_index = coma_index+1
                   ELSE
                      varName = auxName
                   END IF
                   auxName = auxName(start_index:mynames_length)
 
                   local_f = myfield(:,:,:,kk)

                   ncout = nf90_inq_varid(ncid, TRIM(varName), f_id)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                   ncout = nf90_put_var(ncid, f_id, local_f, start = starts, count = counts)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                END DO  

                ncout = nf90_close(ncid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
             END IF
             CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo) 
          END DO

          DEALLOCATE( local_f )
          CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
 
        END SUBROUTINE save_field_R3toRn_ncdfOLD
        !!! OLD END !!!


        !==========================================
        ! SAVE OPTIMIZATION RESULT LIST
        !==========================================
        SUBROUTINE save_to_optimizationResultList(B, J1, iter, optimizationSuccessful)
         USE global_variables
         IMPLICIT NONE

         real(pr), intent(in) :: B, J1
         integer :: iter
         logical, intent(in) :: optimizationSuccessful
         CHARACTER(200) :: filename
         

         filename = HomeDir//"results_q"//lebesgueQTxt//".dat"

         if(.not.optimizationSuccessful) iter = - iter
         
         if(rank==0 .and. verboseOptimization) print*, "saving results to ", filename
         
         if(rank==0) then
            if(B_list_iterator==bIterOffset+1) then
               OPEN(10, FILE = filename, FORM = 'FORMATTED', STATUS = 'REPLACE')
               !WRITE(10,*) "# B  J1 iter"
               WRITE(10, "(4 G20.12)") "B#", "B", "J1", "iter"
               close(10)
            end if

               
            OPEN(10, FILE = filename, FORM = 'FORMATTED', STATUS = 'OLD', POSITION = 'APPEND')
            WRITE(10, "(A20, 2 ES20.12, I20)") bIterTxt, B, J1, Iter
            CLOSE(10)
         end if

       END SUBROUTINE save_to_optimizationResultList

        !==========================================
        !     SAVE OPTIMIZATION DIAGNOSTICS
        !==========================================
        SUBROUTINE save_diagnostics_optim(myOptimType, iter, tau, beta, J, ener, ens, L2div, dEdt_visc, dEdt_NL, dEdt_Heli)   ! Newly modify April 24, 2017
          USE global_variables
          IMPLICIT NONE

          CHARACTER(len=*), INTENT(IN) :: myOptimType
          REAL(pr), DIMENSION(1:3), INTENT(IN) :: ener, ens
          REAL(pr), INTENT(IN) :: tau, beta, J, L2div, dEdt_visc, dEdt_NL, dEdt_Heli
          INTEGER, INTENT(IN) :: iter

          CHARACTER(200) :: filename
          CHARACTER(2) :: K0txt, E0txt, IGtxt
!          CHARACTER(2) :: WEIGHTtxt
          CHARACTER(4) :: Ntxt

          WRITE(K0txt, '(i2.2)') K0_index
          WRITE(E0txt, '(i2.2)') E0_index
          WRITE(IGtxt, '(i2.2)') iguess
!          WRITE(WEIGHTtxt, '(i2.2)') int_WEIGHT
         
!          filename = "/work/yund0050/MultiObjective_095_01/WEIGHT"//WEIGHTtxt//"_E"//E0txt//"_"//myOptimType//"_IG"//IGtxt//"_iter_info.dat"
          !filename = HomeDir//"_E"//E0txt//"_IG"//IGtxt//"_iter_info.dat"   ! Newly added on May 8, 2017
          filename = ConstraintDir//"iteration-info"//"-B"//trim(bIterTxt)//".dat"   ! trim removes the blank space in "end " for "iterend "
          
          IF (iter == 0) THEN
             OPEN(10, FILE = filename, FORM = 'FORMATTED', STATUS = 'REPLACE')
             !WRITE(10,*) "Iter  Tau  Beta  J  Ener  Ens  Div_L2norm  visc_R  NL_R Heli_R"   ! Newly added H_R, means helicity term in the objective function R
             WRITE(10,*) "Iter  Tau  Beta  J"   ! Newly added H_R, means helicity term in the objective function R
             CLOSE(10)
          END IF

          OPEN(10, FILE = filename, FORM = 'FORMATTED', POSITION = 'APPEND')
          !WRITE(10, "(I5.4, 9 ES20.12)") iter, tau, beta, J, SUM(ener), SUM(ens), L2div, dEdt_visc, dEdt_NL, dEdt_Heli
          WRITE(10, "(I5.4, 3 ES20.12)") iter, tau, beta, J
          CLOSE(10)

        END SUBROUTINE save_diagnostics_optim

         !===========================================================
         ! Create constraint B directory
         !===========================================================
         subroutine initializeConstraintDirectory()
            use global_variables
            use mpi
            implicit none

            character(len=:), allocatable :: constParDir


            !create parent directory
            constParDir = HomeDir//"constraintDirs/"

            call createDirectoryIfNonExistent(constParDir)


            if(bIterTxt=="nan") then
               ! create initial directory
               ConstraintDir = constParDir//"B000-initial/"
            else
               !create constraint directory
               ConstraintDir = constParDir//"B"//bIterTxt//"-"//bTxt//"/"
            end if
            call createDirectoryIfNonExistent(ConstraintDir)
            

         end subroutine initializeConstraintDirectory
 
        !===========================================================
        ! create directory
        !===========================================================
         subroutine createDirectoryIfNonExistent(fullDir)
            use global_variables
            use mpi
            implicit none
            character(len=*), intent(in) :: fullDir
            integer :: status
            logical :: dirExists
            logical :: useBarrier


            if(rank==0) then
               !check if exists already
               inquire(file=trim(fullDir), exist=dirExists)
               if(.not. dirExists) then
                  !print*, "creating dir ", fullDir
                  status = system( "mkdir " // fullDir )
                  if(status /= 0) then
                     print*, "mkdir " // fullDir // " failed"
                  else
                     if(verboseOptimization) print*, "mkdir " // fullDir // " successful"
                  end if
               else
                  if(verboseOptimization) print*, "dir ", fullDir, " already exists"
               end if
            end if

            call mpi_barrier(mpi_comm_world, statinfo)

         end subroutine createDirectoryIfNonExistent

         !===========================================================
         ! create directory
         !===========================================================
          subroutine createDirectoryIfNonExistent2(fullDir)
             use global_variables
             use mpi
             implicit none
             character(len=*), intent(in) :: fullDir
             integer :: status
             logical :: dirExists
             logical :: useBarrier
 


             call mpi_barrier(mpi_comm_world, statinfo)
             call sleep(2)
             call mpi_barrier(mpi_comm_world, statinfo)

 
             if(rank==0) then
                !check if exists already
                inquire(file=trim(fullDir), exist=dirExists)
                if(.not. dirExists) then
                   !print*, "creating dir ", fullDir
                   status = system( "mkdir " // fullDir )
                   !print*, "creating dir ", fullDir, " status", status
                   if(status /= 0) then
                      !print*, "mkdir " // fullDir // " failed"
                   end if
                else
                   !if(verboseOptimization) print*, "dir ", fullDir, " already exists"
                end if
             end if
 
             print*, "g", " rank ", rank, " ", fullDir, " statinfo ", statinfo
 
             call mpi_barrier(mpi_comm_world, statinfo)
            

             call sleep(10)

             print*, "h", " rank ", rank, " ", fullDir, " statinfo ", statinfo
 
          end subroutine createDirectoryIfNonExistent2

 
        !===========================================================
        ! SAVE POINTS IN PHYSICAL SPACE THAT MAKE THE RINGS
        !===========================================================
        SUBROUTINE save_ring_data(ringLoc, numPointsRing, myflag)
          USE global_variables
          IMPLICIT NONE

          REAL(pr), DIMENSION(1:numPointsRing,1:3), INTENT(IN) :: ringLoc
          INTEGER, INTENT(IN) :: numPointsRing
          INTEGER, INTENT(IN) :: myflag

          INTEGER :: i
          CHARACTER(200) :: filename
          CHARACTER(2) :: K0txt, E0txt, IGtxt
!          CHARACTER(2) :: WEIGHTtxt

          WRITE(K0txt, '(i2.2)') K0_index
          WRITE(E0txt, '(i2.2)') E0_index
          WRITE(IGtxt, '(i2.2)') iguess
!          WRITE(WEIGHTtxt, '(i2.2)') int_WEIGHT
         
!          filename = "/work/yund0050/MultiObjective_095_01/WEIGHT"//WEIGHTtxt//"_E"//E0txt//"_maxdEdt_K"//K0txt//"_E"//E0txt//"_IG"//IGtxt//"_ringLocation.dat"
          filename = ConstraintDir//"ringLocation"//"-B"//bIterTxt//".dat"   ! Newly added on May 8, 2017

          IF (myflag==1) THEN
             OPEN(10, FILE = filename, FORM = 'FORMATTED', STATUS = 'OLD', POSITION = 'APPEND')
          ELSE
             OPEN(10, FILE = filename, FORM = 'FORMATTED', STATUS = 'REPLACE')
          END IF
          DO i=1,numPointsRing
             WRITE(10, "(3 ES20.12)") ringLoc(i,1), ringLoc(i,2), ringLoc(i,3)
          END DO
          CLOSE(10)

        END SUBROUTINE save_ring_data



        !===========================================================
        !              SAVE ELAPSED TIME
        !===========================================================
        SUBROUTINE save_elapsed_time(mysystem, myflag, t_ini, t_fin)
          USE global_variables
          IMPLICIT NONE

          CHARACTER(len=*), INTENT(IN) :: mysystem
          INTEGER, INTENT(IN) :: myflag
          INTEGER, DIMENSION(8), INTENT(IN) :: t_ini, t_fin
          INTEGER :: elapsed_minutes, elapsed_hours, elapsed_days

          CHARACTER(200) :: filename
          
          elapsed_days = 0
          elapsed_hours = 0
          elapsed_minutes = 0
 
          IF (t_fin(1)==t_ini(1)) THEN
             IF (t_fin(2)==t_ini(2)) THEN
                elapsed_days = t_fin(3) - t_ini(3)
                elapsed_hours = t_fin(5) - t_ini(5) + 24*elapsed_days
                IF (elapsed_hours < 24) THEN
                   elapsed_days = 0
                END IF
                elapsed_minutes = t_fin(6) - t_ini(6) + 60*elapsed_hours + 1440*elapsed_days
                IF (elapsed_minutes < 60) THEN
                   elapsed_hours = 0
                END IF                      
             ELSE
                IF (t_fin(2)==2 .OR. t_fin(2)==4 .OR. t_fin(2)==6 .OR. t_fin(2)==8 .OR. t_fin(2)==9 .OR. t_fin(2)==11) THEN
                   elapsed_days = 31 + t_fin(3) - t_ini(3)
                   elapsed_hours = t_fin(5) - t_ini(5) + 24*elapsed_days
                   IF (elapsed_hours < 24) THEN
                      elapsed_days = 0
                   END IF
                   elapsed_minutes = t_fin(6) - t_ini(6) + 60*elapsed_hours + 1440*elapsed_days
                   IF (elapsed_minutes < 60) THEN
                      elapsed_hours = 0
                   END IF                      
                ELSEIF (t_fin(2)==5 .OR. t_fin(2)==7 .OR. t_fin(2)==10 .OR. t_fin(2)==12) THEN
                   elapsed_days = 30 + t_fin(3) - t_ini(3)
                   elapsed_hours = t_fin(5) - t_ini(5) + 24*elapsed_days
                   IF (elapsed_hours < 24) THEN
                      elapsed_days = 0
                   END IF
                   elapsed_minutes = t_fin(6) - t_ini(6) + 60*elapsed_hours + 1440*elapsed_days
                   IF (elapsed_minutes < 60) THEN
                      elapsed_hours = 0
                   END IF                      
                ELSEIF (t_fin(2)==3) THEN
                   IF (MOD(t_fin(1),4)==0) THEN
                      elapsed_days = 29 + t_fin(3) - t_ini(3)
                      elapsed_hours = t_fin(5) - t_ini(5) + 24*elapsed_days
                      IF (elapsed_hours < 24) THEN
                         elapsed_days = 0
                      END IF
                      elapsed_minutes = t_fin(6) - t_ini(6) + 60*elapsed_hours + 1440*elapsed_days
                      IF (elapsed_minutes < 60) THEN
                         elapsed_hours = 0
                      END IF                      
                   ELSE
                      elapsed_days = 28 + t_fin(3) - t_ini(3)
                      elapsed_hours = t_fin(5) - t_ini(5) + 24*elapsed_days
                      IF (elapsed_hours < 24) THEN
                         elapsed_days = 0
                      END IF
                      elapsed_minutes = t_fin(6) - t_ini(6) + 60*elapsed_hours + 1440*elapsed_days
                      IF (elapsed_minutes < 60) THEN
                         elapsed_hours = 0
                      END IF                      
                   
                   END IF
                END IF
             END IF
          ELSE
             elapsed_days = 31 + t_fin(3) - t_ini(3)
             elapsed_hours = t_fin(5) - t_ini(5) + 24*elapsed_days
             IF (elapsed_hours < 24) THEN
                elapsed_days = 0
             END IF
             elapsed_minutes = t_fin(6) - t_ini(6) + 60*elapsed_hours + 1440*elapsed_days
             IF (elapsed_minutes < 60) THEN
                elapsed_hours = 0
             END IF                      
          END IF   
  
          SELECT CASE (mysystem)

            CASE ("maxdEdt")  
!              filename = "/scratch/yund0050/MultiObjective_095_01/Timing/computing_time_optimIter.dat"
              filename = HomeDir//"Timing/computing_time_optimIter.dat"    ! Newly added on May 8, 2017
              IF (myflag==0) THEN
                  OPEN(10, FILE = filename, FORM = 'FORMATTED', STATUS = 'REPLACE')
              ELSE
                  OPEN(10, FILE = filename, FORM = 'FORMATTED', STATUS = 'OLD', POSITION = 'APPEND')
              END IF
              WRITE(10,*) "E0 = ", E0, ", K0 = ", K0,", N = ", n(1),", N_proc = ",np
              WRITE(10,*) "   ",elapsed_days," d"
              WRITE(10,*) "   ",elapsed_hours," h"
              WRITE(10,*) "   ",elapsed_minutes," m"
              !WRITE(10,*) "   ",elapsed_seconds," s"  
              CLOSE(10)

          END SELECT

        END SUBROUTINE save_elapsed_time

!============================================================================================
!============================================================================================

        !==============================================================
        !                    READ DATA
        !==============================================================
        SUBROUTINE read_data(mydata, size_x, size_y, filename, myformat)
 
          USE global_variables
          IMPLICIT NONE

          REAL(pr), DIMENSION(:,:), INTENT(OUT) :: mydata
          INTEGER, INTENT(IN) :: size_x, size_y
          CHARACTER(len=*), INTENT(IN) :: filename
          CHARACTER(len=*), INTENT(IN) :: myformat 
          INTEGER :: i=1,j=1

          SELECT CASE (myformat)

          CASE (".dat")
                OPEN (10, FILE = filename, FORM = 'FORMATTED', STATUS = 'OLD')
                DO j=1,size_y
                   DO i=1,size_x
                      READ(10, "(F18.15)") mydata(i,j)   
                   END DO
                END DO
                CLOSE(10)

          CASE (".bin")
                OPEN (10, FILE = filename, FORM = 'UNFORMATTED', STATUS = 'OLD')
                DO j=1,size_y
                   DO i=1,size_x
                      READ(10) mydata(i,j)
                   END DO
                END DO
                CLOSE (10)

          END SELECT
        END SUBROUTINE read_data

        !============================================================
        ! READ VORTICITY IN netCDF FORMAT
        !============================================================
        SUBROUTINE read_field_R3toR3_ncdf(myfield, filename, Fx_txt, Fy_txt, Fz_txt, successful)
          USE global_variables
          USE netcdf
          IMPLICIT NONE
          INCLUDE "mpif.h"

          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(OUT) :: myfield
          CHARACTER(len=*), INTENT(IN) :: filename, Fx_txt, Fy_txt, Fz_txt
         
          REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: local_f, global_f
 
          INTEGER :: ncout, ncid, fid, dimids(3)
          INTEGER :: x_dimid, y_dimid, z_dimid
          INTEGER :: fname_len, ii, nx_ncdf, ny_ncdf, nz_ncdf
          INTEGER, DIMENSION(1:3) :: starts, counts
          logical, intent(out) :: successful



          successful = .true.
          if(rank==0) print*, achar(9), achar(9), achar(9), "loading file ", filename

          IF (parallel_data) THEN
             ALLOCATE( local_f(1:n(1),1:n(2),1:local_N) )

             starts = (/ 1, 1, rank*local_N+1 /)
             counts = (/ n(1), n(2), local_N /)            
 
             !--------------------------
             ! START netCDF ROUTINES
             !--------------------------
             DO ii=0,np-1
                IF (rank == ii) THEN 
                   ncout = nf90_open(filename, NF90_NOWRITE, ncid)
                   if (ncout /= NF90_NOERR .and. successful) then
                     CALL ncdf_error_handle(ncout)
                     successful = .false.
                   end if

                   ncout = nf90_inq_dimid(ncid, "x", x_dimid)
                   if (ncout /= NF90_NOERR .and. successful) then
                     CALL ncdf_error_handle(ncout)
                     successful = .false.
                   end if
                   ncout = nf90_inq_dimid(ncid, "y", y_dimid)
                   if (ncout /= NF90_NOERR .and. successful) then
                     CALL ncdf_error_handle(ncout)
                     successful = .false.
                   end if
                   ncout = nf90_inq_dimid(ncid, "z", z_dimid)
                   if (ncout /= NF90_NOERR .and. successful) then
                     CALL ncdf_error_handle(ncout)
                     successful = .false.
                   end if

                   ncout = nf90_inquire_dimension(ncid, x_dimid, len = nx_ncdf)
                   if (ncout /= NF90_NOERR .and. successful) then
                     CALL ncdf_error_handle(ncout)
                     successful = .false.
                   end if
                   ncout = nf90_inquire_dimension(ncid, y_dimid, len = ny_ncdf)
                   if (ncout /= NF90_NOERR .and. successful) then
                     CALL ncdf_error_handle(ncout)
                     successful = .false.
                   end if
                   ncout = nf90_inquire_dimension(ncid, z_dimid, len = nz_ncdf)
                   if (ncout /= NF90_NOERR .and. successful) then
                     CALL ncdf_error_handle(ncout)
                     successful = .false.
                   end if

                   ncout = nf90_inq_varid(ncid, Fx_txt, fid)
                   if (ncout /= NF90_NOERR .and. successful) then
                     CALL ncdf_error_handle(ncout)
                     successful = .false.
                   end if
                   ncout = nf90_get_var(ncid, fid, local_f, start = starts, count = counts)
                   if (ncout /= NF90_NOERR .and. successful) then
                     CALL ncdf_error_handle(ncout)
                     successful = .false.
                   end if
                   myfield(:,:,:,1) = local_f
 
                   ncout = nf90_inq_varid(ncid, Fy_txt, fid)
                   if (ncout /= NF90_NOERR .and. successful) then
                     CALL ncdf_error_handle(ncout)
                     successful = .false.
                   end if
                   ncout = nf90_get_var(ncid, fid, local_f, start = starts, count = counts)
                   if (ncout /= NF90_NOERR .and. successful) then
                     CALL ncdf_error_handle(ncout)
                     successful = .false.
                   end if
                   myfield(:,:,:,2) = local_f
 
                   ncout = nf90_inq_varid(ncid, Fz_txt, fid)
                   if (ncout /= NF90_NOERR .and. successful) then
                     CALL ncdf_error_handle(ncout)
                     successful = .false.
                   end if
                   ncout = nf90_get_var(ncid, fid, local_f, start = starts, count = counts)
                   if (ncout /= NF90_NOERR .and. successful) then
                     CALL ncdf_error_handle(ncout)
                     successful = .false.
                   end if
                   myfield(:,:,:,3) = local_f
 
                   ncout = nf90_close(ncid)
 
                   DEALLOCATE(local_f)
                END IF
                CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
             END DO     
             CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

          ELSE
             IF (rank == 0) THEN
                ALLOCATE( global_f(1:n(1),1:n(2),1:n(3)) )
             END IF
            
             ALLOCATE( local_f(1:n(1),1:n(2),1:local_N) )

             IF (rank == 0) THEN 
                ncout = nf90_open(filename, NF90_NOWRITE, ncid)
                if (ncout /= NF90_NOERR .and. successful) then
                  CALL ncdf_error_handle(ncout)
                  successful = .false.
                end if

                ncout = nf90_inq_dimid(ncid, "x", x_dimid)
                if (ncout /= NF90_NOERR .and. successful) then
                  CALL ncdf_error_handle(ncout)
                  successful = .false.
                end if
                ncout = nf90_inq_dimid(ncid, "y", y_dimid)
                if (ncout /= NF90_NOERR .and. successful) then
                  CALL ncdf_error_handle(ncout)
                  successful = .false.
                end if
                ncout = nf90_inq_dimid(ncid, "z", z_dimid)
                if (ncout /= NF90_NOERR .and. successful) then
                  CALL ncdf_error_handle(ncout)
                  successful = .false.
                end if

                ncout = nf90_inq_varid(ncid, Fx_txt, fid)
                if (ncout /= NF90_NOERR .and. successful) then
                  CALL ncdf_error_handle(ncout)
                  successful = .false.
                end if
                ncout = nf90_get_var(ncid, fid, global_f)
                if (ncout /= NF90_NOERR .and. successful) then
                  CALL ncdf_error_handle(ncout)
                  successful = .false.
                end if
             END IF
             CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
             CALL MPI_SCATTER(global_f, total_local_size, MPI_REAL8, local_f, total_local_size, MPI_REAL8, 0, MPI_COMM_WORLD, Statinfo)
             myfield(:,:,:,1) = local_f
 
             IF (rank == 0) THEN 
                ncout = nf90_inq_varid(ncid, Fy_txt, fid)
                if (ncout /= NF90_NOERR .and. successful) then
                  CALL ncdf_error_handle(ncout)
                  successful = .false.
                end if
                ncout = nf90_get_var(ncid, fid, global_f)
                if (ncout /= NF90_NOERR .and. successful) then
                  CALL ncdf_error_handle(ncout)
                  successful = .false.
                end if
             END IF
             CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
             CALL MPI_SCATTER(global_f, total_local_size, MPI_REAL8, local_f, total_local_size, MPI_REAL8, 0, MPI_COMM_WORLD, Statinfo)
             myfield(:,:,:,2) = local_f
 
             IF (rank == 0) THEN 
                ncout = nf90_inq_varid(ncid, Fz_txt, fid)
                if (ncout /= NF90_NOERR .and. successful) then
                  CALL ncdf_error_handle(ncout)
                  successful = .false.
                end if
                ncout = nf90_get_var(ncid, fid, global_f)
                if (ncout /= NF90_NOERR .and. successful) then
                  CALL ncdf_error_handle(ncout)
                  successful = .false.
                end if
                ncout = nf90_close(ncid)
             END IF
             CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
             CALL MPI_SCATTER(global_f, total_local_size, MPI_REAL8, local_f, total_local_size, MPI_REAL8, 0, MPI_COMM_WORLD, Statinfo)
             myfield(:,:,:,3) = local_f
 
             IF (rank == 0) THEN
                DEALLOCATE( global_f )
             END IF
             DEALLOCATE( local_f )
             CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
          END IF  
        END SUBROUTINE read_field_R3toR3_ncdf


        !============================================
        ! SAVE LINE MINIMIZATION DATA
        !============================================
        SUBROUTINE save_linemin_data(tA, tB, tC, FA, FB, FC, iter, mysystem, mymode)
          USE global_variables
          IMPLICIT NONE
          INCLUDE "mpif.h"
          
          REAL(pr), INTENT(IN) :: tA, tB, tC, FA, FB, FC
          integer, intent(in) :: iter
          CHARACTER(len=*), INTENT(IN) :: mysystem
          CHARACTER(len=*), INTENT(IN) :: mymode

          CHARACTER(200) :: filename
          CHARACTER(2) :: K0txt, E0txt, IGtxt
!          CHARACTER(2) :: WEIGHTtxt

          WRITE(K0txt,'(i2.2)') K0_index
          WRITE(E0txt,'(i2.2)') E0_index
          WRITE(IGtxt,'(i2.2)') iguess

          call createDirectoryIfNonExistent(ConstraintDir//"tau-data")

          IF (rank==0) THEN
!             filename = "/work/yund0050/MultiObjective_095_01/WEIGHT"//WEIGHTtxt//"_E"//E0txt//"_"//mysystem//"_IG"//IGtxt//"_lineMin_info.dat"
             if(iter<0) then
               filename = ConstraintDir//"tau-data/"//"all-J-info-"//optimizationIterationTxt//".dat"
             else
               filename = ConstraintDir//"tau-data/"//"lineMin-info-"//optimizationIterationTxt//".dat"
             end if
             SELECT CASE (mymode)
               CASE ("replace")
                 OPEN(10, FILE = filename, FORM = 'FORMATTED', STATUS = 'REPLACE')
                 WRITE(10, "(G5.4, 6 G20.12)") "iter", "tA", "tB", "tC", "FA", "FB", "FC"
               CASE ("append")
                OPEN(10, FILE = filename, FORM = 'FORMATTED', STATUS = 'OLD', POSITION = 'APPEND')
             END SELECT
             WRITE(10, "(I5.4, 6 ES20.12)") iter, tA, tB, tC, FA, FB, FC
             CLOSE(10)
          END IF 
          CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

        END SUBROUTINE save_linemin_data


        !============================================
        ! Save results from kappa test
        !============================================
        SUBROUTINE save_kappa_test(eps, inner_prod, deltaJ, kappa, kappa_adj, adj_factor, myindex, identifier, useAdjustedKappaTest)
          USE global_variables
          IMPLICIT NONE

          REAL(pr), INTENT(IN) :: eps, kappa_adj, kappa, adj_factor, inner_prod, deltaJ
          INTEGER, INTENT(IN) :: myindex
          CHARACTER(len=*), INTENT(IN) :: identifier
          character(len=:), allocatable :: kappaDir
          character(len=:), allocatable :: filePath
          logical, intent(in) :: useAdjustedKappaTest

          
          
          kappaDir = constraintDir//"kappaTest/"
          call createDirectoryIfNonExistent(kappaDir)
         

          filePath = kappaDir//"kappa"//"_B"//bIterTxt//"_"//identifier//".dat"
          filePath=trim(filePath)

          if(rank==0) then
            if(useAdjustedKappaTest) then
              IF (myindex==1) THEN 
                 OPEN (10, FILE = filePath, FORM = 'FORMATTED', STATUS = 'REPLACE')
                 WRITE(10, "(9 G20.12)") "eps", "deltaJ", "deltaJ/eps", "inner_prod", "kappa_adj", "LOG10|kap_adj-1|", "kappa", "LOG10|kap-1|", "adj_factor"
                 CLOSE(10)
              END IF 
              OPEN (10, FILE = filePath, FORM = 'FORMATTED', STATUS = 'OLD', POSITION = 'APPEND')
              WRITE(10, "(9 ES20.12)") eps, deltaJ, deltaJ/eps, inner_prod, kappa_adj, LOG10(ABS(kappa_adj - 1.0_pr)), kappa, LOG10(ABS(kappa - 1.0_pr)), adj_factor
              CLOSE(10)
            else
              IF (myindex==1) THEN 
                 OPEN (10, FILE = filePath, FORM = 'FORMATTED', STATUS = 'REPLACE')
                 WRITE(10, "(6 G20.12)") "eps", "deltaJ", "deltaJ/eps", "inner_prod", "kappa", "LOG10|kap-1|"
                 WRITE(10, "(6 ES20.12)") eps, deltaJ, deltaJ/eps, inner_prod, kappa, LOG10(ABS(kappa - 1.0_pr))
                 CLOSE(10)
              END IF 
              OPEN (10, FILE = filePath, FORM = 'FORMATTED', STATUS = 'OLD', POSITION = 'APPEND')
              WRITE(10, "(6 ES20.12)") eps, deltaJ, deltaJ/eps, inner_prod, kappa, LOG10(ABS(kappa - 1.0_pr))
              CLOSE(10)
            end if
          end if
 
        END SUBROUTINE save_kappa_test

        !===============================================
        ! NETCDF ERROR HANDLE ROUTINE
        !===============================================
        SUBROUTINE ncdf_error_handle(nerror)
          USE global_variables
          USE netcdf
          IMPLICIT NONE
          INCLUDE "mpif.h"

          INTEGER, INTENT(IN) :: nerror
          CHARACTER(80) :: error_string
          CHARACTER(2) :: K0txt, E0txt, IGtxt

          WRITE(K0txt,'(i2.2)') K0_index
          WRITE(E0txt,'(i2.2)') E0_index
          WRITE(IGtxt,'(i2.2)') iguess
 
          error_string = NF90_STRERROR(nerror)

          print*, " Error reading netCDF file. ", error_string

          IF (rank==0) THEN
             OPEN(10, FILE=constraintDir//"/netcdf.log", POSITION='APPEND')
             WRITE(10,*) " Error reading netCDF file. "//error_string
             CLOSE(10)
          END IF
          CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

        END SUBROUTINE ncdf_error_handle


        !===============================================
        ! OPTIMIZATION ERROR HANDLE ROUTINE
        !===============================================
        SUBROUTINE optim_error_handle(nerror)
          USE global_variables
          IMPLICIT NONE
          INCLUDE "mpif.h"

          INTEGER, INTENT(IN) :: nerror
          CHARACTER(80) :: error_string
          CHARACTER(2) :: K0txt, E0txt, IGtxt

          WRITE(K0txt,'(i2.2)') K0_index
          WRITE(E0txt,'(i2.2)') E0_index
          WRITE(IGtxt,'(i2.2)') iguess

          SELECT CASE (nerror)
            CASE (1)
               error_string = " maxdEdt: Going uphill... Verify gradient!"
            CASE (2)
               error_string = " maxdEdt: Could not bracket minimum."
            CASE (3)
               error_string = " maxdEdt: Decreasing tau..."
            CASE (11)
               error_string = " FixK0E0: Going uphill... Verify gradient!"
            CASE (12)
               error_string = " FixK0E0: Could not bracket minimum."
            CASE (13)
               error_string = " FixK0E0: Decreasing tau..."
            CASE (15)
               error_string = " FixK0E0: Could not bracket minimum. Trying tau = TauMax..."

          END SELECT


          IF (rank==0) THEN   
            OPEN(10, FILE=ConstraintDir//"optimization.log", POSITION='APPEND')
            !OPEN(10, FILE=HomeDir//"/LOGFILE_maxdEdtHeli_E"//E0txt//"_IG"//IGtxt//"_info.log", POSITION='APPEND')
            WRITE(10,*) "      Error during optimization. "//error_string
            CLOSE(10)
          END IF 

        END SUBROUTINE optim_error_handle

        !===============================================
        ! OPTIMIZATION ERROR HANDLE ROUTINE
        !===============================================
        SUBROUTINE optim_msg_handle(nmsg)
          USE global_variables
          IMPLICIT NONE
          INCLUDE "mpif.h"

          INTEGER, INTENT(IN) :: nmsg
          CHARACTER(80) :: msg_string
          CHARACTER(2) :: K0txt, E0txt, IGtxt

          WRITE(K0txt,'(i2.2)') K0_index
          WRITE(E0txt,'(i2.2)') E0_index
          WRITE(IGtxt,'(i2.2)') iguess
 
          SELECT CASE (nmsg)
            CASE (0)
               msg_string = "      Cost functional not increasing."
            CASE (1)
               msg_string = "      Optimization terminated." 
            CASE (10)
               msg_string = "   Starting FixK0E0..."
            CASE (11)
               msg_string = "   FixK0E0 OK!"
            CASE (12)
               msg_string = "   Could not FixK0E0... Stop optimization!"
            CASE (13)
               msg_string = "   Could not move to constraint... Stop optimization!"
            CASE (14)
               msg_string = "   Could not move to constraint... Iterations continue...!" 
            CASE (20)
               msg_string = "      Starting mnbrak..."
            CASE (21)
               msg_string = "      mnbrak OK!"
            CASE (22)
               msg_string = "      mnbrak error... going uphill, verify gradient!"
            CASE (30)
               msg_string = "      Starting Brent method..."
            CASE (31)
               msg_string = "      Brent method OK!"
            CASE (32)
               msg_string = "      Optimal tau is too large!"
            CASE (41)
               msg_string = "      warning tau > tau_max/1000000"
            CASE (42)
               msg_string = "      WARNING abs(inner(gradJproj,normalHs)) > 10e-15, should be 0"
            CASE (43)
               msg_string = "      WARNING ||div(gradJproj)||_2^2 > 10e-15, should be 0"
            CASE (44)
               msg_string = "      WARNING function that should be average free is not"

               

          END SELECT

          IF (rank==0) THEN   
            !OPEN(10, FILE=HomeDir//"/LOGFILE_maxdEdtHeli_E"//E0txt//"_IG"//IGtxt//"_info.log", POSITION='APPEND')
            OPEN(10, FILE=ConstraintDir//"optimization"//"-B"//bIterTxt//".log", POSITION='APPEND')
            WRITE(10,*) msg_string
            CLOSE(10)
          END IF 
          CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

        END SUBROUTINE optim_msg_handle

        !=====================================================
        ! DOUBLE THE RESOLUTION OF A GIVEN FIELD N1= 2*N0
        !=====================================================
        SUBROUTINE interpolate3D(f0, f1, N0, N1)
          USE global_variables
          IMPLICIT NONE
          INCLUDE "mpif.h"

          REAL(pr), DIMENSION(:,:,:), INTENT(IN) :: f0
          REAL(pr), DIMENSION(:,:,:), INTENT(OUT) :: f1
          INTEGER, INTENT(IN) :: N0, N1

          INTEGER :: ii1, jj1, kk1, ii0, jj0, kk0, ii0_aux, jj0_aux, kk0_aux, local_nlast1, stride
          REAL(pr) :: dx1, dy1, dz1, dx0, dy0, dz0, xp, yp, zp
          LOGICAL :: same_plane, same_line, same_point         

          dx0 = 1.0_pr/REAL(N0,pr)
          dy0 = dx0
          dz0 = dx0
          dx1 = 1.0_pr/REAL(N1,pr)
          dy1 = dx1
          dz1 = dx1

          IF (N1>N0) THEN
             DO kk1=1,N1
                IF (MOD(kk1-1,2)==0) THEN
                   same_plane = .TRUE.
                ELSE
                   same_plane = .FALSE.
                END IF
              
                DO jj1=1,N1
                   IF (MOD(jj1-1,2)==0) THEN
                      same_line = .TRUE.
                   ELSE
                      same_line = .FALSE.
                   END IF
                   
                   DO ii1=1,N1
                      IF (MOD(jj1-1,2)==0) THEN
                         same_point = .TRUE.
                      ELSE
                         same_point = .FALSE.
                      END IF
                      
                      IF ( same_plane .AND. same_line .AND. same_point ) THEN
                         ii0 = (ii1+1)/2
                         jj0 = (jj1+1)/2
                         kk0 = (kk1+1)/2
                         f1(ii1,jj1,kk1) = f0(ii0,jj0,kk0)

                      ELSEIF ( same_plane .AND. same_line .AND. .NOT. same_point) THEN
                         ii0 = (ii1)/2
                         jj0 = (jj1+1)/2
                         kk0 = (kk1+1)/2
                         IF (ii0 == N0) THEN
                            ii0_aux = 1
                         ELSE
                            ii0_aux = ii0+1
                         END IF
                         f1(ii1,jj1,kk1) = 0.5_pr*( f0(ii0,jj0,kk0) + f0(ii0_aux,jj0,kk0) )

                      ELSEIF ( same_plane .AND. .NOT. same_line .AND. same_point) THEN
                         ii0 = (ii1+1)/2
                         jj0 = (jj1)/2
                         kk0 = (kk1+1)/2
                         IF (jj0 == N0) THEN
                            jj0_aux = 1
                         ELSE
                            jj0_aux = jj0+1
                         END IF
                         f1(ii1,jj1,kk1) = 0.5_pr*( f0(ii0,jj0,kk0) + f0(ii0,jj0_aux,kk0) )

                      ELSEIF ( same_plane .AND. .NOT. same_line .AND. .NOT. same_point) THEN
                         ii0 = (ii1)/2
                         jj0 = (jj1)/2
                         kk0 = (kk1+1)/2
                         IF (ii0 == N0) THEN
                            ii0_aux = 1
                         ELSE
                            ii0_aux = ii0+1
                         END IF
                         IF (jj0 == N0) THEN
                            jj0_aux = 1
                         ELSE
                            jj0_aux = jj0+1
                         END IF
                         f1(ii1,jj1,kk1) = 0.25_pr*( f0(ii0,jj0,kk0) + f0(ii0_aux,jj0,kk0) + f0(ii0,jj0_aux,kk0) + f0(ii0_aux,jj0_aux,kk0) )

                      ELSEIF ( .NOT. same_plane .AND. same_line .AND. same_point) THEN
                         ii0 = (ii1+1)/2
                         jj0 = (jj1+1)/2
                         kk0 = (kk1)/2
                         IF (kk0 == N0) THEN
                            kk0_aux = 1
                         ELSE
                            kk0_aux = kk0+1
                         END IF
                         f1(ii1,jj1,kk1) = 0.5_pr*( f0(ii0,jj0,kk0) + f0(ii0,jj0,kk0_aux) )

                      ELSEIF ( .NOT. same_plane .AND. same_line .AND. .NOT. same_point) THEN
                         ii0 = (ii1)/2
                         jj0 = (jj1+1)/2
                         kk0 = (kk1)/2
                         IF (ii0 == N0) THEN
                            ii0_aux = 1
                         ELSE
                            ii0_aux = ii0+1
                         END IF
                         IF (kk0 == N0) THEN
                            kk0_aux = 1
                         ELSE
                            kk0_aux = kk0+1
                         END IF
                         f1(ii1,jj1,kk1) = 0.25_pr*( f0(ii0,jj0,kk0) + f0(ii0_aux,jj0,kk0) + f0(ii0,jj0,kk0_aux) + f0(ii0_aux,jj0,kk0_aux) )

                      ELSEIF ( .NOT. same_plane .AND. .NOT. same_line .AND. same_point) THEN
                         ii0 = (ii1+1)/2
                         jj0 = (jj1)/2
                         kk0 = (kk1)/2
                         IF (jj0 == N0) THEN
                            jj0_aux = 1
                         ELSE
                            jj0_aux = jj0+1
                         END IF
                         IF (kk0 == N0) THEN
                            kk0_aux = 1
                         ELSE
                            kk0_aux = kk0+1
                         END IF
                         f1(ii1,jj1,kk1) = 0.25_pr*( f0(ii0,jj0,kk0) + f0(ii0,jj0_aux,kk0) + f0(ii0,jj0,kk0_aux) + f0(ii0,jj0_aux,kk0_aux) )

                      ELSEIF ( .NOT. same_plane .AND. .NOT. same_line .AND. .NOT. same_point) THEN
                         ii0 = (ii1)/2
                         jj0 = (jj1)/2
                         kk0 = (kk1)/2
                         IF (ii0 == N0) THEN
                            ii0_aux = 1
                         ELSE
                            ii0_aux = ii0+1
                         END IF
                         IF (jj0 == N0) THEN
                            jj0_aux = 1
                         ELSE
                            jj0_aux = jj0+1
                         END IF
                         IF (kk0 == N0) THEN
                            kk0_aux = 1
                         ELSE
                            kk0_aux = kk0+1
                         END IF
                         f1(ii1,jj1,kk1) = 0.125_pr*( f0(ii0,jj0,kk0) + f0(ii0_aux,jj0,kk0) + f0(ii0,jj0_aux,kk0) + f0(ii0,jj0,kk0_aux) + &
                                                      f0(ii0,jj0_aux,kk0_aux) + f0(ii0_aux,jj0,kk0_aux) + f0(ii0_aux,jj0_aux,kk0) + f0(ii0_aux,jj0_aux,kk0_aux) )
                      END IF
                   END DO
                END DO
             END DO
          ELSEIF (N0>N1) THEN
             stride = N0/N1
             local_nlast1 = N1/np 

             kk0 = 1 
             DO kk1=1,local_nlast1
                jj0 = 1
                DO jj1=1,N1
                   ii0 = 1 
                   DO ii1=1,N1
                      f1(ii1,jj1,kk1) = f0(ii0,jj0,kk0)
                      ii0 = ii0+stride
                   END DO
                   jj0 = jj0+stride
                END DO 
                kk0 = kk0+stride
             END DO
          END IF

        END SUBROUTINE interpolate3D 


END MODULE 
