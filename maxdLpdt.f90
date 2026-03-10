!-----------------------------------------------------!
! Program used to maximize dEdt                       !
! in the 3D Navier-Stokes system,                     !
! using 2 constraints.                                !
!                                                     !
! Parallel version, Complex-Complex FFT               !   
!                                                     !
! November, 2014.                                     !
!                                                     !
! Author: Diego Ayala                                 !
! Department of Mathematics and Statistics            !
! McMaster University                                 !
!-----------------------------------------------------!
 
   PROGRAM main

      USE global_variables
      USE data_ops
      use fftwfunction
      USE function_ops
      USE optimization
      USE mpi
      IMPLICIT NONE
      !=============================================
      ! Declare parameters
      !=============================================
      
      

      !=============================================
      ! MPI Initialization; and fftw mpi initialization
      !=============================================
      CALL MPI_INIT(Statinfo)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,Statinfo)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,np,Statinfo)
      call fftw_mpi_init()                     ! Newly added in Jan 17, but not sure whether it should be here

      !=============================================
      ! Read in runtime parameters; input manually parameters' values.
      ! Read in code parameters; pass input values to code's parameters. 
      ! Newly added in Jan 23
      !=============================================

      !=============================================
      ! Initial Data
      ! iguess 0 = load FRT_N256E500T017_Uvec_fwdTE0220
      ! iguess 1 = load sines
      ! iguess 2 = load random a
      ! iguess 3 = load random smooth a
      ! iguess 4 = load random expSpec a
      ! iguess 5 = load random polySpec a
      ! iguess 6 = load random k a
      ! iguess 7 = save random
      ! iguess 9 = load temp &loadTempFunctionName
      ! iguess 50 = Arnold-Beltrami-Childress, ...
      !=============================================
      iguess = 9

      if(iguess==9) then
         loadTempFunctionName = "u_result_0303_q9_n1024_B017_iterend.nc"
         standardParams = .true.                         ! use the standard parameters values or is it a test run with strange parameters


         call set_q_resol_bIterOffset_optimIterOffsets(loadTempFunctionName)
         !lebesgueQ = automatically now
         !bIterOffset = automatically now                ! should match the loaded iteration or 0 if new
                                                                  !0 = new optimization
                                                                  !8 if continuing u_result_B009_iter1500.nc
                                                                  !9 if u_result_B009_iterend.nc
         !optimizationIterOffset = automatically now     ! should match the loaded iter or 0 if new
                                                                  !0 = new optimization
                                                                  !1500 if continuing u_result_B009_iter1500.nc
                                                                  !0 if u_result_B009_iterend.nc
                                                                  !just for documentation how many iterations were needed
      else
         lebesgueQ = 9.0_pr
         resol = 256
         standardParams = .false.                         ! use the standard parameters values or is it a test run with strange parameters
         loadTempFunctionName = "iguess00"      ! allocate string resources
         write(loadTempFunctionName, '(A6,I2)') "iguess",iguess
         bIterOffset = 0
         optimizationIterOffset = 0
      end if

      IF (rank==0) print*, "start"

      call setStandardParams()

      if(.not.standardParams)then
         qContinuation = .false.
         if(qContinuation) then
            !!! q continuation !!!

            !!! manual parameters !!!

            if(.true.) then
               !!! 4 to 3 !!!
               scratchPath = "/home/fabianbl/scratch/"
               qContNcFileFolder = scratchPath//"8_1024_production/q4/5/output/ncFiles"
               !qContNcFileFolder = scratchPath//"7_qCont/c_1024_4to3/1/output/ncFiles"
               numberOfqValues = 2
               qStartOffset = 0                 ! offset if wanna start for first file at different q. could be since to continue a sim that ended within a file and different q values
               qContStart = 4.0_pr
               qContEnd = 3.7_pr
               bIterRangeStart = -1             ! specify range of used original b values, negative = all
               bIterRangeEnd = -1
            end if

            if(.false.) then
               !!! 5 to 4 confirmation !!!
               scratchPath = "/home/fabianbl/scratch/"
               !qContNcFileFolder = scratchPath//"5_beta_512_10-5/q5/1/output/ncFiles"
               qContNcFileFolder = scratchPath//"7_qCont/a5_5to4/input"
               numberOfqValues = 5
               qStartOffset = 1                 ! offset if wanna start for first file at different q. could be since to continue a sim that ended within a file and different q values
               qContStart = 5.0_pr
               qContEnd = 4.0_pr
               bIterRangeStart = 19             ! specify range of used original b values, negative = all
               bIterRangeEnd = 21
            end if

            if(.false.) then
               !!! calc enstrophy !!!
               scratchPath = "/home/fabianbl/scratch/"
               qContNcFileFolder = scratchPath//"5_beta_512_10-5/q9/1/output/ncFiles"
               numberOfqValues = 0
               qStartOffset = 0                 ! offset if wanna start for first file at different q. could be since to continue a sim that ended within a file and different q values
               qContStart = 5.0_pr
               qContEnd = 5.0_pr
               bIterRangeStart = 13             ! specify range of used original b values, negative = all
               bIterRangeEnd = 15

               ! skip optimization and kappa test
               MAX_ITER = 0
               kappaTest = .false.
            end if
            
            
            
            
            
            
            
            
            !!! end manual parameters !!!







            if(qContNcFileFolder(len(qContNcFileFolder):len(qContNcFileFolder)) /= "/") then
               qContNcFileFolder = qContNcFileFolder//"/"                                    ! add / if not there in qContNcFileFolder
            end if
            fileList = getListofFilesInDirContainingSearchString(qContNcFileFolder, "iterend")
            call set_q_resol_bIterOffset_optimIterOffsets(fileList(1))
            use_e_u_auto_for_q_less_4 = .true.
            allocate(qContqValues(1:numberOfqValues))
            do B_list_iterator=1,numberOfqValues
               qContqValues(B_list_iterator) = qContStart*((qContEnd/qContStart)**(real(B_list_iterator,pr)/real(numberOfqValues, pr)))
            end do

         else
            !!! normal constraint continuation !!!
            !!! change parameters for test runs !!!

            if(.false.) then
               if(abs(lebesgueQ-3.0_pr)<MACH_EPSILON) then
                  save_uvecEveryXiteration = 10
                  save_dEveryXiteration = 10
                  save_gradL2EveryXiteration = 10
                  save_spectraEveryXiteration = 1
               end if
            end if

            if(.false.) then
               OPTIM_TOL = 1.0e-6_pr
            end if

            if(.true.) then
               !!! overwrite B values !!!
               allocate( B_list_overwrite(1:32+bIterOffset) )
               B_list_overwrite(1) = 1.0_pr
               do B_list_iterator=2,size(B_list_overwrite)
                  constraintB = B_list_overwrite(B_list_iterator-1)*10**(1.0_pr/8.0_pr)
                  if(constraintB>240.0_pr) then
                     constraintB = B_list_overwrite(B_list_iterator-1)*10**(1.0_pr/24.0_pr)
                  end if
                  B_list_overwrite(B_list_iterator) = constraintB
               end do
            end if
            
         end if
      end if

      !=============================================
      !--Initialize parameters values
      !=============================================
      CALL initialize()
      call init_fft()


      !=============================================
      !-- Initialize Data Folders
      !=============================================
      call createDirectoryIfNonExistent(HomeDir)
      call createDirectoryIfNonExistent(ncDir)
      call saveUsedParams(HomeDir)

      !=============================================
      ! Set initial Uvec
      !=============================================
      call initial_guess
      IF (.false.) call save_field_R3toR3_ncdf(Uvec(:,:,:,1), Uvec(:,:,:,2), Uvec(:,:,:,3), "Ux", "Uy", "Uz", ncDir//"initial_u.nc", "netCDF")
      if (.true.) call initializeConstraintDirectory
      if (.true.) call calculateSaveSpectrum(uvec,"initial_u")

      if (rank == 0 .and. ((abs(viscCoefficient-1.0)>MACH_EPSILON))) then
         print*, "WARNING viscCoefficient ",viscCoefficient, "not 1"
      end if
      if (rank == 0 .and. ((abs(pressureCoefficient-1.0)>MACH_EPSILON))) then
         print*, "WARNING pressureCoefficient ", pressureCoefficient, "not 1"
      end if



      if(qContinuation) then
         !=========================================================
         ! q continuation ?
         !=========================================================
         
         fileList = getListofFilesInDirContainingSearchString(qContNcFileFolder, "iterend")
         fileListTemp = fileList

         if(size(fileList)<1) then
            print*, "ERROR NO FILES FOUND IN FOLDER ", qContNcFileFolder, " containing the search string, EXITING"
            call exit(1) 
         end if

         tempInt = 0
         do fileNumber = 1,size(fileList)
            tempLogical = .true.
            B_list_iterator = getBiterfromFileName(fileList(fileNumber))
            if(bIterRangeStart>0) then
               if(B_list_iterator<bIterRangeStart) then
                  tempLogical = .false.
               end if
            end if
            if(bIterRangeEnd>0) then
               if(B_list_iterator>bIterRangeEnd) then
                  tempLogical = .false.
               end if
            end if
            if(tempLogical) then
               tempInt = tempInt + 1
               fileListTemp(tempInt) = fileList(fileNumber)
            end if
         end do
         deallocate(fileList)
         fileList = fileListTemp(1:tempInt)


         if(rank==0) print*, "fileList, n=", size(fileList)
         do fileNumber = 1,size(fileList)
            if(rank==0) print*, achar(9), fileList(fileNumber)
         end do

         fileNameResult = HomeDir//"qContData.dat"

         allocate( qContinuationResults(1:numberOfqValues+1, 1:size(fileList)) )

         !if(rank==0) print*, "number of files in dir ", qContNcFileFolder, " is ", size(fileList)
         do fileNumber = 1,size(fileList)
            write(fileNumberText, '(i3.3)') fileNumber
            B_list_iterator = getBiterfromFileName(fileList(fileNumber))
            write(bIterTxt, '(i3.3)') B_list_iterator
            
            !load function
            inputDir = qContNcFileFolder
            loadTempFunctionName = fileList(fileNumber)
            iguess = 9
            call initial_guess
            call setHostNameAndInputDir   ! resets inputdir for kappa test


            do qNumber=qStartOffset,numberOfqValues

               if(qNumber==0) then
                  lebesgueQ = getQfromFileName(loadTempFunctionName)
               else
                  lebesgueQ = qContqValues(qNumber)
               end if
               call set_e_u_or_just_u()
               call setLebesgueQtext()

               constraintB = calc_global_Lq_norm(uvec)
               write(bTxt, '(ES7.1)') constraintB

               call initializeConstraintDirectory
               call saveUsedParams(ConstraintDir)
               

               !if(qNumber/=0) then
               !   call maxdLqdt
               !else
               !   iter = 0
               !end if

               call maxdLqdt

               qContinuationResults(qNumber+1, fileNumber) = eval_J(Uvec, "maxdLqdt")

               if(rank==0) then
                  fileNameResByB = HomeDir//"results_qCont"//"_B"//bIterTxt//".dat"
                  if(qNumber==qStartOffset) then
                     open(11, FILE = fileNameResByB, FORM = 'FORMATTED', STATUS = 'REPLACE')
                     WRITE(11, "(A20, 2 A20, A20, A20)") "q", "B", "J1", "Iter", "#B"
                     close(11)
                  end if
                  open(11, FILE = fileNameResByB, FORM = 'FORMATTED', STATUS = 'OLD', POSITION = 'APPEND')
                  write(11, "(ES20.12, 2 ES20.12, I20, A20)") lebesgueQ, constraintB, qContinuationResults(qNumber+1, fileNumber), Iter, bIterTxt
                  close(11)
               end if

               !!! by q already done in normal optimization
               !if(rank==0) then
               !   fileNameResByQ = HomeDir//"results_qCont"//"_q"//lebesgueQTxt//".dat"
               !   if(fileNumber==1) then
               !      open(12, FILE = fileNameResByQ, FORM = 'FORMATTED', STATUS = 'REPLACE')
               !      WRITE(12, "(A20, 2 A20, A20, A20)") "#B", "B", "J1", "Iter", "q"
               !      close(12)
               !   end if
               !   open(12, FILE = fileNameResByQ, FORM = 'FORMATTED', STATUS = 'OLD', POSITION = 'APPEND')
               !   write(12, "(A20, 2 ES20.12, I20, ES20.12)") bIterTxt, constraintB, qContinuationResults(qNumber+1, fileNumber), Iter, lebesgueQ
               !   close(12)
               !end if
               if(qStartOffset>0) qStartOffset = 0 ! reset qStartOffset to 0 so that for the other files do all q values 
            end do

            
            if(rank==0) then
               open(10, FILE = fileNameResult, FORM = 'FORMATTED', STATUS = 'REPLACE')
               WRITE(10, "(A20)", advance="no") ""
               do qNumber = 0, numberOfqValues
                  if(qNumber<10) then
                     WRITE(10, "(A19)", advance="no") "q"
                     WRITE(10, "(I1)", advance="no") qNumber
                  else if(qNumber<100) then
                     WRITE(10, "(A18)", advance="no") "q"
                     WRITE(10, "(I2)", advance="no") qNumber
                  else
                     WRITE(10, "(A15)", advance="no") "q"
                     WRITE(10, "(I5)", advance="no") qNumber
                  end if
               end do
               write(10, "(A50)") ""
               WRITE(10, "(A20)", advance="no") "B#"
               WRITE(10, "(ES20.12)", advance="no") qContStart
               do qNumber = 1, numberOfqValues
                  WRITE(10, "(ES20.12)", advance="no") qContqValues(qNumber)
               end do
               write(10, "(A50)") "initial filename"
               close(10)


               open(10, FILE = fileNameResult, FORM = 'FORMATTED', STATUS = 'OLD', POSITION = 'APPEND')
               WRITE(10, "(A20)", advance="no") bIterTxt
               do qNumber = 0, numberOfqValues
                  WRITE(10, "(ES20.12)", advance="no") qContinuationResults(qNumber+1, fileNumber)
               end do
               write(10, "(A50)") fileList(fileNumber)
               close(10)
            end if

         end do


         


      else
         !=========================================================
         ! constraint continuation
         !=========================================================
         
         ! Create B value and result list
         allocate( B_list(1:32+bIterOffset) )
         allocate( optimizationResultList(1:3,0:size(B_list)) )

         !elseif(abs(lebesgueQ-3.0_pr)<MACH_EPSILON) then
         if(abs(lebesgueQ-3.0_pr)<MACH_EPSILON) then
            B_list(1) = 1.0_pr
            do B_list_iterator=2,size(B_list)
               constraintB = B_list(B_list_iterator-1)*10**(1.0_pr/4.0_pr)
               if(constraintB>100.0_pr) constraintB = B_list(B_list_iterator-1)*10**(1.0_pr/16.0_pr)
               if(constraintB>110.0_pr) constraintB = B_list(B_list_iterator-1)*10**(1.0_pr/32.0_pr)
               B_list(B_list_iterator) = constraintB
            end do
         else
            B_list(1) = 1.0_pr
            do B_list_iterator=2,size(B_list)
               constraintB = B_list(B_list_iterator-1)*10**(1.0_pr/8.0_pr)
               B_list(B_list_iterator) = constraintB
            end do
         end if

         if(.not.standardParams) then
            if(allocated(B_list_overwrite)) then
               B_list(:) = B_list_overwrite(:)
            end if
         end if

         do B_list_iterator=1,size(B_list)
            if(rank==0) print*, "B", B_list_iterator, "=", B_list(B_list_iterator)
         end do



         !=========================================================
         ! Loop over different values of constraint B values
         !=========================================================
         do B_list_iterator = 1+bIterOffset,size(B_list)
            constraintB = B_list(B_list_iterator)
            write(bIterTxt, '(i3.3)') B_list_iterator
            write(bTxt, '(ES7.1)') constraintB
            if(rank==0) print*, "constraint B"//bIterTxt, constraintB
            call initializeConstraintDirectory
            !=========================================================
            ! OPTIMIZE !
            !=========================================================
            call maxdLqdt
         END DO
      end if



      call fft_deallocate
      call fftw_mpi_cleanup()          ! Newly Added in Jan 17, but not sure whether it should be here; March 20, 2017
      CALL mpi_finalize(Statinfo)
      
      if (rank==0) print*, "end"
 
   END PROGRAM main
