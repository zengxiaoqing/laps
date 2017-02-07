!-------------------------------------------------------------------          !
! subroutines degrib_nav (to set up array sizes) and degrib_data (to read)    !
! to get all the navigation and parm info that is needed.                     !
!                                                                             !
! Externals:                                                                  !
!    Module TABLE                                                             !
!    Module GRIDINFO                                                          !
!    Module FILELIST                                                          !
!    Subroutine READ_NAMELIST                                                 !
!    Subroutine PARSE_TABLE                                                   !
!    Subroutine CLEAR_STORAGE                                                 !
!    Subroutine RD_GRIB                                                       !
!    Subroutine RD_GRIB2                                                      !
!    Subroutine PUT_STORAGE                                                   !
!    Subroutine CCLOSE                                                        !
!                                                                             !
! Kevin Manning, NCAR/MMM  - original 'pregrid' code, 1998-2001               !
! Paula McCaslin, NOAA
! adapted for SI, 2004                                                        !
! adapted for WPS, 2006                                                       !
! adapted for LAPS, 2007                                                      !
!                                                                             !
!-------------------------------------------------------------------          !

subroutine degrib_nav(gribflnm, vtablefn, nx, ny, nz, &
                      gproj, dx, dy, truelat1, truelat2, stdlon, cgrddef, &
                      cross_dateline, sw1, sw2, ne1, ne2, l_full_run, istatus)

  use table
  use gridinfo
  use storage_module
  use filelist
  use datarray
  use module_debug
  use stringutil
  use map_utils

  implicit none

  character(LEN=256) :: gribflnm
  character(LEN=256) :: vtablefn
  character(LEN=19)  :: hdate
  character(LEN=9)   :: field
  character(LEN=2)   :: gproj
  character(LEN=1)   :: cgrddef
  integer, parameter :: maxbglvl = 52 
  integer, dimension(255) :: iuarr = 0
  integer :: iplvl
  integer :: nlvl
  integer :: ierr
  integer :: nunit1 = 12
  integer :: debug_level = 1000
  integer :: grib_version
  integer :: vtable_columns
  integer :: istatus
  integer :: nx, ny, nz
  real, dimension(maxbglvl) :: plvl
  real    :: level
  real    :: dx, dy  ! Required by Laps, in meters.
  real    :: sw1, sw2, ne1, ne2
  real    :: realI, realJ, diff_lon
  real    :: stdlon, truelat1, truelat2
  type(proj_info) :: proj  ! Declared via "USE map_utils" 
  logical :: val_std = .false.
  logical :: cross_dateline 
  logical :: l_full_run

! -----------------
! Determine GRIB Edition number
  grib_version=0
  call edition_num(nunit1, gribflnm, grib_version, ierr)
  call mprintf((ierr.eq.3),ERROR,"GRIB file problem")
  if (grib_version.eq.2) then
     vtable_columns=11 

#if defined (USE_PNG) && (USE_JPEG2000)
!     call mprintf(val_inf,INFORM, &
!        "Linked in png and jpeg libraries for Grib Edition 2")
#else
     val_std = .true.
     call mprintf(val_std,STDOUT,"WARNING - Grib2 data detected. Code, as built, expected Grib1")
     call mprintf(val_std,STDOUT,"\t- Three external libs are necessary Grib2 decoding version to be compiled: ")
     call mprintf(val_std,STDOUT,"\t- -ljasper, -lpng -lz. If these libs are NOT found with the configuration")
     call mprintf(val_std,STDOUT,"\t- scripts then the build scripts CANNOT build Grib2. See the README file.")
     call mprintf(val_std,STDOUT,"\t*** stopping readgrib ***")
     val_std = .false.
     stop
#endif
  else
     vtable_columns=7 
  endif


! Read Vtable.XXX file; its information in sent to the module 'table'
! and 'maxvar' is incremented with each variable found.

  call parse_table(vtablefn, debug_level,vtable_columns)

  if (.not. l_full_run) then  
      write(6,*)' returning with partial run'
      return
  endif

        nlvl = 0
        plvl = -999.
        call clear_storage

           ! If we need to read a new grib record, then read one.

              if (grib_version.ne.2) then 
                 ! Read one record at a time from GRIB1 (and older Editions) 
                 call rd_grib1(nunit1, gribflnm, level, field, &
                      hdate, ierr, iuarr, debug_level)

                call cclose(iuarr(nunit1), debug_level, ierr)
                iuarr(nunit1) = 0

              else 

                 ! Read one file of records from GRIB2.
                 call rd_grib2(nunit1, gribflnm, hdate, &
                      grib_version, ierr, debug_level)
                 FIELD='NULL'

              endif


        if (map%igrid.eq.0) then ! Lat/Lon grid aka Cylindrical Equidistant
              gproj='LL'
        elseif (map%igrid.eq.5) then ! Polar-Stereographic Grid.
              gproj='PS'
        elseif (map%igrid.eq.3) then ! Lambert Conformal Grid
              gproj='LC'
        elseif(map%igrid.eq.4) then ! Gaussian Grid (we will call it lat/lon)
              gproj='GA'
        endif 

        nx=map%nx
        ny=map%ny
        nz=maxbglvl
        write(6,*)' nz is being set to ',maxbglvl
        dx=map%dx
        dy=abs(map%dx)

        if (map%lat2.gt.10000.) map%lat2=map%lat2/1000000
        if (map%lon2.gt.10000.) map%lon2=map%lon2/1000000

        ! Cross dateline test. (Needs more thought.)
        cross_dateline=.false.
        if(map%lon2.gt.map%lon1) diff_lon=map%lon2-map%lon1
        if(map%lon1.gt.map%lon2) diff_lon=map%lon1-map%lon2
        if (diff_lon.gt.350) cross_dateline=.true.

        if (map%lov .gt.180.) map%lov =map%lov -360
        if (map%lon1.gt.180.) map%lon1=map%lon1-360
        if (map%lon2.gt.180.) map%lon2=map%lon2-360

        stdlon=map%lov
        truelat1=map%truelat1
        truelat2=map%truelat2
        sw1=map%lat1
        sw2=map%lon1
        ne1=map%lat2
        ne2=map%lon2
        cgrddef='N'

        if (gproj.eq.'LC' .or. gproj.eq.'PS') then 
           dx = dx * 1000. ! meters at 25.0 N, for example
           dy = dy * 1000. 
           call map_set(PROJ_LC,sw1,sw2,dx,stdlon,truelat1,truelat2, &
                   nx,ny,proj)
           realI=nx
           realJ=ny
           call ij_to_latlon(proj,realI,realJ,ne1,ne2);
          !call latlon_to_ij(proj,sw1,sw2,realI,realJ); ! call args can be inversed.
        elseif (gproj.eq.'LL') then 
           ! This seems unconventional but follows the pattern used in get_bkgd_mdl_info.f
           stdlon=map%lon1
           truelat1=map%lat1
           truelat2=0.
           if (sw1.gt.ne1) then 
              sw1=map%lat2
              ne1=map%lat1
           endif
        endif

        write(*, *) "---------- "
        write(*, *) "proj ", gproj
        write(*, *) "truelat1, truelat2, stdlon ", truelat1, truelat2, stdlon 
        write(*, *) "dx, dy ", dx, dy
        write(*, *) "nx, ny, nz ", nx, ny, nz
        write(*, *) "sw corner ", sw1, sw2
        write(*, *) "ne corner ", ne1, ne2
        write(*, *) "longitude spans ", diff_lon, " degrees; global ", cross_dateline
        write(*, *) "---------- "


        istatus=1
        !if (ierr.eq.1) print*, "iERR ", ierr !istatus=0

  return
 end subroutine degrib_nav

!-------------------------------------------------------------------

subroutine degrib_data(gribflnm, nx, ny, nz, &
         prbght, htbg, tpbg, shbg, uwbg, vwbg, wwbg, &
         htbg_sfc, tpbg_sfc, shbg_sfc, uwbg_sfc, vwbg_sfc, &
         tdbg_sfc, t_at_sfc, prbg_sfc, mslpbg, pcpbg, crefbg, pwatbg, cwatbg, istatus)

  use table
  use gridinfo
  use storage_module
  use filelist
  use datarray
  use module_debug
  use stringutil

  implicit none

  character(LEN=150) :: gribflnm
  character(LEN=19)  :: hdate
  character(LEN=9)   :: field
  integer, parameter :: maxbglvl = 52  
  integer, dimension(255) :: iuarr = 0
  integer :: nunit1 = 12
  integer :: debug_level = 1000
  integer :: iplvl
  integer :: nlvl
  integer :: itime
  integer :: ntimes
  integer :: ierr
  integer :: grib_version
  integer :: istatus, i, j, k, i4time, ishow_timer 
  integer :: nx, ny, nz, nzsh
  integer :: icn3d, icm(maxbglvl)
  real, dimension(maxbglvl) :: plvl
  real    :: level
  real    :: t_ref, rfill, prsfc, qsfc
  real    :: make_ssh, make_td
  real    :: sumtot, shsum(maxbglvl), shavg, r_bogus_sh
  real    :: it,xe,mrsat,esat
  logical :: readit
  logical :: val_std = .false.
  logical :: isthere_shsfc
  logical :: isthere_tdsfc
  logical :: isthere_prsfc

! *** sfc background arrays.
  real :: prbg_sfc(nx,ny)
  real :: uwbg_sfc(nx,ny)
  real :: vwbg_sfc(nx,ny)
  real :: shbg_sfc(nx,ny)
  real :: tdbg_sfc(nx,ny)
  real :: tpbg_sfc(nx,ny) ! 2m temp
  real :: t_at_sfc(nx,ny) ! skin temp
  real :: htbg_sfc(nx,ny)
  real :: mslpbg(nx,ny)
  real :: pcpbg(nx,ny)
  real :: crefbg(nx,ny)
  real :: pwatbg(nx,ny)
  real :: cwatbg(nx,ny)

! *** 3D background arrays.
  real :: prbght(nx,ny,nz)
  real :: htbg(nx,ny,nz)
  real :: tpbg(nx,ny,nz)
  real :: shbg(nx,ny,nz)
  real :: uwbg(nx,ny,nz)
  real :: vwbg(nx,ny,nz)
  real :: wwbg(nx,ny,nz)

  write(6,*)' Start degrib_data: dims are ',nx,ny,nz

! -----------------
! Determine GRIB Edition number
  grib_version=0
  call edition_num(nunit1, gribflnm, grib_version, ierr)
  call mprintf((ierr.eq.3),ERROR,"GRIB file problem")

! -----------------
! The "Vtable" file was read in degrib_nav (via 'call parse_table'); 
! its information in sent to the module "table".

     ! Set READIT to .TRUE., meaning that we have not read any records yet 
     ! from the file GRIBFLNM.
     readit = .TRUE.  ! i.e., "Yes, we want to read a record."

! LOOP1 reads through the file GRIBFLNM, exiting under two conditions:
!        1) We have hit the end-of-file
!        2) We have read past the ending date HEND.
!
! Condition 2 assumes that the data in file GRIBFLNM are ordered in time.

     LOOP1 : DO
        ! At the beginning of LOOP1, we are at a new time period.
        ! Clear the storage arrays and associated level information.
        nlvl = 0
        plvl = -999.
        call clear_storage

! LOOP0 reads through the file GRIBFLNM, looking for data of the current 
! date.  It exits under the following conditions.
!          1) We have hit the end-of-file
!          2) The GRIBFLNM variable has been updated to a nonexistent file.
!          3) We start reading a new date and the data are assumed to be 
!             ordered by date.
!          4) We have a valid record and the data are not assumed to be 
!             ordered by date.

        LOOP0 : DO

           ! If we need to read a new grib record, then read one.
           if (READIT) then

              if (grib_version.ne.2) then 
                 ! Read one record at a time from GRIB1 (and older Editions) 
                 call rd_grib1(nunit1, gribflnm, level, field, &
                      hdate, ierr, iuarr, debug_level)
              else 

                 ! Read one file of records from GRIB2.
                 call rd_grib2(nunit1, gribflnm, hdate, &
                      grib_version, ierr, debug_level)
                 FIELD='NULL'

              endif

              if (ierr.eq.1) then 
                 ! We have hit the end of a file.  Exit LOOP0 then exit LOOP1
                 exit LOOP0
              endif

           endif

           if (FIELD.EQ.'NULL') then
              ! The data read does not match any fields requested
              ! in the Vtable.  So cycle LOOP0 to read the next GRIB record.
              READIT = .TRUE.
              cycle LOOP0
           endif

! When we have reached this point, we have a data array ARRAY which has 
! some data we want to save, with field name FIELD at pressure level 
! LEVEL (Pa).  Dimensions of this data are map%nx and map%ny.  Put
! this data into storage.

           if (((field == "SST").or.(field == "SKINTEMP")) .and. &
                (level /= 200100.)) level = 200100.
           iplvl = int(level)
           call put_storage(iplvl,field, &
                reshape(rdatarray(1:map%nx*map%ny),(/map%nx, map%ny/)),&
                map%nx,map%ny)
           ! Since we processed the record we just read, we set
           ! READIT to .TRUE. so that LOOP0 will read the next record.
           READIT = .TRUE.

        enddo LOOP0

        i4time = ishow_timer()
!-----
        write(6,*)' call get_lapsbg: dims are ',nx,ny,nz
        call get_lapsbg(nlvl, maxbglvl, plvl, debug_level, nx, ny, nz, &
         prbght, htbg, tpbg, shbg, uwbg, vwbg, wwbg, &
         htbg_sfc, tpbg_sfc, shbg_sfc, uwbg_sfc, vwbg_sfc, &
         tdbg_sfc, t_at_sfc, prbg_sfc, mslpbg, pcpbg, crefbg, pwatbg, cwatbg, istatus)

!-----

        i4time = ishow_timer()

! When we have reached this point, we have either hit the end of file, or 
! hit the end of the current date.  Either way, we want to output
! the data found for this date.  This current date is in HSAVE.
! However, if (HSAVE == 0000-00-00_00:00:00), no output is necessary,
! because that 0000 date is the flag for the very opening of a file.

        ! If we hit the end-of-file, its time to exit LOOP1 so we can
        ! increment the GRIBFLNM to read the next file.
        if (ierr.eq.1) exit LOOP1

     enddo LOOP1

     if (grib_version.ne.2) then
        call cclose(iuarr(nunit1), debug_level, ierr)
        iuarr(nunit1) = 0
     endif 

! ------------- qcmodel sh ----------------

     !nzsh=nz-5
     nzsh=nz
     call qcmodel_sh(nx,ny,1,tdbg_sfc)  !td_sfc is actually = RH for GFS/AVN.
     call qcmodel_sh(nx,ny,nzsh,shbg)   !sh is actually = RH for GFS/AVN.

! ------------- convert rh to q ----------------

! Td or rh liq/ice phase temp thresh
! ---------------------------------

     print*,'convert rh to q - 3D'

     t_ref=-132.0
     rfill = 1e37 ! -99999.

     if(nzsh .gt. maxbglvl)then
        write(6,*)' ERROR: nzsh > maxbglvl ',nzsh,maxbglvl 
        write(6,*)' Increase dimension of maxbglvl'                  
     endif

     do k=1,nzsh
        icm(k)=0
        shsum(k)=0.0

        do j=1,ny
        do i=1,nx

          if(shbg(i,j,k) .gt. 1e-7 .and. shbg(i,j,k) .le.101.00)then
            if(shbg(i,j,k).gt.100.) shbg(i,j,k)=100.
             shbg(i,j,k)=make_ssh(prbght(i,j,k),tpbg(i,j,k)-273.15, &
                                   shbg(i,j,k)/100.,t_ref)*0.001
             shsum(k)=shsum(k)+shbg(i,j,k)
          else
            icm(k)=icm(k)+1
            shbg(i,j,k)=rfill
          endif
 
        enddo
        enddo
 
     enddo

!
!---------- check for and fill missing sh with computed avg for that level
!
     sumtot = 0.0
     icn3d = 0
     r_bogus_sh = 1.0e-6

     do k=1,nz
        if(icm(k).gt.0.and.icm(k).lt.nx*ny)then
           shavg=shsum(k)/(nx*ny-icm(k))
           do j=1,ny
           do i=1,nx
              if(shbg(i,j,k).eq.rfill)shbg(i,j,k)=shavg
           enddo
           enddo
           sumtot=shsum(k)+sumtot
           icn3d=icm(k)+icn3d
        elseif(icm(k).eq.nx*ny)then
           do j=1,ny
           do i=1,nx
              if(shbg(i,j,k).eq.rfill)shbg(i,j,k)=r_bogus_sh
           enddo
           enddo
        endif
     enddo

     if(sumtot.gt.0.0.and.icn3d.ne.nx*ny*nz)then
        print*,'Missing background 3D moisture filled with average of good points'
        print*,'#/% 3D points filled: ',icn3d,icn3d/(nx*ny*nz)
     endif

 
! ------------- convert rh (or sh, if avail) to Td ----------------
! ------------- also calculate sh if unavailable   ----------------

     i4time = ishow_timer()

!    Is this test valid for just GRIB1 or for certain models?
     isthere_shsfc = (is_there(200101,'SH_SFC') .OR. is_there(200100,'SH_SFC')&
                                                .OR. is_there(200102,'SH_SFC'))
     if (isthere_shsfc .eqv. .true.) then
         print*,'SH_SFC is available'
     else
         print*,'SH_SFC is unavailable according to "is_there" test'
     endif

!    Is this test valid for just GRIB1 or for certain models?
     isthere_prsfc = (is_there(200101,'PSFC') .OR. is_there(200100,'PSFC')&
                                              .OR. is_there(200102,'PSFC'))
     if (isthere_prsfc .eqv. .true.) then
         print*,'PSFC is available'
     else
         print*,'PSFC is unavailable according to "is_there" test'
     endif

!    Is this test valid for just GRIB1 or for certain models?
     isthere_tdsfc = (is_there(200101,'TD_SFC') .OR. is_there(200100,'TD_SFC')&
                                                .OR. is_there(200102,'TD_SFC'))

     if(isthere_tdsfc .eqv. .true.)then
       write(6,*)' TD_SFC is available, use it directly...'

     else
       print*,'TD_SFC is unavailable according to "is_there" test'
       print*,'convert rh (or sh, if avail) to Td - sfc'

       k=1
       do j=1,ny
       do i=1,nx

        if(tdbg_sfc(i,j).gt.0.0 .and. tdbg_sfc(i,j).lt.100.001)then
           prsfc=prbg_sfc(i,j)/100.
           if (isthere_shsfc) then ! there is SH_SFC                     
              qsfc=shbg_sfc(i,j)
           else                    ! there is no SH_SFC, using RH_SFC
              qsfc=make_ssh(prsfc,tpbg_sfc(i,j)-273.15,tdbg_sfc(i,j)/100.,t_ref)
              shbg_sfc(i,j) = qsfc
           endif
           tdbg_sfc(i,j)=make_td(prsfc,tpbg_sfc(i,j)-273.15,qsfc,t_ref)+273.15
        else
           !Td is rfill
           tdbg_sfc(i,j)=rfill
        endif

       enddo ! i
       enddo ! j

     endif ! isthere_tdsfc

     write(6,*)' tdbg_sfc range = ',minval(tdbg_sfc),maxval(tdbg_sfc)
     write(6,*)' shbg_sfc range = ',minval(shbg_sfc),maxval(shbg_sfc)
     write(6,*)' mslpbg   range = ',minval(mslpbg)  ,maxval(mslpbg)
     write(6,*)' pcpbg    range = ',minval(pcpbg)   ,maxval(pcpbg)    

! ------------- end convert data ----------------
 
     write(6,*)' returning from degrib_data...'
 
 return

end subroutine degrib_data
