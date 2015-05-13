subroutine output(hdate, nlvl, maxlvl, plvl, interval, iflag, out_format, prefix, debug_level)
!                                                                             !
!*****************************************************************************!
!  Write output to a file.
!                                                                             !
!    hdate:
!    nlvl:  
!    maxlvl:
!    plvl:
!    interval:
!    iflag:  1 = output for ingest into rrpr ; 2 = output for hinterp
!    out_format: requested output format (WPS, SI, or MM5)
!                                                                             !
!*****************************************************************************!

  use table
  use gridinfo
  use storage_module
  use filelist
  use module_debug
  use stringutil

  implicit none

  character(LEN=19) :: hdate
  character(LEN=24) :: hdate_output
  character(LEN=3)  :: out_format
  character(LEN=256)  :: prefix
  integer :: iunit = 13

  real, pointer, dimension(:,:) :: scr2d

  integer :: maxlvl
  integer nlvl, debug_level
  real , dimension(maxlvl) :: plvl
  character (LEN=9) :: field
  real :: level
  integer :: sunit = 14
  integer :: interval
  integer :: iflag
! Local Miscellaneous
  integer :: k, n, m, ilev
  integer :: ii, jj
  real :: maxv, minv
  real :: xplv
  real :: xfcst = 0.
  character (LEN=25) :: units
  character (LEN=46) :: Desc
  logical lopen

! DATELEN:  length of date strings to use for our output file names.
  integer :: datelen

  write(6,*)' Subroutine output (output.F90) - debug_level is ',debug_level

! Decide the length of date strings to use for output file names.  
! DATELEN is 13 for hours, 16 for minutes, and 19 for seconds.
  if (mod(interval,3600) == 0) then
     datelen = 13
  elseif (mod(interval,60) == 0) then
     datelen = 16
  else
     datelen = 19
  endif
 
  call get_plvls(plvl, maxlvl, nlvl)

  if(nlvl .gt. maxlvl)then
      write(6,*)' WARNING: get_plvls returned nlvl > maxlvl ',nlvl,maxlvl
  endif

  if ( debug_level .ge. 0 ) then
  write(*,119) hdate(1:10), hdate(12:19)
119 format(/,79('#'),//,'Inventory for date = ', A10,1x,A8,/)

  write(*,advance='NO', fmt='("PRES", 2x)')
  WRTLOOP : do n = 1, maxvar
     do k = 1, n-1
        if (namvar(k).eq.namvar(n)) cycle WRTLOOP
     enddo
     write(*,advance='NO', fmt='(1x,A8)') namvar(n)
  enddo WRTLOOP
  write(*,advance='YES', fmt='(1x)')

  write(*,FMT='(79("-"))')
  end if
  KLOOP : do k = 1, nlvl
     if ((iflag.eq.2).and.(plvl(k).gt.200100) .and. (plvl(k).lt.200200)) then
        cycle KLOOP
     endif
     ilev = nint(plvl(k))
     if ( debug_level .ge. 0 ) then
     write(*, advance='NO', FMT='(F6.1)') plvl(k)/100.
     end if
     MLOOP : do m = 1, maxvar
        do n = 1, m-1
           if (namvar(m).eq.namvar(n)) cycle MLOOP
        enddo
        if ( debug_level .ge. 0 ) then
        if (is_there(ilev,namvar(m))) then
           write(*, advance='NO', FMT='("  X      ")')
        else
	   if ( plvl(k).gt.200000 ) then
             write(*, advance='NO', FMT='("  O      ")')
	   else
             write(*, advance='NO', FMT='("         ")')
	   endif
        endif
        endif
     enddo MLOOP
     if ( debug_level .ge. 0 ) then
     write(*,advance='YES', fmt='(1x)')
     endif
  enddo KLOOP
  if ( debug_level .ge. 0 ) then
  write(*,FMT='(79("-"))')
  endif
  
  if (iflag.eq.1) then
     if (nfiles.eq.0) then
        open(iunit, file=trim(get_path(prefix))//'PFILE:'//HDATE(1:datelen), form='unformatted', &
             position='REWIND')
        nfiles = nfiles + 1
        filedates(nfiles)(1:datelen) = hdate(1:datelen)
     else
        DOFILES : do k = 1, nfiles
           if (hdate(1:datelen).eq.filedates(k)(1:datelen)) then
              open(iunit, file=trim(get_path(prefix))//'PFILE:'//HDATE(1:datelen), form='unformatted',&
                   position='APPEND')
           endif
        enddo DOFILES
        inquire (iunit, OPENED=LOPEN)
        if (.not. LOPEN) then
           open(iunit, file=trim(get_path(prefix))//'PFILE:'//HDATE(1:datelen), form='unformatted', &
                position='REWIND')
           nfiles = nfiles + 1
           filedates(nfiles)(1:datelen) = hdate(1:datelen)
        endif
     endif
  else if (iflag.eq.2) then
     open(iunit, file=trim(prefix)//':'//HDATE(1:datelen), form='unformatted', &
          position='REWIND')
  endif

  NLOOP : do n = 1, nlvl

     OUTLOOP : do m = 1, maxvar
        field = namvar(m)
        do k = 1, m-1
           if (field.eq.namvar(k)) cycle OUTLOOP
        enddo
        level = plvl(n)
        if ((iflag.eq.2).and.(level.gt.200100) .and. (level.lt.200200)) then
           cycle NLOOP
        endif
        ilev = nint(level)
        desc = ddesc(m)
        if (iflag.eq.2) then
           if (desc.eq.' ') cycle OUTLOOP
        endif
        units = dunits(m)
        if ((iflag.eq.1).or.(iflag.eq.2.and.desc(1:1).ne.' ')) then
          if (is_there(ilev,field)) then 
            call get_dims(ilev, field)

            call refr_storage(ilev, field, scr2d, map%nx, map%ny)

            if (out_format(1:2) .eq. 'SI') then
              write(iunit) 4
              hdate_output = hdate
              write (iunit) hdate_output, xfcst, map%source, field, units, &
                   Desc, level, map%nx, map%ny, map%igrid
              if (map%igrid.eq.3) then ! lamcon
                 write (iunit) map%startloc, map%lat1, map%lon1, map%dx, map%dy, &
                      map%lov, map%truelat1, map%truelat2
              elseif (map%igrid.eq.5) then ! Polar Stereographic
                 write (iunit) map%startloc, map%lat1, map%lon1, map%dx, map%dy, &
                      map%lov, map%truelat1
              elseif (map%igrid.eq.0 .or. map%igrid.eq.4)then ! lat/lon
                 write (iunit) map%startloc, map%lat1, map%lon1, map%dy, map%dx
              elseif (map%igrid.eq.1)then ! Mercator
                 write (iunit) map%startloc, map%lat1, map%lon1, map%dy, map%dx, &
                      map%truelat1
              else
                 call mprintf(.true.,ERROR, &
                "Unrecognized map%%igrid: %i in subroutine output 1",i1=map%igrid)
              endif
              write (iunit) scr2d
	    else if (out_format(1:2) .eq. 'WP') then   
                call mprintf(.true.,DEBUG, &
         "writing in WPS format  iunit = %i, map%%igrid = %i",i1=iunit,i2=map%igrid)
              write(iunit) 5
              hdate_output = hdate
              write (iunit) hdate_output, xfcst, map%source, field, units, &
                   Desc, level, map%nx, map%ny, map%igrid
              if (map%igrid.eq.3) then ! lamcon
                 write (iunit) map%startloc, map%lat1, map%lon1, map%dx, map%dy, &
                      map%lov, map%truelat1, map%truelat2, map%r_earth
              elseif (map%igrid.eq.5) then ! Polar Stereographic
                 write (iunit) map%startloc, map%lat1, map%lon1, map%dx, map%dy, &
                      map%lov, map%truelat1, map%r_earth
              elseif (map%igrid.eq.0 .or. map%igrid.eq.4)then ! lat/lon
                 write (iunit) map%startloc, map%lat1, map%lon1, map%dy, map%dx, &
		      map%r_earth
              elseif (map%igrid.eq.1)then ! Mercator
                 write (iunit) map%startloc, map%lat1, map%lon1, map%dy, map%dx, &
                      map%truelat1, map%r_earth
              else
                 call mprintf(.true.,ERROR, &
                "Unrecognized map%%igrid: %i in subroutine output 1",i1=map%igrid)
              endif
	      write (iunit) map%grid_wind
              write (iunit) scr2d
	    else if (out_format(1:2) .eq. 'MM') then
              write(iunit) 3
              hdate_output = hdate
              write (iunit) hdate_output, xfcst, field, units, Desc, level,&
                   map%nx, map%ny, map%igrid
              if (map%igrid.eq.3) then ! lamcon
                 write (iunit) map%lat1, map%lon1, map%dx, map%dy, map%lov, &
                      map%truelat1, map%truelat2
              elseif (map%igrid.eq.5) then ! Polar Stereographic
                 write (iunit) map%lat1, map%lon1, map%dx, map%dy, map%lov, &
                      map%truelat1
              elseif (map%igrid.eq.0 .or. map%igrid.eq.4)then ! lat/lon
                 write (iunit) map%lat1, map%lon1, map%dy, map%dx
              elseif (map%igrid.eq.1)then ! Mercator
                 write (iunit) map%lat1, map%lon1, map%dy, map%dx, map%truelat1
              else
                 call mprintf(.true.,ERROR, &
                "Unrecognized map%%igrid: %i in subroutine output 1",i1=map%igrid)
              endif
              write (iunit) scr2d
	    endif
              if ( debug_level .gt. 100 ) then
	        call mprintf(.true.,DEBUG, &
	        "hdate = %s,  xfcst = %f ",s1=hdate_output,f1=xfcst)
	        call mprintf(.true.,DEBUG, &
           "map%%source = %s, field = %s, units = %s",s1=map%source,s2=field,s3=units)
	        call mprintf(.true.,DEBUG, &
	            "Desc = %s, level = %f",s1=Desc,f1=level)
	        call mprintf(.true.,DEBUG, &
	            "map%%nx = %i, map%%ny = %i",i1=map%nx,i2=map%ny)
              else if ( debug_level .gt. 0 ) then
	        call mprintf(.true.,STDOUT, &
	      " field = %s, level = %f",s1=field,f1=level)
	        call mprintf(.true.,LOGFILE, &
	      " field = %s, level = %f",s1=field,f1=level)
              end if
              if ( debug_level .gt. 100 ) then
	      maxv = -99999.
	      minv = 999999.
	      do jj = 1, map%ny
	      do ii = 1, map%nx
	        if (scr2d(ii,jj) .gt. maxv) maxv = scr2d(ii,jj)
	        if (scr2d(ii,jj) .lt. minv) minv = scr2d(ii,jj)
	      enddo
	      enddo
	      call mprintf(.true.,DEBUG, &
	         "max value = %f , min value = %f",f1=maxv,f2=minv)
              end if

              nullify(scr2d)

           endif
        endif
     enddo OUTLOOP
  enddo NLOOP

  close(iunit)

end subroutine output


subroutine get_lapsbg(nlvl, maxlvl, plvl, debug_level, nx, ny, nz&
         ,prbght, htbg, tpbg, shbg, uwbg, vwbg, wwbg, htbg_sfc, tpbg_sfc, shbg_sfc& 
         ,uwbg_sfc, vwbg_sfc, tdbg_sfc, t_at_sfc, prbg_sfc, mslpbg, pcpbg, crefbg, pwatbg, istatus)
!                                                                             !
!*****************************************************************************!
!  Write output to a file.
!                                                                             
!    nlvl:  
!    maxlvl:
!    plvl:
!    interval:
!                                                                             !
!*****************************************************************************!

  use table
  use gridinfo
  use storage_module
  use filelist
  use module_debug
  use stringutil

  implicit none

  real, pointer, dimension(:,:) :: scr2d
  integer :: maxlvl
  integer nlvl, debug_level
  real , dimension(maxlvl) :: plvl
  character (LEN=9) :: field
  real :: level
  integer :: sunit = 14
  integer :: interval
! Local Miscellaneous
  integer :: k, n, m, ilev
  integer :: ii, jj
  real :: maxv, minv
  real :: xplv
  real :: xfcst = 0.
  character (LEN=25) :: units
  character (LEN=46) :: Desc
  logical lopen
  integer :: idx
  integer :: istatus 
  integer :: nx, ny, nz 

! *** Background model grid data.
!
!
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
     real :: pratebg(nx,ny)

! *** 3D background arrays.

     real :: prbght(nx,ny,nz)
     real :: htbg(nx,ny,nz)
     real :: tpbg(nx,ny,nz)
     real :: shbg(nx,ny,nz)
     real :: uwbg(nx,ny,nz)
     real :: vwbg(nx,ny,nz)
     real :: wwbg(nx,ny,nz)

     include 'constants.inc' ! for grav

     write(6,*)' Subroutine get_lapsbg (output.F90)...'

     call get_plvls(plvl, maxlvl, nlvl)

     if(nlvl .gt. maxlvl)then
         write(6,*)' WARNING: get_plvls returned nlvl > maxlvl ',nlvl,maxlvl
     else
         write(6,*)' nlvl/maxlvl = ',nlvl,maxlvl
         write(6,*)' plvl = ',plvl
     endif

     NLOOP : do m = 1, maxvar
     idx=0
     OUTLOOP : do n = 1, nlvl
        field = namvar(m)
        level = plvl(n)
        ilev = nint(level)
        desc = ddesc(m)
        if (desc.eq.' ') cycle OUTLOOP
        if (level.lt.200000) idx=idx+1 
        units = dunits(m)
        if (desc(1:1).ne.' ') then

          if(ilev .ge. 200000)write(6,*)'ilev/desc ',ilev,desc

          if (is_there(ilev,field)) then 
              call get_dims(ilev, field)
              call refr_storage(ilev, field, scr2d, map%nx, map%ny)
              if(ilev .ge. 200000)then
                write(6,11)m, idx, field, units(1:5), Desc(1:20), ilev
11              format('  isthere',2i5,a,a,a,i8)
              endif

              if(map%nx .gt. nx .OR. map%ny .gt. ny &
                                .OR. idx    .gt. nz)then      
                  write(6,*)' WARNING: get_lapsbg array bounds issue'
                  write(6,*)' map%nx/map%ny/nx/ny = ' & 
                             ,map%nx,map%ny,nx,ny      
                  write(6,*)' idx/nz = ',idx,nz                            
              endif

              if (field.eq.'HGT') then
                 htbg(:,:,idx) = scr2d
              elseif (field.eq.'GEOPT') then
                 htbg(:,:,idx) = scr2d / grav
              elseif (field.eq.'TT') then
                 tpbg(:,:,idx) = scr2d
              elseif (field.eq.'RH') then
                 shbg(:,:,idx) = scr2d
              elseif (field.eq.'UU') then
                 uwbg(:,:,idx) = scr2d
              elseif (field.eq.'VV') then
                 vwbg(:,:,idx) = scr2d
              elseif (field.eq.'VVEL') then
                 wwbg(:,:,idx) = scr2d
              elseif (field.eq.'HGT_SFC') then
                 htbg_sfc = scr2d
              elseif (field.eq.'GEOPT_SFC') then
                 htbg_sfc = scr2d
              elseif (field.eq.'TT_SKIN') then
                 t_at_sfc = scr2d
              elseif (field.eq.'TD_SFC') then
                 tdbg_sfc = scr2d 
                 write(6,*)' Filling tdbg_sfc with TD_SFC'
              elseif (field.eq.'RH_SFC') then
                 shbg_sfc = scr2d
                 tdbg_sfc = scr2d  !See bgdata/readdgprep.f, line 183.
     		                   !write(*, *) "RH_SFC", scr2d(3,30)
                 write(6,*)' Filling tdbg_sfc with RH_SFC'
                 write(6,*)' Filling shbg_sfc with RH_SFC'
              elseif (field.eq.'SH_SFC') then
                 shbg_sfc = scr2d
                 write(6,*)' Filling shbg_sfc with SH_SFC'
              elseif (field.eq.'UU_SFC') then
                 uwbg_sfc = scr2d
              elseif (field.eq.'VV_SFC') then
                 vwbg_sfc = scr2d
              elseif (field.eq.'TT_SFC') then
                 tpbg_sfc = scr2d
              elseif (field.eq.'PSFC') then
                 prbg_sfc = scr2d
                 write(6,*)' Filling prbg_sfc with PSFC'
              elseif (field.eq.'PMSL') then
                 mslpbg = scr2d
              elseif (field.eq.'APCP') then
                 pcpbg = scr2d
                 write(6,*)' Filling pcpbg with APCP'
              elseif (field.eq.'REFC') then
                 crefbg = scr2d
                 write(6,*)' Filling crefbg with REFC'
              elseif (field.eq.'PWAT') then
                 pwatbg = scr2d
              elseif (field.eq.'PRATE') then
                 pratebg = scr2d
                 write(6,*)' Filling pratebg with PRATE'
              else
                 write(6,*)' field not filled: ',trim(field)
              endif
              nullify(scr2d)

          elseif  (ilev .ge. 200000)  then
              write(6,*)' field is not there: ',ilev,trim(field)

          endif ! if parm exists

          if ((level.gt.200100 .and. level.lt.200200)) then
              cycle OUTLOOP
          endif

       else
          write(6,*)' Description is null ',trim(field)

       endif ! if desc=null
     enddo OUTLOOP
  enddo NLOOP

  idx=0
  do n = 1, nlvl
    level = plvl(n)
    if (level.lt.200000) then
      idx=idx+1
      prbght(:,:,idx) = level/100
    endif
  enddo 

! ------------- fill RH levels if needed ----------------
!
! If upper-air RH is missing, see if we can interpolate from surrounding levels.
! This is a simple vertical interpolation, linear in pressure.
! Currently, this simply fills in one missing level between two present levels. 
! May expand this in the future to fill in additional levels.  May also expand 
! this in the future to vertically interpolate other variables.
!

!     do k = 2, nlvl-1, 1
!        if (plvl(k-1) .lt. 200000.) then
!           if ( (.not. is_there(nint(plvl(k)),'RH')) .and. &
!                ( is_there(nint(plvl(k-1)), 'RH')) .and.&
!                ( is_there(nint(plvl(k+1)), 'RH')) ) then
!              call get_dims(nint(plvl(k+1)), 'RH')
!              call vntrp(plvl, maxlvl, k, "RH      ", map%nx, map%ny)
!              print*,'Missing background 3D moisture interpolate from surrounding levels'
!           endif
!        endif
!     enddo

!----------

!     write(*, *) "OUTPUT htbg(3,30,1)", htbg(3,30,1)
!     write(*, *) "OUTPUT tpbg(3,30,1)", tpbg(3,30,1)
!     write(*, *) "OUTPUT wwbg(3,30,1)", wwbg(3,30,1)
!     write(*, *) "OUTPUT tdbg_sfc(3,30)", tdbg_sfc(3,30)
!     write(*, *) "OUTPUT shbg_sfc(3,30)", shbg_sfc(3,30)
!     do jj = 1, 6 
!        write(*, *) "OUTPUT shbg(3,30,",jj, shbg(3,30,jj)
!     enddo
!     do jj = 1, 26 
!        write(*, *) "OUTPUT pcpbg(30,",jj, pcpbg(30,jj)
!     enddo

      write(*, *) "End of get_lapsbg..."

 return
end subroutine get_lapsbg
