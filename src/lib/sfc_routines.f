cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
c
c
	subroutine write_surface_obs(btime,outfile,n_obs_g,
     &    n_obs_b,wmoid,stations,provider,wx,reptype,autostntype,
     &    store_1,store_2,store_3,store_4,store_5,store_6,store_7,
     &    store_2ea,store_3ea,store_4ea,store_5ea,store_6ea,
     &    store_cldamt,store_cldht,maxsta,jstatus)
c
c*****************************************************************************
c
c	Routine to write the LAPS surface data file.   The data is passed
c       to this routine via the 'store' array.
c
c	Changes:
c		P. Stamus  03-27-98  Original version (from old format
c                                      version of write_surface_obs).
c                          05-01-98  Added soil moisture to 'store_5' & '_5ea'
c                          09-04-98  Final adjustments for operational use.
c
c
c               J. Edwards 09-16-98  moved all format definitions to 
c                                    src/include/lso_formats.inc
c                                    changed 909 definition to allow for 
c                                    missing data
c
c*****************************************************************************
c
	real store_1(maxsta,4), 
     &         store_2(maxsta,3), store_2ea(maxsta,3),
     &         store_3(maxsta,4), store_3ea(maxsta,2),
     &         store_4(maxsta,5), store_4ea(maxsta,2),
     &         store_5(maxsta,4), store_5ea(maxsta,4),
     &         store_6(maxsta,5), store_6ea(maxsta,2),
     &         store_7(maxsta,3),
     &         store_cldht(maxsta,5)
c
	integer jstatus, wmoid(maxsta)
c
	character btime*24, outfile*(*), 
     &         stations(maxsta)*20, provider(maxsta)*11,
     &         wx(maxsta)*25,reptype(maxsta)*6, 
     &         autostntype(maxsta)*6,store_cldamt(maxsta,5)*4
c
c
c.....	Write the file.
c
	open(11,file=outfile,status='unknown')
c
c.....	Write the header.
c
	write(11,900) btime,		! time
     &               n_obs_g,		! # of obs in the laps grid
     &               n_obs_b		! # of obs in the box
c
c.....	Write the station data.
c
	do k=1,n_obs_b
c
           call filter_string(stations(k))
           call filter_string(provider(k))

	   write(11,901) stations(k),           !station id
     &                   wmoid(k),              !WMO id number
     &                   provider(k),           !data provider
     &                   (store_1(k,i),i=1,3),  !lat, lon, elev
     &                   nint(store_1(k,4))     !obs time
c
          call filter_string(reptype(k))
          call filter_string(autostntype(k))
          call filter_string(wx(k))

	  write(11,903)  reptype(k),            !station report type
     &                   autostntype(k),        !station type (manual/auto)
     &                   wx(k)                  !present weather
c
	  write(11,905) store_2(k,1), store_2ea(k,1),   !temp, temp expected accuracy
     &                  store_2(k,2), store_2ea(k,2),   !dew point, dew point exp. accuracy
     &                  store_2(k,3), store_2ea(k,3)    !Rel hum, rh expected accuracy
c
	  write(11,907) store_3(k,1), store_3(k,2),     !wind dir, wind speed
     &                  store_3(k,3), store_3(k,4),     !wind gust dir, wind gust speed
     &                  store_3ea(k,1), store_3ea(k,2)  !dir expected accuracy, spd exp accuracy
c
	  write(11,909) store_4(k,1),                   !altimeter
     &                  store_4(k,2),                   !station pressure
     &                  store_4(k,3),                   !MSL pressure
     &                  nint(store_4(k,4)),             !3-h press change character
     &                  store_4(k,5),                   !3-h pressure change
     &                  store_4ea(k,1), store_4ea(k,2)  !pressure exp accuracy, alt exp accuracy


c
	  write(11,911) store_5(k,1), store_5ea(k,1),   !visibility, vis exp accuracy
     &                  store_5(k,2), store_5ea(k,2),   !solar, solar exp accuracy
     &                  store_5(k,3), store_5ea(k,3),   !soil/water temp, soil/water temp exp accuracy
     &                  store_5(k,4), store_5ea(k,4)    !soil moisture, soil moist exp accuracy
c
	  write(11,913)  store_6(k,1),                  !1-h precipitation
     &                   store_6(k,2),                  !3-h precipitation
     &                   store_6(k,3),                  !6-h precipitation
     &                   store_6(k,4),                  !24-h precipitation
     &                   store_6(k,5),                  !snow depth
     &                   store_6ea(k,1), store_6ea(k,2) !precip and snow exp accuracy
c
	  kkk_s = int(store_7(k,1))
	  write(11,915) kkk_s,                     !num cld layers (store_7(k,1))
     &                  store_7(k,2),              !24-h max temperature
     &                  store_7(k,3)               !24-h min temperature
c
c.....	Write the cloud data if we have any.
c
	  if(kkk_s .gt. 0) then
	    do ii=1,kkk_s
              call filter_string(store_cldamt(k,ii))
  	      write(11,917) store_cldamt(k,ii), store_cldht(k,ii)   !layer cloud amount and height
	    enddo !ii
	  endif
c
	enddo !k
c
	endfile(11)
	close(11)	

c
c..... End of data writing.  Let's go home...
c
	return
	include 'lso_formats.inc'
	end

        subroutine get_sfc_badflag(badflag_out,istatus)

cdoc    Returns "badflag" used in surface code

        include 'laps_sfc.inc'

        badflag_out = badflag

        istatus = 1
        return
        end

        subroutine get_ibadflag(ibadflag,istatus)

cdoc    Returns "ibadflag" used in surface code

        ibadflag = -99
        istatus = 1

        return
        end

         subroutine get_filetime_range(
     1                i4time_ob_b,i4time_ob_a                   ! I
     1               ,i4_contains_early,i4_contains_late        ! I
     1               ,intvl                                     ! I
     1               ,i4time_file_b,i4time_file_a)              ! O

cdoc     Determine the range of needed filetimes, given observation time range
cdoc     and other info about the files.

         integer i4time_ob_b        ! Earliest ob we want
         integer i4time_ob_a        ! Latest ob we want
         integer i4_contains_early  ! Earliest contained ob relative to filetime
         integer i4_contains_late   ! Latest contained ob relative to filetime
         integer intvl              ! Regular time interval of files
        
         character*9 a9_ob_b,a9_ob_a,a9_ft_b,a9_ft_a

!        Range of file times we want to read
         i4time_file_b = i4time_ob_b - i4_contains_late
         i4time_file_a = i4time_ob_a + i4_contains_early

!        Range of filenames at fixed intervals
         i4time_file_b = ( (i4time_file_b + intvl - 1) / intvl)*intvl     
         i4time_file_a = ( (i4time_file_a            ) / intvl)*intvl

         call make_fnam_lp(i4time_ob_b,a9_ob_b,istatus)
         call make_fnam_lp(i4time_ob_a,a9_ob_a,istatus)
         call make_fnam_lp(i4time_file_b,a9_ft_b,istatus)
         call make_fnam_lp(i4time_file_a,a9_ft_a,istatus)

         write(6,1)a9_ob_b,a9_ob_a,a9_ft_b,a9_ft_a                    
1        format(1x,'obs range ',2a10,' filetime range ',2a10)

         return
         end

c
c
      subroutine ck_array_real(var, recNum, filval, badflag)
c
      integer recNum
      real var(recNum), filval, badflag
c
      do i=1,recNum
         if(var(i) .eq. filval) var(i) = badflag
         if(abs(var(i)) .gt. 1e20) var(i) = badflag
      enddo !i

c
      return
      end
c
c
      subroutine ck_array_dble(var, recNum, filval, badflag)
c
      integer recNum
      double precision var(recNum), filval
      real badflag
c
      do i=1,recNum
         if(var(i) .eq. filval) var(i) = badflag
         if(abs(var(i)) .gt. 1e20) var(i) = badflag
      enddo !i
c
      return
      end

        subroutine get_sfc_obtime(int_obtime,i4time_lso ! I
     1                           ,i4time_ob,istatus)    ! O

        character*9 a9time                              ! L

        data icount /0/

        save icount

        if(int_obtime .lt. 0)then ! (e.g. flag value of -100)
            icount = icount + 1
            if(icount .le. 1000)then
                write(6,*)' get_sfc_obtime: int_obtime = ',int_obtime
                go to 900
            endif
        endif

        call make_fnam_lp(i4time_lso,a9time,istatus)
        if(istatus .ne. 1)go to 900

        write(a9time(6:9),1,err=900)int_obtime
 1      format(i4.4)
        
        call i4time_fname_lp(a9time,i4time_ob,istatus)
        if(istatus .ne. 1)return

        if(i4time_ob .lt. (i4time_lso - 43200))then
            i4time_ob = i4time_ob + 86400
        elseif(i4time_ob .gt. (i4time_lso + 43200))then
            i4time_ob = i4time_ob - 86400
        endif            

        if(abs(i4time_ob-i4time_lso) .gt. 43200)then
            go to 900
        endif

        istatus = 1
        return

 900    if(icount .le. 1000)then
            write(6,*)' Error in get_sfc_obtime, unresolved ob time'
        endif
        i4time_ob = i4time_lso ! Assume obtime equals the file time
        istatus = 0
        return

        end

        subroutine sfci4_to_sfchhmm(i4time_ob,            ! I
     1                              hhmm,istatus)         ! O

        ! this routine essentially the reverse operation of 
        ! get_sfc_obtime. that is, from i4time to hhmm.

        integer hhmm

        call CV_I4TIM_INT_LP (i4time_ob,NYEAR,NMONTH,NDAY,NHOUR,
     1                     NMIN,NSEC,ISTATUS)

        if(istatus .eq. 1)then
            hhmm = NHOUR*100 + NMIN
        else
            hhmm = -100 ! flag value
        endif

        return
        end        
