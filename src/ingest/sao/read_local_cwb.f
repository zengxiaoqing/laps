c
        subroutine read_tmeso_data(inpath,maxobs,badflag,ibadflag
     1                            ,i4time_sys,stname       
     1                            ,lats,lons,elev                      ! O
     1                            ,i4time_ob_a,t,td,rh,pcp,sfcp,dd,ff  ! O
     1                            ,num,istatus)                        ! O
c
c======================================================================
c
c     Routine to read the CWB ASCII Mesonet files.
c     
c     Original:  P. Stamus, NOAA/FSL  08 Sep 1999
c     Changes:   S. Albers                         (New format)
c
c======================================================================
c
	real lats(maxobs), lons(maxobs), elev(maxobs)
        real t(maxobs), td(maxobs), rh(maxobs)
        real dd(maxobs), ff(maxobs)
        real sfcp(maxobs), pcp(maxobs)
        integer*4 i4time_ob_a(maxobs)

        integer istart(20), iend(20)
c
        character inpath*(*), stn_id*5, stname(maxobs)*5 
     1           ,a9_to_a8*8, a9time*9, a8time*8, a6time*6, filename*80
     1           ,line*132, hhmm*4, cvt_i4time_wfo_fname13*13
     1           ,a13time*13, a13time_ob*13 

        character*100 c1dum, c2dum, c3dum     
c
c.....  Stuff for the mesonet metadata.
c
	real lat_master(maxobs),lon_master(maxobs),elev_master(maxobs)
c
	character stn_master(maxobs)*5
c
c.....  Get the mesonet metadata (station information).
c
        call read_tmeso_stntbl(inpath,maxobs,badflag,
     &                         stn_master,
     &                         lat_master,lon_master,elev_master,
     &                         num_master,istatus)
	if(istatus .ne. 1)then
            write(6,*)' Error reading mesonet station table'
            return
        endif

c
c.....  Start here.  Fill the output arrays with something, then open
c.....	the file to read.
c
	istatus = 0
	do i=1,maxobs
	   stname(i)(1:5) = '     '
	   i4time_ob_a(i) = ibadflag
	   t(i) = badflag
	   td(i) = badflag
	   rh(i) = badflag
 	   pcp(i) = badflag
 	   sfcp(i) = badflag
	   dd(i) = badflag
	   ff(i) = badflag
	enddo !i

c
        i4time_file = i4time_sys + 8*3600                 ! convert GMT to EAT

        a13time = cvt_i4time_wfo_fname13(i4time_file)

        filename = 'Data.CWB.MSO.'
     1                  //a13time(1:4)//'-'//a13time(5:6) ! yyyy_mm
     1                  //'-'//a13time(7:8)               ! dd
     1                  //'_'//a13time(10:13)             ! hhmm
     1                  //'_m.pri'

        write(6,*)' Mesonet file ',filename

        call s_len(inpath,len_inpath)
        call s_len(filename,len_fname)
 
        num = 0

        open(11,file=inpath(1:len_inpath)//filename(1:len_fname)
     1         ,status='old',err=980)

        if(.true.)goto980 ! for testing
c
c.....  This starts the read loop.  Since we don't know how many 
c.....  stations we have, read until we hit the end of file.
c     
 500    continue
c
        read(11,801,end=550,err=990) line
 801    format(a)

!       Find first dash
        do i = 132,1,-1
            if(line(i:i) .eq. '-')then
                idash = i
            endif
        enddo ! i

        a13time_ob = line(idash-4:idash-1)//line(idash+1:idash+2)    ! yyyymm
     1             //line(idash+4:idash+5)//'_'                      ! dd
     1             //line(idash+7:idash+8)//line(idash+10:idash+11)  ! hhmm

!       Parse the string into contiguous characters
        ivar = 1
        istart(1) = 1

        do i = 1,idash
            if(line(i:i) .eq. ' ' .and. line(i-1:i-1) .ne. ' ')then
                iend(ivar) = i-1
            endif

            if(line(i:i) .eq. ' ' .and. line(i+1:i+1) .ne. ' ')then
                ivar = ivar + 1
                istart(ivar) = i+1
            endif
        enddo

        nvar = ivar - 1

        do ivar = 1,nvar
            write(6,*)ivar,line(istart(ivar):iend(ivar))
        enddo ! i

!       Now we can read the variables
!       read(line,901,err=990)hhmm,stn_id
!901    format(a4,a4)

!       read(line(9:132),*,err=990)rspd,idir,rdum,rdum,rdum
!    1                            ,rdum,rsfcp,rt,rtd
!    1                            ,rdum,slp

        ivar = 1
        read(line(istart(ivar):iend(ivar)),*)rspd

c
c.....  Check for valid date/time...if bad, toss this ob.
c
        write(6,*)'date/time at station: ', stn_id, a13time_ob
c
c.....  Have good date/time...store ob.  Adjust/scale variables while storing.
c
        num = num + 1   !add to count
c
	stname(num)(1:5) = stn_id(1:5)
c
        i4_mm = imin * 60
        i4time_ob_a(num) = i4time_file + i4_mm
c
        if(idir.gt.36 .or. idir.lt.0) then
           dd(num) = badflag
        else
           dd(num) = float(idir * 10)
        endif
c
        if(rspd .lt. 0) then
           ff(num) = badflag
        else
           ff(num) = (rspd * 0.1) * 1.94254 !conv m/s to kt
        endif
c
        if(rt .le. -90) then
           t(num) = badflag
        else
           t(num) = (rt * 0.1) * 9/5 + 32 !conv C to F
        endif
c
        if(rtd .le. -90) then
           td(num) = badflag
        else
           td(num) = (rtd * 0.1) * 9/5 + 32 !conv C to F
        endif
c
        if(irh .lt. 0) then
           rh(num) = badflag
        else
           rh(num) = float(irh)
        endif
c
        if(iprecip .lt. 0) then
           pcp(num) = badflag
        else
           pcp(num) = float(iprecip) * 0.1 * 0.03937 !conv mm to inch
        endif
c
        if(rsfcp .le. 0) then
           sfcp(num) = badflag
        else
           ps = rsfcp * 0.1
           if(ps .lt. 500.) then
              ps = ps + 1000.
           endif
           sfcp(num) = ps
        endif
c
c.....  Go back for the next ob.
c
        go to 500
c
c.....  Hit end of file...that's it.
c
 550    continue
c
c       Match data with metadata for each station, then store the metadata
c       in arrays.
c
        do i=1,num
	   do j=1,num_master
	      if(stname(i)(1:5) .eq. stn_master(j)(1:5)) then
		 lats(i) = lat_master(j)
		 lons(i) = lon_master(j)
		 elev(i) = elev_master(j)
	      endif
	   enddo !j
	enddo !i
c
        print *,' Found ', num, ' mesonet stations.'
	istatus = 1
        return
c     
 980    continue
c
        write(6,*)' Warning: could not open mesonet data file ',filename
	istatus = -1
        return
c     
 990    continue
c
        print *,' ** ERROR reading mesonet data.'
	istatus = -1
        return
c     
        end
c
c
c
        subroutine read_tmeso_stntbl(inpath,maxobs,badflag,stn,
     &                               lat,lon,elev,num,istatus)
c
c======================================================================
c
c     Routine to read station information for the CWB ASCII Mesonet 
c	data.
c     
c     Original:  P. Stamus, NOAA/FSL  08 Sep 1999
c     Changes:   S. Albers                            (New Format)
c
c======================================================================
c
        real lat(maxobs), lon(maxobs), elev(maxobs)
c
	integer stn_id
c
        character inpath*(*), stn_name*5, stn(maxobs)*5
c
c
c.....  Start here.  Fill the output with something, then open the 
c.....	file to read.
c
	do i=1,maxobs
	   lat(i) = badflag
	   lon(i) = badflag
	   elev(i) = badflag
	   stn(i)(1:5) = '     '
	enddo !i
c
        call s_len(inpath,len_inpath)
        open(11,file=inpath(1:len_inpath)//'stn.table',status='old'
     1                                                ,err=990)
c
        num = 0

!       Skip header comments at the top of the file
        do iread = 1,2
            read(11,*,end=550,err=990)
        enddo
c
c.....  This starts the station read loop.  Since we don't know how many 
c.....  stations we have, read until we hit the end of file.
c     
 500    continue
c
        read(11,900,end=550,err=990)stn_id,stn_name
     1                             ,lat_deg,lat_min,lat_sec,alat_sec       
     1                             ,lon_deg,lon_min,lon_sec,alon_sec
     1                             ,elev_m
 900    format(3x,i2,1x,a5,14x                      ! name
     1        ,i2,2x,i2,1x,i2,1x,f3.0,4x            ! lat
     1        ,i3,2x,i2,1x,i2,1x,f3.0               ! lon
     1        ,f12.0)                               ! elevation
c
c.....  Move station info to arrays for sending to calling routine.
c
        alat = float(lat_deg) + float(lat_min)/60. 
     1                        + (float(lat_sec) + alat_sec) /3600.
        alon = float(lon_deg) + float(lon_min)/60. 
     1                        + (float(lon_sec) + alon_sec) /3600.

	num = num + 1
	stn(num) = stn_name(3:5)//'  '
	lat(num) = alat
	lon(num) = alon
	elev(num) = elev_m
c
c.....  Go back for the next ob.
c
        go to 500
c
c.....  Hit end of file...that's it.
c
 550    continue
c
        print *,' Found ', num
     1         , ' mesonet stations in the station table.' 
        istatus = 1
        return
c     
 990    continue
c
        print *,' ** ERROR reading mesonet station table'
        istatus = 0
        return
c     
	end

        subroutine get_box_size(box_size,istatus)

        istatus = 1
        box_size = 1.1

        return
        end
