c
        subroutine read_local_cwb(inpath,maxobs,badflag,ibadflag
     1                            ,i4time_sys,stname       
     1                            ,lats,lons,elev                      ! O
     1                            ,i4time_ob_a,t,td,rh,pcp
     1                            ,stnp,mslp,dd,ff                     ! O
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
        real stnp(maxobs), mslp(maxobs), pcp(maxobs)
        integer*4 i4time_ob_a(maxobs)

        logical l_parse

        integer istart(20), iend(20)
        integer*4 cvt_wfo_fname13_i4time
c
        character*(*)stname(*)

        character inpath*(*), stn_id*3
     1           ,a9_to_a8*8, a9time*9, a8time*8, a6time*6, filename*80
     1           ,line*132, hhmm*4, cvt_i4time_wfo_fname13*13
     1           ,a13time_eat*13, a13time_ob_eat*13, c5_blank*5 


        character*100 c1dum, c2dum, c3dum     
c
c.....  Stuff for the mesonet metadata.
c
	real lat_master(maxobs),lon_master(maxobs),elev_master(maxobs)
c
	character stn_id_master(maxobs)*3
	character stn_name_master(maxobs)*5
c
c.....  Get the mesonet metadata (station information).
c
        call read_tmeso_stntbl(inpath,maxobs,badflag,
     &                         stn_id_master,stn_name_master,
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
        c5_blank = '     '
	do i=1,maxobs
           stname(i) = c5_blank 
           i4time_ob_a(i) = ibadflag
           t(i) = badflag
           td(i) = badflag
           rh(i) = badflag
 	   pcp(i) = badflag
 	   stnp(i) = badflag
           dd(i) = badflag
           ff(i) = badflag
	enddo !i

!       if(.true.)goto980 ! for testing

c
        i4time_file_eat = i4time_sys + 8*3600             ! convert GMT to EAT

        a13time_eat = cvt_i4time_wfo_fname13(i4time_file_eat)

        filename = 'Data.CWB.MSO.'
     1                  //a13time_eat(1:4)//'-'//a13time_eat(5:6) ! yyyy_mm
     1                  //'-'//a13time_eat(7:8)                   ! dd
     1                  //'_'//a13time_eat(10:13)                 ! hhmm
     1                  //'_m.pri'

        write(6,*)' Mesonet file ',filename

        call s_len(inpath,len_inpath)
        call s_len(filename,len_fname)
 
        num = 0

        open(11,file=inpath(1:len_inpath)//filename(1:len_fname)
     1         ,status='old',err=980)

!       if(.true.)goto980 ! for testing
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

        a13time_ob_eat = line(idash-4:idash-1)//line(idash+1:idash+2)   ! yyyymm
     1                 //line(idash+4:idash+5)//'_'                     ! dd
     1                 //line(idash+7:idash+8)//line(idash+10:idash+11) ! hhmm

        i4time_rcvd_eat = cvt_wfo_fname13_i4time(a13time_ob_eat)

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
            if(num .le. 10)write(6,*)ivar,line(istart(ivar):iend(ivar))       
        enddo ! i

        ivar = 1
        read(line(istart(ivar):iend(ivar)),301,err=400)ihr_ob,imin_ob
301     format(2i2)

        ivar = 2
        read(line(istart(ivar):iend(ivar)),*,err=400)stn_id

        ivar = 3
        if(l_parse(line(istart(ivar):iend(ivar)),'/'))then
            rspd = badflag
        else
            read(line(istart(ivar):iend(ivar)),*,err=400)rspd
        endif

        ivar = 4
        if(l_parse(line(istart(ivar):iend(ivar)),'/'))then
            idir = ibadflag
        else
            read(line(istart(ivar):iend(ivar)),*,err=400)idir
        endif

        ivar = 9
        if(l_parse(line(istart(ivar):iend(ivar)),'/'))then
            rstnp = badflag
        else
            read(line(istart(ivar):iend(ivar)),*,err=400)rstnp
        endif

        ivar = 10
        if(l_parse(line(istart(ivar):iend(ivar)),'/'))then
            rt = badflag
        else
            read(line(istart(ivar):iend(ivar)),*,err=400)rt
        endif

        ivar = 11
        if(l_parse(line(istart(ivar):iend(ivar)),'/'))then
            rtd = badflag
        else
            read(line(istart(ivar):iend(ivar)),*,err=400)rtd
        endif

        ivar = 13
        if(l_parse(line(istart(ivar):iend(ivar)),'/'))then
            slp = badflag
        else
            read(line(istart(ivar):iend(ivar)),*,err=400)slp
        endif

        go to 410

 400    write(6,*)' read error in station/variable ',num+1,ivar
        write(6,*)ivar,line(istart(ivar):iend(ivar))

        go to 990

 410    continue
c
c
c.....  Have good date/time...store ob.  Adjust/scale variables while storing.
c
        num = num + 1   ! add to count

        if(num .gt. maxobs)then
            write(6,*)' read_local_cwb error for too many obs: '
     1               ,num,maxobs
            istatus = 0
            return
        endif
c
c       Match data with metadata for this station, then store the metadata
c       in arrays.
c
        imatch = 0
	do j=1,num_master
	   if(stn_id .eq. stn_id_master(j)) then
	       lats(num) = lat_master(j)
	       lons(num) = lon_master(j)
	       elev(num) = elev_master(j)
               imatch=1
	   endif
	enddo !j

        if(imatch .eq. 0)then
            write(6,*)' No station match ',stn_id
        endif
c
 	stname(num) = stn_id//'  '
c
c       Correct received time to arrive at actual observation time
        isecofday_ob_eat = ihr_ob * 3600 + imin_ob * 60
        isecofday_rcvd_eat = i4time_rcvd_eat 
     1                    - (i4time_rcvd_eat/86400)*86400

        idiff_sec = isecofday_ob_eat - isecofday_rcvd_eat
        if(idiff_sec .lt. -43200)idiff_sec = idiff_sec + 86400
        if(idiff_sec .gt. +43200)idiff_sec = idiff_sec - 86400

        i4time_ob_eat = i4time_rcvd_eat + idiff_sec 
        i4time_ob_a(num) = i4time_ob_eat - 8*3600 ! EAT to GMT time zone change
c
        if(num .le. 100)then
            write(6,*)'date/time at station: ', stn_id,' '
     1                ,a13time_ob_eat,idiff_sec
        endif
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
           t(num) = c_to_f(rt) 
        endif
c
        if(rtd .le. -90) then
           td(num) = badflag
        else
           td(num) = c_to_f(rtd) 
        endif
c
!       if(irh .lt. 0) then
!          rh(num) = badflag
!       else
!          rh(num) = float(irh)
!       endif
c
!       if(iprecip .lt. 0) then
!          pcp(num) = badflag
!       else
!          pcp(num) = float(iprecip) * 0.1 * 0.03937 !conv mm to inch
!       endif
c
        if(rstnp .le. 0) then
           stnp(num) = badflag
        else
           stnp(num) = rstnp
        endif
c
        if(slp .le. 800. .or. slp .gt. 1100.)then
           mslp(num) = badflag
        else
           mslp(num) = slp
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
        print *,' Found ', num, ' mesonet stations.'
	istatus = 1
        return
c     
 980    continue
c
        write(6,*)' Warning: could not open mesonet data file ',filename
        num = 0
	istatus = -1
        return
c     
 990    continue
c
        print *,' ** ERROR reading mesonet data.'
        num = 0
	istatus = -1
        return
c     
        end
c
c
c
        subroutine read_tmeso_stntbl(inpath,maxobs,badflag,stn_id,
     &                               stn_name,lat,lon,elev,num,istatus)
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
	character*3 stn_id(maxobs), stn_id_in
        character*5 stn_name(maxobs), stn_name_in
c
        character inpath*(*)
c
c
c.....  Start here.  Fill the output with something, then open the 
c.....	file to read.
c
	do i=1,maxobs
	   lat(i) = badflag
	   lon(i) = badflag
	   elev(i) = badflag
	   stn_id(i) = '   '
	   stn_name(i) = '     '
	enddo !i
c
        call s_len(inpath,len_inpath)
        open(11,file=inpath(1:len_inpath)//'stn-table',status='old'
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
        read(11,900,end=550,err=990)stn_id_in,stn_name_in
     1                             ,lat_deg,lat_min,lat_sec,alat_sec       
     1                             ,lon_deg,lon_min,lon_sec,alon_sec
     1                             ,elev_m
 900    format(2x,a3,1x,a5,14x                      ! name
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
	stn_id(num) = stn_id_in
	stn_name(num) = stn_name_in
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
