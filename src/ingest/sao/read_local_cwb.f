cdis
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis
cdis
cdis   
cdis
      subroutine read_local_cwb ( inpath, maxobs_in, badflag, ibadflag,
     ~           i4time_sys, stname, lats, lons, elev,
     ~           t, t24max, t24min, td, rh, pcp1hr, pcp3hr, pcp6hr,
     ~           pcp24hr, dd, ff, wgdd, wgff, p, mslp, pcc, pc, sr, st,
     ~           num, istatus )           
        
      integer, parameter :: maxobs = 460
      integer, parameter :: maxMso = 40
      integer, parameter :: maxAgr = 20
      integer, parameter :: maxCum = 400
      integer, parameter :: maxShp = 110

      character  stname(maxobs)*5, inpath*100
      integer    pcc(maxobs)
      real  lats(maxobs), lons(maxobs), elev(maxobs)
      real  t(maxobs), t24max(maxobs), t24min(maxobs), td(maxobs)
      real  rh(maxobs), pcp1hr(maxobs), pcp3hr(maxobs), pcp6hr(maxobs)
      real  pcp24hr(maxobs), dd(maxobs), ff(maxobs), wgdd(maxobs)
      real  wgff(maxobs), p(maxobs), mslp(maxobs), pc(maxobs)
      real  sr(maxobs), st(maxobs)

!     inpath= '/pj/laps/dat/loc/'
      badflag= -999.
      ibadflag= -999

      stname = '     '
      pcc    = ibadflag
      lats   = badflag
      lons   = badflag
      elev   = badflag
      t      = badflag
      t24max = badflag
      t24min = badflag
      td     = badflag
      rh     = badflag
      pcp1hr = badflag
      pcp3hr = badflag
      pcp6hr = badflag
      pcp24hr= badflag
      ff     = badflag
      dd     = badflag
      wgff   = badflag
      wgdd   = badflag
      p      = badflag
      mslp   = badflag
      pc     = badflag
      sr     = badflag
      st     = badflag

      call read_meso_cwb (inpath, maxMso, badflag, ibadflag, i4time_sys,
     ~     stname(1:maxMso), lats(1:maxMso), lons(1:maxMso),
     ~     elev(1:maxMso), t(1:maxMso), t24max(1:maxMso), 
     ~     t24min(1:maxMso), td(1:maxMso), rh(1:maxMso), 
     ~     pcp1hr(1:maxMso), pcp3hr(1:maxMso), pcp6hr(1:maxMso),
     ~     pcp24hr(1:maxMso), dd(1:maxMso), ff(1:maxMso), 
     ~     wgdd(1:maxMso), wgff(1:maxMso), p(1:maxMso), mslp(1:maxMso),
     ~     pcc(1:maxMso), pc(1:maxMso), sr(1:maxMso), st(1:maxMso),
     ~     numMso, istatusMso)        

      n= numMso
      np= numMso +1

      call read_agr_cwb (inpath, maxAgr, badflag, i4time_sys,
     ~     stname(np:np+maxAgr), lats(np:np+maxAgr), lons(np:np+maxAgr),
     ~     elev(np:np+maxAgr), t(np:np+maxAgr), t24max(np:np+maxAgr),
     ~     t24min(np:np+maxAgr), td(np:np+maxAgr), rh(np:np+maxAgr),
     ~     pcp1hr(np:np+maxAgr), pcp3hr(np:np+maxAgr),
     ~     pcp6hr(np:np+maxAgr), pcp24hr(np:np+maxAgr), 
     ~     ff(np:np+maxAgr), wgff(np:np+maxAgr), wgdd(np:np+maxAgr),
     ~     p(np:np+maxAgr), sr(np:np+maxAgr), st(np:np+maxAgr),
     ~     numAgr, istatusAgr)     

      np= np +numAgr

      call read_cum_cwb (inpath, maxCum, badflag, i4time_sys,
     ~     stname(np:np+maxCum), lats(np:np+maxCum), lons(np:np+maxCum),
     ~     elev(np:np+maxCum), pcp1hr(np:np+maxCum), 
     ~     pcp3hr(np:np+maxCum), pcp6hr(np:np+maxCum),
     ~     numCum,istatusCum)  

      npC= np
      ntC= numMso +numAgr +numCum
      num= ntC
      np= np +numCum

      call read_shp_cwb (inpath, maxShp,badflag, i4time_sys,
     ~     stname(np:np+maxShp),lats(np:np+maxShp), lons(np:np+maxShp),
     ~     elev(np:np+maxShp), p(np:np+maxShp), t(np:np+maxShp), 
     ~     dd(np:np+maxShp), ff(np:np+maxShp),    
     ~     numShp, istatusShp)                           

      nt= numMso +numAgr +numCum +numShp

      do j= np,nt
      do i= npC,ntC
         if ( stname(i) == stname(j) ) then
            p(i)= p(j)
            t(i)= t(j)
            dd(i)= dd(j)
            ff(i)= dd(j)
         endif
      enddo 
      enddo 

      if ( istatusMso == 1 .and. istatusAgr == 1 .and. istatusCum == 1
     ~                     .and. istatusShp == 1 ) then
         istatus= 1
      else
         istatus= 0
         write(6,*) ' Mso ', istatusMso, ' Agr ', istatusAgr,
     ~              ' Cum ', istatusCum, ' Shp ', istatusShp
      endif

      end
 


      subroutine read_meso_cwb (inpath, maxobs, badflag, ibadflag,
     ~                          i4time_sys, stname, lats, lons, elev,
     ~                          t, t24max, t24min, td, rh, pcp1hr,
     ~                          pcp3hr, pcp6hr, pcp24hr, dd, ff,
     ~                          wgdd, wgff, stnp, mslp, pcc, pc, sr, st,
     ~                          num, istatus)                    
 
c======================================================================
c
c     Routine to read the CWB ASCII Mesonet files.
c     
c======================================================================
 
      real :: lats(maxobs), lons(maxobs), elev(maxobs)
      real :: t(maxobs), t24max(maxobs), t24min(maxobs), td(maxobs)
      real :: rh(maxobs), pcp1hr(maxobs), pcp3hr(maxobs), pcp6hr(maxobs)
      real :: pcp24hr(maxobs), dd(maxobs), ff(maxobs), wgdd(maxobs)
      real :: wgff(maxobs), stnp(maxobs), mslp(maxobs), pcp(maxobs)
      real :: pc(maxobs), sr(maxobs), st(maxobs)
      integer*4 :: pcc(maxobs), i4time_ob_a(maxobs)

      logical :: l_parse

c    larger arrays for istart and iend to read data to make processes smooth
      integer, parameter :: num40 = 40,  num70 = 70
      integer :: istart(num70), iend(num70), hh, flag
      integer*4  cvt_wfo_fname13_i4time
 
      character(len=*) :: stname(maxobs), inpath*(*)
      character(13)  cvt_i4time_wfo_fname13, a13time_eat, a13time_ob_eat
      character  stn_id*3, filename*80, line*320, c5_blank*5
 
c                      Stuff for the mesonet metadata.
      real  lat_master(maxobs), lon_master(maxobs), elev_master(maxobs)
 
      character :: stn_id_master(maxobs)*3, stn_name_master(maxobs)*5
 
c               Get the mesonet metadata (station information).
      call read_tmeso_stntbl (inpath, maxobs, badflag,  
     ~                        stn_id_master, stn_name_master,
     ~                        lat_master, lon_master, elev_master,
     ~                        num_master, istatus)
      if ( istatus /= 1 ) then
         write(6,*) ' Error reading mesonet station table'
         return
      endif

c    Fill the output arrays with something, then open the file to read.
 
      istatus= 0
      c5_blank= '     '
      stname= c5_blank 
      i4time_ob_a = ibadflag
      pcc = ibadflag
      t   = badflag
      td  = badflag
      rh  = badflag
      pcp = badflag
      stnp= badflag
      dd  = badflag
      ff  = badflag
      wgdd= badflag
      wgff= badflag
      pc  = badflag
      sr  = badflag
      st  = badflag

      i4time_file_eat= i4time_sys +8*3600             ! convert GMT to EAT
      a13time_eat= cvt_i4time_wfo_fname13(i4time_file_eat)

      filename= 'Data.CWB.MSO.'
     ~           //a13time_eat(1:4)//'-'//a13time_eat(5:6)            ! yyyy_mm
     ~           //'-'//a13time_eat(7:8) //'_' //'0000' //'_h.pri'    ! dd

      write(6,*) ' Mesonet file ', filename

      call s_len ( inpath, len_inpath )
      call s_len ( filename, len_fname )
 
      num= 0
      num_keep= 0

      open (11,file=inpath(1:len_inpath)//filename(1:len_fname), 
     ~         status='old',err=980)

c.....  This starts the read loop.  Since we don't know how many 
c.....  stations we have, read until we hit the end of file.
      
 500  flag= 0
 
      read (11,'(a)',end=550,err=990) line

c    Find first dash in time portion (two spaces before last dash in string)
      do i= 1,300
         if (line(i:i) == ':')  exit
      enddo 
 
 200  idash= i -9
      a13time_ob_eat= line(idash-4:idash-1) //line(idash+1:idash+2)     ! yyyymm
     ~                //line(idash+4:idash+5) //'_'                     ! dd
     ~                //line(idash+7:idash+8) //line(idash+10:idash+11) ! hhmm

      i4time_rcvd_eat= cvt_wfo_fname13_i4time(a13time_ob_eat)

c                  Parse the string into contiguous characters
      istart= 0
      iend=   0

      ivar= 1
      istart(1)= 1

      do i= 1,idash
         if ( line(i:i) == ' ' .and. line(i-1:i-1) /= ' ' ) then
            iend(ivar)= i-1
         endif

         if ( line(i:i) == ' ' .and. line(i+1:i+1) /= ' ' ) then
            ivar= ivar +1
            istart(ivar)= i+1
         endif
      enddo

      if ( istart(num40) /= 0 )  go to 500

      ivar= 1
      read (line(istart(ivar):iend(ivar)),301,err=400) ihr_ob, imin_ob
 301  format (2i2)

      read (a13time_eat(10:11),'(i2)') hh
      if ( ihr_ob /= hh )  go to 500

      ivar= 2
      read (line(istart(ivar):iend(ivar)),*,err=400) stn_id

      ivar= 3
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rstnp= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=400) rstnp
      endif

      ivar= 4
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         slp= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=400) slp 
      endif

      ivar= 5
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         ipcc= ibadflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=400) ipcc
      endif

      ivar= 6
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rpc= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=400) rpc
      endif

      ivar= 7
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rt= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=400) rt
      endif

      ivar= 10
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rtd= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=400) rtd
      endif

      ivar= 11
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         idir= ibadflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=400) idir
      endif

      ivar= 12
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rspd= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=400) rspd
      endif

      ivar= 13
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rwgff= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=400) rwgff
      endif

      ivar= 14
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         iwgdd= ibadflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=400) iwgdd
      endif

      ivar= 15
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rpcp= badflag
      elseif ( l_parse(line(istart(ivar):iend(ivar)),'000T') ) then
         rpcp= 0
      else
         read (line(istart(ivar):iend(ivar)),*,err=400) rpcp
      endif

      ivar= 21
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rsr= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=400) rsr
      endif

      ivar= 24
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         irh= ibadflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=400) irh
      endif

      ivar= 28
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rst= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=400) rst
      endif

      if ( num == 0 )  go to 410

      do i= 1,num
         if ( stn_id//'  ' == stname(i) ) then 
            num= i
            flag= 1
            go to 530
         else
            cycle
         endif
      enddo
      go to 410

 400  write(6,*) ' read error in station/variable ', num+1, ivar
      write(6,*) ivar, line(istart(ivar):iend(ivar))
      go to 990

c    Have good date/time...store ob.  Adjust/scale variables while storing.
 410  num= num_keep +1                    ! add to count

      if ( num > maxobs ) then
         write(6,*) ' read_local_cwb error for too many obs: ',
     ~              num, maxobs
         istatus= 0
         return
      endif
 
c Match data with metadata for this station, then store the metadata in arrays.
      imatch= 0
      do j= 1,num_master
         if ( stn_id == stn_id_master(j) ) then
            lats(num)= lat_master(j)
            lons(num)= lon_master(j)
            elev(num)= elev_master(j)
            imatch=1
         endif
      enddo 

      if ( imatch == 0 ) then
         write(6,*) ' No station match ',stn_id
      endif
 
      stname(num)= stn_id//'  '
 
c            Correct received time to arrive at actual observation time
      isecofday_ob_eat= ihr_ob *3600 +imin_ob *60
      isecofday_rcvd_eat= i4time_rcvd_eat -(i4time_rcvd_eat/86400)*86400       

      idiff_sec= isecofday_ob_eat - isecofday_rcvd_eat
      if ( idiff_sec < -43200 )  idiff_sec= idiff_sec +86400
      if ( idiff_sec > +43200 )  idiff_sec= idiff_sec - 86400

      i4time_ob_eat= i4time_rcvd_eat +idiff_sec 
      i4time_ob_a(num)= i4time_ob_eat -8*3600   ! EAT to GMT time zone change
 
      if ( num <= 100 ) 
     ~   write(6,*) 'date/time at station: ', stn_id, ' ', 
     ~              a13time_ob_eat, idiff_sec
 
530   if ( rstnp <= 0 ) then
         stnp(num)= badflag
      else
         stnp(num)= rstnp
      endif
 
      if ( slp <= 800. .or. slp > 1100. ) then
         mslp(num)= badflag
      else
         mslp(num)= slp
      endif
 
      pcc(num)= ipcc

      if ( rpc <= 0 ) then
         pc(num)= badflag
      else
         pc(num)= rpc
      endif
 
      if ( rt <= -90 ) then
         t(num)= badflag
      else
         if ( rt > 50. )  rt= - (rt - 50.)
         t(num)= c_to_f(rt) 
      endif
 
      if ( rtd <= -90 ) then
         td(num)= badflag
      else
         if ( rtd > 50. ) rtd= - (rtd - 50.)
         td(num)= c_to_f(rtd) 
      endif
 
      if ( idir > 36 .or. idir < 0 ) then
         dd(num)= badflag
      else
         dd(num)= float(idir * 10)
      endif
 
      if ( rspd < 0 ) then
         ff(num)= badflag
      else
         ff(num)= (rspd * 0.1) * 1.94254           ! conv m/s to kt
      endif
 
      if ( rwgff < 0 ) then
         wgff(num)= badflag
      else
         wgff(num)= (rwgff * 0.1) * 1.94254        ! conv m/s to kt
      endif
 
      if ( iwgdd > 36 .or. iwgdd < 0 ) then
         wgdd(num)= badflag
      else
         wgdd(num)= float(iwgdd * 10)
      endif
 
      if ( rpcp < 0 ) then
         pcp1hr(num)= badflag
      else
         pcp1hr(num)= rpcp * 0.1 * 0.03937         ! conv mm to inch
      endif
 
      if ( rsr < 0 ) then
         sr(num)= badflag
      else
         sr(num)= rsr /1000. /3600.                ! conv mJ/m/m to watt/m/m
      endif
 
      if ( irh < 0 ) then
         rh(num)= badflag
      else
         rh(num)= float(irh)
      endif
 
      if ( rst < 0 ) then
         st(num)= badflag
      else
         st(num)= rst
      endif
 
c                          Go back for the next ob.
      if ( flag == 0 )  num_keep= num 
      num= num_keep
      go to 500
 
 550  call mso_t24_pcp (inpath, filename, stname, t24max, t24min,
     ~                  pcp3hr, pcp6hr, pcp24hr, hh, maxobs, num,
     ~                  badflag)

      do i= 1,maxobs
         t24max(i)= c_to_f(t24max(i))
         t24min(i)= c_to_f(t24min(i))
      enddo

c                        Hit end of file...that's it.
      write(6,*) ' Found ', num, ' mesonet stations.'
      istatus= 1
      return
      
 980  write(6,*) ' Warning: could not open mesonet data file ',filename
      num= 0
      istatus= -1
      return

 990  write(6,*) ' ** ERROR reading mesonet data.'
      num= 0
      istatus= -1
      return
      
      end
 
 
 
      subroutine read_tmeso_stntbl (inpath, maxobs, badflag, stn_id,
     ~                           stn_name, lat, lon, elev, num, istatus)       
 
c======================================================================
c
c     Routine to read station information for the CWB ASCII Mesonet 
c	data.
c     
c======================================================================
 
      real             :: lat(maxobs), lon(maxobs), elev(maxobs)
      character(3)     :: stn_id(maxobs), stn_id_in
      character(5)     :: stn_name(maxobs), stn_name_in
      character(len=*) :: inpath
 
      lat= badflag
      lon= badflag
      elev= badflag
      stn_id= '   '
      stn_name= '     '
 
      call s_len ( inpath, len_inpath )
      open (13,file=inpath(1:len_inpath)//'stn-table',status='old',
     ~                                                err=990)
 
      num= 0

c                 Skip header comments at the top of the file
      do iread= 1,2
         read (13,*,end=550,err=990)
      enddo
 
c.....  This starts the station read loop.  Since we don't know how many 
c.....  stations we have, read until we hit the end of file.

 500  read (13,900,end=550,err=990) stn_id_in,stn_name_in,
     ~                              lat_deg,lat_min,lat_sec,alat_sec,       
     ~                              lon_deg,lon_min,lon_sec,alon_sec,
     ~                              elev_m
 900  format (2x,a3,1x,a5,14x,                     ! name
     ~        i2,2x,i2,1x,i2,1x,f3.0,4x,           ! lat
     ~        i3,2x,i2,1x,i2,1x,f3.0,              ! lon
     ~        f12.0)                               ! elevation
 
c         Move station info to arrays for sending to calling routine.
 
      alat= float(lat_deg) +float(lat_min)/60. 
     ~                     +(float(lat_sec) +alat_sec) /3600.
      alon= float(lon_deg) +float(lon_min)/60. 
     ~                     +(float(lon_sec) +alon_sec) /3600.

      num= num +1
      stn_id(num)= stn_id_in
      stn_name(num)= stn_name_in
      lat(num)= alat
      lon(num)= alon
      elev(num)= elev_m
 
c                         Go back for the next ob.
      go to 500
 
c                        Hit end of file...that's it.
 550  write(6,*) ' Found ', num,
     ~           ' mesonet stations in the station table.' 
      istatus= 1
      return
      
 990  write(6,*) stn_id_in, stn_name_in,
     ~           lat_deg, lat_min, lat_sec, alat_sec,       
     ~           lon_deg, lon_min, lon_sec, alon_sec, elev_m
      write(6,*) ' ** ERROR reading mesonet station table'
      istatus= 0
      return

      end


 
      subroutine mso_t24_pcp (inpath, filename, stname, t24max, t24min,
     ~                        pcp3hr, pcp6hr, pcp24hr, ih, maxobs, num,
     ~                        badflag)

      integer, parameter :: num24 = 24,  num40 = 40,  num70 = 70
      character(len=*) :: stname(maxobs), inpath
      character(2)     :: yy, mm, dd, hh
      character        :: filename*35, fileDummy*35, line*320, stn_id*3
      logical :: l_parse
      integer :: istart(num70), iend(num70), count(num24), d(12) 
      integer :: hr, flag
      real    :: t24max(maxobs), t24min(maxobs)
      real    :: tmax(maxobs,num24), tmin(maxobs,num24)
      real    :: pcp3hr(maxobs), pcp6hr(maxobs), pcp24hr(maxobs)
      real    :: p1hr(maxobs,-4:num24)
      data    d / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /

      count= 0
      stn_id=  '     '
      tmax= badflag
      tmin= badflag
      p1hr=    0
      pcp3hr=  0
      pcp6hr=  0
      pcp24hr= 0
      maxrcd= 3000

      read (filename(14:17),'(i4)') iy
      read (filename(19:20),'(i2)') im
      read (filename(22:23),'(i2)') id

c                 open two files to read data of 24 hours
      do l= 0,1
         id= id -l

         if ( id < 1 ) then
            im= im -1
  
            if ( im < 1 ) then
               im= 12
               iy= iy -1
            endif

            id= d(im)
         endif

         iy= iy -2000
         call i2a ( iy, yy )
         call i2a ( im, mm )
         call i2a ( id, dd )
         fileDummy= 'Data.CWB.MSO.' //'20' //yy //'-' //mm //'-' //dd
     ~                                          //'_' //'0000_h.pri'     
         iy= iy +2000

         call s_len ( inpath, len_inpath )
         call s_len ( fileDummy, len_fname )

         select case ( l ) 
         case ( 0 )
            rewind (11)
         case ( 1 )
            open (11,file=inpath(1:len_inpath)//fileDummy(1:len_fname),
     ~               status='old',err=980)
         end select

         do 200 k= 1,maxrcd
            read (11,'(a)',iostat=istat) line

            if ( istat == -1 .or. istat == -2 )  exit
            if ( istat > 0 )  go to 990
               
c    Find first dash in time portion (two spaces before last dash in string)
            do i= 1,300
               if (line(i:i) == ':')  exit
            enddo

            idash= i -9

c                 Parse the string into contiguous characters
            istart= 0
            iend=   0

            ivar= 1
            istart(1)= 1

            do i= 1,idash
               if ( line(i:i) == ' ' .and. line(i-1:i-1) /= ' ' ) then      
                  iend(ivar)= i-1
               endif

               if ( line(i:i) == ' ' .and. line(i+1:i+1) /= ' ' ) then
                  ivar= ivar +1
                  istart(ivar)= i+1
               endif
            enddo

c                              avoid daily data
            if ( istart(num40) /= 0 )  cycle

            ivar= 1
            read (line(istart(ivar):iend(ivar)),'(2i2)',err=199) hr, mn
            if ( l == 0 .and. hr >  ih )  cycle  
            if ( l == 1 .and. hr <= ih )  cycle  
            if ( hr == 0 )  hr= 24
           
            ivar= 2
            read (line(istart(ivar):iend(ivar)),*,err=200) stn_id

            do i= 1,num
            if ( stname(i) == stn_id//'  ' ) then
               ivar= 8
               if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
                 tmax(i,hr)= badflag
               else
                 read(line(istart(ivar):iend(ivar)),*,err=199)tmax(i,hr)
               endif

               ivar= 9
               if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
                 tmin(i,hr)= badflag
               else
                 read(line(istart(ivar):iend(ivar)),*,err=199)tmin(i,hr)    
               endif

               ivar= 15
               if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
                 p1hr(i,hr)= badflag
               elseif(l_parse(line(istart(ivar):iend(ivar)),'000T'))then
                 p1hr(i,hr)= 0
               else
                 read(line(istart(ivar):iend(ivar)),*,err=199)p1hr(i,hr)
               endif
            endif
            enddo

            cycle

 199        write(6,*) ' read error in station/variable ', j, ivar
            write(6,*) ivar,line(istart(ivar):iend(ivar))
 200     enddo

         if ( ih == 23 )  exit 
 500  enddo 

      do i= 1, num
         t24max(i)= tmax(i,1)
         t24min(i)= tmin(i,1)

         do j= 2,num24
            flag= 0

            if ( tmax(i,j) > 50. .or. tmax(i,j) < -90. ) then
               tmax(i,j)= badflag
               exit
            elseif ( tmin(i,j) > 50. .or. tmin(i,j) < -90. ) then
               tmin(i,j)= badflag
               exit
            endif

            if ( t24max(i) < tmax(i,j) )  t24max(i)= tmax(i,j)
            if ( t24min(i) > tmin(i,j) )  t24min(i)= tmin(i,j)

            flag= 1
         enddo

c  once there is any data missing in tmax or tmin, flag= 0 from prior do loop
         if ( flag /= 1 ) then
            t24max(i)= badflag
            t24min(i)= badflag
            write(6,*) 'too few data to obtain Tmax/Tmin for ',
     ~                 stname(i), ' mesonet station ', j
            cycle
         endif
 900  enddo

c                         calculate accumulated rain gauge
      do i= 1,num
         p1hr(i, 0)= p1hr(i,24)
         p1hr(i,-1)= p1hr(i,23)
         p1hr(i,-2)= p1hr(i,22)
         p1hr(i,-3)= p1hr(i,21)
         p1hr(i,-4)= p1hr(i,20)
      enddo

      do i= 1,num
      do j= ih,ih-2,-1
         if ( p1hr(i,j) == badflag ) then
            pcp3hr(i)= badflag
            write(6,*) 'too few data to estimate pcp3hr for ',
     ~                 stname(i), ' mesonet station ', j
            exit
         else
            pcp3hr(i)= pcp3hr(i) +p1hr(i,j)
         endif
      enddo
      enddo

      do i= 1,num
      do j= ih,ih-5,-1
         if ( p1hr(i,j) == badflag ) then
            pcp6hr(i)= badflag
            write(6,*) 'too few data to estimate pcp6hr for ',
     ~                 stname(i), ' mesonet station ', j
            exit
         else
            pcp6hr(i)= pcp6hr(i) +p1hr(i,j)
         endif
      enddo
      enddo

      do i= 1,num
      do j= 1,num24
         if ( p1hr(i,j) == badflag ) then
            pcp24hr(i)= badflag
            write(6,*) 'too few data to estimate pcp24hr for ',
     ~                 stname(i), ' mesonet station ', j
            exit
         else
            pcp24hr(i)= pcp24hr(i) +p1hr(i,j)
         endif
      enddo
      enddo

      return

 980  write(6,*) ' Warning: could not open mesonet data file ',
     ~           inpath(1:len_inpath)//fileDummy(1:len_fname)
      istatus= -1
      return
        
 990  write(6,*) ' ** ERROR reading mesonet data.'
      num= 0
      istatus= -1
      return

      end
        


      subroutine read_agr_cwb (inpath, maxobs, badflag, i4time_sys,
     ~                         stname, lats, lons, elev,
     ~                         t, t24max, t24min, td, rh,               ! O
     ~                         pcp1hr, pcp3hr, pcp6hr, pcp24hr,         ! O
     ~                         ff, wgff, wgdd, p, sr, st,               ! O
     ~                         num, istatus)                            ! O
 
c======================================================================
c
c     Routine to read the CWB ASCII agricultural station files.
c     
c======================================================================
 
      real  lats(maxobs), lons(maxobs), elev(maxobs)
      real  t(maxobs), t24max(maxobs), t24min(maxobs), td(maxobs)
      real  rh(maxobs), pcp1hr(maxobs), pcp3hr(maxobs), pcp6hr(maxobs)
      real  pcp24hr(maxobs), ff(maxobs), wgdd(maxobs), wgff(maxobs)
      real  p(maxobs), sr(maxobs), st(maxobs)       

      logical  l_parse

      integer  istart(30), iend(30)
 
      character(*)  :: stname(maxobs), inpath
      character(13) :: cvt_i4time_wfo_fname13, a13time_eat
      character     :: filename*80, line*320, c5_blank*5
 
c                        Stuff for the agricultural data.
      real :: lat_master(maxobs),lon_master(maxobs),elev_master(maxobs)
      character(5) :: stn_name_master(maxobs)
 
c               Get the agricultural metadata (station information).
      call read_agr_stntbl (inpath, maxobs, badflag, stn_name_master,
     ~                      lat_master, lon_master, elev_master,
     ~                      num_master, istatus)
      if ( istatus /= 1 ) then
         write(6,*) ' Error reading agricultural station table'
         return
      endif

      istatus=  0
      c5_blank= '     '
      stname=  c5_blank
      t=       badflag
      t24max=  badflag
      t24min=  badflag
      td=      badflag
      rh=      badflag
      pcp1hr=  badflag
      pcp3hr=  badflag
      pcp6hr=  badflag
      pcp24hr= badflag
      ff=      badflag
      wgdd=    badflag
      wgff=    badflag
      p=       badflag
      sr=      badflag
      st=      badflag

      i4time_file_eat= i4time_sys +8*3600             ! convert GMT to EAT
      a13time_eat= cvt_i4time_wfo_fname13(i4time_file_eat)

      filename= 'Data.CWB.AGR.'
     ~           //a13time_eat(1:4) //'-' //a13time_eat(5:6)   ! yyyy_mm
     ~           //'-' //a13time_eat(7:8)                      ! dd
     ~           //'_' //a13time_eat(10:11) //'00.agrsv1'      ! hh

      write(6,*) ' agricultural data file ', filename

      call s_len ( inpath, len_inpath )
      call s_len ( filename, len_fname )
 
      num= 0
      idash= 160

      open (11,file=inpath(1:len_inpath)//filename(1:len_fname),
     ~         status='old',err=980)

      do 300 i= 1,6
 300     read (11,*)

c.....  This starts the read loop.  Since we don't know how many 
c.....  stations we have, read until we hit the end of file.
      
      do 420 j= 1,maxobs

         read (11,'(a)',end=500,err=420) line

c                  Parse the string into contiguous characters
         ivar= 1
         istart(1)= 1

         do i= 1,idash
            if (line(i:i) == ' ' .and. line(i-1:i-1) /= ' ') then
               iend(ivar)= i-1
            endif

            if (line(i:i) == ' ' .and. line(i+1:i+1) /= ' ') then
               ivar= ivar +1
               istart(ivar)= i+1
            endif
         enddo
c        if ( istart(2) == 2 .and. iend(2) == 2 )  go to 500
         if ( istart(2) == 2 .and. iend(2) == 2 )  exit

         ivar= 1
         read (line(istart(ivar):iend(ivar)),'(a5)',err=400) stname(j)

c              ***** means instruments are maintained
c              ignore all observations of the station at that time
         ivar= 3
         if ( l_parse(line(istart(ivar):iend(ivar)),'N/A') ) then
            t(j)= badflag
         elseif (l_parse(line(istart(ivar):iend(ivar)),'*****')) then
            cycle
         else
            read (line(istart(ivar):iend(ivar)),*,err=400) t(j)
         endif

         ivar= 4
         if ( l_parse(line(istart(ivar):iend(ivar)),'N/A') ) then
            t24max(j)= badflag
         else
            read (line(istart(ivar):iend(ivar)),*,err=400) t24max(j)
         endif

         ivar= 5
         if ( l_parse(line(istart(ivar):iend(ivar)),'N/A') ) then
            t24min(j)= badflag
         else
            read (line(istart(ivar):iend(ivar)),*,err=400) t24min(j)
         endif

         ivar= 6
         if ( l_parse(line(istart(ivar):iend(ivar)),'N/A') ) then
            td(j)= badflag
         else
            read (line(istart(ivar):iend(ivar)),*,err=400) td(j)
         endif

         ivar= 7
         if ( l_parse(line(istart(ivar):iend(ivar)),'N/A') ) then
            rh(j)= badflag
         else
            read (line(istart(ivar):iend(ivar)),*,err=400) rh(j)
         endif

         ivar= 8
         if ( l_parse(line(istart(ivar):iend(ivar)),'N/A') ) then
            pcp1hr(j)= badflag
         else
            read (line(istart(ivar):iend(ivar)),*,err=400) pcp1hr(j)
         endif

         ivar= 10
         if ( l_parse(line(istart(ivar):iend(ivar)),'N/A') ) then
            ff(j)= badflag
         else
            read ( line(istart(ivar):iend(ivar)),*,err=400) ff(j)
         endif

         ivar= 11
         if ( l_parse(line(istart(ivar):iend(ivar)),'N/A') ) then
            wgff(j)= badflag
         else
            read (line(istart(ivar):iend(ivar)),*,err=400) wgff(j)
         endif

         ivar= 12
         if ( l_parse(line(istart(ivar):iend(ivar)),'N/A') ) then
            wgdd(j)= badflag
         else
            read (line(istart(ivar):iend(ivar)),*,err=400) wgdd(j)
         endif

         ivar= 14
         if ( l_parse(line(istart(ivar):iend(ivar)),'N/A') ) then
            sr(j)= badflag
         else
            read (line(istart(ivar):iend(ivar)),*,err=400) sr(j)
         endif

         ivar= 17
         if ( l_parse(line(istart(ivar):iend(ivar)),'N/A') ) then
            p(j)= badflag
         else
            read (line(istart(ivar):iend(ivar)),*,err=400) p(j)
         endif

         ivar= 18
         if ( l_parse(line(istart(ivar):iend(ivar)),'N/A') ) then
            st(j)= badflag
         else
            read (line(istart(ivar):iend(ivar)),*,err=400) st(j)
         endif

         go to 410

 400     write(6,*) ' read error in station/variable ', num+1, ivar,
     ~              line(istart(ivar):iend(ivar))
         cycle

 410     num= num +1
 420  enddo
 
c                      obtain Tmax and Tmin in 24 hours
 500  call agr_t24_pcp (inpath, filename, stname, t24max, t24min,
     ~                  pcp1hr, pcp3hr, pcp6hr, pcp24hr,
     ~                  maxobs, num, badflag)
c500  call agrT24maxmin (inpath, filename, stname, t24max, t24min,
c    ~                   maxobs, num, badflag)

c                            quality control
      if ( num > maxobs ) then
         write(6,*) ' read_local_cwb error for too many obs: ',
     ~              num, maxobs
         istatus= 0
         return
      endif
 
      do j= 1, num
c               match data with agricultural data for this station
         imatch= 0
         do i= 1,num_master
            if ( stname(j) == stn_name_master(i) ) then
               stname(j)= stn_name_master(i)
               lats(j)= lat_master(i)
               lons(j)= lon_master(i)
               elev(j)= elev_master(i)
               imatch= 1
            endif
         enddo 

         if ( imatch == 0 ) then
            write(6,*) ' No station match ', stname(j)
         endif
 
         if ( t(j) > 50. .or. t(j) <= -90 ) then
            t(j)= badflag
         else
            t(j)= c_to_f(t(j)) 
         endif
 
         if ( t24max(j) > 50. .or. t24max(j) <= -90. ) then
            t24max(j)= badflag
         else
            t24max(j)= c_to_f(t24max(j)) 
         endif
 
         if ( t24min(j) > 50. .or. t24min(j) <= -90. ) then
            t24min(j)= badflag
         else
            t24min(j)= c_to_f(t24min(j)) 
         endif
 
         if ( td(j) > 50. .or. td(j) <= -90 ) then
            td(j)= badflag
         else
            td(j)= c_to_f(td(j)) 
         endif
 
         if ( ff(j) < 0 ) then
            ff(j)= badflag
         elseif ( ff(j) > wgff(j) ) then
            ff(j)= badflag
            wgff(j)= badflag
         else
            ff(j)= ( ff(j) * 0.1 ) * 1.94254         ! conv m/s to kt
         endif
 
         if ( wgff(j) < 0 ) then
            wgff(j)= badflag
         else
            wgff(j)= ( wgff(j) * 0.1 ) * 1.94254     ! conv m/s to kt
         endif

         if ( pcp1hr(j) < 0 ) then
            pcp1hr(j)= badflag
         else
            pcp1hr(j)= pcp1hr(j) * 0.1 * 0.03937     ! conv mm to inch
         endif

         if ( pcp3hr(j) < 0 ) then
            pcp3hr(j)= badflag
         else
            pcp3hr(j)= pcp3hr(j) * 0.1 * 0.03937     ! conv mm to inch
         endif

         if ( pcp6hr(j) < 0 ) then
            pcp6hr(j)= badflag
         else
            pcp6hr(j)= pcp6hr(j) * 0.1 * 0.03937     ! conv mm to inch
         endif

         if ( pcp24hr(j) < 0 ) then
            pcp24hr(j)= badflag
         else
            pcp24hr(j)= pcp24hr(j) * 0.1 * 0.03937   ! conv mm to inch
         endif
 
         if ( sr(j) < 0 ) then
            sr(j)= badflag
         else
            sr(j)= sr(j) /1000. /3600.               ! conv mJ/m/m to watt/m/m
         endif
 
         if ( wgdd(j) > 360. .or. wgdd(j) < 0 )  wgdd(j)= badflag
         if ( rh(j) < 0 )  rh(j)= badflag
         if ( p(j)  < 0 )  p(j)=  badflag
         if ( st(j) < 0 )  st(j)= badflag
      enddo
 
c                          Hit end of file...that's it.
 550  write(6,*) ' Found ', num, ' agricultural stations.'
      istatus= 1
      return
      
 980  write(6,*) ' Warning: could not open agricultural data file ',
     ~           filename
      num= 0
      istatus= -1
      return

 990  write(6,*) ' ** ERROR reading agricultural data.'
      num= 0
      istatus= -1
      return
      
      end
 
 
 
      subroutine read_agr_stntbl (inpath, maxobs, badflag,
     ~                           stn_name, lat, lon, elev, num, istatus)
 
c======================================================================
c
c Routine to read station information for the CWB ASCII agricultural
c   station data.
c     
c======================================================================
 
      real             :: lat(maxobs), lon(maxobs), elev(maxobs)
      character(5)     :: stn_name(maxobs), stn_name_in
      character(len=*) :: inpath
 
      lat=  badflag
      lon=  badflag
      elev= badflag
      stn_name= '     '
 
      call s_len ( inpath, len_inpath )
      open (13,file=inpath(1:len_inpath)//'agrstn-tbl',status='old'
     ~                                                ,err=990)
 
      num= 0

c.....  This starts the station read loop.  Since we don't know how many 
c.....  stations we have, read until we hit the end of file.

 500  read (13,900,end=550,err=990) stn_name_in, alat, alon, elev_m
 900  format (8x,a5,4x,f7.4,1x,f8.4,f10.1)
 
c           Move station info to arrays for sending to calling routine.
      num= num +1
      stn_name(num)= stn_name_in
      lat(num)= alat
      lon(num)= alon
      elev(num)= elev_m
 
c                         Go back for the next ob.
      go to 500
 
c                        Hit end of file...that's it.
 550  write(6,*) ' Found ', num,
     ~           ' agricultural stations in the station table.' 
      istatus= 1
      return
      
 990  write(6,*) stn_name_in, alat, alon, elev_m
      write(6,*) ' ** ERROR reading agricultural station table'
      istatus= 0
      return

      end

 

      subroutine agr_t24_pcp (inpath, filename, stname, t24max, t24min,
     ~                        pcp1hr, pcp3hr, pcp6hr, pcp24hr,
     ~                        maxobs, num, badflag)

      integer, parameter :: num23= 23 
      character(len=*) :: stname(*)
      character(5)     :: stn(num23,maxobs)
      character(2)     :: yy, mm, dd, hh
      character        :: filename*35, fileDummy*35, inpath*50, line*160
      integer :: istart(30), iend(30), d(12), count(num23), flag
      real :: t24max(maxobs), t24min(maxobs), pcp1hr(maxobs)
      real :: pcp3hr(maxobs), pcp6hr(maxobs), pcp24hr(maxobs)
      real :: tmax(num23,maxobs), tmin(num23,maxobs), p1hr(num23,maxobs)
      logical :: l_parse
      data  d / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /

      stn=    '     '
      tmax=   badflag
      tmin=   badflag
      p1hr=   badflag
      pcp3hr=  pcp1hr
      pcp6hr=  pcp1hr
      pcp24hr= pcp1hr
      count=   0
      idash= 160

      read (filename(14:17),'(i4)') iy
      read (filename(19:20),'(i2)') im
      read (filename(22:23),'(i2)') id
      read (filename(25:26),'(i2)') ih

      k_num23 : do k= 1,num23

         ih= ih -1
         if ( ih < 0 ) then
            ih= ih +24
            id= id -1

            if ( id < 1 ) then
               im= im -1
  
               if ( im < 1 ) then
                  im= 12
                  iy= iy -1
               endif

               id= d(im)
            endif
         endif

         iy= iy -2000
         call i2a ( iy, yy )
         call i2a ( im, mm )
         call i2a ( id, dd )
         call i2a ( ih, hh )
         fileDummy= 'Data.CWB.AGR.' //'20' //yy //'-' //mm
     ~                               //'-' //dd //'_' //hh //'00.agrsv1'     
         iy= iy +2000

         call s_len ( inpath, len_inpath )
         call s_len ( fileDummy, len_fname )

         idash= 160

         filename= inpath(1:len_inpath)//fileDummy(1:len_fname)
         open (24,file=inpath(1:len_inpath)//fileDummy(1:len_fname),
     ~            status='old',err=980)

         do i= 1,6
            read (24,*)
         enddo

         do 200 j= 1,maxobs

            read (24,'(a)',end=200,err=200) line

c                 Parse the string into contiguous characters
            ivar= 1
            istart(1)= 1

            do i= 1,idash
               if ( line(i:i) == ' ' .and. line(i-1:i-1) /= ' ' ) then      
                  iend(ivar)= i-1
               endif

               if ( line(i:i) == ' ' .and. line(i+1:i+1) /= ' ' ) then
                  ivar= ivar +1
                  istart(ivar)= i+1
               endif
            enddo

c               stop reading and open another file when meet ^M
            if ( istart(2) == 2 .and. iend(2) == 2 )  cycle k_num23

            ivar= 1
            read (line(istart(ivar):iend(ivar)),'(a5)',err=399) stn(k,j)
           
            ivar= 4
            if ( l_parse(line(istart(ivar):iend(ivar)),'N/A') ) then
               tmax(k,j)= badflag
            elseif (l_parse(line(istart(ivar):iend(ivar)),'*****')) then
               cycle k_num23
            else
               read (line(istart(ivar):iend(ivar)),*,err=399) tmax(k,j)
            endif

            ivar= 5
            if ( l_parse(line(istart(ivar):iend(ivar)),'N/A') ) then
               tmin(k,j)= badflag
            elseif (l_parse(line(istart(ivar):iend(ivar)),'*****')) then
               cycle k_num23
            else
               read (line(istart(ivar):iend(ivar)),*,err=399) tmin(k,j)       
            endif

            ivar= 8
            if ( l_parse(line(istart(ivar):iend(ivar)),'N/A') ) then
               p1hr(k,j)= badflag
            elseif (l_parse(line(istart(ivar):iend(ivar)),'*****')) then
               cycle k_num23
            else
               read (line(istart(ivar):iend(ivar)),*,err=399) p1hr(k,j)
            endif

            count(k)= count(k) +1
 200     continue
         cycle k_num23

 399     write(6,*) ' read error in station/variable ', j, ivar
         write(6,*) ivar,line(istart(ivar):iend(ivar))
           
      enddo k_num23
c500  continue

      t_stations : do i= 1,num
      t_hours    : do k= 1,num23
         flag= 0

         do j= 1,count(k)
            if ( stname(i) == stn(k,j) ) then
               if ( t24max(i) < tmax(k,j) )  t24max(i)= tmax(k,j)
               if ( t24min(i) > tmin(k,j) )  t24min(i)= tmin(k,j)
               flag= 1
            endif
         enddo

         if ( flag /= 1 ) then
            t24max(i)= badflag
            t24min(i)= badflag
            write(6,*) 'too few data to obtain Tmax/Tmin for ',
     ~                 stname(i), ' agricultral station ', k
            cycle t_stations
         endif
      enddo t_hours
 900  enddo t_stations

      pcp_stations : do i= 1,num
      pcp_hours    : do k= 1,num23
         flag= 0

         do j= 1,count(k)
            if ( stname(i) == stn(k,j) ) then
               if ( p1hr(k,j) /= badflag ) then
                  if ( k < 3 ) then
                     pcp3hr(i)= pcp3hr(i) +p1hr(k,j)
                  endif

                  if ( k < 6 ) then
                     pcp6hr(i)= pcp6hr(i) +p1hr(k,j)
                  endif

                  if ( k < 24 ) then
                     pcp24hr(i)= pcp24hr(i) +p1hr(k,j)
                  endif

                  flag= 1
               else
                  flag= -1
               endif
            endif
         enddo

         if ( flag /= 1 ) then
            if ( k < 3 ) then
               pcp3hr(i)= badflag
               pcp6hr(i)= badflag
               pcp24hr(i)= badflag
            endif

            if ( k < 6 ) then
               pcp6hr(i)= badflag
               pcp24hr(i)= badflag
            endif

            if ( k < 24 ) then
               pcp24hr(i)= badflag
            endif

            write(6,*) 'too few data to obtain accumulated pcp for ',
     ~                 stname(i), ' agricultural station ', k
            cycle pcp_stations
         endif
      enddo pcp_hours
      enddo pcp_stations

      istatus= 1
      return

 980  write(6,*) ' Warning: could not open agricultural data file ',
     ~           inpath(1:len_inpath)//fileDummy(1:len_fname)
      istatus= -1
      return
        
      end
        


      subroutine read_cum_cwb(inpath, maxobs, badflag,
     ~                        i4time_sys, stname,      
     ~                        lats, lons, elev,                     ! O
     ~                        pcp1hr, pcp3hr, pcp6hr,               ! O
     ~                        num, istatus)                         ! O
c
c======================================================================
c
c     Routine to read the CWB ASCII rain gauge files.
c     
c======================================================================
c
      integer, parameter :: loopnum = 400 
      real :: lats(maxobs), lons(maxobs), elev(maxobs)
      real :: pcp1hr(maxobs), pcp3hr(maxobs), pcp6hr(maxobs)

      integer :: dataMin(maxobs)
 
      character(len=*)  :: stname(maxobs), inpath
      character(len=80) :: filename
      character(len=13) :: cvt_i4time_wfo_fname13, a13time_sys

c.....  Stuff for the rain gauge metadata.
 
      real  lat_master(maxobs), lon_master(maxobs), elev_master(maxobs)
 
      character(5)  stn_name_master(maxobs), c5_blank
      character(4)  stn_id_master(maxobs), stn_id(maxobs)
      character(1)  pcp1hrQ(maxobs), pcp3hrQ(maxobs), pcp6hrQ(maxobs)
 
c.....  Get the rain gauge metadata (station information).
 
      call read_cum_stntbl(inpath, maxobs, badflag,
     ~                     stn_id_master, stn_name_master,
     ~                     lat_master, lon_master, elev_master,
     ~                     num_master, istatus)
      if ( istatus /= 1 ) then
         write(6,*) ' Error reading rain gauge station table'
         return
      endif

c   fill the output arrays with something, then open the file to read
 
      istatus= 0
      c5_blank= '     '
      stname= c5_blank 
      pcp1hr= badflag
      pcp3hr= badflag
      pcp6hr= badflag

      a13time_sys= cvt_i4time_wfo_fname13(i4time_sys)

      filename= 'Data.CWB.CUM_UDD.'
     ~             //a13time_sys(1:4)//'-'//a13time_sys(5:6)   ! yyyy_mm
     ~             //'-'//a13time_sys(7:8)                     ! dd
     ~             //'_'//a13time_sys(10:11) //'.cum_hr'       ! hh

      write(6,*) ' Mesonet file ',filename

      call s_len ( inpath, len_inpath )
      call s_len ( filename, len_fname )
 
      num= 0

      open (11,file=inpath(1:len_inpath)//filename(1:len_fname),
     ~         status='old',err=980)

c.....  This starts the read loop.  Since we don't know how many 
c.....  stations we have, read until we hit the end of file.
      
 500  continue

      do i= 1,loopNum
         read (11,511,end=520,err=990) stn_id(i), lats(i), lons(i),
     ~                 pcp1hr(i), pcp3hr(i), pcp6hr(i), dataMin(i),
     ~                 pcp1hrQ(i), pcp3hrQ(i), pcp6hrQ(i)       

         num= num +1
      enddo
 511  format (a4, 11x, 2f8.4, 1x, 3f5.1, i5, 3a1)  

c  Match data with metadata for this station, then store the metadata in arrays.
 520  imatch= 0
 
      do i= 1,num
 
         do j= 1,num_master
	    if ( stn_id(i) == stn_id_master(j) ) then
               stname(i)= stn_name_master(i)
               elev(i)= elev_master(j)
               imatch=1
            endif
         enddo

         if ( imatch == 0 ) then
            write(6,*) ' No station match ', stn_id(i), stn_id_master(j)
         endif

c                               quality control
         if ( dataMin(i) > 60 ) then
            pcp1hr(i)= badflag
            pcp3hr(i)= badflag
            pcp6hr(i)= badflag 
         else
            if ( pcp1hr(i) < 0 .or. pcp1hrQ(i) /= 'G' ) then
               pcp1hr(i)= badflag
            else
               pcp1hr(i)= pcp1hr(i) * 0.1 * 0.03937         ! conv mm to inch
            endif
 
            if ( pcp3hr(i) < 0 .or. pcp3hrQ(i) /= 'G' ) then
               pcp3hr(i)= badflag
            else
               pcp3hr(i)= pcp3hr(i) * 0.1 * 0.03937         ! conv mm to inch
            endif
 
            if ( pcp6hr(i) < 0 .or. pcp6hrQ(i) /= 'G' ) then
               pcp6hr(i)= badflag
            else
               pcp6hr(i)= pcp6hr(i) * 0.1 * 0.03937         ! conv mm to inch
            endif
         endif
      enddo 
 
c                               end of file
 550  write(6,*) ' Found ', num, ' rain gauge stations.'
      istatus= 1
      return
 
 980  write(6,*) ' Warning: could not open rain gauge data file ',
     ~           filename
      num= 0
      istatus= -1
      return
      
 990  write(6,*) ' ** ERROR reading rain gauge data.'
      num= 0
      istatus= -1
      return

      end
 
 
 
      subroutine read_cum_stntbl (inpath, maxobs, badflag, stn_id,
     ~                           stn_name, lat, lon, elev, num, istatus)       
 
c======================================================================
c
c     Routine to read station information for the CWB ASCII rain gauge
c	station data.
c
c======================================================================
 
      real  elev(maxobs)
 
      character(len=*) :: inpath
      character(4)     :: stn_id(maxobs), stn_id_in
      character(5)     :: stn_name(maxobs), stn_name_in
      integer          :: elev_m
 
c.....  Start here.  Fill the output with something, then open the 
c.....	file to read.
 
      elev= badflag
      stn_id= '    '
      stn_name= '     '
 
      call s_len ( inpath, len_inpath )
      open (13,file=inpath(1:len_inpath)//'uddstn-tbl',status='old',
     ~                                                 err=990)
 
      num= 0

c.....  This starts the station read loop.  Since we don't know how many 
c.....  stations we have, read until we hit the end of file.
      
 500  continue
 
      read (13,900,end=550,err=990) stn_id_in, stn_name_in, elev_m
 900  format(a4,2x,a5,23x,i7)

      num= num +1
      stn_id(num)= stn_id_in
      stn_name(num)= stn_name_in
      elev(num)= real(elev_m)
 
c                          Go back for the next ob.
      go to 500
 
c                        Hit end of file...that's it.
 550  continue
 
      write(6,*) ' Found ', num,
     ~           ' rain gauge stations in the station table.' 
      istatus= 1
      return
      
 990  continue
 
      write(6,*) ' ** ERROR reading rain gauge station table'
      write(6,*) stn_id_in, stn_name_in, elev_m
      istatus= 0
      return

      end
 


      subroutine read_shp_cwb (inpath, maxobs, badflag,
     ~                         i4time_sys, stname,      
     ~                         lats, lons, elev,                      ! O
     ~                         p, t, dd, ff,                          ! O
     ~                         num, istatus)                          ! O
 
c======================================================================
c
c     Routine to read the CWB ASCII shp files.
c     
c======================================================================
 
      integer, parameter :: loopnum= 110 
      real lats(maxobs), lons(maxobs), elev(maxobs)
      real p(maxobs), t(maxobs), dd(maxobs), ff(maxobs)

      character(len=*)  :: stname(maxobs), inpath
      character(len=80) :: filename
      character(len=13) :: cvt_i4time_wfo_fname13, a13time_eat

c                         stuff for the shp metadata.
      real  lat_master(maxobs), lon_master(maxobs), elev_master(maxobs)       
      character(5) :: stn_name_master(maxobs), c5_blank
      character(4) :: stn_id_master(maxobs), stn_id(maxobs)
 
c                    Get the shp metadata (station information).
      call read_shp_stntbl (inpath, maxobs, badflag,
     ~                      stn_id_master, stn_name_master,
     ~                      lat_master, lon_master, elev_master,
     ~                      num_master, istatus)
      if ( istatus /= 1 ) then
         write(6,*)' Error reading shp station table'
         return
      endif

c   fill the output arrays with something, then open the file to read
 
      istatus= 0
      c5_blank= '     '
      stname= c5_blank 
      p=  badflag
      t=  badflag
      dd= badflag
      ff= badflag

      i4time_file_eat= i4time_sys +8*3600             ! convert GMT to EAT
      a13time_eat= cvt_i4time_wfo_fname13(i4time_file_eat)

      filename= 'Data.CWB.SHP.'
     ~           //a13time_eat(1:4)//'-'//a13time_eat(5:6)   ! yyyy_mm
     ~           //'-'//a13time_eat(7:8)                     ! dd
     ~           //'_'//a13time_eat(10:11) //'00.shp'        ! hh

      write(6,*) ' Mesonet file ',filename

      call s_len ( inpath, len_inpath )
      call s_len ( filename, len_fname )
 
      num= 0

      open (11,file=inpath(1:len_inpath)//filename(1:len_fname),
     ~         status='old',err=980)

c.....  This starts the read loop.  Since we don't know how many 
c.....  stations we have, read until we hit the end of file.
      
 500  continue

      do i= 1,loopNum
         read (11,511,end=520,err=990) stn_id(i), 
     ~                                 p(i), t(i), ff(i), dd(i) 
         num= num +1
      enddo
 511  format (5x, a4, 14x, 4f6.1)

c  Match data with metadata for this station, then store the metadata in arrays.
 
 520  imatch= 0
 
      do i= 1,num
 
         do j= 1,num_master
            if ( stn_id(i) == stn_id_master(j) ) then
               stname(i)= stn_name_master(i)
               lats(i)= lat_master(j)
               lons(i)= lon_master(j)
               elev(i)= elev_master(j)
               imatch=1
            endif
         enddo

         if ( imatch == 0 ) then
            write(6,*)' No station match ', stn_id(i), stn_id_master(j)       
         endif
 
c                               quality control
         if ( p(i) < 0 )  p(i)= badflag
         if ( t(i) < -50. ) then
            t(i)= badflag
         else
            t(i)= c_to_f(t(i))
         endif

         if ( dd(i) < 0  .or.  dd(i) > 360. )  dd(i)= badflag
         if ( ff(i) < 0 ) then
            ff(i)= badflag
         else
            ff(i)= ( ff(i) * 0.1 ) * 1.94254        ! conv m/s to kt
         endif
 
      enddo
c                               end of file
 550  write(6,*) ' Found ', num, ' shp stations.'
      istatus= 1
      return
 
 980  write(6,*) ' Warning: could not open shp data file ', filename
      num= 0
      istatus= -1
      return
      
 990  write(6,*) ' ** ERROR reading shp data.'
      num= 0
      istatus= -1
      return

      end
 
 
 
      subroutine read_shp_stntbl (inpath, maxobs, badflag, stn_id,
     ~                           stn_name, lat, lon, elev, num, istatus)       
 
c======================================================================
c
c     Routine to read station information for the CWB ASCII shp
c	station data.
c
c======================================================================
 
      real lat(maxobs), lon(maxobs), elev(maxobs)
 
      character(len=*) :: inpath
      character(4)     :: stn_id(maxobs), stn_id_in
      character(5)     :: stn_name(maxobs), stn_name_in
      integer          :: elev_m
 
      lat=  badflag
      lon=  badflag
      elev= badflag
      stn_id=   '    '
      stn_name= '     '
 
      call s_len ( inpath, len_inpath )
      open (13,file=inpath(1:len_inpath)//'shpstn-tbl',
     ~         status='old',err=990)
 
      num= 0

c.....  This starts the station read loop.  Since we don't know how many 
c.....  stations we have, read until we hit the end of file.
      
 500  read (13,900,end=550,err=990) stn_id_in, stn_name_in,
     ~                              lat_deg, lat_min, lat_sec,
     ~                              lon_deg, lon_min, lon_sec, elev_m
 900  format(a4,2x,a5,2x,3i3,2x,2i3,1x,i3,3x,i4)

c.....      Move station info to arrays for sending to calling routine.
      alat= float(lat_deg) +float(lat_min)/60. +float(lat_sec)/3600.
      alon= float(lon_deg) +float(lon_min)/60. +float(lon_sec)/3600.      

      num= num +1
      stn_id(num)= stn_id_in
      stn_name(num)= stn_name_in
      lat(num)= alat
      lon(num)= alon
      elev(num)= real(elev_m)
 
c                           Go back for the next ob.
      go to 500
 
c                         Hit end of file...that's it.
 550  continue
      write(6,*) ' Found ', num, ' shp stations in the station table.'    
      istatus= 1
      return
      
 990  continue
      write(6,*) ' ** ERROR reading shp station table'
      istatus= 0
      return

      end



      subroutine  i2a (ii,aa)

      character(2) :: aa
      integer      :: ii

      if ( ii < 10 ) then
         write (aa,'(a1,i1)') '0', ii
      else
         write (aa,'(i2)') ii
      endif

      return
      end


      subroutine get_box_size (box_size,istatus)

      istatus = 1
      box_size = 1.1

      return
      end

