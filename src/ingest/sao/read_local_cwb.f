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
     ~           i4time_sys, rptTp, stnTp, stname, lats, lons, elev,
     ~           t, t24max, t24min, td, rh, pcp1hr, pcp3hr, pcp6hr,
     ~           pcp24hr, dd, ff, wgdd, wgff, p, mslp, pcc, pc, sr, st,
     ~           num, istatus )
        
      integer, parameter :: maxobs = 1100
      integer, parameter :: maxAgr =  20
      integer, parameter :: maxCum = 800
      integer, parameter :: maxShp = 250

      character(*)   stname(maxobs)
      character(*)  inpath
      character(6)   rptTp(maxobs), stnTp(maxobs)
      integer        pcc(maxobs)
      real  lats(maxobs), lons(maxobs), elev(maxobs)
      real  t(maxobs), t24max(maxobs), t24min(maxobs), td(maxobs)
      real  rh(maxobs), pcp1hr(maxobs), pcp3hr(maxobs), pcp6hr(maxobs)
      real  pcp24hr(maxobs), dd(maxobs), ff(maxobs), wgdd(maxobs)
      real  wgff(maxobs), p(maxobs), mslp(maxobs), pc(maxobs)
      real  sr(maxobs), st(maxobs)

      stname = '     '
      rptTp  = 'LDAD'
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

      numAgr= 0
      numCum= 0
      numShp= 0
      istatusAgr= 0
      istatusCum= 0
      istatusShp= 0

      np= 1
      nq= maxAgr

      call s_len ( inpath, len_inpath )
      inpath= inpath(1:len_inpath)//'agr/'

      call read_agr_cwb (inpath, maxAgr, badflag, i4time_sys,
     ~                   stnTp(np:nq), stname(np:nq), 
     ~                   lats(np:nq), lons(np:nq), elev(np:nq), 
     ~                   t(np:nq), t24max(np:nq), t24min(np:nq), 
     ~                   td(np:nq), rh(np:nq),       
     ~                   pcp1hr(np:nq), pcp3hr(np:nq), 
     ~                   pcp6hr(np:nq), pcp24hr(np:nq),
     ~                   ff(np:nq), wgff(np:nq), wgdd(np:nq), 
     ~                   p(np:nq), sr(np:nq), st(np:nq),
     ~                   numAgr, istatusAgr)

      np= np +numAgr
      nq= np +maxCum
c      inpath= inpath(1:len_inpath)//'cum/'
      inpath= inpath(1:len_inpath)//'new_cum/'
c
c  2013-04 change  cum data format 
c      call read_cum_cwb (inpath, maxCum, badflag, i4time_sys,
      call read_newcum_cwb (inpath, maxCum, badflag, i4time_sys,
     ~                   stnTp(np:nq), stname(np:nq), 
     ~                   lons(np:nq), lats(np:nq), elev(np:nq),
     ~                   pcp1hr(np:nq), pcp3hr(np:nq), 
     ~                   pcp6hr(np:nq), pcp24hr(np:nq),
     ~                   numCum, istatusCum)

      npC= np
      num= numAgr +numCum
      np=  np +numCum
      nq=  np +maxShp
      inpath= inpath(1:len_inpath)//'new_shp/' ! modified by Tin

      call read_shp_cwb (inpath, maxShp, badflag, i4time_sys,
     ~                   stnTp(np:nq), stname(np:nq),
     ~                   lons(np:nq), lats(np:nq), elev(np:nq),
     ~                   p(np:nq),t(np:nq),rh(np:nq),dd(np:nq), 
     ~                   ff(np:nq),numShp, istatusShp)       

      nt= numAgr +numCum +numShp

      do j= np,nt
      do i= npC,num
         if ( stname(i) == stname(j) ) then
            stnTp(i)= stnTp(j)
            p(i)= p(j)
            t(i)= t(j)
            dd(i)= dd(j)
            ff(i)= ff(j)
         endif
      enddo 
      enddo 

      if ( istatusAgr == 1 .and. istatusCum == 1 
     ~                     .and. istatusShp == 1 ) then
         istatus= 1
      else
         istatus= -1
         write(6,*) ' Agr ', istatusAgr, ' Cum ', istatusCum,
     ~              ' Shp ', istatusShp
      endif

      end
 


      subroutine read_agr_cwb (inpath, maxobs, badflag, i4time_sys,
     ~                         stnTp, stname, lats, lons, elev,
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
 
      character(*)  :: stnTp(maxobs), stname(maxobs), inpath
      character(13) :: cvt_i4time_wfo_fname13, a13time_eat,a13time_eat2
      character     :: filename*80, line*180
 
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

      stnTp=   'AGR'
      stname=  '     '
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

      idash= 160
      num=     0
      istatus= 0

      i4time_file_eat= i4time_sys +8*3600             ! convert GMT to EAT
      a13time_eat= cvt_i4time_wfo_fname13(i4time_file_eat)

c
c modified by min-ken.hsieh
c fileanme include mm information
c to prevent ARG data being read at 15,30,45 analysis
c
      filename= 'Data.CWB.AGR.'
     ~           //a13time_eat(1:4) //'-' //a13time_eat(5:6)   ! yyyy_mm
     ~           //'-' //a13time_eat(7:8)                      ! dd
     ~           //'_' //a13time_eat(10:13) //'.agrsv1'      ! hh

      write(6,*) ' agricultural data file ', filename

      call s_len ( inpath, len_inpath )
      call s_len ( filename, len_fname )
 
      open (11,file=inpath(1:len_inpath)//filename(1:len_fname),
     ~         status='old',err=980)

      do i= 1,6
         read (11,*,iostat=istat)
 	 if ( istat < 0 ) then
	    pcp3hr=  badflag
	    pcp6hr=  badflag
	    pcp24hr= badflag
	    write(6,*) ' Warning: an empty agricultural file ',
     ~                 inpath(1:len_inpath)//filename(1:len_fname)
       	    istatus= -1
	    return
	 endif
      enddo
     
      do 450 j= 1,maxobs

         read (11,'(a)',end=500,err=450) line

         call s_len ( line, len_line )
         if ( len_line < 1 ) go to 500

c                  Parse the string into contiguous characters
         istart= 0
         iend=   0
         ivar= 1
         istart(1)= 1
         do i= 1,idash
	    if ( i == 1 )  go to 100

            if (line(i:i) == ' ' .and. line(i-1:i-1) /= ' ') then
               iend(ivar)= i-1
            endif

 100        if (line(i:i) == ' ' .and. line(i+1:i+1) /= ' ') then
               ivar= ivar +1
               istart(ivar)= i+1
            endif
         enddo
         if ( istart(2) == 2 .and. iend(2) == 2 )  exit

         ivar= 1
         read (line(istart(ivar):iend(ivar)),'(a5)',err=399) stname(j)

c              ***** means instruments are maintained
c              ignore all observations of the station at that time
         ivar= 3
         if ( l_parse(line(istart(ivar):iend(ivar)),'N/A') ) then
            t(j)= badflag
         elseif (l_parse(line(istart(ivar):iend(ivar)),'*****')) then
            cycle
         else
            read (line(istart(ivar):iend(ivar)),*,err=399) t(j)
         endif

         ivar= 4
         if ( l_parse(line(istart(ivar):iend(ivar)),'N/A') ) then
            t24max(j)= badflag
         else
            read (line(istart(ivar):iend(ivar)),*,err=399) t24max(j)
         endif

         ivar= 5
         if ( l_parse(line(istart(ivar):iend(ivar)),'N/A') ) then
            t24min(j)= badflag
         else
            read (line(istart(ivar):iend(ivar)),*,err=399) t24min(j)
         endif

         ivar= 6
         if ( l_parse(line(istart(ivar):iend(ivar)),'N/A') ) then
            td(j)= badflag
         else
            read (line(istart(ivar):iend(ivar)),*,err=399) td(j)
         endif

         ivar= 7
         if ( l_parse(line(istart(ivar):iend(ivar)),'N/A') ) then
            rh(j)= badflag
         else
            read (line(istart(ivar):iend(ivar)),*,err=399) rh(j)
         endif

         ivar= 8
         if ( l_parse(line(istart(ivar):iend(ivar)),'N/A') ) then
            pcp1hr(j)= badflag
         else
            read (line(istart(ivar):iend(ivar)),*,err=399) pcp1hr(j)
         endif

         ivar= 10
         if ( l_parse(line(istart(ivar):iend(ivar)),'N/A') ) then
            ff(j)= badflag
         else
            read (line(istart(ivar):iend(ivar)),*,err=399) ff(j)
         endif

         ivar= 11
         if ( l_parse(line(istart(ivar):iend(ivar)),'N/A') ) then
            wgff(j)= badflag
         else
            read (line(istart(ivar):iend(ivar)),*,err=399) wgff(j)
         endif

         ivar= 12
         if ( l_parse(line(istart(ivar):iend(ivar)),'N/A') ) then
            wgdd(j)= badflag
         else
            read (line(istart(ivar):iend(ivar)),*,err=399) wgdd(j)
         endif

         ivar= 14
         if ( l_parse(line(istart(ivar):iend(ivar)),'N/A') ) then
            sr(j)= badflag
         else
            read (line(istart(ivar):iend(ivar)),*,err=399) sr(j)
         endif

         ivar= 17
         if ( l_parse(line(istart(ivar):iend(ivar)),'N/A') ) then
            p(j)= badflag
         else
            read (line(istart(ivar):iend(ivar)),*,err=399) p(j)
         endif

         ivar= 18
         if ( l_parse(line(istart(ivar):iend(ivar)),'N/A') ) then
            st(j)= badflag
         else
            read (line(istart(ivar):iend(ivar)),*,err=399) st(j)
         endif

         go to 400

 399     write(6,*) ' read error in station/variable ', num+1, ivar,
     ~              stname(j), ' ', line(istart(ivar):iend(ivar))
         cycle

 400     num= num +1
 450  enddo
 
c                      obtain Tmax and Tmin in 24 hours
 500  call agr_t24_pcp (inpath, filename, stname, maxobs, badflag,
     ~     num, t24max, t24min, pcp1hr, pcp3hr, pcp6hr, pcp24hr,
     ~     istatus)

      if ( istatus /= 1 ) then
         write(6,*) ' Error executing agr_t24_pcp '
         return
      endif

c                            quality control
      if ( num > maxobs ) then
         write(6,*) ' read_agr_cwb error for too many obs: ',
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
            ff(j)= ff(j) * 1.94254         ! conv m/s to kt
         endif
 
         if ( wgff(j) < 0 ) then
            wgff(j)= badflag
         else
            wgff(j)= wgff(j) * 1.94254     ! conv m/s to kt
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
      write(6,*) ' Found ', num, ' agricultural stations.'
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
 
      real         :: lat(maxobs), lon(maxobs), elev(maxobs)
      character(5) :: stn_name(maxobs), stn_name_in
      character(*) :: inpath
 
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

 

      subroutine agr_t24_pcp (inpath, filename, stname, maxobs, badflag,
     ~           num, t24max, t24min, pcp1hr, pcp3hr, pcp6hr, pcp24hr,
     ~           istatus)

      integer, parameter :: num23= 23 
 
      character(*) :: stname(*), inpath
      character(5) :: stn(num23,maxobs)
      character(2) :: yy, mm, dd, hh
      character    :: filename*35, fileDummy*35, line*180
 
      logical :: l_parse
      integer :: istart(30), iend(30), d(12), count(num23), flag
 
      real :: t24max(maxobs), t24min(maxobs), pcp1hr(maxobs)
      real :: pcp3hr(maxobs), pcp6hr(maxobs), pcp24hr(maxobs)
      real :: tmax(num23,maxobs), tmin(num23,maxobs), p1hr(num23,maxobs)
 
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
      istatus= 0

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

         filename= inpath(1:len_inpath)//fileDummy(1:len_fname)
         open (24,file=inpath(1:len_inpath)//fileDummy(1:len_fname),
     ~            status='old',err=980)

         do i= 1,6
            read (24,*,iostat=istat)
            if ( istat < 0 ) then
               pcp3hr=  badflag
               pcp6hr=  badflag
               pcp24hr= badflag
               write(6,*) ' Warning: an empty agricultural file ',
     ~                    inpath(1:len_inpath)//fileDummy(1:len_fname)
               istatus= -1
               return
            endif
         enddo

         do 200 j= 1,maxobs

            read (24,'(a)',iostat=istat) line
            if     ( istat < 0 ) then
               cycle
            elseif ( istat > 0 ) then
               write(*,*) 'Warning: read error in ',
     ~                    inpath(1:len_inpath)//fileDummy(1:len_fname)  
            endif

c                 Parse the string into contiguous characters
            istart= 0
            iend=   0
            ivar= 1
            istart(1)= 1

            do i= 1,idash
               if ( i == 1 )  go to 100

               if ( line(i:i) == ' ' .and. line(i-1:i-1) /= ' ' ) then      
                  iend(ivar)= i-1
               endif

 100           if ( line(i:i) == ' ' .and. line(i+1:i+1) /= ' ' ) then
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
	 write(6,*) fileDummy(1:len_fname), ' ', stn(k,j)
         write(6,*) ivar,line(istart(ivar):iend(ivar))
           
      enddo k_num23

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
        


      subroutine read_cum_cwb (inpath, maxobs, badflag, i4time_sys,
     ~                         stnTp, stname, lats, lons, elev, 
     ~                         pcp1hr, pcp3hr, pcp6hr, pcp24hr,      ! O
     ~                         num, istatus)                         ! O
c
c======================================================================
c
c     Routine to read the CWB ASCII rain gauge files.
c     
c======================================================================
c
      integer, parameter :: loopnum = 800 
      integer, parameter :: num23= 23

      real :: lats(maxobs), lons(maxobs), elev(maxobs)
      real :: pcp1hr(maxobs), pcp3hr(maxobs), pcp6hr(maxobs)
      real :: pcp24hr(maxobs), p1hr(maxobs,num23)

      integer :: dataMin(maxobs), dMin(maxobs,num23), d(12), flag
 
      character(*)  :: stnTp(maxobs), stname(maxobs), inpath
      character(80) :: filename
      character(13) :: cvt_i4time_wfo_fname13, a13time_sys
      character(4)  :: stn(maxobs,num23)
      character(2)  :: yy, mm, dd, hh
      character(1)  :: pcp1hrQ(maxobs), pcp3hrQ(maxobs), pcp6hrQ(maxobs)
      character(1)  :: p1hrQ(maxobs,num23)

c                    Stuff for the rain gauge metadata.
      real  elev_master(maxobs)
 
      character(5) :: stn_name_master(maxobs)
      character(4) :: stn_id_master(maxobs), stn_id(maxobs)

      data  d / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /
 
c             Get the rain gauge metadata (station information).
      call read_cum_stntbl (inpath, maxobs, badflag,
     ~                      stn_id_master, stn_name_master,
     ~                      elev_master, num_master, istatus)
      if ( istatus /= 1 ) then
         write(6,*) ' Error reading rain gauge station table'
         return
      endif

c   fill the output arrays with something, then open the file to read
      istatus= 0
      flag=    0
      stnTp=   'CUM'
      stn=     '     '
      stname=  '     '
      p1hrQ=   ' '
      pcp1hrQ= ' '
      pcp3hrQ= ' '
      pcp6hrQ= ' '
      p1hr=    badflag
      pcp1hr=  badflag
      pcp3hr=  badflag
      pcp6hr=  badflag
      pcp24hr= badflag

      a13time_sys= cvt_i4time_wfo_fname13(i4time_sys)

      read (a13time_sys( 1: 4),'(i4)') iy
      read (a13time_sys( 5: 6),'(i2)') im
      read (a13time_sys( 7: 8),'(i2)') id
      read (a13time_sys(10:11),'(i2)') ih

      j_num23 : do j= 0,num23
         if ( j > 0 )  ih= ih -1

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
         iy= iy +2000

         filename= 'Data.CWB.CUM_UDD.' //'20' //yy //'-' //mm //'-' 
     ~                                 //dd //'_' //hh //'.cum_hr'      
         call s_len ( inpath, len_inpath )
         call s_len ( filename, len_fname )
 
         open (11,file=inpath(1:len_inpath)//filename(1:len_fname),
     ~            status='old',iostat=istat)
         if ( istat > 0 ) then
            if ( j == 0 ) then
               go to 980
            else
               flag= 1                    ! fail to open rain gauge files
               go to 600
            endif
         endif

         if ( j == 0 ) then
            num= 0

            do i= 1,loopNum
               read (11,501,iostat=istat) stn_id(i), lats(i), lons(i),      
     ~                      pcp1hr(i), pcp3hr(i), pcp6hr(i), dataMin(i),
     ~                      pcp1hrQ(i), pcp3hrQ(i), pcp6hrQ(i)       
 501           format (a4, 11x, 2f8.4, 1x, 3f5.1, i5, 3a1)  

               if ( istat == -1 .or. istat == -2 )  exit
               if ( istat > 0 .or. pcp1hr(i) == badflag ) then
                  pcp1hr= badflag
                  pcp3hr= badflag
                  pcp6hr= badflag
                  go to 990
               endif

               num= num +1
            enddo
         else
            do i= 1,num
               read (11,511,iostat=istat) stn(i,j), p1hr(i,j), 
     ~                                    dMin(i,j), p1hrQ(i,j)
 511           format (a4, 28x, f5.1, 10x, i5, a1)  

               if ( istat == -1 .or. istat == -2 )  exit
               if ( istat > 0 .or. p1hr(1,j) == badflag ) then
                  pcp24hr= badflag
                  go to 990
               endif
            enddo
         endif
      enddo j_num23

c  Match data with metadata for this station, then store the metadata in arrays.
 600  imatch= 0
 
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
c              pcp1hr(i)= pcp1hr(i) *0.1 *0.03937         ! conv mm to inch
c modified by min-ken hsieh
c mm to inch coeff. changed
               pcp1hr(i)= pcp1hr(i) *0.03937         ! conv mm to inch
            endif
 
            if ( pcp3hr(i) < 0 .or. pcp3hrQ(i) /= 'G' ) then
               pcp3hr(i)= badflag
            else
c              pcp3hr(i)= pcp3hr(i) *0.1 *0.03937         ! conv mm to inch
               pcp3hr(i)= pcp3hr(i) *0.03937         ! conv mm to inch
            endif
 
            if ( pcp6hr(i) < 0 .or. pcp6hrQ(i) /= 'G' ) then
               pcp6hr(i)= badflag
            else
c              pcp6hr(i)= pcp6hr(i) *0.1 *0.03937         ! conv mm to inch
               pcp6hr(i)= pcp6hr(i) *0.03937         ! conv mm to inch
            endif
         endif
      enddo 

c               there is no need to calculate 24 hr rain fall any more 
c                     when failing to open rain gauge file 
      if ( flag == 1 )  go to 980

      pcp24hr= pcp1hr

c                      quality control of 24hr precipitation
      pcp_stations : do i= 1,num
      if (pcp24hr(i) .ne. badflag) then	! Check initial 1 hour pcp before accum.
      pcp_hours    : do j= 1,num23
         if ( stn_id(i) == stn(i,j) .and. dMin(i,j) <= 60 .and.
     ~        p1hr(i,j) >= 0        .and. p1hrQ(i,j) == 'G' ) then
c           p1hr(i,j)= p1hr(i,j) *0.1 *0.03937         ! conv mm to inch
            p1hr(i,j)= p1hr(i,j) *0.03937         ! conv mm to inch
            pcp24hr(i)= pcp24hr(i) +p1hr(i,j)    
         else
            p1hr(i,j)=  badflag
            pcp24hr(i)= badflag
            cycle pcp_stations
         endif
      enddo pcp_hours
      endif
      enddo pcp_stations

c                               end of file
      write(6,*) ' Found ', num, ' rain gauge stations.'
      istatus= 1
      return
 
 980  write(6,*) ' Warning: could not open rain gauge file ', filename
      istatus= -1
      return
      
 990  write(6,*) ' ** ERROR reading rain gauge data.'
      num= 0
      istatus= -1
      return

      end
 
 
 
      subroutine read_cum_stntbl (inpath, maxobs, badflag, stn_id,
     ~                           stn_name, elev, num, istatus)       
 
c======================================================================
c
c     Routine to read station information for the CWB ASCII rain gauge
c	station data.
c
c======================================================================
 
      real  elev(maxobs)
 
      character(*) :: inpath
      character(4) :: stn_id(maxobs), stn_id_in
      character(5) :: stn_name(maxobs), stn_name_in
      integer      :: elev_m
 
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
 550  write(6,*) ' Found ', num,
     ~           ' rain gauge stations in the station table.' 
      istatus= 1
      return
      
 990  write(6,*) ' ** ERROR reading rain gauge station table'
      write(6,*) stn_id_in, stn_name_in, elev_m
      istatus= 0
      return

      end
 


      subroutine read_shp_cwb (inpath, maxobs, badflag, i4time_sys,
     ~                         stnTp, stname, lats, lons, elev,       ! O
     ~                         p, t,rh, dd, ff, num, istatus)            ! O
 
c======================================================================
c
c     Routine to read the CWB ASCII shp files.
c     
c======================================================================
 
      integer, parameter :: loopnum= 250
      integer :: imatch,i_invalid
      real lats(maxobs), lons(maxobs), elev(maxobs)
      real p(maxobs),t(maxobs),dd(maxobs),ff(maxobs),rh(maxobs)

      character(*)  :: stnTp(maxobs), stname(maxobs), inpath
      character(80) :: filename
      character(13) :: cvt_i4time_wfo_fname13, a13time_eat,a13time_eat2

c                         stuff for the shp metadata.
      real  lat_master(maxobs), lon_master(maxobs), elev_master(maxobs)       
      character(5) :: stn_name_master(maxobs)
      character(4) :: stn_id_master(maxobs), stn_id(maxobs)
      character(4) :: hhmm(maxobs)
 
c                    Get the shp metadata (station information).

c     ===== Tin added =====
      write(*,*) "***** TIN TEST:START READ SHP DATA *****"

!      call read_shp_stntbl (inpath, maxobs, badflag,
!     ~                      stn_id_master, stn_name_master,
!     ~                      lat_master, lon_master, elev_master,
!     ~                      num_master, istatus)
!      if ( istatus /= 1 ) then
!         write(6,*)' Error reading shp station table'
!         return
!      endif

c   fill the output arrays with something, then open the file to read
 
      istatus= 0
      stnTp=    'SHP'
      stname=   '     '
      p=  badflag
      t=  badflag
      dd= badflag
      ff= badflag

      i4time_file_eat= i4time_sys +8*3600             ! convert GMT to EAT

      a13time_eat= cvt_i4time_wfo_fname13(i4time_file_eat)
      
      a13time_eat2 = cvt_i4time_wfo_fname13(i4time_sys) ! added by Tin (Use UTC)
c
c modified by min-ken.hsieh
c filename include hhmm information to get each 15 min data
c
!      filename= 'Data.CWB.SHP.'
!     ~           //a13time_eat(1:4)//'-'//a13time_eat(5:6)   ! yyyy_mm
!     ~           //'-'//a13time_eat(7:8)                     ! dd
!     ~           //'_'//a13time_eat(10:13) //'.shp'        ! hh

      filename= a13time_eat2(1:8)//a13time_eat2(10:13)//'.' ! yyyymmddhhMM
     ~           //'QPESUMS_STATION.15M.mdf'

      write(6,*) ' SHP: file ',filename

      call s_len ( inpath, len_inpath )
      call s_len ( filename, len_fname )
 
      num= 0

      open (11,file=inpath(1:len_inpath)//filename(1:len_fname),
     ~         status='old',err=980)
      do j=1,3
         read(11,*)
      end do

      do i= 1,loopNum
!         read (11,511,end=600,err=990) stn_id(i), 
!     ~                                 hhmm(i),p(i), t(i), ff(i), dd(i) 
! 511     format (5x, a4, 9x,a4,1x, 4f6.1)
c modified by Tin
         read (11,511,end=600,err=990) stname(i),lats(i),
     ~                                 lons(i),elev(i),dd(i),
     ~                                 ff(i),t(i),rh(i),p(i),hhmm(i)
!HJ W>=D+3 change f4.2 to f5.2
 511     format(a6,21x,f7.4,2x,f8.4,2x,f6.1,2x,f6.1,
     ~          1x,f6.1,1x,f6.1,4x,f5.2,2x,f6.1,90x,a4)
         num= num +1
      enddo

c  Match data with metadata for this station, then store the metadata in arrays.
 600  imatch= 0
c
c modified by min-ken.hseih
c filter out data which labeled wrong time
c and find out how many stns match table
c 
      i_match = 0 
      i_invalid = 0
      do i= 1,num
         if (hhmm(i)/=a13time_eat(10:13)) then
	       i_invalid = i_invalid+1
	       p(i)= badflag
	       t(i)= badflag
               dd(i)= badflag
               ff(i)= badflag
               rh(i)= badflag
         endif

!         do j= 1,num_master
!            if ( stn_id(i) == stn_id_master(j) ) then
!               stname(i)= stn_name_master(i)
!               lats(i)= lat_master(j)
!               lons(i)= lon_master(j)
!               elev(i)= elev_master(j)
!               imatch=1
!               i_match = i_match+1
!            endif
!         enddo

!         if ( imatch == 0 ) then
!            write(6,*)' No station match ', stn_id(i), stn_id_master(j)       
!         endif
 
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
            ff(i)= ff(i) * 1.94254        ! conv m/s to kt
         endif
          rh(i)=rh(i)*100.
         if ( rh(i) < 0  )  rh(i)= badflag
         if ( rh(i) > 100.  )  rh(i)=  100.
      enddo
c                               end of file
      write(6,*) ' Total ', num, ' shp stations.'
      write(6,*) ' Found ', num - i_invalid, ' shp stations.'
!      write(6,*) ' Found ', i_match - i_invalid,
!     +			 ' shp stations match the table.'
      istatus= 1
c     ===== Tin added =====
!      write(*,*) dd
!      write(*,*) i_match
      write(*,*) "***** TIN TEST:END READ SHP DATA *****"
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
 
      real :: lat(maxobs), lon(maxobs), elev(maxobs)
 
      character(*) :: inpath
      character(4) :: stn_id(maxobs), stn_id_in
      character(5) :: stn_name(maxobs), stn_name_in
      integer      :: elev_m
 
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
 550  write(6,*) ' Found ', num, ' shp stations in the station table.'    
      istatus= 1
      return
      
 990  write(6,*) ' ** ERROR reading shp station table'
      istatus= 0
      return

      end




      subroutine get_box_size (box_size,istatus)

      istatus = 1
      box_size = 1.1

      return
      end


      subroutine read_newcum_cwb (inpath, maxobs, badflag, i4time_sys,
     ~                         stnTp, stname, lons, lats, elev, 
     ~                         pcp1hr, pcp3hr, pcp6hr, pcp24hr,      ! O
     ~                         num, istatus)                         ! O
c
c======================================================================
c
c     Routine to read the CWB ASCII rain gauge files.
c     
c======================================================================
c
      integer, parameter :: loopnum = 800 
chj num23 is not defined in this subroution. Added. HJ 10/14/2013
      integer, parameter :: num23= 23

      real :: lats(maxobs), lons(maxobs), elev(maxobs)
      real :: pcp1hr(maxobs), pcp3hr(maxobs), pcp6hr(maxobs)
      real :: pcp10min(maxobs),pcp12hr(maxobs),pcp24hr(maxobs)

      integer :: dataMin(maxobs), dMin(maxobs,num23), d(12), flag
 
      character(*)  :: stnTp(maxobs), stname(maxobs), inpath
      character(80) :: filename
      character(13) :: cvt_i4time_wfo_fname13, a13time_sys
      character(2)  :: yy, mm, dd, hh

c                    Stuff for the rain gauge metadata.

      data  d / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /
 
c   fill the output arrays with something, then open the file to read
      istatus= 0
      flag=    0
      stnTp=   'CUM'
      stn=     '     '
      stname=  '     '
      p1hrQ=   ' '
      pcp1hrQ= ' '
      pcp3hrQ= ' '
      pcp6hrQ= ' '
      p1hr=    badflag
      pcp1hr=  badflag
      pcp3hr=  badflag
      pcp6hr=  badflag
      pcp24hr= badflag

      a13time_sys= cvt_i4time_wfo_fname13(i4time_sys)

      read (a13time_sys( 1: 4),'(i4)') iy
      read (a13time_sys( 5: 6),'(i2)') im
      read (a13time_sys( 7: 8),'(i2)') id
      read (a13time_sys(10:11),'(i2)') ih

         iy= iy -2000
         call i2a ( iy, yy )
         call i2a ( im, mm )
         call i2a ( id, dd )
         call i2a ( ih, hh )
         iy= iy +2000

         filename= '20' //yy //mm //dd//hh //'00.' 
     ~                                 //'QPESUMS_GAUGE.10M.mdf'      
c         filename= 'Data.CWB.CUM_UDD.' //'20' //yy //'-' //mm //'-' 
c     ~                                 //dd //'_' //hh //'.cum_hr'      
         call s_len ( inpath, len_inpath )
         call s_len ( filename, len_fname )
 


         open (11,file=inpath(1:len_inpath)//filename(1:len_fname),
     ~            status='old',err=980,iostat=istat)

         do j=1,3
           read(11,*)
         end do
 
            num=0
         do i= 1,loopnum
               read (11,501,iostat=istat)stname(i),lats(i),lons(i),      
     ~                      elev(i),pcp1hr(i),pcp10min(i),pcp3hr(i),
     ~                      pcp6hr(i),pcp12hr(i),pcp24hr(i)       
c 501           format (8x,a8,9x,f7.4,2x,f8.4,f8.2,6f8.2)  
 501           format (a6,19x,f7.4,2x,f8.4,f8.2,6f8.2)  

               if ( istat == -1 .or. istat == -2 )  exit
               if ( istat > 0 .or. pcp1hr(i) == badflag ) then
                  pcp1hr= badflag
                  pcp3hr= badflag
                  pcp6hr= badflag
                  pcp12hr= badflag
                  pcp24hr= badflag
                  go to 990
               endif

               num= num +1
c                               quality control
            if ( pcp1hr(i) < 0  ) then
               pcp1hr(i)= badflag
            else
               pcp1hr(i)= pcp1hr(i) *0.03937         ! conv mm to inch
            endif
 
            if ( pcp3hr(i) < 0 ) then
               pcp3hr(i)= badflag
            else
               pcp3hr(i)= pcp3hr(i) *0.03937         ! conv mm to inch
            endif
 
            if ( pcp6hr(i) < 0 ) then
               pcp6hr(i)= badflag
            else
               pcp6hr(i)= pcp6hr(i) *0.03937         ! conv mm to inch
            endif

            if ( pcp12hr(i) < 0 ) then
               pcp12hr(i)= badflag
            else
               pcp12hr(i)= pcp6hr(i) *0.03937         ! conv mm to inch
            endif

            if ( pcp24hr(i) < 0 ) then
               pcp24hr(i)= badflag
            else
               pcp24hr(i)= pcp24hr(i) *0.03937         ! conv mm to inch
            endif
      enddo 

c                               end of file
      write(6,*) ' Found ', num, ' rain gauge stations.'
      istatus= 1
      return
 
 980  write(6,*) ' Warning: could not open rain gauge file ', filename
      istatus= -1
      return
      
 990  write(6,*) ' ** ERROR reading rain gauge data.'
      num= 0
      istatus= -1
      return

      end
