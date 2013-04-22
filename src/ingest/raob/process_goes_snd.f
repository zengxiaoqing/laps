
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


      subroutine process_goes_snd (path, path_len, filename, 
     1     file_len,iii,jjj,
     1     i4time_begin,i4time_end, mdf, lun_out, istatus)

c
c	modified for filename change 4/22/2013 DB
c
      implicit none

c     parameter list
      integer istatus
      character*(*) path, filename
      real mdf                  !missing data flag
      integer path_len, file_len, i4time_begin,i4time_end,lun_out
      integer iii,jjj

c     variables for output writing
      integer maxsnd,maxlvl
      parameter (maxsnd = 1)
      parameter (maxlvl = 50)

      integer iwmostanum(maxsnd), nlvl (maxsnd)
      real stalat(maxsnd,maxlvl),stalon(maxsnd,maxlvl),staelev(maxsnd)
      character c5_staid(maxsnd)*5,a9time_ob(maxsnd,maxlvl)*9,
     1     c8_obstype(maxsnd)*8
      real height_m(maxsnd,maxlvl)
      real pressure_mb(maxsnd,maxlvl)
      real temp_c(maxsnd,maxlvl)
      real dewpoint_c(maxsnd,maxlvl)
      real dir_deg(maxsnd,maxlvl)
      real spd_mps(maxsnd,maxlvl)
      integer nsnd ! number of soundings per write call (1)
      integer count_level,ob_counter,id_sat



c     internal variables
      real lat,lon
      integer i4time_record,nbuf_1,nbuf_2
      integer avg_i4time_record ! to create modified time rec
      character*9 c_time_record
      character*9 m_time_record ! modified time record
      character*11 c_time_record_long
      real lat_a(iii,jjj), lon_a(iii,jjj), topo_a (iii,jjj)
      real rnorth,south,east,west

C ?
C ?
C ?  *****  DATA RECORD FORMAT  *****
C ?
C ?    -----  HEADER RECORD (RECORD 1)  -----
C ?
C ?   WORD       DESCRIPTION                         SCALE
C ?  -----       -----------                         -----
C ?    1         DAY (YYDDD)                         1
C ?    2         TIME (HHMMSS)                       1
C ?    3         NUMBER OF RECORDS IN FILE           1
C ?    4         SATID                               1
C ?    5         LAT/LON OF NW CORNER OF DATA        1
C ?    6         LAT/LON OF SE CORNER OF DATA        1
C ?
C ?   ------   DATA (RECORDS 2-N)  ------
C ?   
C ?   WORD       DESCRIPTION                         SCALE
C ?  -----       -----------                         -----
C ?    1         DAY (YYDDD)                         1
C ?    2         TIME (HHMMSS)                       1
C ?    3         SATID                               1
C ?    4         MOD FLAG                            1
C ?    5         LATITUDE (DEG)                      100 
C ?    6         LONGITUDE (DEG)                     100
C ?    7         SAM (#FOVs)                         1
C ?    8         CAPE INDEX                          1
C ?    9-26      OBSERVED BT (K)                     100
C ?    27        SOLAR ZENITH ANGLE (DEG)            100
C ?    28        LOCAL ZENITH ANGLE (DEG)            100
C ?    29        RETRIEVAL TYPE                      1
C ?    30        SKIN TEMP (K)                       100
C ?    31        TOTAL PRECIP WATER (MM)             100
C ?    32        LIFTED INDEX (K)                    100
C ?    33        GUESS TOTAL PRECIP WATER (K)        100
C ?    34        GUESS LIFTED INDEX (K)              100
C ?    35        LAYER PRECIP. WATER (1-.9 SIGMA)    100
C ?    36        LAYER PRECIP. WATER (.9-.7 SIGMA)   100
C ?    37        LAYER PRECIP. WATER (.7-.3 SIGMA)   100
C ?    38        SFC TEMP (K) @ SFC PRESSURE         100
C ?    39-78     TEMP PROFILE (K) (1000-.1MB)        100
C ?    79        SFC DEWPOINT TEMP (K) @ SFC PRESS   100
C ?    80-119    DEWPOINT TEMP (K) (1000-0.1MB)      100
C ?    120       SFC GEOPOTENTIAL HEIGHT (M)         1
C ?    121-160   GEOPOTENTIAL HEIGHT (M)             1
C ?    161       SFC PRESUURE (MB)                   10
C ?    162-201   PRESSURE (1000-0.1MB)               10
C ?    202-210   SPARES                              1
C ?
C ?   ###########################################################

      character*1 csat
      character*5 cyyjd
      character*6 cvar(37)
      character*4 chhmm
      CHARACTER*80 CTEXT
      INTEGER*4 NBUF(210),ISCALE(210)
      REAL*4 RBUF(210)

      DATA ISCALE /4*1,2*100,2*1,20*100,1,90*100,41*1,41*10,9*1/
      data cvar /'DAY   ','TIME  ','SATID ','MOD   ','LAT   ',
     1'LON   ','SAM   ','CA    ','BT01  ','BT02  ','BT03  ','BT04  ',
     2'BT05  ','BT06  ','BT07  ','BT08  ','BT09  ','BT10  ','BT11  ',
     3'BT12  ','BT13  ','BT14  ','BT15  ','BT16  ','BT17  ','BT18  ',
     4'SZA   ','LZA   ','RT    ','TSKIN ','TPW   ','LI    ','GTPW  ',
     5'GLI   ','WV1   ','WV2   ','WV3   '/

      EQUIVALENCE(NBUF(1),RBUF(1))

      integer ilev,jj,js,k,lt,lp,lz,ltd,irec,iend,ibegin,nrec,mm,lll
      integer ilin


c-------------------------------------------------------------------------------------
C

      istatus = 0 ! set initially as bad

c     compute average record time for ingest
      avg_i4time_record = nint((float(i4time_begin)+
     1     float(i4time_end))/2.)
      call make_fnam_lp(avg_i4time_record,m_time_record,istatus)

c     get perimeter values for the domain
      call get_domain_perimeter(iii,jjj,'nest7grid',lat_a,
     1     lon_a, topo_a, 1.0, rnorth,south,east,west,istatus)

c     open file 


       OPEN(70,FILE=path(1:path_len)//filename(1:file_len)
     1     ,ACCESS='DIRECT',RECL=840,
     1 FORM='UNFORMATTED',STATUS='OLD',ERR=9000)
C
      irec = 0
c     read header information

      READ(70,REC=1,ERR=9999) NBUF

C     endian swap if needed
      if (filename(5:6) .lt. '09') then
      do k  = 1, 210
         call endian4 (nbuf(k))
      enddo
      endif


      nrec = nbuf(1)
      write(6,*) 'Header Record Contents:'
      write(6,*) '-----------------------'
      write(6,*) 'Day  of data (YYDDD)     = ',nbuf(1)
      write(6,*) 'Time of data (HHMM)      = ',nbuf(2),'Z'
      write(6,*) 'Number of Records        = ',nbuf(3)
      write(6,*) 'Satellite ID             = ',nbuf(4)
      write(6,*) 'Lat/Lon of NW Corner     = ',nbuf(5)
      write(6,*) 'Lat/Lon of SE Corner     = ',nbuf(6)
      write(6,*)

      id_sat = nbuf(4)

      ibegin = 2
      iend = nbuf(3)  ! number of records

      DO 200 ILIN=ibegin,iend      
        DO 205 LLL=1,210
           NBUF(LLL) = 0
  205   CONTINUE
        READ(70,REC=ILIN,ERR=9999) NBUF 

c     Endian swap if needed
      if (filename(5:6) .lt. '09') then
        do k = 1, 210
           call endian4 (nbuf(k))
        enddo
      endif


        DO 210 MM = 1,210
           RBUF(MM) = (NBUF(MM)/FLOAT(ISCALE(MM)))
  210   CONTINUE
        irec = irec + 1
c        write(6,1000) ilin 
 1000   format(20X,'***  RECORD ',i6,' CONTENTS  ***' /)
c        write(6,1010) (cvar(k),rbuf(k),k=1,4) 
 1010   format(1x,4(a6,' = ',f7.0,1x,'|',2x),/)
c        write(6,1015) (cvar(k),rbuf(k),k=5,7) 
 1015   format(1x,a6,' = ',f6.2,2x,'|',2x,a6,' = ',
     1  f6.2,2x,'|',2x,a6,' = ',f6.0,2x,'|',/)
c        write(6,1020) (cvar(k),rbuf(k),k=8,25)
 1020   format(4(1x,a6,' = ',f8.2,2x,'|',1x),/)
c        write(6,1030) (cvar(k),rbuf(k),k=26,37)
 1030   format(/,4(1x,a6,' = ',f6.2,2x,'|',1x),/)
c        write (6,*) rbuf(37)+rbuf(36)+rbuf(35)
c        write(6,1040) 
 1040   format(/,1x,'PRESSURE',4x,'   T    ',4x,'   TD   ',4x,
     1  '   Z    ')
c        write(6,1050)
 1050   format(1x,'--------',4x,'--------',4x,'--------',4x,'--------')


c     set up time window compare
        nbuf_1=int(rbuf(1))
        nbuf_2=int(rbuf(2))
        write (c_time_record_long,34) nbuf_1,nbuf_2
        c_time_record = c_time_record_long(1:9)
 34     format (i5.5,i6.6)
        call i4time_fname_lp(c_time_record,i4time_record,istatus)
        if (i4time_record.ge.i4time_begin 
     1       .and. 
     1       i4time_record.le.i4time_end) then
           continue
        else
           write (6,*) 'skipping record time bounds ',c_time_record
           go to 200 ! bail on this record
        endif


c     Extract lat and lon

        lat = rbuf(5)
        lon = -1.*rbuf(6) ! our convention is east longitude

c     Test whether point is in the domain

        if (lat.le.rnorth .and. lat .ge. south
     1       .and. lon .ge. west .and. lon .le. east )then
           write (6,*) 'accepting sounding ', lat,lon
           write (6,*) 'routine process_goes_snd'
           continue
        else
           ! outside the perimeter - reject
           write (6,*) 'lat lon reject', lat, lon
           go to 200            !bail on this record

        endif                   !domain test finished

c     good ob
        ob_counter=ob_counter+1

        count_level = 0
        js = 1
        ilev = 0
        write (6,*) lat,lon
        do 170 jj=38,78
          ilev = ilev + 1
          lt   = jj
          ltd  = lt + 41
          lz   = ltd + 41
          lp   = lz + 41
          if(rbuf(lt) .gt. -9.9e+3) then
c            write(6,1060) rbuf(lp),rbuf(lt)-273.15,
c     1            rbuf(ltd)-273.15  ,rbuf(lz)
            count_level=count_level+1

c     Fill variables for calling the routine to write SND
            nsnd = 1
            iwmostanum(1) = 0
            stalat(1,count_level) = lat
            stalon(1,count_level) = lon
            write(c5_staid(1)(2:5),44) ob_counter
 44         format(i4.4)
            c5_staid(1)(1:1) = 'S'
c            a9time_ob (1,count_level) = c_time_record
            a9time_ob (1,count_level) = m_time_record
            c8_obstype(1) = 'GOES    '
            write (c8_obstype(1)(5:6),45)id_sat
 45         format(i2.2)
            height_m (1,count_level) = rbuf(lz)
            pressure_mb (1,count_level) = rbuf(lp)
            temp_c (1,count_level) = rbuf(lt)-273.15
            dewpoint_c (1,count_level) = rbuf(ltd)-273.15
            dir_deg (1,count_level) = mdf
            spd_mps (1,count_level) = mdf
c            staelev(1) = height_m(1,1)
            staelev(1)  = 0.0
          endif
 1060     format(1x,4(f8.2,4x))
  170   continue

c     finish off with stuff for writing 
        nlvl(1) = count_level ! number of levels

c     Call write routine
        write(6,*) 'including satsnd in output file'
        call write_snd (lun_out,
     1       maxsnd,maxlvl,nsnd,
     1       iwmostanum,
     1       stalat,stalon,staelev,
     1       c5_staid,a9time_ob,c8_obstype,
     1       nlvl,
     1       height_m,
     1       pressure_mb,
     1       temp_c,
     1       dewpoint_c,
     1       dir_deg,
     1       spd_mps,
     1       istatus)
C
  200 continue
      go to 9001
 9000 continue
      write (6,*) 'routine process_snd failed on open, 9000'

c     call sdest('open error',0)
C
 9999 CONTINUE
      write (6,*) 'routine process_snd failed on read, 9999'
      istatus   = 0
      CLOSE(70)
      return

 9001 Continue
      write (6,*) 'Routine Process_SND Success on read'
      istatus = 1
      close (70)
      RETURN
      END


      subroutine endian4(byte)
c
      integer*1 byte(4),tmp(4)
c
      do j = 1, 4
        tmp(j) = byte(j)
      enddo  !  j
c
      do j = 1, 4
        j1 = 5-j
        byte(j1) = tmp(j)
      enddo  !  j                      
c
      return
      end

