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
cdis cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
        subroutine get_rad (i4time,pw,plat,plon,npts,istat)

c       $log: get_rad.for,v $
c revision 1.5  1995/09/13  21:35:20  birk
c added disclaimer to files
c
c revision 1.4  1994/09/13  16:32:43  birk
c corrected minor mixed mode problem that didn't port to ibm
c
c revision 1.3  1994/09/13  15:32:03  birk
c updated routine to deal with radiometer data format change
c
c revision 1.2  1994/04/27  16:46:12  birk
c modified routine to fail only in the event no radiometers are avail
c to date it has been failing if denver was not available
c this mod will make the routine more robust
c
c revision 1.1  1994/04/25  15:19:24  birk
c initial revision
c

c       8 jan 1991   -   d. birkenheuer

c this routine is for obtaining the radiometer data from wpl for inclusion
c in the water vapor routines.  it is anticipated that this is only atemporary
c routine at this time. (in order to get through wisp)


c added new data access routines from r. beeler and k. plattner 27 jan 1992
c

        implicit none

        real
     1  pw(4),
     1  plat(4),
     1  plon(4)
        integer
     1  i4time, !note this is satellite i4time
     1  npts,  !npfilepts in calling routine
     1  istat,
     1  istatus


        integer
     1  i4time_d, !desired i4time of the data in question
     1  i4time_f, !computed i4time of profiler data file to open
     1  twenty_seg !the twenty minute segment after the
c                             hour (0,1,2)
        integer
     1  ny_d,
     1  nmo_d,
     1  nd_d,
     1  nh_d,
     1  nm_d,
     1  ns_d

        real
     1  hour_fraction !fractional hour used in determining i4time_f

        character*9
     1  filename
        integer
     1  i,j,k,l,ii,jj,kk,ll
        real
     1  data(10,4)
        integer
     1  weight(5)

        data weight/ -3,12,17,12,-3/

        character*80 header
        character*6 site, type


        integer
     1  iunit,
     1  ierror,
     1  points,
     1  year(30),
     1  month(30),
     1  day(30),
     1  hour(30),
     1  minute(30),
     1  second(30),
     1  tav(30),
     1  pm(30),
     1  md(30),
     1  fh(30),
     1  qc(30)

        real
     1  elevation,
     1  longitude,
     1  latitude,
     1  mrelev(30),
     1  mrazim(30),
     1  vapor(30),
     1  liquid(30),
     1  temp(30),
     1  rh(30),
     1  press(30),
     1  cbt(30)

        character*200 directory
        integer len

        call get_directory('radioprd',directory,len)
c-- init


        istat = 0    ! initial assignment of bad istatus

        do k = 1,4
        do i = 1,10
        data(i,k) = 999.
        enddo
        enddo

        npts = 0
        i = 0
        j = 0
        k = 0
        l = 0
        ii = 0
        jj = 0
        kk = 0
        ll = 0

c -- code


        i4time_d = i4time + 600  ! i4time desired is 10 minutes after the
c                                  start of the satellite scan (i4time).

        call cv_i4tim_int_lp
     1          (i4time_d,ny_d,nmo_d,nd_d,nh_d,nm_d,ns_d,istatus)

        nm_d = nm_d/2   !note integer arithmatic  -- returns even number
        nm_d = nm_d*2

c nm_d now contains the minute of the data desired.

c compute the i4time of the profiler data file to access

        i4time_f = (i4time_d/3600)  * 3600      ! note integer arithmatic
c       i4time_f now contains the "nearest hour"
        hour_fraction = float(i4time_d -i4time_f)/3600.
        twenty_seg = int(hour_fraction*3.)+1  !twenty_seg is the increment
c                                            needed to add to i4time_f
        if(nm_d .eq. 00 .or. nm_d.eq.20 .or. nm_d.eq.40) twenty_seg =
     1  twenty_seg-1
        i4time_f = i4time_f + 1200*twenty_seg  ! 1200 is 20 minutes in seconds.
c end computing the desired file time, now open the appropriate file

        call make_fnam_lp (i4time_f,filename,istatus)
        if( istatus.ne.1) return


c get denver data


        iunit = 24

        open (iunit,file=directory(1:len)//filename//'.dmt',
     1  status='old',
     1  err=72)

        call readhd (iunit, header, type, site, points,
     +     elevation, longitude, latitude, ierror)



        call readmet (iunit, header, points, year, month, day,
     + hour, minute, second, tav, pm, mrelev, mrazim, vapor,
     + liquid, temp, rh, press, cbt, md, fh, qc, ierror)

        close (iunit)

        npts = 0
        if(ierror.eq.0) then

        npts = 1+npts
        plat(npts) = latitude
        plon(npts) = -longitude
        do i = 1,points
        if(nm_d.eq.minute(i)) pw(npts) = vapor(i)
        enddo

        endif

        if ( npts.gt.0 ) then
                if (pw(npts).le.0.0) npts = npts-1
        endif

c end denver

72      continue



c get platteville data


        iunit = 24

        open (iunit,file=directory(1:len)//filename//'.pmt',
     1  status='old',
     1  err=73)

        call readhd (iunit, header, type, site, points,
     +     elevation, longitude, latitude, ierror)

        call readmet (iunit, header, points, year, month, day,
     + hour, minute, second, tav, pm, mrelev, mrazim, vapor,
     + liquid, temp, rh, press, cbt, md, fh, qc, ierror)

        close (iunit)


        if(ierror.eq.0) then

        npts = 1+npts
        plat(npts) = latitude
        plon(npts) = -longitude
        do i = 1,points
        if(nm_d.eq.minute(i)) pw(npts) = vapor(i)
        enddo

        endif

        if (npts.gt.0) then
                if (pw(npts).le.0.0) npts = npts-1
        endif

c end platteville


        if(npts.gt.0) istat = 1

73      return
        end





      subroutine readceil (iunit, header, points, year, month,
     + day, hour, minute, second, tav, detect, cbh, cth, ierror)

      integer i,
     +        iunit,
     +        ierror,
     +        points,
     +        year   (*),
     +        month  (*),
     +        day    (*),
     +        hour   (*),
     +        minute (*),
     +        second (*),
     +        tav    (*),
     +        detect (*),
     +        cbh    (2, *),
     +        cth    (2, *)

c      real    elevation,
c     +        longitude,
c     +        latitude

c      character*6 type, site

      character*(*) header

ccc   read in the header lines.
      read (iunit, 100, err=999, end=99) header

ccc   initialize the data point counter.
      i = 1

ccc   read in the data lines.
10    read (iunit, 200, err=999, end=99)
     +    year(i), month(i), day(i),
     +    hour(i), minute(i), second(i),
     +    tav(i),  detect(i), (cbh(j,i), cth(j,i), j=1,2)
c      write (*, 200)
c     +    year(i), month(i), day(i),
c     +    hour(i), minute(i), second(i),
c     +    tav(i),  detect(i), (cbh(j,i), cth(j,i), j=1,2)

ccc   increment the data point counter.
      i = i + 1
      goto 10

ccc   normal, successful read.
99    points = i - 1
      return

ccc   error reading data file.
999   points = i - 1
      ierror = -1

      return

100   format (a, //)
200   format (i4, 2i2, 1x, 3i2, 1x, i5, 1x, i4,
     +  3x, i6, 2x, i6, 3x, i6, 2x, i6)

      end



      subroutine readhd (iunit, header, type, site, points,
     +     elevation, longitude, latitude, ierror)

      integer iunit,
     +        points,
     +        ierror

      real elevation,
     +     longitude,
     +     latitude

      character*(*) header

      character*6 type, site


c     ... read header information ...
      if (iunit .lt. 0) then
         read (header, 100, err=99, end=99) type, site, points,
     +           elevation, longitude, latitude
      else
         read (iunit, 100, err=99, end=99) type, site, points,
     +           elevation, longitude, latitude
      endif


c     ... normal, successful read ...
999   ierror = 0
      return


c     ... error in read process ...
99    ierror = -1
      return


100   format (a6, 1x, a6, 1x, i6, 1x, f8.3, 1x, f8.3, 1x, f8.3)


      end


      subroutine readmet (iunit, header, points, year, month, day,
     + hour, minute, second, tav, pm, mrelev, mrazim, vapor,
     + liquid, temp, rh, press, cbt, md, fh, qc, ierror)

      integer ymd,hms,
     +        i,
     +        iunit,
     +        ierror,
     +        points,
     +        year   (*),
     +        month  (*),
     +        day    (*),
     +        hour   (*),
     +        minute (*),
     +        second (*),
     +        tav    (*),
     +        pm     (*),
     +        md     (*),
     +        fh     (*),
     +        qc     (*)

c      real    elevation,
c     +        longitude,
c     +        latitude,
      real    mrelev (*),
     +        mrazim (*),
     +        vapor  (*),
     +        liquid (*),
     +        temp   (*),
     +        rh     (*),
     +        press  (*),
     +        cbt    (*)

c      character*6 type, site

      character*(*) header

ccc   read in the header lines.
      read (iunit, 100, err=999, end=99) header

ccc   initialize the data point counter.
      i = 1

ccc   read in the data lines.
c10    read (iunit, 200, err=999, end=99)
10    read (iunit, *, err=999, end=99)
     +    ymd,
     +    hms,
     +    tav(i),  pm(i), mrelev(i), mrazim(i),
     +    vapor(i), liquid(i), temp (i), rh (i),
     +    press (i), cbt (i), md(i), fh(i), qc(i)
c      write (*, 200)
c     +    year(i), month(i), day(i),
c     +    hour(i), minute(i), second(i),
c     +    tav(i),  pm(i), mrelev(i), mrazim(i),
c     +    vapor(i), liquid(i), temp (i), rh (i),
c     +    press (i), cbt (i), md(i), fh(i), qc(i)

c break out year month day

      year(i) = nint(float(ymd)/10000.)
      month(i) = ymd-year(i)*10000
      month(i) = nint(float(month(i))/100.)
      day(i) = ymd-year(i)*10000 - month(i)*100


c break out hour min sec

      hour(i) = nint(float(hms)/10000.)
      minute(i) = hms - hour(i)*10000
      minute(i) = nint(float(minute(i))/100.)
      second(i) = hms - hour(i)*10000 - minute(i)*100

ccc   increment the data point counter.
      i = i + 1
      goto 10

ccc   normal, successful read.
99    points = i - 1
      return

ccc   error reading data file.
999   points = i - 1
      ierror = -1

      return

100   format (a, //)
200   format (i4, 2i2, 1x, 3i2, 1x, i5, 1x, i4, 1x, f6.2, 1x, f6.2,
     +  1x, f8.3, 1x, f8.3, 1x, f6.2, 1x, f7.2, 1x, f9.2, 1x, f7.2,
     +  3x, i1, 3x, i1, 2x, i6)

        end
