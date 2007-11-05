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
        subroutine      get_maps_df (i4time,df,htm,ii,jj,kk,lct,istatus)

c       $log: getmapsdf.for,v $
c revision 1.7  1996/07/22  16:43:35  birk
c extended look back for background field to 30 hours
c
c revision 1.6  1996/03/05  21:09:00  birk
c modified i/o interface to enable subhourly cycle
c
c revision 1.5  1995/09/13  21:35:21  birk
c added disclaimer to files
c
c revision 1.4  1995/09/06  14:25:30  birk
c enabled routine to work with hourly lga files.
c note: the routine used to interpolate, now it does
c not.  furthermore, the routine no longer will look for
c lgf files.  these are no longer generated.
c
c revision 1.3  1994/09/09  16:09:17  birk
c fixed mixed mode problem in pressure_of_level
c
c revision 1.2  1994/07/22  22:01:36  birk
c made code compatible with lapsparms.inc
c
c revision 1.1  1994/04/25  15:08:30  birk
c initial revision
c
c
c ****  this routine reads maps analysis (lma) and forecast (lmf) heights
c          and returns the 8x8 field for the requested time...no interpolation.
c          the routine uses the latest maps run in which analysis/forecast
c          valid times are available for before and after the requested
c          i4time.

c       author j. snook (initial version)
c       adapted by d. birkenheuer  21 feb 91
c       included unix constructs  19 oct 92  db.
c       revamped code to remove the call to get_file_names (highly vms
c       dependent)   2/24/93 db
c
c-------non-executable statements-----------------------------------------------
c
        implicit        none
c
        
        CHARACTER*250 DIRLMA
        CHARACTER*250 DIRLMF
        CHARACTER*31 EXTLMA
        CHARACTER*31 EXTLMF
        integer len
        DATA EXTLMA /'lga'/
        DATA EXTLMF /'lgf'/        

c input parameter

      integer i4time,istatus
      character*4     df               ! desired field (e.g., 'rh ')
      integer ii,jj,kk
      real          htm(ii,jj,kk)
      integer lct !laps cycle time

c  internal variables that depend on dynamic allocation

      integer lapsp (kk)
      character*3     var(kk)
      character*4     lvl_coord(kk)
      character*10    units(kk)
      character*125   comment(kk)

c  internal variables

        integer laps_cycle_time

        integer 
     1          i4time1,
     1          i,k

        integer iter ! iterations reqd to go back 48 hour

        real pressure_of_level ! function call
c
c
c
c
c------------------------------------------------------------------------------
c
        call get_directory(extlma,dirlma,len)
        call get_directory(extlmf,dirlmf,len)
        istatus = 0
        laps_cycle_time = lct

        do k = 1,kk
           lapsp(k) =  nint( pressure_of_level(k)  * .01 )
           var(k) = df
        enddo


        print *, ' '
        print *, ' '
        print *, ' '
        print *, ' '
        print *, ' '
        print *, ' '

c
c ****  new nomenclature, i4time=validtime of data requested (could be
c       forecast time)
c
c       i4time1 = reftime, the 00h forecast time or analysis time.
c

c       cycle through all possible forecasts for current time.

c     add direct call to get laps cycle time so that changes to this can be
c     managed at runtime.
        call get_laps_cycle_time(laps_cycle_time,istatus)
        if (istatus .ne. 1) then
           write(6,*) 'Error obtaining laps_cycle time'
           write (6,*) 'Aborting'
           return
        endif

c     figure iterations to go back 48 hours

        write(6,*) 'laps cycle time is, ', laps_cycle_time
    
        iter = 72

        if(laps_cycle_time.lt.3600) then
           write(6,*) 'will look back ',iter*laps_cycle_time,' seconds'
        else
           write(6,*) 'will look back ',iter*3600,' seconds'
        endif

        write(6,*) 'iterations to attempt = ',iter

        do i = 0,iter

           i4time1 = i4time - i*laps_cycle_time
           if (laps_cycle_time.gt.3600) i4time1 = i4time - i*3600
           


c     attempt reading data for this time pair

           print *, 'attempting to open from lga ', i4time-i4time1

           call read_laps(i4time1,i4time, dirlma,
     1          extlma,ii,jj,
     1          kk,kk,var,lapsp,
     1          lvl_coord,units,comment,htm,istatus)

           if (istatus.eq.1) then !successful read
              write(6,*)  'successful read of input file'
              istatus = 1
              return

           endif

        enddo                   ! i





        write(6,*) 'error finding suitable maps for past'
        write(6,*) ' failure on interpolating to get a good background'
        istatus = 0
        return


c     
        end
