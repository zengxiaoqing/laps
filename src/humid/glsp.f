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
        subroutine glsp (i4time,p,ii,jj,istatus)

c       $log: glsp.for,v $
c revision 1.3  1996/03/13  22:43:30  birk
c changed read_laps_data to read_laps to allow a sub-hourly cycle
c
c revision 1.2  1995/09/13  21:35:22  birk
c added disclaimer to files
c
c revision 1.1  1994/04/25  15:14:19  birk
c initial revision
c

c       gets laps surface pressure

c       input i4time

        implicit none

c        include 'parmtrs.inc' removed dependency 7/29/97 db

c parameter variables

      integer i4time
      integer ii,jj
      real p (ii,jj)
      integer istatus

c variables with dynamic dependence

      real data(ii,jj,1)


c internal variables

        integer i,j
        character*256 dir
        character*31 ext
        character*3 var (1)
        integer lvl (1)
        character*4 lvl_coord  (1)
        character*10 units (1)
        character*125 comment (1)
        integer kmax, ksurf, len

        real numel,sum
        data ext/'lsx'/
 
        call get_directory (ext,dir,len)

        sum = 0.0
        numel = 0.0

        istatus = 0

        lvl (1) = 0


        var(1) = 'ps'
        kmax = 1
        ksurf = 1

        call read_laps (i4time,i4time,dir,ext,ii,jj,
     1  kmax,ksurf,var,lvl,
     1  lvl_coord,units,comment,data,istatus)

        if (istatus.ne.1) then
        write(6,*) 'error in routine glsp'
        endif

        write (6,*) istatus
        write (6,*) lvl_coord,units,comment


        do j = 1,jj
        do i = 1,ii

        numel = numel + 1.
        sum = sum + data(i,j,1)
        p(i,j) = data(i,j,1)


        enddo
        enddo

        sum = sum /numel
        write (6,*) 'average value of field is ',sum


        istatus = 1


        return

        end
