cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps
cdis
cdis    this software and its documentation are in the public domain and
cdis    are furnished "as is."  the united states government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  they assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  if significant modifications or enhancements
cdis    are made to this software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
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

        integer*4 i,j
        character*256 dir
        character*31 ext
        character*3 var (1)
        integer*4 lvl (1)
        character*4 lvl_coord  (1)
        character*10 units (1)
        character*125 comment (1)
        integer*4 kmax, ksurf, len

        real*4 numel,sum
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
