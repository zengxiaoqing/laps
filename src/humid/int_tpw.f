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
        subroutine  int_tpw(data,kstart,qs,ps,plevel,tpw,ii,jj,kk)

c       $log: int_tpw.for,v $
c revision 1.3  1995/09/13  21:35:23  birk
c added disclaimer to files
c
c revision 1.2  1994/07/22  22:02:22  birk
c made code compatible with lapsparms.inc
c
c revision 1.1  1994/04/25  15:01:49  birk
c initial revision
c

c       subroutine to integrate the total precipitable water from the
c       specific humidity data.  consolodating this process, since it is
c       done numerous times in the code

c       birkenheuer  feb 9 1993

        implicit none

        include 'lapsparms.for'
        include 'parmtrs.inc'

c parameter variables

      integer ii,jj,kk
        real data(ii,jj,kk)
        integer kstart (ii,jj)
        real qs(ii,jj)
        real ps(ii,jj)
        real plevel (kk)
        real tpw (ii,jj)


c variables requireing dynamic allocation

c               NONE

c internal variables

        real delta_p
        integer i,j,k

c       set delta_p

        delta_p = plevel(1)-plevel(2)

c       integrate the tpw field

        do j = 1,jj
        do i = 1,ii
c compute the dewpoints
c integrate
        tpw(i,j) = 0.0
        do k = kstart(i,j),kk-1
        tpw(i,j) = tpw(i,j) + ( data(i,j,k)+data(i,j,k+1) )/2. * delta_p
        enddo
        tpw(i,j) = tpw(i,j) *1.e3  !  change to g/kg mb
        tpw(i,j) = tpw(i,j)
     1  + ( qs(i,j) + data(i,j,kstart(i,j)) *1.e3 ) /2.
     1  * (ps(i,j)-plevel(kstart(i,j)) )
        tpw(i,j) = tpw(i,j) / 100. / 9.8

        enddo
        enddo

        return

        end
