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
        subroutine lh3_compress (data,tdata,i4time,lvl,
     1      ii,jj,kk,istatus)

c       $log: lh3_compress.for,v $
c revision 1.5  1996/03/05  21:11:36  birk
c modified i/o interface to enable subhourly cycle
c
c revision 1.4  1995/09/13  21:35:52  birk
c added disclaimer to files
c
c revision 1.3  1995/02/14  15:30:31  birk
c modified the rh variable to rh3 to be compliant with the overall plan give all
c fields unique variable names.
c
c revision 1.2  1994/10/24  17:28:22  birk
c installed variable for missing data flag
c
c revision 1.1  1994/04/25  15:06:16  birk
c initial revision
c

c       birkenheuer
c
c       5/15/91 put in test for -1e30 (removed counter and flags for less than
c                       zero condition
c     1 /4/93  revised to put into normal write laps output
c       note ::  the parameter list changed this iteration including the new
c       variable lvl so the levels will not have to be recomputed.

c       5/14/93: revised to include the improved call to moisture conversion
c       routines to differentiate between ice and liquid phase in the analysis.
c
c     1 0/28/93 put in the additions to enable rh wrt liq.  db


        implicit none

c        include 'lapsparms.for'
c        include 'parmtrs.inc'

c parameter variables

      integer ii,jj,kk

      real data (ii,jj,kk)
      real tdata (ii,jj,kk)
      integer i4time
      integer lvl (kk)
      integer istatus

c variables requiring dynamic allocation

      integer double_lvl(kk*2)
        real rhdata(ii,jj,kk),rhdata_l(ii,jj,kk),
     1  equivalenced_rh(ii,jj,kk*2)
      character var(kk*2)*3,
     1        lvl_coord(kk*2)*4,
     1        units(kk*2)*10,
     1        comment(kk*2)*125
c        data var/kk*'rh3',kk*'rhl'/
c        data lvl_coord/kk*'hpa',kk*'hpa'/
c        data units /kk*'percent',kk*'percent'/

c internal variables

      character extlt1*31,ext*50,rhext*50,extpw*50,ext3*50
        character
     1  dirlt1*50,dir*50,rhdir*50,dirpw*50,dir3*50
      real make_rh   !function
      external make_rh

      integer i,j,k,kbottom, len
      real rmd
      data extpw/'lh1'/
      data ext3/'lh2'/
      data extlt1/'lt1'/
      data ext /'lq3'/
      data rhext /'lh3'/

      istatus = 0               ! bad

      call get_directory(extpw,dirpw,len)
      call get_directory(ext3,dir3,len)
      call get_directory(extlt1,dirlt1,len)
      call get_directory(ext,dir,len)
      call get_directory(rhext,rhdir,len)

c   initialize

c     note that dynamic assignments don't work in data statements

	do k = 1,kk
	   var(k) = 'rh3'
           var(k+kk) = 'rhl'
	   lvl_coord(k) = 'hpa'
	   lvl_coord(k+kk) = 'hpa'
	   units(k) = 'percent'
	   units(k+kk) = 'percent'
	 enddo
											    !k


      call get_r_missing_data(rmd, istatus)

        do k = 1,kk
        double_lvl(k) = lvl(k)
        enddo
        do k = kk+1,kk*2
        double_lvl(k) = lvl(k-kk)
        enddo


c....loop for rh computation

        do k = 1,kk
        do j=1,jj
        do i=1,ii


        if (tdata(i,j,k) .lt. 200. .or. data(i,j,k) .eq. rmd) then
        rhdata(i,j,k) = 0   ! pseudo bad data flag for albers routine
        rhdata_l(i,j,k) = 0   ! pseudo bad data flag for albers routine
        else

c       new addition to call make_rh() note: q must be entered as g/kg
c       as is the convention for the calls thusfar for q.
c       a t_ref of 0.0c is currently default in all routines for the tran-
c       sition temperature for ice and liquid vapor reference.

        rhdata(i,j,k) = 100.* make_rh ( float(lvl(k)),tdata(i,j,k)-273.
     115,
     1  data(i,j,k)*1000., 0.0)
        rhdata_l(i,j,k) = 100.* make_rh ( float(lvl(k)),
     1  tdata(i,j,k)-273.15,
     1  data(i,j,k)*1000., -50.)

        rhdata(i,j,k) = max(0.,rhdata(i,j,k) )
        rhdata(i,j,k) = min(100.,rhdata(i,j,k) )
        rhdata_l(i,j,k) = max(0.,rhdata_l(i,j,k) )
        rhdata_l(i,j,k) = min(100.,rhdata_l(i,j,k) )

        endif

        enddo
        enddo
        enddo

c       put in replication for albers code

c  first pass for regular rh fields

        do i = 1,ii
        do j = 1,jj

        kbottom = kk

        do k = 1,kk
        if(rhdata(i,j,k) .ne. 0.0) then
                kbottom = k
                go to 21
        endif
        enddo
21      continue


        if (kbottom .eq. kk .or. kbottom .eq. 1) then
                continue ! not valid bottom found
        else
                do k = kbottom-1,1, -1
                        rhdata(i,j,k) = rhdata(i,j,kbottom)
                enddo
        endif

        enddo
        enddo


c  second  pass for special rh field

        do i = 1,ii
        do j = 1,jj

        kbottom = kk

        do k = 1,kk
        if(rhdata_l(i,j,k) .ne. 0.0) then
                kbottom = k
                go to 22
        endif
        enddo
22      continue


        if (kbottom .eq. kk .or. kbottom .eq. 1) then
                continue ! not valid bottom found
        else
                do k = kbottom-1,1, -1
                        rhdata_l(i,j,k) = rhdata_l(i,j,kbottom)
                enddo
        endif

        enddo
        enddo

c   replace equivalence statments with do loops

      do k = 1,kk
      do j = 1,jj
      do i = 1,ii

        equivalenced_rh(i,j,k) = rhdata(i,j,k)
        equivalenced_rh(i,j,k+kk) = rhdata_l(i,j,k)

      enddo
      enddo
      enddo

        call write_laps (i4time,i4time,
     1        rhdir,
     1        rhext,
     1        ii,
     1        jj,
     1        kk*2,
     1        kk*2,
     1        var,
     1        double_lvl,
     1        lvl_coord,
     1        units,
     1        comment,
     1        equivalenced_rh,
     1        istatus)

        if(istatus.eq.0) then
        return
        endif



        istatus = 1  ! good

        return

        end
