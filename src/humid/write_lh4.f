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
        subroutine write_lh4 (i4time,tpw,bias_one,ii,jj,istatus)

c       $log: write_lh4.for,v $
c revision 1.4  1996/03/05  21:11:53  birk
c modified i/o interface to enable subhourly cycle
c
c revision 1.3  1995/09/13  21:36:16  birk
c added disclaimer to files
c
c revision 1.2  1994/10/24  17:29:36  birk
c installed variable for missing data flag
c
c revision 1.1  1994/04/25  15:09:27  birk
c initial revision
c

c       modified 4/19/94 moved inc file to after data for cray port db

        implicit none

c  parameter variables

      integer ii,jj
      real tpw (ii,jj)      
      real bias_one
      integer istatus

c internal variables

        integer
     1  i4time,
     1  kmax,
     1  lvl(1)
      real mdf

        integer i,j, len

        character
     1  tpdir*50,
     1  tpext*31,
     1  var(1)*3,
     1  lvl_coord(1)*4,
     1  units(1)*10,
     1  comment(1)*125

        data var/1*'tpw'/
        data lvl_coord/1*'  '/
        data units/1*'m  '/
        data tpext /'lh4'/
c  set missing data flag
      call get_r_missing_data(mdf, istatus)
      if (istatus.ne.1) then
               write(6,*) 'Fatal error in assigning missing data flag'
               write(6,*) 'Assigning default value of 1e+37'
               mdf = 1.e37
               write(6,*) 'Continuing using default'
      endif

      call get_directory(tpext,tpdir,len)

        istatus = 0 ! bad

              if(tpw(1,1).eq.mdf) then ! field is not valid
                     write(6,*) 'TPW field not valid'
                     return
              endif




        write(comment(1),1)  bias_one
1       format('vsn 5.1 tpw via sh integration, bias correction = ',e20.
     110)

        kmax = 1

c       convert from cm to meters for archive

        do i  = 1,ii
        do j  = 1,jj

        tpw(i,j) = tpw(i,j) * 1.e-2
c       place mod for missing data flag
        if (tpw(i,j) .lt. 0.0) tpw (i,j) = mdf

        enddo
        enddo


        call write_laps (i4time,i4time,
     1  tpdir,
     1  tpext,
     1  ii,
     1  jj,
     1  1,
     1  1,
     1  var,
     1  lvl,
     1  lvl_coord,
     1  units,
     1  comment,
     1  tpw,
     1  istatus)

        if(istatus.ne.1) then
        istatus = 134316524
        return
        endif

        istatus = 1 ! success

        print*, 'write laps tpw complete',istatus

        return
        end
