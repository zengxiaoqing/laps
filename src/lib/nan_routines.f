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


       subroutine check_nan (var,istatus)

cdoc  This routine checks a variable for Not a Number conditions (NaN).  It 
cdoc  returns an istatus of 0 if NaNs are detected. Otherwise it returns a 
cdoc  value of 1.

c  author: Dan Birkenheuer
c  date:   12/11/96
c  date:    1/15/97  modified into 4 routines Brandy can't handle 
c                    entry points.
c  date:   2011      modified by Steve Albers to use division test since
c                    inequality test doesn't work with PGF90 unless special
c                    compiler flags are used

       real var
       integer istatus



c  single variable
        if(var .ne. 0.0)then
          if (var / var .NE. 1.0) then ! NaN detected
            istatus = 0
            return
          endif
        endif


        istatus = 1
        return



       end



       subroutine check_nan1 (var1,n,istatus)

cdoc  This routine checks a 1d array for Not a Number conditions (NaN).  It 
cdoc  returns an istatus of 0 if NaNs are detected. Otherwise it returns a 
cdoc  value of 1.

c  author: Dan Birkenheuer
c  date:   12/11/96
c  date:    1/15/97  modified into 4 routines Brandy can't handle 
c                    entry points.
c  date:   2011      modified by Steve Albers to use division test since
c                    inequality test doesn't work with PGF90 unless special
c                    compiler flags are used
       integer istatus
       real var1(n)
       integer i
       integer n




c  single dimension

        do i = 1,n
          if(var1(i) .NE. 0.) then
            if(var1(i) / var1(i) .NE. 1.0) then ! NaN detected
              istatus = 0
              write(6,*)' Nan detected at i: ',i
              return
            endif
          endif
        enddo


        istatus = 1
        return




        end



       subroutine  check_nan2 (var2,n,m,istatus)

cdoc  This routine checks a 2d array for Not a Number conditions (NaN).  It 
cdoc  returns an istatus of 0 if NaNs are detected. Otherwise it returns a 
cdoc  value of 1.

c  author: Dan Birkenheuer
c  date:   12/11/96
c  date:    1/15/97  modified into 4 routines Brandy can't handle 
c                    entry points.
c  date:   2011      modified by Steve Albers to use division test since
c                    inequality test doesn't work with PGF90 unless special
c                    compiler flags are used
       integer istatus
       real var2(n,m)
       integer i,j
       integer m,n



c  double dimension

        do j = 1,m
        do i = 1,n
          if(var2(i,j) .NE. 0.0) then 
            if(var2(i,j) / var2(i,j) .NE. 1.0) then ! NaN detected
              istatus = 0
              write(6,*)' Nan detected at i,j: ',i,j
              return
            endif
          endif
        enddo
        enddo


        istatus = 1
        return

        end


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


       subroutine check_nan3 (var3,n,m,l,istatus)

cdoc  This routine checks a 3d array for Not a Number conditions (NaN).  It 
cdoc  returns an istatus of 0 if NaNs are detected. Otherwise it returns a 
cdoc  value of 1.

c  author: Dan Birkenheuer
c  date:   12/11/96
c  date:    1/15/97  modified into 4 routines Brandy can't handle 
c                    entry points.
c  date:   2011      modified by Steve Albers to use division test since
c                    inequality test doesn't work with PGF90 unless special
c                    compiler flags are used
       integer istatus
       real var3(n,m,l)
       integer i,j,k
       integer l,m,n




c  triple dimension

        do k = 1,l
        do j = 1,m
        do i = 1,n
          if(var3(i,j,k) .NE. 0.) then ! NaN detected
            if(var3(i,j,k) / var3(i,j,k) .NE. 1.0) then ! NaN detected
c LW added for testing
              write(6,*) 'NaN found: ',i,j,k
              istatus = 0
              return
            endif
          endif
        enddo
        enddo
        enddo


        istatus = 1
        return

        end



      function nanf(arg)

cdoc  Function subroutine to check a variable for Nan. 1 means a Nan was 
cdoc  detected, otherwise it returns a 0.

      call check_nan(arg,istatus)

      nanf = 1 - istatus

      return
      end

      

       subroutine  clean_nan2 (var2,n,m,istatus)

       use mem_namelist, ONLY: r_missing_data

       integer istatus
       real var2(n,m)
       integer i,j
       integer m,n

c  double dimension

        do j = 1,m
        do i = 1,n
          if(var2(i,j) .NE. 0.0) then 
            if(var2(i,j) / var2(i,j) .NE. 1.0) then ! NaN detected
              var2(i,j) = r_missing_data
            endif
          endif
        enddo
        enddo


        istatus = 1
        return

        end
