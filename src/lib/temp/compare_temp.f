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



        subroutine compare_temp (
     1                tpass1,cgrid,
     1                ni,nj,nk,r_missing_data,                            ! I
     1                obs_barnes_in,max_obs,ncnt_total_in,                ! I
     1                l_withheld_only,                                    ! I
     1                weight_sfc,                                         ! I  
     1                istatus)                                            ! I/O

C****************************************************************************
C
C  Purpose: Provide a single point out of lapstemp_anal to call
C           diagnostic comparision routines.
C
C
C  Inputs: tpass1
C          istat_radar_vel
C          grid_ra_vel
C          rlat_radar
C          rlon_radar
C          rheight_radar
C          n_radars
C
C  outputs: None
C
C*********************************************************************

C***************** Declarations **************************************
        include 'barnesob.inc'
        type (barnesob) :: obs_barnes_in(max_obs)      
        type (barnesob) :: obs_barnes(max_obs)      

        integer istat_radar_vel
        integer l,n_radars,ni,nj,nk,max_radars

        real lat(ni,nj),lon(ni,nj)

        real tpass1(ni,nj,nk)
        real r_missing_data
        real weight_sfc

        integer max_obstypes
        parameter (max_obstypes=11)

        character*4 cgrid
        character*12 c_obstype_a(max_obstypes)
        logical l_parse, l_point_struct, l_withheld_only, l_compare_ob

C********************************************************************

        write(6,*)' Subroutine compare_temp...',cgrid

        l_point_struct = .true.

!       Copy obs structure into local structure depending on 'l_withheld_only'
        ncnt_total = 0
        do i = 1,ncnt_total_in
            if(l_withheld_only)then
                if(obs_barnes_in(i)%l_withhold)then
                    l_compare_ob = .true.
                    write(6,*)i,obs_barnes_in(i)%i,
     1                        obs_barnes_in(i)%j,obs_barnes_in(i)%type
                else
                    l_compare_ob = .false.
                endif
            else
                l_compare_ob = .true.
            endif

            if(l_compare_ob)then     
                ncnt_total = ncnt_total + 1
                obs_barnes(ncnt_total) = obs_barnes_in(i)
            endif
        enddo ! i      

        write(6,*)'l_withheld_only/ncnt_total_in/ncnt_total=',
     1             l_withheld_only,ncnt_total_in,ncnt_total

        if(ncnt_total .eq. 0)then
            write(6,*)' No obs detected, returning...'
            return
        endif

        call get_temp_obstypes (obs_barnes,max_obs,ncnt_total        ! I
     1                         ,c_obstype_a,max_obstypes,n_obstypes  ! I/O
     1                         ,istatus)                             ! O
        if(istatus .ne. 1)stop

!       n_obstypes = 4
!       c_obstype_a(1) = 'SFC '
!       c_obstype_a(2) = 'PROF'
!       c_obstype_a(3) = 'PIN '
!       c_obstype_a(4) = 'CDW '

        do i_obstype = 1,n_obstypes
            call s_len(c_obstype_a(i_obstype),len_obstype)
            write(6,*)

            if(l_withheld_only)then
                write(6,11)cgrid,c_obstype_a(i_obstype)(1:len_obstype)       
 11             format(1x,'  Comparing ',a,' to ',a4
     1                   ,' Obs (Withheld - prior to QC)')    
            elseif(l_parse(cgrid,'FG'))then
                write(6,12)cgrid,c_obstype_a(i_obstype)(1:len_obstype)       
 12             format(1x,'  Comparing ',a,' to ',a4
     1                   ,' Obs (prior to QC)')    
            else
                write(6,13)cgrid,c_obstype_a(i_obstype)(1:len_obstype)
 13             format(1x,'  Comparing ',a,' to ',a4
     1                   ,' Obs (passing QC)')    
            endif

            call comp_grid_tempobs(tpass1,ni,nj,nk
     1          ,weight_sfc
     1          ,obs_barnes,max_obs,ncnt_total,l_point_struct
     1          ,cgrid,c_obstype_a(i_obstype),r_missing_data,rms)

        enddo

        return
        end



        subroutine comp_grid_tempobs(t_3d,ni,nj,nk
     1  ,weight_ob
     1  ,obs_barnes,max_obs,ncnt_total,l_point_struct
     1  ,c_grid,c_obs,r_missing_data,rms)

        include 'barnesob.inc'
        type (barnesob) obs_barnes(max_obs)      

        real t_3d(ni,nj,nk)

        character*4  c_grid
        character*12 c_ob_type,c_obs,c_obs_left,c_obs_right

        logical l_point_struct

        nobs = 0
        residualu = 0.
        residualv = 0.
        sumu = 0.
        sumv = 0.
        sumsp = 0.

        c_obs_left = c_obs
        call left_justify(c_obs_left)
        call s_len(c_obs_left,len_obstype)

        c_obs_right = c_obs
        call right_justify(c_obs_right)

        write(6,*)'Comparing ',c_obs_left(1:len_obstype)
     1           ,' Temp Obs (passing QC) to ',c_grid,' Grid'
        write(6,2)c_obs_right(4:12),c_grid
2       format(1x,'   i   j   k ',a,' Ob (Struct)    Ob (array)'
     1                           ,a,' Analysis        diff')       

        if(l_point_struct)then

!          Checking the obstype is case insensitive
           call downcase(c_obs_left,c_obs_left)

           do iob = 1,ncnt_total
              c_ob_type = obs_barnes(iob)%type
              call downcase(c_ob_type,c_ob_type)
              call left_justify(c_ob_type)

              if(c_ob_type .eq. c_obs_left)then
                  nobs = nobs + 1
                  il = obs_barnes(iob)%i
                  jl = obs_barnes(iob)%j
                  k = obs_barnes(iob)%k
                  difft = obs_barnes(iob)%value(1) - t_3d(il,jl,k) 

                  sumu = sumu + difft

                  residualu = residualu + difft ** 2

                  if(nobs .le. 200 .OR. nobs .eq. (nobs/10)*10)then
                      write(6,101)il,jl,k
     1                ,obs_barnes(iob)%value(1)
     1                ,t_3d(il,jl,k)
     1                ,difft
101                   format(1x,3i4,3(2x,2f7.1))
                  endif

              endif ! obstype match

           enddo ! iob

        endif ! l_point_struct

        if(nobs .gt. 0)then
            rmsu = sqrt(residualu/nobs)
            bias_t = sumu / float(nobs)

        else
            rmsu = 0.
            bias_t = 0.

        endif

        rms  = sqrt(rmsu**2)

        if(nobs .gt. 0)then
            call upcase(c_obs_left,c_obs_left)
            write(6,102)c_obs_left(1:8),c_grid,nobs,bias_t,rmsu
102         format(' BIAS/RMS between '
     1        ,a,' & ',a,' (n,bias_t,rms) = '
     1        ,i6,6f5.1)
        endif

        return

        end


        subroutine get_temp_obstypes(obs_barnes,max_obs,ncnt_total
     1                              ,c_obstype_a,max_obstypes,n_obstypes
     1                              ,istatus)

        include 'barnesob.inc'
        type (barnesob) obs_barnes(max_obs)      

        character*12 c_obstype_a(max_obstypes)

        logical l_match_found

        n_obstypes = 0

        if(ncnt_total .eq. 0)then
            return
        endif        

        write(6,*)' Subroutine get_temp_obstypes, obstypes found...'

        i = 1
        n_obstypes = 1
        write(6,*)n_obstypes,i,obs_barnes(i)%type
        c_obstype_a(i) = obs_barnes(i)%type

        if(ncnt_total .ge. 2)then
            do i = 1,ncnt_total
                l_match_found = .false.
                do j = 1,n_obstypes
                    if(obs_barnes(i)%type .eq. c_obstype_a(j))then
                        l_match_found = .true.
                    endif
                enddo ! j
                if(.not. l_match_found)then
                    n_obstypes = n_obstypes + 1
                    write(6,*)n_obstypes,i,obs_barnes(i)%type
                    if(n_obstypes .gt. max_obstypes)then
                        write(6,*)' ERROR: too many obstypes'
                        istatus = 0
                        return
                    endif
                    c_obstype_a(n_obstypes) = obs_barnes(i)%type
                endif
            enddo ! i
        endif

        istatus = 1

        return
        end 
