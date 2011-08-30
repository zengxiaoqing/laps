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
        subroutine rdsfc(i4time,heights_3d                              ! I
     1  ,N_SFC,maxstns                                                  ! I
     1  ,lat,lon,r_missing_data                                         ! I
     1  ,n_sfc_obs                                                      ! O
     1  ,grid_laps_wt,grid_laps_u,grid_laps_v                           ! I/O
     1  ,max_obs,obs_point,nobs_point                                   ! I/O
     1  ,ni,nj,nk                                                       ! I
     1  ,istatus)                                                       ! O

!****************************************************************************

        use mem_namelist, ONLY: iwrite_output

!	2006 Feb   Yuanfu Xie	Change rk of obs_point to real rk from
!				sfc_k function and add ri rj values.

        include 'barnesob.inc'
        type (barnesob) :: obs_point(max_obs)                           

!       LAPS Grid Dimensions

        include 'windparms.inc' ! weight_sfc

        real lat(ni,nj)
        real lon(ni,nj)

!       SFC

        integer sfc_i(N_SFC) ! X Sfc coordinates
        integer sfc_j(N_SFC) ! Y Sfc coordinates
        integer sfc_k(N_SFC) ! Z Sfc coordinates
        real    sfc_u(N_SFC) ! u Sfc component
        real    sfc_v(N_SFC) ! v Sfc component

!       Laps Analysis Grids
        real grid_laps_wt(ni,nj,nk)
        real grid_laps_u(ni,nj,nk)
        real grid_laps_v(ni,nj,nk)

!***************************************************************************

        real heights_3d(ni,nj,nk)

        character*13 filename13,c13_fname

        character asc_tim_9*9

        real lat_s(maxstns), lon_s(maxstns), elev_s(maxstns)
        real cover_s(maxstns), hgt_ceil(maxstns), hgt_low(maxstns)
        real t_s(maxstns), td_s(maxstns), pr_s(maxstns), sr_s(maxstns)
        real dd_s(maxstns), ff_s(maxstns), ddg_s(maxstns)  
        real ffg_s(maxstns), vis_s(maxstns)
        real dd_ea(maxstns), ff_ea(maxstns)
c
!       character stations(maxstns)*3, wx_s(maxstns)*8        ! c5_stamus
        character stations(maxstns)*20, provider(maxstns)*11
        character atime*24, infile*270
        character directory*250,ext*31

!       Declarations for new read_surface routine
!       New arrays.f reading in the SFC data from the LSO files
        real   pstn(maxstns),pmsl(maxstns),alt(maxstns)
     1          ,store_hgt(maxstns,5)
        real   ceil(maxstns),lowcld(maxstns),cover_a(maxstns)
     1          ,vis(maxstns),rad(maxstns)

        Integer   obstime(maxstns),kloud(maxstns),idp3(maxstns)

        Character   obstype(maxstns)*8
     1             ,store_emv(maxstns,5)*1,store_amt(maxstns,5)*4


        n_sfc_obs = 0

        ext = 'sag'
        if(iwrite_output .ge. 1)then
            call open_lapsprd_file(32,i4time,ext,istatus)
            if(istatus .ne. 1)go to 888
        endif

        call make_fnam_lp(i4time,asc_tim_9,istatus)

        write(6,*)' Reading SFC Obs: Calling read_sfc ',asc_tim_9
        write(6,*)' N_SFC, maxstns = ',N_SFC, maxstns

        ext = 'lso'
        call get_directory(ext,directory,len_dir) ! Returns top level directory

        c13_fname = filename13(i4time,ext(1:3))
        infile = directory(1:len_dir)//c13_fname

        n_obs_b = 0

        call read_sfc_wind(i4time,'lso',n_obs_g,n_obs_b,obstime
     1                    ,stations,provider,lat_s,lon_s,elev_s
     1                    ,dd_s,ff_s,dd_ea,ff_ea,maxstns,istatus)

100     write(6,*)'n_obs_b=',n_obs_b

        if(n_obs_b .gt. maxstns)then
            write(6,*)' Too many stations',n_obs_b,maxstns
            istatus = 0
            return
        endif

        if(istatus .ne. 1)then
            write(6,*)' Warning: Bad istatus from read_sfc',istatus       
            istatus = 1
            return
        endif

        write(6,12)
12      format(/'             Mapping Sfc Obs'
     1      /'   n  Sta    i   j  k     u      v'
     1      ,'       dd     ff      azi    ran ')

        n_in_domain = 0
        n_overwrite = 0

        do i = 1,n_obs_b 

          if(dd_s(i) .ge. 0.0 .and.
     1       ff_s(i) .ge. 0.0       )then ! Note badflag = -99.9

            n_sfc_obs = n_sfc_obs + 1

            ff_s(i) = ff_s(i) * .518

            call latlon_to_rlapsgrid(lat_s(i),lon_s(i),lat,lon,ni,nj
     1          ,ri,rj,istatus)

            sfc_i(n_sfc_obs) = nint(ri)
            sfc_j(n_sfc_obs) = nint(rj)

            call disp_to_uv(dd_s(i),ff_s(i)
     1                     ,sfc_u(n_sfc_obs),sfc_v(n_sfc_obs))

!           write(6,*)lat_s(i),lon_s(i),dd_s(i),ff_s(i),elev_s(i)

            call left_justify(stations(i))

! ***       Remap SFC observation to LAPS observation grid
!           In bounds?
            if(  sfc_i(n_sfc_obs) .ge. 1 .and. sfc_i(n_sfc_obs) .le. ni
     1     .and. sfc_j(n_sfc_obs) .ge. 1 .and. sfc_j(n_sfc_obs) .le. nj 
     1                                                             )then       

                rk = height_to_zcoord2(elev_s(i)
     1              ,heights_3d,ni,nj,nk
     1              ,sfc_i(n_sfc_obs),sfc_j(n_sfc_obs),istatus)
                if(istatus .ne. 1)then
                    write(6,*)
     1              ' Error: sfc ob is outside range of ht field'        
                    write(6,*)rk,i,elev_s(i),sfc_i(n_sfc_obs),
     1                                       sfc_j(n_sfc_obs),
     1                        (heights_3d(sfc_i(n_sfc_obs)
     1                                   ,sfc_j(n_sfc_obs),k),k=1,nk)       
                    return
                endif

                sfc_k(n_sfc_obs) = nint(rk)

!               Fix for stations that are to low for the grid
                if(sfc_k(n_sfc_obs) .eq. 0)sfc_k(n_sfc_obs) = 1

                if(iwrite_output .ge. 1)then
                    write(32,*)ri,rj,rk,dd_s(i),ff_s(i)
                endif

                k = sfc_k(n_sfc_obs)

                if(grid_laps_wt(sfc_i(n_sfc_obs),sfc_j(n_sfc_obs),k)     
     1                  .ne. r_missing_data)then
                    if(n_overwrite .le. 100)then
                        write(6,*)
     1              ' This next station overwrites another on the grid'
                    endif
                    n_overwrite = n_overwrite + 1
                endif

                grid_laps_u(sfc_i(n_sfc_obs),sfc_j(n_sfc_obs),k)
     1                          = sfc_u(n_sfc_obs)

                grid_laps_v(sfc_i(n_sfc_obs),sfc_j(n_sfc_obs),k)
     1                          = sfc_v(n_sfc_obs)

                grid_laps_wt(sfc_i(n_sfc_obs),sfc_j(n_sfc_obs),k)
     1                          = weight_sfc

                if(n_in_domain .le. 100)then
                    write(6,20)n_sfc_obs,stations(i)(1:5),
     1                     sfc_i(n_sfc_obs),
     1                     sfc_j(n_sfc_obs),
     1                     sfc_k(n_sfc_obs),
     1                     sfc_u(n_sfc_obs),
     1                     sfc_v(n_sfc_obs),
     1                     dd_s(i),ff_s(i)
                endif

                n_in_domain = n_in_domain + 1

!               Add to data structure
                nobs_point = nobs_point + 1
                obs_point(nobs_point)%i = sfc_i(n_sfc_obs)
                obs_point(nobs_point)%j = sfc_j(n_sfc_obs)
                obs_point(nobs_point)%k = sfc_k(n_sfc_obs)
                obs_point(nobs_point)%ri = ri	! Yuanfu
                obs_point(nobs_point)%rj = rj	! Yuanfu
                obs_point(nobs_point)%rk = rk ! sfc_k(n_sfc_obs) Yuanfu
                obs_point(nobs_point)%valuef(1) = sfc_u(n_sfc_obs)
                obs_point(nobs_point)%valuef(2) = sfc_v(n_sfc_obs)
                obs_point(nobs_point)%weight = weight_sfc
                obs_point(nobs_point)%vert_rad_rat = 0.5
                obs_point(nobs_point)%type   = 'sfc'
                obs_point(nobs_point)%file   = 'lso'

                call get_sfc_obtime(obstime(n_sfc_obs),i4time
     1                             ,i4time_ret,istatus)     
                obs_point(nobs_point)%i4time = i4time_ret

            else
                write(6,20)n_sfc_obs,stations(i)(1:5)

            endif ! In horizontal bounds

20          format(i4,1x,a5,2i4,i3,2f7.1,2x,2f7.1,2x,2f7.1,2x,2f7.1)

          endif ! wind is reported

        enddo

        close(32)

        write(6,*)' # of valid sfc stations in domain = ',n_in_domain
        write(6,*)' # of unique grid points with sfc stations = '
     1           ,n_in_domain - n_overwrite

        istatus =1

        return

888     write(6,*)' Open error for SAG file'
        istatus = 0
        close(32)
        return

        end

