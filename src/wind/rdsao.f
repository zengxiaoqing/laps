cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
        subroutine rdsao(i4time,heights_3d
     1  ,N_SAO,maxstns
     1  ,lat,lon,sao_i,sao_j,sao_k,sao_u,sao_v,n_sao_obs
     1  ,grid_laps_wt,grid_laps_u,grid_laps_v
     1          ,ni,nj,nk,istatus)

!       1997 Jun    Ken Dritz     Added N_SAO and maxstns as dummy arguments,
!                                 making arrays dimensioned thereby automatic.
!       1997 Jun    Ken Dritz     Removed include of 'lapsparms.for'.

!****************************************************************************

!       LAPS Grid Dimensions

        include 'windparms.inc' ! weight_sao

        real*4 lat(ni,nj)
        real*4 lon(ni,nj)

!       SAO

        integer sao_i(N_SAO) ! X sao coordinates
        integer sao_j(N_SAO) ! Y sao coordinates
        integer sao_k(N_SAO) ! Z sao coordinates
        real    sao_u(N_SAO) ! u sao component
        real    sao_v(N_SAO) ! v sao component

!       Laps Analysis Grids
        real grid_laps_wt(ni,nj,nk)
        real grid_laps_u(ni,nj,nk)
        real grid_laps_v(ni,nj,nk)

!***************************************************************************

        real*4 heights_3d(ni,nj,nk)

        character*13 filename13,c13_fname

        character asc_tim_9*9

        real*4 lat_s(maxstns), lon_s(maxstns), elev_s(maxstns)
        real*4 cover_s(maxstns), hgt_ceil(maxstns), hgt_low(maxstns)
        real*4 t_s(maxstns), td_s(maxstns), pr_s(maxstns), sr_s(maxstns)
        real*4 dd_s(maxstns), ff_s(maxstns), ddg_s(maxstns)  
        real*4 ffg_s(maxstns), vis_s(maxstns)
c
        character stations(maxstns)*3, wx_s(maxstns)*8        ! c5_stamus
        character atime*24, infile*270
        character directory*250,ext*31

!       Declarations for new read_surface routine
!       New arrays.f reading in the SAO data from the LSO files
        real*4   pstn(maxstns),pmsl(maxstns),alt(maxstns)
     1          ,store_hgt(maxstns,5)
        real*4   ceil(maxstns),lowcld(maxstns),cover_a(maxstns)
     1          ,vis(maxstns),rad(maxstns)

        Integer*4   obstime(maxstns),kloud(maxstns),idp3(maxstns)

        Character   obstype(maxstns)*8
     1             ,store_emv(maxstns,5)*1,store_amt(maxstns,5)*4


        n_sao_obs = 0

        ext = 'sag'
        call open_lapsprd_file(32,i4time,ext,istatus)
        if(istatus .ne. 1)go to 888

        call make_fnam_lp(i4time,asc_tim_9,istatus)

        write(6,*)' Reading SAOs: Calling read_sfc ',asc_tim_9
        ext = 'lso'
        call get_directory(ext,directory,len_dir) ! Returns top level directory

        c13_fname = filename13(i4time,ext(1:3))
        infile = directory(1:len_dir)//c13_fname

!       call read_sfc(infile,maxstns,atime,num_meso,num_saos,
!     &        num_sfc,stations,lat_s,lon_s,elev_s,wx_s,cover_s,hgt_ceil,
!     &        hgt_low,t_s,td_s,dd_s,ff_s,ddg_s,ffg_s,pr_s,sr_s,vis_s,istatus)

        call read_surface_old(infile,maxstns,atime,n_meso_g,n_meso_pos,
     1           n_sao_g,n_sao_pos_g,n_sao_b,n_sao_pos_b,n_obs_g,
     1           n_obs_pos_g,
     1           n_obs_b,n_obs_pos_b,stations,obstype,lat_s,lon_s,
     1           elev_s,wx_s,t_s,td_s,dd_s,ff_s,ddg_s,
     1           ffg_s,pstn,pmsl,alt,kloud,ceil,lowcld,cover_a,rad,idp3,
     1           store_emv,store_amt,store_hgt,vis,obstime,istatus)


100     write(6,*)'num_saos',n_sao_b

        if(n_obs_b .gt. maxstns)then
            write(6,*)' Too many stations',n_obs_b,maxstns
            istatus = 0
            return
        endif

        if(n_sao_b .gt. N_SAO)then
            write(6,*)' Too many SAO stations',n_sao_b,N_SAO
            istatus = 0
            return
        endif

        if(istatus .ne. 1)then
            write(6,*)' Warning: Bad istatus from read_sfc',istatus       
            istatus = 1
            return
        endif

        write(6,12)
12      format(/'             Mapping Sao Obs'
     1      /'   n Sta   i   j  k     u      v'
     1      ,'       dd     ff      azi    ran ')


        do i = 1,n_obs_b 

          if(dd_s(i) .ge. 0.0)then

            n_sao_obs = n_sao_obs + 1

            ff_s(i) = ff_s(i) * .518

            call latlon_to_rlapsgrid(lat_s(i),lon_s(i),lat,lon,ni,nj
     1          ,ri,rj,istatus)

            sao_i(n_sao_obs) = nint(ri)
            sao_j(n_sao_obs) = nint(rj)

            call disp_to_uv(dd_s(i),ff_s(i)
     1                     ,sao_u(n_sao_obs),sao_v(n_sao_obs))

!           write(6,*)lat_s(i),lon_s(i),dd_s(i),ff_s(i),elev_s(i)

! ***       Remap SAO observation to LAPS observation grid
!           In bounds?
            if(  sao_i(n_sao_obs) .ge. 1 .and. sao_i(n_sao_obs) .le. ni
     1     .and. sao_j(n_sao_obs) .ge. 1 .and. sao_j(n_sao_obs) .le. nj 
     1                                                             )then       

                rk = height_to_zcoord2(elev_s(i)
     1              ,heights_3d,ni,nj,nk
     1              ,sao_i(n_sao_obs),sao_j(n_sao_obs),istatus)
                if(istatus .ne. 1)then
                    write(6,*)
     1              ' Error: sfc ob may be outside range of ht field'      
                    return
                endif

                sao_k(n_sao_obs) = nint(rk)

!               Fix for stations that are to low for the grid
                if(sao_k(n_sao_obs) .eq. 0)sao_k(n_sao_obs) = 1

                write(32,*)ri-1.,rj-1.,rk-1.,dd_s(i),ff_s(i)

                k = sao_k(n_sao_obs)

                grid_laps_u(sao_i(n_sao_obs),sao_j(n_sao_obs),k)
     1                          = sao_u(n_sao_obs)

                grid_laps_v(sao_i(n_sao_obs),sao_j(n_sao_obs),k)
     1                          = sao_v(n_sao_obs)

                grid_laps_wt(sao_i(n_sao_obs),sao_j(n_sao_obs),k)
     1                          = weight_sao

                write(6,20)n_sao_obs,stations(i)(1:3),
     1                     sao_i(n_sao_obs),
     1                     sao_j(n_sao_obs),
     1                     sao_k(n_sao_obs),
     1                     sao_u(n_sao_obs),
     1                     sao_v(n_sao_obs),
     1                     dd_s(i),ff_s(i)

            else
                write(6,20)n_sao_obs,stations(i)(1:3)

            endif ! In horizontal bounds

20          format(i4,1x,a3,2i4,i3,2f7.1,2x,2f7.1,2x,2f7.1,2x,2f7.1)

          endif ! wind is reported

        enddo

        close(32)

        istatus =1

        return

888     write(6,*)' Open error for SAG file'
        istatus = 0
        close(32)
        return

        end
