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
        subroutine rdsfc(i4time,heights_3d
     1  ,N_SFC,maxstns
     1  ,lat,lon,n_sfc_obs
     1  ,grid_laps_wt,grid_laps_u,grid_laps_v
     1  ,ni,nj,nk,istatus)

!****************************************************************************

!       LAPS Grid Dimensions

        include 'windparms.inc' ! weight_sfc

        real*4 lat(ni,nj)
        real*4 lon(ni,nj)

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
!       New arrays.f reading in the SFC data from the LSO files
        real*4   pstn(maxstns),pmsl(maxstns),alt(maxstns)
     1          ,store_hgt(maxstns,5)
        real*4   ceil(maxstns),lowcld(maxstns),cover_a(maxstns)
     1          ,vis(maxstns),rad(maxstns)

        Integer*4   obstime(maxstns),kloud(maxstns),idp3(maxstns)

        Character   obstype(maxstns)*8
     1             ,store_emv(maxstns,5)*1,store_amt(maxstns,5)*4


        n_sfc_obs = 0

        ext = 'sag'
        call open_lapsprd_file(32,i4time,ext,istatus)
        if(istatus .ne. 1)go to 888

        call make_fnam_lp(i4time,asc_tim_9,istatus)

        write(6,*)' Reading SFC Obs: Calling read_sfc ',asc_tim_9
        write(6,*)' N_SFC, maxstns = ',N_SFC, maxstns

        ext = 'lso'
        call get_directory(ext,directory,len_dir) ! Returns top level directory

        c13_fname = filename13(i4time,ext(1:3))
        infile = directory(1:len_dir)//c13_fname

        call read_surface_old(infile,maxstns,atime,n_meso_g,n_meso_pos,
     1           n_sao_g,n_sao_pos_g,n_sao_b,n_sao_pos_b,n_obs_g,
     1           n_obs_pos_g,
     1           n_obs_b,n_obs_pos_b,stations,obstype,lat_s,lon_s,
     1           elev_s,wx_s,t_s,td_s,dd_s,ff_s,ddg_s,
     1           ffg_s,pstn,pmsl,alt,kloud,ceil,lowcld,cover_a,rad,idp3,
     1           store_emv,store_amt,store_hgt,vis,obstime,istatus)

100     write(6,*)'n_sao_b=',n_sao_b

        if(n_obs_b .gt. maxstns)then
            write(6,*)' Too many stations',n_obs_b,maxstns
            istatus = 0
            return
        endif

        if(n_sao_b .gt. N_SFC)then
            write(6,*)' Too many SFC stations',n_sao_b,N_SFC
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
     1      /'   n Sta   i   j  k     u      v'
     1      ,'       dd     ff      azi    ran ')


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

                write(32,*)ri-1.,rj-1.,rk-1.,dd_s(i),ff_s(i)

                k = sfc_k(n_sfc_obs)

                grid_laps_u(sfc_i(n_sfc_obs),sfc_j(n_sfc_obs),k)
     1                          = sfc_u(n_sfc_obs)

                grid_laps_v(sfc_i(n_sfc_obs),sfc_j(n_sfc_obs),k)
     1                          = sfc_v(n_sfc_obs)

                grid_laps_wt(sfc_i(n_sfc_obs),sfc_j(n_sfc_obs),k)
     1                          = weight_sfc

                write(6,20)n_sfc_obs,stations(i)(1:3),
     1                     sfc_i(n_sfc_obs),
     1                     sfc_j(n_sfc_obs),
     1                     sfc_k(n_sfc_obs),
     1                     sfc_u(n_sfc_obs),
     1                     sfc_v(n_sfc_obs),
     1                     dd_s(i),ff_s(i)

            else
                write(6,20)n_sfc_obs,stations(i)(1:3)

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
