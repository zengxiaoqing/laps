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
        subroutine rdmeso(i4time,heights_3d
     1  ,N_MESO,maxstns
     1  ,lat,lon,meso_i,meso_j,meso_k,meso_u,meso_v,n_meso_obs
     1  ,grid_laps_wt,grid_laps_u,grid_laps_v
     1          ,ni,nj,nk,istatus)

!       1997 Jun    Ken Dritz      Added N_MESO and maxstns as dummy arguments,
!                                  making arrays with those dimensions
!                                  automatic.
!       1997 Jun    Ken Dritz      Removed include of 'lapsparms.for'.

!**************************************************************************

!       LAPS Grid Dimensions

        include 'windparms.inc' ! weight_meso

        real*4 lat(ni,nj)
        real*4 lon(ni,nj)

!       MESO

        integer meso_i(N_MESO) ! X meso coordinates
        integer meso_j(N_MESO) ! Y meso coordinates
        integer meso_k(N_MESO) ! Z meso coordinates
        real    meso_u(N_MESO) ! u meso component
        real    meso_v(N_MESO) ! v meso component

!       Laps Analysis Grids
        real grid_laps_wt(ni,nj,nk)
        real grid_laps_u(ni,nj,nk)
        real grid_laps_v(ni,nj,nk)

!***************************************************************************

        real*4 heights_3d(ni,nj,nk)

        character*13 filename13,c13_fname
        character*100 c100_fname

        character asc_tim_9*9,asc_tim_vol*24

        real*4 lat_s(maxstns), lon_s(maxstns), elev_s(maxstns)
        real*4 cover_s(maxstns), hgt_ceil(maxstns), hgt_low(maxstns)
        real*4 t_s(maxstns), td_s(maxstns), pr_s(maxstns), sr_s(maxstns)
        real*4 dd_s(maxstns), ff_s(maxstns), ddg_s(maxstns), ffg_s(maxst
     1ns)
        real*4 vis_s(maxstns)
c
        character stations(maxstns)*3, wx_s(maxstns)*8        ! c5_stamus
        character atime*24, infile*70
        character directory*50,ext*31

!       Declarations for new read_surface routine
!       New arrays.f reading in the SAO data from the LSO files
        real*4   pstn(maxstns),pmsl(maxstns),alt(maxstns),store_hgt(maxs
     1tns,5)
        real*4   ceil(maxstns),lowcld(maxstns),cover_a(maxstns),vis(maxs
     1tns)
     1                                          ,rad(maxstns)

        Integer*4   obstime(maxstns),kloud(maxstns),idp3(maxstns)

        Character   obstype(maxstns)*8
     1             ,store_emv(maxstns,5)*1,store_amt(maxstns,5)*4


        ext = 'msg'
        call get_directory(ext,directory,len_dir)

        c13_fname = filename13(i4time,ext(1:3))
        c100_fname = directory(1:len_dir)//c13_fname
        open(32,file=c100_fname,err=888,status='unknown')

        call make_fnam_lp(i4time,asc_tim_9,istatus)

        write(6,*)' Reading mesos: Calling read_sfc ',asc_tim_9
        ext = 'lso'
        call get_directory(ext,directory,len_dir) ! Returns top level directory

        c13_fname = filename13(i4time,ext(1:3))
        infile = directory(1:len_dir)//c13_fname

        write(6,*)infile

        call read_surface_old(infile,maxstns,atime,n_meso_g,n_meso_pos,
     &           n_sao_g,n_sao_pos_g,n_sao_b,n_sao_pos_b,n_obs_g,n_obs_p
     1os_g,
     &           n_obs_b,n_obs_pos_b,stations,obstype,lat_s,lon_s,
     &           elev_s,wx_s,t_s,td_s,dd_s,ff_s,ddg_s,
     &           ffg_s,pstn,pmsl,alt,kloud,ceil,lowcld,cover_a,rad,idp3,
     1store_emv,
     &           store_amt,store_hgt,vis,obstime,istatus)

100     write(6,*)'num_meso',n_meso_g

        if(n_obs_b .gt. maxstns)then
            write(6,*)' Too many stations or bad status',n_obs_b,maxstns
            istatus = 0
            return
        endif

        if(n_meso_g .gt. N_MESO)then
            write(6,*)' Too many mesonet stations',n_meso_g,N_MESO
            istatus = 0
            return
        endif

        if(istatus .ne. 1)then
            write(6,*)' Bad ISTATUS in read_sfc',istatus
            istatus = 0
            return
        endif

        write(6,12)
12      format(/'             Mapping meso Obs'
     1      /'   n Sta  i  j  k   u      v'
     1      ,'       dd     ff    ')

        n_meso_obs = 0

        do i = 1,n_obs_b ! num_meso

          if(  (     obstype(i)(1:4) .eq. 'MESO' 
     1          .or. obstype(i)(1:4) .eq. 'CDOT' )
     1                             .and. dd_s(i) .ge. 0.0)then

            n_meso_obs = n_meso_obs + 1

            ff_s(i) = ff_s(i) * .518

            call latlon_to_rlapsgrid(lat_s(i),lon_s(i),lat,lon,ni,nj
     1          ,ri,rj,istatus)

            meso_i(n_meso_obs) = nint(ri)
            meso_j(n_meso_obs) = nint(rj)

            call disp_to_uv(dd_s(i),ff_s(i)
     1                          ,meso_u(n_meso_obs),meso_v(n_meso_obs))

!           write(6,*)lat_s(i),lon_s(i),dd_s(i),ff_s(i),elev_s(i)

! ***       Remap meso observation to LAPS observation grid
!           In bounds?
            if(      meso_i(n_meso_obs) .ge. 1 
     1         .and. meso_i(n_meso_obs) .le. ni
     1         .and. meso_j(n_meso_obs) .ge. 1 
     1         .and. meso_j(n_meso_obs) .le. nj  )then

                rk = height_to_zcoord2(elev_s(i)
     1          ,heights_3d,ni,nj,nk
     1          ,meso_i(n_meso_obs),meso_j(n_meso_obs),istatus)

                if(istatus .ne. 1)return

                meso_k(n_meso_obs) = nint(rk)

!               Fix for stations that are to low for the grid
                if(meso_k(n_meso_obs) .eq. 0)meso_k(n_meso_obs) = 1

                write(32,*)ri-1.,rj-1.,rk-1.,dd_s(i),ff_s(i)

                k = meso_k(n_meso_obs)

                grid_laps_u(meso_i(n_meso_obs),meso_j(n_meso_obs),k)
     1                          = meso_u(n_meso_obs)

                grid_laps_v(meso_i(n_meso_obs),meso_j(n_meso_obs),k)
     1                          = meso_v(n_meso_obs)

                grid_laps_wt(meso_i(n_meso_obs),meso_j(n_meso_obs),k)
     1                          = weight_meso

                write(6,20)n_meso_obs,stations(i)(1:3),
     1                     meso_i(n_meso_obs),
     1                     meso_j(n_meso_obs),
     1                     meso_k(n_meso_obs),
     1                     meso_u(n_meso_obs),
     1                     meso_v(n_meso_obs),
     1                     dd_s(i),ff_s(i) ! ,azimuth,slant_range/1000.

            else ! Out of bounds
                write(6,20)n_meso_obs,stations(i)(1:3)

            endif ! In horizontal bounds

20          format(i4,1x,a3,3i3,2f7.1,2x,2f7.1,2x,2f7.1,2x,2f7.1)

          endif ! wind is reported

        enddo

        close(32)

        istatus = 1

        return

888     write(6,*)' Open error for MSG file'
        istatus = 0
        close(32)
        return

        end
