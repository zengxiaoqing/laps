
        subroutine get_vis(i4time,solar_alt,l_use_vis,lat          ! Input
     1                    ,i4_sat_window,i4_sat_window_offset      ! Input
     1                    ,cloud_frac_vis_a,albedo,ihist_alb       ! Output
     1                    ,ni,nj,nk,r_missing_data                 ! Input
!    1                    ,istat_vis_a                             ! Output
     1                    ,istatus)                                ! Output

!       Steve Albers 1997

        integer*4 ihist_alb(-10:20)
        integer*4 ihist_frac_sat(-10:20)
        integer*4 istat_vis_a(ni,nj)            ! Cloud mask based on VIS
        real*4 albedo(ni,nj)

        real*4 lat(ni,nj)
        real*4 sfc_albedo(ni,nj)

!       This stuff is for reading VIS data from LVD file
        real*4 solar_alt(ni,nj)
        real*4 cloud_frac_vis_a(ni,nj)
        integer*4 mxstn
        parameter (mxstn = 100)       ! max number of "stations" in data file
        character*9 filename
        character*31 ext

        character*3 lvd_ext
        data lvd_ext /'lvd'/

        character var*3,comment*125,units*10
        logical l_use_vis

!       Initialize histograms
        do i = -10,20
            ihist_alb(i) = 0
            ihist_frac_sat(i) = 0
        enddo ! i

!       Initialize
        do i = 1,ni
        do j = 1,nj
            albedo(i,j) = r_missing_data
            cloud_frac_vis_a(i,j) = r_missing_data
            istat_vis_a(i,j) = 0
        enddo ! j
        enddo ! i

        call get_sfc_albedo(ni,nj,lat,r_missing_data               ! I
     1                     ,sfc_albedo,istat_sfc_albedo)           ! O

        n_missing_albedo = ni*nj

!       Determine whether to use VIS / ALBEDO data
        if(.not. l_use_vis)then
            write(6,*)' Warning: l_use_vis set to not use vis data'
            istatus = 0
            return
        endif

!       Read in albedo data
        write(6,*)' Getting the VIS data from LVD file'
        call make_fnam_lp(i4time,filename,istatus)

        ext = lvd_ext
        var = 'ALB'
        ilevel = 0
        call get_laps_2dvar(i4time+i4_sat_window_offset,i4_sat_window
     1                     ,i4time_nearest,EXT,var,units
     1                     ,comment,ni,nj,albedo,ilevel,istatus)
        write(6,*)' istatus from albedo data = ',istatus
        if(istatus .ne. 1 .and. istatus .ne. -1)then
            write(6,*)' Warning: could not read albedo - '
     1               ,'return from get_vis'
            istatus = 0
            return
        endif

        n_missing_albedo = 0

!       Horizontal array loop
        do i = 1,ni
        do j = 1,nj

!         We will now only use the VIS data if the solar alt exceeds 15 deg
!         11 degrees now used to allow 15 min slack in data availability
          if(solar_alt(i,j) .lt. 11.0)then
              if(albedo(i,j) .ne. r_missing_data)then
                  write(6,*)' Error -  albedo not missing:'
     1                     ,solar_alt(i,j)
                  stop
              endif
!             albedo(i,j) = r_missing_data
          endif

          if(albedo(i,j) .ne. r_missing_data)then

!           Translate the albedo into cloud fraction

!           Store histogram information for satellite data
            iscr_alb  = nint(albedo(i,j)*10.)
            iscr_alb  = min(max(iscr_alb,-10),20)
            ihist_alb(iscr_alb) = ihist_alb(iscr_alb) + 1

            clear_albedo = .2097063
            cloud_frac_vis = albedo_to_cloudfrac2(clear_albedo
     1                                           ,albedo(i,j))

!           Fudge the frac at low solar elevation angles
!           Note the ramp extrapolates down to 9 deg to account for slight
!           errors in determining the solar elevation
            if(.false.)then
!           if(solar_alt(i,j) .lt. 20. .and. solar_alt(i,j) .ge. 9.)then
                frac = (20. - solar_alt(i,j)) / 10.
                term1 = .13 * frac
                term2 = 1. + term1
                cloud_frac_vis = (cloud_frac_vis + term1) * term2
            endif

            iscr_frac_sat = nint(cloud_frac_vis*10.)
            iscr_frac_sat = min(max(iscr_frac_sat,-10),20)
            ihist_frac_sat(iscr_frac_sat) = 
     1      ihist_frac_sat(iscr_frac_sat) + 1

!           Make sure satellite cloud fraction is between 0 and 1
            if(cloud_frac_vis .le. 0.0)cloud_frac_vis = 0.0
            if(cloud_frac_vis .ge. 1.0)cloud_frac_vis = 1.0

            cloud_frac_vis_a(i,j) = cloud_frac_vis

!           Is there enough of a signal from the VIS to say a cloud is present?
            if(       cloud_frac_vis_a(i,j) .gt. 0.5 
     1          .and. sfc_albedo(i,j) .ne. r_missing_data
     1          .and. sfc_albedo(i,j) .le. 0.3              )then
                istat_vis_a(i,j) = 1
            endif

          else
            n_missing_albedo =  n_missing_albedo + 1
            cloud_frac_vis_a(i,j) = r_missing_data

          endif

        enddo ! i
        enddo ! j

        write(6,*)
        write(6,*)' N_MISSING_ALBEDO = ',n_missing_albedo
        write(6,*)

        if(n_missing_albedo .eq. ni*nj)then ! Return with status = 0
            write(6,*)' All albedos were missing - return from get_vis'
            istatus = 0
            return
        endif

        write(6,*)'              HISTOGRAMS'
        write(6,*)' I          ',
     1  ' Albedo  Cld Frac Sat'
        do i = -5,15
            write(6,11)i,ihist_alb(i),ihist_frac_sat(i)
11          format(i4,i12,i12,i12,i12)
        enddo ! i

        write(6,*)
        istatus = 1

        return
        end

        subroutine get_sfc_albedo(ni,nj,lat,r_missing_data           ! I
     1                           ,sfc_albedo,istat_sfc_albedo)       ! O

        character*3 var
        real*4 lat(ni,nj)
        real*4 sfc_albedo(ni,nj), static_albedo(ni,nj)

        var = 'ALB'

        call read_static_grid(ni,nj,var,static_albedo,istatus)

        do i = 1,ni
        do j = 1,nj
            sfc_albedo(i,j) = r_missing_data
        enddo ! j
        enddo ! i

        return
        end

        function albedo_to_cloudfrac2(clear_albedo,albedo)

        cloud_albedo = .4485300

        arg = albedo

        call stretch2(clear_albedo,cloud_albedo,0.,1.,arg)

        albedo_to_cloudfrac2 = arg

        return
        end
