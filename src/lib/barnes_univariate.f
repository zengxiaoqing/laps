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


        subroutine maps_to_laps (maps,laps1,ni,nj,iden_ratio
     1                          ,r_missing_data)

!       Horizontal interpolation from sparse grid to laps grid
!       The input array only has values at the maps grid point locations.
!       The output array has the values filled in.
!       It just so happens the density array is filled in at the MAPS grid locations
!       Revised 1992 Steve Albers

        real*4 maps(ni,nj),laps1(ni,nj)


        DO j_laps=1,nj
        DO i_laps=1,ni

            i_maps=(i_laps-1)/iden_ratio * iden_ratio + 1
            j_maps=(j_laps-1)/iden_ratio * iden_ratio + 1

            if(i_maps .eq. ni)i_maps = i_maps - iden_ratio
            if(j_maps .eq. nj)j_maps = j_maps - iden_ratio

            fraci=(i_laps - i_maps) / float(iden_ratio)
            fracj=(j_laps - j_maps) / float(iden_ratio)

            Z1=maps(i_maps           , j_maps           )
            Z2=maps(i_maps+iden_ratio  , j_maps           )
            Z3=maps(i_maps+iden_ratio  , j_maps+iden_ratio)
            Z4=maps(i_maps           , j_maps+iden_ratio)
C
C INTERPOLATE FROM THE GRID TO THE STATION LOCATION
            if(       z1 .eq. r_missing_data
     1       .or.   z2 .eq. r_missing_data
     1       .or.   z3 .eq. r_missing_data
     1       .or.   z4 .eq. r_missing_data )then
                laps1(i_laps,j_laps) = r_missing_data
            else
                laps1(i_laps,j_laps)=Z1+(Z2-Z1)*fraci+(Z4-Z1)*fracj
     #                             -(Z2+Z4-Z3-Z1)*fraci*fracj
            endif

        enddo ! j_laps
        enddo ! i_laps

!       call smooth(laps1,laps2)

9999    Return
        END


      subroutine barnes_univariate_95(u                     ! Outputs
     1                     ,imax,jmax,kmax,grid_spacing_m   ! Inputs
     1                     ,uo,weight_3d,fnorm,n_fnorm      ! Inputs
     1                     ,l_analyze                       ! Input
     1                     ,density_array_in                ! Local
     1                     ,r0_array_out                    ! Output
     1                     ,n_obs_lvl,istatus)              ! Outputs

      implicit none

      include 'lapsparms.inc'

!     WINDPARMS.FOR
        real*4 hor_radius
        parameter (hor_radius = 35.)

        integer*4 MAX_PR,MAX_PR_LEVELS
        parameter (MAX_PR = 10)
        parameter (MAX_PR_LEVELS = 128)

!       Weights in wind analysis
        real*4 weight_meso
        parameter (weight_meso = 1.03)

        real*4 weight_sao
        parameter (weight_sao = 1.0)

        real*4 weight_pirep
        parameter (weight_pirep = 1.02)

        real*4 weight_prof
        parameter (weight_prof = 1.01)     ! Must be a unique value

        real*4 weight_radar
        parameter (weight_radar = 0.05)

        real*4 density_region_ratio
        parameter (density_region_ratio = 0.85)

        real*4 density_smooth_m
        parameter (density_smooth_m = 60000.)

!       Max Barnes radius (meters). Recommended values: 120000. to 300000.
!                                   Best value: 240000.
        real*4 r0_barnes_max_m
        parameter (r0_barnes_max_m = 240000.)

!       Weight for Model Background. Recommended values: 0. to 1e+30.
!       This will make the output values match the background if far from obs.
!       A value of zero means this parameter is not active.
        real*4 weight_bkg_const
        parameter (weight_bkg_const = 0.)
!       parameter (weight_bkg_const = 1.0e30)

!       1 = NO Successive correction, > 1 indicates successive correction
        integer*4 n_iter_wind
        parameter (n_iter_wind = 2)


!***********************END OF "WINDPARMS.FOR"*********************************


      integer*4 n_fnorm,max_obs
      parameter (max_obs = 4000)        ! max obs per level

      integer*4 imax,jmax,kmax
      integer*4 i,ii,il,j,jj,jl,k,istatus,ncnt,ncnt_total,n,i_subscript
      real*4 r0_norm,r0_norm_sq,r0_value,r0_value_min,ratio_norm_sq
      real*4 sum_u,sum_v,sumwt,weight
      real*4 density,density_max,fnorm_max,density_region_m
      real*4 grid_spacing_m,r0_barnes_max,r0_calc_min

      real*4 uo(imax,jmax,kmax),u(imax,jmax,kmax)
      integer*4 iob(max_obs),job(max_obs)
      real*4 fnorm(0:n_fnorm)
      real*4 weight_3d(imax,jmax,kmax)

      integer*4 n_obs_lvl(kmax),n_cross_in,n_density_smooth

      logical   l_analyze(kmax)
      real*4    r0_array_out(imax,jmax)

      real*4    density_array_in(imax,jmax)
      real*4    wt_lut(-(imax-1):(imax-1),-(jmax-1):(jmax-1))

      istatus = 0

      write(6,1234)
 1234 format(1x,'barnes_new called')


      r0_norm = 10. ! Grid point lengths
      r0_norm_sq = r0_norm**2

!     Minimum radius of influence in grid points (for solid data coverage)
      r0_value_min = 1.7

!     Maximum radius of influence in grid points (for no data coverage)
      r0_barnes_max = r0_barnes_max_m / grid_spacing_m

      n_density_smooth = density_smooth_m / grid_spacing_m

!     Check for possibility of illegal value in fnorm LUT
      ratio_norm_sq = r0_norm_sq / (r0_value_min**2)
      fnorm_max = (float(imax-1)**2 + float(jmax-1)**2) * ratio_norm_sq

      write(6,*)'r0_value_min,fnorm_max,n_fnorm'
     1  ,r0_value_min,fnorm_max,n_fnorm

      write(6,*)'r0_barnes_max',r0_barnes_max
      write(6,*)'n_density_smooth',n_density_smooth

      if(fnorm_max .gt. n_fnorm)then
          write(6,*)
     1           ' Error, fnorm_max is potentially TOO LARGE, '
     1          ,'increase n_fnorm allocation from ',n_fnorm
     1          ,' to > ',fnorm_max
          return
      endif


!     Create a lookup table for wt_lut
      do i = -(imax-1),(imax-1)
      do j = -(jmax-1),(jmax-1)
             wt_lut(i,j) = i*i+j*j
      enddo
      enddo

      density_max = 0.
      ncnt_total=0

!     Loop through vertical levels
      do k=1,kmax
          ncnt=0

!         Initialize density box array
          do j = 1,jmax
          do i = 1,imax
              density_array_in(i,j) = 0.0
!             if(k .eq. 1)density_array_in(i,j) = 1.0
          enddo ! i
          enddo ! j

!         Count obs and determine i,j of obs.
          do j=1,jmax
          do i=1,imax

              if(uo(i,j,k) .ne. r_missing_data)then

                  if(ncnt .ge. max_obs)then
                      write(6,*)' ERROR: Too many obs for barnes',ncnt
                      return
                  else

!                     Increment ob counter
                      ncnt=ncnt+1

                      iob(ncnt)=i
                      job(ncnt)=j
                      density_array_in(i,j) = 1.0

!                     if(k .eq. 13)write(6,71)ncnt,i,j
!       1                               ,uo(i,j,k),weight_3d(i,j,k)
!71                   format(1x,i4,2i3,f11.4,e11.4)


                  endif ! Not too many obs

              endif ! We have an ob

          enddo ! i
          enddo ! j

          n_obs_lvl(k) = ncnt
          ncnt_total = ncnt_total + ncnt

!         Calculate densities and create large density array
          if(l_analyze(k))then
              density_max = 0.

              density_region_m = density_region_ratio * r0_barnes_max_m

              n_cross_in = density_region_m / grid_spacing_m
              n_cross_in = min(n_cross_in,imax,jmax)
              n_cross_in = ((n_cross_in/2)*2) + 1 ! Must be ODD

!             n_cross_in = ((n_cross_in-1) / 2)

!             Calculate average density in big boxes centered on each grid pt
              call smooth_cross_laps(imax,jmax,1,imax,1,jmax
     1                                 ,density_array_in  ! I/O
     1                                 ,r0_array_out      ! Dummy in this case
     1                                 ,n_cross_in)

!             General smooth over of the density field
              do i = 1,n_density_smooth
                  call smooth_cross_laps(imax,jmax,1,imax,1,jmax
     1                                 ,density_array_in  ! I/O
     1                                 ,r0_array_out      ! Dummy in this case
     1                                 ,3)
              enddo

              do i = 1,imax
              do j = 1,jmax

                  density = density_array_in(i,j)

                  if(density .eq. 0.0)then
                      r0_value = r0_barnes_max
                  else
                      r0_value = min(r0_barnes_max,
     1                                   sqrt(1./density) * r0_value_min
     1)
                  endif

                  density_max = max(density,density_max)

!                 Radius value in grid points
                  r0_array_out(i,j) = r0_value

              enddo ! j
              enddo ! i

              if(density_max .gt. 0.)then
                  r0_calc_min = sqrt(1./density_max) * r0_value_min
              else
                  r0_calc_min = 0.
              endif

              write(6,402)k,ncnt,density_max,r0_calc_min,n_cross_in
402           format(' Analyzing: lvl k',i7,'  nobs_lvl',i7
     1                           ,' den/r0/ncross',f10.4,f10.2,i4)

              do j=1,jmax
              do i=1,imax
                  sum_u=0.
                  sumwt=0.

!                 Nominal radius_sq of influence divided by desired radius_sq of influence
                  ratio_norm_sq = r0_norm_sq / r0_array_out(i,j)**2

                  do n=1,ncnt
                      ii=iob(n)
                      jj=job(n)

                      i_subscript = 0.5+(wt_lut(i-ii,j-jj)*ratio_norm_sq
     1)
                      weight = fnorm(i_subscript)
     1                  * weight_3d(ii,jj,k) ! Obs are being weighted

                      sum_u=sum_u+(weight*uo(ii,jj,k))
                      sumwt=sumwt+weight

                      if(i .eq. 23 .and. abs(j-31) .le. 5
     1                   .and. (k .eq. 6 .or. k .eq. 3)   )then
                          write(6,81,err=82)i,j,n,ii,jj,uo(ii,jj,k)
     1          ,weight,sum_u/sumwt,r0_array_out(i,j)
 81                       format(1x,3i4,2i3,f11.4,e11.4,3f8.2)
 82                   endif

                  enddo ! n

                  if (sumwt.eq.0.)then
                      u(i,j,k) = r_missing_data
                  else
                      u(i,j,k)=sum_u/(sumwt + weight_bkg_const)
                  end if

              enddo ! i
              enddo ! j

          else ! if (.not. l_analyze) ...
              write(6,411)k
411           format(1x,i9,' Skipping 2nd Pass')
          endif ! l_analyze

      enddo ! k

      write(6,*)'ncnt_total',ncnt_total

      istatus = 1
      return
      end
