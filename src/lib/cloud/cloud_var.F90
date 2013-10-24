subroutine cloud_var(i4time,lat,lon,ni,nj,nk,kcld,heights_3d,temp_3d,t_gnd_k,cloudcvr,cld_hts,tb8_k,cloud_albedo_vis_a &
                    ,subpoint_lat_clo_vis,subpoint_lon_clo_vis,r_missing_data &    ! I 
                    ,di_dh,dj_dh)                                                  ! O

! This routine evaluates the cost function of the cloud cover field using simple forward models for visible and 11 micron satellite data

include 'trigd.inc'

use cloud_rad
#ifdef CRTM
    use prmtrs_stmas_cloud, only : satellite_obs
    use crtm_kmatrix, only : conf4lvd
!   Declarations for conf4lvd
    TYPE(SATELLITE_OBS) lvd
    character*200 coe_path,moist_unit
#endif


! Note kcld is on the cloud height grid
! Note nk is on the LAPS grid

real cloudcvr(ni,nj,kcld)              ! I/O
real cld_hts(kcld)                     ! I
real tb8_k(ni,nj)                      ! I
real t_gnd_k(ni,nj)                    ! I
real cloud_albedo_vis_a(ni,nj)         ! I

real heights_3d(ni,nj,nk)              ! I
real temp_3d(ni,nj,nk)                 ! I
real r_missing_data                    ! I

real lat(ni,nj)                        ! I
real lon(ni,nj)                        ! I
real subpoint_lat_clo_vis(ni,nj)       ! I
real subpoint_lon_clo_vis(ni,nj)       ! I
real di_dh(ni,nj)                      ! O
real dj_dh(ni,nj)                      ! O

real alt(ni,nj)                        ! L
real azi(ni,nj)                        ! L
real phase(ni,nj)                      ! L
real spec(ni,nj)                       ! L
real dx(ni,nj)                         ! L
real dy(ni,nj)                         ! L
real projrot_laps(ni,nj)               ! L
real cloud_frac_vis_a(ni,nj)           ! L
real t_effective                       ! L
real JCOST

integer max_cld_lyrs                   ! L
parameter (max_cld_lyrs = 100)

! Mainly local use in cvr_to_tb8_effective
real a(max_cld_lyrs)        ! Cloud fractions of layers
real f(max_cld_lyrs)        ! Apparent "cross-section" of cloud layers seen from above
integer ilyr(kcld)          ! Layer index for each cloud lvl (needs KCLOUD)

real di(kcld)               ! Parallax correction in I index for each cloud level
real dj(kcld)               ! Parallax correction in J index for each cloud level
real cldcvr_1d(kcld)        ! 1D cloud cover field, corrected for satellite parallax

! Use for CRTM
integer numgrid(4)

I4_elapsed = ishow_timer()

print*,'start subroutine cloud_var'

JCOST_TB8 = 0.

! Satellite geometry for parallax offset 
write(6,*)' Calling satgeom...'
range_m = 42155680.00
call satgeom(i4time,lat,lon,ni,nj &
            ,subpoint_lat_clo_vis,subpoint_lon_clo_vis,range_m,r_missing_data,Phase &
            ,Spec,alt,azi,istatus)
call get_grid_spacing_array(lat,lon,ni,nj,dx,dy)
call projrot_latlon_2d(lat,lon,ni,nj,projrot_laps,istatus)

! Adjust cloudcvr based on values of tb8 and cloud_frac_vis_a

! Call cvr_to_tb8_effective as a simple forward model (for each grid column)
do i = 1,ni
do j = 1,nj

!   Calculate parallax offset (sat/lvd grid index minus analysis grid index)
    if(alt(i,j) .gt. 0.)then
        ds_dh = tand(90. - alt(i,j))
        azi_grid = azi(i,j) - projrot_laps(i,j)
        do k = 1,kcld
            di_dh(i,j) = (ds_dh / dx(i,j)) * (-sind(azi_grid))
            dj_dh(i,j) = (ds_dh / dy(i,j)) * (-cosd(azi_grid))
            di(k) = cld_hts(k) * di_dh(i,j)
            dj(k) = cld_hts(k) * dj_dh(i,j)
        enddo ! k
    else
        di = 0. 
        dj = 0.
        di_dh(i,j) = 0.
        dj_dh(i,j) = 0.
    endif

    intvl = max(ni/41,1)
    if(j .eq. nj/2 .AND. i .eq. (i/intvl)*intvl)then
        idebugsub = 1
    else
        idebugsub = 0
    endif

    if(i .eq. ni/2 .AND. j .eq. nj/2)then
        idebug = 1
    else
        idebug = 0
    endif

    if(idebugsub .eq. 1)then
        write(6,101)i,j,subpoint_lon_clo_vis(i,j),alt(i,j),azi(i,j),cld_hts(20),di(20),dj(20)
101     format(' sat i/j/sub/alt/azi/ht/di/dj = ',2i4,6f10.3)
    endif

!   Determine 1D cloud cvr in the column corrected for parallax to compare with
!   the satellite data from point i,j
    do k = 1,kcld
        ip = i - nint(di(kcld)) ; ip = max(ip,1); ip = min(ip,ni)
        jp = j - nint(dj(kcld)) ; jp = max(jp,1); jp = min(jp,nj)
        cldcvr_1d(k) = cloudcvr(ip,jp,k)
    enddo ! k

    call cvr_to_tb8_effective(kcld,temp_3d,nk,i,j,ni,nj,a      &
                             ,f,ilyr,cldcvr_1d,cld_hts,t_gnd_k(i,j)  &
                             ,heights_3d,t_effective,nlyr      &
                             ,idebug,istatus)

    if(idebug .eq. 1 .OR. istatus .ne. 1)then
        write(6,*)'tb8 fwd mdl - max cldcvr: ',i,j,maxval(cloudcvr(i,j,:))
        write(6,*)'cldcvr column: ',i,j,cloudcvr(i,j,:)
        write(6,*)'t_effective/t_gnd_k/tb8: ',t_effective,t_gnd_k(i,j),tb8_k(i,j)
        write(6,*)'ilyr: ',ilyr(:)
        write(6,*)'nlyr,f: ',nlyr,(f(l),l=1,nlyr)
    endif

    JCOST_TB8 = JCOST_TB8 + abs(tb8_k(i,j) - t_effective)
enddo
enddo

print*,'Mean value of JCOST_TB8 = ',JCOST_TB8/float(ni*nj)

! Use cloud_frac_vis_a as an upper bound on the cloud fractions

JCOST_VIS = 0.

do i = 1,ni
do j = 1,nj
    cloud_frac_vis_a(i,j) = cloud_albedo_vis_a(i,j)
    if(cloud_frac_vis_a(i,j) .ne. r_missing_data)then
        call albedo_to_clouds(cloud_albedo_vis_a(i,j),cloud_rad_trans,cloud_od,cloud_opacity)

!       Determine 1D cloud cvr in the column corrected for parallax to compare
!       with the satellite data from point i,j
        do k = 1,kcld
            ip = i - nint(di(kcld)) ; ip = max(ip,1); ip = min(ip,ni)
            jp = j - nint(dj(kcld)) ; jp = max(jp,1); jp = min(jp,nj)
            cldcvr_1d(k) = cloudcvr(ip,jp,k)
        enddo ! k

        cvr_max_vert = maxval(cldcvr_1d(:))
        if(cloud_frac_vis_a(i,j) .lt. cvr_max_vert)then
            JCOST_VIS = JCOST_VIS + (cvr_max_vert - cloud_frac_vis_a(i,j))
        endif
    endif
enddo
enddo

print*,'Mean value of JCOST_VIS = ',JCOST_VIS/float(ni*nj)


print*,'finished subroutine cloud_var'

I4_elapsed = ishow_timer()

#ifdef CRTM
    if(.false.)then
        call conf4lvd(lvd,numgrid,'imgr_g13',coe_path,moist_unit,badflag)
    endif
#endif

end subroutine cloud_var

