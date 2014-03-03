
        subroutine hsect_img(i4time_ref,lun,NX_L,NY_L,NZ_L
     1                      ,r_missing_data,laps_cycle_time
     1                      ,plot_parms,namelist_parms,ifield_found)

        real field_2d(NX_L,NY_L)
        character*31 ext
        character*4 var_2d
        character*20 units_2d

        write(6,*)' Subroutine hsect_img'

!       Read in cloud albedo from LIL file
        var_2d = 'CLA'
        vel = 0
        ext = 'lil'
        call get_laps_2dgrid(i4time_ref,86400,i4time_cloud,
     1                       ext,var_2d,units_2d,comment_2d,
     1                       NX_L,NY_L,field_2d,0,istatus)

!       Write as PNG
        write(6,*)' Write all sky cyl ppm file '
!       call writeppm3Matrix(
!    1             isky_rgb_cyl(0,:,:),isky_rgb_cyl(1,:,:)
!    1            ,isky_rgb_cyl(2,:,:),'allsky_rgb_cyl_'//clun)


        return
        end
