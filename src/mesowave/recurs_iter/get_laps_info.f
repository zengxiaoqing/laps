        subroutine get_laps_info(nx,ny,stanlat,
     &   stanlat2,stanlon
     &   ,laps_cycle_time,badflag,maxstations
     &   ,maproj)

        character*6 maproj 
        integer nx,ny,laps_cycle_time,maxstations
        real stanlat,stanlat2,stanlon,badflag
        
c
        call get_grid_dim_xy(nx,ny,istatus)
        if(istatus .ne. 1) print*, 'Error getting nx,ny' 
        call get_standard_latitudes(stanlat,stanlat2,istatus)
        if(istatus .ne. 1) print*, 'Error getting latitudes' 
        call get_standard_longitude(stanlon,istatus)
        if(istatus .ne. 1) print*, 'Error getting longitude'
        call get_laps_cycle_time(laps_cycle_time,istatus)
        if(istatus .ne. 1) print*, 'Error getting cycle time'
        call get_r_missing_data(badflag,istatus) 
        if(istatus .ne. 1) print*, 'Error getting badflag'
        call get_maxstns(maxstations,istatus)
        if(istatus .ne. 1) print*, 'Error getting maxstations'
        call get_c6_maproj(mapproj,istatus)
        if(istatus .ne. 1) print*, 'Error getting map projection'
     
        return
        end