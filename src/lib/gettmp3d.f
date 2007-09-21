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
        subroutine get_temp_3d(i4time_needed,i4time_nearest,iflag,
     1                          imax,jmax,kmax,field_3d,istatus)

cdoc    Returns 3-D temperature (or ht) field within a time window. Has options
cdoc    for theta and/or balanced fields depending on the 'iflag' input.

!       Steve Albers            2000

!       iflag = 0 (Potential Temperature K)
!       iflag = 1 (Temperature K)
!       iflag = 2 (Height m)
!       iflag = 3 (Balanced Pot'l Temperature K)
!       iflag = 4 (Balanced Temperature K)
!       iflag = 5 (Balanced Height m)

        logical ltest_vertical_grid

        character*150 DIRECTORY
        character*31 EXT

        character*125 comment_3d(kmax)
        character*10 units_3d(kmax)
        character*3 var_t
        integer LVL_3d(kmax)
        character*4 LVL_COORD_3d(kmax)

        real field_3d(imax,jmax,kmax)

        character*255 c_filespec

        O_K(T_K,P_PA)   =   O( T_K-273.15 , P_PA/100. )  + 273.15
!       TDA_K(T_K,P_PA) = TDA( T_K-273.15 , P_PA/100. )  + 273.15

        write(6,*)
        write(6,*)' Subroutine get_temp_3d: iflag=',iflag

        EXT = 'lt1'

        if(iflag .lt. 3)then ! 'temp.exe' analysis
            call get_directory(ext,directory,len_dir)
        else
            call get_directory('balance',directory,lend)
            directory=directory(1:lend)//'lt1/'
        endif

        if(iflag .eq. 2 .or. iflag .eq. 5)then
            var_t = 'HT'
        else
            var_t = 'T3'
        endif

        call get_3dgrid_dname(directory
     1           ,i4time_needed,1000000,i4time_nearest
     1           ,ext,var_t,units_3d
     1           ,comment_3d,imax,jmax,kmax,field_3d,istatus)       


!       Test for bad temperatures
        if(var_t .eq. 'T3')then
            if(field_3d(29,29,17) .lt. 173.)istatus = 0
        endif

        IF(ISTATUS .ne. 1)THEN
            write(6,*)' Error Reading in 3d lt1 field ',var_t
        endif

!       Convert from T to Theta
        do i = 1,imax
        do j = 1,jmax

        if(iflag .eq. 0 .or. iflag .eq. 3)then
            do k = 1,kmax
                theta = O_K(field_3d(i,j,k),zcoord_of_level(k))
                field_3d(i,j,k) = theta
            enddo ! k
        endif

        enddo ! j
        enddo ! i

        return

        end

