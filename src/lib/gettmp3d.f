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
     1                          imax,jmax,kmax,temp_3d,istatus)

!       Steve Albers            1990

!       iflag = 0 (Potential Temperature K), iflag = 1 (Temperature K)

        include 'lapsparms.inc'

        character*50 DIRECTORY
        character*31 EXT

        character*125 comment_3d(NZ_L_MAX)
        character*10 units_3d(NZ_L_MAX)
        character*3 var_t(NZ_L_MAX)
        integer*4 LVL_3d(NZ_L_MAX)
        character*4 LVL_COORD_3d(NZ_L_MAX)

        real*4 temp_3d(imax,jmax,kmax)

        real*4 theta(NZ_L_MAX)

!       real*4 temp_std(NZ_L_MAX)/283.15,279.55,275.80,272.75,270.20,
!       1       267.45,263.10,258.65,253.75,248.20,242.85,235.35,
!       1       226.55,223.65,225.35,221.05,221.15/

        character*255 c_filespec

        O_K(T_K,P_PA)   =   O( T_K-273.15 , P_PA/100. )  + 273.15
!       TDA_K(T_K,P_PA) = TDA( T_K-273.15 , P_PA/100. )  + 273.15

        write(6,*)
        write(6,*)' Getting 3_D Temperature Analysis (LT1 or equivalent)
     1'

        EXT = 'lt1'
        call get_directory(ext,directory,len_dir)
        c_filespec = directory(1:len_dir)//'*.'//ext(1:3)

        call get_file_time(c_filespec,i4time_needed,i4time_nearest)

        do k = 1,kmax
            units_3d(k)   = 'K'
            if(vertical_grid .eq. 'HEIGHT')then
                lvl_3d(k) = zcoord_of_level(k)/10
                lvl_coord_3d(k) = 'MSL'
            elseif(vertical_grid .eq. 'PRESSURE')then
              if(i4time_nearest .gt. 963350000)then
                lvl_3d(k) = nint(zcoord_of_level(k))/100
              else ! Read offset old files
                lvl_3d(k) = nint(zcoord_of_level(k+1))/100
              endif

              lvl_coord_3d(k) = 'MB'

            endif

            var_t(k) = 'T3' ! newvar = 'T3', oldvar = 'T'

        enddo ! k

!       write(6,*)' LVL_3d = ',LVL_3d

!       Read in 3d T array
        CALL READ_LAPS_DATA(I4TIME_NEAREST,DIRECTORY,
     1          EXT,imax,jmax,kmax,kmax,
     1          var_t,LVL_3d,LVL_COORD_3d,UNITS_3d,COMMENT_3d,
     1          temp_3d,ISTATUS)

!       ISTATUS = .false. ! Just for testing

!       Test for bad temperatures
        if(temp_3d(29,29,17) .lt. 173.)istatus = 0

        IF(ISTATUS .ne. 1)THEN
            write(6,*)' Error Reading in 3d Temperature Analysis'
!       1            ,', Using Canned Sounding'
!           do k = 1,kmax
!           do j = 1,jmax
!           do i = 1,imax
!               temp_3d(i,j,k) = temp_std(k)
!           enddo
!           enddo
!           enddo
        endif


!       Convert from T to Theta
        do i = 1,imax
        do j = 1,jmax

        if(iflag .eq. 0)then
            do k = 1,kmax
                theta(k) = O_K(temp_3d(i,j,k),zcoord_of_level(k))
                temp_3d(i,j,k) = theta(k)
            enddo ! k
        endif

        enddo ! j
        enddo ! i

        return

        end

