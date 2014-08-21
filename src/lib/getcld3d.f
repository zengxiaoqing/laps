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

        subroutine get_clouds_3dgrid(i4time_needed,i4time_nearest
     1          ,imax,jmax,kcloudin,ext ! ,b_cldcv
     1                          ,clouds_3d,cld_hts,cld_pres,istatus)

!       Steve Albers            1990
!       Steve Albers            1991     (b_cldcv now is an input with adj dim)
!       Steve Albers     30 Oct 1991     cld_pres now has adjustable dims
!       Steve Albers      6 Nov 1991     ext is passed in
!	Linda Wharton    26 Oct 1998	 removed units_2d and comment_2d var not used

!       This is a newer version of 'get_clouds_3d' that reads the heights and
!       approximate Pressures directly from the cloud grid

!       ARGUMENT LIST
!       i4time_needed      Input   Desired i4time
!       i4time_nearest     Output  Actual i4time of nearest file
!       imax,jmax          Input   LAPS horizontal grid dimensions
!       kcloudin           Input   Number of cloud vertical levels (Dense grid)
!                                  This MUST be equal to value in laps_cloud.inc
!                                  if file is before i4time = 969476400
        character*31 ext ! Input   File extension (normally 'LC3')
!       byte b_cldcv(imax,jmax,kcloudin)    ! Output   Clouds in bytes
        real clouds_3d(IMAX,JMAX,KCLOUDIN)! Output 3D Array of Cloud Cover (0-1.2)
        real cld_pres(KCLOUDIN) ! Output Pressures of the vertical grid (Approx)
                                  ! In PASCALS

        include 'laps_cloud.inc'

!       These two arrays deal with old cloud data files only
        real cld_hts_old(KCLOUD)

        data cld_hts_old/1200.,1300.,1400.,1500.,1600.,1700.,1800.,
     11900.,2000.,2200.,2400.,2500.,2600.,2800.,3000.,3200.,
     23400.,3600.,3800.,4000.,4500.,5000.,5500.,6000.,6500.,
     37000.,7500.,8000.,8500.,9000.,9500.,10000.,11000.,12000.,
     413000.,14000.,15000.,16000.,17000.,18000.,19000.,20000./


        real cld_hts_new(KCLOUD)

        data cld_hts_new/1200.,1300.,1400.,1500.,1600.,1700.,1800.,
     11900.,2000.,2100.,2200.,2400.,2600.,2800.,3000.,3200.,
     23400.,3600.,3800.,4000.,4500.,5000.,5500.,6000.,6500.,
     37000.,7500.,8000.,8500.,9000.,9500.,10000.,11000.,12000.,
     413000.,14000.,15000.,16000.,17000.,18000.,19000.,20000./

        character*255 c_filespec
        character*13 filename13
        character*150 directory

!       Stuff for readlapsdata
        character*125 comment_3d(KCLOUD)
        character*10 units_3d(KCLOUD)

        character*3 var_3d(KCLOUD),var_2d
        data var_2d/'LC3'/

        integer LVL_3d(KCLOUD)
        character*4 LVL_COORD_3d(KCLOUD)

        logical l_packed_data
        data l_packed_data /.false./

        common /laps_diag/ no_laps_diag

        if(no_laps_diag .eq. 0)write(6,*)' Getting Cloud Analysis'

        call get_directory(ext,directory,len_dir)
        c_filespec = directory(1:len_dir)//'*.'//ext(1:3)

        call get_file_time(c_filespec,i4time_needed,i4time_nearest)

        if(i4time_nearest .eq. 0)then
            write(6,*)' Error: No cloud (',ext(1:3),') files available'
            istatus = 0
            return
        endif

!       Fill in cloud pressures array with standard atmosphere values
!       This is relavant for old files only

        i4time_switch = 969476400

        if(i4time_nearest .lt. i4time_switch)then
            if(KCLOUD .ne. KCLOUDIN)then
                write(6,*)' Error: Incorrect vertical cloud dimension'
                istatus = 0
                return
            endif
            do i = 1,KCLOUD
                cld_hts(i) = cld_hts_old(i)
            enddo
        else
          if(KCLOUD .eq. KCLOUDIN)then
            do i = 1,KCLOUD
                cld_hts(i) = cld_hts_new(i)
            enddo
          endif
        endif

        if(KCLOUD .eq. KCLOUDIN)then
          do k=1,KCLOUD
            cld_pres(k) = ztopsa(cld_hts(k)) * 100.
          enddo ! k
        endif

!       Read in cloud grid

        if(l_packed_data)then

            if(no_laps_diag .eq. 0)
     1       write(6,*)' Reading in byte array from ',ext(1:3),' directo
     1ry'
            open(21,file=directory(1:len_dir)//filename13(i4time_nearest
     1,ext(1:3))
     1  ,access='sequential',form='unformatted',status='old',err=998)
!           read(21,err=998)(((b_cldcv(i,j,k),i=1,imax),j=1,jmax),k=1,KCLOUDIN)
            read(21,err=100)(cld_hts(k),k=1,KCLOUDIN),(cld_pres(k),k=1,K
     1CLOUDIN)

            if(no_laps_diag .eq. 0)write(6,*)' Updated cloud heights/pre
     1ssures'
            close(21)

            goto90

100         write(6,*)' Using standard atmosphere for cloud pressure arr
     1ay'
            close(21)

!           Convert cloud cover array from byte format
90          do k = 1,KCLOUDIN
            do j = 1,jmax
            do i = 1,imax
!               if(b_cldcv(i,j,k) .ne. -1)then
!                   i_temp = b_cldcv(i,j,k)
!                   clouds_3d(i,j,k) = i_temp / 100.
!               else
!                   clouds_3d(i,j,k) = 0.0
!               endif
            enddo ! i
            enddo ! j
            enddo ! k

            istatus = 1
            return

998         write(6,*)' Error reading cloud (LC3) file'
            istatus = 0
            return

        else ! USE READLAPSDATA

            do k = 1,KCLOUDIN
                units_3d(k) = '    '
                lvl_3d(k) = k
                lvl_coord_3d(k) = '   '
                var_3d(k) = var_2d
            enddo ! k

            CALL READ_LAPS_DATA(I4TIME_NEAREST,DIRECTORY,EXT,imax,jmax,
     1  kcloudin,KCLOUD,VAR_3D,LVL_3D,LVL_COORD_3D,UNITS_3D,
     1                     COMMENT_3D,clouds_3d,ISTATUS)

            if(istatus .ne. 1)then
                write(6,*)' Error reading cloud (LC3) file'
                return
            endif

!           Decode cloud heights and pressures
            do k = 1,KCLOUDIN
                read(comment_3d(k),1,err=999)cld_hts(k),cld_pres(k)
1               format(2e20.8)
            enddo ! k

        endif

        return ! Normal return

999     write(6,*)' Error reading comment field in LC3 file'
        istatus = 0
        return

        end

        subroutine get_clouds_3dpres(i4time_needed,imax,jmax,kmax
     1  ,kcloud,clouds_3d_pres,clouds_3d_height ! ,b_cldcv
     1  ,cld_hts,cld_pres_1d,cld_pres_3d,istatus)

        integer i4time_needed                      ! Input
        integer imax,jmax,kmax                     ! Input (LAPS Dims)
        integer kcloud                             ! Input (normally 42)
        real clouds_3d_pres(imax,jmax,kmax)        ! Output
        real clouds_3d_height(imax,jmax,kcloud)    ! Local Dummy Array
!       byte   b_cldcv(imax,jmax,kcloud)             ! Local Dummy Array
        real cld_hts(kcloud)                       ! Local Dummy Array
        real cld_pres_1d(kcloud)                   ! Local Dummy Array
        real cld_pres_3d(imax,jmax,kcloud)         ! Local Dummy Array
        integer istatus                            ! Output

        character*31 ext

        ext = 'lc3'

        call get_clouds_3dgrid(i4time_needed,i4time_nearest
     1          ,imax,jmax,kcloud,ext ! ,b_cldcv
     1          ,clouds_3d_height,cld_hts,cld_pres_1d,istatus)

        if(istatus .ne. 1)return

        if(i4time_needed .ne. i4time_nearest)then
            istatus = 0
            return
        endif

!       Create cld_pres_3d array
        do k = 1,kcloud
        do j = 1,jmax
        do i = 1,imax
            cld_pres_3d(i,j,k) = cld_pres_1d(k)
        enddo ! i
        enddo ! j
        enddo ! k

        call interp_height_pres(imax,jmax,kmax,kcloud
     1  ,clouds_3d_pres,clouds_3d_height,cld_pres_3d,istatus)

        return
        end


        subroutine interp_height_pres(imax,jmax,kmax,kcloud
     1  ,clouds_3d_pres,clouds_3d_height,cld_pres_3d,istatus)

        integer imax,jmax,kmax                     ! Input (LAPS Dims)
        integer kcloud                             ! Input (normally 42)
        real clouds_3d_pres(imax,jmax,kmax)        ! Output
        real clouds_3d_height(imax,jmax,kcloud)    ! Input
        real cld_pres_3d(imax,jmax,kcloud)         ! Input
        integer istatus                            ! Output

!       Interpolate from height grid to pressure grid (Generating inputs can
!       lead to slowness or inaccuracies)
        do lvl = 1,kmax
            rlvl = pressure_of_level(lvl)

            do j = 1,jmax
            do i = 1,imax

!             Efficiency speedup
              if(i .eq. 1 .or. i .eq. i/5*5)then

!               Default if out of bounds
                kht_low = 1
                kht_high = 1
                frac_low = 1.0
                frac_high = 0.0

                do kht = KCLOUD-1,1,-1
                    if(cld_pres_3d(i,j,kht  ) .ge. rlvl .and.
     1               cld_pres_3d(i,j,kht+1) .le. rlvl           )then
                        kht_low = kht
                        kht_high = kht_low + 1
                        frac_high = (rlvl - cld_pres_3d(i,j,kht)) /
     1                        (cld_pres_3d(i,j,kht+1) - cld_pres_3d(i,j,
     1kht))
                        frac_low  = 1.0 - frac_high
                        goto 10 ! We found the cloud level, kick out of loop
                    endif
                enddo ! kht

10            endif

!d            write(6,1)lvl,kht_low,kht_high,frac_low,frac_high
!d      1            ,cld_pres_3d(i,j,kht_low),rlvl,cld_pres_3d(i,j,kht_high)
!d1           format(1x,' Reading Clouds LVL ',i5,3x,2i5,3x,2f8.3,2x,3f8.2)

              clouds_3d_pres(i,j,lvl)
     1                = clouds_3d_height(i,j,kht_low ) * frac_low
     1                      + clouds_3d_height(i,j,kht_high) * frac_high
            enddo ! i
            enddo ! j

        enddo ! lvl

        return
        end


        subroutine interp_height_pres_fast(imax,jmax,kmax,kcloud
     1  ,clouds_3d_pres,clouds_3d_height,heights_3d,cld_hts,istatus)

        integer imax,jmax,kmax                     ! Input (LAPS Dims)
        integer kcloud                             ! Input (normally 42)
        real clouds_3d_pres(imax,jmax,kmax)        ! Output
        real clouds_3d_height(imax,jmax,kcloud)    ! Input
        real heights_3d(imax,jmax,kmax)            ! Input
        real cld_hts (kcloud)
        integer istatus                            ! Output

!       Interpolate from height grid to pressure grid (This is fast and
!       accurate when heights are supplied)

        do j = 1,jmax
        do i = 1,imax

            kstart = 1

            do lvl = 1,kmax

                do kht = kstart,KCLOUD-1
                    if(cld_hts(kht+1) .ge. heights_3d(i,j,lvl)  )then
                        kht_low  = kht
                        kht_high = kht_low + 1

                        if(cld_hts(kht) .le. heights_3d(i,j,lvl))then
                            frac_high = (heights_3d(i,j,lvl) - cld_hts(k
     1ht)) /
     1                                (cld_hts(kht+1)      - cld_hts(kht
     1))
                            frac_low  = 1.0 - frac_high
                        else ! Take care of lower boundary effects
                            kht_low = 1
                            kht_high = 1
                            frac_low = 1.0
                            frac_high = 0.0
                        endif

                        goto 10 ! We found the cloud level, kick out of loop
                    endif
                enddo ! kht

10            kstart = kht

              clouds_3d_pres(i,j,lvl)
     1                = clouds_3d_height(i,j,kht_low ) * frac_low
     1                + clouds_3d_height(i,j,kht_high) * frac_high

              if(i .eq. imax/2 .and. j .eq. jmax/2)then
                  write(6,1)lvl,kht_low,kht_high,frac_low,frac_high
     1                     ,clouds_3d_height(i,j,kht_low)
     1                     ,clouds_3d_height(i,j,kht_high)
     1                     ,heights_3d(i,j,lvl)
     1                     ,clouds_3d_pres(i,j,lvl)
 1                format(1x,' Reading Clouds LVL ',i5,3x,2i5,3x,2f8.3
     1                                            ,2x,2f8.2,f8.1,f8.2)
              endif

            enddo ! lvl
        enddo ! i
        enddo ! j

        return
        end
