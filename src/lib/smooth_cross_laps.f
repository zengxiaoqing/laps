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

        subroutine smooth_cross_laps(ni,nj,i_l,i_h,j_l,j_h,r4_img
     1                              ,n_cross_in)

        integer n_cross_in    ! ODD
        real r4_img(ni,nj) ! Input/Output
        real r4_buf(ni,nj) ! Local use

        if(n_cross_in .gt. ni .or. n_cross_in .gt. nj)then
            write(6,*)' Error in smooth_cross_laps',n_cross_in,ni,nj
            write(6,*)' Aborting execution'
            stop
        endif

        n_cross = n_cross_in/2

!       Smooth in the I direction

        do j = j_l,j_h

            do i = i_l+n_cross,i_h-n_cross

                isum = 0
                rsum = 0

                do ii = i-n_cross,i+n_cross
                    isum = isum + 1
                    rsum = rsum + r4_img(ii,j)
                enddo ! ii

                r4_buf(i,j) = rsum / float(n_cross_in)

            enddo ! i

            do i = 1,i_l+n_cross-1
                r4_buf(i,j) = r4_buf(i_l+n_cross,j)
            enddo

            do i = i_h-n_cross+1,i_h
                r4_buf(i,j) = r4_buf(i_h-n_cross,j)
            enddo

        enddo ! j

!       Smooth in the J direction

        do i = i_l,i_h

            do j = j_l+n_cross,j_h-n_cross

                isum = 0
                rsum = 0
                do jj = j-n_cross,j+n_cross
                    isum = isum + 1
                    rsum = rsum + r4_buf(i,jj)
                enddo ! ii
                r4_img(i,j) = rsum / float(n_cross_in)

            enddo ! j

            do j = 1,j_l+n_cross-1
                r4_img(i,j) = r4_img(i,j_l+n_cross)
            enddo

            do j = j_h-n_cross+1,j_h
                r4_img(i,j) = r4_img(i,j_h-n_cross)
            enddo

        enddo ! i


        return
        end


