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
        subroutine insert_pireps(i4time,cldcv,cld_hts,wtcldcv
     1        ,cld_snd,wt_snd,i_snd,j_snd,n_cld_snd,max_cld_snd
     1  ,lat,lon,ni,nj,nk,ix_low,ix_high,iy_low,iy_high,max_pireps
     1  ,istatus)

        real*4 lat(ni,nj),lon(ni,nj)

        real*4 cldcv(ni,nj,nk)
        real*4 wtcldcv(ni,nj,nk)
        real*4 cld_hts(nk)

!       Arrays for cloud soundings
        real*4 cld_snd(max_cld_snd,nk)
        real*4 wt_snd(max_cld_snd,nk)
        integer*4 i_snd(max_cld_snd)
        integer*4 j_snd(max_cld_snd)

        character*9  asc9_tim_pirep,asc9_tim_rcvd
        character*80 string
        character*13 filename13

        character*31    ext
        integer*4       len_dir

        logical l_good_pirep

!       1 Feb 1996     Steve Albers        Remove spread( calls

        max_layers = 3  ! Defined in PIN file and NetCDF RAW data

        lun = 11
        ext = 'pin'
        call open_lapsprd_file(lun,i4time,ext,istatus)
        if(istatus .ne. 1)go to 998

        n_cloud = 0
        num_pirep = 0
        num_good_pirep = 0

5       read(11,101,end=900,err=50)string(1:6)
101     format(a6)

50      if(string(2:5) .eq. 'Time')then
!           a9time = string(30:39)
            read(11,151)asc9_tim_pirep,asc9_tim_rcvd
151         format(1x,a9,2x,a9)
            write(6,*)
            write(6,151)asc9_tim_pirep,asc9_tim_rcvd
        endif

        if(string(2:4) .eq. 'Lat')then
            read(11,201)rlat,rlon,ralt
201         format(2(f8.3,2x), f6.0,2i5)
            call latlon_to_rlapsgrid(rlat,rlon,lat,lon,ni,nj,ri,rj,istat
     1us)
!           if(istatus.ne.1)return
            ilaps = nint(ri)
            jlaps = nint(rj)
            write(6,201)rlat,rlon,ralt,ilaps,jlaps
        endif

        if(string(2:5) .eq. 'Clou')then

            if(                           ilaps .ge. ix_low   ! In bounds
     1                      .and. ilaps .le. ix_high  ! In bounds
     1                      .and. jlaps .ge. iy_low   ! In bounds
     1                      .and. jlaps .le. iy_high  ! In bounds
     1                                          )then ! In bounds

              n_cld_snd = n_cld_snd + 1  ! Create a new souding
              num_pirep = num_pirep + 1

              if(num_pirep .gt. max_pireps)then
                  write(6,*)' Too many pireps, check N_PIREP parameter'
                  istatus = 0
                  return
              endif

              l_good_pirep = .false.

              do i = 1,max_layers
                read(11,203,err=500)cbase_ft,ctop_ft,icover
203             format (12x,2f8.0,i5)
                write(6,*,err=500)' Got a cloud report',
     1                  nint(cbase_ft),nint(ctop_ft),icover

                if(icover .gt. 0)then ! Good Report (icover is in eighths)

                  if(cbase_ft .ge. 0. .and. ctop_ft .ge. 0.)then ! Fill Layer

                    cbase_m = cbase_ft / 3.281
                    ctop_m  = ctop_ft  / 3.281

                    cbuf_low = cbase_m - 500.
                    cbuf_high = ctop_m + 500.

                    cover_rpt = float(icover) / 8.0 ! 1.0

                    n_cloud = n_cloud + 1
                    write(6,*)' Good layer report:      '
     1          ,nint(cbase_m),nint(ctop_m),cover_rpt

                    do k=1,nk
                        cover = cover_rpt

!                       Fill in Cloud Layer
                        if(cld_hts(k) .ge. cbase_m .and.
     1                   cld_hts(k) .lt. ctop_m                  )then
                            call spread2(cld_snd,wt_snd,i_snd,j_snd,n_cl
     1d_snd
     1                          ,max_cld_snd,nk,ilaps,jlaps,k,cover,1.)
                            write(6,*)' Fill in k = ',k,cover,cld_hts(k)
                            l_good_pirep = .true.
                        endif

                        cover = .01

!                       Fill in clear buffer under cloud layer
                        if(cld_hts(k) .ge. cbuf_low .and.
     1                   cld_hts(k) .lt. cbase_m                  )then
                            call spread2(cld_snd,wt_snd,i_snd,j_snd,n_cl
     1d_snd
     1                          ,max_cld_snd,nk,ilaps,jlaps,k,cover,1.)
                            write(6,*)' Fill in k = ',k,cover,cld_hts(k)
                            l_good_pirep = .true.
                        endif

!                       Fill in clear buffer above cloud layer
                        if(cld_hts(k) .ge. ctop_m .and.
     1                   cld_hts(k) .lt. cbuf_high                )then
                            call spread2(cld_snd,wt_snd,i_snd,j_snd,n_cl
     1d_snd
     1                          ,max_cld_snd,nk,ilaps,jlaps,k,cover,1.)
                            write(6,*)' Fill in k = ',k,cover,cld_hts(k)
                            l_good_pirep = .true.
                        endif

                    enddo ! K cld_hts

                  elseif(cbase_ft .gt. 0)then
                    cbase_m = cbase_ft / 3.281
                    ctop_m  = cbase_m + 1000.

                    cbuf_low = cbase_m - 500.

                    cover_rpt = float(icover) / 8.0 ! 1.0

                    write(6,*)' Only a base reported:'
     1          ,nint(cbase_m),nint(ctop_m),cover_rpt

                    n_cloud = n_cloud + 1

                    do k=1,nk
                        cover = cover_rpt

!                       Fill in Cloud Layer
                        if(cld_hts(k) .ge. cbase_m .and.
     1                   cld_hts(k) .lt. ctop_m                  )then
                            call spread2(cld_snd,wt_snd,i_snd,j_snd,n_cl
     1d_snd
     1                          ,max_cld_snd,nk,ilaps,jlaps,k,cover,1.)
                            write(6,*)' Fill in k = ',k,cover,cld_hts(k)
                            l_good_pirep = .true.
                        endif

                        cover = .01

!                       Fill in clear buffer under cloud layer
                        if(cld_hts(k) .ge. cbuf_low .and.
     1                   cld_hts(k) .lt. cbase_m                  )then
                            call spread2(cld_snd,wt_snd,i_snd,j_snd,n_cl
     1d_snd
     1                          ,max_cld_snd,nk,ilaps,jlaps,k,cover,1.)
                            write(6,*)' Fill in k = ',k,cover,cld_hts(k)
                            l_good_pirep = .true.
                        endif

                    enddo ! K cld_hts

                  elseif(ctop_ft .gt. 0)then
                    ctop_m  = ctop_ft  / 3.281
                    cbase_m = ctop_m - 500.

                    cbuf_high = ctop_m + 500.

                    cover_rpt = float(icover) / 8.0 ! 1.0

                    write(6,*)' Only a top reported:'
     1          ,nint(cbase_m),nint(ctop_m),cover_rpt

                    n_cloud = n_cloud + 1

                    do k=1,nk
                        cover = cover_rpt

!                       Fill in Cloud Layer
                        if(cld_hts(k) .ge. cbase_m .and.
     1                   cld_hts(k) .lt. ctop_m                  )then
                            call spread2(cld_snd,wt_snd,i_snd,j_snd,n_cl
     1d_snd
     1                          ,max_cld_snd,nk,ilaps,jlaps,k,cover,1.)
                            write(6,*)' Fill in k = ',k,cover,cld_hts(k)
                            l_good_pirep = .true.
                        endif

                        cover = .01

!                       Fill in clear buffer above cloud layer
                        if(cld_hts(k) .ge. ctop_m .and.
     1                   cld_hts(k) .lt. cbuf_high                )then
                            call spread2(cld_snd,wt_snd,i_snd,j_snd,n_cl
     1d_snd
     1                          ,max_cld_snd,nk,ilaps,jlaps,k,cover,1.)
                            write(6,*)' Fill in k = ',k,cover,cld_hts(k)
                            l_good_pirep = .true.
                        endif

                    enddo ! K cld_hts

                  endif ! Base and/or Top reported

                endif ! Cover > 0

              enddo ! i cloud layer

              if(l_good_pirep)num_good_pirep = num_good_pirep + 1

            else ! Out of Bounds
!             Skip over cloud report
              do i = 1,max_layers
                read(11,203,err=500)cbase_ft,ctop_ft,icover
              enddo

              write(6,*)' Out of domain perimeter'

            endif ! In Bounds

        endif ! Cloud Report String

!       if(string(2:5) .eq. 'Sky ')then
!           read(11,204,err=500)isky_cover
!204         format (40x,i4)
!            write(6,*)' sky cover = ',isky_cover
!        endif

500     goto5


900     write(6,*)
        write(6,*)' Completed insertion of pireps'
        write(6,*)' num pireps/num good pireps/cloud layers = '
     1                  ,num_pirep,num_good_pirep,n_cloud
        close(11)

        istatus = 1
        return

998     write(6,*)' Error, could not find the pirep file'
        istatus = 1
        return

        end

