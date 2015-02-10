cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis
cdis
cdis   
cdis
        subroutine insert_pireps(i4time,cld_hts                          ! I
     1        ,default_clear_cover                                       ! I
     1        ,cld_snd,wt_snd,i_snd,j_snd,n_cld_snd,max_cld_snd          ! I/O
     1        ,lat,lon,ni,nj,nk,ix_low,ix_high,iy_low,iy_high,max_pireps ! I
     1        ,istatus)                                                  ! O

        real lat(ni,nj),lon(ni,nj)

        real cld_hts(nk)

!       Arrays for cloud soundings
        real cld_snd(max_cld_snd,nk)
        real wt_snd(max_cld_snd,nk)
        integer i_snd(max_cld_snd)
        integer j_snd(max_cld_snd)

        character*9  asc9_tim_pirep,asc9_tim_rcvd
        character*80 string

        character*31    ext

        logical l_good_pirep

!       1 Feb 1996     Steve Albers        Remove spread( calls

        max_layers = 3  ! Defined in PIN file and NetCDF RAW data

        lun = 11
        ext = 'pin'
        call open_lapsprd_file_read(lun,i4time,ext,istatus)
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
            call latlon_to_rlapsgrid(rlat,rlon,lat,lon,ni,nj,ri,rj
     1                              ,istatus)
!           if(istatus.ne.1)return
            ilaps = nint(ri)
            jlaps = nint(rj)
            write(6,201)rlat,rlon,ralt,ilaps,jlaps
        endif

        if(string(2:5) .eq. 'Clou')then

            if(                   ilaps .ge. ix_low   ! In bounds
     1                      .and. ilaps .le. ix_high  ! In bounds
     1                      .and. jlaps .ge. iy_low   ! In bounds
     1                      .and. jlaps .le. iy_high  ! In bounds
     1                                          )then ! In bounds

              n_cld_snd = n_cld_snd + 1  ! Create a new souding
              num_pirep = num_pirep + 1

              if(num_pirep .gt. max_pireps)then
                  write(6,*)' insert_pireps: Error, too many pireps,'
     1                     ,' check N_PIREP parameter'
                  istatus = 0
                  return
              endif

              l_good_pirep = .false.

!             PIN file input heights are feet MSL
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
     1                       ,nint(cbase_m),nint(ctop_m),cover_rpt

                    do k=1,nk
                        cover = cover_rpt

!                       Fill in Cloud Layer
                        if(cld_hts(k) .ge. cbase_m .and.
     1                     cld_hts(k) .lt. ctop_m                  )then
                            call spread2(cld_snd,wt_snd,i_snd,j_snd
     1                                  ,n_cld_snd,max_cld_snd
     1                                  ,nk,ilaps,jlaps,k,cover,1.)
                            write(6,*)' Fill in k = ',k,cover,cld_hts(k)
                            l_good_pirep = .true.
                        endif

                        cover = default_clear_cover

!                       Fill in clear buffer under cloud layer
                        if(cld_hts(k) .ge. cbuf_low .and.
     1                     cld_hts(k) .lt. cbase_m                 )then       
                            call spread2(cld_snd,wt_snd,i_snd,j_snd
     1                                  ,n_cld_snd,max_cld_snd
     1                                  ,nk,ilaps,jlaps,k,cover,1.)
                            write(6,*)' Fill in k = ',k,cover,cld_hts(k)
                            l_good_pirep = .true.
                        endif

!                       Fill in clear buffer above cloud layer
                        if(cld_hts(k) .ge. ctop_m .and.
     1                     cld_hts(k) .lt. cbuf_high               )then       
                            call spread2(cld_snd,wt_snd,i_snd,j_snd
     1                                  ,n_cld_snd,max_cld_snd
     1                                  ,nk,ilaps,jlaps,k,cover,1.)
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
     1                       ,nint(cbase_m),nint(ctop_m),cover_rpt

                    n_cloud = n_cloud + 1

                    do k=1,nk
                        cover = cover_rpt

!                       Fill in Cloud Layer
                        if(cld_hts(k) .ge. cbase_m .and.
     1                     cld_hts(k) .lt. ctop_m                  )then       
                            call spread2(cld_snd,wt_snd,i_snd,j_snd
     1                                  ,n_cld_snd,max_cld_snd
     1                                  ,nk,ilaps,jlaps,k,cover,1.)
                            write(6,*)' Fill in k = ',k,cover,cld_hts(k)
                            l_good_pirep = .true.
                        endif

                        cover = default_clear_cover

!                       Fill in clear buffer under cloud layer
                        if(cld_hts(k) .ge. cbuf_low .and.
     1                     cld_hts(k) .lt. cbase_m                 )then       
                            call spread2(cld_snd,wt_snd,i_snd,j_snd
     1                                  ,n_cld_snd,max_cld_snd
     1                                  ,nk,ilaps,jlaps,k,cover,1.)
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
     1                       ,nint(cbase_m),nint(ctop_m),cover_rpt

                    n_cloud = n_cloud + 1

                    do k=1,nk
                        cover = cover_rpt

!                       Fill in Cloud Layer
                        if(cld_hts(k) .ge. cbase_m .and.
     1                     cld_hts(k) .lt. ctop_m                  )then
                            call spread2(cld_snd,wt_snd,i_snd,j_snd
     1                                  ,n_cld_snd,max_cld_snd
     1                                  ,nk,ilaps,jlaps,k,cover,1.)
                            write(6,*)' Fill in k = ',k,cover,cld_hts(k)
                            l_good_pirep = .true.
                        endif

                        cover = default_clear_cover

!                       Fill in clear buffer above cloud layer
                        if(cld_hts(k) .ge. ctop_m .and.
     1                     cld_hts(k) .lt. cbuf_high               )then      
                            call spread2(cld_snd,wt_snd,i_snd,j_snd
     1                                  ,n_cld_snd,max_cld_snd
     1                                  ,nk,ilaps,jlaps,k,cover,1.)
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

998     write(6,*)' Warning, could not find the pirep file'
        istatus = 1
        return

        end

