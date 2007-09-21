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
        subroutine steer_grid(i4time_latest,ni,nj,nk
     1    ,xlaps,ylaps,xradar,yradar,grid_ra_vel,grid_ra_ref,max2d_ref
     1    ,max2d_refprv
     1                  ,lat,lon,standard_latdum,standard_londum
     1                  ,iiilut
     1                  ,umean,vmean,steer_u,steer_v,
     1                                                    istatus)

        integer max_storms
        parameter(max_storms=75)

!       Dummy arrays
        real xlaps(ni,nj),ylaps(ni,nj),xradar(ni,nj),yradar(ni,nj)
        real grid_ra_vel(ni,nj,nk),grid_ra_ref(ni,nj,nk)
        real max2d_ref(ni,nj),max2d_refprv(ni,nj)

        real lat(ni,nj),lon(ni,nj)

        real iiilut(-ni:ni,-nj:nj)

        real steer_u(ni,nj),steer_v(ni,nj)
     1       ,umean(ni,nj),vmean(ni,nj)    ! Unmodified mean winds WRT True N
     1       ,storm_u(max_storms),storm_v(max_storms),wt_ob(max_storms)

        integer i4time_latest,istatus,istorm(max_storms),jstorm(max_st
     1orms)

        n_storms = 0
        iswitch_centroids = 1

!       open(86,file='wind.parms',status='old'
!       1                                                               ,err=1)
!       read(86,*,err=1)iswitch_centroids
!       close(86)

!1      write(6,*)' iswitch_centroids = ',iswitch_centroids

        if(iswitch_centroids .eq. 1)then

            call storm_cent_rt(ni,nj,nk,i4time_latest,
     1          xlaps,ylaps,xradar,yradar,grid_ra_vel,grid_ra_ref,max2d_
     1ref,
     1          max2d_refprv,
     1          lat,lon,standard_latdum,standard_londum,
     1    umean,vmean,n_storms,
     1          istorm,jstorm,
     1          storm_u,storm_v,         ! Storm Motions WRT True or Grid N?
     1                                  istatus)

        endif

        write(6,*)
        write(6,*)' Returned from storm_cent_rt'
        write(6,*)'   #             Storm  U/V    Mean U/V'
     1                  ,'  Stm Dir/Spd  Mean Dir/Spd'

        do i = 1,n_storms
            umean_at_storm = umean(istorm(i),jstorm(i))
            vmean_at_storm = vmean(istorm(i),jstorm(i))

            call uv_to_disp(storm_u(i),storm_v(i),dir_storm,spd_storm)

            call uv_to_disp(umean_at_storm,vmean_at_storm
     1                  ,dirmean_at_storm,spdmean_at_storm)

            write(6,101)i,istorm(i),jstorm(i)
     1          ,storm_u(i),storm_v(i)
     1          ,umean_at_storm,vmean_at_storm
     1          ,nint(dirmean_at_storm),spdmean_at_storm
     1          ,nint(dir_storm),spd_storm
101         format(1x,i4,i5,i5,2f6.1,2x,2f6.1,2(2x,i3,'/',f6.1))

        enddo ! i

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    Call OA routine to derive steering wind grid using umean,vmean,storm_u,
c    storm_v

        call analyze_storm_motion(ni,nj,umean,vmean        ! Input
     1    ,n_storms,istorm,jstorm,storm_u,storm_v          ! Input
     1  ,iiilut                                    ! Input (Local)
     1  ,wt_ob                                     ! Input (Local)
     1        ,max_storms                                ! Input
     1  ,steer_u,steer_v                        )  ! Output

        return
        end


        subroutine analyze_storm_motion(ni,nj,umean,vmean ! Input
     1    ,n_storms,istorm,jstorm,storm_u,storm_v         ! Input
     1  ,iiilut                                   ! Input (Local)
     1  ,wt_ob                                    ! Input (Local)
     1        ,max_storms                               ! Input
     1  ,steer_u,steer_v                       )  ! Output

        real iiilut(-ni:ni,-nj:nj)

        real steer_u(ni,nj),steer_v(ni,nj),umean(ni,nj),vmean(ni,nj),
     1       storm_u(max_storms),storm_v(max_storms)

c       integer i4time_latest,istatus,istorm(max_storms),jstorm(max_st
c    1orms)
        integer istorm(max_storms),jstorm(max_storms)

        real wt_ob(max_storms)

        if(n_storms .gt. 0)then
            do i = 1,n_storms
                wt_ob(i) = 1.0
            enddo ! i

            call barnes_r5th(ni,nj,n_storms,max_storms
     1          ,istorm,jstorm,storm_u,wt_ob,umean,iiilut,steer_u)

            call barnes_r5th(ni,nj,n_storms,max_storms
     1          ,istorm,jstorm,storm_v,wt_ob,vmean,iiilut,steer_v)

        else ! Just copy the mean wind field into the output grids
            write(6,*)' No storms, returning mean wind field'
            do j = 1,nj
            do i = 1,ni
                steer_u(i,j) = umean(i,j)
                steer_v(i,j) = vmean(i,j)
            enddo ! i
            enddo ! j

        endif

        return
        end


      subroutine barnes_r5th(ni,nj,ncnt,max_obs
     1                  ,iob,job,obs,wt_ob,background_field
     1                                  ,iiilut,anal)

      real exponent_distance_wt,background_weight
      parameter (exponent_distance_wt = 5.0)
      parameter (background_weight = 1.0)

      integer  n_fnorm
      parameter (n_fnorm = 10000)

      dimension  iob(max_obs),job(max_obs),obs(max_obs),wt_ob(max_obs)
     1  ,fnorm(n_fnorm)
     1        ,anal(ni,nj)

      real background_field(ni,nj)

      real iiilut(-ni:ni,-nj:nj)

      write(6,*)' Barnes_r5th called'

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Error in barnes_r5th: STOP'
          stop
      endif

      do iii = 1,n_fnorm
        fnorm(iii) = (100./float(iii)) ** (exponent_distance_wt / 2.0)
      enddo

      write(6,*)' Ncnt = ',ncnt

C     THIS IS SEt FOR INCREMENTS OF .01
      spcng = 10.
      radm2=1.0/spcng**2
      write(6,*)' radm2*100 = ',radm2*100.

!     Create a lookup table for (iii)
      do i = -ni,ni
      do j = -nj,nj
          rsq=i*i+j*j
          iii=radm2*100.*rsq+1.
          if(iii .gt. n_fnorm)then
              iiilut(i,j) = 0.
          else
              iiilut(i,j) = fnorm(iii)
          endif
!         write(6,1111)i,j,iii,fnorm(iii),iiilut(i,j)
1111      format(2i4,i6,2f10.5)
      enddo
      enddo

!     Note that it is OK if iii exceeds n_fnorm because of the test above
      write(6,*)' Highest value of iii (compared to n_fnorm) = ',iii,n_f
     1norm

      do j=1,nj
      do i=1,ni
          sum=  background_weight * background_field(i,j)
          sumwt=background_weight

          do n=1,ncnt
              ii=iob(n)
              jj=job(n)
              weight = iiilut(i-ii,j-jj) * wt_ob(n) ! Obs are being weighted
              sum=weight*obs(n)+sum
              sumwt=sumwt+weight
          enddo

          if (sumwt.eq.0.)then
              anal(i,j) = r_missing_data
          ELSE
              anal(i,j)=sum/sumwt
          endif

      enddo ! i
      enddo ! j
      return
      end
