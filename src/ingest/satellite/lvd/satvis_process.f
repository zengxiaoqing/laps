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
      subroutine process_vis_satellite(c_sat_id,
     &           c_sat_type,
     &           i4time,
     &           ni,nj,lat,lon,
     &           n_vis_lines,n_vis_elem,
     &           r_grid_ratio,
     &           image_vis,
     &           r_llij_lut_ri,
     &           r_llij_lut_rj,
     &           sublat_d,sublon_d,range_m,
     &           l_national,iskip_bilin,
     &           laps_vis_raw,laps_vis_norm,albedo,
     &           istatus)
c
c**************************************************************************
c
c       Routine to collect satellite data for the LAPS analyses.
c
c       Changes:
c        P.A. Stamus       12-01-92       Original (from 'get_vas_bt' routine).
c                     01-11-93       Write output in LAPS standard format.
c                     02-01-93       Add snooze call for vis data.
c                     09-16-93       Add Bands 3,4,5,12 data.
c       J.R. Smart    03-01-94	     Implement on the Sun.  Adapt to ISPAN
c                             grids. This required removing all references to
c                             GOES mdals/ground station satellite receiving.
c          "          10-26-94       modified for goes-8 data. No need to compute
c                                    radiance and btemp as this conversion is included
c                                    in icnt_lut.
c          "          11-28-95       Extract the visible processing part as a separate
c                                    subroutine.
c          "          10-21-96       Incorporated new normalize brightness procedure to
c				     account for vis data that has been pre-normalized (ie., gvar fsl-conus)
c          "          03-07-97      New normalize brightness routine that accounts for missing data.
c          "          03-12-97      Added c_sat_id to argument list. This directs the normalization procedure.
c
c       Notes:
c       This program processes vis satellite data from ISPAN database
c
c       Variables:
c       laps_vis         RA       O       Visible (raw)
c       vis_norm         RA       O       Visible (normalized)
c       albedo           RA       O       Albedo (0.0 -- 1.0)
c
c
c       Note: For details on the filtering and averaging methods, see the
c       VASDAT2 routine or talk to S. Albers.
c
c****************************************************************************
c
       implicit none
c
       integer ni,nj
c
c..... Grids to put the satellite data on.
c
       real*4 laps_vis_raw(ni,nj)
       real*4 laps_vis_norm(ni,nj)
       real*4 laps_vis_norm_natl(ni,nj)
       real*4 albedo(ni,nj)
       real*4 phase_angle_d(ni,nj)
       real*4 specular_ref_angle_d(ni,nj)
       real*4 rland_frac(ni,nj)
c 
c..... LAPS lat/lon files.
c
       real*4 lat(ni,nj)
       real*4 lon(ni,nj)
c
       integer n_vis_lines,n_vis_elem
       integer iskip_bilin
       integer i_dir
       integer n_missing_albedo
       integer len_dir
       integer isave,jsave
       integer ismax,jsmax
       integer ismin,jsmin

       real*4 r_llij_lut_ri(ni,nj)
       real*4 r_llij_lut_rj(ni,nj)
       real*4 r_grid_ratio
       real*4 image_vis(n_vis_elem,n_vis_lines)
       real*4 albedo_max,albedo_min
       real*4 rspacing_dum
       real*4 r_missing_data
c
       logical l_national

       integer istatus_a
       integer istatus_l
       integer istatus_m
       integer istatus_n
       integer istatus_v
       integer istatus_r
       integer istatus(3)
       integer i,j
c      integer imn,jmn,imx,jmx

       real*4 sublat_d,sublon_d,range_m
c      real*4 difference
c      real*4 diffsum
c      real*4 diffsum_abs
c      real*4 meandiff
c      real*4 meanabsdiff
c      real*4 maxdiff
c      real*4 mindiff
       real*4 ave,adev,sdev,var,skew,curt

       integer*4 icnt
c
       integer i4time,imax,jmax
c
       character*3   var_2d
       character*3   c_sat_type
       character*6   c_sat_id
       character*50  directory
       character*31  ext
       character*10  units_2d
       character*125 comment_2d
c
c Start
c ----------------------
       istatus(1) = 0
       istatus(2) = 0
       istatus(3) = 0

       imax = ni
       jmax = nj

       call get_r_missing_data(r_missing_data,istatus_r)
       if(istatus_r.ne.1)then
          write(6,*)'error getting r_missing_data'
          goto 999
       endif
c
c Call the visible satellite data remapping subroutine.
c
       write(6,*)
       write(6,*)' Processing visible satellite data.'
       write(6,*)' ----------------------------------'
       call satdat2laps_vis(
     &                  r_grid_ratio,
     &                  r_llij_lut_ri,
     &                  r_llij_lut_rj,
     &                  imax,jmax,
     &                  n_vis_lines,n_vis_elem, ! image_vis array dimensions
     &                  image_vis,
     &                  laps_vis_raw,
     &                  istatus_v)
       if(istatus_v .ne. 1) then
          write(*,921)istatus_v
          istatus(1) = -1
       end if
 921   format(1x,' +++ WARNING. Bad istatus_v = ',i3, 'from
     & SATDAT2LAPS_VIS+++')
c
       call check(laps_vis_raw,r_missing_data,istatus_v,imax,jmax)
       if(istatus_v .lt. 0) then
          write(6,*)'Bad visible counts detected'
          write(6,916) istatus_v
          istatus(1) = istatus_v
 916      format(' +++ WARNING. Visible status = ',i8)
       else
          write(6,*)'All laps vis data checked out ok'
       endif
c
c.....       Normalize the VIS data.
c
c For locally produced (ground station ispan look-alike) visible we
c want to reverse the already applied normalization (national scale
c normalization) and then "re-normalized" for the local domain.
c
c
       do j=1,jmax
       do i=1,imax
          laps_vis_norm(i,j)=laps_vis_raw(i,j)
       enddo 
       enddo

       if(c_sat_id.eq.'goes08'.and.c_sat_type.eq.'cdf')then    !this switch is for gvar conus /public data
c 	  						       that is a pre-normalized display-ready image.

          do j=1,jmax
          do i=1,imax
             laps_vis_norm_natl(i,j)=laps_vis_raw(i,j)
          enddo
          enddo

          write(6,*)'These vis data are already normalized'
          i_dir = -1
          l_national = .true.
          write(6,*)'Reverse the normalization'

          call normalize_brightness(i4time,lat,lon,laps_vis_norm,
     &imax,jmax,sublat_d,sublon_d,range_m,l_national,iskip_bilin,
     &r_missing_data,
     &6,i_dir,phase_angle_d,specular_ref_angle_d,istatus_n)
c
          if(istatus_n .ne. 1) then
             write(*,*)'+++WARNING+++ Bad status returned from
     &NORMALIZE LAPS VIS'
          else
             write(*,*)'Visible image successfully un-normalized'
          endif
c
c at this point the counts are goes7 look-a-like.
c So the standard for svs is to save the raw satellite vis counts.
c

          do j=1,jmax
          do i=1,imax
c
c save the unnormalized counts.
c Note: this raw vis data now looks like raw goes-7 counts and we must reverse the stretch.
c
             if(laps_vis_norm(i,j).ne.r_missing_data)then
               call stretch(0., 255., 0., 303.57, laps_vis_norm(i,j))  !reverse the stretch 
               laps_vis_raw(i,j)=laps_vis_norm(i,j)                    !OK, save these. This is svs!
c for goes8 - make it look like goes7
               call stretch(0., 303.57, 0., 255., laps_vis_norm(i,j))
c laps_vis_norm is now ready for local normalization
             endif

          enddo
          enddo
 
          write(6,*)
          i_dir = 1
          l_national = .false.
c
c compute differences between the two difference normalizations
c
c         diffsum = 0
c         diffsum_abs = 0
c         icnt = 0
c         maxdiff = -255
c         mindiff = 255
c         do j=1,jmax
c         do i=1,imax

c            if(laps_vis_norm_natl(i,j).ne.r_missing_data.and.
c    &laps_vis_norm(i,j).ne.r_missing_data)then
c            difference=laps_vis_norm_natl(i,j)-laps_vis_norm(i,j)
c            diffsum = difference + diffsum
c            diffsum_abs = abs(difference) + diffsum_abs
c            icnt = icnt + 1
c            if(difference.gt.maxdiff)then
c               maxdiff=difference
c               imx=i
c               jmx=j
c            endif
c            if(difference.lt.mindiff)then
c               mindiff=difference
c               imn=i
c               jmn=j
c            endif
c            endif
c
c         enddo
c         enddo

c         meandiff = diffsum/icnt
c         meanabsdiff = diffsum_abs/icnt

c         write(6,*)'========================================='
c         write(6,*)'Mean difference (natl-local): ',meandiff
c         write(6,*)'Mean absolute difference: ',meanabsdiff
c         write(6,*)'Max difference: ',maxdiff,' i/j ',imx,jmx
c         write(6,*)'Min difference: ',mindiff,' i/j ',imn,jmn
c         write(6,*)'========================================='
c
c================================================
c WFO Switch
c
       elseif(c_sat_type.eq.'wfo')then !Note: we need to somehow know whether this is goes-9
c 						       or goes-8!
          write(6,*)'These vis data will get normalized - WFO'
          if(c_sat_id.eq.'goes08')then
             write(6,*)'Stretch ',c_sat_id,' to goes7 look-a-like'
             do j=1,jmax
             do i=1,imax
               if(laps_vis_norm(i,j).ne.r_missing_data)then
c for goes8 - make it look like goes7
                 call stretch(0., 303.57, 0., 255., laps_vis_norm(i,j))
               endif
             enddo
             enddo
          elseif(c_sat_id.eq.'goes09')then
             write(6,*)'Stretch ',c_sat_id,' to goes7 look-a-like'
             do j=1,jmax
             do i=1,imax
               if(laps_vis_norm(i,j).ne.r_missing_data)then
c for goes9 - make it look like goes8
                  call stretch(0., 255., 0., 193.,laps_vis_norm(i,j))
c for goes8 - make it look like goes7
                  call stretch(0., 303.57, 0., 255., laps_vis_norm(i,j))
               endif
             enddo
             enddo
          endif
          i_dir = 1

c      else !ispan no longer an option. This would be GVAR ground station. NOT -> !assume this to be ispan data

c         write(6,*)'These vis data will get normalized'

c         do j=1,jmax
c         do i=1,imax
c            laps_vis_norm(i,j)=laps_vis_raw(i,j)
c for goes8 - make it look like goes7
c            call stretch(0., 303.57, 0., 255., laps_vis_norm(i,j))
c         enddo             !i
c         enddo

c         i_dir = 1
c=======================================
c GVAR Switch
c
       elseif(c_sat_type.eq.'gvr'.or.c_sat_type.eq.'gwc')then
          if(c_sat_id.eq.'goes08')then
             write(6,*)'GVAR GOES-8 Vis data'
             write(6,*)'Stretch ',c_sat_id,' to goes7 look-a-like'
             do j=1,jmax
             do i=1,imax
                if(laps_vis_norm(i,j).ne.r_missing_data)then
c for goes8 - make it look like goes7
                   call stretch(0.,303.57,0.,255.,laps_vis_norm(i,j))
                endif
             enddo
             enddo
          elseif(c_sat_id.eq.'goes09')then
             write(6,*)'GVAR GOES-9 Vis data'
             write(6,*)'Stretch ',c_sat_id,' to goes7 look-a-like'
             do j=1,jmax
             do i=1,imax
                if(laps_vis_norm(i,j).ne.r_missing_data)then
c for goes9 - make it look like goes8
                   call stretch(0., 255., 0.,193.,laps_vis_norm(i,j))
c for goes8 - make it look like goes7
                   call stretch(0.,303.57,0.,255.,laps_vis_norm(i,j))
                endif
             enddo
             enddo

          endif
          i_dir = 1

       endif 
c
c ready to normalize vis counts to local domain
c =============================================
       call normalize_brightness(i4time,lat,lon,
     &             laps_vis_norm,imax,jmax,
     &             sublat_d,sublon_d,
     &             range_m,l_national,iskip_bilin,
     &             r_missing_data,6,i_dir,phase_angle_d,
     &             specular_ref_angle_d,
     &             istatus_n)
       if(istatus_n .ne. 1) then
          write(*,*)'+++WARNING+++ Bad status returned from
     &NORMALIZE LAPS VIS'
          istatus(2) = istatus_n
       else
          write(*,*)'Visible image normalized for local domain'
          call check(laps_vis_norm
     &              ,r_missing_data,istatus_n,imax,jmax)
          istatus(2)=istatus_n  
       endif
c =============================================
       ext = 'nest7grid'

! Get the location of the static grid directory
       call get_directory(ext,directory,len_dir)

       var_2d='LDF'
       write(6,*)'read land fraction'

       call rd_laps_static (directory,ext,imax,jmax,1,var_2d,
     1units_2d,comment_2d,
     1rland_frac,rspacing_dum,istatus_l)
       if(istatus_l .ne. 1)then
           write(6,*)' Error reading LAPS static-ldf'
           return
       endif

       istatus_a = 1
       call vis_to_albedo(i4time,c_sat_id,
     &                    laps_vis_norm,
     &                    lat,lon,
     &                    imax,jmax,
     &                    r_missing_data,
     &                    phase_angle_d,
     &                    specular_ref_angle_d,
     &                    rland_frac,
     &                    albedo,
     &                    albedo_min,
     &                    albedo_max,
     &                    n_missing_albedo,
     &                    istatus_a)

       if(istatus_a .ne. 1)then
          write(*,*)'+++WARNING.+++ err status ',istatus_a,' from
     &               VIS_TO_ALBEDO'
          istatus(3) = istatus_a
       else
          write(*,*)'Albedo successfully computed'
       end if
       call check(albedo,r_missing_data,istatus_a,imax,jmax)
       if(istatus_a .lt. 1) then
          print*,' +++ WARNING. Visible status = ',istatus_a
          istatus(3) = istatus_a
       endif
c
       write(*,*)'Successfully remapped vis data'

       call moment(albedo,imax*jmax,
     &             ave,adev,sdev,var,skew,curt,
     &             istatus_m)
       if(istatus_m.ne.0)then
          print*,'Error returned from subroutine moment'
          goto 999
       endif

       icnt=0
       do j=1,jmax
       do i=1,imax
          if(albedo(i,j)-ave.gt.3.0*sdev)then
             icnt=icnt+1
             isave=i
             jsave=j
          endif
          if(albedo(i,j).eq.albedo_max)then
             ismax=i
             jsmax=j
          endif
          if(albedo(i,j).eq.albedo_min)then
             ismin=i
             jsmin=j
          endif
       enddo
       enddo

       i=isave
       j=jsave
       print*,'  ================='
       print*,'  Albedo statistics'
       print*,'  ================='
       if(icnt.gt.0.and.ismax.gt.0.and.jsmax.gt.0.and.
     &ismin.gt.0.and.jsmin.gt.0)then
          print*,'  N > 3 standard dev = ',icnt
          print*,'  Last i/j > 3 stand dev = ',i,j
          print*,'  laps_vis_norm(i,j)= ',laps_vis_norm(i,j)
          print*,'  Average albedo = ', ave
          print*,'  Standard dev = ',sdev
          print*,'  i/j/count of max albedo ',ismax,jsmax,
     &laps_vis_norm(ismax,jsmax)
          print*,'  i/j/count of min albedo ',ismin,jsmin,
     &laps_vis_norm(ismin,jsmin)
       else
          print*,'All albedo < 3 standard dev'
          print*,'No max/min output'
       endif

       goto 999
898    write(6,*)'Error getting r_msng_sat_flag'
c     
 999   return
       end
