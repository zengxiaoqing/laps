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
        subroutine storm_cent_rt(ni,nj,nz,i4time_latest,
     1  xlaps,ylaps,xradar,yradar,grid_ra_vel,grid_ra_ref,max2d_ref,
     1  max2d_refprv,
     1  lat,lon,stdlat,stdlon,umean,vmean,nstorm,
     1  istorm,jstorm,storm_u,storm_v,istatus)

c    Routine is called once per hour (actually at 20 after the hour).
c    Each run, the previous hour centroid file is read and the next
c    three (20 after, 40 after and on the hour) are run.  For example,
c    at 1823, it reads in 1700, and then runs 1720, 1740 and 1800.  At 1900,
c    it reads 1800, and then runs 1820, 1840 and 1900.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                               changes
c
c     7-11-91    mj   Deletes centroids within centroids.
c                mj   Changed threshold for cells within areas to 40 dbz
c                          from 35 dbz.

C INCLUDES AND DECLARATION OF VAR PASSED IN  ===========================

        integer*4 mxs
        parameter(mxs=75)

	integer*4 ni, nj, nz, istatus

        real*4 lat(ni,nj), lon(ni,nj),xlaps(ni,nj),ylaps(ni,nj),
     1       xradar(ni,nj),yradar(ni,nj)

        real*4 grid_ra_vel(ni,nj,nz),grid_ra_ref(ni,nj,nz),
     1       max2d_ref(ni,nj),max2d_refprv(ni,nj)

        integer*4 i4time_latest

        real*4 stdlat
        real*4 stdlon

        real*4 umean(ni,nj)     ! u_2d = previous, umean=current.
        real*4 vmean(ni,nj)     ! v"        "      v"      "

        real*4 storm_u(mxs),storm_v(mxs)
        integer*4 istorm(mxs),jstorm(mxs),nstorm

C LOCAL VARIABLE DECLARATIONS *ALL COMMENTED OUT FOR NOW* ==============

c       include 'segparms.for'

c       integer*4 maxseg,maxsegc
c       parameter(maxseg = 1000, maxsegc = 75)

c       integer*4       ss_normal,rtsys_bad_prod,rtsys_no_data,rtsys_abo
c    1rt_prod
c       parameter      (ss_normal        =1, ! success
c    1          rtsys_bad_prod   =2, ! inappropriate data, insufficient data
c    1          rtsys_no_data    =3, ! no data
c    1          rtsys_abort_prod =4) ! failed to make a prod

c       integer*4 j_status(20),iprod_number(20),i4time_array(20),
c    1            prod_array(10),prod_number

c       character atime*24, filename*9, yesno*1,no_data_type*13,
c    &            outfile*255,atimevil*24,btimevil*19,vilfilename*9,
c    &            jday*3,asc_time*24,jdate*9,timefile*50,
c    &            jdaytime*7,asc_day_hour*7,motion_srce*2,trash*80,
c    &            jdaytimes*7

c       character*1 ichar_new(mxs)

c       real*4 ref_max_seg(mxs,maxsegc),refbox_xmin(mxs),refbox_xmax(mxs
c    1),
c    1       refbox_ymin(mxs),refbox_ymax(mxs),refreg_max(mxs),
c    1       ref_tot_area(mxs),refreg_xcent(mxs),refreg_ycent(mxs),
c    1       rc_x_min(mxs),rc_x_max(mxs),rc_y_min(mxs),rc_y_max(mxs),
c    1       ref_seg_area(mxs,maxsegc),xstart(mxs,maxsegc),xend(mxs,maxs
c    1egc),
c    1       ystart(mxs,maxsegc),yend(mxs,maxsegc),ref_max(maxseg),
c    1       x_start_seg(maxseg),x_end_seg(maxseg),y_start_seg(maxseg),
c    1       y_end_seg(maxseg),refseg_area(mxs,maxsegc),sort_flag(mxs),
c    1       ref_max_y(maxseg),ref_max_x(maxseg),ref_max_segx(mxs,maxseg
c    1c),
c    1       ref_max_segy(mxs,maxsegc),refreg_maxx(mxs),refreg_maxy(mxs)
c    1,
c    1       refreg_xcent_ave(mxs),refreg_ycent_ave(mxs),storm_speed(mxs
c    1),
c    1       storm_dir(mxs),refeq_flag(mxs),refsegeq_flag(mxs),
c    1       ref_change(mxs),
c    1       dist_from(mxs)

c       real*4
c    1       xstart_prv(mxs,maxsegc),xend_prv(mxs,maxsegc),
c    1       ystart_prv(mxs,maxsegc),yend_prv(mxs,maxsegc),
c    1              ref_change_prv(mxs),
c    1              refreg_xcent_prv(mxs),
c    1              refreg_ycent_prv(mxs),
c    1              refbox_xmin_prv(mxs),
c    1              refbox_xmax_prv(mxs),
c    1              refbox_ymin_prv(mxs),
c    1              refbox_ymax_prv(mxs),
c    1              refreg_max_prv(mxs),
c    1              refreg_maxx_prv(mxs),
c    1              refreg_maxy_prv(mxs),
c    1              refreg_xcent_ave_prv(mxs),
c    1              refreg_ycent_ave_prv(mxs),
c    1              storm_speed_prv(mxs),
c    1              storm_dir_prv(mxs)


c       integer*4  ir_cluster_name(maxseg),ir_seg_number(mxs),
c    1           i_merge_flag(mxs),i_valid_rclust(mxs),num_valid_rclust,
c    1           i_seg_num(mxs),i_seg_num_prv(mxs),jgrid(maxseg),
c    1           igrid_start(maxseg),istorm_id_prv(mxs),
c    1           igrid_end(maxseg),j_dist,notalone_flag(maxseg),
c    1           igrid_cent(mxs),jgrid_cent(mxs),istorm_corr(mxs),
c    1           igrid_cent_prv(mxs),jgrid_cent_prv(mxs),
c    1           mw_motion_flag(mxs),i1grid_flag(mxs),prv_max_ipos,
c    1           prv_max_jpos,new_flag(mxs)

c       logical available(maxseg),l_low_fill,l_high_fill

c.....  Stuff for LAPS Surface file (standard).
c
c       real*4 sigecho_flag

c       integer*4 i4time,imax,jmax,i4time_radar,
c    1          i4time_wind_prv

c       character time*35,cm_filename*34,ca_filespec*255

c       character*9 asc9_tim_cent,asc9_htim_wind,asc9_prev,asc9_tim
c       character*125 var_2d
c       character*50  directory
c       character*10  units_2d
c       character*125 comment_2d
c       character*24 asc_time_run
c       character*4  asc_hour
c       character*255 c255_radar_filename

c       character*255 ca_rfilename
c       character*255 ca_cfilename
c       character*255 ca_wfilename

c       character*31  rc_ext
c       character*31  dc_ext
c       character*31  wm_ext

c       data rc_ext /'vrc'/
c       data dc_ext /'rdc'/
c       data wm_ext /'lwm'/
c
c  dimensions for derived laps surface variables
c
c       real*4 ref_thresh,thresh_ovlp,thresh_range
c       data ref_thresh/30.0/,thresh_ovlp/-1.0E-3/,thresh_range/40.0/

c       integer*4 thresh_run
c       data thresh_run/1/

c       real*4 re
c       data re /6371.2293/         !radius of earth

c       real*4 dtr
c       data dtr /0.01745329/       !degrees to radians

c       real*4 r(4),

c       integer*4 i4time_tol,i4time_hour,
c    1          i4time_get,sys$trnlog,lenfil,ISTATUS,
c    1          itimes

C BEGIN SUBROUTINE======================================================

        nstorm = 0
        istatus = 1

        write(6,91)
91      format(' Finished with dummy storm centroid/tracking',/,
     1       ' =====================================')
        write(6,*)
        write(6,*)

        return
        end
