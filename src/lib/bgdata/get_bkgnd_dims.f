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
      subroutine get_bkgnd_dims(bgmodel,bgfname,nx_bg,ny_bg,nz_bg,nrecs,
     +                         cmodel,nf_fid)
      implicit none
      integer bgmodel, nx_bg,ny_bg,nz_bg, nrecs, istatus, nf_fid
      character*(*) bgfname, cmodel
 
      nf_fid=1
      nrecs=1
      if(bgmodel.eq.0) then
         call get_grid_dim_xy(nx_bg,ny_bg,istatus)
         call get_laps_dimensions(nz_bg,istatus)
      else if(bgmodel.eq.1) then
         nx_bg = 81
         ny_bg = 62
         nz_bg = 25        
         cmodel = 'RUC60_NATIVE'   
      else if(bgmodel.eq.2) then
         call get_eta48_dims(bgfname,nx_bg,ny_bg,nz_bg,nrecs,nf_fid)
         cmodel = 'ETA48_CONUS'        
      else if(bgmodel.eq.3) then
         nx_bg = 144
         ny_bg = 73
         nz_bg = 16        
         cmodel = 'NOGAPS (2.5)'            
      else if(bgmodel.eq.4) then
         nx_bg = 93
         ny_bg = 65
         nz_bg = 19        
         cmodel = 'SBN CONUS211'            
      else if(bgmodel.eq.5) then
         nx_bg = 151
         ny_bg = 113
         nz_bg = 40       
         cmodel = 'RUC40_NATIVE'            
      else if(bgmodel.eq.6) then
         nx_bg = 360
         ny_bg = 181
         nz_bg = 26
         cmodel = 'AVN_LL_GRIB'
      else if(bgmodel.eq.7) then
         nx_bg = 185
         ny_bg = 129
         nz_bg = 42
         cmodel = 'ETA48_GRIB'
      else if(bgmodel.eq.8) then
         nx_bg = 360
         ny_bg = 181
         nz_bg = 16        
         cmodel = 'NOGAPS (1.0)'            
      else if(bgmodel.eq.9) then
cc         call get_file_names(bgpath,bg_files,names,max_files,istat)
cc         if (istat .ne. 1) print *,'Error in get_file_names.'
         call get_conus_dims(bgfname,nx_bg,ny_bg,nz_bg)
         cmodel='NWS_CONUS'
      endif
         
      return
      end

