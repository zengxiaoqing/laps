      program dprep_driver
c
      implicit none

      include 'bgdata.inc'
c
c eventually put these parameters into bgdata.inc
      integer mxvars,mxlvls
      parameter(mxvars=10,mxlvls=100)

      integer NX,NY,NZ                !Background grid dimensions

      integer nzbg_ht,nzbg_uv,nzbg_sh,nzbg_ww
c
c cmodel is really only 12 chars but the SBN netcdf carrys 132
c
      integer i,l, filesfound, len
      character*256 nl_file
      integer bgfcnt, max_files, istat, istatus
      parameter (max_files=500)
      character*256 bgfnames(max_files)
      integer  maxtimes
      parameter (maxtimes=24)
      integer inittimes(maxtimes)
      integer bndytimes(maxtimes)
      integer fcstlengths(maxtimes)
      integer bgmodel
      integer oldest_forecast, max_forecast_delta
      integer accepted_files
      logical use_analysis, outdir_defined, bgpath_defined
      character*9 a9
      character*13 cvt_i4time_wfo_fname13
      integer i4time_latest
      integer i4time_now, lendr, itime, ntime
      integer ntbg, idum
      character*256 cfilespec,cfname
      character*256 bgpath, topofname, outdir
      character*256 initbgpaths(maxbgmodels)
      character*256 bndybgpaths(maxbgmodels)
      integer       initbgmodels(maxbgmodels)
      integer       bndybgmodels(maxbgmodels)

c these are dummy variables used in get_bkgd_mdl_info (returns nx,ny,nz)
      character*2 gproj
      real        dlat,dlon,Lat1,Lat2,Lon0,sw(2),ne(2)
      real        centrallat,centrallon
      real        dx,dy
                 
      character*132 initcmodels(maxbgmodels)
      character*132 bndycmodels(maxbgmodels)
      character*132 cmodel
      character*256 reject_files(max_files)
      character*200 cdirstatic
      integer reject_cnt
      data reject_cnt/0/
c     data max_forecast_delta/18/
      data use_analysis/.true./
      data outdir_defined/.false./
      data bgpath_defined/.false./

      namelist /dprep_nl/initbgpaths,initbgmodels,initcmodels
     +                  ,inittimes
     +                  ,bndybgpaths,bndybgmodels,bndycmodels
     +                  ,bndytimes,  fcstlengths
     +                  ,max_forecast_delta,oldest_forecast
     +                  ,use_analysis, outdir

c

c_______________________________________________________________________________
c
c *** Initialize tables.
c
      include 'grid_fname.cmn' 
      call es_ini
c
c *** Read namelist for data sources
c
      call get_config(istatus)
      call get_directory('static',cdirstatic,lendr)
      nl_file = cdirstatic(1:lendr)//'/dprep.nl'
      open(21,FILE=nl_file,status='old',err=900)
      read(21,dprep_nl,err=901)
      close(21)
      call s_len(outdir,len)
      if(len.gt.0) inquire(FILE=outdir,EXIST=outdir_defined)
      if(.not.outdir_defined) then
         call s_len(generic_data_root,lendr)
         outdir=generic_data_root(1:lendr)//'/lapsprd/dprep'
      endif


      accepted_files = 0
      i = 1
      call get_systime(i4time_now,a9,istat)

      read(a9(6:7),'(i2)') ntime
      itime=1

      do while(inittimes(itime).ne.ntime .and. itime.lt.maxtimes)
         itime=itime+1
      enddo   
      
c this section gets model initialization data

      if(inittimes(itime).eq.ntime) then 
         print*, 'Dprep at time ',inittimes(itime),'for forecast of',
     +        fcstlengths(itime), ' hours'

         do while(accepted_files.le.0 .and. i.le.maxbgmodels)

            bgpath_defined=.false.
            bgmodel = initbgmodels(i)
            bgpath = initbgpaths(i)
            cmodel = initcmodels(i)
            call s_len(bgpath,len)
            if(len.gt.0) inquire(FILE=bgpath,EXIST=bgpath_defined)
            if(bgmodel.eq.0.and..not.bgpath_defined) then
               bgpath=generic_data_root(1:lendr)
            elseif(.not.bgpath_defined) then
               print*,'Could not find directory ',bgpath
               stop
            endif

            if(bgmodel .eq. 4)then
               cfilespec=bgpath(1:len)//'/*'
               call get_latest_file_time(cfilespec,i4time_latest)
               cfname = cvt_i4time_wfo_fname13(i4time_latest)

c              get_sbn_dims(cdfname,cmodel,mxvars,mxlvls
c    +,nvars,nxbg,nybg,nzbg_ht,nzbg_sh,nzbg_uv,nzbg_ww
c    +,n_valtimes,istatus)

c              call get_sbn_dims(bgpath,cfname,idum,idum,idum,ntbg)
            endif

            call get_acceptable_files(i4time_now,bgpath,bgmodel
     +          ,bgfnames,max_files,oldest_forecast,0,use_analysis
     +          ,bgfcnt,accepted_files,-1,cmodel,nx,ny,nz
     +          ,reject_files,reject_cnt)

            if(accepted_files .gt. 0)then
               cfname=bgpath(1:len)//bgfnames(accepted_files)


               call get_bkgd_mdl_info(bgmodel,cmodel,cfname
     &,mxvars,mxlvls,nx,ny,nzbg_ht,nzbg_sh,nzbg_uv,nzbg_ww
     &,gproj,dlat,dlon,centrallat,centrallon,dx,dy
     &,Lat1,Lat2,Lon0,sw,ne,istatus)

               if(istatus.ne.1)then
                  print*,'Error getting background model information'
                 stop
               endif

c     print*,nx,ny,nz,i4time_now,bgmodel,bgpath,
c     +         bgfnames,bgfcnt,outdir,filesfound,bgfcnt

               call dprep_sub(nx,ny,nzbg_ht,nzbg_uv,nzbg_sh,nzbg_ww
     +,mxlvls,i4time_now,bgmodel,bgpath,bgfnames,bgfcnt,outdir
     +,accepted_files,cmodel,gproj,dx,dy,istatus)

               if(istatus.ne.0)then
                  print*,'No dprep output due to error'
                  i=maxbgmodels+1
               endif

            else
               i=i+1
            endif

         enddo


         if(accepted_files.le.0) then
            print*,'ERROR: Could not find suitable file to initialize'
c            stop
         endif

      endif
      
      itime=1
      do while(bndytimes(itime).ne.ntime .and. itime.lt.maxtimes)
         itime=itime+1
      enddo


      if(bndytimes(itime).eq.ntime) then

         accepted_files=0
         i=1
         do while(accepted_files.le.0 .and. i.le.maxbgmodels)
            bgmodel = bndybgmodels(i)
            bgpath = bndybgpaths(i)
            cmodel = bndycmodels(i)
            call s_len(bgpath,len)
            if(bgpath(len:len).ne.'/')then
               len=len+1
               bgpath(len:len)='/'
            endif

            if(bgmodel .eq. 4)then
               cfilespec=bgpath(1:len)//'*'
               call get_latest_file_time(cfilespec,i4time_latest)
               if(i4time_latest .ne. 0)then
                  cfname = cvt_i4time_wfo_fname13(i4time_latest)

c                 get_sbn_dims(cdfname,cmodel,mxvars,mxlvls
c    +,nvars,nxbg,nybg,nzbg_ht,nzbg_sh,nzbg_uv,nzbg_ww
c    +,n_valtimes,istatus)

c                 call get_sbn_dims(bgpath,cfname,idum,idum,idum,ntbg)
               else
                  i=maxbgmodels  !6-30-00: Smart: no files (i4time_latest = 0) terminates the loop.
               endif
            endif

            call get_acceptable_files(i4time_now,bgpath,bgmodel
     +          ,bgfnames,max_files,oldest_forecast,max_forecast_delta
     +          ,use_analysis,bgfcnt,accepted_files,fcstlengths(itime)
     +          ,cmodel,nx,ny,nz,reject_files,reject_cnt)

            if(accepted_files .gt. 0)then
               cfname=bgpath(1:len)//bgfnames(accepted_files)

               call get_bkgd_mdl_info(bgmodel,cmodel,cfname
     &,mxvars,mxlvls,nx,ny,nzbg_ht,nzbg_sh,nzbg_uv,nzbg_ww
     &,gproj,dlat,dlon,centrallat,centrallon,dx,dy
     &,Lat1,Lat2,Lon0,sw,ne,istatus)

               print *, nx,ny,nzbg_ht,nzbg_uv,nzbg_sh,nzbg_ww

               call dprep_sub(nx,ny,nzbg_ht,nzbg_uv,nzbg_sh,nzbg_ww
     +,mxlvls,i4time_now,bgmodel,bgpath,bgfnames,accepted_files
     +,outdir,filesfound,cmodel,gproj,dx,dy,istatus)

               if(istatus.ne.0)then
                  print*,'No dprep output due to error'
                  i=maxbgmodels+1
               endif

            else
               i=i+1
            endif
         
         enddo
      endif

      print*, 'Normal completion of dprep'
      if(bgfcnt.le.0 .or. accepted_files.le.0)then
        print*, 'No dprep request for this time'
      endif

      stop

 900  print *, 'ERROR: Could not open file ',nl_file,' to read'
      stop
 901  print *,'ERROR: Could not read namelist'
      write(*,dprep_nl)
      stop
c
      end
c
      subroutine dprep_sub(nx_bg,ny_bg,nzbg_ht,nzbg_uv
     +,nzbg_sh,nzbg_ww,mxlvls,i4time,bgmodel,bgpath
     +,bgfnames,bgfcnt,outdir,filesfound,cmodel,gproj,dx,dy
     +,istatus_dprep)

      implicit none
      include 'netcdf.inc'

      integer nx_bg, ny_bg
      integer nzbg_ht,nzbg_uv,nzbg_sh,nzbg_ww
      integer mxlvls
      integer i4time, i, j, bgtime
c     integer bgmodel, istatus, bgfcnt, k, len, filesfound
      integer bgmodel, istatus, bgfcnt, k, filesfound
      integer istatus_dprep

      character*(*) bgpath,bgfnames(bgfcnt),outdir,cmodel
      character*256 fullname
      character*256 filename
      character*2 gproj
      character*4 ext, af
      character*9 fname
      character*13 fname13,fname9_to_wfo_fname13

c     real La1, La2, Lo1, Lo2, ht( NX,  NY,  NZ), 

c
c *** Background model grid data.
c
c
c *** sfc background arrays.
c
      real, allocatable  :: prbg_sfc(:,:)
      real, allocatable  :: uwbg_sfc(:,:)
      real, allocatable  :: vwbg_sfc(:,:)
      real, allocatable  :: shbg_sfc(:,:)
      real, allocatable  :: tpbg_sfc(:,:)
      real, allocatable  :: htbg_sfc(:,:)
      real, allocatable  :: mslpbg(:,:)
c
c *** 3D background arrays.
c
      real, allocatable  :: prbght(:,:,:)
      real, allocatable  :: prbgsh(:,:,:)
      real, allocatable  :: prbguv(:,:,:)
      real, allocatable  :: prbgww(:,:,:)
      real, allocatable  :: htbg(:,:,:)
      real, allocatable  :: tpbg(:,:,:)
      real, allocatable  :: exbg(:,:,:)
      real, allocatable  :: shbg(:,:,:)
      real, allocatable  :: uwbg(:,:,:)
      real, allocatable  :: vwbg(:,:,:)
      real, allocatable  :: wwbg(:,:,:)

c     real ht( NX,  NY,  NZ), 
c    +   ht_sfc( NX,  NY), p_sfc( NX,  NY), 
c    +   rh( NX,  NY,  NZ), rh_sfc( NX,  NY), 
c    +   th( NX,  NY,  NZ), th_sfc( NX,  NY), 
c    +   uw( NX,  NY,  NZ), uw_sfc( NX,  NY), 
c    +   vw( NX,  NY,  NZ), vw_sfc( NX,  NY),
c    +    p( NX,  NY,  NZ), ww(NX, NY, NZ),
c    +   mslp(NX, NY)

      real sw(2), ne(2)
      real lat1, lat2, lon0
      real dx,dy

      interface
c
         subroutine read_bgdata(mxlvls,nx_bg,ny_bg
     +,nzbg_ht,nzbg_sh,nzbg_uv,nzbg_ww,ctype
     +,bgpath,fname_bg,af_bg,fullname,cmodel,bgmodel
     +,prbght,prbgsh,prbguv,prbgww
     +,htbg,tpbg,uwbg,vwbg,shbg,wwbg
     +,htbg_sfc,prbg_sfc,shbg_sfc,tpbg_sfc
     +,uwbg_sfc,vwbg_sfc,mslpbg,istatus)
c
         real  :: prbg_sfc(:,:)
         real  :: uwbg_sfc(:,:)
         real  :: vwbg_sfc(:,:)
         real  :: shbg_sfc(:,:)
         real  :: tpbg_sfc(:,:)
         real  :: htbg_sfc(:,:)
         real  ::   mslpbg(:,:)
c
         real  :: prbght(:,:,:)
         real  :: prbgsh(:,:,:)
         real  :: prbguv(:,:,:)
         real  :: prbgww(:,:,:)
         real  :: htbg(:,:,:)
         real  :: tpbg(:,:,:)
         real  :: shbg(:,:,:)
         real  :: uwbg(:,:,:)
         real  :: vwbg(:,:,:)
         real  :: wwbg(:,:,:)

         character*4   af_bg
         character*5   ctype
         character*132 cmodel
         character*256 bgpath
         character*256 fname_bg
         character*256 fullname
         integer       bgmodel
         integer       mxlvls
         integer       nx_bg
         integer       ny_bg
         integer       nzbg_ht
         integer       nzbg_sh
         integer       nzbg_uv
         integer       nzbg_ww
         integer       istatus
         end subroutine

         subroutine dprep_output(outdir,fname,cmodel,gproj
     +,sw,ne,nx,ny,nzht,nzuv,nzrh,nzww,prht,pruv,prrh,prww
     +,ht,ex,th,uw,vw,ww,rh,ht_sfc,p_sfc,rh_sfc,th_sfc
     +,uw_sfc,vw_sfc,mslp,lat1,lat2,lon0,dx,dy,istatus)

         real  ::  p_sfc(:,:)
         real  :: uw_sfc(:,:)
         real  :: vw_sfc(:,:)
         real  :: rh_sfc(:,:)
         real  :: th_sfc(:,:)
         real  :: ht_sfc(:,:)
         real  ::   mslp(:,:)
c
         real  :: prht(:,:,:)
         real  :: prrh(:,:,:)
         real  :: pruv(:,:,:)
         real  :: prww(:,:,:)
         real  :: ht(:,:,:)
         real  :: ex(:,:,:)
         real  :: th(:,:,:)
         real  :: rh(:,:,:)
         real  :: uw(:,:,:)
         real  :: vw(:,:,:)
         real  :: ww(:,:,:)

         real  :: sw(2),ne(2)
         real  :: lat1,lat2,lon0
         real  :: dx,dy

         character*132 cmodel
         character*256 outdir
         character*256 bgpath
         character*256 fname
         character*2   gproj
         integer       bgmodel
         integer       mxlvls
         integer       nx
         integer       ny
         integer       nzht
         integer       nzrh
         integer       nzuv
         integer       nzww
         integer       istatus
         end subroutine

      end interface

      allocate (htbg(nx_bg,ny_bg,nzbg_ht))
      allocate (tpbg(nx_bg,ny_bg,nzbg_sh))
      allocate (exbg(nx_bg,ny_bg,nzbg_sh))
      allocate (shbg(nx_bg,ny_bg,nzbg_sh))
      allocate (uwbg(nx_bg,ny_bg,nzbg_uv))
      allocate (vwbg(nx_bg,ny_bg,nzbg_uv))
      allocate (wwbg(nx_bg,ny_bg,nzbg_ww))
      allocate (prbght(nx_bg,ny_bg,nzbg_sh))
      allocate (prbguv(nx_bg,ny_bg,nzbg_uv))
      allocate (prbgsh(nx_bg,ny_bg,nzbg_sh))
      allocate (prbgww(nx_bg,ny_bg,nzbg_ww))

      allocate (htbg_sfc(nx_bg,ny_bg))
      allocate (prbg_sfc(nx_bg,ny_bg))
      allocate (shbg_sfc(nx_bg,ny_bg))
      allocate (uwbg_sfc(nx_bg,ny_bg))
      allocate (vwbg_sfc(nx_bg,ny_bg))
      allocate (tpbg_sfc(nx_bg,ny_bg))
      allocate (mslpbg(nx_bg,ny_bg))

      istatus_dprep=1
      filesfound = 0
      
      
      do k=1,bgfcnt

         call s_len(bgfnames(k),j)
         j=j-13
         fname = bgfnames(k)(j+1:j+9)
         af = bgfnames(k)(j+10:j+13)
         call i4time_fname_lp(fname,bgtime,istatus)

         if(bgmodel.eq.4)then
            bgfnames(k)=fname9_to_wfo_fname13(bgfnames(k))
         endif

         call s_len(bgpath,i)
         fullname = bgpath(1:i)//'/'//bgfnames(k)

         call read_bgdata(mxlvls,nx_bg,ny_bg,
     +	   nzbg_ht,nzbg_sh,nzbg_uv,nzbg_ww,'dprep'
     +    ,bgpath,bgfnames(k),af,fullname,cmodel,bgmodel
     +    ,prbght,prbgsh,prbguv,prbgww
     +    ,htbg, tpbg,uwbg,vwbg,shbg,wwbg
     +    ,htbg_sfc,prbg_sfc,shbg_sfc,tpbg_sfc
     +    ,uwbg_sfc,vwbg_sfc,mslpbg,istatus_dprep)


c        call read_bgdata(nx,ny,nz,bgpath,bgfnames(k),af
c    + ,fullname,cmodel,bgmodel
c    + ,ht,p,th,uw,vw,rh,ww
c    + ,ht_sfc,p_sfc,rh_sfc,th_sfc,uw_sfc,vw_sfc,mslp
c    + ,gproj,lon0,lat1,lat2,istatus_dprep)

         if (istatus_dprep .ne. 0) then

            if (bgmodel .gt. 1 .and. bgmodel .le. 3) then
               fname13=fname(1:9)//af
            elseif (bgmodel .eq. 4) then
               fname13=fname9_to_wfo_fname13(fname)
            endif
            print *,'Error reading background model data for: ',
     .         bgpath(1:i)//'/'//fname13
            print *,'Process aborted for this file.'
            return

         endif


         if(bgmodel.eq.0)then
            ext = '.LPS'
            gproj='PS'
            af='0000'

            call get_laps_data(bgpath
     .,i4time,nx_bg,ny_bg,nzbg_ht,prbght,htbg,tpbg,uwbg
     .,vwbg,shbg,htbg_sfc,prbg_sfc,mslpbg,tpbg_sfc,shbg_sfc
     .,uwbg_sfc,vwbg_sfc,istatus)

            if(istatus.eq.1)then
               call dprep_laps(i4time,nx_bg,ny_bg,nzbg_ht,prbght
     .,htbg,tpbg,uwbg,vwbg,shbg,htbg_sfc,prbg_sfc,mslpbg,tpbg_sfc
     .,shbg_sfc,uwbg_sfc,vwbg_sfc,istatus)
               call get_standard_latitudes(lat1,lat2,istatus)
               call get_standard_longitude(lon0,istatus)
               call get_laps_corners(nx_bg,ny_bg,sw,ne)
               istatus_dprep = 0
            endif

         elseif(bgmodel.eq.2)then
            ext = '.E48'
            call dprep_eta_conusc(nx_bg,ny_bg,nzbg_ht ,htbg
     +,prbght,tpbg,uwbg,vwbg,shbg,htbg_sfc,prbg_sfc, shbg_sfc
     +,tpbg_sfc, uwbg_sfc, vwbg_sfc, mslpbg ,istatus)

         elseif(bgmodel.eq.4)then
            ext='.SBN'
            
         elseif(bgmodel.eq.5)then
            ext='.RUC'
            call dprep_ruc2_pub(nx_bg,ny_bg,nzbg_ht,nzbg_uv
     +,nzbg_sh,nzbg_ww,htbg,prbght,shbg,uwbg,vwbg,tpbg,gproj)
         else
            print*, 'ERROR bgmodel may not be supported ',bgmodel
            stop
         endif

         if(istatus_dprep.eq.0)then
            filename = fname//af(3:4)//'00'//ext
            call dprep_output(outdir,filename
     +,cmodel,gproj,sw,ne,nx_bg,ny_bg,nzbg_ht,nzbg_uv,nzbg_sh
     +,nzbg_ww,prbght,prbguv,prbgsh,prbgww,htbg,exbg,tpbg
     +,uwbg,vwbg,wwbg,shbg,htbg_sfc,prbg_sfc,shbg_sfc,tpbg_sfc
     +,uwbg_sfc,vwbg_sfc,mslpbg,lat1,lat2,lon0,dx,dy,istatus)
            filesfound=filesfound+1
         endif
      enddo

      deallocate (htbg, tpbg, shbg, exbg)
      deallocate (uwbg, vwbg, wwbg)
      deallocate (prbght,prbguv,prbgsh,prbgww)
      deallocate (htbg_sfc,prbg_sfc,shbg_sfc
     +,uwbg_sfc,vwbg_sfc,tpbg_sfc,mslpbg)


      return
      end

      subroutine dprep_output(outdir,fname,cmodel,gproj
     +,sw,ne,nx,ny,nzht,nzuv,nzrh,nzww,prht,pruv,prrh,prww
     +,ht,ex,th,uw,vw,ww,rh
     +,ht_sfc,p_sfc,rh_sfc,th_sfc,uw_sfc,vw_sfc
     +,mslp,lat1,lat2,lon0,dx,dy,istatus)


      implicit none
      character*(*) outdir, fname, cmodel, gproj
      character*256 outfile
      integer nx,ny,nzht,nzuv,nzrh,nzww
      integer nsfcfld, istatus,  n

      real lat1, lat2, lon0, rota
      real sw(2), ne(2)
      real dx,dy

      real, intent(in) ::  ht( :, :, : )
      real, intent(in) ::  ex( :, :, : )
      real, intent(in) ::  rh( :, :, : )
      real, intent(in) ::  th( :, :, : )
      real, intent(in) ::  uw( :, :, : )
      real, intent(in) ::  vw( :, :, : )
      real, intent(in) ::  ww( :, :, : )

      real, intent(in) ::  prht( :, :, : )
      real, intent(in) ::  pruv( :, :, : )
      real, intent(in) ::  prrh( :, :, : )
      real, intent(in) ::  prww( :, :, : )


      real, intent(in) ::  ht_sfc( :, :)
      real, intent(in) ::   p_sfc( :, :)
      real, intent(in) ::  rh_sfc( :, :)
      real, intent(in) ::  th_sfc( :, :)
      real, intent(in) ::  uw_sfc( :, :)
      real, intent(in) ::  vw_sfc( :, :)
      real, intent(in) ::    mslp( :, :)

      real  adum2d( NX,  NY)



c          rh_sfc( NX,  NY), 
c    +     th( NX,  NY,  NZ), th_sfc( NX,  NY), 
c    +     uw( NX,  NY,  NZ), uw_sfc( NX,  NY), 
c    +     vw( NX,  NY,  NZ), vw_sfc( NX,  NY),
c    +     mslp(NX, NY),      adum3d( NX,  NY, NZ)


c     call s_len(model_src,lenm)

c     if(model_src(1:lenm).eq.'mm5')then
c        call output_pregrid_format(cmodel,ex, th, ht
c    & ,uw, vw, rh, mslp
c    & ,adum3d, adum3d, adum3d, adum3d, adum3d, adum2d)
c     endif

      n=index(outdir,' ')-1
      outfile=outdir(1:n)//'/'//fname
      n=index(outfile,' ')-1
      print *,'DPREP data ---> ',outfile(1:n)

      open(1,file=outfile(1:n),status='unknown',form='unformatted')
c
c *** Prepare (Write) header stuff.
c
      nsfcfld=7
      write(1) nx,ny,nzht,nzuv,nzrh,nzww,nx,ny,nsfcfld,gproj

      if(gproj.eq.'PS')then
         rota=0.
         write(1) nx,ny,nzht,lat2,lon0,rota,sw,ne
      elseif(gproj.eq.'LC')then
         write(1) nx,ny,nzht,lat1,lat2,lon0,sw,ne
c      else if(gproj.eq.'CE') then
c         write(1) nx,ny,nzht,lat1,lat2,lon0,sw,ne
      else
         print*,'Unsupported grid projection in dprep.f'
         stop
      endif
c
c *** Write surface data.
c
      write(1) uw_sfc
      write(1) vw_sfc
      write(1) th_sfc
      write(1) ht_sfc
      write(1) rh_sfc
      write(1) p_sfc
      write(1) mslp
c
c *** Write upper air data 
c
      write(1) uw
      write(1) vw
c     write(1) ww
      write(1) th
      write(1) ht
      write(1) rh
      write(1) ex

      print*, (ex(4,4,n),n=1,nzht)
c
      close(1)
c
      return
      end
