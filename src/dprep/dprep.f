      program dprep_driver
c
      implicit none
c
      integer NX,NY,NZ,NRECS          !Background grid dimensions
c
      real*4 esat,es
      common /estab/esat(15000:45000),es(15000:45000)
c      
c cmodel is really only 12 chars but the SBN netcdf carrys 132
c
      character*132 cmodel
      integer i, filesfound, nf_fid, len
      character*180 nl_file, laps_data_root
      integer bgfcnt, max_files, istat
      parameter (max_files=500)
      character*200 bgfnames(max_files)
      integer maxdprepmodels, maxtimes
      parameter (maxdprepmodels=10)
      parameter (maxtimes=24)
      integer inittimes(maxtimes), fcstlengths(maxtimes)
      integer bgmodel, oldest_forecast, max_forecast_delta
      logical use_analysis, outdir_defined, bgpath_defined
      character*9 a9
      integer i4time_now, lendr, itime, ntime
      character*180 bgpath, topofname, outdir
      character*180 initbgpaths(maxdprepmodels)
     +             ,bndybgpaths(maxdprepmodels)
      integer initbgmodels(maxdprepmodels),bndybgmodels(maxdprepmodels)
                 
      character*100 reject_files(max_files)
      integer reject_cnt
      data reject_cnt/0/
      data max_forecast_delta/18/
      data use_analysis/.true./
      data outdir_defined/.false./
      data bgpath_defined/.false./
      namelist /dprep_nl/initbgpaths,bndybgpaths,initbgmodels
     +                  ,bndybgmodels,inittimes, fcstlengths, topofname
     +                  ,outdir,max_forecast_delta,oldest_forecast

c

c_______________________________________________________________________________
c
c *** Initialize tables.
c
      
      call es_ini
c
c *** Read namelist for data sources
c
      do i=1,180
         laps_data_root(i:i) = ' '
      enddo
      call GETENV('LAPS_DATA_ROOT',laps_data_root) 
      call s_len(laps_data_root,lendr)
      if(lendr.le.1) then
         nl_file(1:8) = 'dprep.nl'
      else
         nl_file = laps_data_root(1:lendr)//'/static/dprep.nl'
      endif
      
      open(21,FILE=nl_file,status='old',err=900)
      read(21,dprep_nl,err=901)
      close(21)
      call s_len(outdir,len)
      if(len.gt.0) inquire(FILE=outdir,EXIST=outdir_defined)
      if(.not.outdir_defined) then
         outdir=laps_data_root(1:lendr)//'/lapsprd/dprep'
      endif


      filesfound = 0
      i = 1
      call get_systime(i4time_now,a9,istat)

      read(a9(6:7),'(i2)') ntime
      itime=1
      do while(inittimes(itime).ne.ntime .and. itime.lt.maxtimes)
         itime=itime+1
      enddo   
      

      if(inittimes(itime).eq.ntime) then 
         print*, 'Dprep at time ',inittimes(itime),'for forecast of',
     +        fcstlengths(itime), ' hours'

         do while(filesfound.le.0 .and. i.le.maxdprepmodels)
            bgpath_defined=.false.
            bgmodel = initbgmodels(i)
            bgpath = initbgpaths(i)
            call s_len(bgpath,len)
            if(len.gt.0) inquire(FILE=bgpath,EXIST=bgpath_defined)
            if(bgmodel.eq.0.and..not.bgpath_defined) then
               bgpath=laps_data_root(1:lendr)
            else if(.not.bgpath_defined) then
               print*,'Could not find directory ',bgpath
               stop
            endif




            call get_acceptable_files(i4time_now,bgpath,bgmodel
     +           ,bgfnames,max_files,oldest_forecast,0,use_analysis,
     +           bgfcnt,-1,cmodel,nx,ny,nz,reject_files,reject_cnt)

            i=i+1
            if(bgfcnt.ge.1) then
c     print*,nx,ny,nz,i4time_now,bgmodel,bgpath,
c     +         bgfnames,bgfcnt,outdir,filesfound,bgfcnt
               call dprep_sub(nx,ny,nz,i4time_now,bgmodel,bgpath,
     +              bgfnames,bgfcnt,outdir,filesfound)
            endif

         
         enddo


         if(filesfound.le.0) then
            print*,'ERROR: Could not find suitable file to initialize'
c            stop
         endif

         filesfound=0
         i=1
         do while(filesfound.le.0 .and. i.le.maxdprepmodels)
            bgmodel = bndybgmodels(i)
            bgpath = bndybgpaths(i)



            call get_acceptable_files(i4time_now,bgpath,bgmodel,
     +           bgfnames,max_files,oldest_forecast,max_forecast_delta,
     +           use_analysis,bgfcnt,fcstlengths(itime),cmodel,nx,ny,nz,
     +           reject_files,reject_cnt)


            print *, nx,ny,nz
         
            i=i+1
            if(bgfcnt.gt.0) then
               call dprep_sub(nx,ny,nz,i4time_now,bgmodel
     +              ,bgpath,bgfnames,bgfcnt,outdir,filesfound)
            endif
         
         enddo

      else
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
      subroutine dprep_sub(NX,NY,NZ,i4time,bgmodel,bgpath,bgfnames,
     +     bgfcnt,outdir,filesfound)
      implicit none
      include 'netcdf.inc'
      integer NX, NY, NZ, i4time, i, j, bgtime
      integer bgmodel, istatus, bgfcnt, k, len, filesfound
      real sw(2), ne(2)

      character*(*) bgpath, bgfnames(bgfcnt), outdir
      character*132 fullname
      character*4 ext, gproj, af
      character*9 fname
      character*13 fname13
      integer nxbg, nybg, nzbg(5),ntbg 
      real La1, La2, Lo1, Lo2, ht( NX,  NY,  NZ), 
     +   ht_sfc( NX,  NY), p_sfc( NX,  NY), 
     +   rh( NX,  NY,  NZ), rh_sfc( NX,  NY), 
     +   th( NX,  NY,  NZ), th_sfc( NX,  NY), 
     +   uw( NX,  NY,  NZ), uw_sfc( NX,  NY), 
     +   vw( NX,  NY,  NZ), vw_sfc( NX,  NY),
     +   p(NX,NY,NZ), mslp(NX, NY), ww(nx,ny,nz)

      real lat1, lat2, lon0
      filesfound = 0
      
      do k=1,bgfcnt
         call s_len(bgfnames(k),j)
         j=j-13
         fname = bgfnames(k)(j+1:j+9)
         af = bgfnames(k)(j+10:j+13)
         call i4time_fname_lp(fname,bgtime,istatus)



         call s_len(bgpath,i)
         fullname = bgpath(1:i)//'/'//bgfnames(k)

         if(bgmodel.eq.0) then   
            ext = '.LPS'
            gproj='PS'      
            af='0000'
            call get_laps_data(bgpath,
     .           i4time,nx,ny,nz,p,ht,th,uw,vw,rh,
     .           ht_sfc,p_sfc,mslp,th_sfc,rh_sfc,uw_sfc,vw_sfc,istatus)
            if(istatus.eq.1) then
               call dprep_laps(i4time,nx,ny,nz,p,ht,th,uw,vw,rh,
     .              ht_sfc,p_sfc,mslp,th_sfc,rh_sfc,uw_sfc,vw_sfc,
     .              istatus)
            
               call get_standard_latitudes(lat1,lat2,istatus)
               call get_standard_longitude(lon0,istatus)
               call get_laps_corners(nx,ny,sw,ne)
            endif
         else if(bgmodel.eq.2) then
            ext = '.E48'
            gproj='LC'
            lat1=25.0
            lat2=25.0
            lon0=-95.0
            sw(1)=12.19
            sw(2)=-133.459
            ne(1)=57.29
            ne(2)=-49.3849

            call read_eta_conusc(fullname,nx,ny,nz,
     +           ht, p, th, uw, vw, rh, ht_sfc, p_sfc,
     +           rh_sfc, th_sfc, uw_sfc, vw_sfc, mslp ,istatus)

            call dprep_eta_conusc(nx,ny,nz,
     +           ht, p, th, uw, vw, rh, ht_sfc, p_sfc,
     +           rh_sfc, th_sfc, uw_sfc, vw_sfc, mslp ,istatus)

         else if(bgmodel.eq.4) then
            ext='.SBN'
            lat1=25.0
            lat2=25.0
            lon0=-95.0
            sw(1)=12.19
            sw(2)=-133.459
            ne(1)=57.29
            ne(2)=-49.3849
            
            call get_sbn_dims(bgpath,fname,nxbg,nybg,nzbg,ntbg)
            
            call read_conus_211(bgpath,fname,af,nx,ny,nz,
     .           nxbg,nybg,nzbg,ntbg,
     .           p,ht,th,rh,uw,vw,
     .           p_sfc,uw_sfc,vw_sfc,rh_sfc,th_sfc,
     .           mslp,gproj,2,istatus)

         else if(bgmodel.eq.5) then
            ext='.RUC'
            gproj='LC'
            lat1=25.0
            lat2=25.0
            lon0=-95.0
            sw(1)=16.2810
            sw(2)=-126.1378
            ne(1)=55.4818
            ne(2)=-57.3794

            call read_ruc2_hybb(fullname,nx,ny,nz
     +           ,mslp,ht,p,rh,uw,vw,th,ww
     +           ,istatus)
            if(istatus.gt.0) then
               call dprep_ruc2_pub(nx,ny,nz
     +              ,ht,p,rh,uw,vw,th,gproj)
               

            endif



         else 
            print*, 'ERROR bgmodel may not be supported ',bgmodel
            stop
         endif
         if(istatus.gt.0) then
            call dprep_output(outdir,fname//af(3:4)//'00'//ext,gproj,
     +           sw,ne,nx,ny,nz,ht, p, th, uw, vw, rh, ht_sfc, p_sfc,
     +           rh_sfc, th_sfc, uw_sfc, vw_sfc, mslp ,lat1, lat2, lon0, 
     +           istatus)
            filesfound=filesfound+1
         endif
      enddo

      return
      end

      subroutine dprep_output(outdir,fname,gproj,sw,ne,nx,ny,nz,
     +     ht, ex, th, uw, vw, rh, ht_sfc, p_sfc,rh_sfc, th_sfc, uw_sfc, 
     +     vw_sfc, mslp, lat1, lat2, lon0, istatus)
      implicit none
      character*(*) outdir, fname, gproj
      character*200 outfile
      integer nx,ny,nz, nsfcfld, istatus, ip, jp, n
      real lat1, lat2, lon0, rota
      real sw(2), ne(2)
      real ht( NX,  NY,  NZ), ex(NX,NY,NZ),
     +     ht_sfc( NX,  NY), p_sfc( NX,  NY), 
     +     rh( NX,  NY,  NZ), rh_sfc( NX,  NY), 
     +     th( NX,  NY,  NZ), th_sfc( NX,  NY), 
     +     uw( NX,  NY,  NZ), uw_sfc( NX,  NY), 
     +     vw( NX,  NY,  NZ), vw_sfc( NX,  NY),
     +     mslp(NX, NY) 

      n=index(outdir,' ')-1
      outfile=outdir(1:n)//'/'//fname
      n=index(outfile,' ')-1
      print *,'DPREP data ---> ',outfile(1:n)

      open(1,file=outfile(1:n),status='unknown',form='unformatted')
c
c *** Write header stuff.
c
      ip=1
      jp=1
      nsfcfld=7
      write(1) nx,ny,nz,nx,ny,ip,jp,nsfcfld,gproj

      if(gproj.eq.'PS') then
         rota=0.
         write(1) nx,ny,nz,lat2,lon0,rota,sw,ne
      else if(gproj.eq.'LC') then
         write(1) nx,ny,nz,lat1,lat2,lon0,sw,ne
c      else if(gproj.eq.'CE') then
c         write(1) nx,ny,nz,lat1,lat2,lon0,sw,ne
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
      write(1) th
      write(1) ht
      write(1) rh
      write(1) ex

      print*, (ex(4,4,n),n=1,nz)
c
      close(1)
c
      return
      end


