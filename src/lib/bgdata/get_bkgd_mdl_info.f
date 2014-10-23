      subroutine get_bkgd_mdl_info(bgmodel,cmodel,fullname
     &,nxbg,nybg,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww
     &,gproj,dlat,dlon,centrallat,centrallon,dxbg,dybg
     &,Lat0,Lat1,Lon0,sw,ne,cgrddef,istatus)
c
c JSmart 04-2001
c
c     USE laps_static

      implicit none

      include 'grid_fname.cmn'
      include 'bgdata.inc'

      character*256 fullname
      character*200 fpathname, fpathname_save
      character*200 fpathname_a(maxbgmodels)
      character*200 cfname_internal
      character*200 outdir
      character*200 vtable,headers_file
      character*132 cmodel
      character*30  projname
      character*13  fname13
      character*13  fname9_to_wfo_fname13
      character*4   cf
      character*2   gproj, gproj_save
      character*1   cgrddef,cgrddef_save
      
      integer       i,j,l
      integer       istatus
      integer       nxbg,nybg,nzbg,nxbg_save,nybg_save
      integer       nzbg_ht,nzbg_ht_save
      integer       nzbg_tp
      integer       nzbg_sh
      integer       nzbg_uv
      integer       nzbg_ww
      integer       lenfn,nclen,leng,lenh,lenhf,nheaders,lenp,iheader
      integer       bgmodel
      integer       record
      integer       n_valtimes

      real          Lat0,Lat1,Lat0_save,Lat1_save
      real          Lon0,Lov,Lon0_save
      real          La1in,La2in
      real          Lo1in,Lo2in
      real          La1,Lo1,La2,Lo2
      real          dlat,dlon,dlat_save,dlon_save
      real          sw(2),ne(2),sw_save(2),ne_save(2)
      real          latdxdy,londxdy
      real          rlon00,rlat00
      real          latnxny,lonnxny
      real          centrallat,centrallon
      real          dxbg,dybg
      real          rotation

      logical       cross_dateline,cross_dateline_save, l_valid_header 

      save fpathname_save,nxbg_save, nybg_save, nzbg_ht_save,
     &     gproj_save,dlat_save,dlon_save,Lat0_save,Lat1_save,Lon0_save,       
     &     cgrddef_save,cross_dateline_save,
     &     sw_save,ne_save

      data fpathname_save /'NULL'/

      interface

        subroutine get_eta48_dims(filename,NX,NY,NZ
     &,StdLat1,StdLat2,Lon0,La1,Lo1,La2,Lo2,istatus)
          character*200 filename
          integer NZ, NX, NY 
          integer istatus
          real    StdLat1,StdLat2
          real    Lon0
          real    La1,Lo1
          real    La2,Lo2

        end subroutine

        subroutine get_ruc2_dims(filename,cmodel,NX,NY,NZ
     &,StdLat1,StdLat2,Lon0,La1,Lo1,La2,Lo2,istatus)
          character*200 filename
          character     cmodel*(*)
          integer NZ, NX, NY
          integer istatus
          real    StdLat1,StdLat2
          real    Lon0
          real    La1,Lo1
          real    La2,Lo2
        end subroutine

        subroutine get_attribute_sbn(cdfname,centralLat,centralLon,
     &rlat00,rlon00,latNxNy,lonNxNy,latdxdy,londxdy,dx,dy,nx,ny,
     &rotation,projname,istatus)
          character cdfname*200
          character projname*30
          integer   nx,ny
          real      centralLat
          real      centralLon
          real      rlat00
          real      rlon00
          real      dx,dy
          real      latNxNy
          real      lonNxNy
          real      latdxdy
          real      londxdy
          real      rotation
        end subroutine

        subroutine get_sbn_dims(cdfname,cmodel
     +,nxbg,nybg,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww
     +,n_valtimes,istatus)
          character*132 cmodel
          character*200 cdfname
          integer nxbg,nybg
          integer nzbg_ht
          integer nzbg_tp
          integer nzbg_sh
          integer nzbg_uv
          integer nzbg_ww
          integer n_valtimes 
          integer istatus
        end subroutine

        subroutine readavnpublicdims(fname,x,y,numIsoLevel,record,
     +istatus)
          character*200 fname
          integer       numIsoLevel
          integer       record, x, y
          integer       nf_fid, nf_vid, nf_status
          integer       istatus
        end subroutine

        subroutine read_lapsprd_attr(fullname,
     +Dx, Dy, La1, Lo1, Latin1, Latin2, LoV, 
     +grid_type, La2,Lo2, istatus)
          integer istatus
          real Dx, Dy, La1, La2, Lo1, Lo2, LoV
          real Latin1, Latin2
          character*30 grid_type
          character*(*) fullname
        end subroutine

      end interface


      print*,'Here: get_bkgd_mdl_info'
      call s_len(cmodel,nclen)
      print*,'cmodel = ',TRIM(cmodel)
      print*,'-----------------------'

      istatus=-4	! Yuanfu: change it to -4 for exit if no bgmodel and cmodel match bk data

!     write(*,*) " grib fullname", trim(fullname)
      call s_len(fullname,lenfn)
      write(*,*) " grib fullname(1:lenfn)", fullname(1:lenfn)

      if(bgmodel.eq.0)then 
       if(cmodel(1:nclen).eq.'MODEL_FUA')then
          cfname_internal=fullname(1:lenfn)//".fua"
          call getdims_lapsprd(cfname_internal,nxbg,nybg,nzbg,istatus)
          if(istatus.ne.1)then
             print*,'error returned: getdims_lapsprd'
             return
          endif
          call read_lapsprd_attr(cfname_internal, 
     +     dxbg, dybg, La1, Lo1, La1in, La2in, LoV,
     +     projname, La2,Lo2, istatus)
          if(istatus.ne.1)then
             print*,'error returned: read_lapsprd_attr'
             return
          endif

c this code gets domain info from an existing "static" file
c        call find_domain_name(generic_data_root,grid_fnam_common,
c    &istatus)
c        call get_horiz_grid_spec(generic_data_root)
c        call s_len(grid_type,leng)

c        call s_len(projname,leng)
         
         search_length: do l=30,1,-1
            if(projname(l:l).ne.' ')then
               exit search_length
            endif
         enddo search_length

         if(l.eq.1)then
            print*,'Unable to determine string length'
            print*,'get_mdl_bkgd_info: ',cmodel
            return
         endif

! Hongli Jiang add the tangentail option. 9/4/2012 

         if(projname(1:l).eq. 'polar'.or.
     &projname(1:l).eq.'polar stereographic')then
            gproj='PS'
         elseif(projname(1:l).eq.'lambert')then
            gproj='LC'
         elseif(projname(1:l).eq. 'mercator')then
            gproj='MC'
         elseif(projname(1:l).eq.'secant lambert conformal')then
           gproj='LC'
         elseif(projname(1:l).eq.'tangential lambert conformal')then
           gproj='LC'
         else
           print*,'ERROR: unable to determine gproj setting ',gproj
           print*,'ERROR: in get_bkgd_mdl_info. '
           print*,'ERROR: l/projname ',l,projname(1:l)
           return
         endif

         if(Lo1.gt.180)Lo1=Lo1-360
         if(Lo2.gt.180)Lo2=Lo2-360
         if(LoV.gt.180)LoV=LoV-360
         nzbg_ht=nzbg
         nzbg_tp=nzbg
         nzbg_sh=nzbg
         nzbg_uv=nzbg
         nzbg_ww=nzbg
         sw(1)=La1
         sw(2)=Lo1
         ne(1)=La2
         ne(2)=Lo2
         Lon0=LoV
         Lat0=La1in
         centrallat=La1in
         centrallon=LoV
         Lat1=La2in

       elseif(cmodel(1:nclen).eq.'LAPS_FUA'.or.
     +        cmodel(1:nclen).eq.'LAPS') then

         call get_laps_dimensions(nzbg,istatus)
         call get_grid_dim_xy(nxbg,nybg,istatus)

         nzbg_ht=nzbg
         nzbg_tp=nzbg
         nzbg_sh=nzbg
         nzbg_uv=nzbg
         nzbg_ww=nzbg

         return

       elseif(cmodel(1:nclen).eq.'LAPS')then

         print*,'Error: lga currently not able to use LAPS'
         print*,'analysis for background ... not tested'
         print*,'Modify cmodel variable in background.nl'
         stop

       endif

      endif

c ETA Public
c ----------
      if(bgmodel.eq.2.and.cmodel(1:nclen).eq.'ETA48_CONUS'.or.
     &cmodel(1:nclen).eq.'ORSM_HKO')
     &then
         call get_eta48_dims(fullname,nxbg,nybg,nzbg
     &         ,Lat0,Lat1,Lon0,La1in,Lo1in,La2in,Lo2in,istatus)
         if(istatus.eq.1)then
            if(cmodel(1:nclen).eq.'ETA48_CONUS')then
               gproj='LC'
            elseif(cmodel(1:nclen).eq.'ORSM_HKO')then
               gproj='NP'
               sw(1)=10
               sw(2)=100
               ne(1)=35
               ne(2)=128
               nxbg=113
               nybg=101
            endif
            nzbg_ht=nzbg
            nzbg_tp=nzbg
            nzbg_sh=nzbg
            nzbg_uv=nzbg
            nzbg_ww=nzbg
            centrallon=Lon0
            sw(1)=La1in
            sw(2)=Lo1in
            ne(1)=La2in
            ne(2)=Lo2in
         else
            print*,'Error - get_eta48_dims: ',fullname(1:lenfn)
         endif
      endif

c All SBN grids!
C WNI-BLS ... or Unidata 
c ----------------
      if((bgmodel.eq.4).or.(bgmodel.eq.10))then

         j=lenfn-13
         if(index(fullname(j+1:j+13),'_').eq.0 .and.
     +            fullname(j:j).eq.'/') then
            fname13=fname9_to_wfo_fname13(fullname(j+1:j+9))
c           cf=fullname(j+6:j+9)
c           fullname=fullname(1:j)//fname13//cf

            cfname_internal=fullname(1:j)//fname13            !//cf
         else
            cfname_internal = fullname
         endif

         call s_len(cfname_internal,lenfn)
         IF (bgmodel .eq. 4) THEN
           call get_sbn_dims(cfname_internal,cmodel
     +       ,nxbg,nybg,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww
     +       ,n_valtimes,istatus)
           if(istatus.ne. 1)then
              print*,'Error: get_sbn_dims'
	      print*,bgmodel,cmodel(1:nclen)
	      return
           endif

           print*,'call get_attribute_sbn'
           call get_attribute_sbn(cfname_internal,centralLat
     +       ,centralLon,rlat00,rlon00,latNxNy,lonNxNy,latdxdy
     + ,londxdy,dxbg,dybg,nxbg,nybg,rotation,projname,istatus)

           if(istatus.ne. 1)then
              print*,'Error: get_attribute_sbn'
              print*,bgmodel,cmodel(1:nclen)
              return
           endif
 
c set projection type for gridconv.f

           call s_len(projname,leng)

           if(projname(1:leng).eq.'LAMBERT_CONFORMAL')gproj='LC'
           if(projname(1:leng).eq.'STEREOGRAPHIC')gproj='PS'
           if(projname(1:leng).eq.'CYLINDRICAL_EQUIDISTANT')gproj='LE'

           if(cmodel(1:nclen).eq.'RUC40_NATIVE'.or.
     .        cmodel(1:nclen).eq.'ETA48_CONUS')then

              nzbg_tp=nzbg_tp-1
              nzbg_uv=nzbg_uv-1
              nzbg_sh=nzbg_sh-1
              print*,'Retrieved SBN attributes for ',cmodel(1:nclen)

           elseif(cmodel(1:nclen).eq.'AVN_SBN_CYLEQ')then

c for global AVN, nav code expects grid 1,1 in nw corner
              print*,'set return variables'
              rlat00 =-1.*rlat00
              latNxNY=-1.*latNxNy
              nzbg_tp=nzbg_tp-2
              nzbg_uv=nzbg_uv-2
              nzbg_sh=nzbg_sh-2
              print*,'Retrieved SBN attributes for ',cmodel(1:nclen)
           elseif(cmodel(1:7).eq.'MesoEta')then
              print*,'Retrieved SBN attributes for ',cmodel(1:nclen)
           else 
              print*,'Unknown SBN model type: cmodel = ',cmodel
           endif

           Lon0=centralLon
           Lat0=centralLat
           Lat1=Lat0 ! This has be the second latitude 
                     !(tangent lambert) since no Lat1.
           dlat=dxbg/111.1
           dlon=dybg/111.1
           sw(1)=rlat00 
           sw(2)=rlon00
           ne(1)=latNxNy
           ne(2)=lonNxNy

        elseif(bgmodel.eq.10) THEN  ! WNI-BLS

           call get_unidata_grid(cfname_internal,cmodel
     +       ,nxbg,nybg,Lat0,Lat1,Lon0,La1in,Lo1in,La2in,Lo2in
     +       ,dxbg,dybg,gproj,istatus)
                                           
           if (gproj .EQ. "LL") THEN
             cgrddef = 'S'
             dlat = dybg
             dlon = dxbg
             lat0 = la1in
             lon0 = lo1in
             print *, "UNIDATA LL: dlat/dlon/lat0/lon0 =",
     +          dlat,dlon,lat0,lon0
           endif      
           if(istatus.ne. 1)then
              print*,'Error: get_unidata_grid'
              print*,bgmodel,cmodel(1:nclen)
              return
           endif
           sw(1) = La1in
           sw(2) = Lo1in
           ne(1) = La2in
           ne(2) = Lo2in

           call get_unidata_dims(cfname_internal,cmodel
     +       ,nxbg,nybg,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww
     +       ,n_valtimes,istatus)
           if(istatus.ne. 1)then
              print*,'Error: get_unidata_dims'
              print*,bgmodel,cmodel(1:nclen)
              return
           endif
         print *,'nx / ny / nzbg_ht = ',nxbg ,nybg, nzbg_ht
         print *,'corners  = ', sw(1),sw(2),ne(1),ne(2)
         print *,' gproj = ', gproj

        endif

      endif
c RUC Public
c ----------
      if(bgmodel.eq.5.and.cmodel(1:nclen).eq.'RUC40_NATIVE'
     &                .or.cmodel(1:nclen).eq.'RUC20_NATIVE')then
         call get_ruc2_dims(fullname,cmodel,nxbg,nybg,nzbg
     &              ,Lat0,Lat1,Lon0,La1in,Lo1in,La2in,Lo2in,istatus)
         if(istatus.eq.1)then
            nzbg_ht=nzbg
            nzbg_tp=nzbg
            nzbg_sh=nzbg
            nzbg_uv=nzbg
            nzbg_ww=nzbg
            Lon0=Lon0-360.
            gproj='LC'
            sw(1)=La1in
            sw(2)=-Lo1in    !-360.
            ne(1)=La2in
            ne(2)=Lo2in
         else
            print*,'Error - get_ruc2_dims: ',fullname(1:lenfn)
         endif
      endif
      
c AVN Public
c ----------
      if(bgmodel.eq.6.and.cmodel(1:nclen).eq.'AVN_FSL_NETCDF')
     &then
         call readavnpublicdims(fullname,nxbg,nybg,nzbg,record
     +,istatus)
         if(istatus.eq.1)then
            nzbg_ht=nzbg
            nzbg_tp=nzbg
            nzbg_sh=nzbg
            nzbg_uv=nzbg
            nzbg_ww=nzbg
            gproj='LL'
            Lat0=90.0   !although consistent with AVN Public, doesn't work
            Lon0=0.0    !with gridconv latlon_2_llij
            dlat=1.0 
            dlon=1.0
            cgrddef='N'
         else
            print*,'Error - readavnpublicdims '
         endif

      endif

c AVN binary (at AFWA)
c --------------------
      if(bgmodel.eq.6.and.cmodel(1:nclen).eq.'AVN_AFWA_DEGRIB')
     &then
         nxbg=360
         nybg=181
         nzbg=26
         nzbg_ht=nzbg
         nzbg_tp=nzbg
         nzbg_sh=nzbg
         nzbg_uv=nzbg
         nzbg_ww=nzbg
         gproj='LL'
         Lat0=-90.0
         Lon0=0.0
         dlat=1.0
         dlon=1.0
         cgrddef='S'
      endif

c ECMWF 
c --------------------
      if(bgmodel.eq.12)then
       if(cmodel(1:nclen).eq.'FMI_NETCDF_LL')
     &then
         gproj='LL'
         nxbg=81
         nybg=69
         nzbg=16
         nzbg_ht=nzbg
         nzbg_tp=nzbg
         nzbg_sh=nzbg
         nzbg_uv=nzbg
         nzbg_ww=nzbg
         Lat0=72.0  !NW Corner
         Lon0=15.0   !   "
         sw(1)=54.75
         sw(2)=15.0
         ne(1)=72.
         ne(2)=35.25
         dlat=0.25
         dlon=0.25
         cgrddef='N'
       endif
      endif

c Taiwan FA and NF models: Added NF_15KM, GFS, NF45 and TFS 4-27-04:
c -----------------------------------------------------------------
      if(bgmodel.eq.3)then
       if(cmodel(1:nclen).eq.'CWB_20FA_LAMBERT_RE')then
         nxbg  = 91
         nybg  = 91
         nzbg  = 16
         nzbg_ht=nzbg
         nzbg_tp=nzbg
         nzbg_sh=nzbg
         nzbg_uv=nzbg
         nzbg_ww=nzbg
         gproj='LC'
         Lat0=10.0
         Lat1=40.0
         Lon0=+120.
         sw(1)=15.879
         sw(2)=+112.545
         ne(1)=32.384
         ne(2)=+131.172
         istatus = 1
       elseif(cmodel(1:nclen).eq.'CWB_20FA_LAMBERT_NF')then
         nxbg = 145
         nybg = 139
         nzbg = 11
         nzbg_ht=nzbg
         nzbg_tp=nzbg
         nzbg_sh=nzbg
         nzbg_uv=nzbg
         nzbg_ww=nzbg
         gproj='LC'
         Lat0=10.0
         Lat1=40.0
         Lon0=+120.
         sw(1)=15.80
         sw(2)=+109.24
         ne(1)=34.987
         ne(2)=+131.60
         istatus = 1
       elseif(cmodel(1:nclen).eq.'CWB_20FA_LAMBERT_NF15'.or.
     &        cmodel(1:nclen).eq.'CWB_20FA_LAMBERT_GFS')then
         nxbg = 181
         nybg = 193
         nzbg = 11
         nzbg_ht=nzbg
         nzbg_tp=nzbg
         nzbg_sh=nzbg
         nzbg_uv=nzbg
         nzbg_ww=nzbg
         gproj='LC'
         Lat0=10.0
         Lat1=40.0
         Lon0=+120.
         sw(1)=9.28194
         sw(2)=+109.7727
         ne(1)=35.26665
         ne(2)=+137.7342
         istatus = 1
       elseif(cmodel(1:nclen).eq.'CWB_20FA_LAMBERT_NF45')then
         nxbg = 221
         nybg = 127
         nzbg = 11
         nzbg_ht=nzbg
         nzbg_tp=nzbg
         nzbg_sh=nzbg
         nzbg_uv=nzbg
         nzbg_ww=nzbg
         gproj='LC'
         Lat0=10.0
         Lat1=40.0
         Lon0=+120.
         sw(1)=-5.34068
         sw(2)=+77.9186
         ne(1)=42.92812
         ne(2)=+180.2034
         istatus = 1
       elseif(cmodel(1:nclen).eq.'CWB_20FA_LAMBERT_TFS')then
         nxbg = 229
         nybg = 181
         nzbg = 11
         nzbg_ht=nzbg
         nzbg_tp=nzbg
         nzbg_sh=nzbg
         nzbg_uv=nzbg
         nzbg_ww=nzbg
         gproj='LC'
         Lat0=10.0
         Lat1=40.0
         Lon0=+120.
         sw(1)=-9.902
         sw(2)=+82.854
         ne(1)=+52.219
         ne(2)=+199.610
         istatus = 1
       else
         print*,'bgmodel = 3, but unknown cmodel spec'
         print*,'cmodel = ',cmodel(1:nclen)
         print*,'returning without current bkgd mdl info'
       endif

      endif

c     if (bgmodel .eq. 7) then
c        gproj='LC'
c        nx_lc=nx
c        ny_lc=ny
c        nz_lc=nz
c        lat1=25.0
c        lat1_lc=lat1
c        lat2=25.0
c        lat2_lc=lat2
c        lon0=-95.0
c        lon0_lc=lon0
c        sw(1)=12.19
c        sw(2)=-133.459
c        ne(1)=57.29
c        ne(2)=-49.3849
c     elseif (bgmodel.eq.8)then
c        gproj='LL'
c        nx_ll=nx
c        ny_ll=ny
c        nz_ll=nz
c        lat0=-90.0
c        lon0=0.0
c        lon0_lc=lon0
c        dlat=1.0
c        dlon=1.0
c     endif

c GRIB1 and GRIB2 
c --------------------
      if(bgmodel.eq.13)
     &then

         call get_directory('static',outdir,lenfn)    
         vtable=outdir(1:lenfn)//'Variable_Tables/Vtable.'//cmodel

!        Determine basename and compare to saved and/or read in ones
         call get_directory_length(fullname,lenfn)
         fpathname = fullname(1:lenfn)

         if(fpathname_save .ne. fpathname)then ! header info not saved in memory

!          Read header information                     
!          This is just for testing right now. It's possible that more than
!          just the header info is needed and we may still have to call
!          DEGRIB_NAV at least once to do some internal processing.

           write(6,*)' Checking header file info'
           l_valid_header = .false.
           call s_len(outdir,lenh)
           headers_file=outdir(1:lenh)//'/Variable_Tables/headers.txt'       
           call s_len(headers_file,lenhf)
           write(6,*)' Read headers file = ',headers_file(1:lenhf)
           open(22,file=headers_file(1:lenhf),status='old',err=97)
           read(22,*,err=98,end=98)nheaders
           do iheader = 1,nheaders
             read(22,11,err=98,end=98)fpathname_a(iheader)
11           format(a)
!            write(6,*)' test - header path = ',iheader,  
!    1                 fpathname_a(iheader)             
             read(22,12,err=98,end=98)nxbg,nybg,nzbg_ht,gproj,dlat,dlon
     1         ,Lat0,Lat1,Lon0
     1         ,cgrddef,cross_dateline,sw,ne
12           format(3i6,1x,a2,1x,2f11.2,3f11.5,1x,a1,1x,l1,4f11.5)

             if(cgrddef .ne. 'N' .and. cgrddef .ne. 'S')then
                 write(6,*)' WARNING: cgrddef is ',cgrddef
                 write(6,*)' Rest of header is as follows:'
                 write(6,*)nxbg,nybg,nzbg_ht,gproj,dlat,dlon
     1             ,Lat0,Lat1,Lon0
     1             ,cgrddef,cross_dateline,sw,ne
                 goto 98
             endif

             if(fpathname_a(iheader) .eq. fpathname)then ! valid header
                 l_valid_header = .true. ! can be turned off  
                 goto 99
             else
                 write(6,*)' No header match '
                 write(6,*)fpathname_a(iheader)
                 write(6,*)fpathname
             endif
           enddo ! iheader

           goto 99

 97        write(6,*)' Warning: could not open header info'
           goto 100

 98        write(6,*)' Error reading header info'
           goto 100

 99        write(6,*)' Header info successfully read in'
100        close(22)
        
           if(l_valid_header)then
             write(6,*)' Valid DEGRIB_NAV header info has been read in'
             call s_len(fpathname,lenp)
             write(6,11)fpathname(1:lenp)
             write(6,12)nxbg,nybg,nzbg_ht,gproj,dlat,dlon
     1         ,Lat0,Lat1,Lon0
     1         ,cgrddef,cross_dateline,sw,ne
             istatus = 1
 
             write(*,*) "CALL DEGRIB_NAV:" 
             call degrib_nav(fullname, vtable, nxbg, nybg, nzbg_ht,
     &       gproj,dlat,dlon,Lat0,Lat1,Lon0,cgrddef,cross_dateline,
     &       sw(1),sw(2),ne(1),ne(2),.false.,istatus)

           else ! call degrib_nav to obtain header info
!            fullname, vtable are assumed as inputs, the rest outputs

             write(*,*) "CALL DEGRIB_NAV to obtain header info:" 
             write(*,*) " grib filename", trim(fullname)

             call degrib_nav(fullname, vtable, nxbg, nybg, nzbg_ht,
     &       gproj,dlat,dlon,Lat0,Lat1,Lon0,cgrddef,cross_dateline,
     &       sw(1),sw(2),ne(1),ne(2),.true.,istatus)

             if(cgrddef .ne. 'N' .and. cgrddef .ne. 'S')then
                 write(6,*)' WARNING: undefined cgrddef, setting to N'
                 cgrddef = 'N'
             endif

!            Write header information into a file for later retrieval
             call s_len(outdir,lenh)
             headers_file=outdir(1:lenh)//'/Variable_Tables/headers.txt'       
             call s_len(headers_file,lenhf)
             write(6,*)' Write headers file = ',headers_file(1:lenhf)
             close(22) ! A bug: channel 22 is left open in LAPS; by
Yuanfu Xie
             open(22,file=headers_file(1:lenhf),status='unknown')
             nheaders = 1
             write(22,*)nheaders
             write( 6,*)nheaders
             write(22,11)fpathname
             write( 6,11)fpathname
             write(22,12)nxbg,nybg,nzbg_ht,gproj,dlat,dlon                   
     1                  ,Lat0,Lat1,Lon0,cgrddef,cross_dateline,sw,ne
             write( 6,12)nxbg,nybg,nzbg_ht,gproj,dlat,dlon                   
     1                  ,Lat0,Lat1,Lon0,cgrddef,cross_dateline,sw,ne
             close(22)

           endif ! call degrib_nav

!          Save variables if we need them later
           fpathname_save = fpathname
           nxbg_save = nxbg
           nybg_save = nybg
           nzbg_ht_save = nzbg_ht
           gproj_save = gproj
           dlat_save = dlat
           dlon_save = dlon
           Lat0_save = Lat0
           Lat1_save = Lat1
           Lon0_save = Lon0
           cgrddef_save = cgrddef
           cross_dateline_save = cross_dateline
           sw_save = sw
           ne_save = ne

         else
           write(6,*)' Obtaining DEGRIB_NAV info from saved variables'
           nxbg = nxbg_save
           nybg = nybg_save
           nzbg_ht = nzbg_ht_save
           gproj = gproj_save
           dlat = dlat_save
           dlon = dlon_save
           Lat0 = Lat0_save
           Lat1 = Lat1_save
           Lon0 = Lon0_save
           cgrddef = cgrddef_save
           cross_dateline = cross_dateline_save
           sw = sw_save
           ne = ne_save
           istatus = 1

         endif

         nzbg_tp=nzbg_ht
         nzbg_sh=nzbg_ht
         nzbg_uv=nzbg_ht
         nzbg_ww=nzbg_ht

      endif

      return
      end
