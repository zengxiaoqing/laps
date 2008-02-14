c-----------------------------------------------------------------------
      subroutine readnet_nam(ncid,imx,imy,imz,iblev,
     +                    gh_sfc,gh,rh_2mFH,rh,rh_lbis,t_2mFH,
     +                    t,t_lbis,uw_10mFH,uw,uw_lbis,vw,
     +                    vw_lbis,vw_10mFH,av,pvv,p_sfc,heli,
     +                    cape_sfc,cape_lbis,cin_sfc,cin,
     +                    bli_lbis,pli_lbis,pw,emspMSL,prMSL,
     +                    cp_sfc,tp_sfc,isoLevel,boundryLevel,
     +                    valtime,reftime,origin,model,grid_type,
     +                    j_dim,i_dim,Ni,Nj,La1,La2,Lo1,Lo2,
     +                    Di,Dj,IntLat1,IntLat2,Lon0,start,
     +                    count,vdims,strbuf,num_of_ens)
c-----------------------------------------------------------------------

      include 'netcdf.inc'

c     Define Variables.
c     Variable ids run sequentially from 1 to nvars =   50

      parameter      (nvars=50)          ! number of variables
      parameter      (nrec=1)           ! change this to generalize
      parameter      (ndims=7)          ! number of dimensions
      parameter      (record=1)
      parameter      (namelen=132)
      parameter      (nav=1)  

      integer      rcode                     ! error code
      integer      recdim                    ! record dimension
      integer        imx,imy,imz,iblev
      integer        num_of_ens   

      real         gh_sfc(imx,imy,nrec)
      real         gh(imx,imy,imz,nrec)
      real         rh_2mFH(imx,imy,nrec)
      real         rh(imx,imy,imz,nrec)
      real         rh_lbis(imx,imy,iblev,nrec)
      real         t_2mFH(imx,imy,nrec)
      real         t(imx,imy,imz,nrec)
      real         t_lbis(imx,imy,iblev,nrec)
      real         uw_10mFH(imx,imy,nrec)
      real         uw(imx,imy,imz,nrec)
      real         uw_lbis(imx,imy,iblev,nrec)
      real         vw(imx,imy,imz,nrec)
      real         vw_lbis(imx,imy,iblev,nrec)
      real         vw_10mFH(imx,imy,nrec)
      real         av(imx,imy,imz,nrec)
      real         pvv(imx,imy,imz,nrec)
      real         p_sfc(imx,imy,nrec)
      real         heli(imx,imy,nrec)
      real         cape_sfc(imx,imy,nrec)
      real         cape_lbis(imx,imy,nrec)
      real         cin_sfc(imx,imy,nrec)
      real         cin(imx,imy,nrec)
      real         bli_lbis(imx,imy,nrec)
      real         pli_lbis(imx,imy,nrec)
      real         pw(imx,imy,nrec)
      real         emspMSL(imx,imy,nrec)
      real         prMSL(imx,imy,nrec)
      real         cp_sfc(imx,imy,nrec)
      real         tp_sfc(imx,imy,nrec)

      integer      isoLevel(imz)
      integer      boundryLevel(iblev)

      real*8         valtime(nrec)
      real*8         reftime(nrec)
      character*1    origin(namelen) 
      character*1    model(namelen) 

      character*1    grid_type(namelen,nrec)
      character*1    j_dim(namelen,nrec)
      character*1    i_dim(namelen,nrec)

      integer*2      version
      integer*2      Ni(nav) 
      integer*2      Nj(nav) 

      real         La1(nav) 
      real         La2(nav) 
      real         Lo1(nav) 
      real         Lo2(nav)
      real         Di(nav) 
      real         Dj(nav) 
      real         IntLat1(nav) 
      real         IntLat2(nav) 
      real         Lon0(nav) 

      integer      start(ndims)            ! hyperslab starting index
      integer      count(ndims)            ! hyperslab count from start
      integer        vdims(ndims)            ! max # of var dims
      character*1024 strbuf                  ! string buffer for var
                                             !  and attr names



c     Get info on the record dimension for this file.
      call ncinq(ncid,ndims,nvars,ngatts,recdim,rcode)
      call ncdinq(ncid,recdim,strbuf,nrecs,rcode)
c     nrecs now contains the # of records for this file

c     Retrieve data for gh_sfc variable.
      call ncvinq(ncid,   1,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,   1,start,count,gh_sfc,rcode)

c     Retrieve data for gh variable.
      call ncvinq(ncid,   2,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,   2,start,count,gh,rcode)

c     Retrieve data for rh_2mFH variable.
      call ncvinq(ncid,   3,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,   3,start,count,rh_2mFH,rcode)

c     Retrieve data for rh variable.
      call ncvinq(ncid,   4,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,   4,start,count,rh,rcode)

c     Retrieve data for rh_lbis variable.
      call ncvinq(ncid,   5,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,   5,start,count,rh_lbis,rcode)

c     Retrieve data for t_2mFH variable.
      call ncvinq(ncid,   6,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,   6,start,count,t_2mFH,rcode)

c     Retrieve data for t variable.
      call ncvinq(ncid,   7,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,   7,start,count,t,rcode)

c     Retrieve data for t_lbis variable.
      call ncvinq(ncid,   8,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,   8,start,count,t_lbis,rcode)

c     Retrieve data for uw_10mFH variable.
      call ncvinq(ncid,   9,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,   9,start,count,uw_10mFH,rcode)

c     Retrieve data for uw variable.
      call ncvinq(ncid,  10,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  10,start,count,uw,rcode)

c     Retrieve data for uw_lbis variable.
      call ncvinq(ncid,  11,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  11,start,count,uw_lbis,rcode)

c     Retrieve data for vw variable.
      call ncvinq(ncid,  12,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  12,start,count,vw,rcode)

c     Retrieve data for vw_lbis variable.
      call ncvinq(ncid,  13,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  13,start,count,vw_lbis,rcode)

c     Retrieve data for vw_10mFH variable.
      call ncvinq(ncid,  14,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  14,start,count,vw_10mFH,rcode)

c     Retrieve data for av variable.
      call ncvinq(ncid,  15,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  15,start,count,av,rcode)

c     Retrieve data for pvv variable.
      call ncvinq(ncid,  16,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  16,start,count,pvv,rcode)

c     Retrieve data for p_sfc variable.
      call ncvinq(ncid,  17,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  17,start,count,p_sfc,rcode)

c     Retrieve data for heli variable.
      call ncvinq(ncid,  18,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  18,start,count,heli,rcode)

c     Retrieve data for cape_sfc variable.
      call ncvinq(ncid,  19,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  19,start,count,cape_sfc,rcode)

c     Retrieve data for cape_lbis variable.
      call ncvinq(ncid,  20,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  20,start,count,cape_lbis,rcode)

c     Retrieve data for cin_sfc variable.
      call ncvinq(ncid,  21,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  21,start,count,cin_sfc,rcode)

c     Retrieve data for cin variable.
      call ncvinq(ncid,  22,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  22,start,count,cin,rcode)

c     Retrieve data for bli_lbis variable.
      call ncvinq(ncid,  23,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  23,start,count,bli_lbis,rcode)

c     Retrieve data for pli_lbis variable.
      call ncvinq(ncid,  24,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  24,start,count,pli_lbis,rcode)

c     Retrieve data for pw variable.
      call ncvinq(ncid,  25,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  25,start,count,pw,rcode)

c     Retrieve data for emspMSL variable.
      call ncvinq(ncid,  26,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  26,start,count,emspMSL,rcode)

c     Retrieve data for prMSL variable.
      call ncvinq(ncid,  27,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  27,start,count,prMSL,rcode)

c     Retrieve data for cp_sfc variable.
      call ncvinq(ncid,  28,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  28,start,count,cp_sfc,rcode)

c     Retrieve data for tp_sfc variable.
      call ncvinq(ncid,  29,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  29,start,count,tp_sfc,rcode)

c     Retrieve data for isoLevel variable.
      call ncvinq(ncid,  30,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  30,start,count,isoLevel,rcode)

c     Retrieve data for boundryLevel variable.
      call ncvinq(ncid,  31,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  31,start,count,boundryLevel,rcode)

c     Retrieve data for valtime variable.
      call ncvinq(ncid,  32,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  32,start,count,valtime,rcode)

c     Retrieve data for reftime variable.
      call ncvinq(ncid,  33,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  33,start,count,reftime,rcode)

c     Retrieve data for origin variable.
      call ncvinq(ncid,  34,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgtc(ncid,  34,start,count,origin,lenstr,rcode)

c     Retrieve data for model variable.
      call ncvinq(ncid,  35,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgtc(ncid,  35,start,count,model,lenstr,rcode)

c     Retrieve data for version variable.
      call ncvinq(ncid,  36,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  36,start,count,version,rcode)

c     Retrieve data for grid_type variable.
      call ncvinq(ncid,  37,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgtc(ncid,  37,start,count,grid_type,lenstr,rcode)

c     Retrieve data for j_dim variable.
      call ncvinq(ncid,  38,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgtc(ncid,  38,start,count,j_dim,lenstr,rcode)

c     Retrieve data for i_dim variable.
      call ncvinq(ncid,  39,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgtc(ncid,  39,start,count,i_dim,lenstr,rcode)

c     Retrieve data for Ni variable.
      call ncvinq(ncid,  40,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  40,start,count,Ni,rcode)

c     Retrieve data for Nj variable.
      call ncvinq(ncid,  41,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  41,start,count,Nj,rcode)

c     Retrieve data for La1 variable.
      call ncvinq(ncid,  42,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  42,start,count,La1,rcode)

c     Retrieve data for La2 variable.
      call ncvinq(ncid,  43,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  43,start,count,La2,rcode)

c     Retrieve data for Lo1 variable.
      call ncvinq(ncid,  44,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  44,start,count,Lo1,rcode)

c     Retrieve data for Lo2 variable.
      call ncvinq(ncid,  45,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  45,start,count,Lo2,rcode)

c     Retrieve data for Di variable.
      call ncvinq(ncid,  46,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  46,start,count,Di,rcode)

c     Retrieve data for Dj variable.
      call ncvinq(ncid,  47,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  47,start,count,Dj,rcode)

c     Retrieve data for IntLat1 variable.
      call ncvinq(ncid,  48,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  48,start,count,IntLat1,rcode)

c     Retrieve data for IntLat2 variable.
      call ncvinq(ncid,  49,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  49,start,count,IntLat2,rcode)

c     Retrieve data for Lon0 variable.
      call ncvinq(ncid,  50,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  50,start,count,Lon0,rcode)

c
c     Begin writing statements to use the data.
c

      return
      end
