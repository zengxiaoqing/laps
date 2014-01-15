      program qbalpe_main_sfc
c
      implicit none
c
      integer   nx,ny,nz
      integer   istatus
c_______________________________________________________________________________
c
      call get_grid_dim_xy(nx,ny,istatus)
      if (istatus .ne. 1) then
          write (6,*) 'Error getting horizontal domain dimensions'
          go to 999
      endif
      call get_laps_dimensions(nz,istatus)
      if (istatus .ne. 1) then
          write (6,*) 'Error getting vertical domain dimension'
          go to 999
      endif

      call qbalpe_stag(nx,ny,nz)
c

999   print*,'Done'

1000  end
c
c===============================================================================
c
      subroutine qbalpe_stag(nx,ny,nz)
c  This subroutine reinstates the staggering in the 
c balance package.  This is critical for maximum 
c accuracy and adjustment potential
c
c*********************
      use mem_namelist, ONLY: read_namelist_laps
      use mem_namelist, ONLY: max_pr
c*********************************edtoll
      include 'trigd.inc'
      implicit none
      include 'bgdata.inc'
      include 'grid_fname.cmn'
c
      integer   nx,ny,nz   

c
      real dx(nx,ny),dy(nx,ny),dp(nz)
     .      ,p(nz),pstag(nz),ps(nx,ny),psb(nx,ny)
     .      ,lat(nx,ny),lon(nx,ny),ter(nx,ny)
     .      ,phi(nx,ny,nz),t(nx,ny,nz)
     .      ,u(nx,ny,nz),v(nx,ny,nz),sh(nx,ny,nz)
     .      ,db(nx,ny,nz)
c    .      ,phib(nx,ny,nz),tb(nx,ny,nz)
c    .      ,ub(nx,ny,nz),vb(nx,ny,nz),shb(nx,ny,nz)
c    .      ,phibs(nx,ny,nz),tbs(nx,ny,nz)
c    .      ,ubs(nx,ny,nz),vbs(nx,ny,nz),shbs(nx,ny,nz)
c    .      ,phis(nx,ny,nz),ts(nx,ny,nz)
c    .      ,us(nx,ny,nz),vs(nx,ny,nz),shs(nx,ny,nz)
c    .      ,lapsuo(nx,ny,nz),lapsvo(nx,ny,nz) !t=t0-dt currently not used

      real, allocatable, dimension (:,:,:) ::
     .        lapsu,lapsv,lapssh,lapstemp,lapsphi,omo
     .       ,phib,  tb,  ub,  vb,  shb,  omb
     .       ,phibs, tbs, ubs, vbs, shbs, ombs
     .       ,phis,  ts,  us,  vs,  shs,  oms

      real errt(nx,ny,nz),errw(nx,ny,nz)
     .      ,pd8(nz),pd5(nz),kpd8,kpd5
      real c_errphib,c_delo,c_errphi,c_erru,c_errub

      real om(nx,ny,nz)
c    .      ,omo(nx,ny,nz)
c    .      ,oms(nx,ny,nz)
c    .      ,ombs(nx,ny,nz)
c    .      ,omb(nx,ny,nz)
c    .      ,nu(nx,ny,nz),nv(nx,ny,nz),fu(nx,ny,nz),fv(nx,ny,nz)
c    .      ,wb(nx,ny,nz)
     .      ,re,rdpdg,po,cappa
     .      ,delo,dt
     .      ,gamo
     .      ,obert,oberu,oberw,zter,moderu,moderv,modert
     .      ,erru(nx,ny,nz),errv(nx,ny,nz),errub(nx,ny,nz)
     .      ,errphi(nx,ny,nz),errphib(nx,ny,nz)
     .      ,slastu,slastv,slastt,cor
      real grid_spacing_actual_m
     .      ,grid_spacing_cen_m
     .      ,pdif,dpbl,dpblf
     .      ,u_grid,v_grid
     .      ,u_true,v_true
     .      ,dpp
     .      ,dpht     ! drop height (m)
     .      ,distnf   ! no fly radius around drop zone (m)
     .      ,errwnds
     .      ,errdds
     .      ,errmod
     .      ,errdis
     .      ,o


      real g,sumdt,omsubs,sk,bnd,ff,fo,err,rog,rod
     .      ,sumdz,sumr,sumv2,snxny,sumf,sumt,cl,sl
     .      ,sumtscl,sumkf,sumks,sldata,den,sumom2,sumomt2
     .      ,ffz,sumu,sumv

      real smsng,rdum,dd,ddmin,cx,cy
      real comega_smooth
      real rstats(7)
      real phi3dvar(nz)

c made 2d 2-20-01 JS.
      real  tau(nx,ny)
      real  ro(nx,ny)
      real,   allocatable,dimension (:,:) :: terscl
      integer,allocatable,dimension (:,:) :: ks,kf
      integer ksij,kfij,k8,k5,lmin
c
      integer   itmax,lmax
     .         ,npass
     .         ,i4time_sys
     .         ,i4time_airdrop
     .         ,i,j,k,kk,l,ll,istatus
     .         ,ii,jj,iii,icount
     .         ,io,jo

      integer   itstatus
      integer   init_timer
      integer   ishow_timer
     
      integer   lend
      integer   lends
      integer   lenvg
      integer   lenm
      integer   adv_anal_by_t_min
      integer   k1,k2
      integer   kfirst,klast
      integer   ij
      integer   idum, idist
c     integer   max_pr
      integer   n_snd,nsnd,mxz
c
      logical lrunbal
c 
      logical incl_clom

      logical lrotate/.false./
      logical larray_diag/.false./
      logical frstone,lastone
      logical l_dum
      logical setdelo0 !namelist variable that forces delo=0 if
c                       set = .true.

c************edtoll
      character*255 filename
c**************edtoll
      character*255 staticdir,sfcdir
c     character*255 generic_data_root
      character*125 comment
      character*40  vertical_grid
      character*31  staticext,sfcext
      character*10  units
      character*9   a9_time
      character*9   a9_time_airdrop
      character*8   c8_project
      character*3   cpads_type

c    Added by B. Shaw, 4 Sep 01, modified by Steve Albers
      real, allocatable :: lapsrh(:,:,:)
      real, external :: ssh, make_rh
      real shsat, lapsrh_orig
      real ubias,vbias,urms,vrms,oberr
      logical l_preserve_rh /.true./

c    Arrays for Airdrop application
      real, allocatable, dimension(:,:) :: udrop,vdrop,tdrop,rri,rrj,
     &    rrk,rrit,rrjt,rrkt,rrii,rrjj
      real, allocatable, dimension (:) :: udropc,vdropc,tdropc,rric,rrjc

      real  erf
c_______________________________________________________________________________
c
!     Initialize in case not in namelist
      erf=.01
      itmax=200  !max iterations for relaxation

      call get_balance_nl(lrunbal,adv_anal_by_t_min,cpads_type,
     .                    incl_clom,setdelo0,
     .                    c_erru, c_errub, c_errphi, c_errphib, c_delo,       
     .                    comega_smooth,erf,itmax,
     .                    istatus)
      if(istatus.ne.0)then
         print*,'error getting balance namelist'
         stop
      endif
      print*,'lrotate = ',lrotate
      call downcase(cpads_type,cpads_type)
c
c switch to run balance package or not
c
      if(.not.lrunbal)then
         print*,'Namelist value lrunbal = false '
         print*,'Balance Package not running '
         goto 999
      endif
c
c *** Get times of mass, wind and surface data.
c
      call get_systime(i4time_sys,a9_time,istatus)
      if(istatus .ne. 1)go to 999

      call get_r_missing_data(smsng,istatus)
      if(istatus .ne. 1)go to 999
c
c get pressures and determine pressure intervals.
c
      call get_pres_1d(i4time_sys,nz,p,istatus)
      call get_vertical_grid(vertical_grid,istatus)

      call s_len(vertical_grid,lenvg)
      if(vertical_grid(1:lenvg).eq.'PRESSURE')THEN!Pressure in pa
         do i=2,nz
          dp(i)=p(i-1)-p(i)
          pstag(i-1)=p(i-1)-dp(i)*0.5
         enddo
c        print*,(p(i),dp(i),i=1,nz)
      else
         print*,'vertical grid is not PRESSURE ',vertical_grid
         goto 999
      endif

      call find_domain_name(generic_data_root,staticext,istatus)
      sfcext='lsx'
      call get_directory(staticext,staticdir,lend)
      call get_directory(sfcext,sfcdir,lends)

      re=6371220.
      rdpdg=3.141592654/180.
      cappa=287.053/1004.686
      g=9.80665
      bnd=1.e-30 !value of winds on the terrain face 

      call get_domain_laps(nx,ny,staticext,lat,lon,ter
     1                    ,grid_spacing_cen_m,istatus)

      if (istatus .ne. 1) then
         print *,'Error getting laps lat, lons.'
         stop
      endif
c
      do j=1,ny
      do i=1,nx
         call get_grid_spacing_actual(lat(i,j),lon(i,j)
     1             ,grid_spacing_actual_m,istatus)
         dx(i,j)=grid_spacing_actual_m
         dy(i,j)=dx(i,j)
         
      enddo
      enddo
c
c *** Get background grids
c
      allocate(phib(nx,ny,nz),tb(nx,ny,nz),shb(nx,ny,nz)
     .,ub(nx,ny,nz),vb(nx,ny,nz),omb(nx,ny,nz))

      call get_modelfg_3d(i4time_sys,'U3 ',nx,ny,nz,ub,istatus)
      call get_modelfg_3d(i4time_sys,'V3 ',nx,ny,nz,vb,istatus)
      call get_modelfg_3d(i4time_sys,'T3 ',nx,ny,nz,tb,istatus)
      call get_modelfg_3d(i4time_sys,'HT ',nx,ny,nz,phib,istatus)
      call get_modelfg_3d(i4time_sys,'SH ',nx,ny,nz,shb,istatus)
      call get_modelfg_3d(i4time_sys,'OM ',nx,ny,nz,omb,istatus)

      if(istatus.ne.1)then
         print*,'background model frst guess not obtained'
         return
      endif

      call get_modelfg_2d(i4time_sys,'PSF',nx,ny,psb,istatus)

      if(istatus.ne.1)then
         print*,'background model frst guess not obtained'
         return
      endif

c     check to see that omb is not missing
      where(abs(omb) .gt. 100.) omb=0.
c
c *** Get LAPS 3D analysis grids.
c
      allocate (lapsu(nx,ny,nz),lapsv(nx,ny,nz)   !t=t0
     .         ,lapssh(nx,ny,nz)
     .         ,lapstemp(nx,ny,nz)
     .         ,lapsphi(nx,ny,nz)
     .         ,omo(nx,ny,nz))

c 

      if (.not.incl_clom) then

        call get_laps_3d_analysis_data_isi(i4time_sys,nx,ny,nz
     +,lapsphi,lapstemp,lapsu,lapsv,lapssh,omo,istatus)
        if (istatus .ne. 1) then
         print *,'Error in get_laps_3d_analysis_data_isi...Abort.'
         stop
        endif

      else

        call get_laps_3d_analysis_data(i4time_sys,nx,ny,nz
     +,lapsphi,lapstemp,lapsu,lapsv,lapssh,omo,istatus)
c omo is the cloud vertical motion from lco
        if (istatus .ne. 1) then
         print *,'Error getting LAPS analysis data...Abort.'
         stop
        endif

!       Smooth the cloud omega data
        if (comega_smooth .ne. 0.) then
          if (comega_smooth .eq.-1.) then
            if(grid_spacing_cen_m .le. 4000.)then
              idist = 2
            else
              idist = 1
            endif
          else
            idist = nint(comega_smooth)
          endif

          write(6,*)' Smoothing cloud omega ',comega_smooth,idist
          do k = 1,nz
            where(omo .eq. smsng)omo=0.
            call smooth2 (nx,ny,idist,omo(:,:,k))
            where(omo .eq. 0.)omo = smsng                     
          enddo ! k
        else
            write(6,*)' Skip smoothing cloud omega ',comega_smooth
        endif

      endif

c
c *** Get LAPS 2D surface pressure.
c
      call get_laps_2d(i4time_sys,sfcext,'PS ',units,
     1                  comment,nx,ny,ps,istatus)
c
c *** For Airdrop-LAPS project we want to advance
c *** analyses to the time of payload release as specified
c *** by namelist variable adv_anal_by_t_min. Subr 'advance_grids'
c *** acquires backgrounds at i4time_airdrop and uses them to
c *** advance the "laps" analysis arrays forward in time to
c *** i4time_airdrop. Background arrays are filled with data
c *** at i4time_airdrop.
c *** 
c
      call get_c8_project(c8_project,istatus)
      call upcase(c8_project,c8_project)

      if(c8_project .eq. 'AIRDROP')then

         print*
         print*,' ************************'
         print*,' ******* Airdrop  *******'
         print*,' ************************'
         print*
         print*,' Advance systime by ',
     +adv_anal_by_t_min*60, ' seconds '

         i4time_airdrop=i4time_sys+adv_anal_by_t_min*60

         call advance_grids(i4time_sys,i4time_airdrop
     .,nx,ny,nz,ub,vb,tb,phib,shb,omb,psb
     .,lapsphi,lapstemp,lapsu,lapsv,lapssh,omo,ps,istatus)
         if(istatus.ne.1)then
            print*,'Error returned: advance_grids '
            return
         endif

         call make_fnam_lp(i4time_airdrop,a9_time_airdrop,istatus)

      endif

c *** Not considering non-linear terms for now, so no need to read t0-dt. 
c     call get_laps_wind(winddir,windtime,windext,nx,ny,nz
c    .                  ,lapsuo,lapsvo,istatus)

      if(lrotate)then

         write(6,*)' Rotate u/v to grid north. Input winds are true N'
         write(6,*)' Convert input and background ht to GPM   '

         do k = 1, nz
         do j = 1, ny
         do i = 1, nx

            call uvtrue_to_uvgrid(
     1            lapsu(i,j,k),lapsv(i,j,k)
     1           ,u_grid   ,v_grid
     1           ,lon(i,j)           )
            lapsu(i,j,k) = u_grid                    
            lapsv(i,j,k) = v_grid                    
            call uvtrue_to_uvgrid(
     1             ub(i,j,k),vb(i,j,k)
     1            ,u_grid   ,v_grid
     1            ,lon(i,j)          )
            ub(i,j,k)=u_grid
            vb(i,j,k)=v_grid
            lapsphi(i,j,k)=lapsphi(i,j,k)*g
            phib(i,j,k)=phib(i,j,k)*g
         enddo
         enddo
         enddo

      else

         write(6,*)' Convert input and background ht to GPM   '
         do k = 1, nz
         do j = 1, ny
         do i = 1, nx
            lapsphi(i,j,k)=lapsphi(i,j,k)*g
            phib(i,j,k)=phib(i,j,k)*g
         enddo
         enddo
         enddo

      endif

c *** Get laps surface pressure.
c
c     call get_laps_2d(i4time_sys,sfcext,'PS ',units,
c    1                  comment,nx,ny,ps,istatus)
c
c all pressure is in pascals
c set dynamic weight del using lat and surface pressure
c some comments about analysis constants delo and tau
c delo is 100 x the inverse square of the expected balance residual. This is
c a specified parameter that is constant over the grid
c tau controls the mass distribution of any continuity adjustments.
c See paper mcGinley 1984, Contibutions to AtmosPhysicsvol 57 p527-535
c if tau is large mass adjustment occurs over shallow layers
c tau is based on scaling = griddist**2 Brunt-Vaisala Freq**4*atmheight scale**4
c / ( pressureheight of mtn**2 mean velocity **2 f**2)
c a steep slope will reduce tau so that appropriate terrain induced vertical
c motions will result. The height scale of the atm and brunt vaisala freq
c will control how quickly
c the terrain induced vertical motion will vertically penetrate.
c For the current application we will use the background omega to determine
c an appropriate tau value.
      sumdt=0.
      sumdz=0.
      sumf=0.
      sumt=0.
      den=0.
      sumom2=0.
      sumr=0.
      sumv2=0.

c set the resolvable scale based on a wavelength for scaling 
c purposes. We will use a wavelength 
c resolvable by the grid (4*grid spacing) 
c and a data resolving wavelength based on
c the nyquist interval (4 * data spacing)
c this should be automated to be computed by the actual number of obs/area
c guess number of obs over grid, put in sldata for now.
      sldata=1200.  !number of good upper air obs (sonde quality)
      sldata=sqrt(float((nx-1)*(ny-1))*dx(nx/2,ny/2)**2/sldata)
c set model scale - a mean wave resolvable by the grid
      sl=(4.*dx(nx/2,ny/2)+sldata)/2.
      fo=14.52e-5   !2*omega

      allocate (ks(nx,ny),kf(nx,ny),terscl(nx,ny))

      if(.false.)then
         call terrain_scale_vartau(nx,ny,nz,ter,lapsphi
     &,ks,kf,terscl)
      else

         pd8=p/85000.
         pd5=p/50000.
         kpd8=nz
         kpd5=nz
         do k=1,nz
            if(pd8(k).lt.kpd8.and.pd8(k).ge.1.0)then
               kpd8=pd8(k)
               k8=k
            endif
            if(pd5(k).lt.kpd5.and.pd5(k).ge.1.0)then
               kpd5=pd5(k)
               k5=k
            endif
         enddo

         call terrain_scale(nx,ny,ter,terscl(1,1))
      endif
      ks(1,1)=k8
      kf(1,1)=k5

      print*,'terrain scale = ',terscl(1,1)
      print*,'ks/p = ',ks(1,1),p(ks(1,1))
      print*,'kf/p = ',kf(1,1),p(kf(1,1))

      do i=1,nx
      do j=1,ny

c        terscl(i,j)=terscl(1,1)   !make it constant over domain for now.
         ks(i,j)=ks(1,1)
         kf(i,j)=kf(1,1)
         
         kfij=kf(i,j)
         ksij=ks(i,j)
         sumdt=sumdt+(lapstemp(i,j,kfij)*(100000./p(kfij)
     & )**cappa - lapstemp(i,j,ksij)*(100000./p(ksij))**cappa)
         sumdz=sumdz+(lapsphi(i,j,kfij)-lapsphi(i,j,ksij))/g
         ff = fo*sind(lat(i,j))
         sumf=sumf+ff
         do k=ksij,kfij
           sumt=sumt+(lapstemp(i,j,k)*(100000./p(k))**cappa)
           den =den +p(k)/287.04/lapstemp(i,j,k)
           sumom2=sumom2+omb(i,j,k)**2
           sumv2=sumv2+(lapsu(i,j,k)**2+lapsv(i,j,k)**2)
         enddo
c error terms are the inverse sq error; right now with no
c horizontal stucture. Only vertical error allowed for now.
         do k=1,nz
!           erru(i,j,k)=(0.1*(1.+float(k-1)*.10))**(-2)
!           errub(i,j,k)=(1.5*(1.+float(k-1)*.3))**(-2)
!           errphi(i,j,k)=(15.*(1.+float(k-1)*.1)*g)**(-2)
!           errphib(i,j,k)=(30.*(1.+float(k-1)*.1)*g)**(-2)

            erru(i,j,k)=(c_erru*(1.+float(k-1)*.10))**(-2)
            errub(i,j,k)=(c_errub*(1.+float(k-1)*.3))**(-2)
            errphi(i,j,k)=(c_errphi*(1.+float(k-1)*.1)*g)**(-2)
            errphib(i,j,k)=(c_errphib*(1.+float(k-1)*.1)*g)**(-2)
c
c vertical motions in clear areas come in as the missing data parameter.
c replace missing cloud vv's with background vv's. Unless there is cloud
c we seek to replicate the background vertical motions
            if(abs(omb(i,j,k)).gt.100.)omb(i,j,k)=0. ! check for missing omb
            if(abs(omo(i,j,k)).gt.100.)omo(i,j,k)=omb(i,j,k)
         enddo
      enddo
      enddo

      snxny=float(nx*ny)
      sk=float(kfij-ksij+1)
      sumt=sumt/snxny/sk
      sumdt=sumdt/snxny
      sumdz=sumdz/snxny
      sumf=sumf/snxny
      sumv2=sumv2/snxny/sk
      den=den/snxny/sk
      sumom2=sumom2/snxny/sk
      sumr=sqrt(g*sumdt/sumt/sumdz)
c     delo is scaled as 10% of expected eqn of motion residual ro*U**2/L
      rod=sqrt(sumv2)/(sumf*sldata)
      rog=sqrt(sumv2)/(sumf*sl)
      if(rog.gt.1.) rog=1.
      if(rod.gt.1.) rod=1.
!     delo=100.*sl**2/sumv2**2/rog**2           
      delo=c_delo*sl**2/sumv2**2/rog**2           
      sumomt2=sumv2**3*sumt**2*den**2/(sl**2*rog**2*sumdt**2) 
      do j=1,ny
      do i=1,nx
         kfij=kf(i,j)
         ksij=ks(i,j)
         dpp=p(ksij)-p(kfij)
c scale tau as of 1/(omega)**2 where omega is the background rms omega              
c if background is missing sumom2 will be 0. In that case scale a 
c vertical motion based on a scaled dynamic omega. Utilze the maximum 
c estimated omega for setting the flow blocking parameter tau
        if(sumom2.ne.0.) then
        omsubs=sumom2
         if(sumom2.lt.sumomt2) omsubs=sumomt2
         tau(i,j)= 1./omsubs 
        else
         tau(i,j)=1./sumomt2
        endif

      enddo
      enddo
      deallocate(ks,kf,terscl)

c
c if this is for AIRDROP there is no need to run the balance package, only 
c continuity. set delo=0. In balcon this will skip the balance sequence.
c
      if(c8_project .eq. 'AIRDROP'.or.setdelo0)then
         delo=0.
      endif
c
c print these arrays now.
      print*,'length scales(m): mean wavelength, data wavelength ',sl,
     &sldata 
      print*,'/dthet/thet/dz/den/N/V/f/delo/tau:' 
     &,sumdt,sumt,sumdz,den,sumr,sqrt(sumv2),sumf,delo,tau(1,1)
      print*,'Omegab,Omegadyn,Froude Num ',
     &'DataRossby Num,GrdRosby,aspectP/X:',sqrt(sumom2),sqrt(sumomt2),
     &    (sumr*sumdz/sqrt(sumv2)),rod,rog ,(dpp/dx(nx/2,ny/2))
      
c
c *** Compute non-linear terms (nu,nv) and execute mass/wind balance.
c *** Do for lmax iterations.
c
      lmax=1 

c put lapsphi into phi, lapsu into u, etc
      phi = lapsphi
      u   = lapsu
      v   = lapsv
      t   = lapstemp
      om  = omo
      sh  = lapssh

c
      if(larray_diag)then
         print*
         print*,'Before balstagger - analysis u/v'
         print*
         do k=1,nz
            print*
            print*,'Calling array diagnosis: ',k,' ',p(k)
            print*,'---------------------------------------'
            call array_diagnosis(u(1,1,k),nx,ny,'u-comp    ')
            call array_diagnosis(v(1,1,k),nx,ny,'v-comp    ')
         enddo
      endif
c
c stagger the LAPS grids to prepare for balancing
c
c laps analysis grids first.
      allocate (phis(nx,ny,nz),ts(nx,ny,nz),shs(nx,ny,nz)
     .,us(nx,ny,nz),vs(nx,ny,nz),oms(nx,ny,nz))

      call balstagger(u,v,phi,t,sh,om,us,vs,
     &phis,ts,shs,oms,nx,ny,nz,p,ps,1) 

      if(larray_diag)then
         print*
         print*,'After balstagger - analysis u/v'
         print*
         do k=1,nz
            print*
            print*,'Calling array diagnosis: ',k,' ',p(k)
            print*,'---------------------------------------'
            call array_diagnosis(us(1,1,k),nx,ny,'u-comp    ')
            call array_diagnosis(vs(1,1,k),nx,ny,'v-comp    ')
         enddo
      endif
c
c background grids second.
      allocate (phibs(nx,ny,nz),tbs(nx,ny,nz)
     .,shbs(nx,ny,nz),ombs(nx,ny,nz),ubs(nx,ny,nz)
     .,vbs(nx,ny,nz))
      call balstagger(ub,vb,phib,tb,shb,omb,ubs,vbs,
     &phibs,tbs,shbs,ombs,nx,ny,nz,p,ps,1) 

      deallocate(phib,shb)
   
c put phis into phi, us into u ... ie., put staggered arrays into non-staggered.
      call move_3d(phis,phi,nx,ny,nz)
      call move_3d(us,u,nx,ny,nz)
      call move_3d(vs,v,nx,ny,nz)
      call move_3d(oms,om,nx,ny,nz)
      call move_3d(shs,sh,nx,ny,nz)
      call move_3d(ts,t,nx,ny,nz) 
c
c use these constructs if cloud and backgnd vv's don't exist
c     call zero3d(om,nx,ny,nz)
c     call zero3d(ombs,nx,ny,nz)
c     call zero3d(oms,nx,ny,nz)

c
c ****** Execute mass/wind balance.
c  set maximum relaxation correction for phi in GPM
      err=1.0
c returns staggered grids of full fields u,v,phi

      call balcon(phis,us,vs,oms,phi,u,v,om,phibs,ubs,vbs,ombs
     . ,ts,rod,delo,tau,itmax,err,erru,errphi,errub,errphib
c    . ,nu,nv,fu,fv
     . ,nx,ny,nz,lat,dx,dy,ps,p,dp,lmax,erf)

      if(larray_diag)then
         print*
         print*,'After balcon'
         print*
         do k=1,nz
            print*
            print*,'Calling array diagnosis: ',k,p(k)
            print*,'-------------------------------'
            call array_diagnosis(u(1,1,k),nx,ny,' u-comp   ')
            call array_diagnosis(v(1,1,k),nx,ny,' v-comp   ')
         enddo
      endif
c ******** 10/08 ***** Ed Tollerud
c     call put_sfc_bal(i4time,t_bal,ht_bal,u_bal,v_bal         ! Input
c    1                       ,topo,ni,nj,nk                           ! Input
c    1                       ,istatus                              )  ! Output
c     call put_sfc_bal(i4time,tbs,phibs,ubs,vbs,               ! Input
c    1                       ,ter,nx,ny,nz                           ! Input
c    1                       ,istatus                              )  ! Output
c Question: no balancing on ts? is tbs 

c
c
c *** destagger and Write out new laps fields.
c
c the non-staggered grids must be input with intact boundaries from background 
c  the bs arrays are used as dummy arrays to go from stagger to non staggered
      print*,'Destaggering grids back to LAPS A grid' 
      call balstagger(ubs,vbs,phibs,tbs,
     & shbs,ombs,u,v,phi,t,sh,om,nx,ny,nz,p,ps,-1) 
c   Prior to applying boundary subroutine put non-staggered grids back into
c   u,v,om,t,sh,phi.

c      do k=1,nz/2
c      call diagnose(u,nx,ny,nz,nx/2,ny/2,k,7,'staggered U')
c      call diagnose(v,nx,ny,nz,nx/2,ny/2,k,7,'staggered V')
c      call diagnose(ubs,nx,ny,nz,nx/2,ny/2,k,7,'unstaggd U ')
c      call diagnose(vbs,nx,ny,nz,nx/2,ny/2,k,7,'unstaggd V ')
c      call diagnose(lapsu,nx,ny,nz,nx/2,ny/2,k,7,'original U ')
c      call diagnose(lapsv,nx,ny,nz,nx/2,ny/2,k,7,'original V ')
c      enddo

      if(.false.)then

c     write out input laps vs balanced laps at center grid point
      print*,'Output u and Input u and diff after balance and destagger'
      do k=1,nz
       write (6,1111) ubs(nx/2,ny/2,k),lapsu(nx/2,ny/2,k)
     &        ,ubs(nx/2,ny/2,k)-lapsu(nx/2,ny/2,k)
      enddo
      print*,' Output v and Input v after balance and destagger '
      do k=1,nz
       write (6,1111) vbs(nx/2,ny/2,k),lapsv(nx/2,ny/2,k)
     &        ,vbs(nx/2,ny/2,k)-lapsv(nx/2,ny/2,k)
      enddo
      print*,' Output T and Input T after balance and destagger '
      do k=1,nz
       write (6,1111) tbs(nx/2,ny/2,k),lapstemp(nx/2,ny/2,k)
     &        ,tbs(nx/2,ny/2,k)-lapstemp(nx/2,ny/2,k)
      enddo
      print*,' Output phi and Input phi after balance and destagger '
      do k=1,nz
       write (6,1111) phibs(nx/2,ny/2,k),lapsphi(nx/2,ny/2,k)
     &        ,phibs(nx/2,ny/2,k)-lapsphi(nx/2,ny/2,k)
      enddo
 1111 format(1x,3f9.2)

      endif
      deallocate(lapsphi,lapsu,lapsv,omo)
c temporary revision, ed Tollerud, set npass to work on 64 bit machines, 0709********
      npass=0
c end temporary revision, 0709*********
c adjust surface temps to account for poor phi estimates below ground 
      call sfctempadj(tb,lapstemp,p,ps,nx,ny,nz,npass)  

      call move_3d(phibs,phi,nx,ny,nz)
      call move_3d(ubs,u,nx,ny,nz)
      call move_3d(vbs,v,nx,ny,nz)
      call move_3d(tbs,t,nx,ny,nz)
      call move_3d(ombs,om,nx,ny,nz)
      call move_3d(shbs,sh,nx,ny,nz)

      deallocate (phibs,ubs,vbs,tbs,ombs,shbs)

c JS> commented 8-21-02.  unsure of need for ps at this point?
c     call get_laps_2d(i4time_sys,sfcext,'PS ',units,
c    1                  comment,nx,ny,ps,istatus)

      if(lrotate)then

         print*,'rotate u/v analysis and backgrnd to true north'
         print*,'Put geopotential hgt back to m'
         do k = 1, nz
            do j = 1, ny
            do i = 1, nx
               call uvgrid_to_uvtrue(
     1            u(i,j,k),v(i,j,k)
     1           ,u_true   ,v_true
     1           ,lon(i,j)           )
               u(i,j,k) = u_true
               v(i,j,k) = v_true

               phi(i,j,k)=phi(i,j,k)/g
            enddo
            enddo
         enddo

      else

         print*,'Put geopotential hgt back to m'
         do k = 1, nz
            do j = 1, ny
            do i = 1, nx
               phi(i,j,k)=phi(i,j,k)/g
            enddo
            enddo
         enddo

      endif
c
c  The (dry) dynamic balancing makes small adjustments to the temp
c  field, which changes the saturation vapor pressure, which means
c  rh may not be exactly 100% where there is cloud liquid present.
c  That would cause erroneous evaporation of production of cloud
c  liquid (and associated heating/cooling) in the very first time 
c  step.  No good.


c  This section updates a previous section to adjust the moisture
c  using specific humidity so moisture between the unbalanced and
c  balance analysis is better conserved.  B. Shaw, Sep 01

      allocate (lapsrh(nx,ny,nz))
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
c         
c       Compute the saturation specific humidity

        shsat = ssh(p(k)*0.01,t(i,j,k)-273.15)*0.001

        if(.not. l_preserve_rh)then
c           Ensure the specfic humidity does not exceed
c           the saturation value for this temperature

            lapssh(i,j,k) = MIN(shsat,lapssh(i,j,k))
c
c           Finally, rediagnose RH wrt liquid from the 
c           modified sh field

            lapsrh(i,j,k) = make_rh(p(k)*0.01,t(i,j,k)-273.15
     .                    ,lapssh(i,j,k)*1000.,-132.)*100.         
            lapsrh(i,j,k) = MAX(lapsrh(i,j,k),1.0)

        else ! Keep the RH the same as it was prior to balancing (while
             ! modifying SH)
            lapsrh_orig = make_rh(p(k)*0.01,lapstemp(i,j,k)-273.15  
     .                          ,lapssh(i,j,k)*1000.,-132.)*100.         
            lapsrh(i,j,k) = MAX(lapsrh_orig,1.0)
   
            lapssh(i,j,k) = (lapsrh_orig/100.) * shsat

        endif

      enddo
      enddo
      enddo     

      deallocate(lapstemp)
c
c     New section added by to replicate the values of u/v at 
c     from the lowest p-level still above ground to all levels
c     below ground  (B. Shaw, 12 Apr 02) 

c     do j = 1,ny
c       do i = 1, nx
c         findsfclev:  do k = 1,nz
c           IF (phi(i,j,k) .GT. ter(i,j)) EXIT findsfclev
c         enddo findsfclev
c         IF (k .gt. 1) THEN
c           u(i,j,1:k-1) = u(i,j,k)
c           v(i,j,1:k-1) = v(i,j,k)
c         ENDIF
c       enddo
c     enddo
c
c produce surface balance fields by interpolation to terrain height
c ******** 10/08 ***** Ed Tollerud
c     call put_sfc_bal(i4time,t_bal,ht_bal,u_bal,v_bal         ! Input
c    1                       ,topo,ni,nj,nk                           ! Input
c    1                       ,istatus                              )  ! Output

      call put_sfc_bal(i4time_sys,t,phi,u,v               ! Input
     1                       ,ter,nx,ny,nz                           ! Input
     1                       ,istatus                              )  ! Output
      istatus = ishow_timer()
c
c
c Write balance output (balance/lt1 and balance/lw3).
c ---------------------------------------------------
c
      if(c8_project .eq. 'AIRDROP')then
         i4time_sys = i4time_airdrop
      endif

      call write_bal_laps(i4time_sys,phi,u,v,t,om,lapsrh,lapssh
     .                   ,nx,ny,nz,p,istatus)
      if(istatus.ne.1)then
         write(6,*)'error writing balance fields'
         return
      endif

      if(c8_project .eq. 'AIRDROP')then

c ----------------------------------------------------------------
c -----------    AIRDROP ANALYSIS ERROR SECTION ------------------
c Output required: turbulent compontents of u, v, and w;
c analysis error of u,v,wc t, rh
c compute dropsone wind difference from background
c simple solution is to read from .pig and tmg files, then create a
c truth profile and compute an average analysis error

         write(6,*)' Rotate u/v analysis and background
     1    to true north for AIRDROP analysis'

         do k = 1, nz
            do j = 1, ny
            do i = 1, nx
               call uvgrid_to_uvtrue(
     1            u(i,j,k),v(i,j,k)
     1           ,u_true   ,v_true
     1           ,lon(i,j)           )
               u(i,j,k) = u_true
               v(i,j,k) = v_true
               call uvgrid_to_uvtrue(
     1             ub(i,j,k),vb(i,j,k)
     1            ,u_true   ,v_true
     1            ,lon(i,j)          )
            ub(i,j,k)=u_true
            vb(i,j,k)=v_true
            enddo
            enddo
         enddo

cedtoll         call get_wind_parms(l_dum,l_dum,l_dum
cedtoll     1,rdum,rdum,rdum,rdum,rdum,max_pr,idum,idum,istatus)
c************edtoll
c     call get_directory(staticext,staticdir,lend)
      filename=staticdir(1:lend)//'/wind.nl'
      call read_namelist_laps('wind',filename)
c*************edtoll
c
c ---------------------------------------------------------------------------
c nsnd = the max number of profiles (soundings); used for array dimensioning
c        purposes.
c n_snd= the number of soundings found in the ingest files. For "pin", this
c        is always one while for "snd" 1 <= n_snd <= max_pr, and n_snd is 
c        returned from routine readprg.
c ---------------------------------------------------------------------------
c
         if(cpads_type.eq.'pin')then
            nsnd=1
            n_snd=1
         else
            nsnd=max_pr
            n_snd=4
         endif

         allocate(udrop(nsnd,nz),udropc(nz)
     &           ,vdrop(nsnd,nz),vdropc(nz)
     &           ,tdrop(nsnd,nz),tdropc(nz)
     &           ,rri(nsnd,nz),rrj(nsnd,nz),rrk(nsnd,nz)
     &           ,rrit(nsnd,nz),rrjt(nsnd,nz),rrkt(nsnd,nz)
     &           ,rrii(nsnd,nz),rrjj(nsnd,nz)
     &           ,rric(nz),rrjc(nz))

         if(cpads_type.eq.'pin')then

            call readpig(a9_time,nx,ny,nz,lat,lon
     1,udrop,vdrop,tdrop,rri,rrj,rrit,rrjt,istatus)

         elseif(cpads_type.eq.'snd')then

            call readprg(a9_time,nx,ny,nz,nsnd
     1,udrop,vdrop,tdrop,rri,rrj,rrk,rrit,rrjt,rrkt,n_snd,istatus)

         elseif(cpads_type.eq.'pln')then

            print*,'*********************************************'
            print*,'LAPS PLANNING run for PADS Airdrop mission'
            print*,'*********************************************'
            print*
c
c calculate background density array and convert t to potential temperature
c
            io=nx/2+1
            jo=ny/2+1

            do k=1,nz
            phi3dvar(k)=phi(io,jo,k)
            do j=1,ny
            do i=1,nx
               db(i,j,k)=p(k)/287.04/tb(i,j,k)/(p(k)/100000.)**(2./7.)
                t(i,j,k)=o(t(i,j,k)-273.15,p(k)/100.)
            enddo
            enddo
            enddo
            dt=1000        ! operational dropsonde time cycle in sec
            dpht=5000     ! drop height (m)
            distnf= 20000. ! no fly radius around drop zone (m)
            errwnds=1. ! dropsonde wind error
            errdds=.1 ! model  density error in kg/m**3
            errmod=3. ! model wind error (m/sec)
            errdis=100000. ! model error displacement in m
c -------------------------------------------------------------------------
CC #### PADS Planning tool variance estimator ####
c -------------------------------------------
            call var3d (u,v,db,t ,smsng,grid_spacing_cen_m
     &,dt,dpht,distnf,p,phi3dvar,errwnds,errdds,errmod,errdis,io,jo
     &,nx,ny,nz,lat,lon,ter,a9_time_airdrop)

            goto 998

         else

            print*,'********************************************'
            print*,'   !!! Error: No PADS type indicated !!!'
            print*,'Check static/balance.nl; variable cpads_type'
            print*,'********************************************'
            return

         endif

         if(istatus.eq.-3)then
         print*,'Warning status: ',cpads_type,' istatus=',istatus
            print*,'Dropsonde data not available; using default '
     1    ,'dropsonde (profile in center of grid)'
            go to 99
         elseif(istatus.eq.-1)then
            print*,'Only tmg exists: generate u/v drop profiles
     1 from analysis with gaussian error'
         elseif(istatus.eq.-2)then
             print*,'Only prg exists: generate T drop profile
     1 from analysis with gaussian error'
         endif
c routine to create drop and rr arrays when either prg or tmg
c do not exist
         if(istatus.eq.-1)then
            rri=rrit
            rrj=rrjt
            iii=382983
            slastv=0
            slastu=0
            cor=.5
            do l=1,n_snd
            do k=1,nz
               if(rri(l,k).ne.smsng .and. rrj(l,k).ne.smsng)then
                  ii=rri(l,k)
                  jj=rrj(l,k)
                  slastu= ffz(iii,20,cor,slastu) ! assumed gaussian error 3ms 
                  udrop(l,k)= u(ii,jj,k) + slastu*3.
                  slastv= ffz(iii,20,cor,slastv) ! assumed gaussian error 3ms 
                  vdrop(l,k)= v(ii,jj,k) + slastv*3.   
               endif
            enddo
            enddo
         else ! this is the istatus= -2 option
            rrit=rri
            rrjt=rrj
            iii=382983
            slastt=0
            do l=1,n_snd
            do k=1,nz
               if(rrit(l,k).ne.smsng .and. rrjt(l,k).ne.smsng)then
                  ii=rrit(l,k)
                  jj=rrjt(l,k)
                  slastt=ffz(iii,20,cor,slastt)
                  tdrop(l,k)= t(ii,jj,k) + slastt*2.      
               endif
            enddo
            enddo
         endif
c
c this routine checks all levels for u,v, and t. If missing (either
c no dropsonde or for levels above dropsonde, dropsonde variables 
c will be proxied by the analysis at either the dropsonde points
c (if they exist) otherwise, the center of the grid and
c gaussian noise will be applied: 3 m/s for wind, 1C for temp
c this will allow a variance to be provided above dropsonde levels
c
99    slastu=0
      slastv=0
      slastt=0 
      sumu=0
      sumv=0
      print*,'k rri rrj, udrop, vdrop, u and v at drop point'
      do l=1,n_snd
       icount=0
       do k=1,nz
        if(rri(l,k).eq.smsng.or.rrj(l,k).eq.smsng)then
         icount=icount+1
        endif
       enddo
       if(icount.ne.nz)then
        frstone=.true.
        lastone=.true.
        do k=1,nz
         if(rri(l,k).ne.smsng)then
          if(frstone)then
           frstone=.false.
           kfirst=k
           if(k.gt.1)then
            do kk=1,k-1
             rri(l,kk)=rri(l,k)
             rrj(l,kk)=rrj(l,k)
            enddo
           endif
          endif
         elseif(.not.frstone)then
          if(lastone)then
           lastone=.false.
           klast=k-1
           rri(l,k)=rri(l,klast)
           rrj(l,k)=rrj(l,klast)
          else
           rri(l,k)=rri(l,klast)
           rrj(l,k)=rrj(l,klast)
          endif
         endif
        enddo
       else
        rri(l,1:nz)=float(nx)/2.
        rrj(l,1:nz)=float(ny)/2.
       endif
       do  k=1,nz
        ii=nint(rri(l,k))
        jj=nint(rrj(l,k))
        if(udrop(l,k).eq.smsng.or.vdrop(l,k).eq.smsng)then
         slastu=ffz(iii,20,cor,slastu)
         slastv=ffz(iii,20,cor,slastv)
         udrop(l,k)=u(ii,jj,k)+slastu*3.  
         vdrop(l,k)=v(ii,jj,k)+slastv*3.  
        endif
        if(tdrop(l,k).eq.smsng)then
         slastt=ffz(iii,20,cor,slastt)
         tdrop(l,k)=t(ii,jj,k)+slastt*1.  
        endif
        if(.false.)then
         write(6,1112) k,rri(l,k),rrj(l,k),udrop(l,k),vdrop(l,k),
     &   u(50,16,k),v(50,16,k),u(50,16,k)-udrop(l,k),v(50,16,k)-
     &   vdrop(l,k)
         if(rri(l,k).ne.50) sumu=u(50,16,k)-udrop(l,k)+sumu!only sum over sonde
         if(rri(l,k).ne.50) sumv=v(50,16,k)-vdrop(l,k)+sumv
        endif
       enddo 
1112   format (1x,i3,8f6.2)
       print*,' sonde ',l, 'ballistic u v diff ',sumu,sumv
      enddo ! on l for each sonde

c -----------------------------------------------------------------
c routine to find closest sonde to cneter of grid (drop point) for
c error estimation, when more than one sounding

      if(n_snd.gt.1)then
         cx=nx/2+.5 
         cy=ny/2+.5
         do k=1,nz
          ddmin=1.e30
          lmin=0
          do l=1,n_snd
           if(rri(l,k).ne.smsng.and.rrj(l,k).ne.smsng) then
             dd=(cx-rri(l,k))**2+(cy-rrj(l,k))**2
             if(dd.lt.ddmin) then
                lmin=l
                ddmin=dd
             endif
           endif
          enddo 
          if(lmin.ne.0) then
             udropc(k)=udrop(lmin,k)
             vdropc(k)=vdrop(lmin,k)
             tdropc(k)=tdrop(lmin,k)
             rric(k)=rri(lmin,k)
             rrjc(k)=rrj(lmin,k)
          endif
         enddo ! on k all levels
      else
         rric=rri(1,:)
         rrjc=rrj(1,:)
         udropc=udrop(1,:)
         vdropc=vdrop(1,:)
         tdropc=tdrop(1,:)
      endif
c
c subr profile fills array erru(u),errub(v) with analysis error 
c here we assume the following for dropsonde error: wind 1m/s, temp .5deg
c and model error
      call read_wind3d_wgi(rstats,istatus)
      if(istatus .ne. 1)then
         print*,'Using u/v model errors from climo estimates:'
         moderu=3.
         moderv=3.
         modert=1.0
      else
         print*,'Using u/v model errors from wind analysis:'
         moderu=rstats(6)
         moderv=rstats(7)
         modert=1.0
      endif
      print*,'u/v model error = ',moderu,moderv
      oberu=1.
      obert=.5
      oberw=.05

c  compute TKE
c  the s arrays are used to hold the turbulent components of u,v, and w

      call turb(u,v,om,t,phi,p,us,vs,oms,ts,ter,nx,ny,nz)

      print*
      print*,'Calling subroutine profile. '
      print*

      call profile(udropc,vdropc,tdropc,rric,rrjc,oberu,oberw
     &,obert,moderu,moderv,modert,erru,errv,errw ,errt,u,v,om,
     & lapssh,us,vs,oms,t,phi,ub,vb,omb,tb,p
     &,ter,zter,nx,ny,nz,i4time_sys)

c           call write_errors(a9_time_airdrop,p,erru,errv,errw,errt,
c    1 nx,ny,nz,rri,rrj,zter,alt)

      endif

 998  deallocate(lapsrh)

 999  return
      end

      subroutine sfctempadj(t,to,p,ps,nx,ny,nz,npass)
c routine is designed to correct surface temperatures which have been computed 
c employing the hyposometric equation using geopotential heights of pressure 
c surfaces below ground. 'to' is the laps lt1 which uses surface data near the 
c ground. We want these estimates to be in the 3-D fields rather than  fictional
c hydrostatic estimates from phi. Above the ground balanced temps are good.
c Both the temp mmediately above and below ground are given the lt1. 
c There is a super-adiabatic check in the vertical and adjustments made if 
c necessary.

      integer nx,ny,nz
      real t(nx,ny,nz),to(nx,ny,nz),p(nz),ps(nx,ny)
      cappa=2./7.
      sum=0.
      sum2=0.
      cnt2=0
      cnt=0.
      do j=1,ny  
      do i=1,nx
       do k=2,nz
        if(ps(i,j).ge.p(k)) then
         temp=to(i,j,k)-t(i,j,k)
         sum=temp+sum
         sum2=temp**2+sum2
         cnt=cnt+1
c   put lt1 temps at the levels below the surface
         ksave=k-1
         do l = ksave,1,-1
          t(i,j,l)=to(i,j,l)
         enddo
c   check for superadiabatic layers
c        potemp=t(i,j,k)*(100000./p(k))**cappa   
c        cnt1=0
c        sum1=0
c        do l=k+1,nz
c         potempk=t(i,j,l)*(100000./p(l))**cappa  
c         if(potemp.gt.potempk) then
c          cnt1=cnt1+1.
c          sum1=sum1+(potemp-potempk)
c          cnt2=cnt2+1
c         else
c          lsave=l-1
c          go to 1 ! reached top of adiabatic layer 
c         endif
c        enddo ! on l
c  1     continue
c        potemp=potemp-sum1/(cnt+1.)
c        t(i,j,k)=potemp*(p(k)/100000.)**cappa
c        do ll=k+1,lsave
c         t(i,j,ll)=potemp*(p(ll)/100000.)**cappa
c        enddo ! on ll
c        go to 2 ! go to a new grid point
        endif ! the ps gt p check
       enddo ! on k
c run smoother on temps below ground
       do n=1,npass
        do k=2,ksave
         t(i,j,k)=.5*t(i,j,k)+.25*(t(i,j,k+1)+t(i,j,k-1))
        enddo
       enddo
   2  enddo ! on i
      enddo ! on j
      print*, 'LT1 Temp Adjust Summary '
      print*, 'Bias ',sum/cnt, ' RMS ',sqrt(sum2/cnt) 
c     print*, cnt2 , ' Temps had to be adj for sup-adi lapse rates'
      return
      end 



c===============================================================================
c
      subroutine momres(u,v,phi,nu,nv,fu,fv,wa,delo
     .                 ,nx,ny,nz,lat,dx,dy,ps,p)
c
c *** Momres computes momentum residual for whole domain for the staggered
c grid. Each momentum residual (u component, v component) is computed on the 
c  u and v grids respectively
c
      implicit none
c
      integer nx,ny,nz
     .         ,nxm1,nym1,nzm1
     .         ,i,j,k,ks
c
      real u(nx,ny,nz),v(nx,ny,nz)
     .      ,phi(nx,ny,nz)
     .      ,nu(nx,ny,nz),nv(nx,ny,nz)
     .      ,fu(nx,ny,nz),fv(nx,ny,nz)
     .      ,wa(nx,ny,nz)
     .      ,lat(nx,ny),dx(nx,ny),dy(nx,ny)
     .      ,ps(nx,ny),p(nz)
     .      ,delo,g,rdpdg,fo,bnd,cnt,sum,ang,sin1,f
     .      ,nvv,nuu,errms,fuu,fvv

      real, allocatable, dimension(:,:,:) :: wb,wc
c_______________________________________________________________________________
c
      print *,'momres'
      g=9.80665
      rdpdg=3.141592654/180.
      fo=14.52e-5   !2*omega
      bnd=1.e-30
      nxm1=nx-1
      nym1=ny-1
      nzm1=nz-1
      cnt=0.
      sum=0.

      allocate (wb(nx,ny,nz),wc(nx,ny,nz))

      call zero3d(wa,nx,ny,nz)
      call zero3d(wb,nx,ny,nz)
      call zero3d(wc,nx,ny,nz)
      do j=2,nym1
      do i=2,nxm1
         ang=(lat(i,j)+lat(i,j+1))*rdpdg
         sin1=sin(ang)
         f=fo*sin1
         do k=1,nz
            ks=1
            nvv=(nv(i,j,k)+nv(i,j+1,k)
     .          +nv(i-1,j+1,k)+nv(i-1,j,k))*.25
            nuu=(nu(i,j,k)+nu(i+1,j,k)
     .          +nu(i,j-1,k)+nu(i+1,j-1,k))*.25  
            fvv=(fv(i,j,k)+fv(i,j+1,k)
     .          +fv(i-1,j+1,k)+fv(i-1,j,k))*.25
            fuu=(fu(i,j,k)+fu(i+1,j,k)
     .          +fu(i,j-1,k)+fu(i+1,j-1,k))*.25  
            wb(i,j,k)=0.
            wc(i,j,k)=0.
            if (u(i,j,k) .ne. bnd)
     .          wc(i,j,k)=f*u(i,j,k)+g/dy(i,j)
     .                   *(phi(i,j+1,k)-phi(i,j,k))+nvv-fvv
            if (v(i,j,k) .ne. bnd)
     .          wb(i,j,k)=-f*v(i,j,k)+g/dx(i,j)
     .                   *(phi(i+1,j,k)-phi(i,j,k))+nuu-fuu
         enddo
      enddo
      enddo
      do k=1,nzm1
      do j=2,nym1
      do i=2,nxm1
         if (ps(i,j) .gt. p(k+1)) then !compute rms geostrophic departure
c                                 on the non-vert staggered std grid
            wa(i,j,k)=sqrt((wc(i-1,j-1,k)+wc(i,j-1,k))**2*.25
     .                    +(wb(i-1,j,k)+wb(i-1,j-1,k))**2*.25)
            sum=wa(i,j,k)**2+sum
            cnt=cnt+1.
         endif
      enddo
      enddo
      enddo

      deallocate (wb,wc)

c commented JS 01-16-01
c     write(2) wa,delo
      
      errms=sqrt(sum/cnt)
      write(6,1000) errms,delo
1000  format(1x,'BEFORE/AFTER BALCON...MOMENTUM RESIDUAL FOR DOMAIN'
     &     ,e12.4,' delo= ',e12.4)
      return
      end
      subroutine diagnose(a,nx,ny,nz,ii,jj,kk,ispan,title)
      implicit none
      integer nx,ny,nz
      real a(nx,ny,nz)
      integer ii,jj,kk
      integer ispan,ieast,iwest,jnorth,jsouth
      integer i,j,k
      character*(*) title
      print*,title,' at level ',kk
c     write(9,1002) title,' at level ',kk
 1002 format(1x,a25,a10,i3)
      jnorth=jj+ispan/2+1
      if(jnorth.gt.ny) jnorth=ny
      jsouth=jj-ispan/2+1
      if(jsouth.lt.1) jsouth=1
      iwest=ii-ispan/2+1
      if(iwest.lt.1) ieast=1
      ieast=ii+ispan/2+1
      if(ieast.gt.nx)ieast=nx 
      do j=jnorth,jsouth,-1
c      write(9,1001)  j,iwest,ieast
 1001  format(1x,3i8)
       print*,j,iwest,ieast
       write(6,1000) (a(i,j,kk),i=iwest,ieast)
c      write(9,1000) (a(i,j,kk),i=iwest,ieast)
 1000  format(1x,7e11.5)
      enddo
      return
      end
c     
c===============================================================================
c
c
      subroutine balcon(to,uo,vo,omo,t,u,v,om,tb,ub,vb,omb,tmp,
     .   rod,delo,tau,itmax,err,erru,errph,errub,errphb
c    .,nu,nv,fu,fv
     .,nx,ny,nz,lat,dx,dy,ps,p,dp,lmax,erf_in)
c
c *** Balcon executes the mass/wind balance computations as described
c        mcginley (Meteor and Atmos Phys, 1987) except that
c        this scheme operates on perturbations from background "b"
c        fields. The dynamic constraint is formulated from this 
c        perturbation field. The constraint equation is
c        partial  du'/dt= -ro*nonlin'-d phi'/dx +fv' + friction' 
c        delo determines the magnitude of the residual and is based on
c        scaled non linear term U**2/length* rossby number
c        Both o and b fields arrive staggered .  A perturbation is 
c        computed prior to the dynamic balance.  A new
c        geopotential(t) is computed using relaxation on eqn. (2).  New
c        u, v and omega winds are computed using eqns. (4), (5) and (6)
c        with the new geopotential and neglecting the lagrange multiplier
c        term.  Next the        
c        lagrange multiplier is computed using 3-d relaxation on eqn. (3).
c        U, v and omega are adjusted by adding the lagrange multiplier
c        term with the new lagrange multiplier.
c        the u,v,t are the balanced output arrays
c        The unique aspect of this analysis is that background model error is
c        specified explicitly over the entire grid as determined from 
c        verification stats. The observed error is an array that takes into
c        account both observation error and interpolation error.
c        omo is the cloud consistent vertical motion
c        Change added 2/23/07 to perform mass continity initially with 
c        no terrain (terbnd is not applied) then balance. The final
c        continuity application then zeros out winds on terrain faces
c        so the final adjustment accounts for terrain. This change was
c        applied owing to spurious mass adjustments owing to large shear
c        on terrain faces and large adjustments to phi and thus temperature.
c
      implicit none
c
      integer   nx,ny,nz
     .         ,nxm1,nym1,nzm1
     .         ,itmax,lmax,ks,kf,ittr
     .         ,i,j,k,l,is,ip,js,jp,kp,it,itt
     .         ,icnt,iwpt,istatus
     .         ,ucnt,vcnt,uwpt,vwpt
c
      real t(nx,ny,nz),to(nx,ny,nz),tb(nx,ny,nz)
     .      ,u(nx,ny,nz),uo(nx,ny,nz),ub(nx,ny,nz)
     .      ,v(nx,ny,nz),vo(nx,ny,nz),vb(nx,ny,nz)
     .      ,om(nx,ny,nz),omo(nx,ny,nz),omb(nx,ny,nz)
     .      ,tmp(nx,ny,nz)
c    .,nu(nx,ny,nz),nv(nx,ny,nz),fu(nx,ny,nz),fv(nx,ny,nz)
     .      ,dx(nx,ny),dy(nx,ny),dp(nz)
     .      ,ps(nx,ny),p(nz),pso(nx,ny)
     .      ,lat(nx,ny),ff(nx,ny)
     .      ,erru(nx,ny,nz),errph(nx,ny,nz)
     .      ,errub(nx,ny,nz),errphb(nx,ny,nz)

      real err,rdpdg,bnd,g,fo,r,re,ovr
     .      ,ang,f,cotmax,sin1,fs,cos1,beta
     .      ,a,bb,cortmt,rod
     .      ,dudy,dvdx,dnudx,dnvdy,tt,uot,vot,tot
     .      ,dt2dx2,dt2dy2,slap,force,rest,cot
     .      ,cotma1,cotm5,rho,cotm0,erf_in,erf,dtdx,dtdy,nuu,nvv
     .      ,dldp,dldx,dldy,tsum,r_missing_data
     .      ,usum,vsum,delo,fuu,fvv
     .      ,angu,angv,dyu,dxv

      real term1,term2,term3,term4,term5,term6
     .      ,term7,term8,term9,term10,term11
     .      ,euoay,euoax,snv,snu,dfvdy,dfudx
     .      ,eueub,fob,foax,foay
     .      ,fu2,fv2,fuangu,fvangv,dt

c 2d array now (JS 2-20-01)
      real tau(nx,ny)

c these are used for diagnostics
      integer nf
      parameter (nf=10)
      real  data(nf),fmax,sum
      real  fldmax(nf,nz),fldmin(nf,nz)
      integer fldmxi(nf,nz),fldmxj(nf,nz),ismx,jsmx,ksmx
      integer fldmni(nf,nz),fldmnj(nf,nz)

      integer   itstatus
      integer   init_timer
      integer   ishow_timer

      logical larray_diag/.false./

      real, allocatable, dimension(:,:,:) :: aaa,bbb
      real, allocatable, dimension(:,:,:) :: fu,fv,nu,nv
      real, allocatable, dimension(:,:) :: dxx,dx2,dxs
      real, allocatable, dimension(:,:) :: dyy,dy2,dys
      real, allocatable, dimension(:,:) :: fx,ffx
      real, allocatable, dimension(:,:) :: fy,ffy

c_______________________________________________________________________________
c
      print *,'balcon'
      nxm1=nx-1
      nym1=ny-1
      nzm1=nz-1
      ks=1
      kf=nz
      rdpdg=3.141592654/180.
      bnd=1.e-30
      g=9.80665
      fo=14.52e-5   !2*omega
      r=287.053
      re=6371220.
      ovr=1.        !over-relaxation factor
      ittr=0

      allocate (aaa(nx,ny,nz),bbb(nx,ny,nz))
c
c only need these calculations once
c
c set pso to bottom pressure level
      do j=1,ny
      do i=1,nx
       pso(i,j)=p(1)
      enddo
      enddo 
      do k=1,nz
       do j=1,ny
        js=1
        if(j.eq.ny) js=0
        do i=1,nx
         is=1
         if(i.eq.nx) is=0
         ang=(lat(i+is,j+js)+lat(i,j+js)+lat(i,j)+lat(i+is,j))*rdpdg*.25
         ff(i,j)=fo*sin(ang)
c wind error must be defined at the geopotential stagger points
         aaa(i,j,k)=erru(i,j,k)+errub(i,j,k)+ff(i,j)**2*delo
         bbb(i,j,k)=(erru(i,j,k)+errub(i,j,k))/aaa(i,j,k)
        enddo
       enddo
      enddo

      allocate (dxx(nx,ny),dx2(nx,ny),dxs(nx,ny),
     1          dyy(nx,ny),dy2(nx,ny),dys(nx,ny),
     1 fx(nx,ny),fy(nx,ny),ffx(nx,ny),ffy(nx,ny))

      do j=2,nym1
      do i=2,nxm1
         dxx(i,j)=(dx(i,j)+dx(i,j+1)+dx(i+1,j)+dx(i+1,j+1))
     &*.25
         dx2(i,j)=dxx(i,j)*2.
         dxs(i,j)=dxx(i,j)*dxx(i,j)
         dyy(i,j)=(dy(i,j)+dy(i,j+1)+dy(i+1,j)+dy(i+1,j+1))
     &*.25
         dy2(i,j)=dyy(i,j)*2.
         dys(i,j)=dyy(i,j)*dyy(i,j)
         fx(i,j)=(ff(i+1,j)-ff(i-1,j))/dx2(i,j)
         fy(i,j)=(ff(i,j+1)-ff(i,j-1))/dy2(i,j)
         ffx(i,j)=ff(i,j)*fx(i,j)
         ffy(i,j)=ff(i,j)*fy(i,j)
      enddo
      enddo

c analz with input fields prior to balcon iterations on lmax
      allocate (fu(nx,ny,nz),fv(nx,ny,nz)
     .         ,nu(nx,ny,nz),nv(nx,ny,nz))


c create perturbations
c owing to an artifact of coding the t array is phi
      print*,'initialize timer'
      itstatus=init_timer()
      print*,'-------------------------'

      do l=1,lmax

       itstatus=ishow_timer()

       write(6,*) '|||||||||BALCON ITERATION NUMBER ',l,' ||||||||||'
       print*,'-----------------------------------------------------'
c      write(9,*) '|||||||||BALCON ITERATION NUMBER ',l,' ||||||||||'
c  set convergence error for relaxation of lamda...set erf to desired 
c  accuracy of wind
       erf=erf_in     
       erf=erf*dx(nx/2,ny/2)
c apply continuity to input winds over complete domain with no terrain
       call leib_sub(nx,ny,nz,erf,tau,erru
     .,lat,dx,dy,pso,p,dp,uo,u,vo,v,
     . omo,om,omb,l,lmax)
c
      do k=1,nz/2    
      print*, 'After continuitlevel ',k
      call diagnose(uo,nx,ny,nz,268,87,k,7,'OBSERVED U ')
      call diagnose(vo,nx,ny,nz,268,87,k,7,'OBSERVED V ')
      call diagnose(u,nx,ny,nz,268,87,k,7,'CONTINTY U ')
      call diagnose(v,nx,ny,nz,268,87,k,7,'CONTINTY V ')
      enddo

      itstatus=ishow_timer()

       if(delo.eq.0.) go to 111
c move adjusted fields to observation fields 
       call move_3d(u,uo,nx,ny,nz)
       call move_3d(v,vo,nx,ny,nz)
       call move_3d(om,omo,nx,ny,nz)
c convert input observed fields to perturbations 
       do j=1,ny
        do i=1,nx
         do k=1,nz
          to(i,j,k)=to(i,j,k)-tb(i,j,k)!background and obs are in GPM
          t(i,j,k)=0. !put background grid in solution grid to establish
c                      boundary values(zero perturbation)
          if(ub(i,j,k).ne.bnd)then
             uo(i,j,k)=uo(i,j,k)-ub(i,j,k)
             u(i,j,k)=0.0
          endif
          if(vb(i,j,k).ne.bnd)then
             vo(i,j,k)=vo(i,j,k)-vb(i,j,k)
             v(i,j,k)=0.0
          endif
         enddo
        enddo
       enddo 

      itstatus=ishow_timer()

      call frict(fu,fv,ub,vb,uo,vo,p,ps,tmp
     .                 ,nx,ny,nz,dx,dy,dp,dt,bnd)
      call nonlin(nu,nv,uo,vo,ub,vb,om,omb
     .           ,nx,ny,nz,dx,dy,dp,dt,bnd,rod)
       do k=1,nz/2
       print*, 'Nonlinear terms at lvl ',k, ' at center of  Dennis' 
       call diagnose(nu,nx,ny,nz,nx/2,ny/2,k,7,'Nonlinear U terms')
       call diagnose(nv,nx,ny,nz,nx/2,ny/2,k,7,'Nonlinear V terms')
       enddo

c *** Compute new phi (t array) using relaxation on eqn. (2).
c        beta*dldx term is dropped to eliminate coupling with lambda eqn.
c
       itstatus=ishow_timer()

       fmax=-1.e30
       do 1 it=1,itmax

         cotmax=0.
         sum=0.
         do 2 k=ks,kf
         do 2 j=2,nym1
          do 2 i=2,nxm1

cbeka           do 2 k=ks,kf

               if(k.ne.kf)then
               if(pso(i,j).lt.p(k+1)) then 
                term1=0.
                term2=0.
                term3=0.
                force=0.
                goto 25
               endif
               endif
               eueub=erru(i,j,k)+errub(i,j,k)
               fob=ff(i,j)/bbb(i,j,k) 
               foax=(ff(i+1,j)/aaa(i+1,j,k)-ff(i-1,j)/aaa(i-1,j,k))
     &              /dx2(i,j)
               foay=(ff(i,j+1)/aaa(i,j+1,k)-ff(i,j-1)/aaa(i,j-1,k))
     &              /dy2(i,j)
               euoax=(erru(i+1,j,k)/aaa(i+1,j,k)-erru(i-1,j,k)/
     &                       aaa(i-1,j,k))/dx2(i,j)
               euoay=(erru(i,j+1,k)/aaa(i,j+1,k)-erru(i,j-1,k)/
     &                       aaa(i,j-1,k))/dy2(i,j)
               term1=delo*(-fob*foax-ffx(i,j)/eueub)
               term2=delo*(-fob*foay-ffy(i,j)/eueub)
               term3=-(errph(i,j,k)+errphb(i,j,k))/(bbb(i,j,k)*delo)
               term4=-errph(i,j,k)/(bbb(i,j,k)*delo)
               term5=-fob*euoay-fy(i,j)*erru(i,j,k)/eueub
               term6= fob*euoax+fx(i,j)*erru(i,j,k)/eueub
               term7=ff(i,j)*erru(i,j,k)/eueub
               term8=-term1                       !delo*(fob*foax+ffx/eueub)
               term9=-term2                       !delo*(fob*foay+ffy/eueub)
               term10=term1
               term11=term2
               dudy=(uo(i,j,k)-uo(i,j-1,k))/dyy(i,j)
               dvdx=(vo(i,j,k)-vo(i-1,j,k))/dxx(i,j)
               dfudx=(fu(i+1,j,k)-fu(i-1,j,k)
     .               +fu(i+1,j-1,k)-fu(i-1,j-1,k))*.5/dx2(i,j)
               dfvdy=(fv(i,j+1,k)-fv(i,j-1,k)
     .               +fv(i-1,j+1,k)-fv(i-1,j-1,k))*.5/dy2(i,j)
               dnudx=(nu(i+1,j,k)-nu(i-1,j,k)
     .               +nu(i+1,j-1,k)-nu(i-1,j-1,k))*.5/dx2(i,j)
               dnvdy=(nv(i,j+1,k)-nv(i,j-1,k)
     .               +nv(i-1,j+1,k)-nv(i-1,j-1,k))*.5/dy2(i,j)
               snv=(nv(i-1,j,k)+nv(i,j,k))*.5
               snu=(nu(i,j-1,k)+nu(i,j,k))*.5
               fuu=(fu(i-1,j,k)+fu(i,j,k))*.5
               fvv=(fv(i-1,j,k)+fv(i,j,k))*.5
               dtdy=(t(i,j+1,k)-t(i,j-1,k))/dy2(i,j)
               dtdx=(t(i+1,j,k)-t(i-1,j,k))/dx2(i,j)
               tot=to(i,j,k)
               uot=(uo(i-1,j,k)+uo(i-1,j-1,k))*.5
               vot=(vo(i,j-1,k)+vo(i-1,j-1,k))*.5
               force=term4*tot+term5*uot+term6*vot+term7*(dvdx-dudy)
     &         +(term8*snu+term9*snv)+term10*fuu+term11*fvv
     &         -(dnudx+dnvdy)+dfudx+dfvdy
25              tt=t(i,j,k)
               dt2dx2=(t(i+1,j,k)+t(i-1,j,k)-2.*tt)/dxs(i,j)
               dt2dy2=(t(i,j+1,k)+t(i,j-1,k)-tt*2.)/dys(i,j)
               slap=dt2dx2+dt2dy2+dtdx*term1+dtdy*term2+term3*tt
               rest=slap-force
               cortmt=-2./dxs(i,j)-2./dys(i,j)+term3
               cot=rest/cortmt
               t(i,j,k)=tt-cot*ovr
               cotmax=amax1(cotmax,abs(cot))
               if (it .eq. 1) cotma1=cotmax
               sum=sum+abs(force)
               if(abs(force).gt.fmax) then
                ismx=i
                jsmx=j
                ksmx=k
                fmax=force
               endif
2          continue
c        write(6,1000) it,cotmax,ovr,cotma1
         ittr=ittr+1
         cotm5=cotmax
c
c ****** Recompute over-relaxation factor every fifth iteration.
c
         if (ittr .ne. 5) goto 15
         ittr=0
         rho=(cotm5/cotm0)**.2
         if (rho .le. 1) ovr=2./(1.+sqrt(1.-rho))
         cotm0=cotm5
15       if (cotmax .lt. err) goto 12
         if (it .ne. 1) goto 1
         cotm0=cotm5
1000     format(1x,'PHI SOLVER: it = ',i4,' max correction '
     &          ,' = ',e12.3
     .         ,'ovr =  ',e12.4/1x
     .         ,'first iteration max correction was ',e12.4)

1      continue  !this is the itmax loop


       itstatus=ishow_timer()
       print*,' ---------------------------------------------'
       print*,'Elapsed time after itmax loop (sec): ',itstatus
       print*,' ---------------------------------------------'
12     write(6,1000) it,cotmax,ovr,cotma1
       print*, 'forcing function max '
       print*, 'max ',fmax, ' at ijk ',ismx,jsmx,ksmx
       print*, 'forcing fcn abs average ', sum/float(it*(kf-ks+1)*
     &       (nym1-1)*(nxm1-1)) 

c     write(9,1000) it,cotmax,ovr,cotma1
c     erf=0.
c
c *** Compute new u, v, omega using eqns. (4), (5), (6) with new phi and
c        without the lagrange multiplier terms.
c
       call initmxmn(nf,nz
     &,fldmax,fldmin,fldmxi,fldmxj,fldmni,fldmnj)

       itstatus=ishow_timer()

       do k=1,nz
        do j=1,ny-1
         js = 1
         if (j .eq.1) js=0
         do i=1,nx-1
          is = 1
          if (i.eq.1) is=0
          angu=(lat(i+1,j)+lat(i,j))*.5*rdpdg
          angv=(lat(i,j+1)+lat(i,j))*.5*rdpdg
          dyu=(dy(i+1,j)+dy(i,j))*.5
          dxv=(dx(i+1,j)+dx(i,j))*.5
          fuangu=fo*sin(angu)
          fvangv=fo*sin(angv)
          dtdx=(t(i+1,j,k)-t(i,j,k))/dxv
          dtdy=(t(i,j+1,k)-t(i,j,k))/dyu
          uot=uo(i,j,k)
          vot=vo(i,j,k)
          fvv=(fv(i,j,k)+fv(i,j+1,k)+fv(i-is,j+1,k)+fv(i-is,j,k))*.25
          fuu=(fu(i,j,k)+fu(i+1,j,k)+fu(i+1,j-js,k)+fu(i,j-js,k))*.25
          nvv=(nv(i,j,k)+nv(i,j+1,k)+nv(i-is,j+1,k)+nv(i-is,j,k))*.25
          nuu=(nu(i,j,k)+nu(i+1,j,k)+nu(i+1,j-js,k)+nu(i,j-js,k))*.25
          if (uot .ne. bnd) then
           fu2=(fuangu*fuangu)
           u(i,j,k)=(uot*erru(i,j,k)
     &     -fuangu*delo*(nvv-fvv+dtdy))/aaa(i,j,k)      

           data(1)=nvv
           data(2)=fvv
           data(3)=dtdy
           data(4)=erru(i,j,k)

           call savemxmninfo(1,4,nf,nz,data,i,j,k
     &,fldmax,fldmin,fldmxi,fldmxj,fldmni,fldmnj)

          else 
           u(i,j,k)=bnd
          endif
          if ( vot .ne. bnd) then
           fv2=(fvangv*fvangv)
           v(i,j,k)=(vot*erru(i,j,k)
     &     +fvangv*delo*(nuu-fuu+dtdx))/aaa(i,j,k)       

           data(5)=nuu
           data(6)=fuu
           data(7)=dtdx
           data(8)=erru(i,j,k)
           call savemxmninfo(5,8,nf,nz,data,i,j,k
     &,fldmax,fldmin,fldmxi,fldmxj,fldmni,fldmnj)

          else
           v(i,j,k)=bnd
          endif
         enddo !on i
        enddo ! on j
c   don't change wind on west (u) or south (v) bndry
          do j=1,ny
           u(1,j,k)=uo(1,j,k)
          enddo
          do i=1,nx
           v(i,1,k)=vo(i,1,k)
          enddo
        if(larray_diag)then
           print*
           print*,'Calling array diagnosis: ',k,p(k)
           print*,'-------------------------------'
           call array_diagnosis(u(1,1,k),nx,ny,' u-comp   ')
           call array_diagnosis(v(1,1,k),nx,ny,' v-comp   ')
c          call array_diagnosis(uo(1,1,k),nx,ny,' uo comp  ')
c          call array_diagnosis(vo(1,1,k),nx,ny,' vo comp  ')
           call array_diagnosis(t(1,1,k),nx,ny,'  phi     ')
        endif

       enddo ! on k

       if(larray_diag)then
          call printmxmn(8,nf,nz,p,fldmax,fldmin
     &               ,fldmxi,fldmxj,fldmni,fldmnj)
       endif
c 

       if(larray_diag)then
        do k=1,nz
          print*
          print*,'Calling perturb  array diagnosis: ',k,p(k)
          print*,'-------------------------------'
          call array_diagnosis(u(1,1,k),nx,ny,' u-comp   ')
          call array_diagnosis(v(1,1,k),nx,ny,' v-comp   ')
          call array_diagnosis(t(1,1,k),nx,ny,'  phi     ')
          call array_diagnosis(om(1,1,k),nx,ny,' omega    ')
c          call array_diagnosis(vo(1,1,k),nx,ny,' vo comp  ')
c          call array_diagnosis(uo(1,1,k),nx,ny,' uo comp  ')
        enddo
       endif


c Restore full winds and heights by adding back in background
       do k=1,nz
          do j=1,ny
          do i=1,nx
           if(ub(i,j,k).ne.bnd) u(i,j,k)=u(i,j,k)+ub(i,j,k)
           if(ub(i,j,k).ne.bnd) uo(i,j,k)=uo(i,j,k)+ub(i,j,k)
           if(vb(i,j,k).ne.bnd) v(i,j,k)=v(i,j,k)+vb(i,j,k)
           if(vb(i,j,k).ne.bnd) vo(i,j,k)=vo(i,j,k)+vb(i,j,k)
           t(i,j,k)=t(i,j,k)+tb(i,j,k)
           to(i,j,k)=to(i,j,k)+tb(i,j,k)
         enddo
         enddo
       enddo


       do k=1,nz/2    
       print*, 'After balance ',k
       call diagnose(uo,nx,ny,nz,268,87,k,7,'OBSERVED U ')
       call diagnose(vo,nx,ny,nz,268,87,k,7,'OBSERVED V ')
       call diagnose(u,nx,ny,nz,268,87,k,7,'BALANCED U ')
       call diagnose(v,nx,ny,nz,268,87,k,7,'BALANCED V ')
       call diagnose(to,nx,ny,nz,268,87,k,7,'OBSV HEIGHT')
       call diagnose(t,nx,ny,nz,268,87,k,7,'ADJUSTED HT')
       enddo

       itstatus=ishow_timer()
       print*,'------------------------------------------'
       print*,'Elapsed time (after leib_sub) sec: ',itstatus
       print*,'------------------------------------------'


      enddo ! on lmax

c apply continuity to final winds...done with backgrounds use ub,vb,omb
c at this point introduc terrain and adjust winds over and around terrain
c note that these winds will not be balanced
!     call terbnd(u,v,om,nx,ny,nz,ps,p,bnd)
!     call terbnd(ub,vb,omb,nx,ny,nz,ps,p,bnd)
       call leib_sub(nx,ny,nz,erf,tau,erru
     .,lat,dx,dy,ps,p,dp,u,ub,v,vb,
     . om,omb,omb,l,lmax)
c evaluate dynamic balance and continuity
      call analzo(t,to,ub,uo,vb,vo,omb,omo
     .                ,nu,nv,fu,fv,delo,tau
     .                ,nx,ny,nz
     .                ,lat,dx,dy,ps,p,dp,l,lmax)
c move adjusted fields (ub,vb,omb) to solution fields 
       call move_3d(ub,u,nx,ny,nz)
       call move_3d(vb,v,nx,ny,nz)
       call move_3d(omb,om,nx,ny,nz)
 111  print*,'------------------------------------------------'
      itstatus=ishow_timer()
      print*,'Elapsed time end of balcon loop (sec): ',itstatus
      print*,'------------------------------------------------'

      deallocate (aaa,bbb,fu,fv,nu,nv)
      deallocate (dxx,dx2,dxs,dyy,dy2,dys
     1,fx,fy,ffx,ffy)


      return
      end
c
c ----------------------------------------------------------------
c
      subroutine leib_sub(nx,ny,nz,erf,tau,erru
     .,lat,dx,dy,ps,p,dp,uo,u,vo,v,
     . omo,om,omb,l,lmax)

      implicit none

      integer nx,ny,nz
      integer nxm1,nym1,nzm1
      integer i,j,k,ks,l,lmax
      integer itstatus,ishow_timer

      real      u(nx,ny,nz),uo(nx,ny,nz)
     .      ,v(nx,ny,nz),vo(nx,ny,nz)
     .      ,om(nx,ny,nz),omo(nx,ny,nz),omb(nx,ny,nz)
     .      ,erru(nx,ny,nz)
     .      ,lat(nx,ny),dx(nx,ny),dy(nx,ny)
     .      ,ps(nx,ny),p(nz),dp(nz)

      real dldx,dldy,dldp,sum,sum1,cnt
     .,a,erf,bnd
      real tau(nx,ny)

      real, allocatable, dimension(:,:,:) :: slam,f3,h

      bnd=1.e-30

c
c *** Compute lagrange multiplier (slam) using 3-d relaxtion on eqn. (3).
c
      nxm1=nx-1
      nym1=ny-1
      nzm1=nz-1

      allocate (slam(nx+1,ny+1,nz+1),f3(nx+1,ny+1,nz+1)
     .      ,h(nx+1,ny+1,nz+1))

      call zero3d(slam,nx+1,ny+1,nz+1)
      call zero3d(h,nx+1,ny+1,nz+1)
c
c ****** Compute a/tau (h) term and rhs terms in eqn. (3)
c
      call fthree(f3,uo,vo,omo,omb,h,erru,tau,
     .   nx,ny,nz,lat,dx,dy,dp)

      itstatus=ishow_timer()
c
c ****** Perform 3-d relaxation.
c
      itstatus=ishow_timer()

      call leibp3(slam,f3,200,erf,h
     .  ,nx,ny,nz,dx,dy,ps,p,dp)
 
      itstatus=ishow_timer()
c
c ****** Compute new u, v, omega by adding the lagrange multiplier terms.
c 
      sum=0
      sum1=0
      cnt=0
      do k=1,nz
      ks=1
      if(k.eq.nz) ks=0
      do j=1,nym1
      do i=1,nxm1

         a=2.*erru(i,j,k)
         dldp=(slam(i,j,k)-slam(i,j,k+1))/dp(k+ks)
         dldx=(slam(i+1,j+1,k+1)-slam(i,j+1,k+1))/dx(i,j)
         dldy=(slam(i+1,j+1,k+1)-slam(i+1,j,k+1))/dy(i,j)
         if (u(i,j,k) .ne. bnd) u(i,j,k)=uo(i,j,k)+dldx/a
         if (v(i,j,k) .ne. bnd) v(i,j,k)=vo(i,j,k)+dldy/a
         if (omo(i,j,k).ne.bnd) om(i,j,k)=omo(i,j,k)+
     &     .5*dldp/tau(i,j)
       sum=sum+(u(i,j,k)-uo(i,j,k))**2+(v(i,j,k)-vo(i,j,k))**2
       sum1=sum1+(om(i,j,k)-omo(i,j,k))**2
       cnt=cnt+1.
      enddo
      enddo
      do i=1,nx
         a=2.*erru(i,ny,k)
         dldp=(slam(i,ny,k)-slam(i,ny,k+1))/dp(k+ks)
         dldy=(slam(i+1,ny+1,k+1)-slam(i+1,ny,k+1))/dy(i,ny)
         if (v(i,ny,k) .ne. bnd) v(i,ny,k)=vo(i,ny,k)+dldy/a
         if (omo(i,ny,k).ne.bnd) om(i,ny,k)=omo(i,ny,k)+
     &     .5*dldp/tau(i,ny)
      enddo
      do j=1,ny
         a=2.*erru(nx,j,k)
         dldp=(slam(nx,j,k)-slam(nx,j,k+1))/dp(k+ks)
         dldx=(slam(nx+1,j+1,k+1)-slam(nx,j+1,k+1))/dx(nx,j)
         if (u(nx,j,k) .ne. bnd) u(nx,j,k)=uo(nx,j,k)+dldx/a
         if (omo(nx,j,k).ne.bnd) om(nx,j,k)=omo(nx,j,k)+
     &     .5*dldp/tau(nx,j)
      enddo
      enddo
c   print out rms vector adjustment
      print*, 'RMS Vector(m/s) adjustment after continuity applied ' 
      print*, sqrt(sum/cnt)
      print*, 'RMS omega (Pa/s)adjustment after continuity applied ' 
      print*, sqrt(sum1/cnt)

      deallocate (slam,f3,h)

      return
      end
c
c ---------------------------------------------------------------
c
      subroutine analzo(t,to,u,uo,v,vo,om,omo,nu,nv,fu,fv,
     .               delo,tau ,nx,ny,nz,lat,dx,dy,ps,p,dp,l,lmax)
c
c     analzo is a diagnostic routine that looks at geostrophic residual
c     maxima, continuity residual maxs.  It computes the rms terms in the 
c     variational formalism

      implicit none
c
      integer   nx,ny,nz
     .         ,nzm1,nxm2,nym2
     .         ,isv,jsv,ksv
     .         ,i,j,k,iflag
     .         ,l,lmax,istatus
c
      real t(nx,ny,nz),to(nx,ny,nz)
     .      ,u(nx,ny,nz),uo(nx,ny,nz)
     .      ,v(nx,ny,nz),vo(nx,ny,nz)
     .      ,om(nx,ny,nz),omo(nx,ny,nz)
     .      ,nu(nx,ny,nz),nv(nx,ny,nz)
     .      ,fu(nx,ny,nz),fv(nx,ny,nz)
     .      ,lat(nx,ny),dx(nx,ny),dy(nx,ny)
     .      ,ps(nx,ny),p(nz),dp(nz)
     .      ,delo,rdpdg,fo,bnd
     .      ,conmax,sumom,cont,con,sumu,sumv,sumt,thermu,thermv
     .      ,sumww,resu,resv,tgpu,tgpv,tgpph,tgpc,sumf,sumn
     .      ,ang,sin1,f,nvv,nuu,dtdx,dtdy,dudx,dvdy,domdp
     .      ,uot,vot,tot,uuu,vvv,contm,fuu,fvv,dxx,dyy
     .      ,ures,vres

      real tau(nx,ny)
      real r_missing_data

c_______________________________________________________________________________
c
c     write(9,*) '******ANALZ OUTPUT BY LAPS LAYER***********'
      write(6,*) '******ANALZ OUTPUT BY LAPS LAYER***********'

      nzm1=nz-1
      nxm2=nx-2
      nym2=ny-2
      rdpdg=3.141592654/180.
      fo=14.52e-5   !2*omega
      bnd=1.e-30

      do k=2,nzm1
         write(6,4000) k,p(k+1)/100.
c        write(9,4000) k,p(k+1)/100.
4000     format(1x,'----------Level ',i4,'   ',f5.0,' mb---------')
         conmax=0
         sumom=0
         cont=0
         sumu=0.
         sumv=0.
         sumt=0.
         thermu=0.
         thermv=0.
         sumww=0.
         resu=0.
         resv=0.
         tgpu=0.
         tgpv=0.
         tgpph=0.
         tgpc=0.
         iflag=0
         sumf=0
         sumn=0

         do j=2,nym2
         do i=2,nxm2

         dxx=(dx(i-1,j-1)+dx(i,j-1)+dx(i-1,j)+dx(i,j))*.25
         dyy=(dy(i-1,j-1)+dy(i,j-1)+dy(i-1,j)+dy(i,j))*.25
         if(ps(i,j).gt.p(k+1)) then

            ang=(lat(i,j)+lat(i+1,j))*.5*rdpdg
            sin1=sin(ang)
            f=fo*sin1
c friction is from perturbation winds only so may not be perfect
c for diagnosis near ground
            fvv=(fv(i,j,k)+fv(i,j+1,k)+fv(i-1,j+1,k)+fv(i-1,j,k))/4.
            fuu=(fu(i,j,k)+fu(i+1,j,k)+fu(i+1,j-1,k)+fu(i,j-1,k))/4.
            nvv=(nv(i,j,k)+nv(i,j+1,k)+nv(i-1,j+1,k)+nv(i-1,j,k))/4.
            nuu=(nu(i,j,k)+nu(i+1,j,k)+nu(i+1,j-1,k)+nu(i,j-1,k))/4.
            dtdx=(t(i+1,j,k)-t(i,j,k))/dx(i,j)
            dtdy=(t(i,j+1,k)-t(i,j,k))/dy(i,j)
            dudx=(u(i,j-1,k-1)-u(i-1,j-1,k-1))/dxx
            dvdy=(v(i-1,j,k-1)-v(i-1,j-1,k-1))/dyy    
            domdp=(om(i,j,k-1)-om(i,j,k))/dp(k)
            uot=uo(i,j,k)
            vot=vo(i,j,k)
            tot=to(i,j,k)

            if (uot .ne. bnd) then
               thermu=(u(i,j,k)+dtdy/f)**2+thermu
               sumu=sumu+(u(i,j,k)-uot)**2
               vres=(f*u(i,j,k)+dtdy+nvv-fvv)**2
               resv=vres+resv
               tgpu=tgpu+1.
            endif
            if (vot .ne. bnd) then
               thermv=(v(i,j,k)-dtdx/f)**2+thermv
               sumv=sumv+(v(i,j,k)-vot)**2
               ures=(-f*v(i,j,k)+dtdx+nuu-fuu)**2
               resu=resu+ures
               tgpv=tgpv+1.
            endif

            sumt=sumt+(t(i,j,k)-tot)**2
            tgpph=tgpph+1.

            if (ps(i,j) .gt. p(k+1)) then
               cont=(dudx+dvdy+domdp)**2+cont
               con=dudx+dvdy+domdp
               sumom=(om(i,j,k)-omo(i,j,k))**2+sumom
               uuu=(u(i,j-1,k-1)+u(i-1,j-1,k-1))*.5
               vvv=(v(i-1,j,k-1)+v(i-1,j-1,k-1))*.5
               sumww=uuu**2+vvv**2+sumww
               sumf=fuu**2+fvv**2+sumf
               sumn=nuu**2+nvv**2+sumn
               tgpc=tgpc+1.
               if (abs(con) .ge. conmax) then
                  conmax=abs(con)
                  contm=con
                  isv=i
                  jsv=j
                  ksv=k
                  iflag=1
               endif
            endif
         endif!on below ground check
         enddo
         enddo

         if (tgpu .ne. 0) sumu=sqrt(sumu/tgpu)
         if (tgpv .ne. 0) sumv=sqrt(sumv/tgpv)
         if (tgpph .ne. 0) sumt=sqrt(sumt/tgpph)
         if (tgpu .ne. 0) sumom=sqrt(sumom/tgpu)
         if (tgpu .ne. 0) thermu=sqrt(thermu/tgpu)
         if (tgpv .ne. 0) thermv=sqrt(thermv/tgpv)
         if (tgpc .ne. 0.) cont=sqrt(cont/tgpc)
         if (tgpc .ne. 0.) sumww=sqrt(sumww/tgpc)
         if (tgpc .ne. 0.) sumf =sqrt(sumf/tgpc)
         if (tgpc .ne. 0.) sumn =sqrt(sumn/tgpc)
         if (tgpu .ne. 0) resu=sqrt(resu/tgpu)
         if (tgpv .ne. 0) resv=sqrt(resv/tgpv)
         write(6,1002) sumww
c        write(9,1002) sumww
1002     format(1x,' rms wind speed:',e12.4)
         if (sumww .ne. 0.) sumww=sqrt(thermu**2+thermv**2)/sumww

c took tau out of this write 2-20-01 (JS)
         write(6,1001) sumt,delo,sumu,sumv,
     .             sumom,tau(1,1),cont,thermu,thermv,sumww,resu,resv,
     .             sumf,sumn

c        write(9,1001) sumt,delo,sumu,sumv,
c    .             sumom,tau(1,1),cont,thermu,thermv,sumww,resu,resv

      enddo

1009  format(1x,'relax ',3i4/1x,'con f3 l ',3e12.4/1x,'lijk ',6e12.4/1x,
     .       ' d x y p ',3e12.4/1x,'uuvvoo ',6e12.4)
1001  format(1x,' rms error for each term in functional and del,tau'
     . /1x,'    phi-phio:',e12.5,'   delo:',e12.5,
     . /1x,'    u-uo:    ',e12.5,'   v-vo:',e12.5,
     . /1x,'    om-omo:  ',e12.5,',  !!!tau:',e12.5,
     . /1x,' rms cont eqn error:   ',e12.5
     . /1x,' rms ageostrophic wind:',2e12.5
     . /1x,' implied rossby number:',f6.3 
     . /1x,' momentum resids-u and v eqn:',2e12.4
     . /1x,' rms friction vector        :', e12.4 
     . /1x,' rms nonlinear term value   :', e12.4)
c Write out sample vertical profiles of u, uo, v, vo, phi phio
c in center of grid
      if(.false.)then
      Print*, ' U component in middle of grid: adjusted, non adjusted'
     &          ,' dif output - input'
      Do k=1,nz
       write(6,1111) u(50,16,k), uo(50,16,k)
     &     ,u(50,16,k)-uo(50,16,k)
      enddo
      Print*, ' V component in middle of grid: adjusted, non adjusted'
     &          ,' dif output - input'
      Do k=1,nz
       write(6,1111) v(50,16,k), vo(50,16,k)
     &     ,v(50,16,k)-vo(50,16,k)
      enddo
      Print*, ' Geopotentialin middle of grid: adjusted, non adjusted'
     &          ,' dif output - input'
      Do k=1,nz
       write(6,1111) t(nx/2,ny/2,k), to(nx/2,ny/2,k)
     &     ,t(nx/2,ny/2,k)-to(nx/2,ny/2,k)
      enddo 
 1111 format (1x, 3f8.1)

      endif
c
      return
      end
c
c===============================================================================
      subroutine analz(t,to,u,uo,v,vo,om,omo,slam,f3,nu,nv,fu,fv,
     .               delo,tau ,nx,ny,nz,lat,dx,dy,ps,p,dp,l,lmax)
c
c     analz is a diagnostic routine that looks at geostrophic residual
c     maxima, continuity residual maxs.  It computes the rms terms in the 
c     variational formalism

      implicit none
c
      integer   nx,ny,nz
     .         ,nzm1,nxm2,nym2
     .         ,isv,jsv,ksv
     .         ,i,j,k,iflag
     .         ,l,lmax,istatus
c
      real t(nx,ny,nz),to(nx,ny,nz)
     .      ,u(nx,ny,nz),uo(nx,ny,nz)
     .      ,v(nx,ny,nz),vo(nx,ny,nz)
     .      ,om(nx,ny,nz),omo(nx,ny,nz)
     .      ,slam(nx+1,ny+1,nz+1),f3(nx+1,ny+1,nz+1)
     .      ,nu(nx,ny,nz),nv(nx,ny,nz)
     .      ,fu(nx,ny,nz),fv(nx,ny,nz)
     .      ,lat(nx,ny),dx(nx,ny),dy(nx,ny)
     .      ,ps(nx,ny),p(nz),dp(nz)
     .      ,delo,rdpdg,fo,bnd
     .      ,conmax,sumom,cont,con,sumu,sumv,sumt,thermu,thermv
     .      ,sumww,resu,resv,tgpu,tgpv,tgpph,tgpc
     .      ,ang,sin1,f,nvv,nuu,dtdx,dtdy,dudx,dvdy,domdp
     .      ,uot,vot,tot,uuu,vvv,contm,fuu,fvv,dxx,dyy
     .      ,ures,vres

      real tau(nx,ny)
      real r_missing_data

c_______________________________________________________________________________
c
c     write(9,*) '******ANALZ OUTPUT BY LAPS LAYER***********'
      write(6,*) '******ANALZ OUTPUT BY LAPS LAYER***********'

      nzm1=nz-1
      nxm2=nx-2
      nym2=ny-2
      rdpdg=3.141592654/180.
      fo=14.52e-5   !2*omega
      bnd=1.e-30

      do k=2,nzm1
         write(6,4000) k,p(k+1)/100.
c        write(9,4000) k,p(k+1)/100.
4000     format(1x,'----------Level ',i4,'   ',f5.0,' mb---------')
         conmax=0
         sumom=0
         cont=0
         sumu=0.
         sumv=0.
         sumt=0.
         thermu=0.
         thermv=0.
         sumww=0.
         resu=0.
         resv=0.
         tgpu=0.
         tgpv=0.
         tgpph=0.
         tgpc=0.
         iflag=0

         do j=2,nym2
         do i=2,nxm2

         dxx=(dx(i-1,j-1)+dx(i,j-1)+dx(i-1,j)+dx(i,j))*.25
         dyy=(dy(i-1,j-1)+dy(i,j-1)+dy(i-1,j)+dy(i,j))*.25
         if(ps(i,j).gt.p(k)) then

            ang=(lat(i,j))*rdpdg
            sin1=sin(ang)
            f=fo*sin1
            fvv=(fv(i,j,k)+fv(i,j+1,k)+fv(i-1,j+1,k)+fv(i-1,j,k))/4.
            fuu=(fu(i,j,k)+fu(i+1,j,k)+fu(i+1,j-1,k)+fu(i,j-1,k))/4.
            nvv=(nv(i,j,k)+nv(i,j+1,k)+nv(i-1,j+1,k)+nv(i-1,j,k))/4.
            nuu=(nu(i,j,k)+nu(i+1,j,k)+nu(i+1,j-1,k)+nu(i,j-1,k))/4.
            dtdx=(t(i+1,j,k)-t(i,j,k))/dx(i,j)
            dtdy=(t(i,j+1,k)-t(i,j,k))/dy(i,j)
            dudx=(u(i,j-1,k-1)-u(i-1,j-1,k-1))/dxx
            dvdy=(v(i-1,j,k-1)-v(i-1,j-1,k-1))/dyy    
            domdp=(om(i,j,k-1)-om(i,j,k))/dp(k)
            uot=uo(i,j,k)
            vot=vo(i,j,k)
            tot=to(i,j,k)

            if (uot .ne. bnd) then
               thermu=(u(i,j,k)+dtdy/f)**2+thermu
               sumu=sumu+(u(i,j,k)-uot)**2
               vres=(f*u(i,j,k)+dtdy+nvv-fvv)**2
               resv=vres+resv
               tgpu=tgpu+1.
            endif
            if (vot .ne. bnd) then
               thermv=(v(i,j,k)-dtdx/f)**2+thermv
               sumv=sumv+(v(i,j,k)-vot)**2
               ures=(-f*v(i,j,k)+dtdx+nuu-fuu)**2
               resu=resu+ures
               tgpv=tgpv+1.
            endif

            sumt=sumt+(t(i,j,k)-tot)**2
            tgpph=tgpph+1.

            if (ps(i,j) .gt. p(k)) then
               cont=(dudx+dvdy+domdp)**2+cont
               con=dudx+dvdy+domdp
               sumom=(om(i,j,k)-omo(i,j,k))**2+sumom
               uuu=(u(i,j-1,k-1)+u(i-1,j-1,k-1))*.5
               vvv=(v(i-1,j,k-1)+v(i-1,j-1,k-1))*.5
               sumww=uuu**2+vvv**2+sumww
               tgpc=tgpc+1.
               if (abs(con) .ge. conmax) then
                  conmax=abs(con)
                  contm=con
                  isv=i
                  jsv=j
                  ksv=k
                  iflag=1
               endif
            endif
         endif!on below ground check
         enddo
         enddo

         if (tgpu .ne. 0) sumu=sqrt(sumu/tgpu)
         if (tgpv .ne. 0) sumv=sqrt(sumv/tgpv)
         if (tgpph .ne. 0) sumt=sqrt(sumt/tgpph)
         if (tgpu .ne. 0) sumom=sqrt(sumom/tgpu)
         if (tgpu .ne. 0) thermu=sqrt(thermu/tgpu)
         if (tgpv .ne. 0) thermv=sqrt(thermv/tgpv)
         if (tgpc .ne. 0.) cont=sqrt(cont/tgpc)
         if (tgpc .ne. 0.) sumww=sqrt(sumww/tgpc)
         if (tgpu .ne. 0) resu=sqrt(resu/tgpu)
         if (tgpv .ne. 0) resv=sqrt(resv/tgpv)
         write(6,1002) sumww
c        write(9,1002) sumww
1002     format(1x,' rms wind speed:',e12.4)
         if (sumww .ne. 0.) sumww=sqrt(thermu**2+thermv**2)/sumww

c took tau out of this write 2-20-01 (JS)
         write(6,1001) sumt,delo,sumu,sumv,
     .             sumom,tau(1,1),cont,thermu,thermv,sumww,resu,resv

c        write(9,1001) sumt,delo,sumu,sumv,
c    .             sumom,tau(1,1),cont,thermu,thermv,sumww,resu,resv

         if (iflag.eq. 1) write(6,1009)
     .      isv,jsv,ksv,contm,f3(isv,jsv,ksv),slam(isv,jsv,ksv),
     .      slam(isv+1,jsv,ksv),slam(isv-1,jsv,ksv),slam(isv,jsv+1,ksv),
     .      slam(isv,jsv-1,ksv),slam(isv,jsv,ksv+1),slam(isv,jsv,ksv-1),
     .      dx(isv,jsv),dy(isv,jsv),dp(ksv),
     .      u(isv,jsv-1,ksv-1),u(isv-1,jsv-1,ksv-1),v(isv-1,jsv,ksv-1),
     .      v(isv-1,jsv-1,ksv-1),om(isv,jsv,ksv),om(isv,jsv,ksv-1)
c        if (iflag.eq. 1) write(9,1009)
c    .      isv,jsv,ksv,contm,f3(isv,jsv,ksv),slam(isv,jsv,ksv),
c    .      slam(isv+1,jsv,ksv),slam(isv-1,jsv,ksv),slam(isv,jsv+1,ksv),
c    .      slam(isv,jsv-1,ksv),slam(isv,jsv,ksv+1),slam(isv,jsv,ksv-1),
c    .      dx(isv,jsv),dy(isv,jsv),dp(ksv),
c    .      u(isv,jsv-1,ksv-1),u(isv-1,jsv-1,ksv-1),v(isv-1,jsv,ksv-1),
c    .      v(isv-1,jsv-1,ksv-1),om(isv,jsv,ksv),om(isv,jsv,ksv-1)

      enddo

1009  format(1x,'relax ',3i4/1x,'con f3 l ',3e12.4/1x,'lijk ',6e12.4/1x,
     .       ' d x y p ',3e12.4/1x,'uuvvoo ',6e12.4)
1001  format(1x,' rms error for each term in functional and del,tau'
     . /1x,'    phi-phio:',e12.5,'   delo:',e12.5,
     . /1x,'    u-uo:    ',e12.5,'   v-vo:',e12.5,
     . /1x,'    om-omo:  ',e12.5,',  !!!tau:',e12.5,
     . /1x,' rms cont eqn error:   ',e12.5
     . /1x,' rms ageostrophic wind:',2e12.5
     . /1x,' implied rossby number:',f6.3 
     . /1x,' momentum resids-u and v eqn:',2e12.4)
c
      return
      end
c
c===============================================================================
c
      subroutine terbnd(u,v,om,nx,ny,nz,ps,p,bnd)
c  terrain boundary sets velocity component to 1.e-30 on terrain 
c  faces and in the interior of terrain.  
c  The input variables are on the staggered (mcginley1987) grid

      implicit none
c
      integer   nx,ny,nz,nzm1,i,j,k
      integer   nxm1
      integer   nym1
c
      real u(nx,ny,nz),v(nx,ny,nz)
     .      ,om(nx,ny,nz)
     .      ,ps(nx,ny),p(nz)
     .      ,bnd
c_______________________________________________________________________________
c
      nym1=ny-1
      nxm1=nx-1
      nzm1=nz-1
      do j=2,ny
      do i=2,nx
         om(i,j,1)=bnd
         do k=2,nzm1
            if (ps(i,j) .le. p(k)) then ! point is below ground
                u(i-1,j-1,k-1)=bnd
                u(i,j-1,k-1)=bnd
                v(i-1,j-1,k-1)=bnd
                v(i-1,j,k-1)=bnd
                om(i,j,k)=bnd
                om(i,j,k-1)=bnd
            endif
         enddo
      enddo
      enddo
c
      do k=2,nz
         do j=2,ny
            if(ps(1,j) .le. p(k)) then
               u(1,j-1,k-1)=bnd   
               om(1,j,k)=bnd
               om(1,j,k-1)=bnd
            endif
         enddo

         do i=2,nx
            if(ps(i,1) .le. p(k)) then
               v(i-1,1,k-1)=bnd
               om(i,1,k)=bnd
               om(i,1,k-1)=bnd
            endif
         enddo
      enddo
      return
      end
c
c===============================================================================
c
      subroutine frict(fu,fv,ub,vb,u,v,p,ps,t
     .                 ,nx,ny,nz,dx,dy,dp,dt,bnd)
c
c *** Frict computes frictional acceleration in the first layer
c     near the terrain surface
c     using the parameterization from Haltner Ch. 10.
c     Frict is computed from the perturbed u, v grids
c
      implicit none
c
      integer   nx,ny,nz
     .         ,nxm1,nym1,nzm1
     .         ,i,j,k
     .         ,imx,jmx,imn,jmn
c
      real fu(nx,ny,nz),fv(nx,ny,nz)
     .      ,ub(nx,ny,nz),vb(nx,ny,nz) 
     .      ,u(nx,ny,nz),v(nx,ny,nz),t(nx,ny,nz),bnd  
     .      ,dx(nx,ny),dy(nx,ny),dp(nz),ps(nx,ny),p(nz)
     .      ,ddp(nz),ddp2(nz)
     .      ,dt,dudx,dvdy,dudy,dvdx
     .      ,dudp,dvdp,dudt,dvdt
     .      ,rmx2d,rmn2d,dps

      real gdtvdp_save(nx,ny),gdtudp_save(nx,ny)
      real gdtvdp,gdtudp

      real   daudp,davdp
     .      ,daubdp,davbdp,dubdp,dvbdp
     .      ,d2vbdx,d2vbdy,d2ubdx,d2ubdy
     .      ,d2vdx,d2vdy,d2udx,d2udy,dauvdp,dauvbdp
     .      ,d2udp,d2vdp,d2ubdp,d2vbdp,absuv,absuvb
     .      ,hzdf,eddf
     .      ,grav,r,den,cd

      real, allocatable, dimension(:,:) :: dx2,dy2
c_______________________________________________________________________________
c
c
      print *,'frictn'
c compute difference in frictional component between rawcomponent and 
c background. Fu-Fub, Fv-Fvb
c constnts
      grav=9.808!m/sec2
      r=287.04!gas const
c     hzdf=50000.!m2/sec ! this value in haltiner is suspected to be too large
      hzdf=500.!m2/sec
      eddf=200.!pa2/sec
      cd=.0025 !non dim
      nxm1=nx-1
      nym1=ny-1
      nzm1=nz-1
      call zero3d(fu,nx,ny,nz)
      call zero3d(fv,nx,ny,nz)

      allocate (dx2(nx,ny))
      allocate (dy2(nx,ny))

      do j=1,ny
      do i=1,nx
         dx2(i,j)=dx(i,j)*dx(i,j)
         dy2(i,j)=dy(i,j)*dy(i,j)
         gdtudp_save(i,j)=0.
         gdtvdp_save(i,j)=0.
      enddo
      enddo

      do k=2,nzm1
         ddp(k)=dp(k)+dp(k+1)
         ddp2(k)=ddp(k)*ddp(k)
      enddo

      do k=2,nzm1
      do j=2,nym1
      do i=2,nxm1
         dauvdp=(sqrt(u(i,j,k-1)**2+v(i,j,k-1)**2)-
     &           sqrt(u(i,j,k+1)**2+v(i,j,k+1)**2))/ddp(k)
         absuv=sqrt(u(i,j,k)**2+v(i,j,k)**2)
         if(u(i,j,k).eq.bnd) then
         else
          d2udx=(u(i+1,j,k)+u(i-1,j,k)-2.*u(i,j,k))/dx2(i,j)
          d2udy=(u(i,j+1,k)+u(i,j-1,k)-2.*u(i,j,k))/dy2(i,j)
          d2udp=4.*(u(i,j,k+1)+u(i,j,k-1)-2.*u(i,j,k))/ddp2(k)
          gdtudp=0.
          if(ps(i,j).gt.p(k+1).and.(ps(i,j)-p(k+1)).le.dp(k))then! we're near sfc
           den=p(k+1)/r/t(i,j,k+1)
c          surface friction term is always scaled for a contact layer of 
c          5000. pascals. Smaller dps will give values too high, 
           dps=5000.
           gdtudp=den*cd*grav*abs(u(i,j,k))*u(i,j,k)/dps         
           gdtudp_save(i,j)=gdtudp
          endif
          fu(i,j,k)=hzdf*(d2udx+d2udy)
     &         +eddf*d2udp-gdtudp
         endif
         if(v(i,j,k).eq.bnd) then
         else
          d2vdx=(v(i+1,j,k)+v(i-1,j,k)-2.*v(i,j,k))/dx2(i,j)
          d2vdy=(v(i,j+1,k)+v(i,j-1,k)-2.*v(i,j,k))/dy2(i,j)
          d2vdp=4.*(v(i,j,k+1)+v(i,j,k-1)-2.*v(i,j,k))/ddp2(k)
          gdtvdp=0.
          if(ps(i,j).gt.p(k+1).and.(ps(i,j)-p(k+1)).le.dp(k)) then! we're near sfc
           den=p(k+1)/r/t(i,j,k+1)
           gdtvdp=den*cd*grav*abs(v(i,j,k))*v(i,j,k)/dp(k)      
           gdtvdp_save(i,j)=gdtvdp
          endif
          fv(i,j,k)=hzdf*(d2vdx+d2vdy)
     &              +eddf*d2vdp-gdtvdp
         endif
      enddo
      enddo
      enddo
c fill friction arrays
      do k=2,nzm1
      do j=1,ny
       fu(1,j,k)=fu(2,j,k)
       fu(nx,j,k)=fu(nxm1,j,k)
       fv(1,j,k)=fv(2,j,k)
       fv(nx,j,k)=fv(nxm1,j,k)
      enddo
      do i=1,nx
       fu(i,ny,k)=fu(i,nym1,k)
       fu(i,1,k)=fu(i,2,k)
       fv(i,ny,k)=fv(i,nym1,k)
       fv(i,1,k)=fv(i,2,k)
      enddo
      enddo
c
      deallocate (dx2)
      deallocate (dy2)

      print*,'gdtvdp and gdtudp diagnostics'
      print*,'-----------------------------'

      call get_mxmn_2d(nx,ny,gdtudp_save,rmx2d,rmn2d
     &,imx,jmx,imn,jmn)
      print*,'max/min gdtudp i/j = ',rmx2d,imx,jmx,
     &rmn2d,imn,jmn

      if(rmx2d.gt. 1.0)then
         print*,'WARNING: max gdtupd > 1.0'
      endif

      call get_mxmn_2d(nx,ny,gdtvdp_save,rmx2d,rmn2d
     &,imx,jmx,imn,jmn)
      print*,'max/min gdtvdp i/j = ',rmx2d,imx,jmx,
     &rmn2d,imn,jmn

      if(rmx2d .gt. 1.0)then
         print*,'WARNING: max gdtvpd > 1.0'
      endif



      return
      end
c===============================================================================
c
      subroutine nonlin(nu,nv,u,v,ub,vb,om,omb
     .                 ,nx,ny,nz,dx,dy,dp,dt,bnd,rod)
c
c *** Nonlin computes the non-linear terms (nu,nv) from staggered input.
c     The non linear terms are linearized with a background/perturbation
c     combination. The eularian time terms are assumed to be zero.
c     The non-linear terms nu,nv, are computed on the 
c     u, v grids, respectively. The nonlin term is scaled by the rossby number
c     rod - the length scale is that defined by data resolvability 
c
      implicit none
c
      integer   nx,ny,nz
     .         ,nxm1,nym1,nzm1
     .         ,i,j,k
c
      real nu(nx,ny,nz),nv(nx,ny,nz)
     .      ,u(nx,ny,nz),v(nx,ny,nz)     !time=t
     .      ,ub(nx,ny,nz),vb(nx,ny,nz)     !time=t
     .      ,om(nx,ny,nz),omb(nx,ny,nz)
     .      ,dx(nx,ny),dy(nx,ny),dp(nz),bnd
     .      ,dt,dudx,dvdy,dudy,dvdx
     .      ,dudp,dvdp,dudt,dvdt
     .      ,dvbdy,dubdp,dvbdp,dudby
     .      ,dubdx,dvbdx,dubdy
     .      ,omu,ombu,omv,ombv,vvu,uuv,vvbu,uubv,rod
c_______________________________________________________________________________
c
c
      print *,'nonlin'
      nxm1=nx-1
      nym1=ny-1
      nzm1=nz-1
      call zero3d(nu,nx,ny,nz)
      call zero3d(nv,nx,ny,nz)
c the u,v,om variables are perturbations
      do k=2,nzm1
      do j=2,nym1
      do i=2,nxm1
         if(u(i,j,k).eq.bnd) then
          nu(i,j,k)=0.
         else
          ombu=(omb(i+1,j+1,k+1)+omb(i,j+1,k+1)+omb(i+1,j+1,k)+
     &            omb(i,j+1,k))*.25
          omu=(om(i+1,j+1,k+1)+om(i,j+1,k+1)+om(i+1,j+1,k)+
     &            om(i,j+1,k))*.25
          vvu=(v(i-1,j+1,k)+v(i,j+1,k)+v(i-1,j,k)+v(i,j,k))*.25
          vvbu=(vb(i-1,j+1,k)+vb(i,j+1,k)+vb(i-1,j,k)+vb(i,j,k))*.25
          dubdx=(ub(i+1,j,k)-ub(i-1,j,k))/(2.*dx(i,j))
          dudx=(u(i+1,j,k)-u(i-1,j,k))/(2.*dx(i,j))
          dubdy=(ub(i,j+1,k)-ub(i,j-1,k))/(2.*dy(i,j))
          dudy=(u(i,j+1,k)-u(i,j-1,k))/(2.*dy(i,j))
          dubdp=(ub(i,j,k+1)-ub(i,j,k-1))/(dp(k)+dp(k+1))
          dudp=(u(i,j,k+1)-u(i,j,k-1))/(dp(k)+dp(k+1))
          dudt=0.
          nu(i,j,k)=(dudt+ub(i,j,k)*dudx+vvbu*dudy+ombu*dudp
     &                 +u(i,j,k)*dubdx+vvu*dubdy+omu*dubdp)*rod
c         nu(i,j,k)=(ub(i,j,k)+u(i,j,k))*(dudx+dubdx)+(vvu+vvbu)*
c    &              (dudy+dubdy)+ombu*dudp+omu*dubdp 
         endif
         if(v(i,j,k).eq.bnd) then
          nv(i,j,k)=0.
         else 
          ombv=(omb(i+1,j+1,k+1)+omb(i+1,j,k+1)+omb(i+1,j,k)+
     &            omb(i+1,j+1,k))*.25
          omv=(om(i+1,j+1,k+1)+om(i+1,j,k+1)+om(i+1,j,k)+
     &            om(i+1,j+1,k))*.25
          uuv=(u(i,j,k)+u(i+1,j,k)+u(i,j-1,k)+u(i+1,j-1,k))*.25
          uubv=(ub(i,j,k)+ub(i+1,j,k)+ub(i,j-1,k)+ub(i+1,j-1,k))*.25
          dvbdx=(vb(i+1,j,k)-vb(i-1,j,k))/(2.*dx(i,j))
          dvdx=(v(i+1,j,k)-v(i-1,j,k))/(2.*dx(i,j))
          dvbdy=(vb(i,j+1,k)-vb(i,j-1,k))/(2.*dy(i,j))
          dvdy=(v(i,j+1,k)-v(i,j-1,k))/(2.*dy(i,j))
          dvbdp=(vb(i,j,k+1)-vb(i,j,k-1))/(dp(k)+dp(k+1))
          dvdp=(v(i,j,k+1)-v(i,j,k-1))/(dp(k)+dp(k+1))
          dvdt=0.           
          nv(i,j,k)=(dvdt+uubv*dvdx+vb(i,j,k)*dvdy+ombv*dvdp
     &                 +uuv*dvbdx+v(i,j,k)*dvbdy+omv*dvbdp)*rod
c         nv(i,j,k)=(uubv+uuv)*(dvdx+dvbdx)+(vb(i,j,k)+v(i,j,k))*
c    &                 *(dvdy+dvbdy)+ombv*dvdp+omv*dvbdp
         endif
      enddo
      enddo
      enddo
c fill nonlinear arrays near boundaries, don't want big gradients.
      do k=2,nzm1
      do j=1,ny
       nu(1,j,k)=nu(2,j,k)
       nu(nx,j,k)=nu(nxm1,j,k)
       nv(1,j,k)=nv(2,j,k)
       nv(nx,j,k)=nv(nxm1,j,k)
      enddo
      do i=1,nx
       nu(i,ny,k)=nu(i,nym1,k)
       nu(i,1,k)=nu(i,2,k)
       nv(i,ny,k)=nv(i,nym1,k)
       nv(i,1,k)=nv(i,2,k)
      enddo
      enddo
c
      return
      end
c
c===============================================================================
c
      subroutine fthree(f3,u,v,om,omb,h,erru,tau
     .,nx,ny,nz,lat,dx,dy,dp)
c
c *** Fthree computes a/tau (h) and rhs terms in eqn. (3).
c
      implicit none
c
      integer   nx,ny,nz,nzp1,nxp1,nyp1
     .         ,i,j,k,is,js,ks
c
      real f3(nx+1,ny+1,nz+1),h(nx+1,ny+1,nz+1)
     .      ,erru(nx,ny,nz)
     .      ,u(nx,ny,nz),v(nx,ny,nz)
     .      ,om(nx,ny,nz),lat(nx,ny)
     .      ,omb(nx,ny,nz)
     .      ,dx(nx,ny),dy(nx,ny),dp(nz)
     .      ,dpp,aa,dxx,dyy
     .      ,cont,formax
      real tau(nx,ny)
c_______________________________________________________________________________
c
      print *,'fthree'
      formax=0.
      do k=2,nz
         dpp=dp(k)
         do j=2,ny
         do i=2,nx
            h(i,j,k)=erru(i,j,k)/tau(i,j)
            f3(i,j,k)=-2.*erru(i,j,k)*
     &    ((u(i,j-1,k-1)-u(i-1,j-1,k-1))/dx(i,j)
     .    +(v(i-1,j,k-1)-v(i-1,j-1,k-1))/dy(i,j)
     .    +(om(i,j,k-1)-om(i,j,k))/dpp)
           if (abs(f3(i,j,k)) .ge. formax) then
               formax=abs(f3(i,j,k))
               is=i
               js=j
               ks=k
            endif
         enddo
         enddo
      enddo
      print*, 'f3: Maximum forcing function for lamda ',formax
      print*, 'at ',is,js,ks                                  
c
      return
      end
c     
c===============================================================================
c
      subroutine leibp3(sol,force,itmax,erf,h
     .                 ,nx,ny,nz,dx,dy,ps,p,dp)
c
c *** Leibp3 performs 3-d relaxation.
c
      implicit none
c
      integer   nx,ny,nz
     .         ,nxm1,nym1,nzm1
     .         ,ittr,is,js
     .         ,i,j,k,it,itmax,ia
c
      real sol(nx+1,ny+1,nz+1)
     .      ,force(nx+1,ny+1,nz+1)
     .      ,h(nx+1,ny+1,nz+1)
     .      ,dx(nx,ny),dy(nx,ny)
     .      ,ps(nx,ny),p(nz),dp(nz)
     .      ,erf,si,sj,sk
     .      ,ovr,erb,hh,ertm,ermm,corlm
     .      ,dx2,dx1s,dy2,dy1s,dz,dz2,dz1s
     .      ,aa,cortm,res,cor,corb,corlmm
     .      ,reslm,rho,cor0,cor5

      logical ltest_3d(nx,ny,nz)
c_______________________________________________________________________________
c
c *** Relaxer solver...eqn must be.......                                 
c        sxx+syy+h*szz-force=0   
c
      print *,'start subroutine leibp3'

!     Setup ltest_3d array to help minimize if testing within do loops
      ltest_3d = .false.
      do k=2,nz
        do j=2,ny
          js=1
          if(j.eq.ny)js=0
          do i=2,nx
            is=1
            if(i.eq.nx) is=0
            if (ps(i,j).lt.p(k))then
              ltest_3d(i,j,k) = .true.
            else
              if(ps(i+is,j).lt.p(k)) ltest_3d(i,j,k) = .true.
              if(ps(i-1 ,j).lt.p(k)) ltest_3d(i,j,k) = .true.
              if(ps(i,j+js).lt.p(k)) ltest_3d(i,j,k) = .true.
              if(ps(i,j-1).lt.p(k)) ltest_3d(i,j,k) = .true.
              if(ps(i,j).lt.p(k-1)) ltest_3d(i,j,k) = .true.
            endif
          enddo
        enddo
      enddo

      ovr=1.0
      erb=0.
      si=1.
      sj=1.
      sk=1.
c
c *** First guess here.
c
      nxm1=nx-1
      nym1=ny-1
      nzm1=nz-1
      ittr=0.
      hh=0.
      do it=1,itmax
         ertm=0.
         ermm=0.
         ia=0
         corlm=0.
cbeka changing the loops order
         do k=2,nz
         dz=dp(k)
         dz1s=dz*dz

         do j=2,ny
         js=1
         if(j.eq.ny)js=0
         do i=2,nx
            is=1
            if(i.eq.nx) is=0
            dx1s=dx(i,j)*dx(i,j)
            dy1s=dy(i,j)*dy(i,j)
cbeka           do 2 k=2,nz
cbeka            dz=dp(k)
cbeka            dz1s=dz*dz
            if(ltest_3d(i,j,k))then
               if (ps(i,j).lt.p(k)) go to 2 
               if(ps(i+is,j).lt.p(k)) sol(i+1,j,k)=sol(i,j,k)
               if(ps(i-1 ,j).lt.p(k)) sol(i-1,j,k)=sol(i,j,k)
               if(ps(i,j+js).lt.p(k)) sol(i,j+1,k)=sol(i,j,k)
               if(ps(i,j-1).lt.p(k)) sol(i,j-1,k)=sol(i,j,k)
               if(ps(i,j).lt.p(k-1)) sol(i,j,k-1)=sol(i,j,k)
            endif
            aa=h(i,j,k)
            cortm=-2./dx1s-2./dy1s-aa*2./dz1s
            res= (sol(i+1,j,k)+sol(i-1,j,k))/dx1s+
     .           (sol(i,j+1,k)+sol(i,j-1,k))/dy1s+
     .           aa*(sol(i,j,k+1)+sol(i,j,k-1))/dz1s
     .           +cortm*sol(i,j,k)-force(i,j,k)
            cor=res/cortm
            corb=5.*erf+1.
            if (it .ne. 1) corb=cor/corlmm
            if (abs(corb) .gt. erf) ia=1
            if (abs(cor) .gt. corlm) corlm=abs(cor)
            sol(i,j,k)=sol(i,j,k)-cor*ovr
2            continue
         enddo
         enddo
	enddo
         reslm=corlm*cortm
c         write(6,1001) it,reslm,corlm,corlmm,erb
1200     format(1x,3e12.5)
         erb=amax1(ermm,ertm)
5        ittr=ittr+1
         cor5=corlm
         if (ittr .eq. 5) then
            ittr=0
            rho=(cor5/cor0)**.2
            if (rho .le. 1) ovr=2./(1.+sqrt(1.-rho))
            cor0=cor5
         endif
         if (ia .eq. 1 .and. it .eq. 1) then
            corlmm=corlm
            cor0=corlmm
         endif
         if(corlm.lt.erf) go to 20
      enddo! on it

20    reslm=corlm*cortm

      write(6,1001) it,reslm,corlm,corlmm,erb
      write(6,1002) ovr
c     write(9,1001) it,reslm,corlm,corlmm,erb
c     write(9,1002) ovr
1002       format(1x,'LEIBP3:ovr rlxtn const at fnl ittr = ',e10.4)
1001       format(1x,'iterations= ',i4,' max residual= ',e10.3,
     .    ' max correction= ',e10.3, ' first iter max cor= ',e10.3,
     .    ' max bndry error= ',e10.3)
c
      return
      end
c                           
c===============================================================================
c
      subroutine leib(sol,force,itmax,erf,nx,ny,nz
     .               ,ps,p,a,b,c,d,e,dx,dy,dz,lpress)
c
c *** Leib performs 2-d relaxation.
c
      implicit none
c
      integer   nx,ny,nz
     .         ,nxm1,nym1,nzm1
     .         ,lpress,ittr,itmax,ia
     .         ,i,j,k,kk,it
c
      real sol(nx,ny,nz),force(nx,ny,nz)
     .      ,ps(nx,ny),p(nz)
     .      ,a(nx,ny),b(nx,ny),c(nx,ny)
     .      ,d(nx,ny),e(nx,ny)
     .      ,dx(nx,ny),dy(nx,ny),dz(nz)
     .      ,ovr,reslmm,erf,erb,ertm,ermm,corlm,cortm,res
     .      ,cor,cor5,cor0,corlmm,reslm,rho
     .      ,dx2,dxs,dy2,dys,dz2,dzs,aa,bb,cc,dd
c_______________________________________________________________________________
c
      print *,'leib'
      ovr=1.
      reslmm=0.
      erb=0.
c
c *** First guess here. 
c
      do 6 k=1,nz
      do 6 j=1,ny
      do 6 i=1,nx
         sol(i,j,k)=0.
6     continue
      nxm1=nx-1
      nym1=ny-1
      nzm1=nz-1
      ittr=0.
      do it=1,itmax
         ertm=0.
         ermm=0.
         ia=0
         corlm=0.
         do j=2,nym1
         do i=2,nxm1
            dy2=dy(i,j)*2.
            dys=dy(i,j)*dy(i,j)
            dx2=dx(i,j)*2.
            dxs=dx(i,j)*dx(i,j)
            aa=a(i,j)
            bb=b(i,j)
            cc=c(i,j)
            dd=d(i,j)
            do 2 kk=2,nzm1
               k=nz+1-kk
               dz2=dz(k)*2.
               dzs=dz(k)*dz(k)
               cortm=-2./dxs-2./dys-aa*2./dzs+e(i,j)
c              if(lpress) 20,20,21
c21               sol(i,j,nz)=sol(i,j,nzm1)
c              if(ps(i,j)-p(k)) 25,25,20
c25            if(lpress.eq.1) call bound(sol,i,j,k)
c              goto 2
20             res=(sol(i+1,j,k)+sol(i-1,j,k))/dxs+
     .             (sol(i,j+1,k)+sol(i,j-1,k))/dys+
     .             (sol(i,j,k+1)+sol(i,j,k-1))/dzs*
     .              aa+cortm*sol(i,j,k)-force(i,j,k)+
     .              bb*(sol(i+1,j,k)-sol(i-1,j,k))/dx2+
     .              cc*(sol(i,j+1,k)-sol(i,j-1,k))/dy2+
     .              dd*(sol(i,j,k+1)-sol(i,j,k-1))/dz2
               cor=res/cortm
               if(abs(cor).gt.erf) ia=1
               if(abs(cor).gt.corlm) corlm=abs(cor)
               sol(i,j,k)=sol(i,j,k)-cor*ovr
2           continue
         enddo
         enddo
5        ittr=ittr+1
         cor5=corlm
         if(ittr.eq.5) then
            ittr=0
            rho=(cor5/cor0)**.2
            if(rho.le.1) ovr=2./(1.+sqrt(1.-rho))
            cor0=cor5
         endif
         if (ia .ne. 1) goto 4
         if (it .eq. 1) then
            corlmm=corlm
            cor0=corlmm
         endif
      enddo
4     reslm=corlm*cortm
      write(6,1001) it,reslm,corlm,corlmm,erb
      write(6,1002) ovr
1002       format(1x,'ovr rlxtn const at fnl ittr = ',e10.4)
1001       format(1x,'iterations= ',i4,' max residual= ',e10.3
     .      ,' max correction= ',e10.3, ' first iter max cor= '
     .      ,e10.3,'max bndry error= ',e10.3)
c
      return
      end
c
c===============================================================================
c
      Subroutine balstagger(u,v,phi,t,sh,om,
     &  us,vs,phis,ts,shs,oms,
     &nx,ny,nz,p,ps,idstag)
c
c  This routine takes a standard LAPS field with all variables at
c  each grid point and produces the E-stagger appropriate for applying
c  the dynamic balancing in qbalpe.f   or vice versa
c  idstag > 0 staggers (LAPS -> stagger)
c  idstag < 0 destaggers (stagger -> LAPS). For this latter process
c  the non staggered input grids must be the original gridded laps fields
c  to ensure that the boundaries are consistent with lga backgrounds.  The 
c  destaggering will only process interior grid points on the laps mesh

      real u(nx,ny,nz),v(nx,ny,nz),om(nx,ny,nz),t(nx,ny,nz),  
     &us(nx,ny,nz),vs(nx,ny,nz),oms(nx,ny,nz),ts(nx,ny,nz),  
     &phi(nx,ny,nz),ps(nx,ny)
     &,p(nz),phis(nx,ny,nz),
     &shs(nx,ny,nz),sh(nx,ny,nz)

      real w   ! Mixing ratio, added by B. Shaw, Sep 01
      real, allocatable, dimension(:,:,:,:) :: wr

c wind on terrain face
      bnd=1.e-30
c  gas const
      r=287.04
c  gravity
      g=9.80665
      if(idstag.eq.0) then
         print*, 'no staggering accomplished'
         return
      endif
      if(idstag.gt.0) then

        allocate (wr(nx+1,ny+1,nz,6))

c set vertical stagger first
c first wind level for balcon is second level in LAPS
c omega is shifted one-half in vertical
c t is shifted likewise
c phi is shifted upward one level like winds
c level nz will be the same as laps for all fields
       do j=1,ny
        do i=1,nx
         do k=1,nz-1
          wr(i,j,k,1)=u(i,j,k+1) 
          wr(i,j,k,2)=v(i,j,k+1)
          wr(i,j,k,3)=(om(i,j,k)+om(i,j,k+1))*.5
          wr(i,j,k,4)=(t(i,j,k)+t(i,j,k+1))*.5
          wr(i,j,k,5)=(sh(i,j,k)+sh(i,j,k+1))*.5
          wr(i,j,k,6)=phi(i,j,k+1) 
         enddo
        enddo
       enddo

c horzizontal stagger
c extrapolate to fill north row and east column
c 2nd order taylor series is used
c rest of variables
       do kk=1,6
        do k=1,nz
         do i=1,nx
          wr(i,ny+1,k,kk)=3.*wr(i,ny,k,kk)-3.*wr(i,ny-1,k,kk)
     &                     +wr(i,ny-2,k,kk)
         enddo
         do j=1,ny
          wr(nx+1,j,k,kk)=3.*wr(nx,j,k,kk)-3.*wr(nx-1,j,k,kk)
     &                      +wr(nx-2,j,k,kk)
         enddo
          wr(nx+1,ny+1,k,kk)=3.*wr(nx,ny,k,kk)-3.*wr(nx-1,ny-1,k,kk)
     &                      +wr(nx-2,ny-2,k,kk)
        enddo
       enddo 
c now do horizontal interpolation to staggered grid
       do k=1,nz
        do j=1,ny
         do i=1,nx
          us(i,j,k)=(wr(i+1,j+1,k,1)+wr(i,j+1,k,1))*.5
          vs(i,j,k)=(wr(i+1,j+1,k,2)+wr(i+1,j,k,2))*.5
          ts(i,j,k)=(wr(i,j,k,4)+wr(i+1,j+1,k,4)+wr(i,j+1,k,4)
     %                 +wr(i+1,j,k,4))*.25
          shs(i,j,k)=(wr(i,j,k,5)+wr(i+1,j+1,k,5)+wr(i,j+1,k,5)
     %                 +wr(i+1,j,k,5))*.25
          phis(i,j,k)=(wr(i,j,k,6)+wr(i+1,j+1,k,6)+wr(i,j+1,k,6)
     %                 +wr(i+1,j,k,6))*.25
          oms(i,j,k)=wr(i,j,k,3)
c omega is already on the staggered horizontal mesh
         enddo
        enddo
       enddo
       ts(:,:,nz)=ts(:,:,nz-1)

       deallocate (wr)
       
       return

      else! de-stagger

c horizontal destagger
c re compute hydrostatic virtual t from staggered adjusted phis 
c use a third-order scheme. use local dlnp as a constant for
c scheme for each layer. This avoids precision errors for 
c terms with dlnp**3. Also 2nd deriv terms vanish
       nxx=nx/2
       nyy=ny/2
       do k=2,nz-1
        alpka1k=alog(p(k+1)/p(k))
        t1=1./24./alpka1k 
        do j=1,ny
         do i=1,nx
          for=(phis(i,j,k)-phis(i,j,k-1))/
     &                              alog(p(k)/p(k+1))
          if(k+2.le.nz.and.k.gt.2) then
!          tor= t1*(phis(i,j,k+1)-3.*phis(i,j,k)+3.*phis(i,j,k-1)-
!    &            phis(i,j,k-2))
           tor=0.
          else
           tor=0.
          endif
         ts(i,j,k)=(for-tor)/r          
          if(i.eq.nxx.and.j.eq.nyy)
     &          print*, 'i,j,k, for/r, tor/r temp ',i,j,k,for/r,
     &                  -tor/r,ts(i,j,k)
         enddo
        enddo
       enddo

       tdd=-50.!lowest possible temp for esw calculation

c put background values in stagged temp arrays at top and btm.
       do j=1,ny
        do i=1,nx
          ts(i,j,1)=t(i,j,1)
          ts(i,j,nz)=t(i,j,nz)
        enddo
       enddo
       do k=2,nz
        do j=2,ny
         do i=2,nx
          u(i,j,k)=(us(i,j-1,k-1)+us(i-1,j-1,k-1))*.5
          v(i,j,k)=(vs(i-1,j,k-1)+vs(i-1,j-1,k-1))*.5
          om(i,j,k)=(oms(i,j,k)+oms(i,j,k-1))*.5
          sh(i,j,k)=(shs(i-1,j-1,k-1)+shs(i,j-1,k-1)+
     &               shs(i-1,j,k-1)+shs(i,j,k-1)+
     &               shs(i-1,j-1,k)+shs(i,j-1,k)+
     &               shs(i-1,j,k-1)+shs(i,j,k))*.125
          phi(i,j,k)=(phis(i-1,j-1,k-1)+phis(i,j-1,k-1)+
     &               phis(i-1,j,k-1)+phis(i,j,k-1))*.25
          t(i,j,k)=(ts(i-1,j-1,k-1)+ts(i,j-1,k-1)+
     &               ts(i-1,j,k-1)+ts(i,j,k-1)+
     &               ts(i-1,j-1,k)+ts(i,j-1,k)+
     &               ts(i-1,j,k)+ts(i,j,k))*.125
c  devirtualize temperature

c  New method simply inverts the conventional virtual
c  temperature formula (Tv=T*(1.+0.61w), where
c  w is the mixing ratio) to get T from Tv
c  using the specific humidity.

          w = sh(i,j,k)/(1.-sh(i,j,k))
          t(i,j,k) = t(i,j,k)/(1.+0.61*w)

c         tdm=dwpt(t(i,j,k)-273.16,sh(i,j,k))
c         tvkm=tv(t(i,j,k)-273.16,tdm,p(k)/100.)
c         tvkd=tv(t(i,j,k)-273.16,tdd,p(k)/100.)
c         t(i,j,k)=t(i,j,k)-(tvkm-tvkd)
         enddo
        enddo
       enddo
c final step is bring temps in line with heights on the west and
c south boundaries: use taylor extrapolation
c north and east boundaries for omega
       do k=1,nz
        do j=2,ny
         if(ps(1,j).ge.p(k)) then

          call destagger_y(1,nx,ny,nz,u,j,k,bnd,p)
          call destagger_y(1,nx,ny,nz,v,j,k,bnd,p)
          call destagger_y(1,nx,ny,nz,om,j,k,bnd,p)

         else

          u(1,j,k)=bnd
          v(1,j,k)=bnd
          om(1,j,k)=bnd

         endif

         if(ps(nx,j).ge.p(k)) then

          call destagger_y(nx,nx,ny,nz,om,j,k,bnd,p)

         else

          om(nx,j,k)=bnd

         endif

         call destagger_y(1,nx,ny,nz,t,j,k,bnd,p)
         call destagger_y(1,nx,ny,nz,phi,j,k,bnd,p)
         call destagger_y(1,nx,ny,nz,sh,j,k,bnd,p)

        enddo

        do i=2,nx
         if(ps(i,1).ge.p(k)) then

          call destagger_x(1,nx,ny,nz,u,i,k,bnd,p)
          call destagger_x(1,nx,ny,nz,v,i,k,bnd,p)
          call destagger_x(1,nx,ny,nz,om,i,k,bnd,p)

         else

          u(i,1,k)=bnd
          v(i,1,k)=bnd
          om(i,1,k)=bnd

         endif

         if(ps(i,ny).ge.p(k)) then

          call destagger_x(ny,nx,ny,nz,om,i,k,bnd,p)

         else

          om(i,ny,k)=bnd

         endif

         call destagger_x(1,nx,ny,nz,t,i,k,bnd,p)
         call destagger_x(1,nx,ny,nz,phi,i,k,bnd,p)
         call destagger_x(1,nx,ny,nz,sh,i,k,bnd,p)

        enddo

        if(ps(1,1).ge.p(k)) then

          call destagger_c(nx,ny,nz,u,k,bnd,p)
          call destagger_c(nx,ny,nz,v,k,bnd,p)
          call destagger_c(nx,ny,nz,om,k,bnd,p)

        else

          u(1,1,k)=bnd
          v(1,1,k)=bnd
          om(1,1,k)=bnd

        endif

        call destagger_c(nx,ny,nz,t,k,bnd,p)
        call destagger_c(nx,ny,nz,phi,k,bnd,p)
        call destagger_c(nx,ny,nz,sh,k,bnd,p)

       enddo 
      endif ! destagger

      return
      end
c
c--------------------------------------------------
c
      subroutine destagger_y(i,nx,ny,nz,data,j,k,bnd,p)

      implicit none
      integer nx,ny,nz
      real data(nx,ny,nz)
      real p(nz)
      real bnd,pt1,pt2
      integer i,j,k,jj
      integer i1,i2,i3

      if(i.eq.nx)then
         i1=i-1
         i2=i-2
         i3=i-3
      else
         i1=i+1
         i2=i+2
         i3=i+3
      endif

      if( (data(i1,j,k).ne.bnd).and.
     1    (data(i2,j,k).ne.bnd).and.
     1    (data(i3,j,k).ne.bnd) ) then

        data(i,j,k)=3.*data(i1,j,k)-3.*data(i2,j,k)+data(i3,j,k)

      elseif( (data(i1,j,k).eq.bnd) .and.
     1(data(i2,j,k).ne.bnd) .and. (data(i3,j,k).ne.bnd))then

        data(i,j,k)=(data(i2,j,k)+data(i3,j,k))*.5

      elseif(j.ne.ny)then

        if(data(i,j+1,k).ne.bnd)then

           data(i,j,k)=data(i,j+1,k)

        endif

      else

        pt1=bnd
        pt2=bnd
        do jj=j-1,1,-1
           if(data(i,jj,k).ne.bnd)then
              pt1=data(i,jj,k)
           endif
        enddo
        do jj=j+1,ny
           if(data(i,jj,k).ne.bnd)then
              pt2=data(i,jj,k)
           endif
        enddo
        if(pt1.ne.bnd.and.pt2.ne.bnd)then
           data(i,j,k)=(pt1+pt2)*.5
        endif

      endif

      if(k.ne.nz.and.data(i,j,k).eq.bnd)then

        data(i,j,k)=data(i,j,k+1)+
     1              data(i,j,k+1)*(alog(p(k+1)/p(k)))
      endif

      return
      end

c
c--------------------------------------------------
c
      subroutine destagger_x(j,nx,ny,nz,data,i,k,bnd,p)

      implicit none
      integer nx,ny,nz
      real data(nx,ny,nz)
      real p(nz)
      real bnd,pt1,pt2
      integer i,j,k,ii
      integer j1,j2,j3

      if(j.eq.ny)then
         j1=j-1
         j2=j-2
         j3=j-3
      else
         j1=j+1
         j2=j+2
         j3=j+3
      endif

      if( (data(i,j1,k).ne.bnd).and.
     1    (data(i,j2,k).ne.bnd).and.
     1    (data(i,j3,k).ne.bnd) ) then

        data(i,j,k)=3.*data(i,j1,k)-3.*data(i,j2,k)+data(i,j3,k)

      elseif( (data(i,j1,k).eq.bnd) .and.
     1(data(i,j2,k).ne.bnd) .and. (data(i,j3,k).ne.bnd))then

        data(i,j,k)=(data(i,j2,k)+data(i,j3,k))*.5

      elseif(data(i-1,j,k).ne.bnd)then

        data(i,j,k)=data(i-1,j,k)

      elseif(i.ne.nx)then

        if(data(i+1,j,k).ne.bnd)then

           data(i,j,k)=data(i+1,j,k)

        endif

      else

        pt1=bnd
        pt2=bnd
        do ii=i-1,1,-1
           if(data(ii,j,k).ne.bnd)then
              pt1=data(ii,j,k)
           endif
        enddo
        do ii=i+1,ny
           if(data(ii,j,k).ne.bnd)then
              pt2=data(ii,j,k)
           endif
        enddo
        if(pt1.ne.bnd.and.pt2.ne.bnd)then
           data(i,j,k)=(pt1+pt2)*.5
        endif

      endif

      if(k.ne.nz.and.data(i,j,k).eq.bnd)then

        data(i,j,k)=data(i,j,k+1)+
     1              data(i,j,k+1)*(alog(p(k+1)/p(k)))
      endif

      return
      end
c
c ------------------------------------------------------------
c
      subroutine destagger_c(nx,ny,nz,data,k,bnd,p)

      implicit none
      integer nx,ny,nz
      real data(nx,ny,nz)
      real p(nz)
      real bnd
      integer k

      if( (data(2,2,k).ne.bnd).and.
     1    (data(3,3,k).ne.bnd).and.
     1    (data(4,4,k).ne.bnd) ) then

        data(1,1,k)=3.*data(2,2,k)-3.*data(3,3,k)+data(4,4,k)

      elseif( (data(2,2,k).eq.bnd) .and.
     1(data(3,3,k).ne.bnd) .and. (data(4,4,k).ne.bnd))then

        data(1,1,k)=(data(3,3,k)+data(4,4,k))*.5

      elseif(data(2,1,k).ne.bnd.and.data(1,2,k).ne.bnd)then

        data(1,1,k)=(data(2,1,k)+data(1,2,k))*.5

      elseif(k.ne.nz)then

        data(1,1,k)=data(1,1,k+1)+
     1              data(1,1,k+1)*(alog(p(k+1)/p(k)))
      endif

      return
      end
c--------------------------------------------------
      subroutine savemxmninfo(is,ie,nf,nz,data,in,jn,kn,
     &fldmx,fldmn,fldmxi,fldmxj,fldmni,fldmnj)

      implicit none

      integer in,jn,kn,is,ie,nf,m,nz
      real    fldmx(nf,nz),fldmn(nf,nz)
      real    data(nf)
      integer fldmxi(nf,nz),fldmxj(nf,nz)
      integer fldmni(nf,nz),fldmnj(nf,nz)

      do m=is,ie
         if(data(m).gt.fldmx(m,kn))then
            fldmx(m,kn)=data(m)
            fldmxi(m,kn)=in
            fldmxj(m,kn)=jn
         endif
         if(data(m).lt.fldmn(m,kn))then
            fldmn(m,kn)=data(m)
            fldmni(m,kn)=in
            fldmnj(m,kn)=jn
         endif
      enddo

      return
      end

      subroutine initmxmn(nf,nz
     &,fldmax,fldmin,fldmxi,fldmxj,fldmni,fldmnj)

      implicit none

      integer nx,ny,nz,k,l,nf

      real fldmax(nf,nz)
      real fldmin(nf,nz)
      integer fldmxi(nf,nz),fldmxj(nf,nz)
      integer fldmni(nf,nz),fldmnj(nf,nz)

      integer istatus
      real    r_missing_data

      call get_r_missing_data(r_missing_data,istatus)
      call get_grid_dim_xy(nx,ny,istatus)

      do k=1,nz
      do l=1,nf

         fldmax(l,k)=-r_missing_data-1.
         fldmin(l,k)=r_missing_data
         fldmni(l,k)=nx+1
         fldmnj(l,k)=ny+1
         fldmxi(l,k)=nx+1
         fldmxj(l,k)=ny+1

      enddo
      enddo
      return
      end

      subroutine printmxmn(n,nf,nz,p,fldmax,fldmin
     &,fldmxi,fldmxj,fldmni,fldmnj)

      implicit none

      integer nf, nz
      integer i,k,n
      real    fldmax(nf,nz)
      real    fldmin(nf,nz)
      real    p(nz)
      integer fldmxi(nf,nz)
      integer fldmxj(nf,nz)
      integer fldmni(nf,nz)
      integer fldmnj(nf,nz)


      print*,'max/min info'

      do k=1,nz
         print*
         print*,'  ========================'
         print*,'  level ',k, 'p = ',p(k)
         print*,'  ========================'
         print*,'/   nvv   /   fvv   /',
     &'    dtdy   /   erru   /   nuu   /',
     &'    fuu    /   dtdx   /   erru  /'
         print*,'-------------------------'
     &,'-------------------------------------'
     &,'-------------------------------------'
         print*,'              Maxima'
         print*,'              ======'
         write(6,10)(fldmax(i,k),i=1,n)
         print*,'              Minima'
         print*,'              ======'
         write(6,10)(fldmin(i,k),i=1,n)
      enddo
10    format(1x,8(e10.3,1x))

      return
      end
cc
c------------------------------------------------------------
cc
      subroutine readprg(a9_time,nx,ny,nz,max_pr
     &,udrop,vdrop,tdrop,ri,rj,rk,rit,rjt,rkt,nsnd,istatus)

c this subroutine reads the .prg and .tmg files to recover observed  
c u, v, bnd T profiles, the decimal i,j locations at the nz LAPS levels 

      real udrop(max_pr,nz)
      real vdrop(max_pr,nz)
      real tdrop(max_pr,nz)
      real ri(max_pr,nz)
      real rj(max_pr,nz)
      real rit(max_pr,nz)
      real rjt(max_pr,nz)
      real rk(max_pr,nz)
      real rkt(max_pr,nz)
      real dum2

c variables output:
c     udrop, vdrop tdrop: dropsonde u,v T obs at the nz laps levels
c     ri, rj: real grid coordinates of dropsone position in grid space
c     rit,rjt: real position of temperature sonde in grid space
c
c Note: arrays udrop, vdrop, tdrop
c       must be 2d to allow more than 1.
c
      Character*180 dum
      Character*9 a9_time
      Character*255 dum1
      integer len,istatus,nstar,nblank,ncnt,iflag,ngood
      logical lexist
      real, allocatable, dimension(:,:) :: dd,ff
     1,tt,uu,vv

c     integer max_pr

      integer ngoodlevs(max_pr)

      allocate(dd(max_pr,nz),ff(max_pr,nz),tt(max_pr,nz))
      allocate(uu(max_pr,nz),vv(max_pr,nz))

      call get_r_missing_data(smsng,istatus)
      if(istatus.ne.1)then
         print*,'Error: returned from get_r_missing_data'
         return
      endif

c preset all dropsonde output to missing
      istatus=0
      ri=smsng
      rj=smsng
      rk=smsng
      rkt=smsng
      rit=smsng
      rjt=smsng
      udrop=smsng
      vdrop=smsng
      tdrop=smsng
      dd=smsng
      ff=smsng
      ngoodlevs=0

      pi=4.*atan(1.)
      rdpdg=pi/180. 
      call get_directory('prg',dum,len)
      dum(len+1:len+13)=a9_time//'.prg'
      inquire(file=dum,exist=lexist)
      if(lexist)then
        open (11, file=dum,form='formatted',status='old',err=50)
        iflag=0
        nlev=1
        nsnd=1
        do n=1,nz*max_pr

         read(11,1000,end=1) riin,rjin,rkin,dum1(1:20)
 1000    format(f11.5,2f10.5,a20)

         ri(nsnd,nlev)=riin
         rj(nsnd,nlev)=rjin
         rk(nsnd,nlev)=rkin

         if(dum1(1:1).eq.'*'.and.iflag.eq.0 ) then  !the first level of first sounding

          ri(nsnd,nlev)=smsng
          rj(nsnd,nlev)=smsng
          rk(nsnd,nlev)=smsng
          nlev=nlev+1

         elseif(dum1(1:1).eq.'*'.and.iflag.eq.1)then !could be the end of existing sounding or
c                                                     missing levels within existing sounding.
          if(nlev.gt.nz)then
             nsnd=nsnd+1
             nlev=1
             ri(nsnd,nlev)=smsng
             rj(nsnd,nlev)=smsng
             rk(nsnd,nlev)=smsng
             nlev=nlev+1
             iflag=0
          else
             ri(nsnd,nlev)=smsng
             rj(nsnd,nlev)=smsng
             rk(nsnd,nlev)=smsng
             if(nlev.eq.nz)then
                nlev=1
                nsnd=nsnd+1
             else
                nlev=nlev+1
             endif
          endif

         else

          read(dum1(4:10),100)dd(nsnd,nlev)
100       format(f7.3)
          read(dum1(15:20),101)ff(nsnd,nlev)
101       format(f6.3)
          ri(nsnd,nlev)=ri(nsnd,nlev)+1
          rj(nsnd,nlev)=rj(nsnd,nlev)+1
          ngoodlevs(nsnd)=ngoodlevs(nsnd)+1
          iflag=1
          nlev=nlev+1
         endif
        enddo

1       if(ngoodlevs(nsnd).eq.0)nsnd=nsnd-1

        if(nsnd.lt.1)then
           print*,'No wind data in prg file'
           istatus=-1
           goto 7
        else
           print*,'Found ',nsnd,' Soundings'
           print*,'------------------------'
           do i=1,nsnd
             print*,' Snd#: ',i,' Number of Good Levels = ',ngoodlevs(i)
           enddo
        endif

c convert dd ff to u,v
        do n=1,nsnd 
         do k=1,nz
           if(dd(n,k).ne.smsng)then
           call disp_to_uv(dd(n,k),ff(n,k),udrop(n,k),vdrop(n,k))
           endif
         enddo
        enddo

c if dropsonde is reverse order, flip it.
c       if(rk(1).gt.rk(nsave))then
c          call flip_sonde(mxz,ncnt ,uu,vv,ri,rj,rk)
c       endif
   44   close(11)
      else !on lexist
       print*,'No prg file at this time'
       istatus=-1
      endif

7     continue
c now read the tmg file to get the dropsonde temps
c     call get_directory('tmg',dum,len)
c     dum(len+1:len+13)=a9_time//'.tmg'
c     nn=0
c     inquire(file=dum,exist=lexist)
c     if(lexist)then
c      continue
c      open (11, file=dum,form='formatted',status='old',err=50) 
c      Do n=1,mxz
c       nn=nn+1
c       read(11,*,end=2) aa,bb,cc,ee, dum1           
c       if(dum1.eq.'  ')  go to 2
c       if(dum1.eq.'ACA') then    
c        rit(nn)=aa
c        rjt(nn)=bb
c        rkt(nn)=cc
c        tt(nn)=ee  
c       endif
c      enddo
c since dropsonde is reverse order, flip it.
c 2    nn=nn-1
c      if(nn.lt.1)then
c         istatus = istatus-2
c         deallocate(ri,rj,rk,dd,ff,uu,vv,tt)
c         goto 49
c      endif
c      if(nn.eq.1) then! there is only one ob; put at nearest laps level
c        do k=1,nz! clear arrays
c         tdrop(k)=smsng
c         rit(k)=smsng
c         rjt(k)=smsng
c        enddo
c        k=rk(1)
c        tdrop(k)=tt(1)
c        rit(k)=ri(1)
c        rjt(k)=rj(1) 
c        close (11) ! close file  and return
c        return
c      endif
c      if(rk(1).gt.rk(nn))then
c         call flip_sonde(mxz,nn,uu,vv,rit,rjt,rk)
c      endif
c now interpolate to the laps levels in rk space
c      do k=1,nz
c       rr=k
c       do  n=1,nn-1
c        iflag=0        
c        if (rr.lt.rk(n+1).and.rr.ge.rk(n)) then ! point is between two dropsonde levels
c          aa=(rr-rk(n))/(rk(n+1)-rk(n))
c          tdrop(k)=tt(n)*(1.-aa)+tt(n+1)*aa
c          iflag=1
c          go to 4
c        endif
c      enddo 
c      if(iflag.eq.0) then
c       tdrop(k)=smsng
c      endif
c  4   enddo
c     close 11
c     else ! on lexist
       print*,'No tmg file for this time'
       istatus=istatus-2
c     endif


      deallocate(dd,ff)
      deallocate (uu,vv,tt)

      return

c49    print*,'No temp data in tmg file'
c     return

50    print*,'Error opening file: ',dum(1:len+13)
      return
      end
c-------------------------------------------------------------------
      subroutine readpig(a9_time,nx,ny,nz,lat,lon
     &,udrop,vdrop,tdrop,rri,rrj,rrit,rrjt,istatus)

c this subroutine reads the .pig and .tmg files to recover observed  
c u, v, bnd T profiles, the decimal i,j locations at the nz LAPS levels 
      real udrop(nz),vdrop(nz),tdrop(nz),rri(nz),rrj(nz), 
     & rrit(nz),rrjt(nz),dum2
c variables output:
c     udrop, vdrop tdrop: dropsonde u,v T obs at the nz laps levels
c     rri, rrj: real grid coordinates of dropsone position in grid space 
c
c Note: only one profile is allowed atm. arrays udrop, vdrop, tdrop
c       must be 2d to allow more than 1.
c
      Character*180 dum
      Character*9 a9_time
      Character*3 dum1
      integer len,istatus
      logical lexist
      real, allocatable, dimension(:) :: ri,rj,rk,dd,ff
     1,tt,uu,vv
      real lat(nx,ny),lon(nx,ny)

      integer mxz
      parameter (mxz=500)
      
      allocate(ri(mxz),rj(mxz),rk(mxz))
      allocate(dd(mxz),ff(mxz),tt(mxz))
      allocate(uu(mxz),vv(mxz))

      call get_r_missing_data(smsng,istatus)
      if(istatus.ne.1)then
         print*,'Error: returned from get_r_missing_data'
         return
      endif

c preset all dropsonde output to missing
      istatus=0
      rri=smsng
      rrj=smsng
      udrop=smsng
      vdrop=smsng
      tdrop=smsng

      pi=4.*atan(1.)
      rdpdg=pi/180. 
      call get_directory('pig',dum,len)
      dum(len+1:len+13)=a9_time//'.pig'
      inquire(file=dum,exist=lexist)
      if(lexist)then
        open (11, file=dum,form='formatted',status='old',err=50) 
        Do n=1,mxz
         read(11,*,end=1) ri(n),rj(n),rk(n),dd(n),ff(n),dum1
c change from (0,0,0) origin to (1,1,1) grid origin system.
         ri(n)=ri(n)+1
         rj(n)=rj(n)+1
         rk(n)=rk(n)+1
         if(dum1.eq.'  ') go to 1
        enddo
        print*, 'Suspect read in the pig file'
        istatus=-1
        goto 7

1       nsave=n-1 
        if(nsave.lt.1)then
           print*,'No wind data in pig file'
           istatus=-1
           goto 7
        endif

c convert dd ff to u,v
        do n=1,nsave
           call disp_to_uv(dd(n),ff(n),uu(n),vv(n))
        enddo

c since dropsonde is reverse order, flip it.
        if(rk(1).gt.rk(nsave))then
           call flip_sonde(mxz,nsave,uu,vv,ri,rj,rk)
        endif
c now interpolate to the laps levels in rk space
c if nsave is 1 (one level of data) then assign it to nearest laps level
c first clear out drop winds, make assignment and then close file
        do k=1,nz
         udrop(k)=smsng
         vdrop(k)=smsng
         rri(k)=smsng
         rrj(k)=smsng
        enddo
        if(nsave.eq.1) then
         k=rk(1)
         udrop(k)=uu(1)
         vdrop(k)=vv(1)
         rri(k)=ri(1)
         rrj(k)=rj(1) 
         go to  44 ! close file 
        endif
        do k=1,nz
           rr=k
           do n=1,nsave-1
              iflag=0        
              if( rr.lt.rk(n+1).and.rr.ge.rk(n))then !point laps level is between two dropsonde points.
                 aa=(rr-rk(n))/(rk(n+1)-rk(n))
                 udrop(k)=uu(n)*(1.-aa)+uu(n+1)*aa
                 vdrop(k)=vv(n)*(1.-aa)+vv(n+1)*aa
                 rri(k)=ri(n)*(1.-aa)+ri(n+1)*aa
                 rrj(k)=rj(n)*(1.-aa)+rj(n+1)*aa
             print*, 'k ri rj rkn rkn+1   ',k,ri(n),rj(n),rk(n),rk(n+1)
             print*, 'udrop,uu(n),uu(n+1) ',udrop(k),uu(n),uu(n+1) 
             print*, 'vdrop,vv(n),vv(n+1) ',vdrop(k),vv(n),vv(n+1) 
                 
                 iflag=1
                 go to 3
              endif
           enddo ! on n
           if(iflag.eq.0)then
              udrop(k)=smsng
              vdrop(k)=smsng
           endif
   3    enddo ! on k
   44   close(11)
      else
        print*,'No pig file at this time'
        istatus=-1
      endif

c now read the tmg file to get the dropsonde temps
7     call get_directory('tmg',dum,len)
      dum(len+1:len+13)=a9_time//'.tmg'
      nn=0
      inquire(file=dum,exist=lexist)
      if(lexist)then

       open (11, file=dum,form='formatted',status='old',err=50) 
       Do n=1,mxz
        nn=nn+1
       read(11,*,end=2) aa,bb,cc,ee, dum1           
       if(dum1.eq.'  ')  go to 2
       if(dum1.eq.'ACA') then    
        ri(nn)=aa
        rj(nn)=bb
        rk(nn)=cc
        tt(nn)=ee  
       endif
       enddo
c since dropsonde is reverse order, flip it.
2      nn=nn-1
       if(nn.lt.1)then
          istatus = istatus-2
          deallocate(ri,rj,rk,dd,ff,uu,vv,tt)
          goto 49
       endif
       if(nn.eq.1) then! there is only one ob; put at nearest laps level
         do k=1,nz! clear arrays
          tdrop(k)=smsng
          rrit(k)=smsng
          rrjt(k)=smsng
         enddo
         k=rk(1)
         tdrop(k)=tt(1)
         rrit(k)=ri(1)
         rrjt(k)=rj(1) 
         close (11) ! close file  and return
         return
       endif
       if(rk(1).gt.rk(nn))then
          call flip_sonde(nz,nn,uu,vv,ri,rj,rk)
       endif

c now interpolate to the laps levels in rk space
       do k=1,nz
       rr=k
       do  n=1,nn-1
        iflag=0        
        if (rr.lt.rk(n+1).and.rr.ge.rk(n)) then ! point is between two dropsonde levels
         aa=(rr-rk(n))/(rk(n+1)-rk(n))
         tdrop(k)=tt(n)*(1.-aa)+tt(n+1)*aa
         iflag=1
         go to 4
        endif
       enddo 
       if(iflag.eq.0) then
        tdrop(k)=smsng
       endif
   4   enddo

      else
       print*,'No tmg file for this time'
       istatus=istatus-2
       return
      endif

      istatus = 1  !successful return with both pig and tmg data

      deallocate(ri,rj,rk,dd,ff)
      deallocate (uu,vv,tt)
      return

49    print*,'No temp data in tmg file'
      return
50    print*,'Error opening file: ',dum(1:len+13)
      return
      end
c-------------------------------------------------------------------
      subroutine flip_sonde(mxz,nk,uu,vv,ri,rj,rk)

      implicit none
     
      integer mxz,nk,k
      real uu(mxz),vv(mxz),ri(mxz),rj(mxz),rk(mxz)
      real, allocatable, dimension(:) :: rii,rjj,rkk,u,v
      
      allocate(rii(mxz),rjj(mxz),rkk(mxz),u(mxz),v(mxz))
      do k=1,nk
         rii(nk-k+1)=ri(k)
         rjj(nk-k+1)=rj(k)
         rkk(nk-k+1)=rk(k)
         u(nk-k+1)=uu(k)
         v(nk-k+1)=vv(k)
      enddo
      ri=rii
      rj=rjj
      rk=rkk
      uu=u
      vv=v
      deallocate(rii,rjj,rkk,u,v)
      return
      end
c-------------------------------------------------------------------
      function ffz(ii,n,cor,seed)
c
c*********************************************************************
c
c     Pulls out a random value from a unit normal distribution with standard 
c     deviation = 1.  Input is
c     an integer seed 'ii' and an interation number 'n'. 'n' should
c     be >20 for best results.
c     if iswitch is on (1) then value seed is added to the returned 
c     gaussian number in order to correlate error
c     
c  Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       21 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c       07 Oct 1998  Peter Stamus, NOAA/FSL
c          Change 'ran' to 'ran1' function for portability.
c
c     Notes:
c
c*********************************************************************
c
      sum=0.
      do l=1,n
         sum=sum+ran1(ii)             
      enddo !l
      if(n.ne.1) then
       ffz=sum-float(n)/2.
      else
       ffz=sum
      endif
      ii=ii-1 
      if(ii.ge.0) ii=-ii-1
      ffz=ffz+cor*seed
c
      return
      end
c
c
      function rm(sm,st,mg,mmg,iii)
c this fuction is set up to create a joint pdf for a series of 
c gaussian pdfs. Means are passed thru sm(mg); std dev thru st(mg)
c mmg is actual number of combined distributions. rm is a 
c nonguassian random number times an error
      real sm(mg),st(mg),xo,sum,sum1
c compute overall fuction mean
      fact=1.001 ! ensures that overall function is larger than joint
      sum1=0 
      sum2=0
      sm(mmg+1)=0
      if(mmg.eq.2) then
       sum1=st(1)**2+st(2)**2
       sum2=sm(1)*st(2)**2+sm(2)*st(1)**2
       sm(mmg+1)=sum2/sum1
c      print*,'sm1,st1,sm2,st2 ',sm(1),st(1),sm(2),st(2)
      endif
      if (mmg.eq.3) then
       sum1=(st(1)*st(2))**2+(st(2)*st(3))**2+(st(3)*st(1))**2
       sum2=sm(1)*(st(2)*st(3))**2+sm(2)*(st(1)*st(3))**2+
     1       sm(3)*(st(1)*st(2))**2
       sm(mmg+1)=sqrt(sum2/sum1)
      endif
      if(mmg.ge.4) print*, 'not set up for four+ distributions'
c compute mean st
      if(mmg.eq.2) then
       sum1=st(1)**2+st(2)**2
       sum2=(st(1)*st(2))**2
       st(mmg+1)=sqrt(sum2/sum1)
c      print*,'over function sm st ', sm(mmg+1), st(mmg+1)
      endif
      if (mmg.eq.3) then
       sum1=(st(1)*st(2))**2+(st(2)*st(3))**2+(st(3)*st(1))**2
       sum2=(st(1)*st(2)*st(3))**2
       st(mmg+1)=sqrt(sum2/sum1)
      endif
      if(mmg.ge.4) print*, 'not set up for four+ distributions'

c generate gaussian random number
      icnt=0
      slast=0
      cor=.5
  5   xo=ffz(iii,20,cor,slast)*st(mmg+1) ! offset relative to sm(mmg+1)
      slast=xo
      ep=0.
      do k=1,mmg
       ep=-(xo+sm(mmg+1)-sm(k))**2/(st(k)*st(k)) + ep
      enddo
       ep1=-(xo/fact)**2/(st(mmg+1)*st(mmg+1))
       rat=ep-ep1
       cat=alog(ran1(iii))
c      print*,'xo,ep,ep1,rat,cat,iii', icnt+1,xo,ep,ep1,rat,cat,iii
       icnt=icnt+1
       if(cat.lt.rat) then
        rm=xo+sm(mmg+1)
c       print*,'Normal run rejection method. cnt= ',icnt
        return
       else
        if(icnt.eq.100) then
c        print*,'Problem run rejection method: set to mean. cnt= ',icnt
         rm=sm(mmg+1)
         return
        endif
        go to 5
       endif
       end

c ********************************************************************
      function ran1(idum)              
c
c     Function to generate a random number.  Use this instead of a
c     machine dependent one. idum must be a negative integer
c     
c     Original: Peter Stamus, NOAA/FSL  07 Oct 1998
c     Changes:  John McGinley, NOAA/FSL 20 Apr 00 - changed to ran2
c
c     Notes:
c        From "Numerical Recipes in Fortran", page 272.  There named ran2.
c
c*********************************************************************
c
        integer idum, ia1,ia2,im1,im2,imm1,iq1,iq2,ir1,ir2, ntab,ndiv
        real ran1, am, eps, rnmx
        parameter( ia1= 40014,
     &             ia2=40692,
     &             im1= 2147483563,
     &             im2= 2147483399,
     &             am= 1. / im1,
     &             imm1=im1-1,
     &             iq1= 53668, 
     &             iq2= 52774, 
     &             ir1= 12211,
     &             ir2= 3791,
     &             ntab = 32,
     &             ndiv = 1 + imm1/ntab    ,
     &             eps = 1.2e-7,
     &             rnmx = 1. - eps)
c
        integer idum2,j, k, iv(ntab), iy
        save iv, iy,idum2
        data idum2/123456789/,iv/ntab * 0/, iy/0/
c
c.....  Start here.
c
        if(idum.gt.0) idum=-idum
        if(idum.le.0 ) then  !initialize
           idum = max(-idum,1)           !prevent idum=0
           idum2=idum
           do j=ntab+8,1,-1              !load the shuffle table
              k = idum / iq1
              idum = ia1* (idum - k*iq1) - ir1*k
              if(idum .lt. 0) idum = idum + im1
              if(j .le. ntab) iv(j) = idum
           enddo !j
           iy = iv(1)
        endif
c
        k = idum / iq1                    !start here if not initializing
        idum = ia1* (idum - k*iq1) - ir1*k  
        if(idum .lt. 0) idum = idum + im1
        k = idum2/ iq2                    !start here if not initializing
        idum2= ia2* (idum2- k*iq2) - ir2*k  
        if(idum2.lt. 0) idum2= idum2+ im2
        j = 1 + iy/ndiv
        iy = iv(j)-idum2                  !output prev stored value and 
        iv(j) = idum                      !  refill shuffle table
        if(iy.lt.1)iy=iy+imm1
        ran1 = min(am*iy, rnmx)
c
        return
        end
c
c
      subroutine turb(u,v,om,t,phi,p,us,vs,ws,ts,ter,nx,ny,nz)
c this subroutine estimates turbulent components from the laps 
c grids.
c variables: u,v,om,t,phi,sh are the input laps grids 
c            us,vs,ws are the output turbulent components
      real u(nx,ny,nz),v(nx,ny,nz),om(nx,ny,nz)
     1,t(nx,ny,nz),phi(nx,ny,nz),p(nz),ts(nx,ny,nz)

      real us(nx,ny,nz),vs(nx,ny,nz),ws(nx,ny,nz)
      real ter(nx,ny)

      real, allocatable, dimension(:,:,:) :: tke_3d

      real, allocatable, dimension(:) :: p1d,u1d,v1d,t1d,z1d
     .,dtf


      allocate (p1d(nz),u1d(nz),v1d(nz),t1d(nz)
     .,z1d(nz),dtf(nz))
      allocate (tke_3d(nx,ny,nz))

c prepare input columns for dtf3
        
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C This is the beginning TKE computation
c
c subus called: vertirreg
c functions used: rf, RfKondo
C
      factr = 1.0  ! calibration factor for RUC-20 and RUC-10

      do i=1,nx
      do j=1,ny

      do k=1,nz
      p1d(k) = p(k)
      t1d(k) = t(i,j,k)
      u1d(k) = u(i,j,k)
      v1d(k) = v(i,j,k)
      z1d(k) = phi(i,j,k)
      zter=ter(i,j)
      enddo
c
c Error: arguments differ from those in subroutine
      call compute_dtf3(p1d,t1d,u1d,v1d,z1d,zter,dtf,nz)
C
C The result of the tke computation goes into array tke_3d
C from that we assume isotropy and recover the 3 turbulent components
c in u,v,w
c

      do l=1,nz
      tke3 = dtf(l)/factr
      if(tke3.gt.10.) tke3 = 10.
      tke_3d(i,j,l) = tke3
      us(i,j,l)=sqrt(tke3/3.) 
      vs(i,j,l)=us(i,j,l)
      ws(i,j,l)=vs(i,j,l)
      enddo

      enddo
      enddo

      return
      end

c
c________________________________________________________________________

      subroutine compute_dtf3(p,t,u,v,z,zter,tke_KH,nz)

c
c Adrian Marroquin FSL
c version modified for ITFA
c 11/09/98
c
c Km parameter calibrated to optimize PODn and
c PODy, Km = 75.0 m^2/s
c The following is a table of POD's (PODn = PODy) and
c thresholds for each one of the months from Nov 97
c to June 1998.
c
c
c Nov   Dec   Jan   Feb   Mar   Apr   May   Jun
c
c  66.   66.   61.   62.   68.   65.   64.   60.  GA (I>1, >20kft) all w
c  .58   .71   .72   .86   .73   .73   .42   .47  Thresholds
c  67.   66.   64.   63.   67.   64.   62.   59.  HA (I>1, >w120, >20kft)
c  .65   .67   .63   .71   .70   .68   .4    .4   Thresholds
c  63.   82.   72.   68.   76.   70.   60.   68.  GA (I>4, >20kft) all w
c  .55  1.105  .85   .8    .975  .76   .39   .5   Thresholds
c  62.   82.   72.   66.   75.   68.   59.   67.  HA (I>4, >w120, >20kft)
c  .6    1.15  .95   .91   .93   .78   .39   .55  Thresholds
c
c In the above table GA means General Aviation,
c HA heavy aircraft (commercial)
c I>1 turbulence intensities light or greater,
c >20kft aircraft flying above 20,000 ft,
c >w120 aircraft heavier than 120,000 lbs,
c and I>4 turbulence intensities moderate-to-severe or greater.
c
c WARNING: The above table was generated with TKE from DTF3 verified
c with PIREPs. The model output was from RUC2 (40-km), 40 isentropi
c levels. PIREPs from turbulence related to convection were not
c removed. DTF3 formulation is only applicable to turbulence from
c shear intabilities (clear-air turbulence) specially found in
c upper fronts.
c-------------------------------------------------------------------------
c DTF3 has been tested with 12-15 December 97 case study. In this case
c convective activity was at a minimum in the first half of December 97.
c During 12-15 Dec 97 a quasi-steady front moved across the US accompanied
c with a severe turbulence outbreak. For this case PODy = 93.6% and PODn =
c 76.1% obtained using thresholds for the month of December (see table
c above, PODn = PODy = 82% for December). The difference in PODy's is
c attributed to the inclusion of PIREPs from convection during the second
c half of December 97.
c
c DTF3 works well for moderate-to-severe turbulence or greater
c that affect heavy (commercial aircraft).
c
c-------------------------------------------------------------------------
c Constants from Stull (1988), page 219
c
      parameter(c1=1.44,c2=1.0,c3=1.92,
     *          c13 = c1/c3, c23 = c2/c3, ce = 0.19)
      parameter(alinf = 200.,akarm = 0.35,cr = 0.54)
     *
c
c
      parameter (cepn = 2.5, cepp = 0.76)
c
c pass p, t, u, v, and z
c
        real        p(nz),                     ! pressure in Pa
     1              t(nz),
     1              u(nz),
     1              v(nz),
     1              z(nz)

      real, allocatable, dimension(:):: brnt,shr,ri
     1,epsilon

      real tke_KH(nz)
      real tke
c
      data iepn3/0/
      data akm/75.0/
c
c--------------------------------------------------------------
c compute ri, brnt, and shr
c
      data r/287.04/,rocp/0.286/,g/9.8/
c
c Constants from Stull (1988), page 219
c
      data prands/2.5/
c
c
      if(.not. allocated(brnt))then
         allocate (brnt(nz),ri(nz),shr(nz),
     1            epsilon(nz))
      endif
      klev= nz
      pi1 = (p(1)/100000.)**rocp
      pi2 = (p(2)/100000.)**rocp
      th1 = t(1)/pi1
      th2 = t(2)/pi2

      do k=2,nz-1
      pi3 = (p(k+1)/100000.)**rocp
      th3 = t(k+1)/pi3
c
c vertical derivatives using pressure
c
      brunt = vertirreg(th1,th2,th3,
     *                 p(k-1),p(k),p(k+1))
      shru = vertirreg(u(k-1),u(k),u(k+1),
     *                 p(k-1),p(k),p(k+1))

      shrv = vertirreg(v(k-1),v(k),v(k+1),
     *                 p(k-1),p(k),p(k+1))

c
      beta = g*g*p(k)/(r*pi2*th2*th2)
      brnt(k) = -beta*brunt
      shr(k) = beta*p(k)*(shru*shru+shrv*shrv)/(r*pi2)
      th1 = th2
      th2 = th3
      pi1 = pi2
      pi2 = pi3
      ri(k) = brnt(k)/(shr(k)+1.e-10)
      enddo
c
      brnt(1) = brnt(2)
      shr(1) = shr(2)
      ri(1) = ri(2)
      brnt(nz) = brnt(nz-1)
      shr(nz) = shr(nz-1)

      ri(nz) = ri(nz-1)
c
      do k=1,nz
      if(ri(k).gt.120.) ri(k) = 120.
      enddo
c
c--------------------------------------------------------------
c
c now compute dissipation
c ztop equivalent to cpbl
c
      ztop = zter+3000.
      zsfc = zter
c
      DO K=1,klev
c
      zlev = z(k)
        if(ri(k).gt.0.01) then
        Rff = RfKondo(ri(k))
        else
        Rff = rf(ri(k))
        endif
      epsilon(k) = akm*shr(k)*(c13-c23*Rff)
      if(epsilon(k).lt.0.) epsilon(k) = 0.
c      tke_KH(k) = epsilon(k)
      if(iepn3.eq.0) then                 ! if iepn3 = 1, only epn
      if(brnt(k).le.0.) then
          tke_KH(k) = 0.
      else
         br = sqrt(brnt(k))
         tke_KH(k) = 0.7*epsilon(k)/(ce*br)
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if(zsfc.le.zlev.and.zlev.le.ztop) then                     !
         dz = zlev - zsfc                                           !
         alb = alinf*akarm*dz/(akarm*dz+alinf)                      !
         als = cr*sqrt(tke_KH(k))/br                                !
         all = amin1(alb,als)                                       !
         tke_KH(k) = (all*epsilon(k))**.666666                      !
         endif      

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      endif
      endif
c
      ENDDO    ! end of K-loop (klev)
c
      return
      end


c ============================================================
c
      subroutine get_laps_3d_analysis_data_isi(i4time,nx,ny,nz
     +,phi,t,u,v,sh,omo,istatus)
c
      implicit none

      integer   nx,ny,nz
      integer   i4time
      integer   istatus
      integer   lendlco
      integer   lends
      integer   lendt
      integer   lendw
      integer   lendsh
      integer   i,j,k
      real  r_missing_data
      real  phi(nx,ny,nz),t(nx,ny,nz)
     .       ,u(nx,ny,nz),v(nx,ny,nz),sh(nx,ny,nz)
     .       ,omo(nx,ny,nz)

      real, allocatable :: om(:,:,:)


      character*255 tempdir,winddir,sfcdir,shdir,lcodir
      character*125 comment
      character*31  tempext,windext,sfcext,shext,lcoext
      character*10  units
      logical found_lowest

      call get_r_missing_data(r_missing_data,istatus)
      if (istatus .ne. 1) then
         print *,'Error getting r_missing_data...Abort.'
         return
      endif

      shext='lq3'
      tempext='lt1'
      windext='lw3'
      sfcext='lsx'
      lcoext='lco'

      call get_directory(tempext,tempdir,lendt)
      call get_directory(windext,winddir,lendw)
      call get_directory(sfcext,sfcdir,lends)
      call get_directory(shext,shdir,lendsh)
      call get_directory(lcoext,lcodir,lendlco)

      call get_laps_3d(i4time,nx,ny,nz
     1  ,tempext,'ht ',units,comment,phi,istatus)

      if (istatus .ne. 1) then
         print *,'Error getting LAPS height data...Abort.'
         return
      endif
c
c *** Get laps temps
c
      call get_laps_3d(i4time,nx,ny,nz
     1  ,tempext,'t3 ',units,comment,t,istatus)

      if (istatus .ne. 1) then
         print *,'Error getting LAPS temp data...Abort.'
         return
      endif
c
c *** Get laps spec humidity
c
      call get_laps_3d(i4time,nx,ny,nz
     1  ,shext,'sh ',units,comment,sh,istatus)

      if(istatus .ne. 1)then
         print*,'Error getting LAPS sh  data ... Abort.'
         return
      endif
C    The specific humidity field uses the missing value
c    for below ground points, so we need to fill those in by
c    replicating the lowest valid value downward (upward in array
c    space).

      DO j = 1, ny
        DO i = 1, nx
          k = 1
          found_lowest = .false.

          DO WHILE (.NOT. found_lowest)
            IF (sh(i,j,k) .lt. 1.e37) THEN
              found_lowest = .true.
              sh(i,j,1:k) = sh(i,j,k)
            ELSE
              k = k + 1
              IF (k .ge. nz) THEN
                PRINT *, 'No valid SH found in column!'
                PRINT *, 'I/J = ', i,j
                STOP
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO

c    Make sure we got it right!
c      print *, 'Min/Max/Center values of SH:'
c      DO k = 1, nz
c        print *, minval(sh(:,:,k)),maxval(sh(:,:,k)),
c     +           sh(nx/2,ny/2,k)
c      ENDDO
c
c *** Get laps cloud omega
c


         istatus=0.
      if(istatus .ne. 1)then
         print*,'we are in the _isi version of get_laps_3d' 
         print*,'No LAPS Cld Omega data ....'
         print*,'Initializing omo array with zero'
         call zero3d(omo,nx,ny,nz)
      endif
c
c *** Get laps wind data.
c
      call get_laps_3d(i4time,nx,ny,nz
     1  ,windext,'u3 ',units,comment,u,istatus)

      if (istatus .ne. 1) then
         print *,'Error getting LAPS time 0 u3 data...Abort.'
         return
      endif

      call get_laps_3d(i4time,nx,ny,nz
     1  ,windext,'v3 ',units,comment,v,istatus)

      if (istatus .ne. 1) then
         print *,'Error getting LAPS time 0 v3 data...Abort.'
         return
      endif

      allocate (om(nx,ny,nz))

      call get_laps_3d(i4time,nx,ny,nz
     1  ,windext,'om ',units,comment,om,istatus)

      if (istatus .ne. 1) then
         print *,'Error getting LAPS time 0 v3 data...Abort.'
         return
      endif
      ! BLS commented this out as a test on 11/19/02
      !where(omo .eq. r_missing_data)omo=om

      deallocate(om)
      return
      end
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c---------------------------------------------------------------------
      function RfKondo(ri)
c
      parameter(c0=6.873,c1=7.)
c
c Rfc (critical flux Ri) = 0.143
c
      if(ri.gt.1.) then
      RfKondo = 1./c1
      else
      if(0.01.lt.ri.and.ri.le.1.) then
      d1 = 1.+c0*ri
      ahm = 1./(c0*ri+1./d1)
      RfKondo = ri*ahm
      endif
      endif
c
c for Ri < 0.01 use Rf (Yamada form)
c
      return
      end
c-------------------------------------------------------------------------
      subroutine maxminav2d(aux,mx,my,amin,amax)
c
      dimension aux(mx,my)
c
      anxy = float(mx*my)
c
      amax = -1.e10
      amin = 1.e10
c
      sum = 0.

      do i=1,mx
      do j=1,my
      sum = sum+ aux(i,j)
      if(aux(i,j).ge.amax) amax = aux(i,j)
      if(aux(i,j).le.amin) amin = aux(i,j)
      enddo
      enddo
c
      avfld = sum/anxy
      print*,'  max=',amax,'  min=',amin,' average=',avfld
c
      return
      end
c________________________________________________________________________-
      function vertirreg(f1,f2,f3,x1,x2,x3)
      dx1 = x2-x1
      dx2 = x3-x2
      rat1 = dx1/dx2
      rat2 = 1./rat1
      sdx = 1./(dx1+dx2)
      vertirreg = ((f3-f2)*rat1+(f2-f1)*rat2)*sdx
      return
      end
C_____________________________________________________________
      function prand(ri)
      data b/3.0/, a/6.873/
c
      if(ri.gt.1.) then
      prand = 7.*ri/b
      else
      if(0.01.lt.ri.and.ri.lt.1.) then
      prand = (a*ri*(1.+a*ri)+1.)/(b*(1.+a*ri))
      else
      prand = 1./b
      endif
      endif
c
      return
      end
c____________________________________________________________
      function rifunc(ri2,ri4,ric)
c
      fact = sqrt(4.+ri2*(1+4.*ric)/ric)
      anum = -(2.+ri2/ric)+ri4*fact/sqrt(ric)
      rifunc = anum/(2.*ri2)
c
      return
      end
c--------------------------------------------------------------------------
      function rf(ri)
c
      data c1/.056/,c2/.3/,c3/.3333/
      data a1/.78/,a2/.79/,b1/15./,b2/8./
c
      e1 = b1-6.*a1
      e2 = b1 + 12.*a1*(1.-c2)+3.*b2*(1.-c3)
      e3 = b1*(1.-3.*c1)-6.*a1
      e4 = b1*(1.-3.*c1)+12.*a1*(1.-c2)+9.*a2*(1.-c2)
      e5 = b1+3.*a1*(1.-c2)+3.*b2*(1.-c3)
c
      f1 = 0.5*a2*e5/(a1*e4)
      f2 = a1*e3/(a2*e5)
      f3 = 2.*a1*(e3*e5-2.*e1*e4)/(a2*e5*e5)
      f4 = a1*e3/(a2*e5)
      f42 = f4*f4
c
      rf = f1*(ri+f2-sqrt(ri*ri+f3*ri+f42))
c
      return
      end

