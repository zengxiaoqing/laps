      program qbalpe_main
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
      include 'trigd.inc'
      implicit none
c
      integer   nx,ny,nz
c
      real*4 dx(nx,ny),dy(nx,ny),dp(nz)
     .      ,p(nz),pstag(nz),ps(nx,ny),ter(nx,ny)
     .      ,lat(nx,ny),lon(nx,ny)
     .      ,phi(nx,ny,nz),t(nx,ny,nz)
     .      ,u(nx,ny,nz),v(nx,ny,nz),sh(nx,ny,nz)
     .      ,phib(nx,ny,nz),tb(nx,ny,nz)
     .      ,ub(nx,ny,nz),vb(nx,ny,nz),shb(nx,ny,nz)
     .      ,phibs(nx,ny,nz),tbs(nx,ny,nz)
     .      ,ubs(nx,ny,nz),vbs(nx,ny,nz),shbs(nx,ny,nz)
     .      ,phis(nx,ny,nz),ts(nx,ny,nz)
     .      ,us(nx,ny,nz),vs(nx,ny,nz),shs(nx,ny,nz)
c    .      ,lapsuo(nx,ny,nz),lapsvo(nx,ny,nz) !t=t0-dt currently not used
     .      ,lapsu(nx,ny,nz),lapsv(nx,ny,nz)   !t=t0
     .      ,lapssh(nx,ny,nz),lapslwc(nx,ny,nz)
     .      ,lapstemp(nx,ny,nz)
     .      ,lapsphi(nx,ny,nz)

      real*4 om(nx,ny,nz),omb(nx,ny,nz)
     .      ,omo(nx,ny,nz),oms(nx,ny,nz)
     .      ,ombs(nx,ny,nz)
     .      ,nu(nx,ny,nz),nv(nx,ny,nz)
     .      ,fu(nx,ny,nz),fv(nx,ny,nz)
     .      ,wb(nx,ny,nz)
     .      ,re,rdpdg,po,cappa
     .      ,delo,dt
     .      ,gamo
     .      ,erru(nx,ny,nz),errub(nx,ny,nz)
     .      ,errphi(nx,ny,nz),errphib(nx,ny,nz)

      real*4 grid_spacing_actual_m
     .      ,grid_spacing_cen_m
     .      ,pdif,dpbl,dpblf
     .      ,u_grid,v_grid
     .      ,u_true,v_true
     .      ,dpp

      real*4 g,sumdt,omsubs,sk,bnd,ff,fo,err,rog,rod
     .      ,sumdz,sumr,sumv2,snxny,sumf,sumt,cl,sl
     .      ,sumtscl,sumkf,sumks,sldata,den,sumom2

c made 2d 2-20-01 JS.
      real*4 terscl(nx,ny)
      real*4 tau(nx,ny)
      real*4 ro(nx,ny)
      integer ks(nx,ny)
      integer kf(nx,ny)
      integer ksij,kfij
c
      integer   itmax,lmax
     .         ,masstime,windtime,sfctime,omtime
     .         ,i,j,k,ll,istatus

      integer   itstatus
      integer   init_timer
      integer   ishow_timer
     
      integer   lend
      integer   lends
      integer   lenvg
c
      logical lrunbal
      logical lrotate/.false./
      logical larray_diag/.false./

      character*255 staticdir,sfcdir
      character*255 generic_data_root
      character*125 comment
      character*40  vertical_grid
      character*31  staticext,sfcext
      character*10  units
      character*9   a9_time

c    Added by B. Shaw, 4 Sep 01
      real*4, allocatable :: lapsrh(:,:,:)
      real*4, external :: ssh, make_rh
      real*4 shsat
c_______________________________________________________________________________
c
      call get_balance_nl(lrunbal,istatus)
      if(istatus.ne.0)then
         print*,'error getting balance namelist'
         stop
      endif
      print*,'lrotate = ',lrotate
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
      call get_systime(masstime,a9_time,istatus)
      if(istatus .ne. 1)go to 999
c
c get pressures and determine pressure intervals.
c
      call get_pres_1d(masstime,nz,p,istatus)
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
      itmax=200  !max iterations for relaxation
c
      windtime=(masstime+1800)/3600*3600
      sfctime=windtime
      omtime=windtime/10800*10800
c
c *** Get laps grid lat, lons and grid spacings.
c
c     call get_laps_lat_lon(staticdir,staticext
c    .                     ,nx,ny,lat,lon,istatus)

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
c *** Get LGA fields
      call get_modelfg_3d(masstime,'U3 ',nx,ny,nz,ub,istatus)
      call get_modelfg_3d(masstime,'V3 ',nx,ny,nz,vb,istatus)
      call get_modelfg_3d(masstime,'T3 ',nx,ny,nz,tb,istatus)
      call get_modelfg_3d(masstime,'HT ',nx,ny,nz,phib,istatus)
      call get_modelfg_3d(masstime,'SH ',nx,ny,nz,shb,istatus)
      call get_modelfg_3d(masstime,'OM ',nx,ny,nz,omb,istatus)

c convert background sh to rh
c  COMMENTED OUT...Now do everything in specific humidity
c     do k=1,nz
c        do j=1,ny
c        do i=1,nx
c           rhb(i,j,k)=make_rh(p(k)/100.,tb(i,j,k)-273.15
c    .,rhb(i,j,k)*1000.,-132.)*100.
c        enddo
c        enddo
c     enddo

c
c *** Get laps analysis grids.
c
      call get_laps_analysis_data(masstime,nx,ny,nz
     +,lapsphi,lapstemp,lapsu,lapsv,lapssh,omo,lapslwc,istatus)
c omo is the cloud vertical motion from lco
      if (istatus .ne. 1) then
         print *,'Error getting LAPS analysis data...Abort.'
         stop
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

c     enddo
c     enddo
c
c *** Get laps surface elevations.
c
c     call get_laps_sfc_elev(staticdir,staticext,nx,ny
c    .                      ,ter,istatus)
c
c *** Get laps surface pressure.
c
      call get_laps_2d(masstime,sfcext,'PS ',units,
     1                  comment,nx,ny,ps,istatus)
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

c set the resolvable scale based on data nyquist inteval=4 spacing
c this should be automated to be computed by the actual number of obs/area
c guess number of obs over grid, put in sldata for now.
      sldata=200.!number of good upper air obs
      sldata=4.*sqrt(float((nx-1)*(ny-1))*dx(nx/2,ny/2)**2/sldata)


c set model scale - a low end wave resolvable by the grid
      sl=4.*dx(nx/2,ny/2)

      fo=14.52e-5   !2*omega

      if(.false.)then
         call terrain_scale_vartau(nx,ny,nz,ter,lapsphi
     &,ks,kf,terscl)
      else
         do k=1,nz
            if(p(k).eq.85000)ks(1,1)=k
            if(p(k).eq.50000)kf(1,1)=k
c        ks(1,1)=6
c        kf(1,1)=11
         enddo

         call terrain_scale(nx,ny,ter,terscl(1,1))
      endif

      print*,'terrain scale = ',terscl(1,1)
      print*,'ks/p = ',ks(1,1),p(ks(1,1))
      print*,'kf/p = ',kf(1,1),p(kf(1,1))


      do i=1,nx
      do j=1,ny

         terscl(i,j)=terscl(1,1)   !make it constant over domain for now.
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
            erru(i,j,k)=(1.0*(1.+float(k-1)*.10))**(-2)
c           erru(i,j,k)=(1.5*(1.+float(k-1)*.25))**(-2)
            errub(i,j,k)=(1.5*(1.+float(k-1)*.3))**(-2)
            errphi(i,j,k)=(15.*(1.+float(k-1)*.1)*g)**(-2)
            errphib(i,j,k)=(30.*(1.+float(k-1)*.1)*g)**(-2)
c
c vertical motions in clear areas come in as the missing data parameter.
c replace missing cloud vv's with background vv's. Unless there is cloud
c we seek to replicate the background vertical motions
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
      do j=1,ny
      do i=1,nx
         kfij=kf(i,j)
         ksij=ks(i,j)
         dpp=p(ksij)-p(kfij)
c scale tau as of 1/(omega)**2 where omega is the background rms omega              
         tau(i,j)= 1./sumom2 


      enddo
      enddo
c     delo is scaled as 10% of expected eqn of motion residual ro*U**2/L
      rod=sqrt(sumv2)/(sumf*sldata)
      rog=sqrt(sumv2)/(sumf*sl)
      if(rog.gt.1.) rog=1.
      if(rod.gt.1.) rod=1.
      delo=100.*sl**2/sumv2**2/rod**2           
c
c print these arrays now.
      print*,'/dthet/thet/dz/den/N/V/f/delo/tau:' 
     &,sumdt,sumt,sumdz,den,sumr,sqrt(sumv2),sumf,delo,tau(1,1)
      print*,'Omegab,Froude Num/DataRossby Num,aspectP/X:',sqrt(sumom2),
     &    (sumr*sumdz/sqrt(sumv2)),rod ,(dpp/dx(nx/2,ny/2))
      
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
      call balstagger(ub,vb,phib,tb,shb,omb,ubs,vbs,
     &phibs,tbs,shbs,ombs,nx,ny,nz,p,ps,1) 
   
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

      call terbnd(u,v,om,nx,ny,nz,ps,p,bnd)
      call terbnd(us,vs,oms,nx,ny,nz,ps,p,bnd)
      call terbnd(ubs,vbs,ombs,nx,ny,nz,ps,p,bnd)
c

c
c ****** Execute mass/wind balance.
c  set maximum relaxation correction for phi in GPM
      err=1.0
c returns staggered grids of full fields u,v,phi

      call balcon(phis,us,vs,oms,phi,u,v,om,phibs,ubs,vbs,ombs,
     .       ts,rod,delo,tau,itmax,err,erru,errphi,errub,errphib
     .           ,nu,nv,fu,fv,nx,ny,nz,lat,dx,dy,ps,p,dp,lmax)

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

c
c
c *** destagger and Write out new laps fields.
c
c the non-staggered grids must be input with intact boundaries from background 

      call balstagger(ub,vb,phib,tb,
     & shb,omb,u,v,phi,t,sh,om,nx,ny,nz,p,ps,-1) 
c
c   Prior to applying boundary subroutine put non-staggered grids back into
c   u,v,om,t,sh,phi.

      call move_3d(phib,phi,nx,ny,nz)
      call move_3d(ub,u,nx,ny,nz)
      call move_3d(vb,v,nx,ny,nz)
      call move_3d(tb,t,nx,ny,nz)
      call move_3d(omb,om,nx,ny,nz)
      call move_3d(shb,sh,nx,ny,nz)

      call get_laps_2d(masstime,sfcext,'PS ',units,
     1                  comment,nx,ny,ps,istatus)

      if(lrotate)then

         print*,'rotate u/v back to true north'
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

c       Ensure the specfic humidity does not exceed
c       the saturation value for this temperature

        lapssh(i,j,k) = MIN(shsat,lapssh(i,j,k))
c
c       If cloud water is present, set the sh equal
c       to the saturation value

        if (lapslwc(i,j,k).gt.0.) lapssh(i,j,k)=shsat

c       Finally, rediagnose RH wrt liquid from the 
c       modified sh field

        lapsrh(i,j,k) = make_rh(p(k)*0.01,t(i,j,k)-273.15
     .                    ,lapssh(i,j,k)*1000.,-132.)*100.         
        lapsrh(i,j,k) = MIN(lapsrh(i,j,k),100.)
        lapsrh(i,j,k) = MAX(lapsrh(i,j,k),1.0)
      enddo
      enddo
      enddo     
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

      call write_bal_laps(masstime,phi,u,v,t,om,lapsrh,lapssh
     .                   ,nx,ny,nz
     .                   ,p,istatus)
      if(istatus.ne.1)then
         write(6,*)'error writing balance fields'
         return
      endif
      deallocate(lapsrh)
 999  return
      end

c
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
      integer*4 nx,ny,nz
     .         ,nxm1,nym1,nzm1
     .         ,i,j,k,ks
c
      real*4 u(nx,ny,nz),v(nx,ny,nz)
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
      ieast=ii-ispan/2+1
      if(ieast.lt.1) ieast=1
      iwest=ii+ispan/2+1
      if(iwest.gt.nx)iwest=nx 
      do j=jnorth,jsouth,-1
c      write(9,1001)  j,ieast,iwest
 1001  format(1x,3i8)
       print*,j,ieast,iwest
       write(6,1000) (a(i,j,kk),i=ieast,iwest)
c      write(9,1000) (a(i,j,kk),i=ieast,iwest)
 1000  format(1x,7e11.5)
      enddo
      return
      end
c     
c===============================================================================
c
c
      subroutine balcon(to,uo,vo,omo,t,u,v,om,tb,ub,vb,omb,tmp,
     .   rod,delo,tau,itmax,err,erru,errph,errub,errphb,nu,nv
     .     ,fu,fv,nx,ny,nz,lat,dx,dy,ps,p,dp,lmax)
c
c *** Balcon executes the mass/wind balance computations as described
c        mcginley (Meteor and Appl Phys, 1987) except that
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
      real*4 t(nx,ny,nz),to(nx,ny,nz),tb(nx,ny,nz)
     .      ,u(nx,ny,nz),uo(nx,ny,nz),ub(nx,ny,nz)
     .      ,v(nx,ny,nz),vo(nx,ny,nz),vb(nx,ny,nz)
     .      ,om(nx,ny,nz),omo(nx,ny,nz),omb(nx,ny,nz)
     .      ,nu(nx,ny,nz),nv(nx,ny,nz),tmp(nx,ny,nz)
     .      ,fu(nx,ny,nz),fv(nx,ny,nz)
     .      ,dx(nx,ny),dy(nx,ny),dp(nz)
     .      ,ps(nx,ny),p(nz)
     .      ,lat(nx,ny),ff(nx,ny)
     .      ,erru(nx,ny,nz),errph(nx,ny,nz)
     .      ,errub(nx,ny,nz),errphb(nx,ny,nz)

      real*4 err,rdpdg,bnd,g,fo,r,re,ovr
     .      ,ang,f,cotmax,sin1,fs,cos1,beta
     .      ,a,bb,cortmt,rod
     .      ,dudy,dvdx,dnudx,dnvdy,tt,uot,vot,tot
     .      ,dt2dx2,dt2dy2,slap,force,rest,cot
     .      ,cotma1,cotm5,rho,cotm0,erf,dtdx,dtdy,nuu,nvv
     .      ,dldp,dldx,dldy,tsum,r_missing_data
     .      ,usum,vsum,delo,fuu,fvv
     .      ,angu,angv,dyu,dxv

      real*4 term1,term2,term3,term4,term5,term6
     .      ,term7,term8,term9,term10,term11
     .      ,euoay,euoax,snv,snu,dfvdy,dfudx
     .      ,eueub,fob,foax,foay
     .      ,fu2,fv2,fuangu,fvangv,dt

c 2d arrays now (JS 2-20-01)
      real*4 tau(nx,ny)

c these are used for diagnostics
      integer nf
      parameter (nf=10)
      real*4  data(nf),fmax,sum
      real*4  fldmax(nf,nz),fldmin(nf,nz)
      integer fldmxi(nf,nz),fldmxj(nf,nz),ismx,jsmx,ksmx
      integer fldmni(nf,nz),fldmnj(nf,nz)

      integer   itstatus
      integer   init_timer
      integer   ishow_timer

      logical larray_diag/.false./

      real, allocatable, dimension(:,:,:) :: aaa,bbb
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
      call analzo(to,to,uo,uo,vo,vo,omo,omo
     .                ,nu,nv,fu,fv,delo,tau
     .                ,nx,ny,nz
     .                ,lat,dx,dy,ps,p,dp,l,lmax)

c create perturbations
c owing to an artifact of coding the t array is phi
      print*,'initialize timer'
      itstatus=init_timer()
      print*,'-------------------------'

      do l=1,lmax



       write(6,*) '|||||||||BALCON ITERATION NUMBER ',l,' ||||||||||'
       print*,'-----------------------------------------------------'
c      write(9,*) '|||||||||BALCON ITERATION NUMBER ',l,' ||||||||||'
c  set convergence error for relaxation of lamda...set erf to desired 
c  accuracy of wind
       erf=.01 !m/sec
       erf=erf*dx(nx/2,ny/2)
c apply continuity to input winds
       call leib_sub(nx,ny,nz,erf,tau,erru
     .,lat,dx,dy,ps,p,dp,uo,u,vo,v,
     . omo,om,omb,l,lmax)

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

      call frict(fu,fv,ub,vb,uo,vo,p,ps,tmp
     .                 ,nx,ny,nz,dx,dy,dp,dt,bnd)
      call nonlin(nu,nv,uo,vo,ub,vb,om,omb
     .           ,nx,ny,nz,dx,dy,dp,dt,bnd,rod)
c *** Compute new phi (t array) using relaxation on eqn. (2).
c        beta*dldx term is dropped to eliminate coupling with lambda eqn.
c
       fmax=-1.e30
       do 1 it=1,itmax

         cotmax=0.
         sum=0.
         do 2 j=2,nym1
          do 2 i=2,nxm1
            

           do 2 k=ks,kf

               if(k.ne.kf)then
               if(ps(i,j).lt.p(k+1)) then 
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
   
      call diagnose(to,nx,ny,nz,ismx,jsmx,ksmx,7,'INPUT GEOPOTENTIALS')
      call diagnose(uo,nx,ny,nz,ismx,jsmx,ksmx,7,'INPUT U-COMPONENT  ')
      call diagnose(vo,nx,ny,nz,ismx,jsmx,ksmx,7,'INPUT V-COMPONENT  ')
      call diagnose(omo,nx,ny,nz,ismx,jsmx,ksmx,7,'INPUT OMEGA')

c     write(9,1000) it,cotmax,ovr,cotma1
c     erf=0.
c
c *** Compute new u, v, omega using eqns. (4), (5), (6) with new phi and
c        without the lagrange multiplier terms.
c
       call initmxmn(nf,nz
     &,fldmax,fldmin,fldmxi,fldmxj,fldmni,fldmnj)

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
       call diagnose(t,nx,ny,nz,ismx,jsmx,ksmx,7,'BALANCED PERT PHIS')
       call diagnose(u,nx,ny,nz,ismx,jsmx,ksmx,7,'BALANCED U-PERTUR ')
       call diagnose(v,nx,ny,nz,ismx,jsmx,ksmx,7,'BALANCED V-PERTUR ')
       call diagnose(omo,nx,ny,nz,ismx,jsmx,ksmx,7,'BALANCED OMEGA')       

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


       itstatus=ishow_timer()
       print*,'------------------------------------------'
       print*,'Elapsed time (after leib_sub) sec: ',itstatus
       print*,'------------------------------------------'


      enddo ! on lmax

c apply continuity to final winds...done with backgrounds use ub,vb,omb
       call leib_sub(nx,ny,nz,erf,tau,erru
     .,lat,dx,dy,ps,p,dp,u,ub,v,vb,
     . om,omb,omb,l,lmax)
c evaluate dynamic balance and continuity
      call analzo(t,to,ub,uo,vb,vo,omb,omo
     .                ,nu,nv,fu,fv,delo,tau
     .                ,nx,ny,nz
     .                ,lat,dx,dy,ps,p,dp,l,lmax)
c move adjusted fields (uo,vo,omo) to solution fields 
       call move_3d(ub,u,nx,ny,nz)
       call move_3d(vb,v,nx,ny,nz)
       call move_3d(omb,om,nx,ny,nz)
      print*,'------------------------------------------------'
      itstatus=ishow_timer()
      print*,'Elapsed time end of balcon loop (sec): ',itstatus
      print*,'------------------------------------------------'

      deallocate (aaa,bbb)
      deallocate (dxx,dx2,dxs,dyy,dy2,dys
     1,fx,fy,ffx,ffy)


      return
      end
c
c ---------------------------------------------------------------
c
      subroutine leib_sub(nx,ny,nz,erf,tau,erru
     .,lat,dx,dy,ps,p,dp,uo,u,vo,v,
     . omo,om,omb,l,lmax)

      implicit none

      integer nx,ny,nz
      integer nxm1,nym1,nzm1
      integer i,j,k,ks,l,lmax

      real*4      u(nx,ny,nz),uo(nx,ny,nz)
     .      ,v(nx,ny,nz),vo(nx,ny,nz)
     .      ,om(nx,ny,nz),omo(nx,ny,nz),omb(nx,ny,nz)
     .      ,erru(nx,ny,nz)
     .      ,lat(nx,ny),dx(nx,ny),dy(nx,ny)
     .      ,ps(nx,ny),p(nz),dp(nz)

      real*4 dldx,dldy,dldp,sum,cnt
     .,a,erf,bnd
      real*4 tau(nx,ny)

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
c
c ****** Perform 3-d relaxation.
c
      call leibp3(slam,f3,200,erf,h
     .  ,nx,ny,nz,dx,dy,ps,p,dp)
c
c ****** Compute new u, v, omega by adding the lagrange multiplier terms.
co
      sum=0
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
         if (u(i,j,k) .ne. bnd) u(i,j,k)=uo(i,j,k)+.5*dldx/a
         if (v(i,j,k) .ne. bnd) v(i,j,k)=vo(i,j,k)+.5*dldy/a
         if (omo(i,j,k).ne.bnd) om(i,j,k)=omo(i,j,k)+
     &     .5*dldp/tau(i,j)
       sum=sum+(u(i,j,k)-uo(i,j,k))**2+(v(i,j,k)-vo(i,j,k))**2
       cnt=cnt+1.
      enddo
      enddo
      do i=1,nx
         a=2.*erru(i,ny,k)
         dldp=(slam(i,ny,k)-slam(i,ny,k+1))/dp(k+ks)
         dldy=(slam(i+1,ny+1,k+1)-slam(i+1,ny,k+1))/dy(i,ny)
         if (v(i,ny,k) .ne. bnd) v(i,ny,k)=vo(i,ny,k)+.5*dldy/a
         if (omo(i,ny,k).ne.bnd) om(i,ny,k)=omo(i,ny,k)+
     &     .5*dldp/tau(i,ny)
      enddo
      do j=1,ny
         a=2.*erru(nx,j,k)
         dldp=(slam(nx,j,k)-slam(nx,j,k+1))/dp(k+ks)
         dldx=(slam(nx+1,j+1,k+1)-slam(nx,j+1,k+1))/dx(nx,j)
         if (u(nx,j,k) .ne. bnd) u(nx,j,k)=uo(nx,j,k)+.5*dldx/a
         if (omo(nx,j,k).ne.bnd) om(nx,j,k)=omo(nx,j,k)+
     &     .5*dldp/tau(nx,j)
      enddo
      enddo
c   print out rms vector adjustment
      print*, 'RMS Vector adjustment after continuity applied ' 
      print*, sqrt(sum/cnt)

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
      real*4 t(nx,ny,nz),to(nx,ny,nz)
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

      real*4 tau(nx,ny)
      real*4 r_missing_data

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
     . /1x,' rms friction vector        :',2e12.4 
     . /1x,' rms nonlinear term value   :',2e12.4)
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
      real*4 t(nx,ny,nz),to(nx,ny,nz)
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

      real*4 tau(nx,ny)
      real*4 r_missing_data

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
      real*4 u(nx,ny,nz),v(nx,ny,nz)
     .      ,om(nx,ny,nz)
     .      ,ps(nx,ny),p(nz)
     .      ,bnd
c_______________________________________________________________________________
c
      nym1=ny-1
      nxm1=nx-1
      nzm1=nz-1
      do j=2,nym1
      do i=2,nxm1
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
         do j=1,ny
            if(ps(1,j) .le. p(k)) then
c              u(1,j-1,k-1)=bnd   !commented out JS 01-17-01
               u(1,j,k-1)=bnd     !added this
               om(1,j,k)=bnd
               om(1,j,k-1)=bnd
            endif
            if(ps(nx,j) .le. p(k)) then
c              u(nx,j-1,k-1)=bnd
c              u(nxm1,j-1,k-1)=bnd !commented JS 01-17-01
               u(nx,j,k-1)=bnd
               u(nxm1,j,k-1)=bnd
               om(nx,j,k)=bnd
               om(nx,j,k-1)=bnd
            endif
         enddo

         do i=1,nx
            if(ps(i,1) .le. p(k)) then
c              v(i-1,1,k-1)=bnd
               v(i,1,k-1)=bnd
               om(i,1,k)=bnd
               om(i,1,k-1)=bnd
            endif
            if(ps(i,ny) .le. p(k)) then
c              v(i-1,ny,k-1)=bnd    !commented JS 01-17-01
c              v(i-1,nym1,k-1)=bnd
               v(i,ny,k-1)=bnd
               v(i,nym1,k-1)=bnd
               om(i,ny,k)=bnd
               om(i,ny,k-1)=bnd
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
      real*4 fu(nx,ny,nz),fv(nx,ny,nz)
     .      ,ub(nx,ny,nz),vb(nx,ny,nz) 
     .      ,u(nx,ny,nz),v(nx,ny,nz),t(nx,ny,nz),bnd  
     .      ,dx(nx,ny),dy(nx,ny),dp(nz),ps(nx,ny),p(nz)
     .      ,ddp(nz),ddp2(nz)
     .      ,dt,dudx,dvdy,dudy,dvdx
     .      ,dudp,dvdp,dudt,dvdt
     .      ,rmx2d,rmn2d

      real*4 gdtvdp_save(nx,ny),gdtudp_save(nx,ny)
      real*4 gdtvdp,gdtudp

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
           gdtudp=den*cd*grav*abs(u(i,j,k))*u(i,j,k)/dp(k)       
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
      real*4 nu(nx,ny,nz),nv(nx,ny,nz)
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
      real*4 f3(nx+1,ny+1,nz+1),h(nx+1,ny+1,nz+1)
     .      ,erru(nx,ny,nz)
     .      ,u(nx,ny,nz),v(nx,ny,nz)
     .      ,om(nx,ny,nz),lat(nx,ny)
     .      ,omb(nx,ny,nz)
     .      ,dx(nx,ny),dy(nx,ny),dp(nz)
     .      ,dpp,aa,dxx,dyy
     .      ,cont,formax
      real*4 tau(nx,ny)
c_______________________________________________________________________________
c
      print *,'fthree'
      formax=0.
      do k=2,nz
         dpp=dp(k)
         do j=2,ny
         do i=2,nx
            js=1
            if (j .eq. 2) js=0
            dyy=dy(i,j)*float(js+1)
            is=1
            if (i .eq. 2) is=0
            dxx=dx(i,j)*float(is+1)
            h(i,j,k)=erru(i,j,k)/(tau(i,j))
            f3(i,j,k)=-2.*erru(i,j,k)*
     &    ((u(i,j-1,k-1)-u(i-1,j-1,k-1))/dx(i,j)
     .    +(v(i-1,j,k-1)-v(i-1,j-1,k-1))/dy(i,j))
     .    -h(i,j,k)*(om(i,j,k-1)-om(i,j,k))/dpp
           if (abs(f3(i,j,k)) .ge. formax) then
               formax=abs(f3(i,j,k))
            endif
         enddo
         enddo
      enddo
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
      real*4 sol(nx+1,ny+1,nz+1)
     .      ,force(nx+1,ny+1,nz+1)
     .      ,h(nx+1,ny+1,nz+1)
     .      ,dx(nx,ny),dy(nx,ny)
     .      ,ps(nx,ny),p(nz),dp(nz)
     .      ,erf,si,sj,sk
     .      ,ovr,erb,hh,ertm,ermm,corlm
     .      ,dx2,dx1s,dy2,dy1s,dz,dz2,dz1s
     .      ,aa,cortm,res,cor,corb,corlmm
     .      ,reslm,rho,cor0,cor5
c_______________________________________________________________________________
c
c *** Relaxer solver...eqn must be.......                                 
c        sxx+syy+h*szz-force=0   
c
      print *,'leibp3'
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
         do j=2,ny
         js=1
         if(j.eq.ny)js=0
         do i=2,nx
            is=1
            if(i.eq.nx) is=0
            dx1s=dx(i,j)*dx(i,j)
            dy1s=dy(i,j)*dy(i,j)
            do 2 k=2,nz
            dz=dp(k)
            dz1s=dz*dz
            if (ps(i,j).lt.p(k)) go to 2 
            aa=h(i,j,k)
            if(ps(i+is,j).lt.p(k)) sol(i+1,j,k)=sol(i,j,k)
            if(ps(i-1 ,j).lt.p(k)) sol(i-1,j,k)=sol(i,j,k)
            if(ps(i,j+js).lt.p(k)) sol(i,j+1,k)=sol(i,j,k)
            if(ps(i,j-1).lt.p(k)) sol(i,j-1,k)=sol(i,j,k)
            if(ps(i,j).lt.p(k-1)) sol(i,j,k-1)=sol(i,j,k)
            cortm=-2./dx1s-2./dy1s-aa*2./dz1s
                  res=(sol(i+1,j,k)+sol(i-1,j,k))/dx1s+
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
1002       format(1x,'LIEBP3:ovr rlxtn const at fnl ittr = ',e10.4)
1001       format(1x,'iterations= ',i4,' max residual= ',e10.3,
     .    ' max correction= ',e10.3, ' first iter max cor= ',e10.3,
     .    'max bndry error= ',e10.3)
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
      real*4 sol(nx,ny,nz),force(nx,ny,nz)
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

       deallocate (wr)
       
       return

      else! de-stagger

c horizontal destagger
c re compute hydrostatic virtual t from staggered adjusted phis 
c use a second-order scheme
       do k=2,nz-1
        do j=1,ny
         do i=1,nx
          fo=(phis(i,j,k)-phis(i,j,k-1))/
     &                              alog(p(k)/p(k+1))
          so=0.
          if(k+2.le.nz.and.k.ge.3) 
     &    so=.5*((phis(i,j,k+1)-phis(i,j,k-1))/alog(p(k)/p(k+2))-
     &         (phis(i,j,k)-phis(i,j,k-2))/alog(p(k-1)/p(k+1)))
         ts(i,j,k)=(fo-so)/r          
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
