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
     .      ,u(nx,ny,nz),v(nx,ny,nz),rh(nx,ny,nz)
     .      ,phib(nx,ny,nz),tb(nx,ny,nz)
     .      ,ub(nx,ny,nz),vb(nx,ny,nz),rhb(nx,ny,nz)
     .      ,phibs(nx,ny,nz),tbs(nx,ny,nz)
     .      ,ubs(nx,ny,nz),vbs(nx,ny,nz),rhbs(nx,ny,nz)
     .      ,phis(nx,ny,nz),ts(nx,ny,nz)
     .      ,us(nx,ny,nz),vs(nx,ny,nz),rhs(nx,ny,nz)
c    .      ,lapsuo(nx,ny,nz),lapsvo(nx,ny,nz) !t=t0-dt currently not used
     .      ,lapsu(nx,ny,nz),lapsv(nx,ny,nz)   !t=t0
     .      ,lapsrh(nx,ny,nz),lapslwc(nx,ny,nz)
     .      ,lapstemp(nx,ny,nz)
     .      ,lapsphi(nx,ny,nz)

      real*4 om(nx,ny,nz),omb(nx,ny,nz)
     .      ,omo(nx,ny,nz),oms(nx,ny,nz)
     .      ,ombs(nx,ny,nz)
     .      ,nu(nx,ny,nz),nv(nx,ny,nz)
     .      ,fu(nx,ny,nz),fv(nx,ny,nz)
     .      ,wb(nx,ny,nz)
     .      ,re,rdpdg,po,cappa
     .      ,delo,tau,dt
     .      ,gamo
     .      ,erru(nx,ny,nz),errub(nx,ny,nz)
     .      ,errphi(nx,ny,nz),errphib(nx,ny,nz)

      real*4 grid_spacing_actual_m
     .      ,pdif,dpbl,dpblf
     .      ,u_grid,v_grid
     .      ,u_true,v_true

      real*4 g,sumdt,omsubs,sk,terscl,bnd,ff,fo,err
     .      ,sumdz,sumr,sumv2,snxny,sumf,sumt,cl,sl,ro
c
      integer   itmax,lmax
     .         ,masstime,windtime,sfctime,omtime
     .         ,i,j,k,ll,istatus
     .         ,ks,kf

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
c
c_______________________________________________________________________________
c
      call get_balance_nl(lrunbal,gamo,delo,istatus)
      if(istatus.ne.0)then
         print*,'error getting balance namelist'
         stop
      endif
      print*,'lrunbal = ',lrunbal
      print*,'gamo = ',gamo
      print*,'delo = ',delo
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
      call get_laps_lat_lon(staticdir,staticext
     .                     ,nx,ny,lat,lon,istatus)
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
      call get_modelfg_3d(masstime,'SH ',nx,ny,nz,rhb,istatus)
      call get_modelfg_3d(masstime,'OM ',nx,ny,nz,omb,istatus)
c
c *** Get laps analysis grids.
c
      call get_laps_analysis_data(masstime,nx,ny,nz
     +,lapsphi,lapstemp,lapsu,lapsv,lapsrh,omo,lapslwc,istatus)
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
      call get_laps_sfc_elev(staticdir,staticext,nx,ny
     .                      ,ter,istatus)
c
c *** Get laps surface pressure.
c
      call get_laps_2d(masstime,sfcext,'PS ',units,
     1                  comment,nx,ny,ps,istatus)

      sl=1200000. ! four times the average data spacing
      fo=14.52e-5   !2*omega
      terscl=3000.
      ks=6
      kf=14
      do j=1,ny
      do i=1,nx
c all pressure is in pascals
c set dynamic weight del using lat and surface pressure
c some comments about analysis constants delo and tau
c delo is the inverse square of the expected balance residual. This is
c a specified parameter that is constant over the grid
c tau controls the mass distribution of any continuity adjustments.
c if tau is large mass adjustment occurs over shallow layers
c tau is based on scaling = scale ht**2 Brunt-Vaisala Freq**4/
c   (density**2 gravity**2 mean velocity **2 f**2)
c scale ht is based on cloud depth for cloud adjustment/and or terrain ht
c a best guess is 3000m which works for terrain and clouds in most places
c To represent the  layer of the atm we are considering 850 to 500mb
         
           sumdt=sumdt+(lapstemp(i,j,kf)*(100000./p(kf))**cappa 
     &                - lapstemp(i,j,ks)*(100000./p(ks))**cappa)
           sumdz=sumdz+(lapsphi(i,j,kf)-lapsphi(i,j,ks))/g
           ff = fo*sind(lat(i,j))
           sumf=sumf+ff
         do k=ks,kf
           sumt=sumt+(lapstemp(i,j,k)*(100000./p(k))**cappa)
           sumr=sumr+p(k)/287.04/lapstemp(i,j,k)
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
c simulated cloud
c           sk=(i-nx/2)**2+(k-nz/2)**2+(j-ny/2)**2       
c           omo(i,j,k)=-1.*exp(-sk/100.)
c
c replace missing cloud vv's with background vv's.
            if(abs(omo(i,j,k)).gt.100.)omo(i,j,k)=omb(i,j,k)
         enddo
      enddo
      enddo
      snxny=float(nx*ny)
      sk=float(kf-ks+1)
      sumt=sumt/snxny/sk
      sumdt=sumdt/snxny
      sumdz=sumdz/snxny
      sumf=sumf/snxny
      sumv2=sumv2/snxny/sk
      sumr=sumr/snxny/sk
      ro=sqrt(sumv2)/(sumf*sl)  ! rossby number for dynamic adjustment
      tau=terscl**2*(g*sumdt/sumt/sumdz)**2/(sumr**2*g**2*sumv2*sumf**2)
      print*,'dthet/thet/dz/den/V/f/tau,ro: ',sumdt,sumt,sumdz,sumr,
     &     sqrt(sumv2),sumf,tau,ro
      if (ro.gt.1) ro=1. 
c
c *** Compute non-linear terms (nu,nv) and execute mass/wind balance.
c *** Do for lmax iterations.
c
      lmax=3 
c put lapsphi into phi, lapsu into u, etc
      call move_3d(lapsphi,phi,nx,ny,nz)
      call move_3d(lapsu,u,nx,ny,nz)     
      call move_3d(lapsv,v,nx,ny,nz)    
      call move_3d(lapstemp,t,nx,ny,nz)
      call move_3d(omo,om,nx,ny,nz)
      call move_3d(lapsrh,rh,nx,ny,nz)
c
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
      call balstagger(u,v,phi,t,rh,om,us,vs,
     &phis,ts,rhs,oms,nx,ny,nz,p,ps,1) 

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
      call balstagger(ub,vb,phib,tb,rhb,omb,ubs,vbs,
     &phibs,tbs,rhbs,ombs,nx,ny,nz,p,ps,1) 
   
c put phis into phi, us into u ... ie., put staggered arrays into non-staggered.
      call move_3d(phis,phi,nx,ny,nz)
      call move_3d(us,u,nx,ny,nz)
      call move_3d(vs,v,nx,ny,nz)
      call move_3d(oms,om,nx,ny,nz)
      call move_3d(rhs,rh,nx,ny,nz)
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
c ****** Compute non-linear terms (nu,nv) from observed field.
c
      call nonlin(nu,nv,us,vs,ubs,vbs,oms,ombs
     .           ,nx,ny,nz,dx,dy,dp,dt,bnd)
      call frict(fu,fv,ubs,vbs,us,vs,p,ps,ts
     .                 ,nx,ny,nz,dx,dy,dp,dt,bnd)

c
c ****** Compute momentum residual from obs fields over whole domain.

      call momres(us,vs,phis,nu,nv,fu,fv,wb,delo
     .           ,nx,ny,nz,lat,dx,dy,ps,p)
c
c ****** Execute mass/wind balance.
c
      err=1.0
c returns staggered grids of full fields u,v,phi

      call balcon(phis,us,vs,oms,phi,u,v,om,phibs,ubs,vbs,ombs,
     .         ro,delo,tau,itmax,err,erru,errphi,errub,errphib
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
      call momres(u,v,phi,nu,nv,fu,fv,wb,delo,nx,ny,nz,lat,dx,dy,ps,p)
c
c *** destagger and Write out new laps fields.
c
c the non-staggered grids must be input with intact boundaries from background 

      call balstagger(ub,vb,phib,tb,
     & rhb,omb,u,v,phi,t,rh,om,nx,ny,nz,p,ps,-1) 
c
c   Prior to applying boundary subroutine put non-staggered grids back into
c   u,v,om,t,rh,phi.

      call move_3d(phib,phi,nx,ny,nz)
      call move_3d(ub,u,nx,ny,nz)
      call move_3d(vb,v,nx,ny,nz)
      call move_3d(tb,t,nx,ny,nz)
      call move_3d(omb,om,nx,ny,nz)
      call move_3d(rhb,rh,nx,ny,nz)

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


      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
        if(lapslwc(i,j,k).gt.0. .AND.
     .        lapslwc(i,j,k).lt.100.)
     .                rh(i,j,k)=100.
      enddo
      enddo
      enddo

c
      call write_bal_laps(masstime,phi,u,v,t,om,rh,nx,ny,nz
     .                   ,p,istatus)
      if(istatus.ne.1)then
         write(6,*)'error writing balance fields'
         return
      endif

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
1000  format(1x,'BEFORE/AFTER BALCON.....MOMENTUM RESIDUAL FOR DOMAIN '
     &     ,e12.4,' delo= ',e12.4)
      return
      end
      subroutine diagnose(a,nx,ny,nz,ii,jj,kk,ispan,title)
      implicit none
      real a(nx,ny,nz)
      integer nx,ny,nz
      integer ii,jj,kk
      integer ispan,ieast,iwest,jnorth,jsouth
      integer i,j,k
      character*(*) title
      print*,title,' at level ',kk
      write(9,1002) title,' at level ',kk
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
       write(9,1001)  j,ieast,iwest
 1001  format(1x,3i8)
       print*,j,ieast,iwest
       write(6,1000) (a(i,j,kk),i=ieast,iwest)
       write(9,1000) (a(i,j,kk),i=ieast,iwest)
 1000  format(1x,7e11.5)
      enddo
      return
      end
c     
c===============================================================================
c
c
      subroutine balcon(to,uo,vo,omo,t,u,v,om,tb,ub,vb,omb,
     .   ro,delo,tau,itmax,err,erru,errph,errub,errphb,nu,nv
     .     ,fu,fv,nx,ny,nz,lat,dx,dy,ps,p,dp,lmax)
c
c *** Balcon executes the mass/wind balance computations as described
c        mcginley (Meteor and Appl Phys, 1987) except that
c        this scheme operates on perturbations from background "b"
c        fields. The dynamic constraint is formulated from this 
c        perturbation field. The constraint equation is
c          du'/dt= -ro*nonlin'-d phi'/dx +fv' + friction' 
c        ro is a measure of how well the observation field determines
c        non linear structure. ro is similar to the rossby number
c        V/fL where L is the resolved wavelength and V is a representative
c        velocity.  If obs were everywhere ro =1;
c        very sparse obs ro ~ 0. If data doesn't support it the 
c        adjustment perturbation will be geostrophic.
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
     .      ,nu(nx,ny,nz),nv(nx,ny,nz)
     .      ,fu(nx,ny,nz),fv(nx,ny,nz)
     .      ,dx(nx,ny),dy(nx,ny),dp(nz)
     .      ,ps(nx,ny),p(nz)
     .      ,lat(nx,ny),ff(nx,ny)
     .      ,erru(nx,ny,nz),errph(nx,ny,nz)
     .      ,errub(nx,ny,nz),errphb(nx,ny,nz)

      real*4 tau,err,rdpdg,bnd,g,fo,r,re,ovr
     .      ,ang,f,cotmax,sin1,fs,cos1,beta
     .      ,a,bb,cortmt,ro
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
     .      ,fu2,fv2,fuangu,fvangv

c these are used for diagnostics
      integer nf
      parameter (nf=10)
      real*4  data(nf)
      real*4  fldmax(nf,nz),fldmin(nf,nz)
      integer fldmxi(nf,nz),fldmxj(nf,nz)
      integer fldmni(nf,nz),fldmnj(nf,nz)

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

        call diagnose(to,nx,ny,nz,26,3,7,7,'INPUT    GEOPOTENTIALS')
        call diagnose(uo,nx,ny,nz,26,3,7,7,'INPUT    U-COMPONENT  ')
        call diagnose(vo,nx,ny,nz,26,3,7,7,'INPUT    V-COMPONENT  ')
c create perturbations
c owing to an artifact of coding the t array is phi
      do l=1,lmax
       write(6,*) '|||||||||BALCON ITERATION NUMBER ',l,' ||||||||||'
       print*,'-----------------------------------------------------'
c      write(9,*) '|||||||||BALCON ITERATION NUMBER ',l,' ||||||||||'
       do j=1,ny
        do i=1,nx
         do k=1,nz
          to(i,j,k)=to(i,j,k)-tb(i,j,k)!background an obs are in GPM
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

c *** Compute new phi (t array) using relaxation on eqn. (2).
c        beta*dldx term is dropped to eliminate coupling with lambda eqn.
c
       do 1 it=1,itmax

         cotmax=0.
         do 2 j=2,nym1
          do 2 i=2,nxm1
            
c commented by JS 01-16-01 (calcs only needed once)
c          dxx(i,j)=(dx(i,j)+dx(i,j+1)+dx(i+1,j)+dx(i+1,j+1))*.25
c          dx2(i,j)=dxx(i,j)*2.
c          dxs(i,j)=dxx(i,j)*dxx(i,j)
c          dyy(i,j)=(dy(i,j)+dy(i,j+1)+dy(i+1,j)+dy(i+1,j+1))*.25
c          dy2(i,j)=dyy(i,j)*2.
c          dys(i,j)=dyy(i,j)*dyy(i,j) 
c          fx(i,j)=(ff(i+1,j)-ff(i-1,j))/dx2
c          fy(i,j)=(ff(i,j+1)-ff(i,j-1))/dy2
c          ffx(i,j)=ff(i,j)*fx(i,j)
c          ffy(i,j)=ff(i,j)*fy(i,j)

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
     &            +ro*(term8*snu+term9*snv)+term10*fuu+term11*fvv
     &                -ro*(dnudx+dnvdy)+dfudx+dfvdy
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

2          continue
c        write(6,1000) it,itt,cotmax,ovr,cotma1
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
1000     format(1x,'PHI SOLVER: it = ',i4,' max correction for '
     &          ,a4,' = ',e12.3
     .         ,'ovr =  ',e12.4/1x
     .         ,'first iteration max correction was ',e12.4)

1      continue  !this is the itmax loop

12     write(6,1000) it,itt,cotmax,ovr,cotma1

c     write(9,1000) it,itt,cotmax,ovr,cotma1
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
           u(i,j,k)=(uot*erru(i,j,k)-fuangu*delo*(ro*nvv-fvv+dtdy))
     &                        /aaa(i,j,k)      

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
           v(i,j,k)=(vot*erru(i,j,k)+fvangv*delo*(ro*nuu-fuu+dtdx))
     &                              /aaa(i,j,k)       

           data(5)=nuu
           data(6)=fuu
           data(7)=dtdx
           data(8)=erru(i,j,k)
           call savemxmninfo(5,8,nf,nz,data,i,j,k
     &,fldmax,fldmin,fldmxi,fldmxj,fldmni,fldmnj)

          else
           v(i,j,k)=bnd
          endif
         enddo
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

       enddo

       call printmxmn(8,nf,nz,p,fldmax,fldmin
     &,fldmxi,fldmxj,fldmni,fldmnj)
c 
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

        call diagnose(t,nx,ny,nz,26,3,7,7,'BALANCED GEOPOTENTIALS')
        call diagnose(u,nx,ny,nz,26,3,7,7,'BALANCED U-COMPONENT  ')
        call diagnose(v,nx,ny,nz,26,3,7,7,'BALANCED V-COMPONENT  ')
       if(.true.)then
        do k=1,nz
          print*
          print*,'Calling array diagnosis: ',k,p(k)
          print*,'-------------------------------'
          call array_diagnosis(u(1,1,k),nx,ny,' u-comp   ')
          call array_diagnosis(v(1,1,k),nx,ny,' v-comp   ')
c          call array_diagnosis(uo(1,1,k),nx,ny,' uo comp  ')
c          call array_diagnosis(vo(1,1,k),nx,ny,' vo comp  ')
          call array_diagnosis(t(1,1,k),nx,ny,'  phi     ')
        enddo
       endif


       erf=100.

       call leib_sub(nx,ny,nz,erf,tau,delo
     .,lat,dx,dy,ps,p,dp,t,to,uo,u,vo,v,om,omo,nu,nv,fu,fv)

c move adjusted fields to observation-driven fields for next iteration
       call move_3d(t,to,nx,ny,nz)
       call move_3d(u,uo,nx,ny,nz)
       call move_3d(v,vo,nx,ny,nz)
       call move_3d(om,omo,nx,ny,nz)

      enddo ! on lmax

      deallocate (aaa,bbb)
      deallocate (dxx,dx2,dxs,dyy,dy2,dys
     1,fx,fy,ffx,ffy)


      return
      end
c
c ---------------------------------------------------------------
c
      subroutine leib_sub(nx,ny,nz,erf,tau,delo
     .,lat,dx,dy,ps,p,dp,t,to,uo,u,vo,v,om,omo,nu,nv,fu,fv)

      implicit none

      integer nx,ny,nz
      integer nxm1,nym1,nzm1
      integer i,j,k,ks

      real*4 t(nx,ny,nz),to(nx,ny,nz)
     .      ,u(nx,ny,nz),uo(nx,ny,nz)
     .      ,v(nx,ny,nz),vo(nx,ny,nz)
     .      ,om(nx,ny,nz),omo(nx,ny,nz)
     .      ,nu(nx,ny,nz),nv(nx,ny,nz)
     .      ,fu(nx,ny,nz),fv(nx,ny,nz)
     .      ,lat(nx,ny),dx(nx,ny),dy(nx,ny)
     .      ,ps(nx,ny),p(nz),dp(nz)

      real*4 ang,rdpdg,sin1,dldx,dldy,dldp
     .,a,f,fo,fs,erf,tau,delo,bnd

      real, allocatable, dimension(:,:,:) :: slam,f3,h

      fo=14.52e-5   !2*omega
      rdpdg=3.141592654/180.
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
      call fthree(f3,u,v,omo,delo,nu,nv,h,tau,
     .   nx,ny,nz,lat,dx,dy,dp)
c
c ****** Perform 3-d relaxation.
c
      erf=.1
      call leibp3(slam,f3,200,erf,h
     .  ,nx,ny,nz,dx,dy,ps,p,dp)
c
c ****** Compute new u, v, omega by adding the lagrange multiplier terms.
c
      do k=1,nz
      ks=1
      if(k.eq.nz) ks=0
      do j=1,nym1
      do i=1,nxm1

c        ang=(lat(i,j))*rdpdg
c        sin1=sin(ang)
c        f=fo*sin1
c        fs=f*f
c        a=1.+fs*delo
c commemnted 01-02-01
         a=1.0
         dldp=(slam(i,j,k)-slam(i,j,k+1))/dp(k+ks)
         dldx=(slam(i+1,j+1,k+1)-slam(i,j+1,k+1))/dx(i,j)
         dldy=(slam(i+1,j+1,k+1)-slam(i+1,j,k+1))/dy(i,j)
         if (u(i,j,k) .ne. bnd) u(i,j,k)=u(i,j,k)+.5*dldx/a
         if (v(i,j,k) .ne. bnd) v(i,j,k)=v(i,j,k)+.5*dldy/a
         if (omo(i,j,k).ne.bnd) om(i,j,k)=omo(i,j,k)+.5*dldp/tau
      enddo
      enddo
      enddo
      call diagnose(slam,nx+1,ny+1,nz+1,27,5,8,8,'LAMDA ')
      call diagnose(u,nx,ny,nz,26,3,7,7,'CONTIN U-COMPONENT  ')
      call diagnose(v,nx,ny,nz,26,3,7,7,'CONTIN V-COMPONENT  ')
      call diagnose(om,nx,ny,nz,27,4,7,7,'CONTIN OMEGA    ')
      call diagnose(om,nx,ny,nz,27,4,8,7,'CONTIN OMEGA    ')
      call analz(t,to,u,uo,v,vo,om,omo
     .                ,slam,f3,nu,nv,fu,fv,delo,tau
     .                ,nx,ny,nz
     .                ,lat,dx,dy,ps,p,dp)


      deallocate (slam,f3,h)

      return
      end
c
c ---------------------------------------------------------------
c
      subroutine analz(t,to,u,uo,v,vo,om,omo,slam,f3,nu,nv,fu,fv,
     .               delo,tau ,nx,ny,nz,lat,dx,dy,ps,p,dp)
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
     .      ,delo,tau,rdpdg,fo,bnd
     .      ,conmax,sumom,cont,con,sumu,sumv,sumt,thermu,thermv
     .      ,sumww,resu,resv,tgpu,tgpv,tgpph,tgpc
     .      ,ang,sin1,f,nvv,nuu,dtdx,dtdy,dudx,dvdy,domdp
     .      ,uot,vot,tot,uuu,vvv,contm,fuu,fvv,dxx,dyy
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
               resv=(f*u(i,j,k)+dtdy+nvv-fvv)**2+resv
               tgpu=tgpu+1.
            endif
            if (vot .ne. bnd) then
               thermv=(v(i,j,k)-dtdx/f)**2+thermv
               sumv=sumv+(v(i,j,k)-vot)**2
               resu=resu+(-f*v(i,j,k)+dtdx+nuu-fuu)**2
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

         write(6,1001) sumt,delo,sumu,sumv,
     .             sumom,tau,cont,thermu,thermv,sumww,resu,resv

c        write(9,1001) sumt,delo,sumu,sumv,
c    .             sumom,tau,cont,thermu,thermv,sumww,resu,resv

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
     . /1x,'    om-omo:  ',e12.5,'   tau:',e12.5,
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
c *** Frict computes frictional acceleration in the near the terrain surface
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
      hzdf=50000.!m2/sec
      eddf=20.!pa2/sec
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
         dauvbdp=(sqrt(ub(i,j,k-1)**2+vb(i,j,k-1)**2)-
     &           sqrt(ub(i,j,k+1)**2+vb(i,j,k+1)**2))/ddp(k)
         absuv=sqrt(u(i,j,k)**2+v(i,j,k)**2)
         absuvb=sqrt(ub(i,j,k)**2+vb(i,j,k)**2)
         if(u(i,j,k).eq.bnd) then
          fu(i,j,k)=0.
         else
          d2udx=(u(i+1,j,k)+u(i-1,j,k)-2.*u(i,j,k))/dx2(i,j)
          d2udy=(u(i,j+1,k)+u(i,j-1,k)-2.*u(i,j,k))/dy2(i,j)
          d2udp=4.*(u(i,j,k+1)+u(i,j,k-1)-2.*u(i,j,k))/ddp2(k)
          d2ubdx=(ub(i+1,j,k)+ub(i-1,j,k)-2.*ub(i,j,k))/dx2(i,j)
          d2ubdy=(ub(i,j+1,k)+ub(i,j-1,k)-2.*ub(i,j,k))/dy2(i,j)
          d2ubdp=4.*(ub(i,j,k+1)+ub(i,j,k-1)-2.*ub(i,j,k))/ddp2(k)
          gdtudp=0.
          if(ps(i,j).gt.p(k+1).and.(ps(i,j)-p(k+1)).le.dp(k))then! we're near sfc
           den=p(k+1)/r/t(i,j,k+1)
           dudp=(u(i,j,k-1)-u(i,j,k+1))/ddp(k)
           dubdp=(ub(i,j,k-1)-ub(i,j,k+1))/ddp(k)
           gdtudp=den*cd*grav*(dudp*absuv+dauvdp*u(i,j,k)
     &        -dubdp*absuvb-dauvbdp*ub(i,j,k))
           gdtudp_save(i,j)=gdtudp
          endif
          fu(i,j,k)=hzdf*(d2udx+d2udy-d2ubdx-d2ubdy)
     &         +eddf*(d2udp-d2ubdp)-gdtudp
         endif
         if(v(i,j,k).eq.bnd) then
          fv(i,j,k)=0.
         else
          d2vdx=(v(i+1,j,k)+v(i-1,j,k)-2.*v(i,j,k))/dx2(i,j)
          d2vdy=(v(i,j+1,k)+v(i,j-1,k)-2.*v(i,j,k))/dy2(i,j)
          d2vdp=4.*(v(i,j,k+1)+v(i,j,k-1)-2.*v(i,j,k))/ddp2(k)
          d2vbdx=(vb(i+1,j,k)+vb(i-1,j,k)-2.*vb(i,j,k))/dx2(i,j)
          d2vbdy=(vb(i,j+1,k)+vb(i,j-1,k)-2.*vb(i,j,k))/dy2(i,j)
          d2vbdp=4.*(vb(i,j,k+1)+vb(i,j,k-1)-2.*vb(i,j,k))/ddp2(k)
          gdtvdp=0.
          if(ps(i,j).gt.p(k+1).and.(ps(i,j)-p(k+1)).le.dp(k)) then! we're near sfc
           den=p(k+1)/r/t(i,j,k+1)
           dvdp=(v(i,j,k-1)-v(i,j,k+1))/ddp(k)
           davdp=(abs(v(i,j,k-1))-abs(v(i,j,k+1)))/ddp(k)
           dvbdp=(vb(i,j,k-1)-vb(i,j,k+1))/ddp(k)
           davbdp=(abs(vb(i,j,k-1))-abs(vb(i,j,k+1)))/ddp(k)
           gdtvdp=den*cd*grav*(dvdp*absuv+dauvdp*v(i,j,k)
     &             -dvbdp*absuvb-dauvbdp*vb(i,j,k))
           gdtvdp_save(i,j)=gdtvdp
           if(abs(gdtvdp).gt.0.5)then
              print*,'abs(gdtvdp) > 0.5 :', gdtvdp
              print*,'den/dvdp/absuv/dauvdp/absuvb/dvbdp/dauvbdp'
              print*,den,dvdp,absuv,dauvdp,absuvb,dvbdp,dauvbdp
           endif
          endif
          fv(i,j,k)=hzdf*(d2vdx+d2vdy-d2vbdx-d2vbdy)
     &              +eddf*(d2vdp-d2vbdp)-gdtvdp
         endif
      enddo
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

c     do j=1,ny,2
c        print*,'j ',j
c        print*,'i/gdtudp/gdtvdp ',(i,' ',gdtudp_save(i,j)
c    &,' ',gdtvdp_save(i,j),': ',i=1,nx,2)
c     enddo

c     print*,'j = 1'
c     print*,'i/gdtudp/gdtvdp ',(i,' ',gdtudp_save(i,1)
c    &,' ',gdtvdp_save(i,1),': ',i=1,nx,2)
c     print*
c     print*,'j = ',ny
c     print*,'i/gdtudp/gdtvdp ',(i,' ',gdtudp_save(i,ny)
c    &,' ',gdtvdp_save(i,ny),': ',i=1,nx,2)
c     print*
c     print*,'i = 1'
c     print*,'j/gdtudp/gdtvdp ',(j,' ',gdtudp_save(1,j)
c    &,' ',gdtvdp_save(1,j),': ',j=1,ny,2)
c     print*
c     print*,'i = ',nx
c     print*,'j/gdtudp/gdtvdp ',(j,' ',gdtudp_save(nx,j)
c    &,' ',gdtvdp_save(nx,j),': ',j=1,ny,2)


      return
      end
c===============================================================================
c
      subroutine nonlin(nu,nv,u,v,ub,vb,om,omb
     .                 ,nx,ny,nz,dx,dy,dp,dt,bnd)
c
c *** Nonlin computes the non-linear terms (nu,nv) from staggered input.
c     The non linear terms are linearized with a background/perturbation
c     combination. The eularian time terms are assumed to be zero.
c     The non-linear terms nu,nv, are computed on the 
c     u, v grids, respectively
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
     .      ,omu,ombu,omv,ombv,vvu,uuv,vvbu,uubv
c_______________________________________________________________________________
c
c
      print *,'nonlin'
      nxm1=nx-1
      nym1=ny-1
      nzm1=nz-1
      call zero3d(nu,nx,ny,nz)
      call zero3d(nv,nx,ny,nz)
c create perturbations in the u,v,om variables
      do k=1,nz
       do j=1,ny
        do i=1,nx
         if(u(i,j,k).ne.bnd) u(i,j,k)=u(i,j,k)-ub(i,j,k)
         if(v(i,j,k).ne.bnd) v(i,j,k)=v(i,j,k)-vb(i,j,k)
         if(om(i,j,k).ne.bnd) om(i,j,k)=om(i,j,k)-omb(i,j,k)
        enddo
       enddo  
      enddo 
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
          nu(i,j,k)=dudt+ub(i,j,k)*dudx+vvbu*dudy+ombu*dudp
     &                 +u(i,j,k)*dubdx+vvu*dubdy+omu*dubdp
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
          nv(i,j,k)=dvdt+uubv*dvdx+vb(i,j,k)*dvdy+ombv*dvdp
     &                 +uuv*dvbdx+v(i,j,k)*dvbdy+omv*dvbdp
         endif
      enddo
      enddo
      enddo
c make winds whole again
      do k=1,nz
       do j=1,ny
        do i=1,nx
         if(u(i,j,k).ne.bnd) u(i,j,k)=u(i,j,k)+ub(i,j,k)
         if(v(i,j,k).ne.bnd) v(i,j,k)=v(i,j,k)+vb(i,j,k)
         if(om(i,j,k).ne.bnd) om(i,j,k)=om(i,j,k)+omb(i,j,k)
        enddo
       enddo  
      enddo 
c
      return
      end
c
c===============================================================================
c
      subroutine fthree(f3,u,v,om,delo,nu,nv,h,tau
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
     .      ,u(nx,ny,nz),v(nx,ny,nz)
     .      ,om(nx,ny,nz),lat(nx,ny)
     .      ,nu(nx,ny,nz),nv(nx,ny,nz)
     .      ,dx(nx,ny),dy(nx,ny),dp(nz)
     .      ,delo,tau,fo,rdpdg,dpp,f,aa,dxx,dyy
     .      ,dnudy,dnvdx,cont,formax
c_______________________________________________________________________________
c
      print *,'fthree'
      fo=14.52e-5   !2*omega
      rdpdg=3.141592654/180.
      formax=0.
      do k=2,nz
         dpp=dp(k)
         do j=2,ny
         do i=2,nx
c           f=sin(lat(i,j)*rdpdg)*fo
c           aa=1.+f*f*delo 
            aa=1.0
            js=1
            if (j .eq. 2) js=0
            dyy=dy(i,j)*float(js+1)
            h(i,j,k)=aa/tau
            is=1
            if (i .eq. 2) is=0

c           dxx=dx(i,j)*float(is+1)
c           dnudy=(nu(i,j,k-1)+nu(i-1,j,k-1)
c    .            -nu(i,j-js-1,k-1)-nu(i-1,j-js-1,k-1))/dyy*.25
c           dnvdx=(nv(i,j,k-1)+nv(i,j-1,k-1)
c    .            -nv(i-is-1,j,k-1)-nv(i-is-1,j-1,k-1))/dxx*.25

            cont=((u(i,j-1,k-1)-u(i-1,j-1,k-1))/dx(i,j)
     .           +(v(i-1,j,k-1)-v(i-1,j-1,k-1))/dy(i,j)
     .           +(om(i,j,k-1)-om(i,j,k))/dpp)*aa
            f3(i,j,k)=-2.*cont! remove this term test +2.*f*delo*(dnvdx-dnudy)
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
      Subroutine balstagger(u,v,phi,t,rh,om,
     &  us,vs,phis,ts,rhs,oms,
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
     &rhs(nx,ny,nz),rh(nx,ny,nz)

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
          wr(i,j,k,5)=(rh(i,j,k)+rh(i,j,k+1))*.5
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
          rhs(i,j,k)=(wr(i,j,k,5)+wr(i+1,j+1,k,5)+wr(i,j+1,k,5)
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
c re compute hydrostatic virtual t from adjusted phis 
       do k=2,nz-1
        do j=1,ny
         do i=1,nx
          ts(i,j,k)=(phis(i,j,k)-phis(i,j,k-1))/
     &                              (alog(p(k)/p(k+1))*r)
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
          rh(i,j,k)=(rhs(i-1,j-1,k-1)+rhs(i,j-1,k-1)+
     &               rhs(i-1,j,k-1)+rhs(i,j,k-1)+
     &               rhs(i-1,j-1,k)+rhs(i,j-1,k)+
     &               rhs(i-1,j,k-1)+rhs(i,j,k))*.125
          phi(i,j,k)=(phis(i-1,j-1,k-1)+phis(i,j-1,k-1)+
     &               phis(i-1,j,k-1)+phis(i,j,k-1))*.25
          t(i,j,k)=(ts(i-1,j-1,k-1)+ts(i,j-1,k-1)+
     &               ts(i-1,j,k-1)+ts(i,j,k-1)+
     &               ts(i-1,j-1,k)+ts(i,j-1,k)+
     &               ts(i-1,j,k)+ts(i,j,k))*.125
c  devirtualize temperature
          tdm=dwpt(t(i,j,k)-273.16,rh(i,j,k))
          tvkm=tv(t(i,j,k)-273.16,tdm,p(k)/100.)
          tvkd=tv(t(i,j,k)-273.16,tdd,p(k)/100.)
          t(i,j,k)=t(i,j,k)-(tvkm-tvkd)
         enddo
        enddo
       enddo
c final step is bring temps in line with heights on the west and
c south boundaries: use taylor extrapolation
c north and east boundaries for omega
       do k=1,nz
        do j=2,ny
         if(ps(1,j).ge.p(k)) then
          u(1,j,k)=3.*u(2,j,k)-3.*u(3,j,k)+u(4,j,k)
          v(1,j,k)=3.*v(2,j,k)-3.*v(3,j,k)+v(4,j,k)
          om(1,j,k)=3.*om(2,j,k)-3.*om(3,j,k)+om(4,j,k)
         else
          u(1,j,k)=bnd
          v(1,j,k)=bnd
          om(1,j,k)=bnd
         endif
         if(ps(nx,j).ge.p(k)) then
          om(nx,j,k)=3.*om(nx-1,j,k)-3.*om(nx-2,j,k)+om(nx-3,j,k)
         else
          om(nx,j,k)=bnd
         endif
         t(1,j,k)=3.*t(2,j,k)-3.*t(3,j,k)+t(4,j,k)
         phi(1,j,k)=3.*phi(2,j,k)-3.*phi(3,j,k)+phi(4,j,k)
         rh(1,j,k)=3.*rh(2,j,k)-3.*rh(3,j,k)+rh(4,j,k)
        enddo
        do i=2,nx
         if(ps(i,1).ge.p(k)) then
          u(i,1,k)=3.*u(i,2,k)-3.*u(i,3,k)+u(i,4,k)
          v(i,1,k)=3.*v(i,2,k)-3.*v(i,3,k)+v(i,4,k)
          om(i,1,k)=3.*om(i,2,k)-3.*om(i,3,k)+om(i,4,k)
         else
          u(i,1,k)=bnd
          v(i,1,k)=bnd
          om(i,1,k)=bnd
         endif
         if(ps(i,ny).ge.p(k)) then
          om(i,ny,k)=3.*om(i,ny-1,k)-3.*om(i,ny-2,k)+om(i,ny-3,k)
         else
          om(i,ny,k)=bnd
         endif
         t(i,1,k)=3.*t(i,2,k)-3.*t(i,3,k)+t(i,4,k)
         phi(i,1,k)=3.*phi(i,2,k)-3.*phi(i,3,k)+phi(i,4,k)
         rh(i,1,k)=3.*rh(i,2,k)-3.*rh(i,3,k)+rh(i,4,k)
        enddo
        if(ps(1,1).ge.p(k)) then
          u(1,1,k)=3.*u(2,2,k)-3.*u(3,3,k)+u(4,4,k)
          v(1,1,k)=3.*v(2,2,k)-3.*v(3,3,k)+v(4,4,k)
          om(1,1,k)=3.*om(2,2,k)-3.*om(3,3,k)+om(4,4,k)
        else
          u(1,1,k)=bnd
          v(1,1,k)=bnd
          om(1,1,k)=bnd
        endif
        t(1,1,k)=3.*t(2,2,k)-3.*t(3,3,k)+t(4,4,k)
        phi(1,1,k)=3.*phi(2,2,k)-3.*phi(3,3,k)+phi(4,4,k)
        rh(1,1,k)=3.*rh(2,2,k)-3.*rh(3,3,k)+rh(4,4,k)
       enddo 
      endif ! destagger

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

