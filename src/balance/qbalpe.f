      program qbalpe_main
c
      implicit none
c
      integer   nx,ny,nz
      integer   istatus
c_______________________________________________________________________________
c

      call get_laps_config('nest7grid',istatus)

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

      call qbalpe(nx,ny,nz)
c

999   print*,'Done'

1000  end
c
c===============================================================================
c
      subroutine qbalpe(nx,ny,nz)
c
      include 'trigd.inc'
      implicit none
c
      integer*4 nx,ny,nz
c
      real*4 dx(nx,ny),dy(nx,ny),dp(nz)
     .      ,p(nz),ps(nx,ny),ter(nx,ny)
     .      ,lat(nx,ny),lon(nx,ny)
     .      ,phi(nx,ny,nz),t(nx,ny,nz)
     .      ,u(nx,ny,nz),v(nx,ny,nz)
     .      ,lapsuo(nx,ny,nz),lapsvo(nx,ny,nz) !t=t0
     .      ,lapsu(nx,ny,nz),lapsv(nx,ny,nz)   !t=t0-dt
     .      ,dir(nx,ny),spd(nx,ny)
     .      ,temp(nx,ny,nz)
     .      ,lapsphi(nx,ny,nz)
     .      ,laps3d(nx,ny,nz,2)
     .      ,om(nx,ny,nz)
     .      ,omo(nx,ny,nz),oms(nx,ny,nz)
     .      ,nu(nx,ny,nz),nv(nx,ny,nz)
     .      ,uo(nx,ny,nz)
     .      ,vo(nx,ny,nz)
     .      ,xx(nx,ny,nz),xxx(nx,ny,nz)
     .      ,wb(nx,ny,nz)
     .      ,wr1(nx,ny,nz),wr2(nx,ny,nz),wr3(nx,ny,nz)
     .      ,re,rdpdg,po,cappa,gam(nz),dppp(nz)
     .      ,delo,gamo,del(nx,ny,nz),tau,dt
     .      ,err,grid_spacing_actual_m
     .      ,pdif,dpbl,dpblf
     .      ,u_grid,v_grid
     .      ,u_true,v_true
c
      integer*4 ip(nz),icon,itmax,lmax
     .         ,masstime,windtime,sfctime,omtime
     .         ,i,j,k,ll,istatus

      integer*4 lend
      integer*4 lends
      integer*4 lendt
      integer*4 lendw
c
      logical lrunbal
      logical lstagger
      logical lnon_linear

      character*255 staticdir,tempdir,winddir,sfcdir
      character*125 comment
      character*31  staticext,tempext,windext,sfcext
      character*10  units
      character*9   a9_time
c
c     data p/1100.,1050.,1000., 950., 900., 850., 800.
c    .      , 750., 700., 650., 600., 550., 500., 450.
c    .      , 400., 350., 300., 250., 200., 150., 100./
c     data dp/0.,20*50./
c_______________________________________________________________________________
c
      include 'lapsparms.cmn'
      
      call get_balance_nl(lrunbal,lstagger,icon,gamo,delo,tau
     1,lnon_linear,istatus)
      if(istatus.ne.0)then
         print*,'error getting balance namelist'
         stop
      endif

      print*,'lrunbal = ',lrunbal
      print*,'lstagger = ',lstagger
      print*,'icon = ',icon
      print*,'gamo = ',gamo
      print*,'delo = ',delo
      print*,'tau  = ',tau
      print*,'lnon_linear = ',lnon_linear
c
c switch to run balance package or not
c
      if(.not.lrunbal)then
         print*,'Namelist value lrunbal = false '
         print*,'Balance Package not running '
         goto 999
      endif

      if(vertical_grid.eq.'PRESSURE')THEN!Pressure in mb
         do i=1,nz
          p(i)=pressure_bottom_l/100.-((i-1)*(pressure_interval_l/100.))
          if(i.gt.1)dp(i)=pressure_interval_l/100.
         enddo
      else
         print*,'vertical grid is not PRESSURE ',vertical_grid
         goto 999
      endif

      staticext='nest7grid'
      tempext='lt1'
      windext='lw3'
      sfcext='lsx'
      call get_directory(staticext,staticdir,lend)
      call get_directory(tempext,tempdir,lendt)
      call get_directory(windext,winddir,lendw)
      call get_directory(sfcext,sfcdir,lends)

      re=6371220.
      rdpdg=3.141592654/180.
      po=p(1)
      cappa=287.053/1004.686
      itmax=200  !max iterations for relaxation
c
c some words on gamo.  gamo is the mid-atmospheric ratio between the
c rms wind adjustment and rms geopotential adjustment.  As you go up
c in height you would want greater potential adjustments of phi relative
c to wind, simply because error in Phi increases upward.    The standard
c ratio for equal adjustment of both winds and height is 1. x 10**-2
c For 10 times more adjustment in phi use .01 x 10**-2.  These numbers are
c related to expected observaional and analysis error of the lt1 and lw3 
c fields. 
c
      do k=1,nz
         gam(k)=gamo
         dppp(k)=dp(k)*100.
         ip(k)=int(p(k))
      enddo
c
c *** Get times of mass, wind and surface data.
c
c     call gettime(masstime)

      call get_systime(masstime,a9_time,istatus)
      if(istatus .ne. 1)go to 999
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
      rdpdg=3.14159/180.
c
      do j=1,ny
      do i=1,nx
         call get_grid_spacing_actual(lat(i,j),lon(i,j)
     1             ,grid_spacing_actual_m,istatus)
         dx(i,j)=grid_spacing_actual_m
         dy(i,j)=dx(i,j)
         
         do k=1,nz
c reduce geostrophic weight to 10% within dpbl of surface 
           dpbl=100.
           dpblf=dpbl*.9
           pdif=(ps(i,j)-p(k))
           if (pdif.le.dpbl) pdif=dpbl  
           del(i,j,k)=abs(delo*(sind(lat(i,j))*(1.- dpblf/pdif)))
         enddo
      enddo
      enddo
c
c *** Get laps heights.
c
      call get_laps_3d(masstime,nx,ny,nz
     1  ,tempext,'ht',units,comment,lapsphi,istatus)

      if (istatus .ne. 1) then
         print *,'Error getting LAPS height data...Abort.'
         stop
      endif
c
c *** Get laps temps
c
      call get_laps_3d(masstime,nx,ny,nz
     1  ,tempext,'t3',units,comment,temp,istatus)

      if (istatus .ne. 1) then
         print *,'Error getting LAPS temp data...Abort.'
         stop
      endif
c
c *** Get laps wind data.
c        Read t=t0 first, then read t=t0-dt.
c *** Not considering non-linear terms for now, so no need to read t0-dt. 
c
c     call get_laps_wind(winddir,windtime,windext,nx,ny,nz
c    .                  ,lapsuo,lapsvo,istatus)

      call get_laps_3d(masstime,nx,ny,nz
     1  ,windext,'u3',units,comment,lapsu,istatus)

      if (istatus .ne. 1) then
         print *,'Error getting LAPS time 0 u3 data...Abort.'
         stop
      endif

      call get_laps_3d(masstime,nx,ny,nz
     1  ,windext,'v3',units,comment,lapsv,istatus)

      if (istatus .ne. 1) then
         print *,'Error getting LAPS time 0 v3 data...Abort.'
         stop
      endif

      write(6,*)' Rotating u/v to grid north'
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         call uvtrue_to_uvgrid(
     1            lapsu(i,j,k),lapsv(i,j,k)
     1           ,u_grid   ,v_grid
     1           ,lon(i,j)           )
         lapsu(i,j,k) = u_grid                    
         lapsv(i,j,k) = v_grid                    
      enddo
      enddo
      enddo

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
      call get_laps_2d_sfc(sfcdir,sfctime,sfcext,nx,ny
     .                ,1,13,ps,istatus)
      do j=1,ny
      do i=1,nx
         ps(i,j)=ps(i,j)*0.01
      enddo
      enddo
c
c *** Read maps qg omega.
c
c     call get_maps_qgom(omtime,omo,istatus)
c
c *** Compute non-linear terms (nu,nv) and execute mass/wind balance.
c *** Do for lmax iterations.
c
      lmax=1
      if(icon.eq.1)lmax=2
      do ll=1,lmax
c
         do k=1,nz
         do j=1,ny
         do i=1,nx
            phi(i,j,k)=lapsphi(i,j,k)
            u(i,j,k)=lapsu(i,j,k)   !t=t0-dt
            v(i,j,k)=lapsv(i,j,k)   !t=t0-dt
            uo(i,j,k)=u(i,j,k)      !t=t0
            vo(i,j,k)=v(i,j,k)      !t=t0
            xx(i,j,k)=phi(i,j,k)
            xxx(i,j,k)=u(i,j,k)
            t(i,j,k)=v(i,j,k)
         enddo
         enddo
         enddo
c stagger the LAPS grids to prepare for balancing
c
         if(lstagger) call balstagger(xxx,t,omo,xx,temp,
     &nx,ny,nz,wr1,wr2,wr3,p,1)

         call terbnd(xxx,t,omo,nx,ny,nz,ps,p)
         if (ll .eq. 1) then
            do k=1,nz
            do j=1,ny
            do i=1,nx
               oms(i,j,k)=omo(i,j,k)
            enddo
            enddo
            enddo
         endif
c
c ****** Compute non-linear terms (nu,nv).
c
         if(lnon_linear) call nonlin(nu,nv,uo,vo,u,v,oms
     .              ,nx,ny,nz,dx,dy,dppp,dt)
c
c ****** Compute momentum residual over whole domain.
c
         if (ll .eq. 1) call momres(xxx,t,phi,nu,nv,wb,0.
     .                             ,nx,ny,nz,lat,dx,dy,ps,p)
c
c ****** Execute mass/wind balance.
c
         err=.1
         call balcon(xx,xxx,t,omo,phi,u,v,oms,del,tau,itmax,err
     .              ,nu,nv,icon,nx,ny,nz,lat,dx,dy,ps,p,dp,gam)
c
      enddo

      call momres(u,v,phi,nu,nv,wb,delo,nx,ny,nz,lat,dx,dy,ps,p)
c
c *** destagger and Write out new laps fields.
c
      if(lstagger) then
         call balstagger(u,v,oms,phi,temp,
     &                   nx,ny,nz,wr1,wr2,wr3,p,-1)
        else
c will give non-staggered temps from new balanced phis
         call phigns(phi,temp,nx,ny,nz,0,p,-1)
      endif

      write(6,*)' Rotating balanced u/v to true north'
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         call uvgrid_to_uvtrue(
     1            u(i,j,k),v(i,j,k)
     1           ,u_true   ,v_true
     1           ,lon(i,j)           )
         u(i,j,k) = u_true
         v(i,j,k) = v_true
      enddo
      enddo
      enddo


c
      call write_bal_laps(masstime,phi,u,v,temp,oms,nx,ny,nz
     .                   ,ip,istatus)
      if(istatus.ne.1)then
         write(6,*)'error writing balance fields'
         return
      endif

999   return
      end

c
c===============================================================================
c
      subroutine momres(u,v,phi,nu,nv,wa,del
     .                 ,nx,ny,nz,lat,dx,dy,ps,p)
c
c *** Momres computes momentum residual for whole domain.
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
     .      ,wa(nx,ny,nz),wb(nx,ny,nz)
     .      ,wc(nx,ny,nz),wd(nx,ny)
     .      ,lat(nx,ny),dx(nx,ny),dy(nx,ny)
     .      ,ps(nx,ny),p(nz)
     .      ,del,g,rdpdg,fo,bnd,cnt,sum,ang,sin1,f
     .      ,nvv,nuu,errms
c_______________________________________________________________________________
c
      print *,'momres'
      g=9.80665
      rdpdg=3.141592654/360.
      fo=14.52e-5   !2*omega
      bnd=1.e-30
      nxm1=nx-1
      nym1=ny-1
      nzm1=nz-1
      cnt=0.
      sum=0.
      call zero(wa,nx,ny,nz)
      call zero(wb,nx,ny,nz)
      call zero(wc,nx,ny,nz)
      do j=1,ny
      do i=1,nx
         wd(i,j)=0
      enddo
      enddo
      do j=2,nym1
      do i=2,nxm1
         ang=(lat(i,j)+lat(i,j+1))*rdpdg
         sin1=sin(ang)
         f=fo*sin1
         do k=1,nz
            ks=1
            if (k .eq. nz) ks=0
            nvv=(nv(i,j,k)+nv(i,j+1,k)
     .          +nv(i-1,j+1,k)+nv(i-1,j,k))*.25
            nuu=(nu(i,j,k)+nu(i+1,j,k)
     .          +nu(i,j-1,k)+nu(i+1,j-1,k))*.25  
c           if (k .eq. 13) wd(i,j)=sqrt(nvv**2+nuu**2)
            if (u(i,j,k) .ne. bnd)
     .          wc(i,j,k)=f*u(i,j,k)+g/dy(i,j)
     .                   *(phi(i,j+1,k)-phi(i,j,k))+nvv
            if (v(i,j,k) .ne. bnd)
     .          wb(i,j,k)=-f*v(i,j,k)+g/dx(i,j)
     .                   *(phi(i+1,j,k)-phi(i,j,k))+nuu
         enddo
      enddo
      enddo
      do k=1,nzm1
      do j=2,nym1
      do i=2,nxm1
         if (ps(i,j) .gt. p(k+ks)) then
            wa(i,j,k)=sqrt((wc(i-1,j-1,k)+wc(i,j-1,k))**2*.25
     .                    +(wb(i-1,j,k)+wb(i-1,j-1,k))**2*.25)
            sum=wa(i,j,k)**2+sum
            cnt=cnt+1.
         endif
      enddo
      enddo
      enddo
c     call prt(wd,sc,ci,nx,ny,1,nz,4hnon ,4hlin )                 
c     do 6 k=3,nz,8                                                   
c     call prt(wa,sc,ci,nx,ny,k,nz,4hmom ,4hres )                 
c    6 continue                                                          
      errms=sqrt(sum/cnt)
      write(6,1000) errms,del
1000  format(1x,'momentum residual for domain ',e12.4,' del= ',e12.4)
      return
      end
c
c===============================================================================
c
      subroutine omega(om,u,v,t,vort,absvort
     .                ,dx,dy,lat,p,dppp,ps,nx,ny,nz)
c
      implicit none
c
      integer*4 nx,ny,nz,nzm1,i,j,k
c
      real*4 dx(nx,ny),dy(nx,ny),dppp(nz)
     .      ,lat(nx,ny),p(nz),ps(nx,ny)
     .      ,u(nx,ny,nz),v(nx,ny,nz)
     .      ,t(nx,ny,nz)
     .      ,vort(nx,ny,nz),absvort(nx,ny,nz)
     .      ,om(nx,ny,nz)
     .      ,a(nx,ny),b(nx,ny)
     .      ,c(nx,ny),f(nx,ny)
     .      ,xxx(nx,ny,nz),xxxx(nx,ny,nz)
     .      ,xxxxx(nx,ny,nz)
     .      ,re,r,cappa,g,rdpdg,fo,stab,delp
     .      ,ang,cos1,sin1,tbar,denbar
     .      ,potm5,potm15
     .      ,beta,dd,twolam,twoa,cee
c_______________________________________________________________________________
c                                           
      print *,'omega'
      re=6371000.
      r=287.053
      cappa=287.053/1004.686
      g=9.80665
      rdpdg=3.141592654/180.
      fo=14.52e-5   !2*omega
      stab=0.
      nzm1=nz-1
      call zero(om,nx,ny,nz)
      call zero(xxxxx,nx,ny,nz)
      call zero(xxxx,nx,ny,nz)
      call zero(xxx,nx,ny,nz)
      call zero(vort,nx,ny,nz)
      delp=(p(5)-p(13))*100.
      do j=1,ny
      do i=1,nx
         ang=lat(i,j)*rdpdg
         cos1=cos(ang)
         sin1=sin(ang)
         tbar=(t(i,j,13)+t(i,j,5))/2.
         denbar=100.*(p(5)+p(13))/r/tbar/2.
         potm5=t(i,j,5)*(1000./p(5))**cappa
         potm15=t(i,j,13)*(1000./p(13))**cappa
         stab=(alog(potm15)-alog(potm5))/denbar/delp+stab
         f(i,j)=(fo*sin1)                             
         b(i,j)=0.
         c(i,j)=-sin1/cos1/re
         do k=1,nz
            if(i.eq.1.or.j.eq.1) then
               om(i,j,k)=0.
               absvort(i,j,k)=f(i,j)
            else
               xxxx(i,j,k)=(-u(i-1,j-1,k)+u(i,j-1,k))/dx(i,j)
     .                    -(v(i-1,j,k)-v(i-1,j-1,k))/dy(i,j)
               xxx(i,j,k)=(v(i,j,k)-v(i-1,j,k))/dx(i,j)
     .                   +(u(i,j,k)-u(i,j-1,k))/dy(i,j)
               vort(i,j,k)=(v(i,j,k)-v(i-1,j,k))/dx(i,j)
     .                    -(u(i,j,k)-u(i,j-1,k))/dy(i,j)
               absvort(i,j,k)=f(i,j)+vort(i,j,k)
            endif
         enddo
      enddo
      enddo
      stab=stab/(float(nx)*float(ny))
      write(6,1001) stab
1001  format(1x,' stability is ',e12.4)
      do j=1,ny
      do i=1,nx
         a(i,j)=f(i,j)**2/stab
      enddo
      enddo
      do k=2,nzm1
      do j=2,ny
      do i=2,nx
         ang=lat(i,j)*rdpdg
         cos1=cos(ang)
         beta=fo*cos1*rdpdg/dy(i,j)
         dd=sqrt(dx(i,j)**2+dy(i,j)**2)
         twolam=(xxxx(i,j,k)+xxxx(i,j,k-1))/2.*(xxx(i,j,k-1)+
     .           xxx(i,j-1,k-1)+xxx(i-1,j-1,k-1)+xxx(i-1,j,k-1)-
     .           (xxx(i,j,k)+xxx(i,j-1,k)+xxx(i-1,j-1,k)+xxx(i-1,j,k)))/
     .           4./dppp(k)+               
     .           (xxxx(i,j,k)-xxxx(i,j,k-1))/2.*(xxx(i,j,k-1)+
     .           xxx(i,j-1,k-1)+xxx(i-1,j-1,k-1)+xxx(i-1,j,k-1)+
     .           (xxx(i,j,k)+xxx(i,j-1,k)+xxx(i-1,j-1,k)+xxx(i-1,j,k)))/
     .           4./dppp(k)
         twoa=r/f(i,j)/((p(k)+p(k-1))/2.)/100.*
     .         ((t(i-1,j,k)-t(i,j-1,k))/dd*
     .         (-vort(i-1,j-1,k-1)-vort(i-1,j-1,k)+
     .         vort(i,j,k)+vort(i,j,k-1))/dd-(t(i,j,k)-t(i-1,j-1,k))/dd*
     .         (vort(i-1,j,k-1)+vort(i-1,j,k)-vort(i,j-1,k)-
     .         vort(i,j-1,k-1))/dd)
         cee=(-v(i-1,j,k)+v(i-1,j,k-1)-v(i-1,j-1,k)+v(i,j,k-1))/
     .         dppp(k)/2.*beta
         xxxxx(i,j,k)=f(i,j)/stab*(twoa-twolam+cee)
         if (i .eq. 14 .and. j .eq. 11) 
     .      write(6,1999) xxxxx(i,j,k),f(i,j),twoa,twolam,cee
1999     format(1x,5e12.4)
      enddo
      enddo
      enddo
      call leib(om,xxxxx,40,1.e-5,nx,ny,nz,ps,p
     .         ,a,b,c,b,b,dx,dy,dppp,0)
      return
      end
c
c===============================================================================
c
      subroutine balcon(to,uo,vo,omo,t,u,v,om,del,tau,itmax,err,nu,nv
     .                 ,icon,nx,ny,nz,lat,dx,dy,ps,p,dp,gam)
c
c *** Balcon executes the mass/wind balance computations as described
c        mcginley (Meteor and Appl Phys, 1987).  A new
c        geopotential(t) is computed using relaxation on eqn. (2).  New
c        u, v and omega winds are computed using eqns. (4), (5) and (6)
c        with the new geopotential and neglecting the lagrange multiplier
c        term.  If icon is zero, balcon returns.  Otherwise, the
c        lagrange multiplier is computed using 3-d relaxation on eqn. (3).
c        U, v and omega are adjusted by adding the lagrange multiplier
c        term with the new lagrange multiplier.
c
      implicit none
c
      integer   nx,ny,nz
     .         ,nxm1,nym1,nzm1
     .         ,itmax,icon,ks,kf,ittr
     .         ,ibnd(nx,ny,nz)
     .         ,i,j,k,is,ip,js,jp,kp,it,itt
     .         ,icnt,iwpt,istatus
     .         ,ucnt,vcnt,uwpt,vwpt
c
      real*4 t(nx,ny,nz),to(nx,ny,nz)
     .      ,u(nx,ny,nz),uo(nx,ny,nz)
     .      ,v(nx,ny,nz),vo(nx,ny,nz)
     .      ,om(nx,ny,nz),omo(nx,ny,nz)
     .      ,nu(nx,ny,nz),nv(nx,ny,nz)
     .      ,dx(nx,ny),dy(nx,ny),dp(nz)
     .      ,ps(nx,ny),p(nz),gam(nz)
     .      ,lat(nx,ny),slam(nx,ny,nz)
     .      ,h(nx,ny,nz),f3(nx,ny,nz)
     .      ,bndx(nx,ny,nz),bndy(nx,ny,nz)
     .      ,del(nx,ny,nz),tau,err,rdpdg,bnd,g,fo,r,re,ovr
     .      ,dy2,dys,ang,f,cotmax,sin1,fs,cos1,beta
     .      ,a,bb,term1,term2,term3,term4,dx2,dxs,cortmt
     .      ,dudy,dvdx,dnudx,dnvdy,snv,tt,uot,vot,tot
     .      ,dt2dx2,dt2dy2,slap,force,rest,cot
     .      ,cotma1,cotm5,rho,cotm0,erf,dtdx,dtdy,nuu,nvv
     .      ,dldp,dldx,dldy,tsum,r_missing_data
     .      ,usum,vsum

      real*4 ttemp(nx,ny,nz)
      real*4 utemp(nx,ny,nz)
      real*4 vtemp(nx,ny,nz)

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
      do k=1,nz
      do j=1,ny
      do i=1,nx
         ibnd(i,j,k)=0
      enddo
      enddo
      enddo
      do k=1,nz
         kp=1
         if (k .eq. nz) kp=0
         do j=1,ny
            js=1
            jp=1
            if (j .eq. 1) js=0
            if (j .eq. ny) jp=0
            do i=1,nx
               is=1
               ip=1
               if (i .eq. 1) is=0
               if (i .eq. nx) ip=0
               to(i,j,k)=to(i,j,k)*g
               t(i,j,k)=t(i,j,k)*g
               if (ps(i,j) .le. p(k+kp)) then!turn on terrain switch
                  ibnd(i,j,k)=1
                  ibnd(i+ip,j+jp,k)=1
                  ibnd(i,j+jp,k)=1
                  ibnd(i+ip,j,k)=1
               endif
            enddo
         enddo
      enddo
c *** Compute new phi (t array) using relaxation on eqn. (2).
c        beta*dldx term is dropped to eliminate coupling with lambda eqn.
c
      do 1 it=1,itmax
         cotmax=0.
c        if(icon.eq.1) call prt(t,sc,ci,nx,ny,4,nz,4hgepo,4h    )
         do 2 j=2,nym1
         do 2 i=2,nxm1
            ang=(lat(i,j-1)+lat(i,j))/2.*rdpdg
            sin1=sin(ang)
            f=fo*sin1
            fs=f*f
            cos1=cos(ang)
            beta=fo*cos1/re
            dx2=dx(i,j)*2.
            dxs=dx(i,j)*dx(i,j)
            dy2=dy(i,j)*2.
            dys=dy(i,j)*dy(i,j)
            do 2 k=ks,kf
               a=(1.+fs*del(i,j,k))
               bb=(2.*f*beta*del(i,j,k))/a
               term1=beta
               term2=-f
               term3=-tan(ang)/re-bb
               term4=-gam(k)/del(i,j,k)*a
               cortmt=-2./dxs-2./dys+term4
               dudy=(uo(i,j,k)-uo(i,j-1,k))/dy(i,j)
               dvdx=(vo(i,j,k)-vo(i-1,j,k))/dx(i,j)
               dnudx=(nu(i+1,j,k)-nu(i-1,j,k)
     .               +nu(i+1,j-1,k)-nu(i-1,j-1,k))*.5/dx2
               dnvdy=(nv(i,j+1,k)-nv(i,j-1,k)
     .               +nv(i-1,j+1,k)-nv(i-1,j-1,k))*.5/dy2
               snv=(nv(i-1,j,k)+nv(i,j,k))*.5
               tt=t(i,j,k)
               dtdy=(t(i,j+1,k)-t(i,j-1,k))/dy2
               uot=(uo(i,j,k)+uo(i,j-1,k))*.5
               tot=to(i,j,k)
               dt2dx2=(t(i+1,j,k)+t(i-1,j,k)-2.*tt)/dxs
               dt2dy2=(t(i,j+1,k)+t(i,j-1,k)-tt*2.)/dys
               slap=dt2dx2+dt2dy2+dtdy*term3+term4*(tt)
               force=-(term4*(-tot)+uot*term1+
     .                term2*(dvdx-dudy)+dnudx+dnvdy-bb*snv)
               rest=slap-force
               cot=rest/cortmt
               t(i,j,k)=tt-cot*ovr
               cotmax=amax1(cotmax,abs(cot))
               if (it .eq. 1) cotma1=cotmax
2          continue
c        write(6,1000) it,itt,cotmax,ovr,cotma1
5        ittr=ittr+1
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
1000     format(1x,'it = ',i4,' max correction for ',a4,' = ',e12.3
     .         ,'ovr =  ',e12.4/1x
     .         ,'first iteration max correction was ',e12.4)
1     continue
12    write(6,1000) it,itt,cotmax,ovr,cotma1
c     erf=0.
c     if(icon.eq.1)call prt(t,sc,ci,nx,ny,4,nz,4hgeo,4hpot)
c69      call zero(slam,nx,ny,nz)
c
c *** Compute new u, v, omega using eqns. (4), (5), (6) with new phi and
c        without the lagrange multiplier terms.
c
      do k=1,nzm1
      do j=1,nym1
      do i=1,nxm1                                                  
         js=1
         if (j .eq. 1) js=0
         ang=(lat(i,j))*rdpdg
         sin1=sin(ang)
         f=fo*sin1
         fs=f*f
         a=1.+fs*del(i,j,k) 
         is=1
         if (i .eq. 1) is=0
         dtdx=(t(i+1,j,k)-t(i,j,k))/dx(i,j)
         dtdy=(t(i,j+1,k)-t(i,j,k))/dy(i,j)
         uot=uo(i,j,k)
         vot=vo(i,j,k)
         nvv=(nv(i,j,k)+nv(i,j+1,k)+nv(i-is,j+1,k)+nv(i-is,j,k))/4.
         nuu=(nu(i,j,k)+nu(i+1,j,k)+nu(i+1,j-js,k)+nu(i,j-js,k))/4.
c the solution wind will be from lw3 if the wind is near the bndry
c since u and v were from lw3 on entry to balcon
         if (uot .ne. bnd) then
            u(i,j,k)=(uot-f*del(i,j,k)*(nvv+dtdy))/a
          else 
            if(ps(i,j).lt.p(k)) u(i,j,k)=bnd
         endif
         if (vot .ne. bnd) then
            v(i,j,k)=(vot+f*del(i,j,k)*(nuu+dtdx))/a
          else
            if (ps(i,j).lt.p(k)) v(i,j,k)=bnd
         endif
      enddo
      enddo
      enddo
c
c *** Compute lagrange multiplier (slam) using 3-d relaxtion on eqn. (3).
c
      if (icon .ne. 0) then
         call zero(slam,nx,ny,nz)
         call zero(h,nx,ny,nz)
c
c ****** Compute a/tau (h) term and rhs terms in eqn. (3)
c
         call fthree(f3,u,v,omo,del,nu,nv,h,tau
     .              ,nx,ny,nz,lat,dx,dy,dp)
         if (icon .eq. 1) call analz(to,to,uo,uo,vo,vo,omo,omo
     .                   ,slam,f3,nu,nv,ibnd,del,tau
     .                   ,nx,ny,nz
     .                   ,lat,dx,dy,ps,p,dp,gam)
c
c ****** Perform 3-d relaxation.
c
         erf=.1
         call leibp3(slam,f3,200,erf,h
     .              ,nx,ny,nz,dx,dy,ps,p,dp)
c        call prt(f3,sc,ci,nx,ny,5,nz,4hforc ,4hng  )
c        call prt(slam,sc,ci,nx,ny,5,nz,4hlamd ,4ha   )
c
c ****** Compute new u, v, omega by adding the lagrange multiplier terms.
c
         do k=1,nzm1
         do j=1,nym1
         do i=1,nxm1
            ang=(lat(i,j))*rdpdg                                  
            sin1=sin(ang)
            f=fo*sin1
            fs=f*f
            dldp=(slam(i,j,k)-slam(i,j,k+1))/dp(k+1)*.01
            dldx=(slam(i+1,j+1,k+1)-slam(i,j+1,k+1))/dx(i,j)
            dldy=(slam(i+1,j+1,k+1)-slam(i+1,j,k+1))/dy(i,j)
            if (u(i,j,k) .ne. bnd) u(i,j,k)=u(i,j,k)+.5*dldx/a
            if (v(i,j,k) .ne. bnd) v(i,j,k)=v(i,j,k)+.5*dldy/a
            if (omo(i,j,k).ne.bnd) om(i,j,k)=omo(i,j,k)+.5*dldp/tau
         enddo
         enddo
         enddo
      endif
      if (icon .eq. 1) call analz(t,to,u,uo,v,vo,om,omo
     .                ,slam,f3,nu,nv,ibnd,del,tau
     .                ,nx,ny,nz
     .                ,lat,dx,dy,ps,p,dp,gam)
c     if (icon .eq. 1) call prt(om,sc,ci,nx,ny,4,nz,4h om ,4h    )   
c     if (icon .eq. 1) call prt(om,sc,ci,nx,ny,8,nz,4h om ,4h    )   
      do k=1,nz
      do j=1,ny
      do i=1,nx
         t(i,j,k)=t(i,j,k)/g
         ttemp(i,j,k)=t(i,j,k)
         utemp(i,j,k)=u(i,j,k)
         vtemp(i,j,k)=v(i,j,k)
      enddo
      enddo
      enddo
c
c insure that we don't have boundary pressure gradients due to
c boundary terrain/pressure intersections. Also check u/v components
c on the boundary.
c
      call get_r_missing_data(r_missing_data,istatus)

      do k=1,nz-1
c w/e sides
         do j=1,ny
            if(ps(1,j).le.p(k+1)) then 
               t(1,j,k)=r_missing_data
               u(1,j,k)=r_missing_data
               v(1,j,k)=r_missing_data
            endif
            if(ps(nx,j).le.p(k+1))then
               t(nx,j,k)=r_missing_data
               u(nx,j,k)=r_missing_data
               v(nx,j,k)=r_missing_data
            endif
         enddo
c s/n sides
         do i=1,nx
            if(ps(i,1).le.p(k+1)) then
               t(i,1,k)=r_missing_data
               u(i,1,k)=r_missing_data
               v(i,1,k)=r_missing_data
            endif
            if(ps(i,ny).le.p(k+1))then
               t(i,ny,k)=r_missing_data
               u(i,ny,k)=r_missing_data
               v(i,ny,k)=r_missing_data
            endif
         enddo
      enddo
c
      do k=1,nz-1
c western boundary
         iwpt=0
         do j=1,ny
            if(t(1,j,k).eq.r_missing_data)then
               tsum=0.0
               usum=0.0
               vsum=0.0
               icnt=0
               ucnt=0
               vcnt=0
               jp=1
               if(j.eq.ny.or.j.eq.1)jp=0
               if(ps(1,j+jp).gt.p(k))then
                  if(t(1,j+jp,k).ne.r_missing_data)then
                     tsum=t(1,j+jp,k)+tsum
                     icnt=icnt+1
                  endif
                  if(u(1,j+jp,k).ne.r_missing_data)then
                     usum=u(1,j+jp,k)+usum
                     ucnt=ucnt+1
                  endif
                  if(v(1,j+jp,k).ne.r_missing_data)then
                     vsum=v(1,j+jp,k)+vsum
                     vcnt=vcnt+1
                  endif
               endif
               if(ps(1,j-jp).gt.p(k))then
                  if(t(1,j-jp,k).ne.r_missing_data)then
                     tsum=t(1,j-jp,k)+tsum
                     icnt=icnt+1
                  endif
                  if(u(1,j-jp,k).ne.r_missing_data)then
                     usum=u(1,j-jp,k)+usum
                     ucnt=ucnt+1
                  endif
                  if(v(1,j-jp,k).ne.r_missing_data)then
                     vsum=v(1,j-jp,k)+vsum
                     vcnt=vcnt+1
                  endif
               endif
               if(ps(2,j).gt.p(k))then
                  if(t(2,j,k).ne.r_missing_data)then
                     tsum=t(2,j,k)+tsum
                     icnt=icnt+1
                  endif
                  if(u(2,j,k).ne.r_missing_data)then
                     usum=u(2,j,k)+usum
                     ucnt=ucnt+1
                  endif
                  if(v(2,j,k).ne.r_missing_data)then
                     vsum=v(2,j,k)+vsum
                     vcnt=vcnt+1
                  endif
               endif

               if(icnt.gt.0)then
                  t(1,j,k)=tsum/float(icnt)
               else
                  t(1,j,k)=ttemp(1,j,k)
                  iwpt=iwpt+1
               endif
               if(ucnt.gt.0)then
                  u(1,j,k)=usum/float(ucnt)
               else
                  u(1,j,k)=utemp(1,j,k)
                  uwpt=uwpt+1
               endif
               if(vcnt.gt.0)then
                  v(1,j,k)=vsum/float(vcnt)
               else
                  v(1,j,k)=vtemp(1,j,k)
                  vwpt=vwpt+1
               endif
            endif
         enddo
c        print*,'no hgt pts on western bndry for avg = ',iwpt
c        print*,'no u pts on western bndry for avg = ',uwpt
c        print*,'no v pts on western bndry for avg = ',vwpt
c eastern boundary
         iwpt=0
         uwpt=0
         vwpt=0
         do j=1,ny
            if(t(nx,j,k).eq.r_missing_data)then
               tsum=0.0
               usum=0.0
               vsum=0.0
               icnt=0
               ucnt=0
               vcnt=0
               jp=1
               if(j.eq.ny.or.j.eq.1)jp=0
               if(ps(nx,j+jp).gt.p(k))then
                  if(t(nx,j+jp,k).ne.r_missing_data)then
                     tsum=t(nx,j+jp,k)+tsum
                     icnt=icnt+1
                  endif
                  if(u(nx,j+jp,k).ne.r_missing_data)then
                     usum=u(nx,j+jp,k)+usum
                     ucnt=ucnt+1
                  endif
                  if(v(nx,j+jp,k).ne.r_missing_data)then
                     vsum=v(nx,j+jp,k)+vsum
                     vcnt=vcnt+1
                  endif
               endif
               if(ps(nx,j-jp).gt.p(k))then
                  if(t(nx,j-jp,k).ne.r_missing_data)then
                     tsum=t(nx,j-jp,k)+tsum
                     icnt=icnt+1
                  endif
                  if(u(nx,j-jp,k).ne.r_missing_data)then
                     usum=u(nx,j-jp,k)+usum
                     ucnt=ucnt+1
                  endif
                  if(v(nx,j-jp,k).ne.r_missing_data)then
                     vsum=v(nx,j-jp,k)+vsum
                     vcnt=vcnt+1
                  endif
               endif
               if(ps(nx-1,j).gt.p(k))then
                  if(t(nx-1,j,k).ne.r_missing_data)then
                     tsum=t(nx-1,j,k)+tsum
                     icnt=icnt+1
                  endif
                  if(u(nx-1,j,k).ne.r_missing_data)then
                     usum=u(nx-1,j,k)+usum
                     ucnt=ucnt+1
                  endif
                  if(v(nx-1,j,k).ne.r_missing_data)then
                     vsum=v(nx-1,j,k)+vsum
                     vcnt=vcnt+1
                  endif
               endif

               if(icnt.gt.0)then
                  t(nx,j,k)=tsum/float(icnt)
               else
                  t(nx,j,k)=ttemp(nx,j,k)
                  iwpt=iwpt+1
               endif
               if(ucnt.gt.0)then
                  u(nx,j,k)=usum/float(ucnt)
               else
                  u(nx,j,k)=utemp(nx,j,k)
                  uwpt=uwpt+1
               endif
               if(vcnt.gt.0)then
                  v(nx,j,k)=vsum/float(vcnt)
               else
                  v(nx,j,k)=vtemp(nx,j,k)
                  vwpt=vwpt+1
               endif
            endif
         enddo
c        print*,'no pts on eastern bndry for avg = ',iwpt
c        print*,'no u pts on eastern bndry for avg = ',uwpt
c        print*,'no v pts on eastern bndry for avg = ',vwpt
c southern boundary
         iwpt=0
         do i=1,nx
            if(t(i,1,k).eq.r_missing_data)then
               tsum=0.0
               usum=0.0
               vsum=0.0
               icnt=0
               ucnt=0
               vcnt=0
               ip=1
               if(i.eq.nx.or.i.eq.1)ip=0
               if(ps(i+ip,1).gt.p(k))then
                  if(t(i+ip,1,k).ne.r_missing_data)then
                     tsum=t(i+ip,1,k)+tsum
                     icnt=icnt+1
                  endif
                  if(u(i+ip,1,k).ne.r_missing_data)then
                     usum=u(i+ip,1,k)+usum
                     ucnt=ucnt+1
                  endif
                  if(v(i+ip,1,k).ne.r_missing_data)then
                     vsum=v(i+ip,1,k)+vsum
                     vcnt=vcnt+1
                  endif
               endif
               if(ps(i-ip,1).gt.p(k))then
                  if(t(i-ip,1,k).ne.r_missing_data)then
                     tsum=t(i-ip,1,k)+tsum
                     icnt=icnt+1
                  endif
                  if(u(i-ip,1,k).ne.r_missing_data)then
                     usum=u(i-ip,1,k)+usum
                     ucnt=ucnt+1
                  endif
                  if(v(i-ip,1,k).ne.r_missing_data)then
                     vsum=v(i-ip,1,k)+vsum
                     vcnt=vcnt+1
                  endif
               endif
               if(ps(i,2).gt.p(k))then
                  if(t(i,2,k).ne.r_missing_data)then
                     tsum=t(i,2,k)+tsum
                     icnt=icnt+1
                  endif
                  if(u(i,2,k).ne.r_missing_data)then
                     usum=u(i,2,k)+usum
                     ucnt=ucnt+1
                  endif
                  if(v(i,2,k).ne.r_missing_data)then
                     vsum=v(i,2,k)+vsum
                     vcnt=vcnt+1
                  endif
               endif

               if(icnt.gt.0)then
                  t(i,1,k)=tsum/float(icnt)
               else
                  t(i,1,k)=ttemp(i,1,k)
                  iwpt=iwpt+1
               endif
               if(ucnt.gt.0)then
                  u(i,1,k)=usum/float(ucnt)
               else
                  u(i,1,k)=utemp(i,1,k)
                  uwpt=uwpt+1
               endif
               if(vcnt.gt.0)then
                  v(i,1,k)=vsum/float(vcnt)
               else
                  v(i,1,k)=vtemp(i,1,k)
                  vwpt=vwpt+1
               endif
            endif
         enddo
c        print*,'no hgt pts on southern bndry for avg = ',iwpt
c        print*,'no u pts on southern bndry for avg = ',uwpt
c        print*,'no v pts on southern bndry for avg = ',vwpt
c northern boundary
         iwpt=0
         uwpt=0
         vwpt=0
         do i=1,nx
            if(t(i,ny,k).eq.r_missing_data)then
               tsum=0.0
               usum=0.0
               vsum=0.0
               icnt=0
               ucnt=0
               vcnt=0
               ip=1
               if(i.eq.nx.or.i.eq.1)ip=0
               if(ps(i+ip,ny).ge.p(k))then
                  if(t(i+ip,ny,k).ne.r_missing_data)then
                     tsum=t(i+ip,ny,k)+tsum
                     icnt=icnt+1
                  endif
                  if(u(i+ip,ny,k).ne.r_missing_data)then
                     usum=u(i+ip,ny,k)+usum
                     ucnt=ucnt+1
                  endif
                  if(v(i+ip,ny,k).ne.r_missing_data)then
                     vsum=v(i+ip,ny,k)+vsum
                     vcnt=vcnt+1
                  endif
               endif
               if(ps(i-ip,ny).ge.p(k))then
                  if(t(i-ip,1,k).ne.r_missing_data)then
                     tsum=t(i-ip,ny,k)+tsum
                     icnt=icnt+1
                  endif
                  if(u(i-ip,1,k).ne.r_missing_data)then
                     usum=t(i-ip,ny,k)+usum
                     ucnt=ucnt+1
                  endif
                  if(v(i-ip,1,k).ne.r_missing_data)then
                     vsum=t(i-ip,ny,k)+vsum
                     vcnt=vcnt+1
                  endif
               endif
               if(ps(i,ny-1).ge.p(k))then
                  if(t(i,ny-1,k).ne.r_missing_data)then
                     tsum=t(i,ny-1,k)+tsum
                     icnt=icnt+1
                  endif
                  if(u(i,ny-1,k).ne.r_missing_data)then
                     usum=u(i,ny-1,k)+usum
                     ucnt=ucnt+1
                  endif
                  if(v(i,ny-1,k).ne.r_missing_data)then
                     vsum=v(i,ny-1,k)+vsum
                     vcnt=vcnt+1
                  endif
               endif
               if(icnt.gt.0)then
                  t(i,ny,k)=tsum/float(icnt)
               else
                  t(i,ny,k)=ttemp(i,ny,k)
                  iwpt=iwpt+1
               endif
               if(ucnt.gt.0)then
                  u(i,ny,k)=usum/float(ucnt)
               else
                  u(i,ny,k)=utemp(i,ny,k)
                  uwpt=uwpt+1
               endif
               if(vcnt.gt.0)then
                  v(i,ny,k)=vsum/float(vcnt)
               else
                  v(i,ny,k)=vtemp(i,ny,k)
                  vwpt=vwpt+1
               endif
            endif
         enddo
c        print*,'no hgt pts on northern bndry for avg = ',iwpt
c        print*,'no u pts on northern bndry for avg = ',uwpt
c        print*,'no v pts on northern bndry for avg = ',vwpt
c        print*,'--------------------------------------'
c        print*
      enddo
c
      return
      end
c
c===============================================================================
c
      subroutine analz(t,to,u,uo,v,vo,om,omo,slam,f3,nu,nv,ibnd,del,tau
     .                ,nx,ny,nz,lat,dx,dy,ps,p,dp,gam)
c
      implicit none
c
      integer   nx,ny,nz
     .         ,ibnd(nx,ny,nz)
     .         ,nzm1,nxm2,nym2
     .         ,isv,jsv,ksv
     .         ,i,j,k
c
      real*4 t(nx,ny,nz),to(nx,ny,nz)
     .      ,u(nx,ny,nz),uo(nx,ny,nz)
     .      ,v(nx,ny,nz),vo(nx,ny,nz)
     .      ,om(nx,ny,nz),omo(nx,ny,nz)
     .      ,slam(nx,ny,nz),f3(nx,ny,nz)
     .      ,nu(nx,ny,nz),nv(nx,ny,nz)
     .      ,lat(nx,ny),dx(nx,ny),dy(nx,ny)
     .      ,ps(nx,ny),p(nz),dp(nz),gam(nz)
     .      ,del(nx,ny,nz),tau,rdpdg,fo,bnd
     .      ,conmax,sumom,cont,con,sumu,sumv,sumt,thermu,thermv
     .      ,sumww,resu,resv,tgpu,tgpv,tgpph,tgpc
     .      ,ang,sin1,f,nvv,nuu,dtdx,dtdy,dudx,dvdy,domdp
     .      ,uot,vot,tot,uuu,vvv,contm
c_______________________________________________________________________________
c
      nzm1=nz-1
      nxm2=nx-2
      nym2=ny-2
      rdpdg=3.141592654/180.
      fo=14.52e-5   !2*omega
      bnd=1.e-30
      do k=1,nzm1
         write(6,4000) k,p(k+1)
4000     format(1x,'Level ',i4,'   ',f5.0,' mb')
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
         do j=2,nym2
         do i=2,nxm2
            ang=(lat(i,j))*rdpdg
            sin1=sin(ang)
            f=fo*sin1
            nvv=(nv(i,j,k)+nv(i,j+1,k)+nv(i-1,j+1,k)+nv(i-1,j,k))/4.
            nuu=(nu(i,j,k)+nu(i+1,j,k)+nu(i+1,j-1,k)+nu(i,j-1,k))/4.
            dtdx=(t(i+1,j,k)-t(i,j,k))/dx(i,j)
            dtdy=(t(i,j+1,k)-t(i,j,k))/dy(i,j)
            if (k .ne. 1) then
               dudx=(u(i,j-1,k-1)-u(i-1,j-1,k-1))/dx(i,j)
               dvdy=(v(i-1,j,k-1)-v(i-1,j-1,k-1))/dy(i,j)
               domdp=(om(i,j,k-1)-om(i,j,k))/dp(k)*.01
            endif
            uot=uo(i,j,k)
            vot=vo(i,j,k)
            tot=to(i,j,k)
            if (uot .ne. bnd) then
               thermu=(u(i,j,k)+dtdy/f)**2+thermu
               sumu=sumu+(u(i,j,k)-uot)**2
               resv=(f*u(i,j,k)+dtdy+nvv)**2+resv
               tgpu=tgpu+1.
            endif
            if (vot .ne. bnd) then
               thermv=(v(i,j,k)-dtdx/f)**2+thermv
               sumv=sumv+(v(i,j,k)-vot)**2
               resu=resu+(-f*v(i,j,k)+dtdx+nuu)**2
               tgpv=tgpv+1.
            endif
            if (ibnd(i,j,k) .ne. 1) then
               sumt=sumt+(t(i,j,k)-tot)**2
               tgpph=tgpph+1.
            endif
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
               endif
            endif
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
1002     format(1x,' rms wind speed:',e12.4)
         if (sumww .ne. 0.) sumww=sqrt(thermu**2+thermv**2)/sumww
c
c commented out write statment (J.Smart). This has to do with parameter del which is
c                                         now a 3-d array.

c        write(6,1001) sumt,del,sumu,gam(k),sumv,
c    .             tau,sumom,cont,thermu,thermv,sumww,resu,resv
c        if (k .ne. 1) write(6,1009)
c    .      isv,jsv,ksv,contm,f3(isv,jsv,ksv),slam(isv,jsv,ksv),
c    .      slam(isv+1,jsv,ksv),slam(isv-1,jsv,ksv),slam(isv,jsv+1,ksv),
c    .      slam(isv,jsv-1,ksv),slam(isv,jsv,ksv+1),slam(isv,jsv,ksv-1),
c    .      dx(isv,jsv),dy(isv,jsv),dp(ksv),
c    .      u(isv,jsv-1,ksv-1),u(isv-1,jsv-1,ksv-1),v(isv-1,jsv,ksv-1),
c    .      v(isv-1,jsv-1,ksv-1),om(isv,jsv,ksv),om(isv,jsv,ksv-1)
      enddo
1009  format(1x,'relax ',3i4/1x,'con f3 l ',3e12.4/1x,'lijk ',6e12.4/1x,
     .       ' d x y p ',3e12.4/1x,'uuvvoo ',6e12.4)
1001  format(1x,' rms error for each term in functional and del,gam,tau'
     . /1x,'    phi-phio:',e12.5,'   del:',e12.5,
     . /1x,'    u-uo:    ',e12.5,'   gam:',e12.5,
     . /1x,'    v-vo:    ',e12.5,'   tau:',e12.5,
     . /1x,'    om-omo:  ',e12.5
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
      subroutine terbnd(u,v,om,nx,ny,nz,ps,p)
c  terrain boundary sets velocity component to 1.e-30 on terrain 
c  faces and in the interior of terrain.  this rountine has
c been modified for the AFWA implementation.  It does not use 
c the original stagger associated with the mcginley(1987) paper
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
      bnd=1.e-30
      nym1=ny-1
      nxm1=nx-1
      nzm1=nz-1
      do j=2,nym1
      do i=2,nxm1
         om(i,j,1)=bnd
         do k=2,nzm1
            if (ps(i,j) .le. p(k)) then
c                u(i-1,j-1,k-1)=bnd
c                u(i,j-1,k-1)=bnd
c                v(i-1,j-1,k-1)=bnd
c                v(i-1,j,k-1)=bnd
c  These are AFWA changes
                 u(i,j,k)=bnd
                 u(i+1,j,k)=bnd
                 v(i,j,k)=bnd
                 v(i,j+1,k)=bnd
                 om(i,j,k-1)=bnd
                 om(i,j,k)=bnd
            endif
         enddo
      enddo
      enddo
c
      do k=2,nz
         do j=1,ny
            if(ps(1,j) .le. p(k)) then
c  These are AFWA changes
               u(1,j,k)=bnd
               v(1,j,k)=bnd
               om(1,j,k)=bnd
               om(1,j,k-1)=bnd
            endif
            if(ps(nx,j) .le. p(k)) then
               u(nx,j,k)=bnd
               v(nx,j,k)=bnd
               om(nx,j,k)=bnd
               om(nx,j,k-1)=bnd
            endif
         enddo

         do i=1,nx
            if(ps(i,1) .le. p(k)) then
c  These are AFWA changes
               u(i,1,k)=bnd
               v(i,1,k)=bnd
               om(i,1,k)=bnd
               om(i,1,k-1)=bnd
            endif
            if(ps(i,ny) .le. p(k)) then
               u(i,ny,k)=bnd
               v(i,ny,k)=bnd
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
      subroutine nonlin(nu,nv,uo,vo,u,v,om
     .                 ,nx,ny,nz,dx,dy,dp,dt)
c
c *** Nonlin computes the non-linear terms (nu,nv).
c *** For now, nonlin assumes a non-staggered grid and uses a backwards
c        time difference.
c
      implicit none
c
      integer   nx,ny,nz
     .         ,nxm1,nym1,nzm1
     .         ,i,j,k
c
      real*4 nu(nx,ny,nz),nv(nx,ny,nz)
     .      ,uo(nx,ny,nz),vo(nx,ny,nz)   !time=t
     .      ,u(nx,ny,nz),v(nx,ny,nz)     !time=t-dt
     .      ,om(nx,ny,nz)
     .      ,dx(nx,ny),dy(nx,ny),dp(nz)
     .      ,dt,dudx,dvdy,dudy,dvdx
     .      ,dudp,dvdp,dudt,dvdt
c_______________________________________________________________________________
c
c *** A 1 2 1 time filter is applied over all gradient terms (not used now...).
c
      print *,'nonlin'
      nxm1=nx-1
      nym1=ny-1
      nzm1=nz-1
      call zero(nu,nx,ny,nz)
      call zero(nv,nx,ny,nz)
      do k=2,nzm1
      do j=2,nym1
      do i=2,nxm1
         dudx=(uo(i+1,j,k)-uo(i-1,j,k))/(2.*dx(i,j))
         dvdx=(vo(i+1,j,k)-vo(i-1,j,k))/(2.*dx(i,j))
         dudy=(uo(i,j+1,k)-uo(i,j-1,k))/(2.*dy(i,j))
         dvdy=(vo(i,j+1,k)-vo(i,j-1,k))/(2.*dy(i,j))
         dudp=(uo(i,j,k+1)-uo(i,j,k-1))/(dp(k)+dp(k+1))
         dvdp=(vo(i,j,k+1)-vo(i,j,k-1))/(dp(k)+dp(k+1))
         dudt=(uo(i,j,k)-u(i,j,k))/dt
         dvdt=(vo(i,j,k)-v(i,j,k))/dt
         nu(i,j,k)=dudt+uo(i,j,k)*dudx+vo(i,j,k)*dudy+om(i,j,k)*dudp
         nv(i,j,k)=dvdt+uo(i,j,k)*dvdx+vo(i,j,k)*dvdy+om(i,j,k)*dvdp
         if (i .eq. 7 .and. j .eq. 12 .and. k .eq. 12) write(6,1000)
     .       dudt,dudx,dudy,dudp,dvdt,dvdx,dvdy,dvdp
1000      format(1x,'nonlin ',7e12.4/1x,7e12.4)
      enddo
      enddo
      enddo
c
      return
      end
c
c===============================================================================
c
      subroutine fthree(f3,u,v,om,del,nu,nv,h,tau
     .                 ,nx,ny,nz,lat,dx,dy,dp)
c
c *** Fthree computes a/tau (h) and rhs terms in eqn. (3).
c
      implicit none
c
      integer   nx,ny,nz
     .         ,i,j,k,is,js,ks
c
      real*4 f3(nx,ny,nz),h(nx,ny,nz)
     .      ,u(nx,ny,nz),v(nx,ny,nz)
     .      ,om(nx,ny,nz),lat(nx,ny)
     .      ,nu(nx,ny,nz),nv(nx,ny,nz)
     .      ,dx(nx,ny),dy(nx,ny),dp(nz)
     .      ,del(nx,ny,nz),tau,fo,rdpdg,dpp,f,aa,dxx,dyy
     .      ,dnudy,dnvdx,cont,formax
c_______________________________________________________________________________
c
      print *,'fthree'
      fo=14.52e-5   !2*omega
      rdpdg=3.141592654/180.
      formax=0.
      do k=2,nz
         ks=1
         if (k .eq. 1) ks=0
         dpp=dp(k)*100.
         do j=2,ny
         do i=2,nx
            f=sin(lat(i,j)*rdpdg)*fo
            aa=1.+f*f*del(i,j,k) 
            js=1
            if (j .eq. 2) js=0
            dyy=dy(i,j)*float(js+1)
            h(i,j,k)=aa/tau
            is=1
            if (i .eq. 2) is=0
            dxx=dx(i,j)*float(is+1)
            dnudy=(nu(i,j,k-1)+nu(i-1,j,k-1)
     .            -nu(i,j-js-1,k-1)-nu(i-1,j-js-1,k-1))/dy(i,j)*.25
            dnvdx=(nv(i,j,k-1)+nv(i,j-1,k-1)
     .            -nv(i-is-1,j,k-1)-nv(i-is-1,j-1,k-1))/dx(i,j)*.25
            cont=((u(i,j-1,k-1)-u(i-1,j-1,k-1))/dx(i,j)
     .           +(v(i-1,j,k-1)-v(i-1,j-1,k-1))/dy(i,j)
     .           +(om(i,j,k-1)-om(i,j,k))/dpp)*aa
            f3(i,j,k)=-2.*cont
            if (abs(f3(i,j,k)) .ge. formax) then
c              write(6,3004) i,j,k,cont,dnvdx,dnudy,
c    .                       dxx,dyy,nu(i,j,k),nv(i,j,k)  
3004           format(1x,3i4,8e12.4)
c              write(6,3004) i,j,k,nu(i,j,k-ks)
c    .                      ,nu(i-is,j,k-ks),nu(i,j-js-1,k-ks)
c    .                      ,nu(i-is,j-js-1,k-ks),nv(i,j,k-ks)
c    .                      ,nv(i,j-js,k-ks),nv(i-is-1,j,k-ks)
c    .                      ,nv(i-is-1,j-js,k-ks)
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
     .         ,ittr
     .         ,i,j,k,it,itmax,ia
c
      real*4 sol(nx,ny,nz),force(nx,ny,nz)
     .      ,h(nx,ny,nz)
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
c        sxx+syy+a*szz+b*sxz+c*syz+d*sxy+e*sx+f*sy+g*sz+h*s-force=0   
c
      print *,'leibp3'
      ovr=1.
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
         do j=2,nym1
         do i=2,nxm1
            dx2=dx(i,j)*2.
            dx1s=dx(i,j)*dx(i,j)
            dy2=dy(i,j)*2.
            dy1s=dy(i,j)*dy(i,j)
            do 2 k=2,nzm1
            dz=dp(k)*100.
            dz2=dz*2.
            dz1s=dz*dz
            if (ps(i,j)-p(k)) 2,2,20
20          aa=h(i,j,k)
            cortm=-2.*si/dx1s-2.*sj/dy1s-aa*2.*sk/dz1s+hh
                  res=(sol(i+1,j,k)+sol(i-1,j,k))/dx1s*si+
     .           (sol(i,j+1,k)+sol(i,j-1,k))/dy1s*sj+
     .           (sol(i,j,k+1)+sol(i,j,k-1))/dz1s*sk*
     .           aa+cortm*sol(i,j,k)-force(i,j,k)
            if (ps(i+1,j)-p(k)) 100,101,101
100             res=res+(-sol(i+1,j,k)+sol(i,j,k))/dx1s
            cortm=cortm+1./dx1s
101             if (ps(i,j+1)-p(k)) 102,103,103
102             res=res+(-sol(i,j+1,k)+sol(i,j,k))/dy1s
            cortm=cortm+1./dy1s
103             if (ps(i-1,j)-p(k)) 104,105,105
104             res=res+(-sol(i-1,j,k)+sol(i,j,k))/dx1s
            cortm=cortm+1./dx1s
105             if (ps(i,j-1)-p(k)) 106,107,107
106             res=res+(-sol(i,j-1,k)+sol(i,j,k))/dy1s
            cortm=cortm+1./dy1s
107             if (ps(i,j)-p(k-1)) 108,109,109
108             res=res-(sol(i,j,k-1)-sol(i,j,k))*aa/dz1s
            cortm=cortm+aa/dz1s
109             continue
            cor=res/cortm
            corb=5.*erf+1.
            if (it .ne. 1) corb=cor/corlmm
            if (abs(corb) .gt. erf) ia=1
            if (abs(cor) .gt. corlm) corlm=abs(cor)
c            if(i.eq.14.and.k.eq.4) write(6,3004) 
c      1         j,sol(i,j,k),force(i,j,k),res,cor,aa,cortm,
c      2         sol(i+1,j,k),sol(i-1,j,k),sol(i,j+1,k),sol(i,j-1,k),
c      3         sol(i,j,k+1),sol(i,j,k-1)
3004             format(1x,'relax ',i3,6e12.4/1x,6e12.4)
            sol(i,j,k)=sol(i,j,k)-cor*ovr
2            continue
         enddo
         enddo
c        write(6,1200) ((sol(i,11,22-k),i=6,8),k=1,21)
         reslm=corlm*cortm
c        write(6,1001) it,reslm,corlm,corlmm,erb
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
      enddo
      reslm=corlm*cortm
      write(6,1001) it,reslm,corlm,corlmm,erb
      write(6,1002) ovr
1002       format(1x,'ovr rlxtn const at fnl ittr = ',e10.4)
1001       format(1x,'iterations= ',i4,' max residual= ',e10.3,
     .    ' max correction= ',e10.3, ' first iter max cor= ',e10.3,
     .    'max bndry error= ',e10.3)
c
      return
      end
c
c===============================================================================
c
      subroutine zero(a,nx,ny,nz)
c
      implicit none
c
      integer   nx,ny,nz,i,j,k
c
      real*4 a(nx,ny,nz)
c_______________________________________________________________________________
c
      do k=1,nz
      do j=1,ny
      do i=1,nx
         a(i,j,k)=0.
      enddo
      enddo
      enddo
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
      Subroutine balstagger(u,v,om,phi,t,nx,ny,nz,
     &                      wr1,wr2,wr3,p,idstag)
c
c  This routine takes a standard LAPS field with all variables at
c  each grid point and produces the E-stagger appropriate for applying
c  the dynamic balancing in qbalpe.f  idstag > 0 staggers (LAPS -> stagger)
c  idstag < 0 destaggers (stagger -> LAPS)

      real u(nx,ny,nz),v(nx,ny,nz),om(nx,ny,nz),t(nx,ny,nz), 
     &phi(nx,ny,nz),wr1(nx,ny,nz),wr2(nx,ny,nz),wr3(nx,ny,nz)
     &,p(nz)

c  gas const
      r=287.04
c  gravity
      g=9.80665
      if(idstag.eq.0) then
         print*, 'no staggering accomplished'
         return
      endif
      if(idstag.gt.0) then

c set vertical stagger first
c first wind level for balcon is second level in LAPS
c omega is shifted one-half in vertical
c t is shifted likewise
c phi is shifted upward one level like winds
c level nz will be the same as laps for all fields

       do k=1,nz-1
        do j=1,ny
         do i=1,nx
          u(i,j,k)=u(i,j,k+1)
          v(i,j,k)=v(i,j,k+1)
          t(i,j,k)=(t(i,j,k)+t(i,j,k+1))*.5
          om(i,j,k)=(om(i,j,k)+om(i,j,k+1))*.5
         enddo
        enddo
       enddo

c horzizontal stagger

       do k=1,nz
        do j=1,ny-1
         do i=1,nx-1
          u(i,j,k)=(u(i+1,j+1,k)+u(i,j+1,k))*.5
          v(i,j,k)=(v(i+1,j+1,k)+v(i+1,j,k))*.5
          t(i,j,k)=(t(i,j,k)+t(i+1,j+1,k)+t(i,j+1,k)+t(i+1,j,k))*.25
c omega is already on the standard horizontal mesh
         enddo
        enddo
       enddo
c re integrate phi hydrostatically from staggered t using phi 1 from LAPS
       do k=1,nz-1
        do j=1,ny
         do i=1,nx
          phi(i,j,k)=phi(i,j,k)+r*t(i,j,k)*alog(p(k)/p(k+1))/g
         enddo
        enddo
       enddo
c for level nz 
       do j=1,ny
        do i=1,nx
          phi(i,j,nz)=phi(i,j,nz)+r*t(i,j,nz)*alog(p(nz-1)/p(nz))/g
        enddo
       enddo
       
       return

      else! de-stagger

c we begin by writing level 1 staggered into the laps level 1; use work arrays
c winds first
       do j=1,ny
        do i=1,nx
         wr1(i,j,1)=u(i,j,1)
         wr2(i,j,1)=v(i,j,1)
         wr3(i,j,1)=om(i,j,1)
        enddo
       enddo
c horizontal destagger

       do k=2,nz
        do j=2,ny
         do i=2,nx
          wr1(i,j,k)=(u(i,j-1,k-1)+u(i-1,j-1,k-1))*.5
          wr2(i,j,k)=(v(i-1,j,k-1)+v(i-1,j-1,k-1))*.5
          wr3(i,j,k)=(om(i,j,k)+om(i,j,k-1))*.5
         enddo
        enddo
       enddo

c now some of the boundary rows and columns

       do k=2,nz
        wr1(1,1,k)=u(1,1,k-1)
        wr2(1,1,k)=v(1,1,k-1)
        wr3(1,1,k)=(om(1,1,k-1)+om(1,1,k))*.5
        do i=2,nx
         wr1(i,1,k)=wr1(i,2,k)
         wr2(i,1,k)=v(i-1,1,k-1)
         wr3(i,1,k)=(om(i,1,k)+om(i,1,k-1))*.5
        enddo
        do j=2,ny
         wr1(1,j,k)=u(1,j-1,k-1)
         wr2(1,j,k)=wr2(2,j,k)
         wr3(1,j,k)=(om(1,j,k)+om(1,j,k-1))*.5
        enddo
       enddo ! on k

c transfer over to u,v,om arrays

       do k=1,nz
        do j=1,ny
         do i=1,nx
          u(i,j,k)=wr1(i,j,k)
          v(i,j,k)=wr2(i,j,k)
          om(i,j,k)=wr3(i,j,k)
         enddo
        enddo
       enddo
 
c now geopotential

       do k=2,nz
        do j=2,ny
         do i=2,nx
          wr3(i,j,k)=(phi(i-1,j-1,k-1)+phi(i,j-1,k-1)
     &               +phi(i-1,j,k-1)+phi(i,j,k-1))*.25
          wr2(i,j,k)=(t(i-1,j-1,k-1)+t(i,j-1,k-1)
     &               +t(i-1,j,k-1)+t(i,j,k))*.25
         enddo
        enddo
       enddo
       do k=2,nz
        wr3(1,1,k)=phi(1,1,k-1)
        wr2(1,1,k)=t(1,1,k-1)
       enddo
       do k=2,nz
        do j=2,ny
         wr3(1,j,k)=(phi(1,j,k-1)+phi(1,j-1,k-1))*.5
         wr2(1,j,k)=(t(i,1,k-1)+t(i,1,k-1))*.5
        enddo
        do i=2,nx
         wr3(i,1,k)=(phi(i,1,k-1)+phi(i,1,k-1))*.5
         wr2(i,1,k)=(t(i,1,k-1)+t(i,1,k-1))*.5
        enddo
       enddo
       do j=1,ny
        do i=1,nx
         wr3(i,j,1)=wr3(i,j,2)+r*wr2(i,j,2)*alog(p(2)/p(1))/g
         wr2(i,j,1)=wr2(i,j,2)
        enddo
       enddo
c call phig to compute temps from heights on the native laps grid
c temps will sill be staggered and must be destaggered below
       call phig(wr3,wr2,nx,ny,nz,0,p,-1)
       do k=2,nz
        do j=1,ny
         do i=1,nx
          wr2(i,j,k)=(wr2(i,j,k)+wr2(i,j,k-1))*.5
         enddo
        enddo
       enddo

       do k=1,nz
        do j=1,ny
         do i=1,nx
          phi(i,j,k)=wr3(i,j,k)
          t(i,j,k)=wr2(i,j,k)
         enddo
        enddo
       enddo

      endif

      return
      end

c
      subroutine phig(phi,t,nx,ny,nz,itshif,p,ittop)
c
c *** phig computes heights from temps or vice versa for
c        ittop equal to 1 and -1 respectively.
c
      implicit none
c
      integer*4 nx,ny,nz
     .         ,itshif,ittop
     .         ,i,j,k
c
      real*4 phi(nx,ny,nz),t(nx,ny,nz)
     .      ,p(nz)
     .      ,z(50),r,g,rog,gor,ddz,ddzi
c
      parameter (r=287.053,g=9.80665,rog=r/g,gor=g/r)
c_______________________________________________________________________________
c
      do k=1,nz
         z(k)=alog(p(1)/p(k))
      enddo
c
c *** Only for data that has no temps at lvl 1 (1025mb).
c 
      if (itshif .ne. 0) then
         do j=1,ny
         do i=1,nx
            t(i,j,1)=t(i,j,2)*2.-t(i,j,3)
         enddo
         enddo
      endif
c
c *** Scheme assumes 1000mb ht is correct.
c
      if (ittop .eq. 1) then
         do k=2,nz
            ddz=z(k)-z(k-1)
            do j=1,ny
            do i=1,nx
               phi(i,j,k)=(t(i,j,k))*rog*ddz+phi(i,j,k-1)
            enddo
            enddo
         enddo
      elseif (ittop .eq. -1) then
         do k=2,nz
            ddzi=1./(z(k)-z(k-1))
            do j=1,ny
            do i=1,nx
               t(i,j,k)=gor*(phi(i,j,k)-phi(i,j,k-1))*ddzi
            enddo
            enddo
         enddo
      endif
c
      return
      end
c
      subroutine phigns(phi,t,nx,ny,nz,itshif,p,ittop)
cThis is the version of phig for the AFWA implementation. 
c
c *** phig computes heights from temps or vice versa for
c        ittop equal to 1 and -1 respectively.
c
      implicit none
c
      integer*4 nx,ny,nz
     .         ,itshif,ittop
     .         ,i,j,k
c
      real*4 phi(nx,ny,nz),t(nx,ny,nz)
     .      ,p(nz)
     .      ,z(50),r,g,rog,gor,ddz,ddzi,ddzj
c
      parameter (r=287.053,g=9.80665,rog=r/g,gor=g/r)
c_______________________________________________________________________________
c
      do k=1,nz
         z(k)=alog(p(1)/p(k))
      enddo
c
c *** Only for data that has no temps at lvl 1 (1025mb).
c 
      if (itshif .ne. 0) then
         do j=1,ny
         do i=1,nx
            t(i,j,1)=t(i,j,2)*2.-t(i,j,3)
         enddo
         enddo
      endif
c
c *** Scheme assumes 1100mb ht is correct.
c
      if (ittop .eq. 1) then
         do k=2,nz
            ddz=z(k)-z(k-1)
            do j=1,ny
            do i=1,nx
               phi(i,j,k)=.5*(t(i,j,k)+t(i,j,k-1))*rog*ddz+phi(i,j,k-1)
            enddo
            enddo
         enddo
      elseif (ittop .eq. -1) then
c scheme returns input temps at k=1 and nz
         do k=2,nz-1 
            ddzi=1./(z(k)-z(k-1))
            ddzj=1./(z(k+1)-z(k))
            do j=1,ny
            do i=1,nx
               t(i,j,k)=.5*gor*(phi(i,j,k)-phi(i,j,k-1))*ddzi+
     &                   .5*gor*(phi(i,j,k+1)-phi(i,j,k))*ddzj
            enddo
            enddo
         enddo
      endif
c
      return
      end
c
     
