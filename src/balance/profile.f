      subroutine profile(udrop,vdrop,tdrop,rri,rrj,oberu,oberw,obert
     1,moderu,moderv,modert,erru, errv, errw,errt, u,v,om,sh,us,vs,oms,
     1 t,phi,ub,vb,omb,tb,p,ter,zter, nx,ny,nz,i4time)
c
c  this routine takes the measured fit between the obs and background
c  and constructs analysis error estimates based on an estimated
c  observation error, a projected truth based on gaussian random error.
c  arrays erru errv errw errt are the output error estimates  
c  variables are stored in each. Here is the key for each k profile:
c
c  erru(1,1,k) background u bias; (2,1,k) backgound variance; (1,5,k) u profile
c  (1,6,k) turb u comp; (1,2,k) anal bias; (2,2,k) anal var;
c
c  errv(1,1,k) background v bias; (2,1,k) backgound variance; (1,5,k) v profile
c  (1,6,k) turb v comp; (1,2,k) anal bias; (2,2,k) anal var;
c
c  errw(1,1,k) background w bias; (2,1,k) backgound variance; (1,5,k) w profile
c  (1,6,k) turb w comp; (1,2,k) anal bias; (2,2,k) anal var;
c
c  errt(1,1,k) background den bias; (2,1,k) backgound variance; (1,5,k)  phi profile
c  (1,6,k) T profile; (1,7,k) density profile               
c  (1,2,k) anal bias; (2,2,k) anal var; 
c  (1,3,k) press anal bias; (2,3,k) press anal var
c  (1,9,k) dropsonde profile;(1,10,k) background profiles
c  errv(1,1,k) 
c  us,vs,oms are the turbulent components here to evaluate covariance

      real erru(nx,ny,nz),errv(nx,ny,nz),errw(nx,ny,nz),errt(nx,ny,nz)
     1   ,u(nx,ny,nz), ub(nx,ny,nz), us(nx,ny,nz)
     1   ,v(nx,ny,nz), vb(nx,ny,nz),vs(nx,ny,nz)
     1   ,t(nx,ny,nz), tb(nx,ny,nz),oms(nx,ny,nz)
     1   ,phi(nx,ny,nz),om(nx,ny,nz),omb(nx,ny,nz)
     1   ,sh(nx,ny,nz)
     1   ,rri(nz),     rrj(nz), ter(nx,ny)
      real p(nz),zter,aa,bb
      real ubias,vbias,urms,vrms
      real udrop(nz),vdrop(nz),tdrop(nz)  
      real sm(10),st(10) 
c                                        
      real ffz,rm,slastu,slastv,slastt,cor
      character*9 dum,a9_time

c  variables: udrop vdrop tdrop are the u,v,T from the dropsonde
c             u,v are the LAPS analyzed velocities
c             ub,vb are the background grids
c             us,vs,oms are the turbulent velocity components
c  create a "truth " profile using nongaussian distribution of 
c dropsonde and model error in combined pdf
c We now save each truth profile so arrays are 2D

      real ut(nz,nx),vt(nz,nx),Tt(nz,nx),wt(nz,nx),pht(nz)
      logical lfndref,l_interp
      integer alt

      nsave=0
      l_interp=.true.
      iii=i4time
      call make_fnam_lp(i4time,a9_time,istatus)
      call zero3d(erru,nx,ny,nz)
      call zero3d(errv,nx,ny,nz)
      call zero3d(errt,nx,ny,nz)
      call zero3d(errw,nx,ny,nz)

      call get_r_missing_data(smsng,istatus)
      r=287.04
      g=9.808
      cnt=0

c we need to use non guassian distributions for truth estimates
c this routine creates two functions for the rejection method: one
c gaussian, one non gaussian. This sets up function rm
c recover variable profiles and put them in the error arrays 

      do alt=0,2

       if(l_interp)then

        l_interp=.false.  !interpolations are only required first time

        do k=1,nz

         if(rri(k).ne.smsng .and. rrj(k).ne.smsng)then
           ii=nint(rri(k))
           jj=nint(rrj(k))
           if(ii.ge.nx)then
              ii=nx-1
              print*,'Warning: dropsonde point on S edge of domain'
           endif
           if(jj.ge.ny)then
              jj=ny-1
              print*,'Warning: dropsonde point on N edge of domain'
           endif
           aa=rri(k)-float(ii)
           bb=rrj(k)-float(jj)
c terrain
           if(nsave .eq. 0)then
              zter=ter(ii,jj)*(1.-aa)*(1.-bb)
     1            +ter(ii+1,jj)*aa*(1.-bb)
     1            +ter(ii,jj+1)*bb*(1.-aa)
     1            +ter(ii+1,jj+1)*aa*bb

           endif
           nsave=k
         else
           if(nsave.eq.0)then
              zter=ter(ii,jj)*(1.-aa)*(1.-bb)
     1            +ter(ii+1,jj)*aa*(1.-bb)
     1            +ter(ii,jj+1)*bb*(1.-aa)
     1            +ter(ii+1,jj+1)*aa*bb
           endif
           ii=nx/2
           jj=ny/2
           aa=0.
           bb=0.
         endif

c these are interpolated (to dropsonde point)
c  profiles of turbulence and analysis
         erru(1,6,k)=us(ii,jj,k)*(1.-aa)*(1.-bb)
     1   +us(ii+1,jj,k)*aa*(1.-bb)+
     1   us(ii,jj+1,k)*bb*(1.-aa)+us(ii+1,jj+1,k)*aa*bb
         erru(1,5,k)=u(ii,jj,k)*(1.-aa)*(1.-bb)
     1+u(ii+1,jj,k)*aa*(1.-bb)+
     1 u(ii,jj+1,k)*bb*(1.-aa)+u(ii+1,jj+1,k)*aa*bb
         errv(1,6,k)=vs(ii,jj,k)*(1.-aa)*(1.-bb)
     1+vs(ii+1,jj,k)*aa*(1.-bb)+
     1 vs(ii,jj+1,k)*bb*(1.-aa)+vs(ii+1,jj+1,k)*aa*bb
         errv(1,5,k)=v(ii,jj,k)*(1.-aa)*(1.-bb)
     1+v(ii+1,jj,k)*aa*(1.-bb)+
     1 v(ii,jj+1,k)*bb*(1.-aa)+v(ii+1,jj+1,k)*aa*bb
         errt(1,6,k)=t(ii,jj,k)*(1.-aa)*(1.-bb)
     1+t(ii+1,jj,k)*aa*(1.-bb)+
     1 t(ii,jj+1,k)*bb*(1.-aa)+t(ii+1,jj+1,k)*aa*bb
c density
         errt(1,7,k)=p(k)/r/errt(1,6,k)          
c specific humidity
         errt(1,8,k)=sh(ii,jj,k)*(1.-aa)*(1.-bb)
     1+sh(ii+1,jj,k)*aa*(1.-bb)+
     1 sh(ii,jj+1,k)*bb*(1.-aa)+sh(ii+1,jj+1,k)*aa*bb
c height
         errt(1,5,k)=phi(ii,jj,k)*(1.-aa)*(1.-bb)
     1+phi(ii+1,jj,k)*aa*(1.-bb)+
     1 phi(ii,jj+1,k)*bb*(1.-aa)+phi(ii+1,jj+1,k)*aa*bb
c dropsonde
         erru(1,9,k)=udrop(k)
         errv(1,9,k)=vdrop(k)
         errt(1,9,k)=tdrop(k)
c background
         erru(1,10,k)=ub(ii,jj,k)*(1.-aa)*(1.-bb)
     1+ub(ii+1,jj,k)*aa*(1.-bb)+
     1 ub(ii,jj+1,k)*bb*(1.-aa)+ub(ii+1,jj+1,k)*aa*bb
         errv(1,10,k)=vb(ii,jj,k)*(1.-aa)*(1.-bb)
     1+vb(ii+1,jj,k)*aa*(1.-bb)+
     1 vb(ii,jj+1,k)*bb*(1.-aa)+vb(ii+1,jj+1,k)*aa*bb
         errt(1,10,k)=tb(ii,jj,k)*(1.-aa)*(1.-bb)
     1+tb(ii+1,jj,k)*aa*(1.-bb)+
     1 tb(ii,jj+1,k)*bb*(1.-aa)+tb(ii+1,jj+1,k)*aa*bb
         errw(1,10,k)=-r*errt(1,10,k)/p(k)/g
     1*(omb(ii,jj,k)*(1.-aa)*(1.-bb)
     1+ omb(ii+1,jj,k)*aa*(1.-bb)+
     1  omb(ii,jj+1,k)*bb*(1.-aa)+omb(ii+1,jj+1,k)*aa*bb)

c vertical motion w
         errw(1,5,k)=-r*errt(1,6,k)/p(k)/g*(om(ii,jj,k)*(1.-aa)*(1.-bb)
     1+om(ii+1,jj,k)*aa*(1.-bb)+
     1 om(ii,jj+1,k)*bb*(1.-aa)+om(ii+1,jj+1,k)*aa*bb)
c turb vv
         errw(1,6,k)=(oms(ii,jj,k)*(1.-aa)*(1.-bb)
     1+oms(ii+1,jj,k)*aa*(1.-bb)+
     1 oms(ii,jj+1,k)*bb*(1.-aa)+oms(ii+1,jj+1,k)*aa*bb)

        enddo ! on k

       endif  ! interpolations
c
c --------------------------------------------------------------------
c ---------------- now create truth profiles -------------------------
c --------------------------------------------------------------------
c
c Alternative 2 is based on use of the mean truth profile. Thus, set
c the analysis err arrays equal to the mean truth.

       if(alt.eq.2)then
c        do k=1,nz
            erru(1,5,:)=ut(:,1)
            errv(1,5,:)=vt(:,1)
            errt(1,6,:)=Tt(:,1)
            errw(1,5,:)=Wt(:,1) 
c        enddo
       endif

       do n=4,nx ! these are the random estimates for the statistics
c this creates nx-3 unique truth profiles: nx is used because it is the 
c first dimension of the err arrays 
c  do background errors first
c set up rejection method non-gaussian variable function
        mmg = 2
        mg=10
        do k=1,nz
c        print*,'K ',k
         sm(1)=udrop(k)
         sm(2)=erru(1,10,k)   !ub(ii,jj,k)
         st(1)=oberu
         st(2)=moderu  
         if(udrop(k).ne.smsng)ut(k,n)=rm(sm,st,mg,mmg,iii)
         sm(1)=vdrop(k)
         sm(2)=errv(1,10,k)   !vb(ii,jj,k)
         st(1)=oberu
         st(2)=moderv  
         if(vdrop(k).ne.smsng)vt(k,n)=rm(sm,st,mg,mmg,iii)      
         sm(1)=tdrop(k)
         sm(2)=errt(1,10,k)   !tb(ii,jj,k)
         st(1)=obert
         st(2)=modert
         if(tdrop(k).ne.smsng)Tt(k,n)=rm(sm,st,mg,mmg,iii)      
         Wt(k,n)=oberw*ffz(iii,20,cor,0.)

c for each truth profile compute error and store in the err arrays by
c height level (k) and truth number (n)

         berru=erru(1,10,k)-ut(k,n)
         berrv=errv(1,10,k)-vt(k,n)
         berrt=errt(1,10,k)  !tb(ii,jj,k)*(1.-aa)*(1.-bb)+tb(ii+1,jj,k)*aa*(1.-bb)+
c                         1   tb(ii,jj+1,k)*bb*(1.-aa)+tb(ii+1,jj+1,k)*aa*bb
         errw(n,1,k)=errw(1,10,k)-wt(k,n)
         berrt=p(k)/r*(1./berrt-1./Tt(k,n))
         erru(n,1,k)=berru
         errv(n,1,k)=berrv
         errt(n,1,k)=berrt
        enddo! on k
       enddo! on n
c
c now analysis errors with a different ensemble of truth estimates
c ----------------------------------------------------------------
       do n=4,nx
        do k=1,nz
c create truth profiles using non-guassian random numbers
c        print*,'K ',k
         sm(1)=udrop(k)
         sm(2)=erru(1,5,k)  !ub(ii,jj,k)
         st(1)=oberu
         st(2)=moderu
         if(udrop(k).ne.smsng)ut(k,n)=rm(sm,st,mg,mmg,iii)
         sm(1)=vdrop(k)
         sm(2)=errv(1,5,k)  !vb(ii,jj,k)
         st(1)=oberu
         st(2)=moderv  
         if(vdrop(k).ne.smsng)vt(k,n)=rm(sm,st,mg,mmg,iii)      
         sm(1)=tdrop(k)
         sm(2)=errt(1,6,k)  !errt(1,5,k)  !tb(ii,jj,k)
         st(1)=obert
         st(2)=modert  
         if(tdrop(k).ne.smsng)Tt(k,n)=rm(sm,st,mg,mmg,iii)      
         Wt(k,n)=oberw*ffz(iii,20,cor,0.)
c integrate hypsometric eqn to recover heights
         if(zter.lt.errt(1,5,k)) then! begin to integrate truth ht
          tave=.5*(Tt(k,n)+errt(1,8,k)/6.+Tt(k-1,n)+errt(1,8,k)/6.)! virt temp
          pht(k)=pht(k-1)+r*tave/g*alog(p(k-1)/p(k))
         else
          kstart=k+1
          pht(k)=errt(1,5,k)
         endif

         aerru=erru(1,5,k)-ut(k,n)
         aerrv=errv(1,5,k)-vt(k,n)
         aerrtt=errt(1,6,k)

c        t(ii,jj,k)*(1.-aa)*(1.-bb)+t(ii+1,jj,k)*aa*(1.-bb)+
c    1   t(ii,jj+1,k)*bb*(1.-aa)+t(ii+1,jj+1,k)*aa*bb ! temp at drop
c convert omega errors to w errors
c        errw(n,2,k)=-r*aerrtt/p(k)/g*(om(ii,jj,k)*(1.-aa)*(1.-bb)
c    1   +om(ii+1,jj,k)*aa*(1.-bb)+
c    1   om(ii,jj+1,k)*bb*(1.-aa)+om(ii+1,jj+1,k)*aa*bb)-wt(k)

         errw(n,2,k)=errw(1,5,k)-wt(k,n)
         aerrt=p(k)/r*(1./aerrtt-1./Tt(k,n))!density err
         erru(n,2,k)=aerru
         errv(n,2,k)=aerrv
         errt(n,2,k)=aerrt
         if(errt(1,5,k).gt.zter) then
           errt(n,4,k)=(errt(1,5,k)-pht(k))  ! ht error
         endif
         errt(n,3,k)=g*p(k)/r/aerrtt*errt(n,4,k) ! pressure error from hts               
        enddo! on k
       enddo ! on n
c
c ----------------------------------------------------------------
c now compute bias error by summing over all nx-4 scenarios
       sumbiasu=0
       sumbiasv=0
       cntt=0
       do k=1,nz   
        cnt=0
        sumbeu=0
        sumbev=0
        sumbet=0
        sumaeu=0
        sumaev=0
        sumaet=0
        sumbew=0
        sumaew=0
        sumapp=0
        sumph=0.
        sumut=0.
        sumvt=0.
        sumTt=0.
        sumWt=0.
        if(alt.eq.0.or.alt.eq.2)then  !original method (alt=0: non-zero mean bias)
                                      !mean truth method (alt=2: zero mean bias).
           do n=4,nx
             cnt=cnt+1.
             sumbeu=erru(n,1,k)+sumbeu
             sumbev=errv(n,1,k)+sumbev
             sumbet=errt(n,1,k)+sumbet
             sumaeu=erru(n,2,k)+sumaeu
             sumaev=errv(n,2,k)+sumaev
             sumaet=errt(n,2,k)+sumaet
             sumbew=errw(n,1,k)+sumbew
             sumaew=errw(n,2,k)+sumaew
             sumapp=errt(n,3,k)+sumapp
             sumph=errt(n,4,k)+sumph ! ht error
             if(alt.eq.0)then
                sumut=ut(k,n)+sumut
                sumvt=vt(k,n)+sumvt
                sumTt=Tt(k,n)+sumTt
                sumWt=Wt(k,n)+sumWt
             endif
          enddo ! on n

        else     ! alternative 1: zero-mean bias
          cnt=1
        endif

c     bias estimates and mean truth profile
        if(cnt.gt.0)then
          erru(1,1,k)=sumbeu/cnt
          errv(1,1,k)=sumbev/cnt
          errt(1,1,k)=sumbet/cnt
          erru(1,2,k)=sumaeu/cnt
          errv(1,2,k)=sumaev/cnt
          errt(1,2,k)=sumaet/cnt !density
          errw(1,1,k)=sumbew/cnt
          errw(1,2,k)=sumaew/cnt
          errt(1,3,k)=sumapp/cnt! pressure
          errt(1,4,k)=sumph/cnt !height
          if(alt.eq.0)then
             ut(k,1)=sumut/cnt  !save mean profiles in element (*,1)
             vt(k,1)=sumvt/cnt
             Tt(k,1)=sumTt/cnt
             Wt(k,1)=sumWt/cnt
          endif
        endif
        if(k.le.nsave.and.k.ge.kstart) then !sum from kstart to nsave            
          cntt=cntt+1
          sumbiasu=erru(1,2,k)+sumbiasu
          sumbiasv=errv(1,2,k)+sumbiasv
        endif
       enddo ! on k
       print*, 'Total u,v bias ', sumbiasu, sumbiasv
       if(cntt.gt.0)then
          print*, 'Layer u,v bias ',sumbiasu/cntt, sumbiasv/cntt
       endif
c these are mean bias errors for column
       erru(1,11,1)=sumbiasu
       errv(1,11,1)=sumbiasv
      
c now compute variance
       do k=1,nz
        sumbeu=0
        sumbev=0
        sumbet=0
        sumbew=0
        sumaeu=0
        sumaev=0
        sumaet=0
        sumaew=0
        sumapp=0.
        sumph=0
        cnt=0.
       
        do n=4,nx ! for all truth senarios
         cnt=cnt+1.
         sumbeu=sumbeu+(erru(n,1,k)-erru(1,1,k))**2 !backgrnd error variance
         sumaeu=sumaeu+(erru(n,2,k)-erru(1,2,k))**2 !analysis error variance
         sumbev=sumbev+(errv(n,1,k)-errv(1,1,k))**2 
         sumaev=sumaev+(errv(n,2,k)-errv(1,2,k))**2
         sumbet=sumbet+(errt(n,1,k)-errt(1,1,k))**2 
         sumaet=sumaet+(errt(n,2,k)-errt(1,2,k))**2
         sumbew=sumbew+(errw(n,1,k)-errw(1,1,k))**2
         sumaew=sumaew+(errw(n,2,k)-errw(1,2,k))**2
         sumapp=sumapp+(errt(n,3,k)-errt(1,3,k))**2
         sumph =sumph+(errt(n,4,k)-errt(1,4,k))**2
c covariance calc err ( ,6, ) is the turb component
         erru(3,2,k)=(erru(n,2,k)-erru(1,2,k))*(erru(n,6,k)
     1    *ffz(iii,20,cor,0.))
     1    +erru(3,2,k)
         errv(3,2,k)=(errv(n,2,k)-errv(1,2,k))*(errv(n,6,k)
     1    *ffz(iii,20,cor,0.))
     1    +errv(3,2,k)
         errw(3,2,k)=(errw(n,2,k)-errw(1,2,k))*(errw(n,6,k)
     1    *ffz(iii,20,cor,0.))
     1    +errw(3,2,k)
        enddo! on n
c   variance estimates
        if(cnt.gt.0)then
          erru(2,1,k)=(sumbeu/cnt)
          errv(2,1,k)=(sumbev/cnt)
          errt(2,1,k)=(sumbet/cnt)
          erru(2,2,k)=(sumaeu/cnt)
          errv(2,2,k)=(sumaev/cnt)
          errt(2,2,k)=(sumaet/cnt)! density
          errw(2,1,k)=(sumbew/cnt)
          errw(2,2,k)=(sumaew/cnt)
          errt(2,3,k)=(sumapp/cnt)! pressure
          erru(3,2,k)=(erru(3,2,k)/cnt)
          errv(3,2,k)=(errv(3,2,k)/cnt)
          errw(3,2,k)=(errw(3,2,k)/cnt)
          errt(2,4,k)=(sumph/cnt)! height
        endif
       enddo ! on k

       call write_errors(a9_time,p,erru,errv,errw,errt
     1,nx,ny,nz,rri,rrj,zter,alt)

      enddo  ! on number of alternatives

      return
      end 
c
c------------------------------------------------
c
      Subroutine write_errors(a9_time,p,erru,errv,errw,errt
     1,nx,ny,nz,rri,rrj,zter,alt)
c this routine writes an ascii file to the output log summarizing errors
      integer lun
      integer alt
      real p(nz)
      real erru(nx,ny,nz),errv(nx,ny,nz),errt(nx,ny,nz)
      real errw(nx,ny,nz)
      real rri(nz),rrj(nz)
      character*200 cfname_out
      character*9   a9_time
c
c direct output to lapsprd/balance/air/yyjjjhhmm.air -- ascii output
c
      call get_directory('balance',cfname_out,lend)

      if(alt==0)then
         cfname_out=cfname_out(1:lend)//'air/'//a9_time//'.air'
         call s_len(cfname_out,lend)
         print*,'Airdrop output filename: ',cfname_out(1:lend)
         lun = 20
         open(lun,file=cfname_out,form='formatted',status='unknown'
     +,err=909)

      elseif(alt==1)then
         cfname_out=cfname_out(1:lend)//'air/'//a9_time//'.airp1'
         call s_len(cfname_out,lend)
         print*,'Airdrop output filename: ',cfname_out(1:lend)
         lun = 20
         open(lun,file=cfname_out,form='formatted',status='unknown'
     +,err=909)

      elseif(alt==2)then
         cfname_out=cfname_out(1:lend)//'air/'//a9_time//'.airp2'
         call s_len(cfname_out,lend)
         print*,'Airdrop output filename: ',cfname_out(1:lend)
         lun = 20
         open(lun,file=cfname_out,form='formatted',status='unknown'
     +,err=909)

      endif

      r=287.04
      sum=0
      sum1=0
      sum2=0

      write(lun,1011)
 1011 format(6x,'************AIRDROP OUTPUT**************')
      write(lun,1012)
 1012 format(6x,'AIRDROP PROFILES FROM LAPS ANALYSIS')
      write(lun,1013)
 1013 format(4x,'P',6x,'HT',6x,'U',8x,'V',8x,'W',8x,'T',7x,'Den')
      write(lun,1014)
 1014 format(3x,'mb',7x,'m',5x,'m/sec',4x,'m/sec',4x,'m/sec',6x,'K'
     +,6x,'kg/m3')
      do k=1,nz
       ii=rri(k)
       jj=rrj(k)
       if(zter.lt.errt(1,5,k)) then
        write(lun,1005) p(k)/100., errt(1,5,k),erru(1,5,k),
     1  errv(1,5,k),errw(1,5,k),errt(1,6,k),errt(1,7,k)
       endif
      enddo
      write(lun,1015) 
       write(lun,1022)
 1022 format(6x,'DROPSONDE AND BACKGROUND WIND PROFILES  ')
      write(lun,1023)
 1023 format(4x,'P',6x,'HT',6x,'UD',7x,'VD',7x,'UB',7x,'VB'
     & ,7x,'TD',7x,'TB')
      write(lun,1024)
 1024 format(3x,'mb',7x,'m',5x,'m/sec',4x,'m/sec',4x,'m/sec',4x,'m/sec'
     +,6x,'K',6x,'K')
     
      do k=1,nz
       ii=rri(k)
       jj=rrj(k)
       if(zter.lt.errt(1,5,k)) then
        write(lun,1025) p(k)/100., errt(1,5,k),erru(1,9,k),
     1  errv(1,9,k),erru(1,10,k),errv(1,10,k),errt(1,9,k),errt(1,10,k)          
       endif
      enddo
      write(lun,1015) 
 1015 format(/)
      write(lun,1016)
 1016 format(8x, 'BIAS ERROR ESTIMATES')
      write(lun,1082)
      write(lun,1083)
      write(lun,1084)
      do k=1,nz
       if(zter.lt.errt(1,5,k)) then
        write(lun,1086) p(k)/100., errt(1,5,k),erru(1,1,k),
     1  errv(1,1,k),errw(1,1,k),
     1  erru(1,2,k),errv(1,2,k),errw(1,2,k),sum,sum,sum
       endif
      enddo
      write(lun,1015)
      write(lun,1016)
      write(lun,1017)
 1017 format(1x,'P LVL  HT LVL ')
      write(lun,1009)
      write(lun,1010)
      do k=1,nz
       if(zter.lt.errt(1,5,k)) then
       write(lun,1008) p(k)/100., errt(1,5,k),errt(1,3,k)/100.
     &,errt(1,4,k),  errt(1,1,k),errt(1,2,k)
     
       endif
      enddo
   
      write(lun,1015)
      write(lun,1028)
 1028 format(10x,'TOTAL COLUMN WIND BIAS')
      write(lun,1029) erru(1,11,1), errv(1,11,1)
 1029 format(1x,'U BIAS ',f7.2,' V BIAS ',f7.2)
      write(lun,1015)
      write(lun,1018)
 1018 format(10x,'VARIANCE ESTIMATES')
      write(lun,1002)
      write(lun,1003)
      write(lun,1004)

      do k=1,nz
       ii=rri(k)
       jj=rrj(k)
       if(zter.lt.errt(1,5,k)) then
        write(lun,1000) p(k)/100., errt(1,5,k),erru(2,1,k),
     1  errv(2,1,k),errw(2,1,k),
     1  erru(2,2,k),errv(2,2,k),errw(2,2,k),erru(1,6,k)**2,
     1  errv(1,6,k)**2,
     1  errw(1,6,k)**2,errt(2,1,k),errt(2,2,k)
       endif
      enddo

      write(lun,1015)
      write(lun,1019)
 1019 format(3x,'TOTAL VARIANCE ESTIMATES')
      write(lun,1020)
 1020 format( '  P     HT     Analysis + Turbulence ')
      write(lun,1021)
 1021 format( '                    (m2/sec2)  ' )
      write(lun,1007)
 1007 format(19x,'U',8x,'V',9x,'W',6x,'Den',7x,'P',9x,'HT') 
 1009 format(19x,'P',8x,'HT',8x,'Density')
 1010 format(19x,'mb',8x,'m',8x,'Back',4x,'Anal')

      do k=1,nz
       ii=rri(k)
       jj=rrj(k)
c combined variance formula - err(3,2,k) contains covariance
       sum=erru(2,2,k)+erru(1,6,k)**2+2.*erru(3,2,k)
       sum1=errv(2,2,k)+errv(1,6,k)**2+2.*errv(3,2,k)
       sum2=errw(2,2,k)+errw(1,6,k)**2+2.*errw(3,2,k)
       if(zter.lt.errt(1,5,k)) then
        write(lun,1001) p(k)/100., errt(1,5,k),
     1 sum,sum1,sum2
     1,errt(2,2,k),errt(2,3,k)/10000.,errt(2,4,k)
       endif
      enddo

      close(lun)

      return

 1000 format(1x,f5.0,f8.0,2(f5.1,f5.1,f5.3),3f5.2,2f8.6)
 1006 format(1x,f5.0,f8.0,2(f5.1,f5.1,f5.3),3f5.2,2f7.4)
 1086 format(1x,f5.0,f8.0,2(f7.1,f7.1,f7.3),3f7.2)
 1002 format(18x,'Background',5x,'Analysis',7x,'Turbulence',5x,
     1'Density')
 1082 format(18x,'Background',12x,'Analysis',12x,'Turbulence')
 1083 format(20x,3('(m/sec)',13x))
 1003 format(20x,3('(m2/sec2)',5x),'(kg2/m6)')
 1004 format(3x,'P',7x,'HT  ',3('  U    V    W  '),'  Bgnd     Anal')
 1084 format(3x,'P',7x,'HT  ',3('   U      V      W    '))
 1001 format(1x,f5.0,f8.0,3f9.3,2f9.6,f9.2)
 1005 format(1x,f5.0,f8.0,5f9.3)
 1025 format(1x,f5.0,f8.0,6f9.3)
 1008 format(1x,f5.0,f8.0,4f9.3)

  909 print*,'Error opening airdrop output file'
      return
      end
