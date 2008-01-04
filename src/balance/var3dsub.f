c     call var3d (u,v,d,t,up,vp,dp,smsng,dx,dt,dpht,distnf,p,ht
c    &,errwnds,errdds,errmod,io,jo,nx,ny,nz,slat,slon,ter,a9_time)
c
c *** subroutine extracted from test var3d program 1-4-06 JRS.
c *** removed up,vp,dp from argument list as these are not needed
c *** as return variables in calling routine.
c

c     subroutine var3d (u,v,d,t,up,vp,dp,smsng,dx,dt,dpht,distnf,
c    &p,ht,errwnds,errdds,errmod,io,jo,nx,ny,nz,slat,slon,ter)

      subroutine var3d (u,v,d,t,smsng,dx,dt,dpht,distnf,p,ht
     &,errwnds,errdds,errmod,errdis,io,jo,nx,ny,nz,slat,slon,
     &ter,a9_time)

      implicit none
     
c this subroutine takes in fields of u,v, and density and estimates a
c variance over the grid by determining a weighted difference between
c each grid point and a central profile up,vp,dp. The weighting function is
c a skewed gaussian based on the mean wind
c input  variables
c             u(nx,ny,nz) 3D grid of winds
c             v(nx,ny,nz)     "
c             t(nx,ny,nz) 3D grid of potl temperatures
c             d(nx,ny,nz) density
c             up(nz) U-wind profile at drop point
c             vp(nz) V-wind profile
c             dp(nz) density   profile
c output fields
c             upv(nz) estimated profile variance
c             vpv(nz)
c             dpv(nz)
c internal fields
c             wt(nx,ny)
c             wt3(nx,ny,nz)        
 
       integer ns
       parameter (ns=4)
       integer isave(ns),jsave(ns)
       integer nse
       parameter (nse=2**ns-1)
       integer index(nse),inxx(nse)	
       integer inx(nse)
       character*9 a9_time
       character*24 lable,title
       integer nx, ny, nz
       real, intent(in) :: u(nx,ny,nz),v(nx,ny,nz),d(nx,ny,nz)
     &,t(nx,ny,nz)
       real, intent(in) :: slat(nx,ny),slon(nx,ny),ter(nx,ny)

       real upi(nse),vpi(nse),dpi(nse),twi(nse)
       real, allocatable, dimension(:,:) :: wt
       real, allocatable, dimension(:,:,:) :: wt3
       real wtsave(ns),drpx(ns),drpy(ns)
       real alph,beta,sumu,sumv,sumd
       real p(nz),ht(nz)
       real ud(nz),vd(nz),ubb(nz),vbb(nz),wbb(nz),uab(nz),vab(nz)
       real up(nz),vp(nz),dp(nz),tp(nz),dd(nz)
       real upv(nz),vpv(nz),tpv(nz),wp(nz),dpv(nz),
     &      wab(nz),dbb(nz),dab(nz),ubv(nz),vbv(nz),wbv(nz),
     &      uav(nz),vav(nz),wav(nz),utv(nz),vtv(nz),wtv(nz),
     &      dav(nz),dbv(nz),pav(nz),htav(nz),tav(nz)

       real eponexy,gamma,sumw,ca,ddd2,ct,epone
       real ss2,sa,sumt,tpvs,sn2,sumww,sumdd,zter,smax,sumvv,sumuu
       real errmod,dt,dx,smsng,dpht,errdds,errwnds,distnf,st,ubar
       real vbar,vbtot,ubtot,errdis
       real ri,rj,ter_hgt_drop,droplat,droplon

       integer i,j,k,kk,l,mm,lsave,n,lmax,m,in,jo,io, ncnt,ix,jy
       integer ktop,kbot,nofly2,idist2
       integer user_sonde_loc(2)
       integer nofly_area(nx,ny)
       integer istatus

       data index /1,2,3,4,12,13,14,23,24,34,123,124,134,234,1234/

       print*,'Module var3d'
       print*,'nx/ny/nz:      ',nx,ny,nz
       print*,'grid center:   ',io,jo
       print*,'nofly radius:  ',distnf
       print*,'drop height:   ',dpht
       print*,'delta t (sec): ',dt
        
c find kbot and ktop

c Note: This assumes the grid center of the domain is the drop
c       location for the PADS mission. Ie., drop location is specified
c       by vars grid_cen_lat_cmn & grid_cen_lon_cmn in nest7grid.parms

       call get_grid_center(droplat,droplon,istatus)

       call latlon_to_rlapsgrid(droplat,droplon,slat,slon,nx,ny,ri,rj
     1                                ,istatus)
       call bilinear_laps(ri,rj,nx,ny,ter,ter_hgt_drop)

c   compute ktop
       do k=1,nz
        if (dpht.lt.ht(k)) then
         ktop=k
         go to 5
        endif
       enddo
c   compute kbot
5      do k=1,nz
        if (ht(k).ge.ter_hgt_drop) then
         kbot=k
         go to 6
        endif
       enddo

       if(kbot.gt.ktop)then
          print*,'Error: kbot > ktop in var3dsub '
          print*,'Terminating ', kbot,ktop
          stop
       endif

c determine mean wind and compute profiles

6      sumu=0.
       sumv=0.
       ncnt=0
       do k=kbot,ktop 
        do j=1,ny
         do i=1,nx
          if (u(i,j,k).ne.smsng) then
           sumu=sumu+u(i,j,k)
           sumv=sumv+v(i,j,k)
           ncnt=ncnt+1
          endif
         enddo
        enddo
       enddo
       ubar=sumu/float(ncnt)
       vbar=sumv/float(ncnt)
       if(ubar.lt.1.and.ubar.ge.0.) ubar=1
       if(ubar.gt.-1.and.ubar.lt.0.) ubar=-1
       if(vbar.lt.1.and.vbar.ge.0.) vbar=1
       if(ubar.gt.-1.and.ubar.lt.0.) vbar=-1
       print*,'ubar,vbar ',ubar,vbar

c profiles
       ubtot=0.
       vbtot=0.
       print*, 'ktop ', ktop
       print*, 'kbot ', kbot
       print*
       print*, 'Height profile in drop layer'
       do k=kbot,ktop
        print*,'k, ht(m) ',k,ht(k)
        up(k)=u(io,jo,k)
        vp(k)=v(io,jo,k)
        wp(k)=0.
        dp(k)=d(io,jo,k)
        tp(k)=t(io,jo,k)
        ud(k)=up(k)
        vd(k)=vp(k)
        dd(k)=dp(k)
        ubb(k)=0.
        vbb(k)=0.
        wbb(k)=0.
        uab(k)=0.
        vab(k)=0.
        wab(k)=0
        dab(k)=0.
        dbb(k)=0.
        ubv(k)=errmod**2
        vbv(k)=errmod**2
        wbv(k)=errmod**2/10000.
        utv(k)=(.05*up(k))**2
        vtv(k)=(.05*vp(k))**2
        wtv(k)=utv(k)+vtv(k)
        dav(k)=(dp(k)*.01)**2
        pav(k)=(p(k)*.0001)**2
        htav(k)=pav(k)/dp(k)*9.8
       enddo
c compute weight function
       allocate (wt(nx,ny))
c   alpha is the along flow weight tuning parameter based on expected model 
c   phase error in meters (errdis), the grid distance dx (m).D=errdis/dx 
c   alpha is 1/D**2 making the weight 1/e for distance D (in I,J space)
c   beta is the cross flow fall off value and is set to 3*alpha  

       alph=(dx/errdis)**2     
       beta=3.*alph
       ix=ubar*dt/dx+.5
       jy=vbar*dt/dx+.5     
       print*, 'full wind offset i,j', ix,jy
       st=vbar/sqrt(ubar**2+vbar**2)
       ct=ubar/sqrt(ubar**2+vbar**2)
       do j=1,ny
        do i=1,nx
         ddd2=float(i-io+ix)**2+float(j-jo+jy)**2 
         if(ddd2.eq.0.) then
          epone=0.
          go to 33
         endif
         sa=float(j-jo+jy)*st  
         ca=float(i-io+ix)*ct 
         ss2=(sa+ca)**2              
         sn2=ddd2-ss2                 
         epone= alph*ss2+beta*sn2         
  33     wt(i,j)=exp(-epone)
        enddo
       enddo

c compute variance estimate
       allocate (wt3(nx,ny,nz))
       tpvs=0.
       do k=kbot,ktop
        sumu=0.
        sumv=0.
        sumt=0.
        sumw=0.
        do j=1,ny
         do i=1,nx
          if (u(i,j,k).ne.smsng) then
           wt3(i,j,nz)=(wt(i,j)*(u(i,j,k)-up(k)))**2
           sumu=sumu+wt3(i,j,nz)
           wt3(i,j,nz-1)=(wt(i,j)*(v(i,j,k)-vp(k)))**2
           sumv=sumv+wt3(i,j,nz-1)
           wt3(i,j,nz-2)=wt3(i,j,nz)+wt3(i,j,nz-1)
           sumd=sumd+(wt(i,j)*(d(i,j,k)-dp(k)))**2
           sumw=sumw+wt(i,j)
          endif
         enddo 
        enddo
c variances for forecast quality
        upv(k)=sumu/sumw !u variance due to model error at drop location
        vpv(k)=sumv/sumw !v variance due to model error at drop location
        tpv(k)=upv(k)+vpv(k) ! total variance 
        if(k.le.ktop)tpvs=tpvs+tpv(k)
        dpv(k)=sumd/float(ncnt)
        write(6,1000) 'uv,vv,tot,dv ',upv(k),vpv(k),tpv(k),dpv(k)
 1000   format(1x,a13,4f9.5)
       enddo 
       print*,'total variance', tpvs

c sonde emmulator
c recompute weights

c assume every grid location is a sonde
c determine impact on drop location
c recompute ubar and vbar below ktop
       sumu=0.
       sumv=0.
       ncnt=0
       do k=kbot,ktop 
        do j=1,ny
         do i=1,nx
          if (u(i,j,k).ne.smsng) then
           sumu=sumu+u(i,j,k)
           sumv=sumv+v(i,j,k)
           ncnt=ncnt+1
          endif
         enddo
        enddo
       enddo
       ubar=sumu/float(ncnt) ! mean wind through drop layer
       vbar=sumv/float(ncnt)
       if(ubar.lt.1.and.ubar.ge.0.) ubar=1
       if(ubar.gt.-1.and.ubar.lt.0.) ubar=-1
       if(vbar.lt.1.and.vbar.ge.0.) vbar=1
       if(ubar.gt.-1.and.ubar.lt.0.) vbar=-1
       print*,'ubar,vbar below ktop ',ubar,vbar
c the tuning parameters alph,beta, and gamma now mimic an isentropic
c elliptical weighting function so values are different from previous 
c using correlation distance as a scaling parameter
       alph=(2.*dx/errdis)**2
       beta=2.*alph
       gamma=0.001
       ix=ubar*dt/dx+.5
       jy=vbar*dt/dx+.5     
       print*, 'offset i,j in drop layer', ix,jy
       st=vbar/sqrt(ubar**2+vbar**2)
       ct=ubar/sqrt(ubar**2+vbar**2)
       do j=1,ny
        do i=1,nx
         ddd2=float(i+ix-io)**2+float(j+jy-jo)**2 
         if(ddd2.eq.0.) then
          eponexy=0.
          go to 333
         endif
         sa=float(j+jy-jo)*st 
         ca=float(i+ix-io)*ct  
         ss2=(sa+ca)**2                
         sn2=ddd2-ss2                 
         eponexy= alph*ss2+beta*sn2         
 333     do k=kbot,nz     
          do kk=ktop,kbot,-1
             if(u(i,j,k).ne.smsng) then
                epone=eponexy+gamma*(t(i,j,k)-tp(kk))**2
                wt3(i,j,k)=exp(-epone)
             else
                wt3(i,j,k)=smsng
             endif
          enddo
         enddo
        enddo
       enddo
c now add up all the weights below ktop
       do j=1,ny
        do i=1,nx
         wt(i,j)=0.
         ncnt=0
         do k=kbot,ktop
          if(wt3(i,j,k).ne.smsng) then
           wt(i,j)=wt(i,j)+wt3(i,j,k)
           ncnt=ncnt+1
          endif
         enddo
         wt(i,j)=wt(i,j)/float(ncnt)
        enddo
       enddo
c now find maximum wt for area
c read PADS input of nofly area file and user sonde location information
      call read_PADS_nofly_and_sonde(nx,ny
     1,nofly_area,user_sonde_loc,istatus)

      idist2=(6.0*dt/dx)**2
      nofly2=(distnf/dx)**2
      do l=1,ns
       smax=0.
       do j=1,ny
        do i=1,nx
         if(wt(i,j).gt.smax) then
          if(((i-io)**2+(j-jo)**2).lt.nofly2) go to 1
          do m=1,l-1
           if(((isave(m)-i)**2+(jsave(m)-j)**2).lt.idist2) go to 1
          enddo
          smax=wt(i,j)
          isave(l)=i
          jsave(l)=j
          wtsave(l)=smax
         endif
    1   enddo
       enddo
       write(6,1001) 'l,is,js,wt ',l,isave(l),jsave(l),wtsave(l)
       drpx(l)=isave(l)
       drpy(l)=jsave(l)
 1001  format(1x,a11, 3i6,f9.5)
      enddo
c compute variance for each sonde set or scenario
      lmax=2**ns-1
      do l=1,lmax   
       call decode_index(l,index,in,inx,lmax,ns)
       sumu=0.
       sumv=0.
       sumd=0.
       sumw=0.
       do k=kbot,ktop
        sumuu=0.
        sumvv=0.
        sumdd=0.
        sumww=0.
        do m=1,in
         mm=inx(m)
         sumuu=sumuu+wt3(isave(mm),jsave(mm),k)*errwnds
         sumvv=sumvv+wt3(isave(mm),jsave(mm),k)*errwnds
         sumdd=sumdd+wt3(isave(mm),jsave(mm),k)*errdds              
         sumww=sumww+wt3(isave(mm),jsave(mm),k)
        enddo
c sum over ktop levels
        sumw=sumw+1.
        uav(k)=((sumuu+sqrt(upv(k)))/(sumww+1.))**2
        vav(k)=((sumvv+sqrt(vpv(k)))/(sumww+1.))**2
        tav(k)=uav(k)+vav(k)
        dav(k)=((sumdd+sqrt(dpv(k)))/(sumww+1.))**2
        sumu=sumu+uav(k)
        sumv=sumv+vav(k)
        sumd=sumd+dav(k)
c       write(6,1002) 'Anvar sen,u,v,tot,den ',index(l),uav(k),
c    &                 vav(k),tav(k),dav(k)

 1002   format(1x,a22,i6,4f9.5)
       enddo


       upi(l)=(sumu/sumw)
       vpi(l)=(sumv/sumw)
       twi(l)=upi(l)+vpi(l)
       dpi(l)=(sumd/sumw)
       write(6,1003)'scenario ',index(l),' uvar ',upi(l),' vvar ',vpi(l)
     & ,' tot  ',twi(l),' denvar ',dpi(l),' pct  ',twi(l)/tpvs*100. 
 1003  format(1x,a9,i6,3(a6,f6.3),a8,f7.5,a6,f6.3)
c determine variance minimum 
      call write_plan_errors(a9_time,p,ht,up,vp,wp,tp,dp,ud,vd,
     & ubb,vbb,wbb,uab,vab,wab,dbb,dab,vbtot,vbtot,ubv,vbv,wbv
     & ,uav,vav,wav,utv,vtv,wtv,dav,dbv,pav,htav,io,jo,nz,zter
     & ,index(l),99)
      enddo !on l 

      deallocate (wt,wt3)

c rank scenarios as a funcion of variance reduction
      do l=1,lmax
       inxx(l)=0 
      enddo
      do l=1,lmax
       smax=1.e30
       do m=1,lmax
        do n=1,l
         if(m.eq.inxx(n)) go to 34
        enddo
        if (twi(m).le.smax) then
         lsave =m
         smax=twi(m)
        endif
  34   enddo
       inxx(l)=lsave
      enddo
      do l=1,lmax
       m=inxx(l)
       write(6,1004) l, ' scenario ',index(m),' totl vari ',twi(m),
     &   ' pct reduction ', (1.-twi(m)/tpvs)*100.
 1004 format(1x,i2,a10,i6,a11,f9.5,a15,f7.3)
      enddo 

      return
      end

      subroutine decode_index(l,index,in,inx,lmax,ns)
c this routine takes in an integer index (l) and breaks it down into 
c individual members in array inx, counts the members..in
c example input 123, inx(l) is 3, 2, 1, and in=3   
      integer index(lmax),in,inx(lmax)
      do i100000=0,9 
       do i10000=0,9
        do i1000=0,9
         do i100=0,9
          do i10=0,9
           do i1=0,9
            iii=i100000*100000+i10000*10000+i1000*1000+
     &          i100*100+i10*10+i1
            if(index(l).eq.iii) then
             do ie=1,6
              if(iii.lt.10**ie) then
               in=ie
               go to 1
              endif
             enddo
    1        inx(1)=i1
             inx(2)=i10
             inx(3)=i100 
             inx(4)=i1000
             inx(5)=i10000
             inx(6)=i100000
             return
            endif
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo 
       print*, index(l), "error in index...sonde scenario not found" 
       return
       end


c
      subroutine write_plan_errors(a9_time,p,ht,up,vp,wp,tp,dp,ud,vd,
     & ubb,vbb,wbb,uab,vab,wab,dbb,dab,ubtot,vbtot,ubv,vbv,wbv
     & ,uav,vav,wav,utv,vtv,wtv,dav,dbv,pav,htav,io,jo,nz,zter,indx,alt)
c this routine writes an ascii file to the output log summarizing errors
      integer lun
      integer alt,li,lend,indx
      real p(nz),ht(nz)
       real up(nz),vp(nz),wp(nz),dp(nz),tp(nz)
       real ud(nz),vd(nz),ubb(nz),vbb(nz),wbb(nz),uab(nz),vab(nz),
     &      wab(nz),dbb(nz),dab(nz),ubv(nz),vbv(nz),wbv(nz),
     &      uav(nz),vav(nz),wav(nz),utv(nz),vtv(nz),wtv(nz),
     &      dav(nz),dbv(nz),pav(nz),htav(nz)  
      real rri(nz),rrj(nz)
      character*200 cfname_out
      character*9   a9_time
      character*4   cindx
c
c direct output to lapsprd/balance/pln/yyjjjhhmm.pln"scenario" -- ascii output
c
      call get_directory('balance',cfname_out,lend)

      if(alt.eq.0)then
         cfname_out=cfname_out(1:lend)//'air/'//a9_time//'.air'
c        call s_len(cfname_out,lend)
         print*,'Airdrop output filename: ',cfname_out(1:lend)
         lun = 20
         open(lun,file=cfname_out,form='formatted',status='unknown'
     +,err=909)

      elseif(alt.eq.1)then
         cfname_out=cfname_out(1:lend)//'air/'//a9_time//'.airp1'
c        call s_len(cfname_out,lend)
         print*,'Airdrop output filename: ',cfname_out(1:lend)
         lun = 20
         open(lun,file=cfname_out,form='formatted',status='unknown'
     +,err=909)

      elseif(alt.eq.2)then
         cfname_out=cfname_out(1:lend)//'air/'//a9_time//'.airp2'
c        call s_len(cfname_out,lend)
         print*,'Airdrop output filename: ',cfname_out(1:lend)
         lun = 20
         open(lun,file=cfname_out,form='formatted',status='unknown'
     +,err=909)

      elseif(alt.eq.99)then
         write(cindx,'(i4)')indx
         search_len: do li=4,1,-1
              if(cindx(li:li).eq.' ')then
                exit search_len
              endif
         enddo search_len
         li=li+1
         if(li.le.1)li=1
         call s_len(cfname_out,lend)
         cfname_out=cfname_out(1:lend)//'pln/'//a9_time//
     +'.pln'//cindx(li:4)
c        print*,'Airdrop planning output filename: ',cfname_out(1:lend)
         lun = 20
         open(lun,file=cfname_out,form='formatted',status='unknown'
     +,err=909)

      endif

  11  r=287.04
      sum=0
      sum1=0
      sum2=0

      if(alt.ne.99)then
         write(lun,1011)
      else
         write(lun,1050)
      endif
 1011 format(6x,'************AIRDROP OUTPUT**************')
 1050 format(6x,'*******AIRDROP PLANNING OUTPUT**********')
      if(alt.ne.99)then
         write(lun,1012)
      else
         write(lun,1052)
      endif 
 1012 format(6x,'AIRDROP PROFILES FROM LAPS ANALYSIS')
 1052 format(6x,'PLANNING PROFILES FROM LAPS ANALYSIS')
      write(lun,1013)
 1013 format(4x,'P',6x,'HT',6x,'U',8x,'V',8x,'W',8x,'T',7x,'Den')
      write(lun,1014)
 1014 format(3x,'mb',7x,'m',5x,'m/sec',4x,'m/sec',4x,'m/sec',6x,'K'
     +,6x,'kg/m3')
       ii=io     
       jj=jo    
      do k=1,nz
       if(zter.lt.ht(k)) then
        write(lun,1005) p(k)/100., ht(k),up(k),
     1  vp(k),wp(k),tp(k),dp(k)
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
     
       ii=io    
       jj=jo    
      do k=1,nz
       if(zter.lt.ht(k)) then
        write(lun,1025) p(k)/100., ht(k),ud(k),
     1  vd(k),ud(k),vd(k),tp(k),tp(k)          
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
       if(zter.lt.ht(k)) then
        write(lun,1086) p(k)/100., ht(k),ubb(k),
     1  vbb(k),wbb(k),
     1  uab(k),vab(k),wab(k),sum,sum,sum
       endif
      enddo
      write(lun,1015)
      write(lun,1016)
      write(lun,1017)
 1017 format(1x,'P LVL  HT LVL ')
      write(lun,1009)
      write(lun,1010)
      do k=1,nz
       if(zter.lt.ht(k)) then
       write(lun,1008) p(k)/100., ht(k),sum,sum,sum,sum
     
       endif
      enddo
   
      write(lun,1015)
      write(lun,1028)
 1028 format(10x,'TOTAL COLUMN WIND BIAS')
      write(lun,1029) ubtot,vbtot                  
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
       if(zter.lt.ht(k)      ) then
        write(lun,1000) p(k)/100., ht(k),ubv(k),
     1  vbv(k),wbv(k),
     1  uav(k),vav(k),wav(k),utv(k),vtv(k),wtv(k),dbv(k),dav(k)
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
       sum=uav(k)+utv(k)
       sum1=vav(k)+vtv(k)
       sum2=wav(k)+wtv(k)                            
       if(zter.lt.ht(k)) then
        write(lun,1001) p(k)/100., ht(k),       
     1 sum,sum1,sum2
     1,dav(k),pav(k),htav(k)                       
       endif
      enddo
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
C
C -----------
C
      subroutine read_PADS_nofly_and_sonde(nx,ny
     1,nofly_area,user_sonde_loc,istatus)

      implicit none

      character*255 directory
      character*255 cfilespec
      integer     nx,ny
      integer     i,j,idum,jdum
      integer     len_dir
      integer     len_cfspec
      integer     nofly_area(nx,ny)  !output
      integer     user_sonde_loc(2)  !output
      integer     istatus
      logical       L1

      istatus = 0
c first is no-fly-area file
      call get_directory('static',directory,len_dir)
      cfilespec=directory(1:len_dir)//'/'//'no-fly-area.txt'
      call s_len(cfilespec,len_cfspec)
      inquire(file=cfilespec(1:len_cfspec),EXIST=L1)
      if(.not.L1)then   !the file does not exist

         print*,"User file no-fly-area.txt does not exist"
         print*,"file spec: ",cfilespec(1:len_cfspec)
         nofly_area=1       !No restricted fly zones anywhere in domain

      else

         print*,"Reading user provided file no-fly-area.txt"
         open(10,FILE=cfilespec(1:len_cfspec),STATUS="old",ERR=90)
         read(10,*)        !    50   50   : nRows nCols (1,1 = SW corner)
         read(10,*)        !    Row   Col : 1=Fly, 0=No Fly
         do i=1,nx
          do j=1,ny
            read(10,*,err=91,end=92)idum,jdum,nofly_area(i,j)
          enddo
         enddo
         goto 95
92       print*,"!! Warning: End-of-file. ",cfilespec(1:len_cfspec)
95       close(10)

      endif
c second is user provided sonde location file
      cfilespec=directory(1:len_dir)//'/'//'sonde_location.txt'
      call s_len(cfilespec,len_cfspec)
      inquire(file=cfilespec(1:len_cfspec),EXIST=L1)
      user_sonde_loc(1)=-99
      user_sonde_loc(2)=-99
      if(.not.L1)then   !the file does not exist

         print*,"User file sonde_location.txt does not exist"
         print*,"file spec: ",cfilespec(1:len_cfspec)
         print*,"Return missing data (-99) for sonde location"

      else

         print*,"Reading user provided file sonde_location.txt"
         open(11,FILE=cfilespec(1:len_cfspec)
     1,STATUS="old",ERR=90)
         read(11,*,err=91,end=93)user_sonde_loc(1),user_sonde_loc(2)
         goto 96
93       print*,"!! Warning: End-of-file. ",cfilespec(1:len_cfspec)
96       close(11)

      endif
      istatus = 1
      return

 90   print*,"Error opening existing file: "
      print*,"  filename: ",cfilespec(1:len_cfspec)
      return
 91   print*,"Error reading file: ",cfilespec(1:len_cfspec)
      return
      end
