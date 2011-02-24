      subroutine tcbogus(jx,ix,nz,ht,tp,rh_c,uw,vw,
     +                   pss,tps,rhs,uws,vws,mslp,
     +                   lat1,lat2,lon0,sw,ne,dskm,glat,glon,
     +                   pressures_pa,filename,bgmodel,cwb_type)
!
!     Bogusing balanced tropical cyclones based on Rankin Vortex
!     Input file : tcbogus.nl
!======================================================================
!        field      unit          description
!         ht         m        geopotential height (upper)
!         tp         k        temperature (upper)
!         rh_c       %        relative humidity (upper)
!         uw        m/s       u-component wind (upper)
!         vw        m/s       v-component wind (upper)
!         tps        k        surface temperature
!         rhs        %        surface relative humidity
!         uws       m/s       surface u-component wind
!         vws       m/s       surface v-component wind
!         mslp      pa        mean sea level pressure
!      pressure_pa  pa        model vertical level
!       filename              background model full path
!        bgmodel              model number (support CWB/NFS bgmodel=3 now)
!         pss       hPa       surface pressure (not necessary!)
!
!      Guo-Ji Jian (CWB,Taiwan)       June, 2002
!      wen-ho Wang modify slightly     Sep, 2004
!      Steve Albers (more general inputs)   2010

      include 'trigd.inc'

      parameter(max_no=3)
      parameter(kb=7,k300=8,k850=4)

      character filename*255,nest7grid*150,cwb_type*8
      integer bgmodel,ix,jx,nz,itc_no,icount,cen_i,cen_j,l,
     +        y1,y2,y3,y4,m1,m2,d1,d2,a1,a2,b1,b2,istatus
      real ht(jx,ix,nz),tp(jx,ix,nz),rh_c(jx,ix,nz),uw(jx,ix,nz),
     +     vw(jx,ix,nz),pss(jx,ix),tps(jx,ix),rhs(jx,ix),uws(jx,ix),
     +     vws(jx,ix),mslp(jx,ix),pressures_pa(nz),ds
      real to(ix,jx,nz+1),uo(ix,jx,nz+1),vo(ix,jx,nz+1),ho(ix,jx,nz+1),
     +     rh(ix,jx,nz+1),xmapc(ix,jx),xmapd(ix,jx),pseald(ix,jx),
     +     psealc(ix,jx),f(ix,jx),p(nz+1),gas,g,mdl_hr
      real vort(ix,jx,nz+1),psi(ix,jx,nz+1),psio(ix,jx,nz+1),
     +     psi1(ix,jx,nz+1),psi2(ix,jx,nz+1),chi(ix,jx),rd(ix,jx),
     +     ff(ix,jx),outo(ix,jx),out1(ix,jx),fip(ix,jx,nz+1)
      real u1(ix,jx,nz+1),v1(ix,jx,nz+1),t1(ix,jx,nz+1),h1(ix,jx,nz+1),
     +     vorttc(ix,jx,nz+1),utcc(ix,jx,nz+1),vtcc(ix,jx,nz+1),
     +     wpsio(ix,jx),wpsi1(ix,jx),dd(ix,jx,nz+1),pd(ix,jx),
     +     pr(ix,jx),dot(ix,jx),hw(ix,jx,nz+1)
      real glat(jx,ix),glon(jx,ix),sw(2),ne(2),x,y,s,cone,xmin,ymin,
     +     dx,dy,lat1,lat2,lon0,dskm,init_hr,omg2,alpha_i(max_no),
     +     vmax(max_no),cal_lat(max_no),cal_lon(max_no)
      real speed(max_no),direc(max_no),cen_lat(max_no),
     +     cen_lon(max_no),rmax(max_no),r_7deg(max_no),alphax,vs,
     +     rs,rout,xix,xjx,dealbox,beta(kb),r,ang
      real year_m,month_m,day_m,hour_m,tothour_m,year,month,day,hour,
     +     tothour
      integer i4time_closest,i4time_sys,dt,i4time1,iread
      real dt1
      character*9 name1,a9
      logical bogus
      include 'grid_fname.cmn'
      data beta/0.90,1.00,0.90,0.80,0.60,0.45,0.35/
      NAMELIST /tcbogus_nl/bogus,itc_no,dskm,year,month,day,hour,
     +                     cen_lat,cen_lon,alpha_i,vmax,rmax,
     +                     r_7deg,speed,direc


! Bogusing namelist and driver
      call get_systime(i4time_sys,a9,istatus)
      call get_directory('tcbogus',nest7grid,len_dir)
      call get_file_time(nest7grid,i4time_sys,i4time_closest)
      call make_fnam_lp(i4time_closest,name1,istatus)
      nest7grid=nest7grid(1:len_dir)//name1//'_tcbogus.nl'
!      print*,'--- the table of tcbogus is following ---'
!      print *,nest7grid
      open(99,file=nest7grid,status='old',form='formatted',iostat=iread,
     + err=101)
      read(99,tcbogus_nl)
      close(99)
 101  if ( iread == 0 ) then 
       print*,'--- open  tcbogus file successful ---'
      else
       print*,'--- No   tcbogus file  ---'
      endif  

      
!  -- set i4time_sys be laps_time ---
!      i4time1=i4time_sys_gg()
      call get_systime(i4time_sys,a9,istatus)
      dt=i4time_sys - i4time_closest
      dt1=dt/3600.
        print*,'---check the file time:(sec) ---'
        print*,'-- dt,laps_time,latest ---',dt,i4time_sys,i4time_closest
!       if ( dt .ge. 10800 ) then 
       if ( dt .ge. 12000 ) then 
        print*,'---- the data is too old to use( about 3.3 hrs )---',dt1
           bogus=.false.
       else if ( dt .le. -2700) then  ! at most new 45 mins
        print*,'---- the data is too new to use(sec)---',dt
           bogus=.false.
       else 
        print*,'---- the data should be used (? hrs)---',dt1
           bogus=.true.
       end if
!
      if (bogus) then
       print *,' '
       print *,'Do Tropical cyclone bogusing'
      else
       print *,'No Tropical cyclone bogusing'
       return
      endif

      call s_len(filename,l)
       print *,' test cyclone bogusing',filename,l,jx,ix,nz

!      if (bgmodel.eq.3 .and. filename(l-15:l-14).eq.'nf') then
      if (bgmodel.eq.3 .and. cwb_type.eq.'nf') then
! OLD CWB NFS_15KM model
       print *,'----- Using OLD CWB ---- '
       lat1=10.
       lat2=40.
       lon0=120.
       sw(1)=15.80
       sw(2)=+109.24
       ne(1)=34.987
       ne(2)=+131.60
       dskm=15.
      elseif (bgmodel.eq.3 .and. cwb_type.eq.'nf15' .or. 
     +        cwb_type.eq.'gfs') then
! NEW CWB NFS_15KM & GFS_180 (interpolated into 15KM) model
       print *,'----- Using NEW CWB ---- '
       lat1=10.
       lat2=40.
       lon0=120.
       sw(1)=9.28194
       sw(2)=+109.7727
       ne(1)=35.26665
       ne(2)=+137.7342
       dskm=15.
      elseif (bgmodel.eq.3 .and. cwb_type.eq.'nf45') then
! NEW CWB NFS_45km  model
       print *,'----- Using NEW CWB 45 ---- '
       lat1=10.
       lat2=40.
       lon0=120.
       sw(1)=-5.34068
       sw(2)=+77.9186
       ne(1)=42.9281
       ne(2)=+180.2034
       dskm=45.
      elseif (bgmodel.eq.3 .and. cwb_type.eq.'tfs') then
! NEW CWB TFS_45KM  model
       print *,'----- Using NEW CWB TFS  ---- '
       lat1=10.
       lat2=40.
       lon0=120.
       sw(1)=-9.902
       sw(2)=+82.854
       ne(1)=52.219
       ne(2)=+199.610
       dskm=45.
      else
       print *,'Other model'
       return
      endif
!!
       print *,'----- Using cwb_type is ---- ',cwb_type

       call lc_param11(s,cone,xmin,ymin,dx,dy,jx,ix,
     +                 nz,lat1,lat2,lon0,sw,ne)

       do i=1,jx
        do j=1,ix
         x=(i-1)*dx+xmin
         y=(j-1)*dy+ymin
         glon(i,j)=lon0+atand(-s*x/y)/cone
         glat(i,j)=(90.-
     +        2.*atand((x/sind(cone*(glon(i,j)-lon0)))**(1./cone)))/s
        enddo
       enddo

       kx=nz+1
       p(1)=1001.
       do k=2,kx
        p(k)=pressures_pa(k-1)/100.
       enddo

       y1=ichar(filename(l-13:l-13))-48
       y2=ichar(filename(l-12:l-12))-48
       y3=ichar(filename(l-11:l-11))-48
       y4=ichar(filename(l-10:l-10))-48
       m1=ichar(filename(l-09:l-09))-48
       m2=ichar(filename(l-08:l-08))-48
       d1=ichar(filename(l-07:l-07))-48
       d2=ichar(filename(l-06:l-06))-48
       a1=ichar(filename(l-05:l-05))-48
       a2=ichar(filename(l-04:l-04))-48
       b1=ichar(filename(l-02:l-02))-48
       b2=ichar(filename(l-01:l-01))-48
       year_m=1000.*y1+100.*y2+10.*y3+y4
       month_m=10.*m1+m2
       day_m=10.*d1+d2
       hour_m=10.*(a1+b1)+(a2+b2)
 

! Tropical cyclone information
      call hourcalc(year,month,day,hour,tothour)
      call hourcalc(year_m,month_m,day_m,hour_m,tothour_m)
      do n=1,itc_no
       if (direc(n).le.90.) then
        cal_lat(n)=cen_lat(n)+(tothour_m-tothour)*
     +             speed(n)*cosd(direc(n))/111.
        cal_lon(n)=cen_lon(n)+(tothour_m-tothour)*
     +             speed(n)*sind(direc(n))/111.
       elseif (direc(n).gt.90. .and. direc(n).le.180.) then
        cal_lat(n)=cen_lat(n)-(tothour_m-tothour)*
     +             speed(n)*sind(direc(n)-90.)/111.
        cal_lon(n)=cen_lon(n)+(tothour_m-tothour)*
     +             speed(n)*cosd(direc(n)-90.)/111.
       elseif (direc(n).gt.180. .and. direc(n).le.270.) then
        cal_lat(n)=cen_lat(n)-(tothour_m-tothour)*
     +             speed(n)*cosd(direc(n)-180.)/111.
        cal_lon(n)=cen_lon(n)-(tothour_m-tothour)*
     +             speed(n)*sind(direc(n)-180.)/111.
       else
        cal_lat(n)=cen_lat(n)+(tothour_m-tothour)*
     +             speed(n)*sind(direc(n)-270.)/111.
        cal_lon(n)=cen_lon(n)-(tothour_m-tothour)*
     +             speed(n)*cosd(direc(n)-270.)/111.
       endif
      enddo

      gas=287.
      g=9.81
      ds=dskm*1000.
      omg2=2.*7.292/100000.

! Tropical cyclone numbers : itc_no (max_no=3)
      do icount=1,itc_no

       call cen_posi(cal_lat(icount),cal_lon(icount),
     +               glon,glat,cen_i,cen_j,jx,ix)
       if (cen_j.le.3 .or. cen_j.ge.jx-3 .or.
     +     cen_i.le.3 .or. cen_i.ge.ix-3) then
        print *,'Out or near model boundary. Ignore bogusing '
        return
       endif

       alphax=alpha_i(icount)
       vs=vmax(icount)*1.2
       rs=rmax(icount)
       rout=r_7deg(icount)
       dealbox=int(r_7deg(icount)/dskm)-1.0
       xix=cen_i+0.2
       xjx=cen_j+0.2
       print *, 'Tropical cyclone location (lon,lat) ',
     +           cal_lon(icount),cal_lat(icount)
       print *, '                          (  i,  j) ',
     +           cen_j,cen_i
       print *,vs,rout/dskm,dealbox

! Background model data
      do kk=1,kx
       do ii=1,ix
        do jj=1,jx
         if (kk.eq.1) then
          to(ii,jj,kk)=tps(jj,ii)
          uo(ii,jj,kk)=uws(jj,ii)
          vo(ii,jj,kk)=vws(jj,ii)
          rh(ii,jj,kk)=rhs(jj,ii)
         else
          to(ii,jj,kk)=tp(jj,ii,kk-1)
          uo(ii,jj,kk)=uw(jj,ii,kk-1)
          vo(ii,jj,kk)=vw(jj,ii,kk-1)
          ho(ii,jj,kk)=ht(jj,ii,kk-1)
          rh(ii,jj,kk)=rh_c(jj,ii,kk-1)
         endif
        enddo
       enddo
      enddo

      do ii=1,ix
       do jj=1,jx
        pss(ii,jj)=pss(ii,jj)
        psealc(ii,jj)=mslp(jj,ii)/100.
        pseald(ii,jj)=mslp(jj,ii)/100.
        f(ii,jj)=omg2*sind(glat(jj,ii))
        xmapc(ii,jj)=1.
        xmapd(ii,jj)=1.
       enddo
      enddo

! Rankin vortex profile
      do k=1,kb
       do i=1,ix
        do j=1,jx
         r=sqrt((i-xix)*(i-xix)+(j-xjx)*(j-xjx))*dskm
         if (r.gt.rs) then
          vr=vs*(rs/r)**(alphax)
         else
          vr=vs*(r/rs)
         endif
         ang=atan(abs((i-xix)/(j-xjx)))
         uux=beta(k)*vr*sin(ang)
         vvx=beta(k)*vr*cos(ang)
         if(float(i).gt.xix .and. float(j).lt.xjx) then
          utcc(i,j,k+1)=-uux
          vtcc(i,j,k+1)=-vvx
         elseif(float(i).lt.xix .and. float(j).lt.xjx) then
          utcc(i,j,k+1)=+uux
          vtcc(i,j,k+1)=-vvx
         elseif(float(i).lt.xix .and. float(j).gt.xjx) then
          utcc(i,j,k+1)=+uux
          vtcc(i,j,k+1)=+vvx
         else
          utcc(i,j,k+1)=-uux
          vtcc(i,j,k+1)=+vvx
         endif
        enddo
       enddo
      enddo


! dew point temp. depression
      do k=1,kx
       do i=1,ix
        do j=1,jx
         rh(i,j,k)=amax1(rh(i,j,k),5.)
         dd(i,j,k)=to(i,j,k)-1./(1./to(i,j,k)-
     +             (1./5418.12)*alog(rh(i,j,k)/100.0))
        enddo
       enddo
      enddo

! Replace u,v,t,h inside the box
      ix1=int(xix-0.5)-int(dealbox)
      ix2=int(xix+0.5)+int(dealbox)
      jx1=int(xjx-0.5)-int(dealbox)
      jx2=int(xjx+0.5)+int(dealbox)

      call vor(uo,vo,xmapd,xmapc,ix,jx,kx,ds,vort)
      do k=2,k300+1
       do i=1,ix
        do j=1,jx
         ff(i,j)=vort(i,j,k)
         chi(i,j)=0.
        enddo
       enddo
       call relax(chi,ff,rd,ix,jx,ds,p(k),istatus)
       if (istatus.ne.1) then
        print *,'bogusing failed !'
        return
       endif

       do i=1,ix
        do j=1,jx
         psio(i,j,k)=chi(i,j)
        enddo
       enddo
      enddo

      do k=2,k300+1
       do i=ix1,ix2
        do j=jx1,jx2
         r=sqrt((i-xix)*(i-xix)+(j-xjx)*(j-xjx))*dskm
         if (r.le.rout) vort(i,j,k)=vort(i,j,k)*0.
        enddo
       enddo
      enddo

      do k=2,k300+1
       do i=1,ix
        do j=1,jx
         ff(i,j)=vort(i,j,k)
         chi(i,j)=0.
        enddo
       enddo
       call relax(chi,ff,rd,ix,jx,ds,p(k),istatus)
       if (istatus.ne.1) then
        print *,'bogusing failed !'
        return
       endif

       do i=1,ix
        do j=1,jx
         psi1(i,j,k)=chi(i,j)
        enddo
       enddo
      enddo

! Calculate perturbation

      do k=2,k300+1
       do i=1,ix
        do j=1,jx
         psi(i,j,k)=psio(i,j,k)-psi1(i,j,k)
        enddo
       enddo
      enddo

      do k=2,k300+1
       do i=2,ix-1
        do j=2,jx-1
         upp=-((psi(i  ,j  ,k)+psi(i  ,j-1,k))-
     +         (psi(i-1,j-1,k)+psi(i-1,j  ,k)))/(2.*ds)
         vpp=+((psi(i  ,j  ,k)+psi(i-1,j  ,k))-
     +         (psi(i-1,j-1,k)+psi(i  ,j-1,k)))/(2.*ds)
         u1(i,j,k)=uo(i,j,k)-upp
         v1(i,j,k)=vo(i,j,k)-vpp
        enddo
       enddo
      enddo
      call fillit(u1,ix,jx,kx,ix,jx,2,ix-1,2,jx-1)
      call fillit(v1,ix,jx,kx,ix,jx,2,ix-1,2,jx-1)

      do k=2,k300+1
       do i=1,ix
        do j=1,jx
         wpsio(i,j)=psio(i,j,k)
         wpsi1(i,j)=psi1(i,j,k)
        enddo
       enddo
       call balance(f,wpsio,ix,jx,ds,outo)
       call balance(f,wpsi1,ix,jx,ds,out1)
       do i=1,ix
        do j=1,jx
         ff(i,j)=outo(i,j)-out1(i,j)
         chi(i,j)=0.
        enddo
       enddo
       call relax(chi,ff,rd,ix,jx,ds,p(k),istatus)
       if (istatus.ne.1) then
        print *,'bogusing failed !'
        return
       endif

       call crs2dot(chi,dot,ix,jx,ix,jx)

       do i=1,ix
        do j=1,jx
         fip(i,j,k)=chi(i,j)
         hw(i,j,k)=dot(i,j)/G
        enddo
       enddo
      enddo

! Find lower level minimum
      hminc=0.
      hmind=0.
      do i=1,ix
       do j=1,jx
        if (fip(i,j,2).le.hminc) then
         hminc=fip(i,j,2)
         imin1=i
         jmin1=j
        endif
        if (hw(i,j,2).le.hmind) then
         hmind=hw(i,j,2)*g
         imin2=i
         jmin2=j
        endif
       enddo
      enddo

! Temperature perturbation and background temperature
      do i=1,ix
       do j=1,jx
        do k=2,2
         fip(i,j,k)=fip(i,j,k850)
         hw(i,j,k)=hw(i,j,k850)
        enddo
       enddo
      enddo

      do k=2,k300+1
       do i=1,ix
        do j=1,jx
         h1(i,j,k)=ho(i,j,k)-hw(i,j,k)
        enddo
       enddo
      enddo

      do k=2,k300+1
       do i=1,ix
        do j=1,jx
         if (k.eq.2) then
          t1(i,j,k)=to(i,j,k)-(-1./gas)*(fip(i,j,k+1)-
     +              fip(i,j,k))/log(p(k+1)/p(k))
         elseif (k.eq.k300+1) then
          t1(i,j,k)=to(i,j,k)-(-1./gas)*(fip(i,j,k)-
     +              fip(i,j,k-1))/log(p(k)/p(k-1))
         else
          t1(i,j,k)=to(i,j,k)-(-1./gas)*(fip(i,j,k+1)-
     +              fip(i,j,k-1))/log(p(k+1)/p(k-1))
         endif
        enddo
       enddo
      enddo

! Calculate vorticity with u1 v1
      call vor(u1,v1,xmapd,xmapc,ix,jx,kx,ds,vort)
      call vor(utcc,vtcc,xmapd,xmapc,ix,jx,kx,ds,vorttc)

! Replace vorticity using Rankin vortex vorticity
      do k=2,k300+1
       do i=ix1,ix2
        do j=jx1,jx2
         r=sqrt((i-xix)*(i-xix)+(j-xjx)*(j-xjx))*dskm
         if (r.le.rout) vort(i,j,k)=vorttc(i,j,k)
        enddo
       enddo
      enddo

      do k=2,k300+1
       do i=1,ix
        do j=1,jx
         ff(i,j)=vort(i,j,k)
         chi(i,j)=0.
        enddo
       enddo
       call relax(chi,ff,rd,ix,jx,ds,p(k),istatus)
       if (istatus.ne.1) then
        print *,'bogusing failed !'
        return
       endif

       do i=1,ix
        do j=1,jx
         psi2(i,j,k)=chi(i,j)
        enddo
       enddo
      enddo

      do k=2,k300+1
       do i=1,ix
        do j=1,jx
         psi(i,j,k)=psi2(i,j,k)-psi1(i,j,k)
        enddo
       enddo
      enddo

      do k=2,k300+1
       do i=2,ix-1
        do j=2,jx-1
         upp=-((psi(i  ,j  ,k)+psi(i  ,j-1,k))-
     +         (psi(i-1,j-1,k)+psi(i-1,j  ,k)))/(2.*ds)
         vpp=+((psi(i  ,j  ,k)+psi(i-1,j  ,k))-
     +         (psi(i-1,j-1,k)+psi(i  ,j-1,k)))/(2.*ds)
         uo(i,j,k)=u1(i,j,k)+upp
         vo(i,j,k)=v1(i,j,k)+vpp
        enddo
       enddo
      enddo
      
      do k=2,k300+1
       do i=1,ix
        do j=1,jx
         wpsio(i,j)=psi2(i,j,k)
         wpsi1(i,j)=psi1(i,j,k)
        enddo
       enddo
       call balance(f,wpsio,ix,jx,ds,outo)
       call balance(f,wpsi1,ix,jx,ds,out1)
       do i=1,ix
        do j=1,jx
         ff(i,j)=outo(i,j)-out1(i,j)
         chi(i,j)=0.
        enddo
       enddo
       call relax(chi,ff,rd,ix,jx,ds,p(k),istatus)
       if (istatus.ne.1) then
        print *,'bogusing failed !'
        return
       endif

       call crs2dot(chi,dot,ix,jx,ix,jx)

       do i=1,ix
        do j=1,jx
         fip(i,j,k)=chi(i,j)
         hw(i,j,k)=dot(i,j)/g
        enddo
       enddo
      enddo

! Find lower level minimum
      hminc=0.
      hmind=0.
      do i=1,ix
       do j=1,jx
        if (fip(i,j,2).le.hminc) then
         hminc=fip(i,j,2)
         imin1=i
         jmin1=j
        endif
        if (hw(i,j,2).le.hmind) then
         hmind=hw(i,j,2)*g
         imin2=i
         jmin2=j
        endif
       enddo
      enddo

! Calculate new Temperature 
      do i=1,ix
       do j=1,jx
        do k=2,2
         fip(i,j,k)=fip(i,j,k850)
         hw(i,j,k)=hw(i,j,k850)
        enddo
       enddo
      enddo

      do k=2,k300+1
       do i=1,ix
        do j=1,jx
         ho(i,j,k)=h1(i,j,k)+hw(i,j,k)
        enddo
       enddo
      enddo

      do k=2,k300+1
       do i=1,ix
        do j=1,jx
         if (k.eq.2) then
          to(i,j,k)=t1(i,j,k)+(-1./gas)*(fip(i,j,k+1)-
     +              fip(i,j,k))/log(p(k+1)/p(k))
         elseif (k.eq.k300+1) then
          to(i,j,k)=t1(i,j,k)+(-1./gas)*(fip(i,j,k)-
     +              fip(i,j,k-1))/log(p(k)/p(k-1))
         else
          to(i,j,k)=t1(i,j,k)+(-1./gas)*(fip(i,j,k+1)-
     +              fip(i,j,k-1))/log(p(k+1)/p(k-1))
         endif
        enddo
       enddo
      enddo

      imaxt=1
      jmaxt=1
      tmax=to(1,1,k850)
      do i=ix1,ix2
       do j=jx1,jx2
        if (to(i,j,k850).gt.tmax) then
         tmax=to(i,j,k850)
         imaxt=i
         jmaxt=j
        endif
       enddo
      enddo

! RH

      do k=2,k300+1
       do i=1,ix-1
        do j=1,jx-1
         rh(i,j,k)=100.*exp(5418.12*(1./to(I,J,K)-
     +                      1./(to(I,J,K)-dd(I,J,K))))
        enddo
       enddo
      enddo
 
! Calculate slp
      poo=1000.
      density=1.225
      do i=1,ix
       do j=1,jx
        pr(i,j)=density*fip(i,j,k850)/100.
       enddo
      enddo

      do i=1,ix-1
       do j=1,jx-1
        pd(i,j)=psealc(i,j)-poo
       enddo
      enddo

      call delta(pr,ix,jx,ds,out1)
      call delta(pd,ix,jx,ds,outo)
      do i=1,ix
       do j=1,jx
        ff(i,j)=outo(i,j)
        chi(i,j)=pd(i,j)
       enddo
      enddo

      do i=ix1,ix2
       do j=jx1,jx2
        r=sqrt((i-xix)*(i-xix)+(j-xjx)*(j-xjx))*dskm
        if (r.le.rout) ff(i,j)=outo(i,j)*0.
       enddo
      enddo

      do i=ix1,ix2
       do j=jx1,jx2
        r=sqrt((i-xix)*(i-xix)+(j-xjx)*(j-xjx))*dskm
        if (r.le.rout) then
         ff(i,j)=out1(i,j)
         chi(i,j)=pr(i,j)
        endif
       enddo
      enddo
      call relax(chi,ff,rd,ix,jx,ds,p(1),istatus)
      if (istatus.ne.1) then
       print *,'bogusing failed !'
       return
      endif

      do i=1,ix
       do j=1,jx
        psealc(i,j)=chi(i,j)+poo
       enddo
      enddo
      call crs2dot(psealc,dot,ix,jx,ix,jx)
      do i=1,ix
       do j=1,jx
        pseald(i,j)=dot(i,j)
       enddo
      enddo

! Where is the min slp of bogusing tropical cyclone ?
      pc1=psealc(ix1,jx1)
      pc2=pseald(ix1,jx1)
      iminpc=1
      jminpc=1
      iminpd=1
      jminpd=1
      do i=ix1,ix2
       do j=jx1,jx2
        if (psealc(i,j).lt.pc1) then
         pc1=psealc(i,j)
         iminpc=i
         jminpc=j
        endif 
        if (pseald(i,j).lt.pc2) then
         pc2=pseald(i,j)
         iminpd=i
         jminpd=j
        endif
       enddo
      enddo

! Surface wind
      do i=1,ix
       do j=1,jx
        utcc(i,j,1)=uo(i,j,2)
        vtcc(i,j,1)=vo(i,j,2)
       enddo
      enddo
      call vor(uo,vo,xmapd,xmapc,ix,jx,kx,ds,vort)

      do i=1,ix
       do j=1,jx
        ff(i,j)=vort(i,j,1)
        chi(i,j)=0.
       enddo
      enddo
      call relax(chi,ff,rd,ix,jx,ds,p(1),istatus)
      if (istatus.ne.1) then
       print *,'bogusing failed !'
       return
      endif

      do i=1,ix
       do j=1,jx
        psio(i,j,1)=chi(i,j)
       enddo
      enddo

      do i=ix1,ix2
       do j=jx1,jx2
        r=sqrt((i-xix)*(i-xix)+(j-xjx)*(j-xjx))*dskm
        if (r.le.rout) vort(i,j,1)=0.
       enddo
      enddo

      do i=1,ix
       do j=1,jx
        ff(i,j)=vort(i,j,1)
        chi(i,j)=0.
       enddo
      enddo
      call relax(chi,ff,rd,ix,jx,ds,p(1),istatus)
      if (istatus.ne.1) then
       print *,'bogusing failed !'
       return
      endif

      do i=1,ix
       do j=1,jx
        psi1(i,j,1)=chi(i,j)
       enddo
      enddo

      do i=1,ix
       do j=1,jx
        psi(i,j,1)=psio(i,j,1)-psi1(i,j,1)
       enddo
      enddo

      do i=2,ix-1
       do j=2,jx-1
        upp=-((psi(i  ,j  ,1)+psi(i  ,j-1,1))-
     +        (psi(i-1,j-1,1)+psi(i-1,j  ,1)))/(2.*ds)
        vpp=+((psi(i  ,j  ,1)+psi(i-1,j  ,1))-
     +        (psi(i-1,j-1,1)+psi(i  ,j-1,1)))/(2.*ds)
        u1(i,j,1)=uo(i,j,1)-upp
        v1(i,j,1)=vo(i,j,1)-vpp
       enddo
      enddo
      call fillit(u1,ix,jx,kx,ix,jx,2,ix-1,2,jx-1)
      call fillit(v1,ix,jx,kx,ix,jx,2,ix-1,2,jx-1)

      call vor(u1,v1,xmapd,xmapc,ix,jx,kx,ds,vort)
      call vor(utcc,vtcc,xmapd,xmapc,ix,jx,kx,ds,vorttc)

      do i=ix1,ix2
       do j=jx1,jx2
        r=sqrt((i-xix)*(i-xix)+(j-xjx)*(j-xjx))*dskm
        if (r.le.rout) vort(i,j,1)=vorttc(i,j,1)
       enddo
      enddo

      do i=1,ix
       do j=1,jx
        ff(i,j)=vort(i,j,1)
        chi(i,j)=0.
       enddo
      enddo
      call relax(chi,ff,rd,ix,jx,ds,p(1),istatus)
      if (istatus.ne.1) then
       print *,'bogusing failed !'
       return
      endif

      do i=1,ix
       do j=1,jx
        psi2(i,j,1)=chi(i,j)
       enddo
      enddo

      do i=1,ix
       do j=1,jx
        psi(i,j,1)=psi2(i,j,1)-psi1(i,j,1)
       enddo
      enddo

      do i=2,ix-1
       do j=2,jx-1
        upp=-((psi(i  ,j  ,1)+psi(i  ,j-1,1))-
     +        (psi(i-1,j-1,1)+psi(i-1,j  ,1)))/(2.*ds)
        vpp=+((psi(i  ,j  ,1)+psi(i-1,j  ,1))-
     +        (psi(i-1,j-1,1)+psi(i  ,j-1,1)))/(2.*ds)
        uo(i,j,1)=u1(i,j,1)+upp
        vo(i,j,1)=v1(i,j,1)+vpp
       enddo
      enddo

! Find the max. wind
      vmax0=uo(cen_i,cen_j,1)
      vmax1=uo(cen_i,cen_j,2)
      do i=ix1,ix2
       do j=jx1,jx2
        v00=sqrt(uo(i,j,1)*uo(i,j,1)+vo(i,j,1)*vo(i,j,1))
        v11=sqrt(uo(i,j,2)*uo(i,j,2)+vo(i,j,2)*vo(i,j,2))
        if (v00.gt.vmax0) then
         vmax0=v00
         imax0=i
         jmax0=j
        endif
        if (v11.gt.vmax1) then
         vmax1=v11
         imax1=i
         jmax1=j
        endif
       enddo
      enddo

      print *,' '
      print *,'Tropical cyclone bogusing result : '
      print 2000,pc1,glat(jminpc,iminpc),glon(jminpc,iminpc)
      print 2001,vmax0,glat(jmax0,imax0),glon(jmax0,imax0)
      print 2002,vmax1,glat(jmax1,imax1),glon(jmax1,imax1)
      print 2003,tmax,glat(jmaxt,imaxt),glon(jmaxt,imaxt)
      print *,' '
2000  format(' Min sea level pressure    : ',f7.2,' hPa  at',
     +        f7.2,'N',f7.2,'E')
2001  format(' Max wind speed at surface : ',f7.2,' m/s  at',
     +        f7.2,'N',f7.2,'E')
2002  format(' Max wind speed at 1000hPa : ',f7.2,' m/s  at',
     +        f7.2,'N',f7.2,'E')
2003  format(' 850hPa max temperaturea   : ',f7.2,' K    at',
     +        f7.2,'N',f7.2,'E')

! Surface RH
      do i=1,ix-1
       do j=1,jx-1
        rh(i,j,1)=100.*exp(5418.12*(1./to(I,J,1)-
     +                     1./(to(I,J,1)-dd(I,J,1))))
       enddo
      enddo

!**********************penny**************************
      rovcp=287./1004.
      print *, 'Penny test!'
      do i=1,ix-1
       do j=1,jx-1
        ps=psealc(i,j)
!        to(i,j,1)=to(i,j,2)*((ps/1000.)**rovcp)
        to(i,j,1)=to(i,j,3)*((ps/925.)**rovcp)
        if (to(i,j,1) .lt. 283.15) then
        print *, to(i,j,1)-273.15,ps,i,j,'sfct'
        endif
        to(i,j,2)=to(i,j,3)*((1000./925.)**rovcp)
        if (to(i,j,2) .lt. 283.15) then
        print *, to(i,j,2)-273.15,p1000,i,j,'1000t'
        endif
       enddo
      enddo
********************************************************

! After bogusing
       do kk=1,kx
        do ii=1,ix
         do jj=1,jx
          if (kk.eq.1) then
           tps(jj,ii)=to(ii,jj,kk)
           uws(jj,ii)=uo(ii,jj,kk)
           vws(jj,ii)=vo(ii,jj,kk)
           rhs(jj,ii)=rh(ii,jj,kk)
          else
           tp(jj,ii,kk-1)=to(ii,jj,kk)
           uw(jj,ii,kk-1)=uo(ii,jj,kk)
           vw(jj,ii,kk-1)=vo(ii,jj,kk)
           ht(jj,ii,kk-1)=ho(ii,jj,kk)
           rh_c(jj,ii,kk-1)=rh(ii,jj,kk)
          endif
         enddo
        enddo
       enddo

       do ii=1,ix
        do jj=1,jx
         pss(ii,jj)=pss(ii,jj)
         mslp(jj,ii)=psealc(ii,jj)*100.
        enddo
       enddo

      enddo
      return
      end

      subroutine lc_param11(s,cone,xmin,ymin,dx,dy,
     *                    nx,ny,nz,lat1,lat2,lon0,sw,ne)

      include 'trigd.inc'

      real s,cone,r,
     .       xmin,xmax,ymin,ymax,
     .       dx,dy
c
      real lat1,lat2,lon0,       !Lambert-conformal std lat1, lat2, lon
     .       sw(2),ne(2)           !SW lat, lon, NE lat, lon
      integer nx,ny,nz           !No. of LC domain grid points
c_______________________________________________________________________________
c
      if (lat1 .ge. 0.) then
         s=1.
      else
         s=-1.
      endif
      if (lat1 .ne. lat2) then
         cone=alog(cosd(lat1)/cosd(lat2))/
     .        alog(tand(45.-s*lat1/2.)/tand(45.-s*lat2/2.))
      else
         cone=cosd(90.-s*lat1)
      endif
c
      r=(tand(45.-s*sw(1)/2.))**cone
      xmin=r*sind(cone*(sw(2)-lon0))
      ymin=-s*r*cosd(cone*(sw(2)-lon0))
      r=(tand(45.-s*ne(1)/2.))**cone
      xmax=r*sind(cone*(ne(2)-lon0))
      ymax=-s*r*cosd(cone*(ne(2)-lon0))
      dx=(xmax-xmin)/float(nx-1)
      dy=(ymax-ymin)/float(ny-1)
c
      return
      end

      subroutine hourcalc(year,month,day,hour,tothours)

      parameter(iydim=1000,imdim=12,begin_year=1950)
      real year,month,day,hour,tothours,yrdays,mndays
      real dayarray(iydim,imdim),yrmark(iydim)

      yrdays=0.
      mndays=0.
      tothours=0.
      do 20 i=1,iydim
20     yrmark(i)=i+begin_year-1

      do i=1,iydim
       do j=1,imdim
        if ((j.eq.1).or.(j.eq.3).or.(j.eq.5).or.(j.eq.7).or.
     +      (j.eq.8).or.(j.eq.10).or.(j.eq.12)) then
            dayarray(i,j)=31.
        elseif ((j.eq.4).or.(j.eq.6).or.(j.eq.9).or.(j.eq.11)) then
            dayarray(i,j)=30.
        elseif (j.eq.2) then
         if (mod(int(yrmark(i)),4).eq.0) then
            dayarray(i,j)=29.
         else
            dayarray(i,j)=28.
         endif
        endif
       enddo
      enddo

      do i=1,int(year-begin_year)
       do j=1,imdim
        yrdays=yrdays+dayarray(i,j)
       enddo
      enddo

      do j=1,int(month-1)
       mndays=mndays+dayarray(int(year-begin_year+1),j)
      enddo

      tothours=(yrdays+mndays+(day-1))*24+hour
      return
      end

      subroutine cen_posi(lat,lon,glon,glat,cen_I,cen_J,JX,IX)
      REAL lat,lon,glon(JX,IX),glat(JX,IX),dis2,min
      INTEGER cen_I,cen_J

      min=1.E+20
      do i=1,JX
       do j=1,IX
        dis2=(lat-glat(i,j))**2+(lon-glon(i,j))**2
        if (dis2.lt.min) then
          cen_I=j
          cen_J=i
          min=dis2
        endif
       enddo
      enddo

      return
      end

      subroutine vor(u,v,dmf,xmf,ix,jx,kx,ds,vort)
      integer ix,jx,kx
      real u(ix,jx,kx),v(ix,jx,kx),vort(ix,jx,kx)
      real dmf(ix,jx),xmf(ix,jx),ds,u1,u2,u3,u4,
     +     v1,v2,v3,v4

      do k=1,kx
       do j=1,jx-1
        do i=1,ix-1
         u1=u(i,j,k)/dmf(i,j)
         u2=u(i+1,j,k)/dmf(i+1,j)
         u3=u(i,j+1,k)/dmf(i,j+1)
         u4=u(i+1,j+1,k)/dmf(i+1,j+1)

         v1=v(i,j,k)/dmf(i,j)
         v2=v(i+1,j,k)/dmf(i+1,j)
         v3=v(i,j+1,k)/dmf(i,j+1)
         v4=v(i+1,j+1,k)/dmf(i+1,j+1)

         vort(i,j,k)=((v4-v2+v3-v1)-(u2-u1+u4-u3))
     +               /(2.*ds)*(xmf(i,j)*xmf(i,j))
        enddo
       enddo
      enddo
      call fillit(vort,ix,jx,kx,ix,jx,1,ix-1,1,jx-1)

      return
      end

      subroutine relax(chi,ff,rd,ix,jx,ds,level,istatus)
      integer ix,jx,mm,icou,istatus
      real chi(ix,jx),ff(ix,jx),rd(ix,jx),fac,
     +     chimx,rdmax,eeps,epx,alpha,level

      istatus=0
      eeps=1.E-6
      alpha=0.45
      mm=5000
      fac=ds*ds*1.

      do i=1,ix
       do j=1,jx
        rd(i,j)=0.
       enddo
      enddo

      do j=1,jx-1
       do i=1,ix-1
        ff(i,j)=ff(i,j)*fac
       enddo
      enddo

      do icou=1,mm
       chimx=0.
       do j=2,jx-2
        do i=2,ix-2
         chimx=amax1(abs(chi(i,j)),chimx)
        enddo
       enddo

       epx=chimx*eeps/alpha
       do j=2,jx-2
        do i=2,ix-2
         rd(i,j)=chi(i,j+1)+chi(i,j-1)+chi(i+1,j)+
     +           chi(i-1,j)-4.*chi(i,j)-ff(i,j)
         chi(i,j)=chi(i,j)+rd(i,j)*alpha
        enddo
       enddo

       rdmax=0.
       do j=2,jx-2
        do i=2,ix-2
         rdmax=amax1(abs(rd(i,j)),rdmax)
        enddo
       enddo
       if (rdmax.le.epx) then
        if (level.eq.1001.) then
         print 2004,icou
        else
         print 2005,level,icou
        endif
        istatus=1
        call fillit(chi,ix,jx,1,ix,jx,2,ix-2,2,jx-2)
        return
       endif
      enddo

      print *,'Error ! relaxation not converge.  bogusing abort ! '
2004  format(' Surface  done, iterations = ',i4)
2005  format(1x,f5.0,' hPa  done, iterations = ',i4)
      return
      end

      subroutine fillit(f,ix,jx,kx,imx,jmx,ifirst,ilast,jfirst,jlast)
      integer ix,jx,kx,imx,jmx,ifirst,ilast,jfirst,jlast
      real f(ix,jx,kx)

      do k=1,kx
       do j=jfirst,jlast
        do i=1,ifirst-1
         f(i,j,k)=f(ifirst,j,k)
        enddo
        do i=ilast+1,imx
         f(i,j,k)=f(ilast,j,k)
        enddo
       enddo

       do i=1,imx
        do j=1,jfirst-1
         f(i,j,k)=f(i,jfirst,k)
        enddo
        do j=jlast+1,jmx
         f(i,j,k)=f(i,jlast,k)
        enddo
       enddo
      enddo

      return
      end

      subroutine balance(f,psi,ix,jx,ds,out)
      integer ix,jx
      real f(ix,jx),psi(ix,jx),out(ix,jx),ds,
     +     psixx,psiyy,psiy,psixy

      do i=2,ix-2
       do j=2,jx-2
        psixx=(psi(i,j+1)+psi(i,j-1)-2.*psi(i,j))/(ds*ds)
        psiyy=(psi(i+1,j)+psi(i-1,j)-2.*psi(i,j))/(ds*ds)
        psiy =(psi(i+1,j)-psi(i-1,j))/(2.*ds)
        psixy=(psi(i+1,j+1)+psi(i-1,j-1)-
     +         psi(i-1,j+1)-psi(i+1,j-1))/(4.*ds*ds)
        out(i,j)=0.25*(f(i,j)+f(i+1,j)+f(i+1,j+1)+
     +                 f(i,j+1))*(psixx+psiyy)+psiy*
     +                (f(i+1,j+1)+f(i+1,j)-f(i,j)-f(i,j+1))/
     +                (2.*ds)-2.*(psixy*psixy-psixx*psiyy)
       enddo
      enddo
      call fillit(out,ix,jx,1,ix,jx,2,ix-2,2,jx-2)

      return
      end

      subroutine delta(psi,ix,jx,ds,out)
      integer ix,jx
      real psi(ix,jx),out(ix,jx),ds
      real psixx,psiyy

      do i=2,ix-2
       do j=2,jx-2
        psixx=(psi(i,j+1)+psi(i,j-1)-2.*psi(i,j))/(ds*ds)
        psiyy=(psi(i+1,j)+psi(i-1,j)-2.*psi(i,j))/(ds*ds)
        out(i,j)=psixx+psiyy
       enddo
      enddo
      call fillit(out,ix,jx,1,ix,jx,2,ix-2,2,jx-2)

      return
      end

      subroutine crs2dot(crs,dot,ix,jx,iend,jend)
      integer ix,jx,iend,jend
      real crs(ix,jx),dot(ix,jx)

      do j=2,jend-1
       do i=2,iend-1
        dot(i,j)=(crs(i,j)+crs(i-1,j)+crs(i,j-1)+crs(i-1,j-1))/4.
       enddo
      enddo

      do j=2,jend-1
       dot(1,j)=((crs(1,j)+crs(1,j-1))*1.5-
     +           (crs(2,j)+crs(2,j-1))*0.5)/2.
       dot(iend,j)=((crs(iend-1,j)+crs(iend-1,j-1))*1.5-
     +              (crs(iend-2,j)+crs(iend-2,j-1))*0.5)/2.
      enddo

      do i=2,iend-1
       dot(i,1)=((crs(i,1)+crs(i-1,1))*1.5-
     +           (crs(i,2)+crs(i-1,2))*0.5)/2.
       dot(i,jend)=((crs(i,jend-1)+crs(i-1,jend-1))*1.5-
     +              (crs(i,jend-2)+crs(i-1,jend-2))*0.5)/2.
      enddo

      dot(1,1)=4.*crs(1,1)-(dot(1,2)+dot(2,1)+dot(2,2))
      dot(iend,1)=4.*crs(iend-1,1)-(dot(iend-1,1)+
     +               dot(iend-1,2)+dot(iend,2))
      dot(1,jend)=4.*crs(1,jend-1)-(dot(1,jend-1)+
     +               dot(2,jend-1)+dot(2,jend))
      dot(iend,jend)=4.*crs(iend-1,jend-1)-(dot(iend,jend-1)+
     +                  dot(iend-1,jend-1)+dot(iend-1,jend))

      return
      end
