c
c
        Subroutine model(dta,var,ua,va,theta,mwt,lat,lon,elev,
     &                   imax,m,mm,time,cycle,amp,amp1,peak,peak1)
c
c*********************************************************************
c
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       24 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c       14 Dec 1999  John McGinley and Peter Stamus, NOAA/FSL
c          Completely new version.
c
c     Notes:
c          var is the variable for advective term
c
c*********************************************************************
c
        parameter (badflag=-99.9)
        common tab(10000),pi,re,rdpdg,reorpd
        real dta(m),mwt(m,m),lat(m),lon(m),elev(m)
        real advect(mm)  ,var(m) 
        real kappa,ua(m),va(m),theta(m) 
        integer time(m), time_i
c
c.....  Start here.
c
        kappa=2./7.
        stdlps=-6.5/1000.
        Tstd=287.16
        pstd=1013.2
        g=9.808
        time_i = time(1)
c
c.....  Here would be reading model grids and interpolating to stns
c.....  set offset for heating: max temp 00Z; min temp 12Z
c
        ddiur=-amp*(2.*pi/(24.*cycle)) *         !correcting waves
     &         sin(2.*pi/24.*(time_i-peak))*cycle
     &        -amp1*2.*pi/12.*sin(2.*pi/12.*float(time_i)-peak1)
c
        do i=1,imax
           sum=0.
           cnt=0.    
           do j=1,imax
              if(i.eq.j) go to 1 
              ang=(lat(i)+lat(j))*.5*rdpdg
              dth=var(j)-var(i)
              dx=re*(lon(j)-lon(i))*cos(ang)*rdpdg
              dy=re*(lat(j)-lat(i))*rdpdg
              rr=(dx*dx+dy*dy) 
              if(ua(i).ne.badflag.and.ua(j).ne.badflag) then
                 uav=(ua(i)+ua(j))*.5
              else
                 uav=0.
              endif
              if(va(i).ne.badflag.and.va(j).ne.badflag) then
                 vav=(va(i)+va(j))*.5
              else
                 vav=0.
              endif
c
c corrects for wind component along station/station vector
              sum=-(uav*dth*dx/rr+vav*dth*dy/rr)*cycle+sum
              cnt=cnt+1.
 1         enddo !on j
           advect(i)=sum/cnt
        enddo !on i
c
c.....  mwt is a model spatial correlation to input what the expected
c.....  spatial error is likely to be with the model
c
        do i=1,imax
           dta(i)=ddiur+advect(i)
        enddo !on i
c
        return
        end
c
c
	Subroutine perturb(ta,tb,dta,imax,m,offset,onoff)     
c       
c*********************************************************************
c
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       24 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c
c     Notes:
c
c*********************************************************************
c
	parameter(badflag=-99.9)
	real ta(m),tb(m),dta(m)
	integer onoff
c
	if (onoff.eq.1) then
	   do i=1,imax
	      if(ta(i).ne.badflag.and.tb(i).ne.badflag) then
		 dta(i)=ta(i)-tb(i)+offset
	      else
		 dta(i)=badflag
	      endif
	   enddo !i
	else
	   do i=1,imax
	      ta(i)=tb(i)+(dta(i)-offset)
	   enddo !i
	endif
c
	return
	end
c
c        
	Subroutine weights(mwt,dwt,ua,va,theta,lat,lon,icycle,imax,m)
c
c*********************************************************************
c
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       24 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c       14 Dec 1999  Peter Stamus, NOAA/FSL
c          New version (minor changes from old).
c
c     Notes:
c
c*********************************************************************
c
	common tab(10000),pi,re,rdpdg,reorpd
	real mwt(m,m),dwt(m,m),ua(m),va(m),theta(m),lat(m),lon(m)
        integer icycle
c       
        cycle = float( icycle )
	expsclr=1000.
	hcutoff= 300000.  !dist in m for wt e**-1
	hcutoffm=300000.  
	vcutoff=5.	 !deg C for vertical wt to go to e**-1
	C3=2.		 !wt on wind factor
	c1=sqrt(expsclr)/hcutoff**2
	c1m=sqrt(expsclr)/hcutoffm**2
	c2=hcutoff**2/vcutoff**2
	c2m=hcutoffm**2/vcutoff**2
c
c corrects for wind component along station/station vector
c
	do i=1,imax
	do j=1,imax
	   ang=(lat(i)+lat(j))*.5*rdpdg
	   dy=(lat(j)-lat(i))*reorpd
	   dx=(lon(j)-lon(i))*reorpd*cos(ang)
	   rr=(dx*dx+dy*dy)+1.e-20
	   r=sqrt(rr)
	   uav=dx*(ua(i)+ua(j))*.5/r
	   vav=dy*(va(i)+va(i))*.5/r
	   dxx=dx+uav*cycle/c3
	   dyy=dy+vav*cycle/c3
	   dz=theta(j)-theta(i)
	   iii=int(c1*((dx*dx+dy*dy) +c2*dz*dz))+1
	   iiii=int(c1m*((dxx*dxx+dyy*dyy) +c2m*dz*dz))+1
	   if(iii.gt.10000) iii=10000
	   if(iiii.gt.10000) iiii=10000 
	   dwt(i,j)=tab(iii)
	   mwt(i,j)=tab(iiii)
	enddo !j
	enddo !i
c
c normalize dwt,mwt
c
	do i=1,imax
	   sum=0.
	   summ=0.
	   do j=1,imax
	      sum=dwt(i,j)+sum
	      summ=mwt(i,j)+summ
	   enddo !on j
	   do j=1,imax
	      dwt(i,j)=dwt(i,j)/sum
	      mwt(i,j)=mwt(i,j)/summ
	   enddo !on j
	enddo !on i
c
	return
	end
c
c
        subroutine project(ord,y,monster,nv,nvar,maxsta,m,ncycles,
     &                     nn,atime,it,icyc,oberr,badthr)     
c
c*********************************************************************
c
c     This routine uses Taylor series of order 'ord' to forecast 
c     the variable 'nv' one cycle.  
c
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       24 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c       14 Dec 1999  John McGinley and Peter Stamus, NOAA/FSL
c          Completely new version.
c
c     Notes:
c          For the time being (14 Dec 1999) 'ord' is limited to 2.
c
c*********************************************************************
c
        parameter(badflag=-99.9)
        common tab(10000),pi,re,rdpdg,reorpd
        integer ord
        real y(m), monster(m,ncycles,nvar)
        character*24 atime
c
c.....  Apply a filter to each derivative
c
        em1=exp(-1.)
        em1s=em1*em1
        do i=1,maxsta
           a=1.
           b=1.
           if(monster(i,2,nv).eq.badflag) a=0.
           if(monster(i,3,nv).eq.badflag) b=0.
           sum0 = monster(i,1,nv)*(1.+a*em1+b*em1s*.5) + 
     &            a*monster(i,2,nv)*(-a*em1-b*em1s) + 
     &            monster(i,3,nv)*(b*em1s*.5)    
           y(i)=sum0
        enddo !on i
c
        return 
        end

