c
c
        Subroutine model(F,dta,var,ua,va,theta,mwt,lat,lon,elev
     1   ,imax,m,time,cycle,amp,amp1,peak,peak1)
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
c var is the variable for advective term
        parameter (badflag=-99.9)
        parameter(im=500,jm=500)
        common tab(10000),pi,re,rdpdg,reorpd
        real F(m,m),dta(m),mwt(m,m),lat(m),lon(m),elev(m)
        real advect(im)  ,var(m) 
        real kappa,ua(m),va(m),theta(m) 
c
        kappa=2./7.
        stdlps=-6.5/1000.
        Tstd=287.16
        pstd=1013.2
        g=9.808
c
c set offset for heating: max temp 00Z; min temp 12Z
c
         ddiur=-amp*(2.*pi/(24.*cycle))
     &      *sin(2.*pi/24.*(time-peak))*cycle
c
c      correcting waves
c
     &  -amp1*2.*pi/12.*sin(2.*pi/12.*float(time-peak1))
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
c
              sum=-(uav*dth*dx/rr+vav*dth*dy/rr)*cycle+sum
              cnt=cnt+1.
 1         enddo !on j
           advect(i)=sum/cnt
        enddo !on i
c
c mwt is a model spatial correlation to input what the expected
c spatial error is likely to be with the model
c
        do j=1,imax
        do i=1,imax
           F(i,j)=(1.+(ddiur+advect(i))/dta(i) )*mwt(i,j)
        enddo!on i 
        enddo!on j
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
      Subroutine weights(mwt,dwt,ua,va,theta,lat,lon,
     & cycle,imax,m)
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
      common tab(10000),pi,re,rdpdg,reorpd
      real mwt(m,m),dwt(m,m),ua(m),va(m),theta(m),lat(m),lon(m)
c
      expsclr=1000.
      hcutoff= 50000. !dist in m for wt e**-1
      hcutoffm=50000.
      vcutoff=5. !deg C for vertical wt to go to e**-1
      C3=2. !wt on wind factor
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
         iii=int(c1*((dxx*dxx+dyy*dyy) +c2*dz*dz))+1
         iiii=int(c1m*((dxx*dxx+dyy*dyy) +c2m*dz*dz))+1
         if(iii.gt.10000) iii=10000
         if(iiii.gt.10000) iiii=10000 
         dwt(i,j)=tab(iii)
         mwt(i,j)=tab(iiii)
      enddo !on j
      enddo !on i
c normalize dwt,mwt
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
        subroutine fouranal(monster,fcf,maxsta,m,maxvar,nvar,nn,icyc)
c
c*********************************************************************
c
c     Routine to perform fourier analysis on nvar variables.
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
        parameter       (badflag = -99.9)
        parameter (im=500,jm=500)
        real monster (m,m,nvar), fcf(m,m,nvar)
        real data(im)
c
        step=float(icyc)/float(nn)   
        do k=1,maxvar
           do i=1,maxsta
              do kk=1,nn
                 nn2=2*kk
                 x=step*float(kk)
                 do ll=2,icyc
                    dx=x-float(ll)   
                    if(dx.le.0.) then
                       data(nn2-1)=(1.+dx)*monster(i,icyc-ll+1,k)
     &                      -dx*monster(i,icyc-ll+2,k)
                       data(nn2)=0
                       if(monster(i,icyc-ll+2,k).eq.badflag) then
                          data(nn2-1)=badflag
                          fcf(i,1,k)=badflag
                          go to 2
                       endif 
                       go to 1
                    endif
                 enddo !on ll
 1            enddo !on kk
              call four1(data,nn,1)
              do l=1,2*nn,2
                 fcf(i,l,k)=data(l)/float(nn)
                 fcf(i,l+1,k)=data(l+1)/float(nn)
              enddo !l
 2         enddo !on i
        enddo !on k
c
        return
        end
c
c
       subroutine four1(data,nn,isign)
c
c*********************************************************************
c
c     Routine to find the complex fourier coefs for the complex 
c     data array data(nn).  
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
       integer isign,nn
       real data(2*nn)
       integer i,istep,j,m,mmax,n
       real tempi, tempr
       double precision theta,wi,wpi,wpr,wr,wtemp
c
       n=2*nn
       j=1
       do i=1,n,2
          if(j.gt.i) then
             tempr=data(j)
             tempi=data(j+1)
             data(j)=data(i)
             data(j+1)=data(i+1)
             data(i)=tempr
             data(i+1)=tempi
          endif
          m=n/2
 1        if((m.ge.2).and.(j.gt.m)) then
             j=j-m
             m=m/2
             go to 1
          endif
          j=j+m
       enddo !on i
       mmax=2
  2    if(n.gt.mmax) then
          istep=2*mmax
          theta=6.28318530717959d0/(isign*mmax)
          wpr=-2.d0*sin(0.5d0*theta)**2
          wpi=sin(theta)
          wr=1.d0
          wi=0.d0
          do m=1,mmax,2
             do i=m,n,istep
                j=i+mmax
                tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
                tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
                data(j)=data(i)-tempr
                data(j+1)=data(i+1)-tempi
                data(i)=data(i)+tempr
                data(i+1)=data(i+1)+tempi
             enddo !i
             wtemp=wr
             wr=wr*wpr-wi*wpi+wr
             wi=wi*wpr+wtemp*wpi+wi
          enddo !m
          mmax=istep
          go to 2
       endif
c
       return
       end
c
c
         subroutine project(fcf,y,monster,nv,nvar,
     &                maxsta,m,nn,atime,it,icyc,oberr,badthr)     
c
c*********************************************************************
c
c     Routine uses fourier series to forecast the variable
c     nv one cycle.
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
         common tab(10000),pi,re,rdpdg,reorpd
         real fcf(m,m,nvar),y(m),monster(m,m,nvar),fff(200)
         character*24 atime
c     
         step=float(icyc)/float(nn)  
         thresh=3.*oberr
         em1=1.-.36789
         em1s=em1*em1
         nn2=2*nn
         open(9,file='/scratch/localhost/mcginley/fort.9',form=
     &        'formatted',access='sequential',status='old')
         do i=1,maxsta
c     if(it.lt.104)go to 1
            if(fcf(i,1,nv).eq.badflag) go to 1
            k=icyc+1
            t=float(k-1)/step
            sum=0
            do l=1,nn-1,2
               f=(float(l+1)/2.-1.)/float(nn)
               sum=sum+fcf(i,l,nv)*cos(2.*pi*f*t)
     &              +fcf(i,l+1,nv)*sin(2.*pi*f*t)
            enddo !l
c
c filtering
c
            f=.5
            l=nn+1
            sum=sum+fcf(i,l,nv)*cos(2.*pi*f*t)
     &           +fcf(i,l+1,nv)*sin(2.*pi*f*t)
            sum=sum+fcf(i,l,nv)*cos(2.*pi*(-f)*t)
     &           +fcf(i,l+1,nv)*sin(2.*pi*(-f)*t)
            do l=nn+3,nn2-1,2
               f=-(float(nn2+3-l)/2.-1.)/float(nn)
               sum=sum+fcf(i,l,nv)*cos(2.*pi*f*t)
     &              +fcf(i,l+1,nv)*sin(2.*pi*f*t)
            enddo !l
            if(it.gt.66.and.it.lt.84)
     &           write(9,*) 'STN',i,'VARIABLE ',nv,'ITERATION',it
            fff(1)=fcf(i,1,nv)
            ll=1
c     do l=3,nn+1,2
            do l=3,9,2          !filtering applied - cut off > 1/3hr freq
               ll=ll+1
               f=(float(l+1)/2.-1.)/float(nn)
               fm=-f
               fff(ll)=fff(ll-1)+fcf(i,l,nv)*cos(2.*pi*f*t)
     &              +fcf(i,l+1,nv)*sin(2.*pi*f*t)
     &              +fcf(i,nn2-l+2,nv)*cos(2.*pi*fm*t)
     &              +fcf(i,nn2-l+1,nv)*sin(2.*pi*fm*t)
               sumf=fff(ll)
            enddo !l
c           if(it.lt.90.and.it.gt.66) then
c           do l=1,nn
c           write(9,*) 'l,mon,fff',l,monster(i,l,nv),fff(l)
c           fff(l)=0.
c           enddo
c           write(9,*)  'fourier series- unfil, fil',sum,sumf 
c           endif
 1          a=1.
            b=1.
            if(monster(i,2,nv).eq.badflag) a=0.
            if(monster(i,3,nv).eq.badflag) b=0
            sum0=monster(i,1,nv)*(1.+a*em1+b*em1s*.5)+a*monster(i,2,nv)*
     &           (-a*em1-b*em1s)+monster(i,3,nv)*(b*em1s*.5)    
            y(i)=sum0
            if(it.lt.26) go to 2
c     if(abs(sum-sum0).lt.thresh)then 
c           write(8,*) 'fft est stn ',i, ' of ',sum,' acptd ovr ',
c    &       '2nd order interp  ',sum0,atime(1:17)
c           y(i)= sum   
c           else
c           y(i)=sum0
c           write(8,*) 'fft est stn ',i, ' of ',sum,' RJCTD FOR ',
c    &       '2nd ordr  ',y(i),'  ',atime(1:17)
c           endif
            write(9,*) 'i,fourier unfil, filt est',i,sum,sumf
            write(9,*) 'taylor series',sum0,monster(i,1,nv),
     &           monster(i,2,nv),
     &           monster(i,3,nv)
            y(i)=sum0           !ob guess set to taylor for the time being
 2          continue
         enddo !on i
c
         return 
         end

