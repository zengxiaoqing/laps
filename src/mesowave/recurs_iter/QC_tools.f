c
c
        subroutine score(TT,T,BB,XB,XT,TA,imax,scf,sca,scb,sco,
     &                   thresh,qcstat,m,badflag)
c
c*********************************************************************
c
c     This routine computes QC efficiency.
c     
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       24 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c       15 Dec 1999  John McGinley and Peter Stamus, NOAA/FSL
c          New version.
c
c     Notes:
c
c*********************************************************************
c     
        real    TT(m),XT(m),TA(m),T(m),XB(m),BB(m),thresh
        integer scf(2,2),sca(2,2),scb(2,2),sco(2,2)
        integer i,ii,jj,jjj,qcstat(m),iii 
c
        do i=1,imax
           if(qcstat(i).eq.1) go to 1
           if(abs(TT(i)-T(i)).gt.thresh ) then !ob is bad
              jj=2
           else
              jj=1
           endif
           if(abs(T(i)-XT(i)).gt.thresh) then !model tosses
              ii=2
              qcstat(i)=qcstat(i)+10000
           else
              ii=1
           endif
           if(abs(T(i)-XB(i)).gt.thresh) then !buddy check tosses ob
              iii=2
              qcstat(i)=qcstat(i)+1000
           else
              iii=1
           endif
           if(abs(T(i)-TA(i)).gt.thresh) then !kalman optimum tosses 
              jjj=2
              qcstat(i)=qcstat(i)+100000
           else
              jjj=1
           endif
           if(abs(T(i)-BB(i)).gt.thresh) then !obtend tosses 
              iiii=2
              qcstat(i)=qcstat(i)+100
           else
              iiii=1
           endif
           scf(ii,jj)=scf(ii,jj)+1
           sca(jjj,jj)=sca(jjj,jj)+1
           scb(iii,jj)=scb(iii,jj)+1
           sco(iiii,jj)=sco(iiii,jj)+1
 1      enddo
c
        return
        end 
c
c
      function ffz(ii,n)
c
c*********************************************************************
c
c     Pulls out a random value from a normal distribution.  Input is
c     an integer seed 'ii' and an interation number 'n'. 'n' should
c     be >20 for best results.
c     
c     Original: John McGinley, NOAA/FSL  Spring 1998
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
         imix=9.*ran1(ii)+1
         do i=1,imix
            ii=ran1(ii)*20000000.-10000000.
         enddo !i
         sum=sum+ran1(ii)
      enddo !l
      ffz=sum-float(n)/2.
c
      return
      end
c
c
      subroutine stats_qc(sc,m,n)
c
c*********************************************************************
c
c     
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       21 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c
c     Notes:
c
c*********************************************************************
c
      real pgot,tss,pbok,pbot,pgok
      integer sc(m,n)
c
      if((sc(1,2)+sc(2,2)).eq.0) then
         pbok=-99.9
         pbot=-99.9
      else
         pbok=float(sc(1,2))/float(sc(1,2)+sc(2,2))
         pbot=1.-pbok
      endif
      if((sc(1,1)+sc(2,1)).eq.0) then
         pgot=-99.9
         pgok=-99.9
      else
         pgot=float(sc(2,1))/float(sc(1,1)+sc(2,1))
         pgok=1.-pgot
      endif
      if ((sc(1,2)+sc(2,1)+sc(2,2)).eq.0.) then
         csi=-99.9
      else
         csi=float(sc(2,2))/float(sc(1,2)+sc(2,1)+sc(2,2))
      endif
      if((sc(2,2)+sc(1,2)).eq.0..or.(sc(1,1)+sc(2,1)).eq.0.) then
         tss=-99.9
      else
         tss=float(sc(2,2))/float(sc(2,2)+sc(1,2))-float(sc(2,1))
     1        /float(sc(1,1)+sc(2,1))
      endif 
      write(6,1000) pbot,pbok,pgot,pgok,csi,tss
 1000 format(/1x,'Fraction of bad obs tossed/kept ',2f6.3/
     1     1x,'Fraction of good obs tossed/kept ',2f6.3/
     2     1x,'CSI = ',f6.3,'  TSS = ',f6.3)
c
      return
      end
c
c
        Subroutine thet(t,theta,maxstn,elev,m)         
c
c*********************************************************************
c
c     The theta used here is based on height only with an estimate of 
c     and average temperature which is valid for the domain: theta is C
c     
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       24 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c       15 Dec 1999  Peter Stamus, NOAA/FSL
c          Housekeeping...change 'if' structure in do loop.
c
c     Notes:
c
c*********************************************************************
c
        parameter (badflag=-99.9)
        real t(m),theta(m),elev(m)
        real kappa
c
        kappa=2./7.
        stdlps=-6.5/1000.
        Tstd=287.16
        pstd=1013.2
        g=9.808
        rr=287.04
        do l=1,maxstn
           if(t(l).eq.badflag) then
              theta(l)= badflag
           else
              tbar=2.*Tstd+stdlps*elev(l)
              ppp=pstd*exp(-elev(l)*g/rr/tbar)
              theta(l)=((t(l)-32.)*.5555556+273.16) * 
     &                 (1000./ppp)**kappa - 273.16
           endif
        enddo !on l
c
        return
        end
c
c
        Subroutine thet2T(t,theta,maxstn,elev,m)          
c
c*********************************************************************
c
c     This subroutine converts from theta (C) to Temp(F).                    
c     
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       24 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c       15 Dec 1999  Peter Stamus, NOAA/FSL
c          Housekeeping...restructure 'if' block in do loop.
c
c     Notes:
c
c*********************************************************************
c
        parameter(badflag=-99.9)
        real t(m),theta(m),elev(m)
        real kappa
c
        kappa=2./7.
        stdlps=-6.5/1000.
        Tstd=287.16
        pstd=1013.2
        g=9.808
        rr=287.04
        do l=1,maxstn
           if(theta(l).eq.badflag) then
              t(l)=badflag
           else 
              tbar=2.*Tstd+stdlps*elev(l)
              ppp=pstd*exp(-elev(l)*g/rr/tbar)
              t(l) = ( (ppp/1000.)**(kappa) * 
     &               (theta(l)+273.16) - 273.16 ) / .5555556 + 32.
           endif
        enddo !on l
c
        return
        end
c
c
         Subroutine convuv(dd,ff,u,v,maxsta,m,badflag)
c
c*********************************************************************
c
c     
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       24 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c       15 Dec 1999  Peter Stamus, NOAA/FSL
c          Housekeeping...restructure 'if' block in do loop.
c
c     Notes:
c
c*********************************************************************
c
         real dd(m),ff(m),u(m),v(m)
c
         rdpdg=3.14159122/180.
         do l=1,maxsta
            if(dd(l).le.badflag.or.ff(l).le.badflag) then
               u(l)=badflag
               v(l)=badflag
            else
               ang=dd(l)-270.
               ang=rdpdg*ang
               u(l)=cos(ang)*ff(l)
               v(l)=-sin(ang)*ff(l)
            endif
         enddo !on l
c
         return
         end 
c
c
         Subroutine reorder(ta,tda,dda,ffa,lata,lona,eleva,pstna,
     &     pmsla,alta,stna,providera,reptypea,indexa,
     &     tb,tdb,ddb,ffb,latb,lonb,elevb,pstnb,pmslb,altb,stnb,
     &     providerb,reptypeb,indexb,
     &     maxstaa,maxstab,m,badflag)
c
c*********************************************************************
c
c     
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       24 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c       07 Oct 1998  Peter Stamus, NOAA/FSL
c          Added provider and report type variables.
c       22 Oct 1998  Peter Stamus, NOAA/FSL
c          Added index variables to track ob location in arrays.
c       12 Mar 1999  Peter Stamus, NOAA/FSL
c          Fixed bug with dup station handling.
c
c     Notes:
c
c*********************************************************************
c
         real ta(m),tda(m),dda(m),ffa(m),lata(m),lona(m),eleva(m)
         real pstna(m),pmsla(m),alta(m)
         real tb(m),tdb(m),ddb(m),ffb(m),latb(m),lonb(m),elevb(m)
         real pstnb(m),pmslb(m),altb(m)
         integer indexa(m), indexb(m), iholdx, ix
         character stna(m)*5,stnb(m)*5,holdnam*5, ch*1
         character providera(m)*11, providerb(m)*11, holdprov*11
         character reptypea(m)*6, reptypeb(m)*6, holdrep*6
c
         holdnam = '     '
         max=maxstaa
         if(maxstab.gt.max)max=maxstab
c
c first check the master list to find duplicates
c
         do k=1,maxstab
         do l=1,maxstab
            if(k.ne.l.and.stnb(k).eq.stnb(l)) then
               if(latb(k).ne.latb(l).or.lonb(k).ne.lonb(l)) then
c
c     we've found a duplicate name
c
                  holdnam = stnb(l)
                  stnb(l)(5:5)='A'
                  write(*,1021) holdnam, stnb(l)
 1021             format(
     &      1x,'Station in master list renamed - old/new: ',a5,1x,a5)
               endif 
            endif
         enddo !l
         enddo !k
c         do k=1,max
c         write(*,1000) k,stna(k),lata(k),lona(k),stnb(k),latb(k),
c    1    lonb(k)
c         enddo
 1000    format(1x,i3,1x,a5,1x,f8.3,1x,f8.3,1x,a5,1x,f8.3,1x,f8.3)
c
         holdnam = '     '
         do k=1,maxstaa
 3          continue
            iflag=0
            do l=1,maxstab
               if(stnb(l)(1:5).eq.stna(k)(1:5)) then
                  if(lata(k).eq.latb(l).and.lona(k).eq.lonb(l)) then 
                     iflag=1    !the station is on the master list
                  else  
                     holdnam = stna(k) !diff stn has a duplicate name
                     ix = ichar(stna(k)(5:5))
                     if(ix .lt. 64) ix = 64
                     ix = ix + 1
                     if(ix .gt. 90) then
                        print *,
     &            ' WARNING. Over 26 dup stn names for ', stna(k)
                        print *,' Will try to continue...'
                     endif
                     ch = char(ix)
                     stna(k)(5:5) = ch
                     write(*,1020) holdnam,stna(k)
 1020                format(1x,
     &  'Master/newlist duplicate station renamed - old/new: ',a5,1x,a5)
                     go to 3 
                  endif
               endif
            enddo !l
            if(iflag.eq.0) then !the station needs to be added to the list
               maxstab=maxstab+1
               indexb(maxstab) = indexa(k)
               stnb(maxstab)=stna(k)
               latb(maxstab)=lata(k)
               lonb(maxstab)=lona(k)
               elevb(maxstab)=eleva(k)
               providerb(maxstab) = providera(k)
               reptypeb(maxstab) = reptypea(k)
               tb(maxstab)=badflag
               tdb(maxstab)=badflag
               ddb(maxstab)=badflag
               ffb(maxstab)=badflag
               pstnb(maxstab)=badflag
               pmslb(maxstab)=badflag
               altb(maxstab)=badflag
            endif
         enddo !on k
c
c now get both lists in the same order with the same stations
c
         holdnam = '     '
         do k=1,maxstab
            iflag=0
            do l=1,maxstaa
               if (stna(l).eq.stnb(k)) then
                  iflag=1
                  iholdx = indexa(k)
                  holdt=ta(k)
                  holdtd=tda(k)
                  holddd=dda(k)
                  holdff=ffa(k)
                  holdlt=lata(k)
                  holdln=lona(k)
                  holdel=eleva(k)
                  holdps=pstna(k)
                  holdpm=pmsla(k)
                  holdal=alta(k)
                  holdnam=stna(k)
                  holdprov = providera(k)
                  holdrep = reptypea(k)
                  indexa(k) = indexa(l)
                  ta(k)=ta(l)
                  tda(k)=tda(l)
                  dda(k)=dda(l)
                  ffa(k)=ffa(l)
                  lata(k)=lata(l)
                  lona(k)=lona(l)
                  eleva(k)=eleva(l)
                  pstna(k)=pstna(l)
                  pmsla(k)=pmsla(l)
                  alta(k)=alta(l)
                  stna(k)=stna(l)
                  providera(k) = providera(l)
                  reptypea(k) = reptypea(l)
                  indexa(l) = iholdx
                  ta(l)=holdt
                  tda(l)=holdtd
                  dda(l)=holddd
                  ffa(l)=holdff
                  lata(l)=holdlt
                  lona(l)=holdln
                  eleva(l)=holdel
                  pstna(l)=holdps
                  pmsla(l)=holdpm
                  alta(l)=holdal
                  stna(l)=holdnam
                  providera(l) = holdprov
                  reptypea(l) = holdrep
                  go to 1
               endif
            enddo !l
            if(iflag.eq.0) then ! there is no match for an existing stn
c make a place for the master list station in putit in the new obs
               maxstaa=maxstaa+1
               do l =maxstaa,k+1,-1
                  indexa(l) = indexa(l-1)
                  stna(l)=stna(l-1)
                  providera(l) = providera(l-1)
                  reptypea(l) = reptypea(l-1)
                  tda(l)=tda(l-1)
                  ta(l)=ta(l-1)
                  dda(l)=dda(l-1)
                  ffa(l)=ffa(l-1)
                  lata(l)=lata(l-1)
                  lona(l)=lona(l-1) 
                  eleva(l)=eleva(l-1)
                  pstna(l)=pstna(l-1)
                  pmsla(l)=pmsla(l-1)
                  alta(l)=alta(l-1)
               enddo !on l
               indexa(k) = maxstaa
               stna(k)=stnb(k)
               providera(k) = providerb(k)
               reptypea(k) = reptypeb(k)
               ta(k)=badflag 
               tda(k)=badflag 
               dda(k)=badflag 
               ffa(k)=badflag 
               lata(k)=latb(k) 
               lona(k)=lonb(k) 
               eleva(k)=elevb(k)
               pstna(k)=badflag 
               pmsla(k)=badflag
               alta(k)=badflag 
            endif
 1       enddo !k
c
         return
         end
c
c
        subroutine fill(stn,t,lat,lon,elev,theta,
     &     lapse,maxsta,m,tmx,tmn)
c
c*********************************************************************
c
c     Routine fills array t with an interpolated value if missing
c     uses a pass for theta analysis then the variable
c     tmx and tmn are the extremes of the variable lapse rate in unit/m)
c     
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       24 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c       15 Dec 1999  John McGinley and Peter Stamus, NOAA/FSL
c          Add to write code 'if' block (1000 format).
c
c     Notes:
c
c*********************************************************************
c
        parameter(badflag=-99.9)
        common tab(10000),pi,re,rdpdg,reorpd 
        real t(m),lat(m),lon(m),elev(m),theta(m)
        real lapse
        character stn(m)*5
c
        expsclr=1000.
        hcutoff=200000.         !dist in m for wt e**-1
        vcutoff=5.              !ht(m) in vertical wt to go to e**-1
        c1=sqrt(expsclr)/hcutoff**2
        c2=hcutoff**2/vcutoff**2
        thmx=8./1000.           ! C/m
        thmn=-1./1000.
        call regress(elev,theta,maxsta,m,thlapse,sint,thmx,thmn)
        call regress(elev,t,maxsta,m,lapse,sintt,tmx,tmn)
c
c compute avg variable 
        tbar=0.
        ebar=0.
        cnt=0.
        do i=1,maxsta
           if (t(i).ne.badflag) then
              cnt=cnt+1.
              tbar=tbar+t(i)
              ebar=ebar+elev(i)
           endif
        enddo !i
        tbar=tbar/cnt
        ebar=ebar/cnt
        do i=1,maxsta 
           if(t(i).eq.badflag) then
              sumwt=0.
              sum=0.
              do j = 1,maxsta
                 if(t(j).eq.badflag) go to 3
                 ang=(lat(i)+lat(j))*.5*rdpdg
                 dy=(lat(j)-lat(i))*reorpd
                 dx=(lon(j)-lon(i))*reorpd*cos(ang)
                 dz=theta(j)-theta(i)
                 if(theta(j).eq.badflag.or.theta(i).eq.badflag) then
                    dz=(elev(j)-elev(j))*thlapse   
                 endif
                 iii=int(c1*((dx*dx+dy*dy) +c2*dz*dz))+1
                 if(iii.gt.10000) iii=10000
                 sum=tab(iii)*(t(j)-(tbar+lapse*(elev(j)-ebar)))+sum     
                 sumwt=tab(iii)+sumwt
 3            enddo !on j
              if(sumwt .ne. 0.) then
                 t(i)=sum/sumwt+tbar+(elev(i)-ebar)*lapse
                 write(6,1000) i, stn(i), t(i)
 1000            format(1x,' Replaced missing ob at stn ',i5,1x,a5,
     &                  ' with ',F8.3) 
              else
                 write(6,*) ' No fill possible for station ', i
              endif
           endif 
        enddo !i
c
        return
        end
c
c
        subroutine fillone(stn,t,lat,lon,elev,theta,
     &    lapse,maxsta,m,tmx,tmn)
c
c*********************************************************************
c
c     Returns a buddy check value for a station location not using 
c     the station value.  Instead, it uses a pass for theta analysis 
c     then the variable tmx and tmn are the extremes of the variable 
c     lapse rate in unit/m .
c     
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       24 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c       04 Aug 1999  Add missing calls to 'regress'.
c
c     Notes:
c
c*********************************************************************
c
        parameter(badflag=-99.9)
        common tab(10000),pi,re,rdpdg,reorpd 
        real t(m),tfl(m),lat(m),lon(m),elev(m),theta(m)
        real lapse
        character stn(m)*5
c     
        expsclr=1000.
        hcutoff=200000.         !dist in m for wt e**-1
        vcutoff=5.              !ht(m) in vertical wt to go to e**-1
        c1=sqrt(expsclr)/hcutoff**2
        c2=hcutoff**2/vcutoff**2
        thmx=8./1000.           ! C/m
        thmn=-1./1000.
        call regress(elev,theta,maxsta,m,thlapse,sint,thmx,thmn)
        call regress(elev,t,maxsta,m,lapse,sintt,tmx,tmn)
c
c     compute avg variable 
        tbar=0.
        ebar=0.
        cnt=0.
        do i=1,maxsta
           if (t(i).ne.badflag) then
              cnt=cnt+1.
              tbar=tbar+t(i)
              ebar=ebar+elev(i)
           endif
        enddo !i
        tbar=tbar/cnt
        ebar=ebar/cnt
        do i=1,maxsta
           sumwt=0.
           sum=0.
           do j = 1,maxsta
              if (j.eq.i) go to 3
              if(t(j).eq.badflag) go to 3
              ang=(lat(i)+lat(j))*.5*rdpdg
              dy=(lat(j)-lat(i))*reorpd
              dx=(lon(j)-lon(i))*reorpd*cos(ang)
              dz=theta(j)-theta(i)
              if(theta(j).eq.badflag.or.theta(i).eq.badflag) then
                 dz=(elev(j)-elev(j))*thlapse   
              endif
              iii=int(c1*((dx*dx+dy*dy) +c2*dz*dz))+1
              if(iii.gt.10000) iii=10000
              sum=tab(iii)*(t(j)-(tbar+lapse*(elev(j)-ebar)))+sum     
              sumwt=tab(iii)+sumwt
 3         enddo !on j
           if(sumwt.ne.0) tfl(i) = sum/sumwt + tbar + 
     &                             (elev(i) - ebar) * lapse
        enddo !on i
        do i=1,maxsta 
           t(i)=tfl(i)
        enddo !i
c
        return
        end
c
c
        subroutine regress(elev,t,imax,m,slp,sint,smx,smn)
c     
c*********************************************************************
c
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
        parameter (badflag=-99.9)
        real t(m),elev(m)
c     
        sum=0.
        sum1=0.
        sum2=0.
        sum12=0.
        sum11=0.
        do i=1,imax
           if(t(i).eq.badflag) go to 1
           sum12=t(i)*elev(i)+sum12
           sum1=elev(i)+sum1
           sum2=t(i)+sum2
           sum11=elev(i)*elev(i)+sum11
           sum=sum+1.
 1      enddo !i
        slp=(sum*sum12-sum1*sum2)/(sum*sum11-sum1*sum1)    
        sint=(sum2-slp*sum1)/sum 
        if (slp.lt.smn) then
           write(6,1001) slp,smn
           slp=smn
        endif
        if(slp.gt.smx) then
           slp = smx
c       write(6,1001) slp,smx
c           write(6,1001) slp,sint
        endif
 1001   format(1x,' slope of ',e12.3,' being replaced by set limit ',
     &      e12.3)
c
        return
        end
c
c
        subroutine qcset(t,ta,ncm,qcstat,m,maxsta,badflag,gross)
c
c*********************************************************************
c
c     qcstat is 1 if ob is missing; 10 if failed gross err;
c     100 if ob test fails; 1000 if model test fails; 
c     10000 opt check fails
c     
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       24 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c       15 Dec 1999  John McGinley and Peter Stamus, NOAA/FSL
c          New version.
c
c     Notes:
c
c*********************************************************************
c
        real t(m), ta(m), ncm(m)
        integer qcstat(m)
c
        do i=1,maxsta
           qcstat(i) = 0
           if(t(i) .eq. badflag) then
              qcstat(i) = 1
              ncm(i) = ncm(i) + 1.
           else
              if(ta(i) .ne. badflag) then
                 if(abs(t(i)-ta(i)) .gt. gross) then
                    qcstat(i) = 10
                    write(6,*) ' qcset:stn ',i,t(i),
     &                         ' fails gross error ck relative to ',
     &                         'estimate ',ta(i),' :ob replaced.'
                    t(i) = ta(i)
                    ncm(i) = ncm(i) + 1.
                 else
                    ncm(i) = 0.
                 endif
              endif
           endif
        enddo !i
c
        return
        end
c
c
        subroutine replace(a,b,imax,jmax,m,n)
c
c*********************************************************************
c
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
        real a(m,n),b(m,n)    
c
        do i=1,imax
        do j=1,jmax
           b(i,j)=a(i,j)
        enddo !j
        enddo !i
c
        return
        end 
c
c
        subroutine replacei(ia,ib,imax,m)
c
c*********************************************************************
c
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
        integer ia(m),ib(m)
c
        do i=1,imax
           ib(i)=ia(i)
        enddo !i
c
        return
        end 
c
c
        Subroutine errorproc(dt,y,by,md,xt,x,wr,vr,ar,imax,m,qcstat,
     &                  oberr,wm,wo,wb,length,badthr,atime,icnt,nvar)
c
c*********************************************************************
c
c     
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       24 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c       07 Oct 1998  Peter Stamus, NOAA/FSL
c          Change 'ran' to 'ran1' function for portability.
c       15 Dec 1999  John McGinley and Peter Stamus, NOAA/FSL
c          New version.
c
c     Notes:
c
c*********************************************************************
c
        parameter(badflag=-99.9)    
        real dt(m),y(m),xt(m),by(m),md(m),x(m),wr(m),vr(m),ar(m)
        real me,oe,ce
        real wm(m),wo(m),wb(m)
        real tdt(m),b(m),tmb(m),tob(m),tcb(m),tmr(m),tor(m),tcr(m)
        real totr(m),totb(m)
        integer sca(2,2),scf(2,2),scb(2,2),sco(2,2),imax,qcstat(m)
        integer scat(2,2),scft(2,2),scbt(2,2),scot(2,2),on,off
c
c.....  Truth routine and missing ob replacement
c
        icnt=icnt+1 
        iiiii=icnt
        iiiii=ran1(iiiii)*20000000.-10000000.
        on=1
        off=0
        if(totr(nvar).ne.0) then
           me=sqrt(tmr(nvar)/totr(nvar))
           oe=sqrt(tor(nvar)/totr(nvar))
           ce=sqrt(tcr(nvar)/totr(nvar))
        endif
        thresh=oberr*badthr 
        do i=1,imax
           iiiii=ran1(iiiii)*2000000.-1000000.
           if (dt(i).ne.badflag) then
              tdt(i)=dt(i)-oberr*ffz(iiiii,20)
           else
              tdt(i) =
     &           .3333*(xt(i)-me*ffz(iiiii,20)) +
     &           .3333*( x(i)-ce*ffz(iiiii,20)) +
     &           .3333*( y(i)-oe*ffz(iiiii,20))
           endif
        enddo !i
c
c.....  Compute performance stats
c
        stdwr=0.
        stdar=0.
        stdvr=0.
        stdfo=0.
        stdco=0.
        stdoo=0.
        sumfo=0.
        sumco=0.
        sumoo=0.
        sumvr=0.
        sumwr=0.
        sumar=0.
        sumtot=0.
        sumtotb=0.
        sum1=0.
        sum2=0.
        sum3=0.
        do i=1,imax
           sumtot=sumtot+1.
           WR(i)=xt(i)-tdt(i)
           sumwr=sumwr+WR(i)
           stdwr=WR(i)*WR(i)+stdwr
           if(dt(i).ne.badflag) then 
              VR(i)=dt(i)-tdt(i)
           else 
              VR(i)=oberr*(-1.)**(int(10.*ran1(i)+1.))
           endif 
           sumvr=sumvr+VR(i)
           stdvr=VR(i)*VR(i)+stdvr
           AR(i)=X(i)-TdT(i)
           sumar=sumar+AR(i)
           stdar=AR(i)*AR(i)+stdar
c
c  do kalmod errors
           wo(i)=(y(i)-tdt(i))**2+float(length-1)*wo(i)**2
           wo(i)=sqrt(wo(i)/float(length))
           wb(i)=(by(i)-tdt(i))**2+float(length-1)*wb(i)**2
           wb(i)=sqrt(wb(i)/float(length))
           wm(i)=(md(i)-tdt(i))**2+float(length-1)*wm(i)**2
           wm(i)=sqrt(wm(i)/float(length))
           sum1=sum1+wo(i)
           sum2=sum2+wb(i)
           sum3=sum3+wm(i)
c     write(6,*) 'KALMOD ERRORS - OB - BUD - MOD :',i,wo(i),wb(i),wm(i)
           if(dt(i).eq.badflag) go to 5
           sumco=X(i)-dT(i)+sumco
           sumfo=XT(i)-dT(i)+sumfo
           stdfo=(XT(i)-dT(i))**2+stdfo
           stdco=(X(i)-dT(i))**2+stdco
           sumtotb=sumtotb+1.
           b(i)=dt(i)
 5         continue   
        enddo !i
c
        print*,'OVERALL FRACTION OB,BUD,NWP '
        Write(6,*) 'OVERALL FRACTION OB,BUD,NWP '
        print*,.5*(sum2+sum3)/(sum1+sum2+sum3),.5*(sum1+sum3)/
     1       (sum1+sum2+sum3),.5*(sum1+sum2)/(sum1+sum2+sum3) 
        Write(6,*) .5*(sum2+sum3)/(sum1+sum2+sum3),.5*(sum1+sum3)/
     1       (sum1+sum2+sum3),.5*(sum1+sum2)/(sum1+sum2+sum3) 
        tmb(nvar)=tmb(nvar)+sumwr
        tob(nvar)=tob(nvar)+sumvr 
        tcb(nvar)=tcb(nvar)+sumar
        tmr(nvar)=tmr(nvar)+stdwr
        tor(nvar)=tor(nvar)+stdvr
        tcr(nvar)=tcr(nvar)+stdar
        totr(nvar)=totr(nvar)+sumtot
        totb(nvar)=totb(nvar)+sumtot 
        if(sumtot.eq.0.) then
           stdar=badflag
           stdvr=badflag
           stdwr=badflag
           sumvr=badflag
           sumar=badflag
           sumwr=badflag
        else
           stdar=sqrt(stdar/sumtot)
           stdvr=sqrt(stdvr/sumtot)
           stdwr=sqrt(stdwr/sumtot)
           sumvr=sumvr/sumtot
           sumar=sumar/sumtot
           sumwr=sumwr/sumtot
        endif
        if(sumtotb.eq.0.) then
           stdco=badflag
           stdfo=badflag
           stdoo=badflag
           sumco=badflag
           sumoo=badflag
           sumfo=badflag
        else
           stdco=sqrt(stdco/sumtotb)
           stdfo=sqrt(stdfo/sumtotb)
           stdoo=sqrt(stdoo/sumtotb)
           sumco=sumco/sumtotb
           sumfo=sumfo/sumtotb
           sumoo=sumoo/sumtotb
        endif
        write(6,1061) stdwr,sumwr,stdvr,sumvr,stdar,sumar,stdfo,
     &       sumfo,stdco,sumco 
 1061   format(1x,'RMS/Bias errors for truth and observation estimates'/
     &       1x,'F X = XT    - truth RMS ',f8.3,' Bias ',f8.3/
     &       1x,'Observ      - truth RMS ',f8.3,' Bias ',f8.3/
     &       1x,'Kalman X    - truth RMS ',f8.3,' Bias ',f8.3/
     &       1x,'F X = XT    - obser RMS ',f8.3,' Bias ',f8.3/
     &       1x,'Kalman X    - obser RMS ',f8.3,' Bias ',f8.3)
        write(6,*) 'RUNNING BIAS AND RMS BY VARIABLE' 
        write(6,*) 'FX=XT ',(tmb(nvar)/totb(nvar)),sqrt(tmr(nvar)/totr(
     &                            nvar))
        write(6,*) 'OBSER ',(tob(nvar)/totb(nvar)),sqrt(tor(nvar)/totr(
     &                            nvar))
        write(6,*) 'KAL X ',(tcb(nvar)/totb(nvar)),sqrt(tcr(nvar)/totr(
     &                            nvar))
 1000   format(1x,'TRUE VALUES FOR TIME ',i4)
c
c zero out scoring arraYs
        do i=1,2
        do j=1,2
           scf(i,j)=0
           scb(i,j)=0
           sca(i,j)=0 
           sco(i,j)=0
        enddo !j
        enddo !i
        call score(tdt,b,y,by,md,xt,imax,scf,sca,scb,sco,
     1       thresh,qcstat,m,badflag)
c     write(6,*) 'CONTINGECY TABLES FOR ERROR THRESHOLD OF ',thresh 
c     write(6,1011) sco(1,1), sco(1,2), sco(2,1),sco(2,2)
 1011   format(//1x,'QC EFFICIENCY FOR OB TND'/14x,' GOOD OB    BAD OB'
     1 //1x,' KEPT     ',2I10//1x,' TOSSED   ',2I10)
c     call stats_qc(sco,2,2)
c     write(6,1010) scb(1,1), scb(1,2), scb(2,1),scb(2,2)
 1010   format(//1x,'QC EFFICIENCY FOR BUDDY'/14x,' GOOD OB    BAD OB'
     1 //1x,' KEPT     ',2I10//1x,' TOSSED   ',2I10)
c     call stats_qc(scb,2,2)
c     write(6,1008) scf(1,1), scf(1,2), scf(2,1),scf(2,2)
 1008   format(//1x,'QC EFFICIENCY FOR NWP    '/14x,' GOOD OB    BAD OB'
     1 //1x,' KEPT     ',2I10//1x,' TOSSED   ',2I10)
c     call stats_qc(scf,2,2)
c     write(6,1009) sca(1,1), sca(1,2), sca(2,1),sca(2,2)
 1009   format(//1x,'QC EFFICIENCY FOR KAL XT '/14x,' GOOD OB    BAD OB'
     1 //1x,' KEPT     ',2I10//1x,' TOSSED   ',2I10)
        do i=1,2
        do j=1,2
           scft(i,j)=scft(i,j)+scf(i,j)
           scbt(i,j)=scbt(i,j)+scb(i,j)
           scat(i,j)=scat(i,j)+sca(i,j)
           scot(i,j)=scot(i,j)+sco(i,j)
        enddo !j
        enddo !i
c     write(6,1021) scot(1,1), scot(1,2), scot(2,1),scot(2,2)
 1021   format(//1x,'CUM QC EFFNCY FOR OB TND'/14x,' GOOD OB    BAD OB'
     1 //1x,' KEPT     ',2I10//1x,' TOSSED   ',2I10)
c     call stats_qc(scot,2,2)
c     write(6,1020) scbt(1,1), scbt(1,2), scbt(2,1),scbt(2,2)
 1020   format(//1x,'CUM QC EFFNCY FOR BUDDY'/14x,' GOOD OB    BAD OB'
     1 //1x,' KEPT     ',2I10//1x,' TOSSED   ',2I10)
c     call stats_qc(scbt,2,2)
c     write(6,1018) scft(1,1), scft(1,2), scft(2,1),scft(2,2)
 1018   format(//1x,'CUM QC EFFNCY FOR NWP'/14x,' GOOD OB    BAD OB'
     1 //1x,' KEPT     ',2I10//1x,' TOSSED   ',2I10)
c     call stats_qc(scft,2,2)
c     write(6,1019) scat(1,1), scat(1,2), scat(2,1),scat(2,2)
 1019   format(//1x,'CUM QC EFFNCY FOR XT'/14x,' GOOD OB    BAD OB'
     1 //1x,' KEPT     ',2I10//1x,' TOSSED   ',2I10)
c     call stats_qc(scat,2,2)
c
        return
        end
c
c    
        subroutine flag(qcstat,qcdat,m,ntest)
c
c*********************************************************************
c
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
        integer qcstat(m),qcdat(m)
c
        isum=0
        do i=ntest-1,1,-1
           ib=10**i
           isum=isum+ib
           ia=qcstat(i)-isum
           if(ia.ge.0) then
              qcdat(i+1)=1
           else
              qcdat(i+1)=0
           endif
        enddo !i
c
        return
        end
c
c
        Function pythag(a,b)
c
c*********************************************************************
c
c     Function computes (a**2+b**2)**.5 without over or underflow.
c     
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       07 Oct 1998  Peter Stamus, NOAA/FSL
c          Housekeeping changes for use in LAPS.
c
c     Notes:
c
c*********************************************************************
c
        real a,b,pythag
        real absa,absb
c
        absa=abs(a)
        absb=abs(b)
c
        if(absa.gt.absb) then
           pythag=absa*sqrt(1.+(absb/absa)**2)
        else
           if(absb.eq.0.) then
              pythag=0.
           else
              pythag=absb*sqrt(1.+(absa/absb)**2)
           endif
        endif
c
        return
        end
c
c
        function ran1(idum)
c
c*********************************************************************
c
c     Function to generated a random number.  Use this instead of a
c     machine dependent one.
c     
c     Original: Peter Stamus, NOAA/FSL  07 Oct 1998
c     Changes:
c
c     Notes:
c        From "Numerical Recipes in Fortran", page 271.
c
c*********************************************************************
c
        integer idum, ia, im, iq, ir, ntab, ndiv
        real ran1, am, eps, rnmx
        parameter( ia = 16807,
     &             im = 2147483647,
     &             am = 1. / im,
     &             iq = 127773,
     &             ir = 2836,
     &             ntab = 32,
     &             ndiv = 1 + (im-1) / ntab,
     &             eps = 1.2e-7,
     &             rnmx = 1. - eps)
c
        integer j, k, iv(ntab), iy
        save iv, iy
        data iv/ntab * 0/, iy/0/
c
c.....  Start here.
c
        if(idum.le.0 .or. iy.eq.0) then  !initialize
           idum = max(-idum,1)           !prevent idum=0
           do j=ntab+8,1,-1              !load the shuffle table
              k = idum / iq
              idum = ia * (idum - k*iq) - ir*k
              if(idum .lt. 0) idum = idum + im
              if(j .le. ntab) iv(j) = idum
           enddo !j
           iy = iv(1)
        endif
c
        k = idum / iq                     !start here if not initializing
        idum = ia * (idum - k*iq) - ir*k  
        if(idum .lt. 0) idum = idum + im
        j = 1 + iy/ndiv
        iy = iv(j)                        !output prev stored value and 
        iv(j) = idum                      !  refill shuffle table
        ran1 = min(am*iy, rnmx)
c
        return
        end
c
c
        subroutine subob(a,b,c,ncm,h,imax,m,badflag)
c
c*********************************************************************
c
c     Original: John McGinley, NOAA/FSL  December 1999
c     Changes:
c       15 Dec 1999  Peter Stamus, NOAA/FSL
c          Housekeeping changes for use in LAPS.
c
c     Notes:
c
c*********************************************************************
c
        real a(m),b(m),c(m),ncm(m),nsq,h,hsq    
c
        hsq=h*h
        do i=1,imax
           nsq=ncm(i)*ncm(i) 
           if(a(i).eq.badflag) then
              a(i)=(nsq*b(i)+hsq*c(i))/(nsq+hsq)
           endif
        enddo !i
c
        return
        end 
c
c
        subroutine replacemsng(a,b,imax,jmax,m,n,badflag)
c
c*********************************************************************
c
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       15 Dec 1999  Peter Stamus, NOAA/FSL
c          Housekeeping changes for use in LAPS.
c
c     Notes:
c
c*********************************************************************
c
        real a(m,n),b(m,n)    
        do j=1,jmax
        do i=1,imax
           if(b(i,j).eq.badflag) b(i,j)=a(i,j)
        enddo !i
        enddo !j
c
        return
        end 
