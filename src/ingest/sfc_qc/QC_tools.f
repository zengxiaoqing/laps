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
         sum=sum+ran1(ii)             
      enddo !l
      if(n.ne.1) then
       ffz=sum-float(n)/2.
      else
       ffz=sum
      endif
      ii=ii-1 
      if(ii.ge.0) ii=-ii-1
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
         Subroutine trimlist(ta,tda,dda,ffa,lata,lona,eleva,pstna,
     &     pmsla,alta,stna,providera,reptypea,indexa,maxstaa,
     &     m,badflag)
c This routine is designed to reduce the size of a station list by
c removing stations with all missing data (t, td, dd, ff, pmsl, alt)
c This should remove precip onlystations and those stations that
c have perpetually missing data. This will also remove stations
c that have all blanks for a name as this will screw up subr reorder.
         real ta(m),tda(m),dda(m),ffa(m),lata(m),lona(m),eleva(m)
         real pstna(m),pmsla(m),alta(m)
         integer indexa(m), iholdx, ix
         character stna(m)*5,holdnam*5, ch*1
         character providera(m)*11, holdprov*11
         character reptypea(m)*6,  holdrep*6
         holdnam='     '
c for all stations
          k=0
    3     k=k+1
    2     if(k.gt.maxstaa) return
          if(stna(k).eq.holdnam) then
           print*, 'Stn ',k,' has no name...removing data'
           go to 1
          endif 
          if(ta(k).ne.badflag) go to 3
          if(tda(k).ne.badflag) go to 3
          if(dda(k).ne.badflag) go to 3
          if(alta(k).ne.badflag) go to 3
          if(pmsla(k).ne.badflag) go to 3
           print*,stna(k),' has all missing data...removed from list'
    1      do l=k,maxstaa
            ta(l)=ta(l+1) 
            tda(l)=tda(l+1)
            dda(l)=dda(l+1)
            ffa(l)=ffa(l+1)
            pmsla(l)=pmsla(l+1)
            alta(l)=alta(l+1)
            stna(l)=stna(l+1)
            providera(l)=providera(l+1)
            indexa(l)=indexa(l+1)
            reptypea(l)=reptypea(l+1)
            lata(l)=lata(l+1)
            lona(l)=lona(l+1)
            eleva(l)=eleva(l+1)
            pstna(l)=pstna(l+1)
           enddo
           maxstaa=maxstaa-1
          go to 2
         end


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
         max=maxstaa
         if(maxstab.gt.max)max=maxstab
         if(maxstaa.eq.0) return
c
c first check the master list to find duplicates
c
         do k=1,maxstab
         holdnam = '     '
         do l=1,maxstab
            if(k.ne.l.and.stnb(k).eq.stnb(l)) then
               if(latb(k).ne.latb(l).or.lonb(k).ne.lonb(l)) then
                  print*, k,stnb(k),l,stnb(l)
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
c check the new list and see if that stn is in the old list
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
c make a place for the master list station in put it in the new ob list
c at the same position as in the master list.
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
        subroutine grosserr(t,tf,stn,maxsta,m,badflag,badthr,groserr
     &          ,qcstat,onoff)
c this routine performs a gross error check using a station average with
c +/- groserr added as limits. If onoff is 1 then it uses tf as an
c estimate with +/- groserr added as limits
        real t(m),tf(m)
        integer qcstat(m),onoff
        character stn(m)*5
        sum=0
        sum2=0
        cnt=0.
        do i=1,maxsta
         qcstat(i)=0.
         if (t(i).ne.badflag) then
          sum=t(i)+sum
          cnt=1.+cnt 
         endif
        enddo
        if(cnt.ne.0) then
         sum=sum/cnt
        else
         print*,'gross error check:entire data set is missing'
         return
        endif
        cnt=0.
        do i=1,maxsta
         if(t(i).ne.badflag) then
          sum2=(t(i)-sum)**2+sum2
          cnt=cnt+1.
         endif
        enddo
        sum2=sqrt(sum2/cnt)
        topg=sum+badthr*sum2
        botg=sum-badthr*sum2
        if(onoff.eq.0) then
         print*,'avg,stddev,low limit,high limit ',sum,sum2,botg,topg
        endif
        do i=1,maxsta
         if(t(i).eq.badflag) go to 1 
         if(tf(i).eq.badflag) go to 3
         if(onoff.eq.0) go to 3
           top=tf(i)+groserr
           bot=tf(i)-groserr
         if (t(i).lt.bot) go to 3
         if (t(i).gt.top) go to 3
         go to 1
   3     if(onoff.eq.1.and.tf(i).ne.badflag) then
          print*,t(i),'at stn ',i,'  ',stn(i),' fails gross kalman ',
     &     'check against ',tf(i),' ...setting to badflag'
           t(i)=badflag
           qcstat(i)=10
         else
          if (abs(t(i)-sum).lt.(badthr*sum2)) go to 1
            print*,t(i),'at stn ',i,'  ',stn(i),
     &    ' fails gross error check ...setting to badflag'
            t(i)=badflag
            qcstat(i)=10
         endif
   1    enddo
        return
        end

        subroutine fillqc(stn,t,lat,lon,elev,theta,
     &     lapse,maxsta,m,tmx,tmn,hcutoff,vcutoff,grosserr)
c
c*********************************************************************
c
c     Routine quality contols with buddy values and fills
c     array t with an interpolated value if missing
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
        real lapse,w(m)
        character stn(m)*5
c
        expsclr=1000.
c       hcutoff=specify         !dist in m for wt e**-1
c       vcutoff=specify         !ht(m) in vertical wt to go to e**-1
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
           w(i)=badflag
           if (t(i).ne.badflag) then
              cnt=cnt+1.
              tbar=tbar+t(i)
              ebar=ebar+elev(i)
           endif
        enddo !i
        if(cnt.eq.0.) then
         write(6,*) 'Fillqc: all missing data - returning'
         return
        endif
        tbar=tbar/cnt
        ebar=ebar/cnt
        nct=0
        do i=1,maxsta 
           c11=c1
           c22=c2
  10          sumwt=0.
              sum=0.
              do j = 1,maxsta
                 if(i.eq.j) go to 3
                 if(t(j).eq.badflag) go to 3
                 ang=(lat(i)+lat(j))*.5*rdpdg
                 dy=(lat(j)-lat(i))*reorpd
                 dx=(lon(j)-lon(i))*reorpd*cos(ang)
                 dz=theta(j)-theta(i)
                 if(theta(j).eq.badflag.or.theta(i).eq.badflag) then
                    dz=(elev(j)-elev(j))*thlapse   
                 endif
                 iii=int(c11*((dx*dx+dy*dy) +c22*dz*dz))+1
                 if(iii.gt.10000) iii=10000
                 sum=tab(iii)*(t(j)-(tbar+lapse*(elev(j)-ebar)))+sum     
                 sumwt=tab(iii)+sumwt
 3            enddo !on j
              if(sumwt .ne. 0.) then
                 sum=sum/sumwt+tbar+(elev(i)-ebar)*lapse
                 if(t(i).eq.badflag) then
                  w(i)=sum
                  nct=nct+1
                 else
                  if(abs(t(i)-sum).gt.grosserr) then
                   print*,'Sub fillqc replaces bad value of ',t(i),
     &              'at stn ',stn(i),' ',i,' with buddy value ',sum            
                   t(i)=badflag
                   w(i)=sum
                  endif
c                write(6,1000) i, stn(i), t(i)
 1000            format(1x,' Replaced missing ob at stn ',i5,1x,a5,
     &                  ' with ',F8.3) 
                 endif
              else
                 write(6,*) 'No fill possible for ',stn(i),',station ',i
                 write(6,*) '..          '
                 w(i)=t(i)    
              endif
        enddo !i
        do i=1,maxsta
         if(w(i).ne.badflag) t(i)=w(i)
        enddo
        write(6,*) 'Total number of missing obs replaced is ',nct
c
        return
        end
c
        subroutine fill(stn,t,lat,lon,elev,theta,
     &     lapse,maxsta,m,tmx,tmn,hcutoff,vcutoff)
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
        real lapse,w(m)
        character stn(m)*5
c
        expsclr=1000.
c       hcutoff=specify         !dist in m for wt e**-1
c       vcutoff=specify         !ht(m) in vertical wt to go to e**-1
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
           w(i)=badflag
           if (t(i).ne.badflag) then
              cnt=cnt+1.
              tbar=tbar+t(i)
              ebar=ebar+elev(i)
           endif
        enddo !i
        if(cnt.eq.0.) then
         write(6,*) 'Fill: All missing data - returning'
         return
        endif
        tbar=tbar/cnt
        ebar=ebar/cnt
        nct=0
        do i=1,maxsta 
           c11=c1
           c22=c2
  10       if(t(i).eq.badflag) then
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
                 iii=int(c11*((dx*dx+dy*dy) +c22*dz*dz))+1
                 if(iii.gt.10000) iii=10000
                 sum=tab(iii)*(t(j)-(tbar+lapse*(elev(j)-ebar)))+sum     
                 sumwt=tab(iii)+sumwt
 3            enddo !on j
              if(sumwt .ne. 0.) then
                 w(i)=sum/sumwt+tbar+(elev(i)-ebar)*lapse
                 nct=nct+1
c                write(6,1000) i, stn(i), w(i)
 1000            format(1x,' Replaced missing ob at stn ',i5,1x,a5,
     &                  ' with ',F8.3) 
              else
                 write(6,*) 'No fill possible for ',stn(i),',station ',i
                 write(6,*) '..setting parameter to missing           '
                 w(i)=t(i)   
              endif
           endif 
        enddo !i
        do i=1,maxsta
         if(w(i).ne.badflag) t(i)=w(i)
        enddo
        write(6,*) 'Total number of missing obs replaced is ',nct
c
        return
        end
c
c
        subroutine fillone(stn,t,lat,lon,elev,theta,
     &    lapse,maxsta,m,tmx,tmn,hcutoff,vcutoff)
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
c       hcutoff=specify         !dist in m for wt e**-1
c       vcutoff=specify         !ht(m) in vertical wt to go to e**-1
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
c          write(6,1001) slp,smn
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
           if(qcstat(i).eq.10) then
                    write(6,*) ' qcset:stn ',i,
     &                     ' failed gross error ck, replaced with ',
     &                        'Kalman estimate ',ta(i)
                    t(i) = ta(i)
           else
             if(t(i) .eq. badflag) then
               qcstat(i) = 1
               ncm(i) = ncm(i) + 1.
             else
               qcstat(i)=0.
               ncm(i) = 0.
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
     &        oberr,wk,wm,wo,wb,length,badthr,atime,icnt,nvar,n,i4time
     &            ,tor,tcr,tmr,tob,tcb,tmb,totr,totb,bmonster,ihr)
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
        real me,oe,ce,bmonster(m,24,nvar)
        real wm(m),wo(m),wb(m),wk(m)
        real tdt(m),b(m),tmb(nvar),tob(nvar),tcb(nvar),tmr(nvar)
        real totr(nvar),totb(nvar),tor(nvar),tcr(nvar),oberr2
        integer sca(2,2),scf(2,2),scb(2,2),sco(2,2),imax,qcstat(m)
        integer scat(2,2),scft(2,2),scbt(2,2),scot(2,2),on,off,ihr
c
c.....  Truth routine and missing ob replacement
c
c set period for averaging the bias error ~ 4 days for testing
        biaslen=4 !days for averaging bias
        icnt=icnt+1
        iiiii=-i4time+icnt
        on=1
        off=0
        if(totr(nvar).ne.0) then
        endif
        thresh=oberr*badthr 
        oberr2=oberr*oberr
        do i=1,imax
           me=wm(i)**2
           oe=wo(i)**2
           ce=wb(i)**2
c maximum likelyhood estimate of truth from all predictors
           tdt(i)= 
     &     (by(i)/ce+md(i)/me+y(i)/oe+dt(i)/oberr2)
     &            /(1./oe+1./me+1./ce+1./oberr2)
           fact=1.
c insert random observation error but skew it toward above truth estimate
           delt=oberr*ffz(iiiii,20)
           if(tdt(i).gt.dt(i)) then
             delt=abs(delt) 
           else
             delt=-abs(delt)
           endif
           tdt(i)=dt(i)+delt
        enddo !i

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
        sum4=0.
        do i=1,imax
           sumtot=sumtot+1.
           WR(i)=xt(i)-tdt(i)
           sumwr=sumwr+WR(i)
           stdwr=WR(i)*WR(i)+stdwr
           if(dt(i).ne.badflag) then 
              VR(i)=dt(i)-tdt(i)
           else 
              VR(i)=oberr*ffz(iiiii,20)
           endif 
           sumvr=sumvr+VR(i)
           stdvr=VR(i)*VR(i)+stdvr
           AR(i)=X(i)-TdT(i)
           sumar=sumar+AR(i)
           stdar=AR(i)*AR(i)+stdar
c
c  do kalmod errors
c length is the averaging period in cycles for the kalmods and 
c days for the bias correction. Make bias length 8 days
           wo(i)=(y(i)-tdt(i))**2+float(length-1)*wo(i)**2
           wo(i)=sqrt(wo(i)/float(length))
           wb(i)=(by(i)-tdt(i))**2+float(length-1)*wb(i)**2
           wb(i)=sqrt(wb(i)/float(length))
           wm(i)=(md(i)-tdt(i))**2+float(length-1)*wm(i)**2
           wm(i)=sqrt(wm(i)/float(length))
           wk(i)=(ar(i)-tdt(i))**2+float(length-1)*wk(i)**2
           wk(i)=sqrt(wk(i)/float(length))
           bmonster(i,ihr,n)=(y(i)-tdt(i)+(biaslen-1.)
     &                       *bmonster(i,ihr,n))/biaslen    
           sum1=sum1+wo(i)
           sum2=sum2+wb(i)
           sum3=sum3+wm(i)
           sum4=sum4+wk(i)
      write(6,9) 'KALMOD ERRORS - OB - BUD - MOD :',i,wo(i),wb(i),wm(i)
      write(6,9) 'KALMOD PERCENTS-OB - BUD - MOD :',i,.5*(wb(i)+wm(i))/
     &  (wo(i)+wb(i)+wm(i)),.5*(wo(i)+wm(i))/(wo(i)+wb(i)+wm(i)),
     &  .5*(wo(i)+wb(i))/(wo(i)+wb(i)+wm(i))
  9   format(1x,a32,i3,3f8.3)
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
        Write(6,*) 'OVERALL FRACTION OB,BUD,NWP '
        Write(6,*) .5*(sum2+sum3)/(sum1+sum2+sum3),.5*(sum1+sum3)/
     1       (sum1+sum2+sum3),.5*(sum1+sum2)/(sum1+sum2+sum3) 
        tmb(n)=tmb(n)+sumwr
        tob(n)=tob(n)+sumvr 
        tcb(n)=tcb(n)+sumar
        tmr(n)=tmr(n)+stdwr
        tor(n)=tor(n)+stdvr
        tcr(n)=tcr(n)+stdar
        totr(n)=totr(n)+sumtot
        totb(n)=totb(n)+sumtot 
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
        write(6,*) 'FX=XT ',(tmb(n)/totb(n)),sqrt(tmr(n)/totr(
     &                            n))
        write(6,*) 'OBSER ',(tob(n)/totb(n)),sqrt(tor(n)/totr(
     &                            n))
        write(6,*) 'KAL X ',(tcb(n)/totb(n)),sqrt(tcr(n)/totr(
     &                            n))
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
      
        subroutine checknan_pt(x,nan_flag)
c
c       Routine to check a single value for NaN's.
c
        real x
c
        nan_flag = 1
c
        if( nan( x ) .eq. 1) then
           print *,' ** ERROR. Found a NaN '
           nan_flag = -1
           return
        endif
c
        return
        end
c
c
        subroutine checknan_1d(x,ni,nan_flag)
c
c       Routine to check a real 1-d array for NaN's.
c
        integer ni
        real x(ni)
c
        nan_flag = 1
c
        do i=1,ni
           if( nan( x(i) ) .eq. 1) then
              print *,' ** ERROR. Found a NaN at ', i
              nan_flag = -1
              return
           endif
        enddo !i
c
        return
        end
c
c
        subroutine checknan_2d(x,ni,nj,nan_flag)
c
c       Routine to check a real 2-d array for NaN's.
c
        integer ni,nj
        real x(ni,nj)
c
        nan_flag = 1
c
        do j=1,nj
        do i=1,ni
           if( nan( x(i,j) ) .eq. 1) then
              print *,' ** ERROR. Found a NaN at ', i, j
              nan_flag = -1
              return
           endif
        enddo !i
        enddo !j
c
        return
        end
c

c
c
	subroutine windconvert(uwind,vwind,direction,speed,
     &                         ni,nj,badflag)
c
c======================================================================
c
c       Given wind components, calculate the corresponding speed and 
c       direction.  Hacked up from the windcnvrt_gm program.
c
c
c       Argument     I/O   Type       Description
c      --------	     ---   ----   -----------------------------------
c       UWind         I    R*4A    U-component of wind
c       VWind         I	   R*4A    V-component of wind
c       Direction     O    R*4A    Wind direction (meteoro. degrees)
c       Speed         O    R*4A    Wind speed (same units as input)
c       ni,nj         I    I       Grid dimensions
c       badflag       I    R*4     Bad flag value
c
c       Notes:
c       1.  If magnitude of UWind or VWind > 500, set the speed and 
c           direction set to the badflag value.
c
c       2.  Units are not changed in this routine.
c
c======================================================================
c
	real  uwind(ni,nj), vwind(ni,nj)
	real  direction(ni,nj), speed(ni,nj)
c
	do j=1,nj
	do i=1,ni
	   if(abs(uwind(i,j)).gt.500. .or. 
     &                           abs(vwind(i,j)).gt.500.) then
	      speed(i,j) = badflag
	      direction(i,j) = badflag
c
	   elseif(uwind(i,j).eq.0.0 .and. vwind(i,j).eq.0.0) then
	      speed(i,j) = 0.0
	      direction(i,j) = 0.0			!Undefined
c
	   else
	      speed(i,j) = 
     &          sqrt(uwind(i,j)*uwind(i,j) + vwind(i,j)*vwind(i,j))  !speed
	      direction(i,j) = 
     &          57.2957795 * (atan2(uwind(i,j),vwind(i,j))) + 180.   !dir
	   endif
	enddo !i
	enddo !j
c
	return
	end
c
c

	subroutine reduce_p(temp,dewp,pres,elev,lapse_t,
     &                          lapse_td,redpres,ref_lvl,badflag)
c
c
c================================================================================
c   This routine is designed to reduce the mesonet plains stations' pressure
C   reports to the elevation of the Boulder mesonet station, namely 1612 m.  The
C   hydrostatic equation is used to perform the reduction, with the assumption
C   that the mean virtual temperature in the layer between the station in ques-
C   tion and Boulder can approximated by the station virtual temperature.  This
C   is a sufficient approximation for the mesonet plains stations below about
C   7000 ft.  For the mountain and higher foothill stations (Estes Park,
C   Rollinsville, Ward, Squaw Mountain, and Elbert), a different technique needs
C   to be utilized in order to decently approximate the mean virtual temperature
C   in the layer between the station and Boulder. Ideas to do this are 1) use
C   the free air temperature from a Denver sounding, or 2) use the data from
C   each higher station to construct a vertical profile of the data and iterate
C   downward to Boulder.

C	D. Baker	 2 Sep 83  Original version.
C	J. Wakefield	 8 Jun 87  Changed ZBOU from 1609 to 1612 m.
c	P. Stamus 	27 Jul 88  Change ranges for good data tests.
c	P. Stamus	05 Dec 88  Added lapse rate for better computation
c	 				of mean virtual temps.
c			19 Jan 89  Fixed error with lapse rates (sheeze!)
c			19 Dec 89  Change reduction to 1500 m.
c			20 Jan 93  Version with variable reference level.
c                       25 Aug 97  Changes for dynamic LAPS.
c       P. Stamus       15 Nov 99  Change checks of incoming values so the
c                                    AF can run over Mt. Everest and Antarctica.
c
c       Notes:  This routine may or may not be giving reasonable results over
c               extreme areas (Tibet, Antarctica).  As noted above, there are
c               questions about use of the std lapse rate to do these reductions.
c               15 Nov 99
c
c================================================================================
c
	implicit none
	real lapse_t, lapse_td, temp, dewp, pres, elev, redpres, ref_lvl
	real badflag
	real gor, ctv
	parameter(gor=0.03414158,ctv=0.37803)
	real dz, dz2, t_mean, td_mean, td, e, tkel, tv, esw
!	DATA GOR,ZBOU,CTV/.03414158,1612.,.37803/
!	data gor,ctv/.03414158,.37803/
		!GOR= acceleration due to gravity divided by the dry air gas
		!     constant (9.8/287.04)
		!F2M= conversion from feet to meters
		! *** zbou is now the standard (reduction) level 12-19-89 ***
		!CTV= 1-EPS where EPS is the ratio of the molecular weight of
		!     water to that of dry air.


C** Check input values......good T, Td, & P needed to perform the reduction.
cc	if(dewp.gt.temp .or. pres.le.620. .or. pres.gt.1080. .or.
cc     &      temp.lt.-30. .or. temp.gt.120. .or. dewp.lt.-35. .or.
cc     &      dewp.gt.90.) then
	if(dewp.gt.temp .or. pres.le.275. .or. pres.gt.1150. .or.
     &      temp.lt.-130. .or. temp.gt.150. .or. dewp.lt.-135. .or.
     &      dewp.gt.100.) then
	   print *,' Warning. Bad input to reduce_p routine.'
	   redpres = badflag	!FLAG VALUE RETURNED FOR BAD INPUT
	   return
	endif

	dz= elev - ref_lvl	!thickness (m) between station & reference lvl
	dz2 = 0.5 * dz		! midway point in thickness (m)
	t_mean = temp - (lapse_t * dz2)	! temp at midpoint (F)
	td_mean = dewp - (lapse_td * dz2)	! dewpt at midpoint (F)
	TD= 0.55556 * (td_mean - 32.)		! convert F to C
	e= esw(td)		!saturation vapor pressure
	tkel= 0.55556 * (t_mean - 32.) + 273.15	! convert F to K
	tv= tkel/(1.-ctv*e/pres)	!virtual temperature

	redpres= pres*exp(gor*(dz/tv))	!corrected pressure
c
	return
	end


        Subroutine recover(stn,wx,index,b,br,bb,ncm,qcsta,
     &   tr,tf,tdr,tdf,ddr,ffr,
     &   uf,vf,pmslr,pmslf,badflag,grosserr,maxsta,m)

c       This subroutine checks a variable that has a suspected gross
c       error failure, by determining if the gross failure was due to
c       a local weather event, a discontinuity, or due to the reappearance
c       of the ob after a long string of missings. If criteria are met
c       the ob passes the gross error check

        real b(m),br(m),tr(m),tf(m),tdr(m),tdf(m),ddr(m),ffr(m)
        real uf(m),vf(m),pmslr(m),pmslf(m),ncm(m),bb(m)
        integer qcsta(m),index(m)
        character stn(m)*5,wx(m)*25,blank*25
 
        blank='                         '
   
        do i=1,maxsta
         ii=index(i)
         if(qcsta(i).eq.10) then 
          if(ncm(i).gt.0.) then! failed because it reappeared aft it was msng
           if(abs(b(i)-bb(i)).lt.grosserr) then ! it is cnstnt with buds
            write(6,*) 'QC reversal due to reappearance of ob at stn ',
     &      i,' ',stn(i),' after ',ncm(i),' missing cycles'
            write(6,*) '...reset to raw value of ', br(i)
            go to 1
           else
            go to 2 ! gross error failure stands
           endif
          endif
          if(wx(ii).ne.blank) then ! likely failed due to local weather
           write(6,*) 'QC reversal due to local weather at stn ' 
     &     ,i,' ',stn(i),'...reset to raw value of ', br(i),
     &     ' weather is ',wx(ii)
           go to 1
          endif 
          ic=0
          if(tr(i).ne.badflag.and.abs(tr(i)-tf(i)).gt.10.) ic=ic+1
          if(tdr(i).ne.badflag.and.abs(tdr(i)-tdf(i)).gt.10.) ic=ic+1

          spd = sqrt(uf(i)*uf(i) + vf(i)*vf(i) )   !speed
          dir = 57.2957795 * (atan2(uf(i),vf(i))) + 180.   !dir
       
          if(ffr(i).ne.badflag.and.abs(ffr(i)-spd).gt.15.) ic=ic+1
          if(ddr(i).ne.badflag.and.abs(ddr(i)-dir).gt.60.) ic=ic+1
          if(pmslr(i).ne.badflag.and.abs(pmslr(i)-pmslf(i)).gt.10.)
     &             ic=ic+1
          if(ic.ge.3) then
           write(6,*) 'QC reversal due to ',ic,' variables undergoing '
     &     , 'sig change at stn ',i,' ',stn(i),'...reset to ',br(i)
           go to 1
          endif
          go to 2
   1      b(i)=br(i)
          qcsta(i)=0
         endif
   2    enddo
        return
        end
 
           
