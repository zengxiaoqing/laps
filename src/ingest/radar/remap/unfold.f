       subroutine unfold(nogate,nobeam,dat,beam,velny,bmiss)
c***********************************************************************
c description      : To unfold Vr PPI data for LAPS system.
c                                                                       
c I/O parameters   :                                                    
c  I/O/W   name,            type,      description
c    I     nogate           integer    gate number in one ray.
c    I     nobeam           integer    total ray number of PPI sweep scan.
c   I/O  dat(nogate,nobeam) real array radial velocity in unit of m/s.
c    W datav(nobeam,nogate) real array working array.
c    I     velny            real       Nquist velocity in unit of m/s.
c    I     bmiss            real       missing value for radial velocity.
c                                      if( dat(i,j).eq.bmiss )then
c                                          dat(i,j) is missing data. 
c
c Date :
c   Feb. 25, 2004 (S.-M. Deng)
c***********************************************************************

       parameter( jb=420 )
       dimension dat(nogate,nobeam),beam(jb)
       real datav(nobeam,nogate)
       real ref(1000),ref1(1000),ref2(1000)
       real refb(1000)
       dimension iref(1000)
       dimension th(jb),vr(jb),cc(7)

       scale1=0.8
       scale2=0.7
       scale=1.0
       dummy=-999.
       rmiss=abs(bmiss)-1.
       do i=1,nobeam
       do j=1,nogate
          datav(i,j)=dat(j,i)
          if( abs(dat(j,i)).ge.rmiss )datav(i,j)=dummy
       enddo
       enddo

       call chkvel(nobeam,nogate,datav,beam,velny,dummy)

c -----initialize array data
       do k=1,nogate
          ref(k)=dummy
          ref1(k)=dummy
          ref2(k)=dummy
          iref(k)=0
       enddo

c ---- find the minimum radial velocity beam
       knref=40
 101   jzero=0
       fmean=0.5*velny
       kntol=0
       do j=1,nobeam
          knum=0
          dmean=0.
          do k=5,nogate
             if( datav(j,k).gt.dummy+1. )then
                dmean=dmean+abs(datav(j,k))
                knum=knum+1
             endif
          enddo
          if( knum.gt.knref )then
              kntol=kntol+1
              dmean=dmean/float(knum)
              if( dmean.lt.fmean )then
                  jzero=j
                  fmean=dmean
              endif
          endif
       enddo
       if( kntol.lt.30 .or. fmean.gt.10. )then
          knref=knref-10
          if( knref.ge.10 )go to 101
          if( knref.eq.0 .and. fmean.gt.10. )then
              knref=5
              go to 101
          endif
       endif
       if( jzero.eq.0 )go to 99        
       iodd=0
       ieven=0 
       do k=5,nogate
          if( datav(jzero,k).gt.dummy )then
              if( datav(jzero,k).lt.0. )then
                  ieven=ieven+1
              elseif( datav(jzero,k).gt.0. )then 
                  iodd=iodd+1
              endif
          endif
       enddo
       if( iodd.lt.ieven )then
           fmean=-1.*fmean
       elseif( iodd.eq.ieven )then
               fmean=0.
       endif 

c ---- unfold in this minimum velocity beam
       do k=5,nogate
          ref(k)=0.
       enddo
       ibeam=0
       call shearb1s(ibeam,fmean,jzero,scale,ref,iref,datav,nobeam
     +              ,nogate,dummy,velny)
       do k=5,nogate
          refb(k)=ref(k)
          ref(k)=datav(jzero,k)
       enddo

c ---- patch the reference velocity data in beam
       call bpatch(ref,nobeam,nogate,dummy,velny)
       call bpatch(ref,nobeam,nogate,dummy,velny)
       call bpatch(ref,nobeam,nogate,dummy,velny)
       call bpatch1(ref,nobeam,nogate,dummy,velny) 
       do k=5,nogate
          ref1(k)=ref(k)
       enddo
       jend=jzero-nobeam/2.
       if( jend.le.0 )jend=jend+nobeam

c ---- for clockwise direction
       do k=5,nogate
          ref2(k)=ref(k)
          iref(k)=0
       enddo
       if( jend.gt.jzero )then
           bmean=fmean 
           ibeam=0
           do j=jzero,jend
              call shearb1s(ibeam,bmean,j,scale,ref,iref
     1                     ,datav,nobeam,nogate,dummy,velny)
              call chgbeam(ibeam,bmean,j,ref,ref2,iref,datav
     1                    ,nobeam,nogate,dummy,velny)
           enddo
       else
           bmean=fmean 
           ibeam=0
           do j=jzero,nobeam
              call shearb1s(ibeam,bmean,j,scale,ref,iref
     1                     ,datav,nobeam,nogate,dummy,velny)
              call chgbeam(ibeam,bmean,j,ref,ref2,iref,datav
     1                    ,nobeam,nogate,dummy,velny)
           enddo
           do j=1,jend
              call shearb1s(ibeam,bmean,j,scale,ref,iref
     1                     ,datav,nobeam,nogate,dummy,velny)
              call chgbeam(ibeam,bmean,j,ref,ref2,iref,datav
     1                    ,nobeam,nogate,dummy,velny)
           enddo
       endif

c ---- for counter-clockwise direction
       do k=5,nogate
          ref(k)=ref1(k)
          ref2(k)=ref1(k)
          iref(k)=0
       enddo
       if( jend.gt.jzero )then
           bmean=fmean 
           ibeam=0
           do j=jzero,1,-1
              call shearb1s(ibeam,bmean,j,scale,ref,iref
     1                     ,datav,nobeam,nogate,dummy,velny)
              call chgbeam(ibeam,bmean,j,ref,ref2,iref,datav
     1                    ,nobeam,nogate,dummy,velny)
           enddo
           do j=nobeam,jend,-1 
              call shearb1s(ibeam,bmean,j,scale,ref,iref,datav
     1                     ,nobeam,nogate,dummy,velny)
              call chgbeam(ibeam,bmean,j,ref,ref2,iref,datav
     1                    ,nobeam,nogate,dummy,velny)
           enddo
       else
           bmean=fmean 
           ibeam=0
           do j=jzero,jend,-1
              call shearb1s(ibeam,bmean,j,scale,ref,iref
     1                     ,datav,nobeam,nogate,dummy,velny)
              call chgbeam(ibeam,bmean,j,ref,ref2,iref,datav
     1                    ,nobeam,nogate,dummy,velny)
           enddo
       endif
 99    continue

       pi=acos(-1.)
       do j=1,nogate
          n=0
          do i=1,nobeam
             if( datav(i,j).gt.-900. )then
                 n=n+1
                 th(n)=beam(i)
                 vr(n)=datav(i,j)
             endif
          enddo
          if( n.gt.180 )then
              call sov_coe7(n,th,vr,cc,ier)
              if( ier.ne.0 )then
                  do i=1,nobeam
                     if( datav(i,j).gt.-900. )then
                         dat(j,i)=datav(i,j)
                     endif
                  enddo
              else
                  do i=1,nobeam
                     if( datav(i,j).gt.-900. )then
                         dat(j,i)=datav(i,j)
                     else
                         if( abs(dat(j,i)).lt.rmiss )then
                             deg=pi*beam(i)/180.
                          val=cc(1)+cc(2)*cos(deg)+cc(3)*sin(deg)
     1                      +cc(4)*cos(2.*deg)+cc(5)*sin(2.*deg)
     2                      +cc(6)*cos(3.*deg)+cc(7)*sin(3.*deg)
                            datref=dat(j,i)
                           
                            call rfold(datref,val,velny
     1                                ,dummy,1.0)
                            if( datref.gt.-900. )then
                                dat(j,i)=datref
                            endif
                         endif
                     endif
                  enddo
              endif
          else
              do i=1,nobeam
                 if( datav(i,j).gt.-900. )then
                     dat(j,i)=datav(i,j)
                 endif
              enddo
          endif
       enddo

       velny2=2.0*velny
       veck=1.6*velny
       do 30 j=1,nobeam
       do 10 i=2,nogate-1,1
         if((abs(dat(i,j)).gt.rmiss).or.(datav(j,i-1).lt.-900.))go to 10
          if( dat(i,j)*dat(i-1,j).gt.0. )go to 10
          ck=dat(i,j)-dat(i-1,j)
          if( (abs(ck).gt.veck).and.(abs(dat(i-1,j)).gt.veck) )then
              if( ck.gt.0. )then
                  dat(i,j)=dat(i,j)-velny2
              else
                  dat(i,j)=dat(i,j)+velny2
              endif 
              datav(j,i)=dat(i,j)
          endif
 10    continue
       do 20 i=nogate-1,2,-1
         if((abs(dat(i,j)).gt.rmiss).or.(datav(j,i+1).lt.-900.))go to 20
          if( dat(i,j)*dat(i+1,j).gt.0. )go to 20
          ck=dat(i,j)-dat(i+1,j) 
          if( (abs(ck).gt.veck).and.(abs(dat(i+1,j)).gt.veck) )then
              if( ck.gt.0. )then
                  dat(i,j)=dat(i,j)-velny2
              else
                  dat(i,j)=dat(i,j)+velny2
              endif
              datav(j,i)=dat(i,j)
          endif
 20    continue
 30    continue

       do 60 i=2,nogate-1,1
       do 40 j=2,nobeam,1
         if((abs(dat(i,j)).gt.rmiss).or.(datav(j-1,i).lt.-900.))go to 40
          if( dat(i,j)*dat(i,j-1).gt.0. )go to 40
          ck=dat(i,j)-dat(i,j-1)
          if( (abs(ck).gt.veck).and.(abs(dat(i,j-1)).gt.veck) )then
              if( ck.gt.0. )then
                  dat(i,j)=dat(i,j)-velny2
              else
                  dat(i,j)=dat(i,j)+velny2
              endif 
              datav(j,i)=dat(i,j)
          endif
 40    continue
       do 50 j=nobeam-1,1,-1
         if((abs(dat(i,j)).gt.rmiss).or.(datav(j+1,i).lt.-900.))go to 50
          if( dat(i,j)*dat(i,j+1).gt.0. )go to 50
          ck=dat(i,j)-dat(i,j+1)
          if( (abs(ck).gt.veck).and.(abs(dat(i,j+1)).gt.veck) )then
              if( ck.gt.0. )then
                  dat(i,j)=dat(i,j)-velny2
              else
                  dat(i,j)=dat(i,j)+velny2
              endif 
              datav(j,i)=dat(i,j)
          endif
 50    continue
 60    continue

       do 90 j=1,nobeam
       do 70 i=2,nogate-1,1
         if((abs(dat(i,j)).gt.rmiss).or.(datav(j,i-1).lt.-900.))go to 70
          if( dat(i,j)*dat(i-1,j).gt.0. )go to 70
          ck=dat(i,j)-dat(i-1,j)
          if( (abs(ck).gt.veck).and.(abs(dat(i-1,j)).gt.veck) )then
              if( ck.gt.0. )then
                  dat(i,j)=dat(i,j)-velny2
              else
                  dat(i,j)=dat(i,j)+velny2
              endif 
              datav(j,i)=dat(i,j)
          endif
 70    continue
       do 80 i=nogate-1,2,-1
         if((abs(dat(i,j)).gt.rmiss).or.(datav(j,i+1).lt.-900.))go to 80
          if( dat(i,j)*dat(i-1,j).gt.0. )go to 80
          ck=dat(i,j)-dat(i+1,j) 
          if( (abs(ck).gt.veck).and.(abs(dat(i+1,j)).gt.veck) )then
              if( ck.gt.0. )then
                  dat(i,j)=dat(i,j)-velny2
              else
                  dat(i,j)=dat(i,j)+velny2
              endif
              datav(j,i)=dat(i,j)
          endif
 80    continue
 90    continue

       return
       end

       subroutine chkvel(nobeam,nogate,datav,beam,velny,dummy)
       parameter( jb=420 )
       real datav(nobeam,nogate),beam(jb)

       do j=1,nobeam
       do k=5,nogate
          if( datav(j,k).gt.dummy )then
              fmean=0.
              idigit=0
              itotal=0
              js=j-3
              je=j+3
              do j1=max(j-3,1),min(j+3,nobeam)
                 difbeam=abs(beam(j)-beam(j1))
                 if(difbeam .lt. 3. .or.
     +              difbeam .gt. 357.) then 
                    do k1=max(k-3,5),min(k+3,nogate)
                       itotal=itotal+1
                       if( datav(j1,k1) .gt. dummy) then
                           fmean=fmean+datav(j1,k1)
                           idigit=idigit+1
                       endif
                    enddo
                 endif
              enddo
              if( js.le.0 )then
                  do j1=j-3+nobeam,nobeam
                     difbeam=abs(beam(j)-beam(j1))
                     if( difbeam.lt.3. .or.
     +                   difbeam .gt. 357.)then
                         do k1=max(k-3,5),min(k+3,nogate)
                            itotal=itotal+1
                            if( datav(j1,k1).gt.dummy )then
                                fmean=fmean+datav(j1,k1)
                               idigit=idigit+1
                            endif
                         enddo
                     endif
                  enddo
              endif
              if( je.gt.nobeam )then
                  do j1=1,je-nobeam
                     difbeam=abs(beam(j)-beam(j1))
                     if( difbeam.lt.3. .or.
     +                   difbeam .gt. 357.) then
                         do k1=max(k-3,5),min(k+3,nogate)
                            if( datav(j1,k1).gt.dummy )then
                                fmean=fmean+datav(j1,k1)
                                idigit=idigit+1
                            endif
                         enddo
                     endif
                  enddo
              endif
              if( itotal.gt.10 )then
                  if( idigit.le.itotal/5+1 )then
                      datav(j,k)=dummy
                      go to 801
                  endif
              endif
              fmean=fmean/float(idigit)
              iunf=0
 3            diff=abs(datav(j,k)-fmean)
              if( diff.gt.(1.0*velny) )then
                  iunf=iunf+1
                  if( iunf.gt.3 )then
                      datav(j,k)=dummy
                      go to 801
                  endif
                  if( datav(j,k).gt.fmean )then
                      datav(j,k)=datav(j,k)-2.*velny
                  else
                      datav(j,k)=datav(j,k)+2.*velny
                  endif
                  go to 3
              endif
          endif

 801      continue
       enddo
       enddo
       return
       end

       subroutine chgbeam(ibeam,bmean,j,ref,ref2,iref,datav
     1                   ,nobeam,nogate,dummy,velny)
       parameter( jb=420 )
       real datav(nobeam,nogate),beam(jb)
       real ref(1000),ref2(1000)
       dimension iref(1000)

       do k=5,nogate
          ref(k)=datav(j,k)
       enddo

c ---- patch the reference velocity data in beam
       call bpatch(ref,nobeam,nogate,dummy,velny)
       call bpatch(ref,nobeam,nogate,dummy,velny)
       call bpatch(ref,nobeam,nogate,dummy,velny)
       call bpatch1(ref,nobeam,nogate,dummy,velny)
       kk2=0
       tmean=0.
       do k=5,nogate
          if( datav(j,k).gt.dummy+1. )then
              kk2=kk2+1
              tmean=tmean+datav(j,k)
          endif
       enddo
       if( kk2 .ge. 10 )then
           bmean=tmean/float(kk2)
           ibeam=0
       else
           ibeam=ibeam+1
           if( ibeam.gt.50 )bmean=0.
       endif

       do k=5,nogate
          if( ref(k).le.dummy+1. )then
            if( iref(k).le.30 .and. ref2(k).gt.dummy+1. )then
              if( ref(k-1).gt.dummy .or. ref(k+1).gt.dummy )then
                  if( abs(ref(k-1)-ref2(k)).le.10. )then
                      ref(k)=ref2(k)
                      iref(k)=iref(k)+1
                  elseif( abs(ref(k+1)-ref2(k)).le.10. )then
                          ref(k)=ref2(k)
                          iref(k)=iref(k)+1
                      else
                          ref(k)=dummy
                          ref2(k)=dummy
                          iref(k)=31
                  endif
              else
                  ref(k)=ref2(k)
                  iref(k)=iref(k)+1
              endif
            else
              ref(k)=dummy
              ref2(k)=dummy
              iref(k)=31
            endif
          else
            ref2(k)=ref(k)
            iref(k)=0
          endif
       enddo
       return
       end

       subroutine shearb1s(ibeam,bmean,j,scale,ref,iref,datav
     +                    ,nobeam,nogate,dummy,velny)
       parameter( jb=420 )
       real datav(nobeam,nogate),beam(jb)
       real ref(1000)
       dimension iref(1000)
       dimension ka(4),kb(4),chose(4)
       dimension rev0(1000),rev1(1000)

       do k=1,nogate
          rev0(k)=datav(j,k)
       enddo
       kreverse=0
       ireverse=0 
       kfirst=5
       kend=nogate
       kstep=1
 101   continue
       kkk=0
       ifirst=0
       do k=5,nogate
          if( datav(j,k).gt.dummy )then
              if( ref(k).gt.dummy )then
                  call rfold(datav(j,k),ref(k),velny
     1	                    ,dummy,1.0*scale)
              endif
          endif
       enddo
       do k=kfirst,kend,kstep
          if( (datav(j,k).gt.dummy).and.(ref(k).gt.dummy)
     +        .and.(abs(datav(j,k)-ref(k)).lt.10.).and.
     +        (iref(k).le.10) )then 
              ifirst=k
              kkk=kkk+1
              if(kkk .eq. 15) go to 801
          endif
       enddo
       if( kkk.lt.5 )then
 805       continue
           kk3=3
           ifirst1=0
           do k=10,nogate-5
              if( datav(j,k).gt.dummy )then 
                  kk2=0
                  do kk=k-5,k+5
                     if( datav(j,kk).gt.dummy )then
                         kk2=kk2+1
                     endif
                  enddo
                  if( kk2.gt.kk3 )then
                      kk3=kk2
                      ifirst1=k
                      if( kk3.gt.10 )go to 802
                  endif
              endif
           enddo
           if( ifirst1.eq.0 )then
               do k=5,nogate
                  datav(j,k)=dummy
               enddo
               return
           endif
 802       continue
           do k=1,4
              ka(k)=0
              kb(k)=0
              chose(k)=dummy
           enddo
           ka(1)=1
           chose(1)=datav(j,ifirst1)
           kb(1)=ifirst1
           km=1
           do k=ifirst1-5,ifirst1+5
              if( datav(j,k).gt.dummy )then
                  do kk=1,km
                if(abs(datav(j,k)-chose(kk)).le.10.)go to 803
                  enddo
                  km=km+1
                  kb(km)=k
                  chose(km)=datav(j,k)
              endif
 803          continue
           enddo
           do k=ifirst1-5,ifirst1+5 
              if( datav(j,k).gt.dummy )then 
                  if( ref(k).gt.dummy .and. 
     1                abs(datav(j,k)-ref(k)).lt.10. )then
                      ifirst=k
                      go to 801
                  endif
                  do kk=1,km
                 if(abs(datav(j,k)-chose(kk)).le.10.)then
                        ka(kk)=ka(kk)+1
                     endif
                  enddo
              endif
           enddo

           ifirst=kb(1)
           kt=ka(1)
           if( km.gt.1 )then
               do k=2,km
                  if( ka(k).gt.kt )then
                      kt=ka(k)
                      ifirst=kb(k)
                  elseif( ka(k).eq.kt )then
                      if( ibeam.lt.50 )then
                          if( abs(datav(j,kb(k))-bmean).lt.
     +                       abs(datav(j,ifirst)-bmean) )then
                             kt=ka(k)
                             ifirst=kb(k)
                          endif
                      else
                          if( abs(datav(j,kb(k))).lt. 
     +                        abs(datav(j,ifirst)) )then
                              kt=ka(k)
                              ifirst=kb(k)
                          endif
                      endif
                  endif
               enddo
           endif

           if( abs(datav(j,ifirst)).gt.10. )then
               if( abs(datav(j,ifirst)-bmean).gt.velny .and.
     +             abs(bmean).ge.10. )then
                   call rfold(datav(j,ifirst),bmean
     +                       ,velny,dummy,1.0*scale)

                   if( abs(datav(j,ifirst)).gt.10. .and.
     +                 abs(datav(j,ifirst)-bmean ).gt. 
     +         abs(datav(j,ifirst)) )datav(j,ifirst)=dummy
                   elseif( abs(datav(j,ifirst)-bmean).gt.velny 
     1         .and.  abs(bmean).lt.10. .and. ibeam.lt.10 )then
                          call rfold(datav(j,ifirst),bmean
     1                              ,velny,dummy,1.0*scale)
                          if( abs(datav(j,ifirst)).gt.10. .and.
     +       abs(datav(j,ifirst)-bmean) .gt.
     +       abs(datav(j,ifirst))) datav(j,ifirst)=dummy  
                   endif
                   if( datav(j,ifirst).le.dummy+1. )go to 805 
               endif
           endif
         
 801       continue

           idigi=ifirst 
           iend=5
           istep=-1
           call shearb2s(kreverse,j,ifirst,iend,istep,idigi,scale
     1            ,bmean,ref,iref,nobeam,nogate,dummy,datav,velny)
           idigi=ifirst
           iend=nogate
           istep=1
           call shearb2s(kreverse,j,ifirst,iend,istep,idigi,scale
     1            ,bmean,ref,iref,nobeam,nogate,dummy,datav,velny)
           if( ireverse.eq.0 )then
               if( kreverse.eq.0 )return
               do k=1,nogate
                  rev1(k)=datav(j,k)
               enddo
               do k=1,nogate
                  datav(j,k)=rev0(k)
               enddo
               kfirst=nogate
               kend=5
               kstep=-1
               ireverse=1
               go to 101
           else
               do k=1,nogate
                  rev0(k)=datav(j,k)
               enddo

               do k=1,nogate
                  if( abs(rev0(k)-rev1(k)).le.1. )then
                     if( abs(rev0(k)-ref(k)).gt.velny .and.
     +               abs(rev0(k)) .gt. 2.*velny .and. ref(k)
     +               .gt. dummy .and. iref(k) .le. 10) then
                         datav(j,k)=dummy
                     else
                         datav(j,k)=rev0(k)
                     endif
                  else
                    if( ref(k).gt.dummy+1. .and. 
     1                  iref(k).le.3 )then
                      if( abs(rev0(k)-ref(k)).le.10.)then
                        if( abs(rev0(k)) .le. 2.*velny )then
                              datav(j,k)=rev0(k)
                        else 
                              datav(j,k)=dummy
                        endif
                      elseif( abs(rev1(k)-ref(k)).le.10.)then
                           if( abs(rev1(k)).le.2.*velny )then
                               datav(j,k)=rev1(k)
                           else 
                               datav(j,k)=dummy
                           endif
                        else
                           datav(j,k)=dummy
                      endif
                    else
                      datav(j,k)=dummy
                    endif
                  endif
               enddo
               ireverse=0
           endif
  
       return
       end

       subroutine shearb2s(kreverse,j,ifirst,iend,istep,idigi,scale
     1           ,bmean,ref,iref,nobeam,nogate,dummy,datav,velny)
       parameter( jb=420 )
       real datav(nobeam,nogate),beam(jb)
       dimension ref(1000)
       dimension iref(1000)

       mfold=0
       do k=ifirst+istep,iend,istep
         if( datav(j,k).gt.dummy )then
           diff1=abs(datav(j,k)-datav(j,idigi))
           if( diff1.gt.10. )then 
             if( ref(k).gt.dummy )then 
               if( abs(ref(k)-datav(j,k)).ge.10. )then
                 if( diff1.gt.velny )then
                     call rfold(datav(j,k),datav(j,idigi)
     +                         ,velny,dummy,1.0*scale)  
                     if(abs(datav(j,k)-datav(j,idigi)).ge.10.
     +              .and. abs(ref(k)-datav(j,k)).ge.10. )then
                        datav(j,k)=dummy
                     endif
                     if( abs(datav(j,k)-ref(k)) .gt. 
     1                   1.0*velny )then
                         datav(j,k)=dummy
                     endif
                 else
                     datav(j,k)=dummy
                 endif
               else
                 if( diff1.gt.velny )then
                     if( abs(k-idigi).le.10 )then
                         call rfold(datav(j,k),datav(j,idigi)
     +                             ,velny,dummy,1.0*scale)
                      if(abs(datav(j,k)-datav(j,idigi)).ge.10.
     +               .and. abs(ref(k)-datav(j,k)).ge.10. )then
                          datav(j,k)=dummy
                      endif
                      if( abs(datav(j,k)-ref(k)).gt.1.0*velny
     +                   .and. iref(k).le.10 )then
                          datav(j,k)=dummy
                          kreverse=1
                      endif
                     else
                        if( abs(datav(j,k)-bmean) .gt. 
     1                      abs(datav(j,k)) .or. 
     2                   abs(datav(j,k)-bmean).gt.velny )then  
                          call rfold(datav(j,k),datav(j,idigi)
     +                              ,velny,dummy,1.0*scale)
                          if( abs(datav(j,k)-datav(j,idigi)) 
     1                        .ge. 10.)then
                              datav(j,k)=dummy
                          endif
                        else
                            mfold=1
                        endif
                     endif
                 else
                     if( abs(k-idigi).le.5 )then
                         datav(j,k)=dummy
                     endif
                 endif
               endif

             else

               if( diff1.gt.velny )then
                   call rfold(datav(j,k),datav(j,idigi)
     +                       ,velny,dummy,1.0*scale)
                   if( abs(datav(j,k)-datav(j,idigi)).ge.10. 
     +                 .and. abs(k-idigi) .le. 5) then
                       datav(j,k)=dummy
                   elseif(abs(k-idigi) .gt. 50 .and. 
     +                abs(datav(j,k)-datav(j,idigi)).gt.15.
     +              .and. abs(datav(j,k)-bmean).gt.15. .and.
     +              abs(datav(j,k)) .gt. abs(datav(j,idigi)))
     +              then
                         datav(j,k)=dummy
                       elseif( abs(datav(j,k)).gt.velny .and. 
     +                         datav(j,k).gt.dummy )then
                           if((abs(datav(j,k)-bmean).gt.velny) 
     1                      .and. (abs(datav(j,k))-10. .gt. 
     2                           abs(datav(j,idigi))) )then
                                 call rfold(datav(j,k),bmean
     1                                 ,velny,dummy,1.0*scale)
                               if(abs(datav(j,k)).ge.10. .and.
     +                          datav(j,k)*bmean .lt. 0.) then
                                  datav(j,k)=dummy
               endif
               if(abs(datav(j,k)-datav(j,idigi)) 
     +           .ge. 15.) then
                datav(j,k)=dummy
               endif
              endif
             endif
            else
             if(abs(k-idigi) .le. 5) then
              datav(j,k)=dummy
             else
              kright=0
              kfault=0
              if(istep .eq. 1) then
               kkend=min(k+10,nogate)
              else
               kkend=max(k-10,5)
              endif
              do kk=k+istep,kkend,istep
               if(datav(j,kk) .gt. dummy) then
                if(abs(datav(j,kk)-datav(j,idigi)) 
     +             .lt. 10.) then 
                 kfault=kfault+1
                elseif(abs(datav(j,kk)-datav(j,k))
     +             .lt. 10.) then
                 kright=kright+1
                endif
               endif
              enddo

              if( (kright+kfault).lt.2 )then
               datav(j,k)=dummy
              else
               if( (kfault-kright).gt.2 )then
                datav(j,k)=dummy
               else
                if( abs(datav(j,k)).ge.10. )then
                 if( abs(datav(j,k)-bmean) .gt. 
     +              abs(datav(j,idigi)-bmean) .and.
     +              abs(datav(j,k)-bmean) .ge. 10.) then
                  datav(j,k)=dummy
                 endif
                endif
               endif
              endif
             endif
            endif
           endif
          endif
          if( datav(j,k).gt.dummy )idigi=k
         endif
       enddo
       if( mfold.eq.0 )go to 103 
       idigi=ifirst
       k=ifirst+istep

 101   continue 
       if( k.eq.iend)go to 103
       if( datav(j,k).gt.dummy ) then
         if( (abs(datav(j,k)-datav(j,idigi)) 
     +      .gt. velny )
     +      .and. (abs(k-idigi) .le. 50)) then
          idigi1=k
          do kk=k+istep,iend,istep
           if( datav(j,kk).gt.dummy )then
            if(abs(datav(j,kk)-datav(j,idigi1))
     +          .lt. velny) then
             if(abs(datav(j,kk)-datav(j,idigi))
     +          .lt. velny) then
              k=kk+istep
              idigi=kk
              go to 101
             endif 
             idigi1=kk
             go to 102
            else
             if((abs(datav(j,kk)-datav(j,idigi))
     +          .lt. velny) .or. (datav(j,k)-
     +          datav(j,idigi))*(datav(j,kk)-
     +          datav(j,idigi1)) .lt. 0.) then
              do kkk=k,kk-istep,istep
               datav(j,kkk)=dummy
              enddo
             endif
             k=kk+istep
             idigi=kk
             go to 101 
            endif
           endif
 102       continue
          enddo 
         endif
         idigi=k
       endif
       k=k+istep
       go to 101
 103   continue
       return
       end

       subroutine bpatch1(ref,nsamp,dummy,velny)
       dimension ref(1000)

       do k=6,nsamp-1
          if( ref(k).le.dummy+1. .and. ref(k-1).gt.dummy )then
              k1=k-1
              do kk=k+1,nsamp
                 if( ref(kk).gt.dummy )then
                     k2=kk
                     go to 201
                 endif
              enddo
 201          continue
              if( abs(ref(k1)-ref(k2)).le.10. .and.
     +            k2-k1.le.80 )then
                  do kk=k1+1,k2-1
                     ref(kk)=ref(k1)+float(kk-k1)*
     +                       (ref(k2)-ref(k1))/float(k2-k1)
                  enddo
              endif
          endif
       enddo

       do k=6,nsamp
          if( ref(k-1).gt.dummy .and. ref(k).gt.dummy 
     +       .and. abs(ref(k-1)-ref(k)).gt.1.5*velny )then
              do kk=k+1,min(k+20,nsamp)
                 if( ref(kk).lt.dummy+1 )go to 202 
                 if( abs(ref(kk)-ref(k-1)).lt.velny )go to 203
              enddo
              go to 202
 203          do kkk=k,kk-1
                 ref(kkk)=ref(k-1)+float(kkk-k+1)*
     +                    (ref(kk)-ref(k-1))/float(kk-k+1)
              enddo
 202          continue
          endif
       enddo
       return

       end 
           
       subroutine bpatch(ref,nsamp,dummy,velny)
       dimension ref(1000)

       do k=6,nsamp-1
          if( ref(k) .le. dummy+1. )then
              kk1=0
              kk2=0
              do kk=k,max(k-5,5),-1
                 if( ref(kk).gt.dummy+1. )then
                     kk1=k-kk
                     go to 11
                 endif
              enddo
 11           continue
              do kk=k,min(k+5,nsamp)
                 if( ref(kk).gt.dummy+1. )then
                     kk2=kk-k
                     go to 12
                 endif
              enddo
 12           continue
              if( kk1.gt.0 .and. kk2.gt.0 )then
                  call rfold(ref(k+kk2),ref(k-kk1),velny
     1                      ,dummy,1.4)
                  if( ref(k+kk2).gt.dummy+1. )then
                      ref(k)=(float(kk2)*ref(k-kk1)+float(kk1)*
     +                ref(k+kk2))/float(kk1+kk2)
                  endif
              endif
          endif
       enddo
       return
       end

       subroutine rfold(datav,ref,velny,dummy,scale)
       datav1=datav
       iunf=0
 3     diff=abs(datav-ref)
       if( diff.gt.(scale*velny) )then
           iunf=iunf+1
           if( iunf.gt.3 )then
               datav=dummy
               return
           endif
           if( datav.gt.ref )then
               datav=datav-2.*velny
           else
               datav=datav+2.*velny
           endif
           go to 3
       endif
       return
       end

      subroutine sov_coe7(n,th,vr,cc,ier)
c***********************************************************************
c description      : To compute Fourier coefficient of Vr by curve-     
c                    fitting method.                                    
c                                                                       
c I/O parameters   :                                                    
c  input:  name,      type,       description                           
c          n          integer     the total number.                     
c          th(n)      real array  the input angle in unit of degrees.   
c          vr(n)      real array  the input Vr value.                   
c  output: name,      type,       description                           
c          cc(7)      real array  the Fourier coefficient of Vr.        
c          ier        integer     the error message.                    
c                                 =0, no error message.                 
c                                 =1, cannot compute the coefficient.   
c                                                                       
c called fun./sub. : sov_eq7                                            
c                                                                       
c note:                                                                 
c         Vr(th)=cc(1)+cc(2)*cos(th)+cc(3)*sin(th)+                     
c                     +cc(4)*cos(2*th)+cc(5)*sin(2*th)+                 
c                     +cc(6)*cos(3*th)+cc(7)*sin(3*th)                  
c                                                                       
c auther           : Deng Shiung Ming                                   
c                    - Institude for Information Industry               
c                                                                       
c create date      : Aug 24, 1995                                       
c***********************************************************************
      parameter( max=420 )
      dimension th(n),vr(n),cc(7),a(7,7),y(7)
      dimension cs1(max),sn1(max),cs2(max),sn2(max),cs3(max),sn3(max)

      if( n.gt.max )then
          ier=1
          print*,'Input total number too large'
          print*,'total number: ',n
          print*,'the max. number: ',max
          print*,'Please modified input number or max. parameter'
          return
      endif
c-----------------------------------------------------------------------
c* To compute the 7-D maxtrix.

      pi=acos(-1.)
      do i=1,n
         thi=pi*th(i)/180.
         cs1(i)=cos(thi)
         sn1(i)=sin(thi)
         cs2(i)=cos(2.*thi)
         sn2(i)=sin(2.*thi)
         cs3(i)=cos(3.*thi)
         sn3(i)=sin(3.*thi)
      enddo

      do i=1,7
      do j=1,7
         a(i,j)=0.
      enddo
      enddo

      a(1,1)=n
      do k=1,n
         a(2,1)=a(2,1)+cs1(k)
         a(3,1)=a(3,1)+sn1(k)
         a(4,1)=a(4,1)+cs2(k)
         a(5,1)=a(5,1)+sn2(k)
         a(6,1)=a(6,1)+cs3(k)
         a(7,1)=a(7,1)+sn3(k)
         a(2,2)=a(2,2)+cs1(k)**2
         a(3,2)=a(3,2)+sn1(k)*cs1(k)
         a(4,2)=a(4,2)+cs2(k)*cs1(k)
         a(5,2)=a(5,2)+sn2(k)*cs1(k)
         a(6,2)=a(6,2)+cs3(k)*cs1(k)
         a(7,2)=a(7,2)+sn3(k)*cs1(k)
         a(3,3)=a(3,3)+sn1(k)**2
         a(4,3)=a(4,3)+cs2(k)*sn1(k)
         a(5,3)=a(5,3)+sn2(k)*sn1(k)
         a(6,3)=a(6,3)+cs3(k)*sn1(k)
         a(7,3)=a(7,3)+sn3(k)*sn1(k)
         a(4,4)=a(4,4)+cs2(k)**2
         a(5,4)=a(5,4)+sn2(k)*cs2(k)
         a(6,4)=a(6,4)+cs3(k)*cs2(k)
         a(7,4)=a(7,4)+sn3(k)*cs2(k)
         a(5,5)=a(5,5)+sn2(k)**2
         a(6,5)=a(6,5)+cs3(k)*sn2(k)
         a(7,5)=a(7,5)+sn3(k)*sn2(k)
         a(6,6)=a(6,6)+cs3(k)**2
         a(7,6)=a(7,6)+sn3(k)*cs3(k)
         a(7,7)=a(7,7)+sn3(k)**2
      enddo
      a(1,2)=a(2,1)
      a(1,3)=a(3,1)
      a(1,4)=a(4,1)
      a(1,5)=a(5,1)
      a(1,6)=a(6,1)
      a(1,7)=a(7,1)
      a(2,3)=a(3,2) 
      a(2,4)=a(4,2) 
      a(2,5)=a(5,2) 
      a(2,6)=a(6,2) 
      a(2,7)=a(7,2) 
      a(3,4)=a(4,3)
      a(3,5)=a(5,3)
      a(3,6)=a(6,3)
      a(3,7)=a(7,3)
      a(4,5)=a(5,4)
      a(4,6)=a(6,4)
      a(4,7)=a(7,4)
      a(5,6)=a(6,5)
      a(5,7)=a(7,5)
      a(6,7)=a(7,6)

c-----------------------------------------------------------------------
c* To compute the Y.

      do k=1,7
         y(k)=0.
      enddo
      do i=1,n
         y(1)=y(1)+vr(i)
         y(2)=y(2)+vr(i)*cs1(i)
         y(3)=y(3)+vr(i)*sn1(i)
         y(4)=y(4)+vr(i)*cs2(i)
         y(5)=y(5)+vr(i)*sn2(i)
         y(6)=y(6)+vr(i)*cs3(i)
         y(7)=y(7)+vr(i)*sn3(i)
      enddo


c-----------------------------------------------------------------------
c* To compute the Fourier coefficient of Vr.

      call sov_eq7(a,y,cc,ier)

      return
      end

      subroutine sov_eq7(a,y,x,ier)
c***********************************************************************
c description      : To sove 7-linear equation A*X=Y, by  X=(A**-1)*Y.  
c                                                                       
c I/O parameters   :                                                    
c  input:  name,      type,       description                           
c          a(7,7)     real array  the coefficient matrix.               
c          y(7)       real array  the input matrix.                     
c  output: name,      type,       description                           
c          x(7)       real array  the solved matrix.                    
c          ier        integer     the error message.                    
c                                 =0, no error message.                 
c                                 =1, cannot sove the matrix equation.  
c                                                                       
c called fun./sub. : invers7                                            
c                                                                       
c auther           : Deng Shiung Ming                                   
c                    - Institude for Information Industry               
c                                                                       
c create date      : Aug 24, 1995                                       
c***********************************************************************
      dimension a(7,7),x(7),y(7),b(7,7)

c-----------------------------------------------------------------------
c* To compute the inverse of matrix A.

      call invers7(a,b,ier)
      if( ier.eq.1 )return

      do k=1,7
         x(k)=0.
         do i=1,7
            x(k)=x(k)+b(i,k)*y(i)
         enddo
      enddo

      return
      end

      subroutine invers7(a,b,ier)
c***********************************************************************
c description      : To compute the inverse 7-D matrix.                 
c                                                                       
c I/O parameters   :                                                    
c  input:  name,      type,       description                           
c          a(7,7)     real array  7-D matrix.                           
c  output: name,      type,       description                           
c          b(7,7)     real array  inverse 7-D matrix of a.              
c          ier        integer     the error message.                    
c                                 =0, no error message.                 
c                                 =1, cannot compute the inverse matrix 
c                                                                       
c called fun./sub. : det6,     det7                                     
c                                                                       
c auther           : Deng Shiung Ming                                   
c                    - Institude for Information Industry               
c                                                                       
c create date      : Aug 24, 1995                                       
c***********************************************************************
      dimension a(7,7),b(7,7),c(6,6)

c-----------------------------------------------------------------------
c* To compute the determinant of 7-D matrix.

      ier=0
      call det7(a,det)
      if( abs(det).lt.0.00001 )then
          ier=1
          print*,'Cannot compute the inverse of 7-D matrix.'
          return
      endif

c-----------------------------------------------------------------------
c* To compute b(1,j),j=1,7

      do j=1,6
      do i=1,6
         i1=i+1
         j1=j+1
         c(i,j)=a(i1,j1)
      enddo
      enddo
      call det6(c,b(1,1))

      c(1,1)=a(1,2)
      c(1,2)=a(1,3)
      c(1,3)=a(1,4)
      c(1,4)=a(1,5)
      c(1,5)=a(1,6)
      c(1,6)=a(1,7)
      call det6(c,b(1,2))

      c(2,1)=a(2,2)
      c(2,2)=a(2,3)
      c(2,3)=a(2,4)
      c(2,4)=a(2,5)
      c(2,5)=a(2,6)
      c(2,6)=a(2,7)
      call det6(c,b(1,3))

      c(3,1)=a(3,2)
      c(3,2)=a(3,3)
      c(3,3)=a(3,4)
      c(3,4)=a(3,5)
      c(3,5)=a(3,6)
      c(3,6)=a(3,7)
      call det6(c,b(1,4))

      c(4,1)=a(4,2)
      c(4,2)=a(4,3)
      c(4,3)=a(4,4)
      c(4,4)=a(4,5)
      c(4,5)=a(4,6)
      c(4,6)=a(4,7)
      call det6(c,b(1,5))

      c(5,1)=a(5,2)
      c(5,2)=a(5,3)
      c(5,3)=a(5,4)
      c(5,4)=a(5,5)
      c(5,5)=a(5,6)
      c(5,6)=a(5,7)
      call det6(c,b(1,6))

      c(6,1)=a(6,2)
      c(6,2)=a(6,3)
      c(6,3)=a(6,4)
      c(6,4)=a(6,5)
      c(6,5)=a(6,6)
      c(6,6)=a(6,7)
      call det6(c,b(1,7))

c-----------------------------------------------------------------------
c* To compute b(2,j),j=1,7

      do i=1,6
         i1=i+1
         c(i,1)=a(i1,1)
      enddo
      do j=2,6
      do i=1,6
         i1=i+1
         j1=j+1
         c(i,j)=a(i1,j1)
      enddo
      enddo
      call det6(c,b(2,1))

      c(1,1)=a(1,1)
      c(1,2)=a(1,3)
      c(1,3)=a(1,4)
      c(1,4)=a(1,5)
      c(1,5)=a(1,6)
      c(1,6)=a(1,7)
      call det6(c,b(2,2))

      c(2,1)=a(2,1)
      c(2,2)=a(2,3)
      c(2,3)=a(2,4)
      c(2,4)=a(2,5)
      c(2,5)=a(2,6)
      c(2,6)=a(2,7)
      call det6(c,b(2,3))

      c(3,1)=a(3,1)
      c(3,2)=a(3,3)
      c(3,3)=a(3,4)
      c(3,4)=a(3,5)
      c(3,5)=a(3,6)
      c(3,6)=a(3,7)
      call det6(c,b(2,4))

      c(4,1)=a(4,1)
      c(4,2)=a(4,3)
      c(4,3)=a(4,4)
      c(4,4)=a(4,5)
      c(4,5)=a(4,6)
      c(4,6)=a(4,7)
      call det6(c,b(2,5))

      c(5,1)=a(5,1)
      c(5,2)=a(5,3)
      c(5,3)=a(5,4)
      c(5,4)=a(5,5)
      c(5,5)=a(5,6)
      c(5,6)=a(5,7)
      call det6(c,b(2,6))

      c(6,1)=a(6,1)
      c(6,2)=a(6,3)
      c(6,3)=a(6,4)
      c(6,4)=a(6,5)
      c(6,5)=a(6,6)
      c(6,6)=a(6,7)
      call det6(c,b(2,7))

c-----------------------------------------------------------------------
c* To compute b(3,j),j=1,7

      do i=1,6
         i1=i+1
         c(i,1)=a(i1,1)
         c(i,2)=a(i1,2)
      enddo
      do j=3,6
      do i=1,6
         i1=i+1
         j1=j+1
         c(i,j)=a(i1,j1)
      enddo
      enddo
      call det6(c,b(3,1))

      c(1,1)=a(1,1)
      c(1,2)=a(1,2)
      c(1,3)=a(1,4)
      c(1,4)=a(1,5)
      c(1,5)=a(1,6)
      c(1,6)=a(1,7)
      call det6(c,b(3,2))

      c(2,1)=a(2,1)
      c(2,2)=a(2,2)
      c(2,3)=a(2,4)
      c(2,4)=a(2,5)
      c(2,5)=a(2,6)
      c(2,6)=a(2,7)
      call det6(c,b(3,3))

      c(3,1)=a(3,1)
      c(3,2)=a(3,2)
      c(3,3)=a(3,4)
      c(3,4)=a(3,5)
      c(3,5)=a(3,6)
      c(3,6)=a(3,7)
      call det6(c,b(3,4))

      c(4,1)=a(4,1)
      c(4,2)=a(4,2)
      c(4,3)=a(4,4)
      c(4,4)=a(4,5)
      c(4,5)=a(4,6)
      c(4,6)=a(4,7)
      call det6(c,b(3,5))

      c(5,1)=a(5,1)
      c(5,2)=a(5,2)
      c(5,3)=a(5,4)
      c(5,4)=a(5,5)
      c(5,5)=a(5,6)
      c(5,6)=a(5,7)
      call det6(c,b(3,6))

      c(6,1)=a(6,1)
      c(6,2)=a(6,2)
      c(6,3)=a(6,4)
      c(6,4)=a(6,5)
      c(6,5)=a(6,6)
      c(6,6)=a(6,7)
      call det6(c,b(3,7))

c-----------------------------------------------------------------------
c* To compute b(4,j),j=1,7

      do i=1,6
         i1=i+1
         c(i,1)=a(i1,1)
         c(i,2)=a(i1,2)
         c(i,3)=a(i1,3)
      enddo
      do j=4,6
      do i=1,6
         i1=i+1
         j1=j+1
         c(i,j)=a(i1,j1)
      enddo
      enddo
      call det6(c,b(4,1))

      c(1,1)=a(1,1)
      c(1,2)=a(1,2)
      c(1,3)=a(1,3)
      c(1,4)=a(1,5)
      c(1,5)=a(1,6)
      c(1,6)=a(1,7)
      call det6(c,b(4,2))

      c(2,1)=a(2,1)
      c(2,2)=a(2,2)
      c(2,3)=a(2,3)
      c(2,4)=a(2,5)
      c(2,5)=a(2,6)
      c(2,6)=a(2,7)
      call det6(c,b(4,3))

      c(3,1)=a(3,1)
      c(3,2)=a(3,2)
      c(3,3)=a(3,3)
      c(3,4)=a(3,5)
      c(3,5)=a(3,6)
      c(3,6)=a(3,7)
      call det6(c,b(4,4))

      c(4,1)=a(4,1)
      c(4,2)=a(4,2)
      c(4,3)=a(4,3)
      c(4,4)=a(4,5)
      c(4,5)=a(4,6)
      c(4,6)=a(4,7)
      call det6(c,b(4,5))

      c(5,1)=a(5,1)
      c(5,2)=a(5,2)
      c(5,3)=a(5,3)
      c(5,4)=a(5,5)
      c(5,5)=a(5,6)
      c(5,6)=a(5,7)
      call det6(c,b(4,6))

      c(6,1)=a(6,1)
      c(6,2)=a(6,2)
      c(6,3)=a(6,3)
      c(6,4)=a(6,5)
      c(6,5)=a(6,6)
      c(6,6)=a(6,7)
      call det6(c,b(4,7))

c-----------------------------------------------------------------------
c* To compute b(5,j),j=1,7

      do i=1,6
         i1=i+1
         c(i,1)=a(i1,1)
         c(i,2)=a(i1,2)
         c(i,3)=a(i1,3)
         c(i,4)=a(i1,4)
         c(i,5)=a(i1,6)
         c(i,6)=a(i1,7)
      enddo
      call det6(c,b(5,1))

      c(1,1)=a(1,1)
      c(1,2)=a(1,2)
      c(1,3)=a(1,3)
      c(1,4)=a(1,4)
      c(1,5)=a(1,6)
      c(1,6)=a(1,7)
      call det6(c,b(5,2))

      c(2,1)=a(2,1)
      c(2,2)=a(2,2)
      c(2,3)=a(2,3)
      c(2,4)=a(2,4)
      c(2,5)=a(2,6)
      c(2,6)=a(2,7)
      call det6(c,b(5,3))

      c(3,1)=a(3,1)
      c(3,2)=a(3,2)
      c(3,3)=a(3,3)
      c(3,4)=a(3,4)
      c(3,5)=a(3,6)
      c(3,6)=a(3,7)
      call det6(c,b(5,4))

      c(4,1)=a(4,1)
      c(4,2)=a(4,2)
      c(4,3)=a(4,3)
      c(4,4)=a(4,4)
      c(4,5)=a(4,6)
      c(4,6)=a(4,7)
      call det6(c,b(5,5))

      c(5,1)=a(5,1)
      c(5,2)=a(5,2)
      c(5,3)=a(5,3)
      c(5,4)=a(5,4)
      c(5,5)=a(5,6)
      c(5,6)=a(5,7)
      call det6(c,b(5,6))

      c(6,1)=a(6,1)
      c(6,2)=a(6,2)
      c(6,3)=a(6,3)
      c(6,4)=a(6,4)
      c(6,5)=a(6,6)
      c(6,6)=a(6,7)
      call det6(c,b(5,7))

c-----------------------------------------------------------------------
c* To compute b(6,j),j=1,7

      do i=1,6
         i1=i+1
         c(i,1)=a(i1,1)
         c(i,2)=a(i1,2)
         c(i,3)=a(i1,3)
         c(i,4)=a(i1,4)
         c(i,5)=a(i1,5)
         c(i,6)=a(i1,7)
      enddo
      call det6(c,b(6,1))

      c(1,1)=a(1,1)
      c(1,2)=a(1,2)
      c(1,3)=a(1,3)
      c(1,4)=a(1,4)
      c(1,5)=a(1,5)
      c(1,6)=a(1,7)
      call det6(c,b(6,2))

      c(2,1)=a(2,1)
      c(2,2)=a(2,2)
      c(2,3)=a(2,3)
      c(2,4)=a(2,4)
      c(2,5)=a(2,5)
      c(2,6)=a(2,7)
      call det6(c,b(6,3))

      c(3,1)=a(3,1)
      c(3,2)=a(3,2)
      c(3,3)=a(3,3)
      c(3,4)=a(3,4)
      c(3,5)=a(3,5)
      c(3,6)=a(3,7)
      call det6(c,b(6,4))

      c(4,1)=a(4,1)
      c(4,2)=a(4,2)
      c(4,3)=a(4,3)
      c(4,4)=a(4,4)
      c(4,5)=a(4,5)
      c(4,6)=a(4,7)
      call det6(c,b(6,5))

      c(5,1)=a(5,1)
      c(5,2)=a(5,2)
      c(5,3)=a(5,3)
      c(5,4)=a(5,4)
      c(5,5)=a(5,5)
      c(5,6)=a(5,7)
      call det6(c,b(6,6))

      c(6,1)=a(6,1)
      c(6,2)=a(6,2)
      c(6,3)=a(6,3)
      c(6,4)=a(6,4)
      c(6,5)=a(6,5)
      c(6,6)=a(6,7)
      call det6(c,b(6,7))

c-----------------------------------------------------------------------
c* To compute b(7,j),j=1,7

      do i=1,6
         i1=i+1
         c(i,1)=a(i1,1)
         c(i,2)=a(i1,2)
         c(i,3)=a(i1,3)
         c(i,4)=a(i1,4)
         c(i,5)=a(i1,5)
         c(i,6)=a(i1,6)
      enddo
      call det6(c,b(7,1))

      c(1,1)=a(1,1)
      c(1,2)=a(1,2)
      c(1,3)=a(1,3)
      c(1,4)=a(1,4)
      c(1,5)=a(1,5)
      c(1,6)=a(1,6)
      call det6(c,b(7,2))

      c(2,1)=a(2,1)
      c(2,2)=a(2,2)
      c(2,3)=a(2,3)
      c(2,4)=a(2,4)
      c(2,5)=a(2,5)
      c(2,6)=a(2,6)
      call det6(c,b(7,3))

      c(3,1)=a(3,1)
      c(3,2)=a(3,2)
      c(3,3)=a(3,3)
      c(3,4)=a(3,4)
      c(3,5)=a(3,5)
      c(3,6)=a(3,6)
      call det6(c,b(7,4))

      c(4,1)=a(4,1)
      c(4,2)=a(4,2)
      c(4,3)=a(4,3)
      c(4,4)=a(4,4)
      c(4,5)=a(4,5)
      c(4,6)=a(4,6)
      call det6(c,b(7,5))

      c(5,1)=a(5,1)
      c(5,2)=a(5,2)
      c(5,3)=a(5,3)
      c(5,4)=a(5,4)
      c(5,5)=a(5,5)
      c(5,6)=a(5,6)
      call det6(c,b(7,6))

      c(6,1)=a(6,1)
      c(6,2)=a(6,2)
      c(6,3)=a(6,3)
      c(6,4)=a(6,4)
      c(6,5)=a(6,5)
      c(6,6)=a(6,6)
      call det6(c,b(7,7))

c-----------------------------------------------------------------------
c* To compute inverse of 7-D matrix.

      do i=1,7
      do j=1,7
         ind=(i+j)/2
         isn=(i+j)-2*ind
         if( isn.eq.0 )then
             isign=1
         else
             isign=-1
         endif
         b(i,j)=isign*b(i,j)/det
      enddo
      enddo

      return
      end

      subroutine det7(a,det)
c***********************************************************************
c description      : To compute the determinant of 7-D matrix.          
c                                                                       
c I/O parameters   :                                                    
c  input:  name,      type,       description                           
c          a(7,7)     real array  7-D matrix.                           
c  output: name,      type,       description                           
c          det        real        determinant of 6-D matrix.            
c                                                                       
c called fun./sub. : det6                                               
c                                                                       
c auther           : Deng Shiung Ming                                   
c                    - Institude for Information Industry               
c                                                                       
c create date      : Aug 24, 1995                                       
c***********************************************************************
      dimension a(7,7),b(6,6)

c-----------------------------------------------------------------------
c* To assign 6-D matrix.

      do j=1,6
      do i=1,6
         i1=i+1
         j1=j+1
         b(i,j)=a(i1,j1)
      enddo
      enddo

c-----------------------------------------------------------------------
c* To compute 6-D determinant.

      call det6(b,d11)

      b(1,1)=a(1,2)
      b(1,2)=a(1,3)
      b(1,3)=a(1,4)
      b(1,4)=a(1,5)
      b(1,5)=a(1,6)
      b(1,6)=a(1,7)
      call det6(b,d21)

      b(2,1)=a(2,2)
      b(2,2)=a(2,3)
      b(2,3)=a(2,4)
      b(2,4)=a(2,5)
      b(2,5)=a(2,6)
      b(2,6)=a(2,7)
      call det6(b,d31)

      b(3,1)=a(3,2)
      b(3,2)=a(3,3)
      b(3,3)=a(3,4)
      b(3,4)=a(3,5)
      b(3,5)=a(3,6)
      b(3,6)=a(3,7)
      call det6(b,d41)

      b(4,1)=a(4,2)
      b(4,2)=a(4,3)
      b(4,3)=a(4,4)
      b(4,4)=a(4,5)
      b(4,5)=a(4,6)
      b(4,6)=a(4,7)
      call det6(b,d51)

      b(5,1)=a(5,2)
      b(5,2)=a(5,3)
      b(5,3)=a(5,4)
      b(5,4)=a(5,5)
      b(5,5)=a(5,6)
      b(5,6)=a(5,7)
      call det6(b,d61)

      b(6,1)=a(6,2)
      b(6,2)=a(6,3)
      b(6,3)=a(6,4)
      b(6,4)=a(6,5)
      b(6,5)=a(6,6)
      b(6,6)=a(6,7)
      call det6(b,d71)

c-----------------------------------------------------------------------
c* To compute 7-D determinant.

      det=a(1,1)*d11-a(2,1)*d21+a(3,1)*d31-a(4,1)*d41+a(5,1)*d51
     *   -a(6,1)*d61+a(7,1)*d71

      return
      end

      subroutine det6(a,det)
c***********************************************************************
c description      : To compute the determinant of 6-D matrix.          
c                                                                       
c I/O parameters   :                                                    
c  input:  name,      type,       description                           
c          a(6,6)     real array  6-D matrix.                           
c  output: name,      type,       description                           
c          det        real        determinant of 6-D matrix.            
c                                                                       
c called fun./sub. : det5                                               
c                                                                       
c auther           : Deng Shiung Ming                                   
c                    - Institude for Information Industry               
c                                                                       
c create date      : Aug 24, 1995                                       
c***********************************************************************
      dimension a(6,6),b(5,5)

c-----------------------------------------------------------------------
c* To assign 5-D matrix.

      do j=1,5
      do i=1,5
         i1=i+1
         j1=j+1
         b(i,j)=a(i1,j1)
      enddo
      enddo

c-----------------------------------------------------------------------
c* To compute 5-D determinant.

      call det5(b,d11)

      b(1,1)=a(1,2)
      b(1,2)=a(1,3)
      b(1,3)=a(1,4)
      b(1,4)=a(1,5)
      b(1,5)=a(1,6)
      call det5(b,d21)

      b(2,1)=a(2,2)
      b(2,2)=a(2,3)
      b(2,3)=a(2,4)
      b(2,4)=a(2,5)
      b(2,5)=a(2,6)
      call det5(b,d31)

      b(3,1)=a(3,2)
      b(3,2)=a(3,3)
      b(3,3)=a(3,4)
      b(3,4)=a(3,5)
      b(3,5)=a(3,6)
      call det5(b,d41)

      b(4,1)=a(4,2)
      b(4,2)=a(4,3)
      b(4,3)=a(4,4)
      b(4,4)=a(4,5)
      b(4,5)=a(4,6)
      call det5(b,d51)

      b(5,1)=a(5,2)
      b(5,2)=a(5,3)
      b(5,3)=a(5,4)
      b(5,4)=a(5,5)
      b(5,5)=a(5,6)
      call det5(b,d61)

c-----------------------------------------------------------------------
c* To compute 6-D determinant.

      det=a(1,1)*d11-a(2,1)*d21+a(3,1)*d31-a(4,1)*d41+a(5,1)*d51
     *   -a(6,1)*d61

      return
      end

      subroutine det5(a,det)
c***********************************************************************
c description      : To compute the determinant of 5-D matrix.          
c                                                                       
c I/O parameters   :                                                    
c  input:  name,      type,       description                           
c          a(5,5)     real array  5-D matrix.                           
c  output: name,      type,       description                           
c          det        real        determinant of 5-D matrix.            
c                                                                       
c called fun./sub. : det4                                               
c                                                                       
c auther           : Deng Shiung Ming                                   
c                    - Institude for Information Industry               
c                                                                       
c create date      : Aug 24, 1995                                       
c***********************************************************************
      dimension a(5,5),b(4,4)

c-----------------------------------------------------------------------
c* To assign 4-D matrix.

      do j=1,4
      do i=1,4
         i1=i+1
         j1=j+1
         b(i,j)=a(i1,j1)
      enddo
      enddo

c-----------------------------------------------------------------------
c* To compute 4-D determinant.

      call det4(b,d11)

      b(1,1)=a(1,2)
      b(1,2)=a(1,3)
      b(1,3)=a(1,4)
      b(1,4)=a(1,5)
      call det4(b,d21)

      b(2,1)=a(2,2)
      b(2,2)=a(2,3)
      b(2,3)=a(2,4)
      b(2,4)=a(2,5)
      call det4(b,d31)

      b(3,1)=a(3,2)
      b(3,2)=a(3,3)
      b(3,3)=a(3,4)
      b(3,4)=a(3,5)
      call det4(b,d41)

      b(4,1)=a(4,2)
      b(4,2)=a(4,3)
      b(4,3)=a(4,4)
      b(4,4)=a(4,5)
      call det4(b,d51)

c-----------------------------------------------------------------------
c* To compute 5-D determinant.

      det=a(1,1)*d11-a(2,1)*d21+a(3,1)*d31-a(4,1)*d41+a(5,1)*d51

      return
      end

      subroutine det4(a,det)
c***********************************************************************
c description      : To compute the determinant of 4-D matrix.          
c                                                                       
c I/O parameters   :                                                    
c  input:  name,      type,       description                           
c          a(4,4)     real array  4-D matrix.                           
c  output: name,      type,       description                           
c          det        real        determinant of 4-D matrix.            
c                                                                       
c note:                                                                 
c          det=  a11*(a22*a33*a44+a23*a34*a42+a24*a32*a43               
c                    -a24*a33*a42-a22*a34*a43-a23*a32*a44)              
c               -a21*(a12*a33*a44+a13*a34*a42+a14*a32*a43               
c                    -a14*a33*a42-a12*a34*a43-a13*a32*a44)              
c               +a31*(a12*a23*a44+a13*a24*a42+a14*a22*a43               
c                    -a14*a23*a42-a12*a24*a43-a13*a22*a44)              
c               -a41*(a12*a23*a34+a13*a24*a32+a14*a22*a33               
c                    -a14*a23*a32-a12*a24*a33-a13*a22*a34)              
c                                                                       
c auther           : Deng Shiung Ming                                   
c                    - Institude for Information Industry               
c                                                                       
c create date      : Aug 24, 1995                                       
c***********************************************************************
      dimension a(4,4)

c-----------------------------------------------------------------------
c* To compute determine of 4-D matrix
 
      a11=a(2,2)*a(3,3)*a(4,4)+a(2,3)*a(3,4)*a(4,2)+a(2,4)*a(3,2)*a(4,3)
     *   -a(2,4)*a(3,3)*a(4,2)-a(2,2)*a(3,4)*a(4,3)-a(2,3)*a(3,2)*a(4,4)

      a21=a(1,2)*a(3,3)*a(4,4)+a(1,3)*a(3,4)*a(4,2)+a(1,4)*a(3,2)*a(4,3)
     *   -a(1,4)*a(3,3)*a(4,2)-a(1,2)*a(3,4)*a(4,3)-a(1,3)*a(3,2)*a(4,4)

      a31=a(1,2)*a(2,3)*a(4,4)+a(1,3)*a(2,4)*a(4,2)+a(1,4)*a(2,2)*a(4,3)
     *   -a(1,4)*a(2,3)*a(4,2)-a(1,2)*a(2,4)*a(4,3)-a(1,3)*a(2,2)*a(4,4)

      a41=a(1,2)*a(2,3)*a(3,4)+a(1,3)*a(2,4)*a(3,2)+a(1,4)*a(2,2)*a(3,3)
     *   -a(1,4)*a(2,3)*a(3,2)-a(1,2)*a(2,4)*a(3,3)-a(1,3)*a(2,2)*a(3,4)

      det=a(1,1)*a11-a(2,1)*a21+a(3,1)*a31-a(4,1)*a41

      return
      end
