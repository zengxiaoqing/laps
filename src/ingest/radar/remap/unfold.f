       subroutine unfold(nogate,nobeam,dat,beam,velny,bmiss)
c***********************************************************************
c description      : To unfold Vr PPI data for LAPS system.
c                                                                       
c I/O parameters   :                                                    
c  I/O/W   name,            type,      description
c    I     nogate           integer    gate number in one ray.
c    I     nobeam           integer    total ray number of PPI sweep scan.
c   I/O  dat(nogate,nobeam) real array radial velocity in unit of m/s.
c    I     beam(nobeam)     real array azimthual angle. (degrees)
c    I     velny            real       Nyquist velocity in unit of m/s.
c    I     bmiss            real       missing value for radial velocity.
c                                      if( dat(i,j).eq.bmiss )then
c                                          dat(i,j) is missing data. 
c Called Function:
c   get_azimuindex,  del_isolated,  simple_unf,  get_vrgroup,  gpunfold
c
c Date :
c   Jul. 10, 2007 (S.-M. Deng)
c***********************************************************************

       implicit none
       integer nogate,nobeam
       real dat(nogate,nobeam),beam(nobeam)
       real velny,bmiss

       integer nog,nob
       parameter( nog=940,nob=430 )
       integer idat(nog*nob)
       integer ind(nog*nob)    !=0, no group data; =-1, cannot group data;
                               !=-2, missing data.

       integer jend,ivad,ier,jend4
       integer i,j,j1 

c -----to check Nyquist velocity. If velny < 10 m/s, donot do unfold.

       if( velny.lt.12. )return

c -----to get the last index along azimthual direction.

       call get_azimuindex(nobeam,beam,jend,ivad)
       if( ivad.ne.0 )then
           print*,'Input Vr PPI data is not periodic azimuthal angle.'
c          print*,'We unfold periodic azimuthal angle data, '
c          print*,' so this input Vr data isnot be unfolded.'
c          return
       endif

c -----to delete isolated data.

       call del_isolated(nogate,nobeam,jend,dat,velny,bmiss,ivad)

c ---- to simple unfold by 3 x 3 cell.

       call simple_unf(nogate,nobeam,jend,dat,velny,bmiss,ivad)

c ---- to unfold radial velocity data by group method.

       if( ivad.eq.0 )then
           jend4=jend+4
       else
           jend4=nobeam
       endif
       call gpunfold(nogate,nobeam,jend,jend4,dat,velny,bmiss,ivad
     1              ,idat,ind)
       if( nobeam.gt.jend )then
           do j=jend+1,nobeam
              j1=j-jend
              beam(j)=beam(j1)
              do i=1,nogate
                 dat(i,j)=dat(i,j1)
              enddo
           enddo
       endif

       return
       end

       subroutine gpunfold(nogate,nobeam,jend,jend4,dat,velny
     1                    ,bmiss,ivad,idat,ind)
c***********************************************************************
c description      : To unfold radial velocity data group method.
c                                                                       
c I/O parameters   :                                                    
c  I/O/W   name,            type,      description
c    I     nogate           integer    gate number in one ray.
c    I     nobeam           integer    total ray number of PPI sweep scan.
c    I     jend             integer    the last index along azimthual direction.
c    I     jend4            integer    =jend+4
c   I/O  dat(nogate,nobeam) real array radial velocity in unit of m/s.
c    I     velny            real       Nquist velocity in unit of m/s.
c    I     bmiss            real       missing value for radial velocity.
c    I     ivad             integer    ppi scanning message.
c                                      =0, is 360 degrees scan;  =1, isnot
c    W   idat(nogate,jend4) int array  radial velocity in unit of 100 * m/s.
c    O    ind(nogate,jend4) int array  group index.
c
c Called Function: 
c    get_ishr,     subgroup,     merge_bound,  get_boundary, sortgp,
c    mergegp,      get_zeroline, near_unf,     azim_unf,     ray_unf,       
c    ray_unf1
c***********************************************************************
       implicit none
       integer nogate,nobeam,jend,jend4,ivad
       real dat(nogate,nobeam),velny,bmiss
       integer idat(nogate,jend4)
       integer ind(nogate,jend4)   !=0, no group data; =-1, cannot group data;
                                   !=-2, missing data.

       integer shear
       parameter( shear=250 ) !=2.5 m/s
       integer miss,shearp,shearn,nquist
       integer i,j,j1,k,jend1,jend2,jend3,ii,jj,ij
       integer i1,i2,j2,kk

       integer gst,gst1,ist,jst,ier,n,ngst,nmax
       integer numzero,zmax,kst,idst,istatus

       integer nog,nob
       parameter( nog=940,nob=430 )
       integer ishr(8*nog*nob)  !=1, small shear index.
       integer ibd(nog*nob)     !=0, isnot boundary; =1, boundary data.
       integer idgp(nog*nob),numgp(nog*nob)
       integer iwkd(nog*nob)
       integer iwzi(nog),jwz1(nog),jwz2(nog)
       integer numz,num,ivalue
       real value

       nquist=100*velny
       shearp=1.5*nquist
       shearn=-shearp 
       miss=-32767
       if( ivad.eq.0 )then
           jend1=jend+1
           jend2=jend+2
           jend3=jend+3
       else
           jend3=jend-1
           jend2=jend-2
       endif
      
       if( ivad.eq.0 )then
           do j=1,jend
           do i=1,nogate
              j1=j+2
              if( abs(dat(i,j)).lt.99. )then
                  idat(i,j1)=100*dat(i,j)
              else
                  idat(i,j1)=miss
              endif
           enddo
           enddo
           do i=1,nogate
              idat(i,1)=idat(i,jend1)
              idat(i,2)=idat(i,jend2)
              idat(i,jend3)=idat(i,3)
              idat(i,jend4)=idat(i,4)
           enddo
       else
           do j=1,nobeam
           do i=1,nogate
              if( abs(dat(i,j)).lt.99. )then
                  idat(i,j)=100*dat(i,j)
              else
                  idat(i,j)=miss
              endif
           enddo
           enddo
       endif

c -----to fill blank data.

       do j=1,jend4
       do i=2,nogate-1,1
          if( (idat(i,j).eq.miss).and.(idat(i-1,j).ne.miss).and.
     1        (idat(i+1,j).ne.miss) )then
              idat(i,j)=idat(i-1,j)
          endif 
       enddo
       enddo

       do j=2,jend3
       do i=1,nogate
          if( (idat(i,j).eq.miss).and.(idat(i,j-1).ne.miss).and.
     1        (idat(i,j+1).ne.miss) )then
              idat(i,j)=idat(i,j-1)
          endif
       enddo
       enddo

c -----to get shear index.

       call get_ishr(nogate,jend4,miss,shear,idat,ishr)

       print*,' '
       print*,'starting group data.'
c -----to group data.

       do j=1,jend4
       do i=1,nogate
          if( idat(i,j).ne.miss )then
              ind(i,j)=0
          else
              ind(i,j)=-2
          endif
       enddo
       enddo

       gst1=1
       do j=3,jend2,1
       do i=3,nogate-2,1
          if( ind(i,j).eq.0 )then
              ist=i
              jst=j
              call subgroup(nogate,jend4,ishr,ind,ist,jst,gst1,ier)
              if( ier.ne.0 )then
                  ind(ist,jst)=-1
              else
                  gst=gst1
                  gst1=gst1+1
              endif
          endif
       enddo
       enddo
       do j=2,jend3,1
       do i=2,nogate-1,1
          if( (ind(i,j).eq.-1).or.(ind(i,j).eq.0) )then
              do jj=j-1,j+1,1
              do ii=i-1,i+1,1
                 if( ind(ii,jj).gt.0 )then
                     if( abs(idat(ii,jj)-idat(i,j)).lt.shear )then
                         ind(i,j)=ind(ii,jj)
                         go to 5
                     endif
                     if( (idat(i,j)-idat(ii,jj)).gt.shearp )then
                         ind(i,j)=ind(ii,jj)
                         idat(i,j)=idat(i,j)-2*nquist
                         go to 5
                     endif
                     if( (idat(i,j)-idat(ii,jj)).lt.shearn )then
                         ind(i,j)=ind(ii,jj)
                         idat(i,j)=idat(i,j)+2*nquist
                         go to 5
                     endif
                 endif
              enddo
              enddo
          endif
 5        continue
       enddo
       enddo
       print*,'ending group data, gst: ',gst
       print*,' '

c -----to merge group index for boundary area.

       if( gst.le.1 )go to 59
       if( ivad.eq.0 )then 
           do k=1,10 
              call merge_bound(nogate,jend4,gst,ind,istatus)
              if( istatus.eq.0 )go to 6
           enddo
 6         continue
       endif
 
c -----to merge sub-group data.

       print*,'merge sub-group data.'
       do k=1,gst
          n=0
          do j=1,jend4
          do i=1,nogate
             if( ind(i,j).eq.k )n=n+1
          enddo
          enddo
          idgp(k)=k
          numgp(k)=n
       enddo
       call sortgp(gst,idgp,numgp)
       call get_boundary(nogate,jend4,ind,ibd)
       call mergegp(nogate,jend4,idat,ibd,ind,shear,gst,idgp,numgp)
       call sortgp(gst,idgp,numgp)
       ngst=gst
       do i=2,gst
          if( numgp(i).eq.0 )then
              ngst=i-1
              go to 10
          endif
       enddo
 10    continue
       if( ngst.le.1 )go to 59

       print*,'gst: ',gst,' ngst: ',ngst 
       print*,' '
      
c -----to get starting group.

       if( ivad.eq.0 )then
           ij=0
           do j=1,jend
           do i=1,nogate
              j1=j+2
              ij=ij+1
              iwkd(ij)=idat(i,j1)
           enddo
           enddo 
           call get_zeroline(nogate,jend,iwkd,miss,numz,iwzi
     1                      ,jwz1,jwz2,ier)
       endif
       nmax=10
       if( nmax.gt.ngst )nmax=ngst
       zmax=0
       idst=idgp(1)
       if( (ier.eq.0).and.(ivad.eq.0) )then
           do k=1,nmax
              numzero=0
              do kk=1,numz
                 i1=iwzi(kk)-5
                 if( i1.lt.1 )i1=1
                 i2=iwzi(kk)+5
                 if( i2.gt.nogate )i2=nogate
                 j1=jwz1(kk)+2-5
                 if( j1.lt.1 )j1=1
                 j2=jwz1(kk)+2+5
                 if( j2.gt.jend4 )j2=jend4
                 do j=j1,j2,1
                 do i=i1,i2,1
                    if( ind(j,j).eq.idgp(k) )then
                        if( abs(idat(i,j)).lt.200 )then
                            numzero=numzero+1
                        endif
                    endif
                 enddo
                 enddo 
                 if( jwz2(kk).gt.0 )then
                     j1=jwz2(kk)+2-5
                     if( j1.lt.1 )j1=1
                     j2=jwz2(kk)+2+5
                     if( j2.gt.jend4 )j2=jend4
                     do j=j1,j2,1
                     do i=i1,i2,1
                        if( ind(j,j).eq.idgp(k) )then
                            if( abs(idat(i,j)).lt.200 )then
                                numzero=numzero+1
                            endif
                        endif
                     enddo
                     enddo 
                 endif
              enddo
              if( numzero.gt.zmax )then
                  zmax=numzero
                  idst=idgp(k)
              endif
           enddo
       else
           do k=1,nmax
              numzero=0
              do j=3,jend3,1
              do i=1,nogate
                 if( ind(i,j).eq.idgp(k) )then
                     if( abs(idat(i,j)).lt.100 )then
                         numzero=numzero+1
                     endif
                 endif
              enddo
              enddo
              if( numzero.gt.zmax )then
                  zmax=numzero
                  idst=idgp(k)
              endif
           enddo
       endif

c -----to unfold data by near method.
       
       print*,'unfold data by near method'
       do 30 k=1,100
          call get_boundary(nogate,jend4,ind,ibd)
          call near_unf(nogate,jend4,ivad,idat,ibd,ind,nquist,shear
     1                 ,ngst,idgp,numgp,idst,istatus)
          if( istatus.eq.0 )go to 40
          call sortgp(ngst,idgp,numgp)
          gst=ngst
          do i=2,gst
             if( numgp(i).eq.0 )then
                 ngst=i-1
                 go to 20
             endif
          enddo
 20       continue
          print*,k,'gst: ',gst,' ngst: ',ngst
 30    continue
 40    continue
       if( ngst.eq.1 )go to 59

       do j=2,jend3,1
       do i=2,nogate-1,1
          if( (ind(i,j).eq.-1).or.(ind(i,j).eq.0) )then
              do jj=j-1,j+1,1
              do ii=i-1,i+1,1
                 if( ind(ii,jj).gt.0 )then
                     if( abs(idat(ii,jj)-idat(i,j)).lt.shear )then
                         ind(i,j)=ind(ii,jj)
                         go to 45
                     endif
                     if( (idat(i,j)-idat(ii,jj)).gt.shearp )then
                         ind(i,j)=ind(ii,jj)
                         idat(i,j)=idat(i,j)-2*nquist
                         go to 45
                     endif
                     if( (idat(i,j)-idat(ii,jj)).lt.shearn )then
                         ind(i,j)=ind(ii,jj)
                         idat(i,j)=idat(i,j)+2*nquist
                         go to 45
                     endif
                 endif
              enddo
              enddo
          endif
 45       continue
       enddo
       enddo

c -----to unfold data by ray direction mothod.

       print*,' '
       print*,'unfold along ray direction'
       call ray_unf(nogate,jend4,idat,ind,nquist,shear
     1             ,ngst,idgp,numgp,idst,istatus)
       if( istatus.eq.1 )then 
           call sortgp(ngst,idgp,numgp)
           gst=ngst
           do i=2,gst
              if( numgp(i).eq.0 )then
                  ngst=i-1
                  go to 50
              endif
           enddo
 50        continue
       endif
       if( ngst.le.1 )go to 59
       print*,'gst: ',gst,' ngst: ',ngst

       if( (ivad.eq.0).and.(numgp(2).gt.20) )then

c -----to unfold data by azimithal direction mothod.

           print*,' '
           print*,'unfold along azimuthal direction'
           call azim_unf(nogate,jend4,idat,ind,nquist,shear
     1                  ,ngst,idgp,numgp,idst,istatus)
           if( istatus.eq.1 )then 
               call sortgp(ngst,idgp,numgp)
               gst=ngst
               do i=2,gst
                  if( numgp(i).eq.0 )then
                      ngst=i-1
                      go to 51
                  endif
               enddo
 51            continue
           endif
           if( ngst.le.1 )go to 59
           print*,'gst: ',gst,' ngst: ',ngst

c -----to unfold data by ray direction mothod again.

           if( numgp(2).gt.20 )then
               print*,' '
               print*,'unfold along ray direction again'
               call ray_unf(nogate,jend4,idat,ind,nquist,shear
     1                     ,ngst,idgp,numgp,idst,istatus)
               if( istatus.eq.1 )then 
                   call sortgp(ngst,idgp,numgp)
                   gst=ngst
                   do i=2,gst
                      if( numgp(i).eq.0 )then
                          ngst=i-1
                          go to 52
                      endif
                   enddo
 52                continue
               else
                   go to 59
               endif
               if( ngst.le.1 )go to 59
               print*,'gst: ',gst,' ngst: ',ngst
           endif

       endif

c -----to unfold data by ray direction mothod re-again.

       if( numgp(2).gt.20 )then
           print*,' '
           print*,'unfold along ray direction re-again'
           call ray_unf1(nogate,jend4,idat,ind,nquist,shear
     1                  ,ngst,idgp,numgp,idst,istatus)
           if( istatus.eq.1 )then 
               call sortgp(ngst,idgp,numgp)
               gst=ngst
               do i=2,gst
                  if( numgp(i).eq.0 )then
                      ngst=i-1
                      go to 53
                  endif
               enddo
 53            continue
           else
               go to 59
           endif
           print*,'gst: ',gst,' ngst: ',ngst
       endif

       if( (ivad.eq.0).and.(numgp(2).gt.20) )then

c -----to unfold data by azimithal direction mothod again.

           print*,' '
           print*,'unfold along azimuthal direction again'
           call azim_unf(nogate,jend4,idat,ind,nquist,shear
     1                  ,ngst,idgp,numgp,idst,istatus)
           if( istatus.eq.1 )then 
               call sortgp(ngst,idgp,numgp)
               gst=ngst
               do i=2,gst
                  if( numgp(i).eq.0 )then
                      ngst=i-1
                      go to 54
                  endif
               enddo
 54            continue
           endif
           if( ngst.le.1 )go to 59
           print*,'gst: ',gst,' ngst: ',ngst

c -----to unfold data by ray direction mothod again.

           if( numgp(2).gt.20 )then
               print*,' '
               print*,'unfold along ray direction again'
               call ray_unf1(nogate,jend4,idat,ind,nquist,shear
     1                      ,ngst,idgp,numgp,idst,istatus)
               if( istatus.eq.1 )then 
                   call sortgp(ngst,idgp,numgp)
                   gst=ngst
                   do i=2,gst
                      if( numgp(i).eq.0 )then
                          ngst=i-1
                          go to 55
                      endif
                   enddo
 55                continue
               else
                   go to 59
               endif
               if( ngst.le.1 )go to 59
               print*,'gst: ',gst,' ngst: ',ngst
           endif

       endif
 59    continue

c -----to check unfold data again.

       if( ivad.eq.0 )then
           ij=0
           do j=1,jend
           do i=1,nogate
              j1=j+2
              ij=ij+1
              iwkd(ij)=idat(i,j1)
           enddo
           enddo 
           call get_zeroline(nogate,jend,iwkd,miss,numz,iwzi
     1                      ,jwz1,jwz2,ier)
           if( ier.ne.0 )go to 60
           num=0
           value=0.
           do j=1,numz
              j1=jwz1(j)+jend/2
              if( j1.gt.jend )j1=j1-jend
              j1=j1+2
              ii=iwzi(j)
              i1=iwzi(j)+1
              if( i1.lt.1 )i1=1
              i2=iwzi(j)-1 
              if( i2.lt.nogate )i2=nogate
              if( idat(ii,j1).ne.miss )then
                  num=num+1
                  value=value+idat(ii,j1)
              endif
              if( idat(i1,j1).ne.miss )then
                  num=num+1
                  value=value+idat(i1,j1)
              endif
              if( idat(i2,j1).ne.miss )then
                  num=num+1
                  value=value+idat(i2,j1)
              endif
           enddo
           if( num.lt.20 )go to 60
           ivalue=value/float(num)
           if( ivalue.lt.(-1.2*nquist) )then
               do j=1,jend4
               do i=1,nogate
                  if( ind(i,j).eq.idst )then
                      idat(i,j)=idat(i,j)+2*nquist
                  endif
               enddo
               enddo
           endif
           if( ivalue.gt.(1.2*nquist) )then
               do j=1,jend4
               do i=1,nogate
                  if( ind(i,j).eq.idst )then
                      idat(i,j)=idat(i,j)-2*nquist
                  endif
               enddo
               enddo
           endif
       endif

 60    continue

       do j=1,jend
       do i=1,nogate
          j1=j+2
          if( abs(dat(i,j)).lt.99.0 )then
              dat(i,j)=0.01*idat(i,j1)
          endif
       enddo
       enddo

       return
       end

       subroutine get_zeroline(ng,nb,idat,miss,num,indi,jwz1,jwz2,ier)
c***********************************************************************
c Subroutine/Function : get_zeroline
c
c Usage :
c    call get_zeroline(ng,nb,idat,miss,num,indi,jwz1,jwz2,ier)
c
c Description      : To get zero line j-index from zero value data.
c
c Arguments :
c  I/O/W   name,      type,       description
c    I     ng         integer     first dimension of input data. 
c    I     nb         integer     second dimension of input data.
c    I    idat(ng,nb) int array   radial velocity.
c    I     miss       integer     missing or bad value.
c    O     num        integer     number of output data.
c    O     indi(ng)   integer     i-index.
c    O     jwz1(ng)   integer     first zero j-index.
c    O     jwz2(ng)   integer     second zero j-index.
c    O     ier        integer     error message. (=0, success; =1, failure)
c
c Modules Called : find_zero,  find_zero1
c***********************************************************************
       implicit none
       integer ng,nb,miss,num,ier
       integer idat(ng,nb),indi(ng),jwz1(ng),jwz2(ng)
       integer nog,nob
       parameter( nog=940,nob=430 )
       integer iwk(nob),iw1(nog),iw2(nog),iwx(nog),iwd(nog)
       integer jaa(4,nog),jbb(4,nog)
       integer idev(4),idv1,idv2

       integer nmax
       parameter( nmax=10 )

       integer i,j,jj,k,n,kk,nn
       integer jz1,jz2,nb2,nb23,nb6,ierr
       real sx,sx2,sy,sxy,det
       real a,b

       ier=1
       num=0
       nb2=nb/2
       nb23=2*nb/3
       nb6=nb/6
       n=0
       do i=1,ng
          jwz1(i)=0
          jwz2(i)=0
       enddo
       do i=1,ng
          do j=1,nb
             iwk(j)=idat(i,j)
          enddo
          call find_zero(nb,iwk,miss,jz1,jz2,ierr)
          if( ierr.eq.0 )then
              n=n+1
              indi(n)=i
              jwz1(n)=jz1
              jwz2(n)=jz2
          endif
       enddo
       if( n.lt.nmax )go to 100

       do k=1,4
          jj=nb*(k-1)/8+1
          idev(k)=0.
          do j=1,n
             idv1=abs(jwz1(j)-jj)
             if( idv1.gt.nb23 )idv1=nb-idv1
             if( jwz2(j).gt.0 )then
                 idv2=abs(jwz2(j)-jj)
                 if( idv2.gt.nb23 )idv2=nb-idv2
                 if( idv1.lt.idv2 )then
                     jaa(k,j)=jwz1(j)
                     jbb(k,j)=jwz2(j)
                     idev(k)=idev(k)+idv1
                 else
                     jbb(k,j)=jwz1(j)
                     jaa(k,j)=jwz2(j)
                     idev(k)=idev(k)+idv2
                 endif
             else
                 jaa(k,j)=jwz1(j)
                 jbb(k,j)=0
                 idev(k)=idev(k)+idv1
             endif
          enddo
       enddo
       kk=1
       idv1=idev(1)
       do k=2,4
          if( idv1.gt.idev(k) )then
              kk=k
              idv1=idev(k)
          endif
       enddo
       jj=nb*(kk-1)/8+1
       nn=0
       do j=1,n
          idv1=jaa(kk,j)-jj
          if( idv1.gt.nb23 )idv1=idv1-nb
          if( abs(idv1).lt.nb6 )then
              nn=nn+1
              iw1(nn)=jaa(kk,j)
              iw2(nn)=jbb(kk,j)
              iwx(nn)=indi(j)
              iwd(nn)=idv1
          endif
       enddo

       if( nn.lt.nmax )go to 100
       sx=0.
       sx2=0.
       sy=0.
       sxy=0.
       do k=1,nn
          sx=sx+iwx(k)
          sx2=sx2+iwx(k)*iwx(k)
          sy=sy+iwd(k)
          sxy=sxy+iwx(k)*iwd(k)
       enddo
       det=nn*sx2-sx*sx
       if( abs(det).lt.1. )return
       a=(sy*sx2-sxy*sx)/det
       b=(nn*sxy-sx*sy)/det
       num=0
       do k=1,nn
          idv1=abs(a+b*iwx(k)-iwd(k))
          if( idv1.lt.20 )then
              num=num+1
              indi(num)=iwx(k)
              jwz1(num)=iw1(k)
              jwz2(num)=iw2(k)
              if( iw2(k).gt.0 )then
                  idv2=abs(abs(iw1(k)-iw2(k))-nb2)
                  if( idv2.gt.30 )then
                      jwz2(num)=0
                  endif
              endif
          endif
       enddo
       if( num.ge.nmax )then
           ier=0
           return
       endif
 
 100   continue
       do i=1,ng
          jwz1(i)=0
          jwz2(i)=0
       enddo
       n=0
       do i=1,ng
          do j=1,nb
             iwk(j)=idat(i,j)
          enddo
          call find_zero1(nb,iwk,miss,jz1,ierr)
          if( ierr.eq.0 )then
              n=n+1
              indi(n)=i
              jwz1(n)=jz1
          endif
       enddo
       if( n.lt.nmax )return

       do k=1,4
          jj=nb*(k-1)/8+1
          idev(k)=0.
          do j=1,n
             idv1=abs(jwz1(j)-jj)
             if( idv1.gt.nb23 )idv1=nb-idv1
             jaa(k,j)=jwz1(j)
             idev(k)=idev(k)+idv1 
          enddo
       enddo
       kk=1
       idv1=idev(1)
       do k=2,4
          if( idv1.gt.idev(k) )then
              kk=k
              idv1=idev(k)
          endif
       enddo
       jj=nb*(kk-1)/8+1
       nn=0
       do j=1,n
          idv1=jaa(kk,j)-jj
          if( idv1.gt.nb23 )idv1=idv1-nb
          if( abs(idv1).lt.nb6 )then
              nn=nn+1
              iw1(nn)=jaa(kk,j)
              iwx(nn)=indi(j)
              iwd(nn)=idv1
          endif
       enddo

       if( nn.lt.nmax )return
       sx=0.
       sx2=0.
       sy=0.
       sxy=0.
       do k=1,nn
          sx=sx+iwx(k)
          sx2=sx2+iwx(k)*iwx(k)
          sy=sy+iwd(k)
          sxy=sxy+iwx(k)*iwd(k)
       enddo
       det=nn*sx2-sx*sx
       if( abs(det).lt.1. )return
       a=(sy*sx2-sxy*sx)/det
       b=(nn*sxy-sx*sy)/det
       num=0
       do k=1,nn
          idv1=abs(a+b*iwx(k)-iwd(k))
          if( idv1.lt.20 )then
              num=num+1
              indi(num)=iwx(k)
              jwz1(num)=iw1(k)
              jwz2(num)=0
          endif
       enddo
       if( num.ge.nmax )then
           ier=0
           return
       endif
 
       return
       end

       subroutine find_zero1(nb,idat,miss,jz1,ier)
c***********************************************************************
c Subroutine/Function : find_zero1
c
c Usage :
c    call find_zero1(nb,idat,miss,jz1,ier)
c
c Description      : To find the azimuthal position of zero value.
c
c Arguments :
c  I/O/W   name,      type,       description
c    I     nb         integer     total ray number of input sweep scan layer.
c    I     idat(nb)   int array   radial velocity.
c    I     miss       integer     missing or bad value.
c    O     jz1        integer     the first index of zero value.
c                                 =0, if no zero value
c    O     ier        integer     error message. (=0, success; =1, failure)
c
c Modules Called : none
c***********************************************************************
       implicit none
       integer nb,miss,jz1,ier
       integer idat(nb)

       integer loop
       parameter( loop=10 )
       integer ixx(21),iyy(21),jwork(loop),jz(loop)
       integer jconst,idev,i,ii,k,kk,nn,i1,i2,n,ix1,ix2,iy1,iy2,ixy,idet
       integer isum,ivmax
       real a,b,sum

       jz1=0
       ier=1
       jconst=-40000
       idev=200
       ivmax=600
       do i=1,loop
          jz(i)=0
          jwork(i)=0
       enddo

       kk=nb/10
       nn=0
       do 10 k=1,kk
          i1=10*(k-1)-10
          i2=10*(k-1)+10
          n=0
          do i=i1,i2,1
             ii=i
             if( ii.lt.1 )ii=ii+nb
             if( ii.gt.nb )ii=ii-nb
             if( (idat(ii).ne.miss).and.(abs(idat(ii)).lt.ivmax) )then
                 n=n+1
                 ixx(n)=i-i1+1
                 iyy(n)=idat(ii)
             endif
          enddo
          if( n.lt.8 )go to 10
          ix1=0
          ix2=0
          iy1=0
          ixy=0
          do i=1,n
             ix1=ix1+ixx(i)
             ix2=ix2+ixx(i)*ixx(i)
             iy1=iy1+iyy(i)
             ixy=ixy+ixx(i)*iyy(i)
          enddo
          idet=n*ix2-ix1*ix1
          if( idet.eq.0 )go to 10
          a=(iy1*ix2-ix1*ixy)/float(idet)
          b=(n*ixy-ix1*iy1)/float(idet)
          iy1=a+b*ixx(1)
          iy2=a+b*ixx(n)
          if( iy1*iy2.ge.jconst )go to 10
          sum=0.
          do i=1,n
             sum=sum+(a+b*ixx(i)-iyy(i))**2
          enddo
          isum=sqrt(sum)/n
          if( isum.gt.idev )go to 10

          nn=nn+1
          if( nn.gt.loop )go to 20
          jwork(nn)=i1-(n-1)*iy1/(iy2-iy1)
          if( jwork(nn).gt.nb )jwork(nn)=jwork(nn)-nb
          if( jwork(nn).lt.1 )jwork(nn)=jwork(nn)+nb
 10    continue
 20    continue
       if( nn.eq.0 )return
       n=1
       jz(n)=jwork(1)
       if( nn.ge.2 )then
           do i=2,nn
              if( abs(jwork(i-1)-jwork(i)).lt.9 )then
                  jz(n)=0.5*(jwork(i-1)+jwork(i))
              else
                 n=n+1
                 jz(n)=jwork(i)
              endif
           enddo
       endif
       jz1=jz(1)
       ier=0

       return
       end

       subroutine find_zero(nb,idat,miss,jz1,jz2,ier)
c***********************************************************************
c Subroutine/Function : find_zero
c
c Usage :
c    call find_zero(nb,idat,miss,jz1,jz2,ier)
c
c Description      : To find the azimuthal position of zero value.
c
c Arguments :
c  I/O/W   name,      type,       description
c    I     nb         integer     total ray number of input sweep scan layer.
c    I     idat(nb)   int array   radial velocity.
c    I     miss       integer     missing or bad value.
c    O     jz1        integer     the first index of zero value.
c                                 =0, if no zero value
c    O     jz2        integer     the second index of zero value.
c                                 =0, if no zero value
c    O     ier        integer     error message. (=0, success; =1, failure)
c
c Modules Called : none
c***********************************************************************
       implicit none
       integer nb,miss,jz1,jz2,ier
       integer idat(nb)

       integer loop
       parameter( loop=10 )
       integer ixx(21),iyy(21),jwork(loop),jz(loop)
       integer jconst,idev,i,ii,k,kk,nn,i1,i2,n,ix1,ix2,iy1,iy2,ixy,idet
       integer isum,ivmax,nb2,nb23,nb49
       real a,b,sum

       jz1=0
       jz2=0
       ier=1
       jconst=-40000
       idev=200
       ivmax=600
       nb2=nb/2
       nb23=2*nb/3
       nb49=4*nb/9
       do i=1,loop
          jz(i)=0
          jwork(i)=0
       enddo

       kk=nb/10
       nn=0
       do 10 k=1,kk
          i1=10*(k-1)-10
          i2=10*(k-1)+10
          n=0
          do i=i1,i2,1
             ii=i
             if( ii.lt.1 )ii=ii+nb
             if( ii.gt.nb )ii=ii-nb
             if( (idat(ii).ne.miss).and.(abs(idat(ii)).lt.ivmax) )then
                 n=n+1
                 ixx(n)=i-i1+1
                 iyy(n)=idat(ii)
             endif
          enddo
          if( n.lt.8 )go to 10
          ix1=0
          ix2=0
          iy1=0
          ixy=0
          do i=1,n
             ix1=ix1+ixx(i)
             ix2=ix2+ixx(i)*ixx(i)
             iy1=iy1+iyy(i)
             ixy=ixy+ixx(i)*iyy(i)
          enddo
          idet=n*ix2-ix1*ix1
          if( idet.eq.0 )go to 10
          a=(iy1*ix2-ix1*ixy)/float(idet)
          b=(n*ixy-ix1*iy1)/float(idet)
          iy1=a+b*ixx(1)
          iy2=a+b*ixx(n)
          if( iy1*iy2.ge.jconst )go to 10
          sum=0.
          do i=1,n
             sum=sum+(a+b*ixx(i)-iyy(i))**2
          enddo
          isum=sqrt(sum)/n
          if( isum.gt.idev )go to 10

          nn=nn+1
          if( nn.gt.loop )go to 20
          jwork(nn)=i1-(n-1)*iy1/(iy2-iy1)
          if( jwork(nn).gt.nb )jwork(nn)=jwork(nn)-nb
          if( jwork(nn).lt.1 )jwork(nn)=jwork(nn)+nb
 10    continue
 20    continue
       if( nn.eq.0 )return
       n=1
       jz(n)=jwork(1)
       if( nn.ge.2 )then
           do i=2,nn
              if( abs(jwork(i-1)-jwork(i)).lt.9 )then
                  jz(n)=0.5*(jwork(i-1)+jwork(i))
              else
                 n=n+1
                 jz(n)=jwork(i)
              endif
           enddo
       endif
       if( n.eq.1 )return
       if( n.eq.2 )then
           idev=abs(jz(1)-jz(2))
           if( idev.gt.nb23 )idev=nb-idev
           if( idev.gt.nb49 )then
               jz1=jz(1)
               jz2=jz(2)
               ier=0 
           endif
           return
       endif
       do k=1,n-1
          kk=k+1
          idev=abs(jz(k)-jz(kk))
          if( idev.gt.nb23 )idev=nb-idev
          isum=abs(idev-nb2)
          do i=k+1,n
             idev=abs(jz(k)-jz(i))
             if( idev.gt.nb23 )idev=nb-idev
             if( isum.gt.(abs(idev-nb2)) )then
                 isum=abs(idev-nb2)
                 kk=i
             endif 
          enddo
          idev=abs(jz(k)-jz(kk))
          if( idev.gt.nb23 )idev=nb-idev
          if( idev.gt.nb49 )then
              jz1=jz(k)
              jz2=jz(kk)
              ier=0
              return
          endif
       enddo
       return
       end

       subroutine ray_unf1(ng,nb,idat,ind,nquist,shear,num,id,ngp,idst
     1                    ,istatus)
c***********************************************************************
c description      : To unfold data by along ray direction mothod.
c                                                                       
c I/O parameters   :                                                    
c  I/O/W   name,            type,      description
c    I     ng               integer    first dimension of input data array.
c    I     nb               integer    second dimension of input data array.
c   I/O    idat(ng,nb)    int array  input data array
c   I/O    ind(ng,nb)     int array  group index.
c    I     nquist           integer    nquist velocity.
c    I     shear            integer    shear velocity.
c    I     num              integer    input dimension.
c    I     id(num)          int array  group index.
c   I/O    ngp(num)         int array  size of group index.
c    I     idst             integer    starting group index.
c    O     istatus          integer    =1, unfold data.
c
c Called Function: curvefitting
c***********************************************************************
       implicit none
       integer ng,nb
       integer idat(ng,nb),ind(ng,nb)
       integer nquist,shear,num
       integer id(num),ngp(num)
       integer idst,istatus

       integer k,kst
       integer nquist2,shearr
       real shearp,shearn,shearp3,shearn3
       integer nog,nob
       parameter( nog=940,nob=430 )
       real dat(nog),wkdat(nog,nob)

       integer i,j,n
       integer nmax
       parameter( nmax=5 )
       real value,dev,base

       istatus=0 
       kst=1
       do k=1,num
          if( id(k).eq.idst )kst=k
       enddo 
       nquist2=2*nquist
       shearp=0.012*nquist
       shearn=-shearp
       shearp3=0.032*nquist
       shearn3=-shearp3

       shearr=2*shear
       do j=1,nb
       do i=1,ng
          wkdat(i,j)=-999.0
       enddo
       enddo

       do j=1,nb
          n=0
          do i=1,ng
             dat(i)=-999.0
             if( ind(i,j).eq.idst )then
                 n=n+1
                 dat(n)=0.01*idat(i,j)
             endif
          enddo
          if( n.gt.nmax )then
              call curvefitting(ng,dat)
              do i=1,ng
                 wkdat(i,j)=dat(i)
              enddo
          endif
       enddo

       do 10 k=1,num
          if( id(k).eq.idst )go to 10
          n=0
          value=0.
          base=0.
          do i=1,ng
          do j=3,nb-2,1
             if( (ind(i,j).eq.id(k)).and.
     1           (wkdat(i,j).gt.-300.) )then
                 n=n+1
                 value=value+0.01*idat(i,j)
                 base=base+wkdat(i,j)
             endif 
          enddo
          enddo
          if( n.lt.nmax )go to 10
          value=value/float(n)
          base=base/float(n)
          dev=value-base
          if( (dev.gt.shearp).and.(value.gt.2.) )then
              istatus=1
              if( dev.gt.shearp3 )then
                  do j=1,nb
                  do i=1,ng
                     if( ind(i,j).eq.id(k) )then
                         idat(i,j)=idat(i,j)-2*nquist2
                         ind(i,j)=idst
                     endif
                  enddo
                  enddo
              else 
                  do j=1,nb
                  do i=1,ng
                     if( ind(i,j).eq.id(k) )then
                         idat(i,j)=idat(i,j)-nquist2
                         ind(i,j)=idst
                     endif
                  enddo
                  enddo
              endif
              ngp(kst)=ngp(kst)+ngp(k)
              ngp(k)=0
              go to 10
          endif
          if( (dev.lt.shearn).and.(value.lt.-2.) )then
              istatus=1
              if( dev.lt.shearn3 )then
                  do j=1,nb
                  do i=1,ng
                     if( ind(i,j).eq.id(k) )then
                         idat(i,j)=idat(i,j)+2*nquist2
                         ind(i,j)=idst
                     endif
                  enddo
                  enddo
              else
                  do j=1,nb
                  do i=1,ng
                     if( ind(i,j).eq.id(k) )then
                         idat(i,j)=idat(i,j)+nquist2
                         ind(i,j)=idst
                     endif
                  enddo
                  enddo
              endif
              ngp(kst)=ngp(kst)+ngp(k)
              ngp(k)=0
              go to 10 
          endif 
          if( abs(100.*dev).lt.shearr )then
              istatus=1
              do j=1,nb
              do i=1,ng
                 if( ind(i,j).eq.id(k) )then
                     ind(i,j)=idst
                 endif
              enddo
              enddo
              ngp(kst)=ngp(kst)+ngp(k)
              ngp(k)=0
          endif
 10    continue

       do j=3,nb-2
       do i=1,ng
          if( (ind(i,j).ne.-2).and.(ind(i,j).ne.idst).and.
     1        (wkdat(i,j).gt.-300.) )then
              dev=0.01*idat(i,j)-wkdat(i,j)
              if( (dev.gt.shearp).and.(idat(i,j).gt.200) )then
                  if( dev.gt.shearp3 )then
                      idat(i,j)=idat(i,j)-2*nquist2
                  else
                      idat(i,j)=idat(i,j)-nquist2
                  endif
                  ind(i,j)=idst
                  ngp(kst)=ngp(kst)+1
              endif
              if( (dev.lt.shearn).and.(idat(i,j).lt.-200) )then
                  if( dev.lt.shearn3 )then
                      idat(i,j)=idat(i,j)+2*nquist2
                  else
                      idat(i,j)=idat(i,j)+nquist2
                  endif
                  ind(i,j)=idst
                  ngp(kst)=ngp(kst)+1
              endif
              if( abs(100*dev).lt.shearr )then
                  ind(i,j)=idst
                  ngp(kst)=ngp(kst)+1
              endif
          endif
       enddo
       enddo

       return
       end

       subroutine ray_unf(ng,nb,idat,ind,nquist,shear,num,id,ngp,idst
     1                   ,istatus)
c***********************************************************************
c description      : To unfold data by along ray direction mothod.
c                                                                       
c I/O parameters   :                                                    
c  I/O/W   name,            type,      description
c    I     ng               integer    first dimension of input data array.
c    I     nb               integer    second dimension of input data array.
c   I/O    idat(ng,nb)    int array  input data array
c   I/O    ind(ng,nb)     int array  group index.
c    I     nquist           integer    nquist velocity.
c    I     shear            integer    shear velocity.
c    I     num              integer    input dimension.
c    I     id(num)          int array  group index.
c   I/O    ngp(num)         int array  size of group index.
c    I     idst             integer    starting group index.
c    O     istatus          integer    =1, unfold data.
c
c Called Function: curvefitting
c***********************************************************************
       implicit none
       integer ng,nb
       integer idat(ng,nb),ind(ng,nb)
       integer nquist,shear,num
       integer id(num),ngp(num)
       integer idst,istatus

       integer k,kst
       integer nquist2,shearr
       real shearp,shearn,shearp3,shearn3
       integer nog,nob
       parameter( nog=940,nob=430 )
       real dat(nog),wkdat(nog,nob)

       integer i,j,n
       integer nmax
       parameter( nmax=5 )
       real value,dev,base

       istatus=0 
       kst=1
       do k=1,num
          if( id(k).eq.idst )kst=k
       enddo 
       nquist2=2*nquist
       shearp=0.012*nquist
       shearn=-shearp
       shearp3=0.032*nquist
       shearn3=-shearp3
       shearr=2*shear
       do j=1,nb
       do i=1,ng
          wkdat(i,j)=-999.0
       enddo
       enddo

       do j=1,nb
          n=0
          do i=1,ng
             dat(i)=-999.0
             if( ind(i,j).eq.idst )then
                 n=n+1
                 dat(n)=0.01*idat(i,j)
             endif
          enddo
          if( n.gt.nmax )then
              call curvefitting(ng,dat)
              do i=1,ng
                 wkdat(i,j)=dat(i)
              enddo
          endif
       enddo

       do 10 k=1,num
          if( id(k).eq.idst )go to 10
          n=0
          value=0.
          base=0.
          do i=1,ng
          do j=3,nb-2,1
             if( (ind(i,j).eq.id(k)).and.
     1           (wkdat(i,j).gt.-300.) )then
                 n=n+1
                 value=value+0.01*idat(i,j)
                 base=base+wkdat(i,j)
             endif 
          enddo
          enddo
          if( n.lt.nmax )go to 10
          value=value/float(n)
          base=base/float(n)
          dev=value-base
          if( (dev.gt.shearp).and.(value.gt.2.) )then
              istatus=1
              if( dev.gt.shearp3 )then
                  do j=1,nb
                  do i=1,ng
                     if( ind(i,j).eq.id(k) )then
                         idat(i,j)=idat(i,j)-2*nquist2
                         ind(i,j)=idst
                     endif
                  enddo
                  enddo
              else  
                  do j=1,nb
                  do i=1,ng
                     if( ind(i,j).eq.id(k) )then
                         idat(i,j)=idat(i,j)-nquist2
                         ind(i,j)=idst
                     endif
                  enddo
                  enddo
              endif
              ngp(kst)=ngp(kst)+ngp(k)
              ngp(k)=0
              go to 10
          endif
          if( (dev.lt.shearn).and.(value.lt.-2.) )then
              istatus=1
              if( dev.lt.shearn3 )then
                  do j=1,nb
                  do i=1,ng
                     if( ind(i,j).eq.id(k) )then
                         idat(i,j)=idat(i,j)+2*nquist2
                         ind(i,j)=idst
                     endif
                  enddo
                  enddo
              else
                  do j=1,nb
                  do i=1,ng
                     if( ind(i,j).eq.id(k) )then
                         idat(i,j)=idat(i,j)+nquist2
                         ind(i,j)=idst
                     endif
                  enddo
                  enddo
              endif
              ngp(kst)=ngp(kst)+ngp(k)
              ngp(k)=0
              go to 10 
          endif 
          if( abs(100.*dev).lt.shearr )then
              istatus=1
              do j=1,nb
              do i=1,ng
                 if( ind(i,j).eq.id(k) )then
                     ind(i,j)=idst
                 endif
              enddo
              enddo
              ngp(kst)=ngp(kst)+ngp(k)
              ngp(k)=0
          endif
 10    continue

       return
       end

       subroutine curvefitting(ng,dat)
c***********************************************************************
c description      : To get data by curving fitting.
c                                                                       
c I/O parameters   :                                                    
c  I/O/W   name,            type,      description
c    I     ng               integer    input dimension.
c   I/O    dat(ng)          real array input data array.
c
c Called Function: none
c***********************************************************************
       implicit none
       integer ng
       real dat(ng)
       integer i,ia,ib,n,ii
       real sumv,va,vb,xa,xb,xx

       do i=1,ng
          if( dat(i).gt.-900. )then
              ia=i
              go to 10
          endif
       enddo
 10    continue
       do i=ng,1,-1
          if( dat(i).gt.-900. )then
              ib=i
              go to 20
          endif
       enddo
 20    continue
       if( ia.ge.ib )return

       n=0
       sumv=0.
       do i=ia,ib,1
          if( n.gt.4 )then
              va=sumv/float(n)
              go to 30
          endif
          if( dat(i).gt.-900. )then
              n=n+1
              sumv=sumv+dat(i)
          endif
       enddo
       return
 30    continue
       do i=1,ia-1
          dat(i)=va
       enddo
       n=0
       sumv=0.
       do i=ib,ia,-1
          if( n.gt.4 )then
              vb=sumv/float(n)
              go to 40
          endif
          if( dat(i).gt.-900. )then
              n=n+1
              sumv=sumv+dat(i)
          endif
       enddo
       return
 40    continue
       do i=ib+1,ng
          dat(i)=vb
       enddo

       do 70 i=ia+1,ib-1,1
          if( dat(i).lt.-900. )then
              do ii=i-1,ia,-1
                 if( dat(ii).gt.-900. )then
                     va=dat(ii)
                     xa=ii
                     go to 50
                 endif
              enddo
              go to 70
 50           continue
              do ii=i+1,ib,1
                 if( dat(ii).gt.-900. )then
                     vb=dat(ii)
                     xb=ii
                     go to 60
                 endif
              enddo
              go to 70
 60           continue 
              xx=i
              dat(i)=va+(xx-xa)*(vb-va)/(xb-xa)
          endif
 70    continue

       return
       end
 
       subroutine azim_unf(ng,nb,idat,ind,nquist,shear,num,id,ngp,idst
     1                    ,istatus)
c***********************************************************************
c description      : To unfold data by along azimuthal direction mothod.
c                                                                       
c I/O parameters   :                                                    
c  I/O/W   name,            type,      description
c    I     ng               integer    first dimension of input data array.
c    I     nb               integer    second dimension of input data array.
c   I/O    idat(ng,nb)    int array  input data array
c   I/O    ind(ng,nb)     int array  group index.
c    I     nquist           integer    nquist velocity.
c    I     shear            integer    shear velocity.
c    I     num              integer    input dimension.
c    I     id(num)          int array  group index.
c   I/O    ngp(num)         int array  size of group index.
c    I     idst             integer    starting group index.
c    O     istatus          integer    =1, unfold data.
c
c Called Function: sov_coe7
c***********************************************************************
       implicit none
       integer ng,nb
       integer idat(ng,nb),ind(ng,nb)
       integer nquist,shear,num
       integer id(num),ngp(num)
       integer idst,istatus

       integer k,kst
       integer nquist2,shearr
       real shearp,shearn,shearp3,shearn3
       real pi
       parameter( pi=3.1415926 ) 
       integer nog,nob
       parameter( nog=940,nob=430 )
       real ang(nob),dat(nob),th(nob),cc(7)
       real wkdat(nog,nob)

       integer i,j,n,m,nlimit,mlimit,ier
       parameter( nlimit=10,mlimit=60 )
       integer nmax
       parameter( nmax=5 )
       real az,value,dev,base

       istatus=0
       kst=1
       do k=1,num
          if( id(k).eq.idst )kst=k
       enddo 
       nquist2=2*nquist
       shearp=0.012*nquist
       shearn=-shearp
       shearp3=0.032*nquist
       shearn3=-shearp3
       shearr=2*shear
       do j=1,nb
          ang(j)=(j-3)*360./float(nb-4)
       enddo
       do j=1,nb
       do i=1,ng
          wkdat(i,j)=-999.0
       enddo
       enddo
       do i=1,ng
          n=0
          m=0
          do j=3,nb-2,1
             if( m.gt.mlimit )go to 10
             if( ind(i,j).eq.idst )then
                 n=n+1
                 dat(n)=0.01*idat(i,j)
                 th(n)=ang(j)
                 m=0
             else
                 m=m+1
             endif
          enddo
          if( n.lt.nlimit )go to 10
          call sov_coe7(n,th,dat,cc,ier)
          if( ier.eq.1 )go to 10
          do j=3,nb-2,1
             if( ind(i,j).eq.idst )then
                 wkdat(i,j)=0.01*idat(i,j)
             else
                 if( ind(i,j).ne.-2 )then
                     az=ang(j)*pi/180.
                     wkdat(i,j)=cc(1)+cc(2)*cos(az)+cc(3)*sin(az)
     1                         +cc(4)*cos(2.*az)+cc(5)*sin(2.*az)
     2                         +cc(6)*cos(3.*az)+cc(7)*sin(3.*az)
                 endif
             endif
          enddo
 10       continue
       enddo

       do 20 k=1,num
          if( id(k).eq.idst )go to 20
          n=0
          value=0.
          base=0.
          do i=1,ng
          do j=3,nb-2,1
             if( (ind(i,j).eq.id(k)).and.
     1           (wkdat(i,j).gt.-300.) )then
                 n=n+1
                 value=value+0.01*idat(i,j)
                 base=base+wkdat(i,j)
             endif 
          enddo
          enddo
          if( n.lt.nmax )go to 20
          value=value/float(n)
          base=base/float(n)
          dev=value-base
          if( dev.gt.shearp )then
              istatus=1
              if( dev.gt.shearp3 )then
                  do j=1,nb
                  do i=1,ng
                     if( ind(i,j).eq.id(k) )then
                         idat(i,j)=idat(i,j)-2*nquist2
                         ind(i,j)=idst
                     endif
                  enddo
                  enddo
              else
                  do j=1,nb
                  do i=1,ng
                     if( ind(i,j).eq.id(k) )then
                         idat(i,j)=idat(i,j)-nquist2
                         ind(i,j)=idst
                     endif
                  enddo
                  enddo
              endif
              ngp(kst)=ngp(kst)+ngp(k)
              ngp(k)=0
              go to 20
          endif
          if( dev.lt.shearn )then
              istatus=1
              if( dev.lt.shearn3 )then
                  do j=1,nb
                  do i=1,ng
                     if( ind(i,j).eq.id(k) )then
                         idat(i,j)=idat(i,j)+2*nquist2
                         ind(i,j)=idst
                     endif
                  enddo
                  enddo
              else
                  do j=1,nb
                  do i=1,ng
                     if( ind(i,j).eq.id(k) )then
                         idat(i,j)=idat(i,j)+nquist2
                         ind(i,j)=idst
                     endif
                  enddo
                  enddo
              endif
              ngp(kst)=ngp(kst)+ngp(k)
              ngp(k)=0
              go to 20 
          endif 
          if( abs(100.*dev).lt.shearr )then
              istatus=1
              do j=1,nb
              do i=1,ng
                 if( ind(i,j).eq.id(k) )then
                     ind(i,j)=idst
                 endif
              enddo
              enddo
              ngp(kst)=ngp(kst)+ngp(k)
              ngp(k)=0
          endif
 20    continue

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
c***********************************************************************
      implicit none
      integer n,ier,max
      parameter( max=430 )
      real th(n),vr(n),cc(7),a(7,7),y(7)
      real cs1(max),sn1(max),cs2(max),sn2(max),cs3(max),sn3(max)
      integer i,j,k
      real pi,thi

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
c called fun./sub. : det7
c***********************************************************************
        implicit none
        integer ier
        real a(7,7),y(7),x(7),b(7,7),d(7)
        integer i,j
        real det

        call det7(a,det)
        ier=0 
        if( abs(det).lt.0.0000001 )then
            ier=1
            return
        endif
        do j=1,7
        do i=1,7
           b(i,j)=a(i,j)
        enddo
        enddo 
        do i=1,7
           b(1,i)=y(i)
        enddo
        call det7(b,d(1))
        do i=1,7
           b(1,i)=a(1,i)
           b(2,i)=y(i)
        enddo
        call det7(b,d(2))
        do i=1,7
           b(2,i)=a(2,i)
           b(3,i)=y(i)
        enddo
        call det7(b,d(3))
        do i=1,7
           b(3,i)=a(3,i)
           b(4,i)=y(i)
        enddo
        call det7(b,d(4))
        do i=1,7
           b(4,i)=a(4,i)
           b(5,i)=y(i)
        enddo
        call det7(b,d(5))
        do i=1,7
           b(5,i)=a(5,i)
           b(6,i)=y(i)
        enddo
        call det7(b,d(6))
        do i=1,7
           b(6,i)=a(6,i)
           b(7,i)=y(i)
        enddo
        call det7(b,d(7))
        do i=1,7
           x(i)=d(i)/det
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
c***********************************************************************
      implicit none
      real a(7,7),det,b(6,6)
      integer i,j,i1,j1
      real d11,d21,d31,d41,d51,d61,d71

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
c***********************************************************************
      implicit none
      real a(6,6),det,b(5,5)
      integer i,j,i1,j1
      real d11,d21,d31,d41,d51,d61

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
c***********************************************************************
      implicit none
      real a(5,5),det,b(4,4)
      integer i,j,i1,j1
      real d11,d21,d31,d41,d51

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
c***********************************************************************
      real a(4,4),det
      real a11,a21,a31,a41

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

       subroutine near_unf(ng,nb,ivad,idat,ibd,ind,nquist,shear
     1                    ,num,id,ngp,idst,istatus)
c***********************************************************************
c description      : To unfold data by near mothod.
c                                                                       
c I/O parameters   :                                                    
c  I/O/W   name,            type,      description
c    I     ng               integer    first dimension of input data array.
c    I     nb               integer    second dimension of input data array.
c    I     ivad             integer    ppi scanning message.
c                                      =0, is 360 degrees scan;  =1, isnot
c   I/O    idat(ng,nb)      int array  input data array
c    I     ibd(ng,nb)       int array  boundary index.
c   I/O    ind(ng,nb)       int array  group index.
c    I     nquist           integer    nquist velocity.
c    I     shear            integer    shear velocity.
c    I     num              integer    input dimension.
c    I     id(num)          int array  group index.
c   I/O    ngp(num)         int array  size of group index.
c    I     idst             integer    starting group index.
c    O     istatus          integer    =0, no unfold dat; =1, unfold data.
c
c Called Function: none
c***********************************************************************
       implicit none
       integer ng,nb,ivad
       integer idat(ng,nb),ibd(ng,nb),ind(ng,nb)
       integer nquist,num,shear
       integer id(num),ngp(num)
       integer idst
       integer istatus

       integer kst
       integer nmax
       parameter( nmax=10000 )
       integer numsh(nmax),numpo(nmax),numne(nmax)
       integer numpo3(nmax),numne3(nmax)
       integer numid
       integer shearp,shearn,shearp3,shearn3,ivalue
       integer i,j,ii,jj,k,n
      
       istatus=0
       if( num.le.1 )return 
       numid=id(1)
       kst=1
       do k=2,num
          if( numid.lt.id(k) )numid=id(k)
          if( id(k).eq.idst )kst=k
       enddo
       if( numid.gt.nmax )then
           print*,'Please modify parameter: nmax  in get_unfindx.'
           print*,'numid: ',numid,' nmax: ',nmax
           return
       endif

       shearp=1.6*nquist+0.5
       shearp3=3.2*nquist 
       shearn=-shearp
       shearn3=-shearp3
       do i=1,numid
          numsh(i)=0
          numpo(i)=0
          numne(i)=0
          numpo3(i)=0
          numne3(i)=0
       enddo
       do 10 j=2,nb-1,1
       do 10 i=2,ng-1,1
          if( (ibd(i,j).ne.1).or.(ind(i,j).ne.idst) )go to 10
          if( (ibd(i+1,j).eq.1).and.(ind(i,j).ne.ind(i+1,j)) )then
              ii=ind(i+1,j)
              ivalue=idat(i+1,j)-idat(i,j)
              if( ivalue.gt.shearp )then
                  if( ivalue.gt.shearp3 )then 
                      numpo3(ii)=numpo3(ii)+1
                  else
                      numpo(ii)=numpo(ii)+1
                  endif
              else
                  if( ivalue.lt.shearn )then
                      if( ivalue.lt.shearn3 )then
                          numne3(ii)=numne3(ii)+1
                      else
                          numne(ii)=numne(ii)+1
                      endif
                  else
                      if( abs(ivalue).lt.shear )numsh(ii)=numsh(ii)+1
                  endif
              endif
          endif
          if( (ibd(i-1,j).eq.1).and.(ind(i,j).ne.ind(i-1,j)) )then
              ii=ind(i-1,j)
              ivalue=idat(i-1,j)-idat(i,j)
              if( ivalue.gt.shearp )then
                  if( ivalue.gt.shearp3 )then
                      numpo3(ii)=numpo3(ii)+1
                  else
                      numpo(ii)=numpo(ii)+1
                  endif
              else
                  if( ivalue.lt.shearn )then
                      if( ivalue.lt.shearn3 )then
                          numne3(ii)=numne3(ii)+1
                      else
                          numne(ii)=numne(ii)+1
                      endif
                  else
                      if( abs(ivalue).lt.shear )numsh(ii)=numsh(ii)+1
                  endif
              endif
          endif
          if( (ibd(i,j+1).eq.1).and.(ind(i,j).ne.ind(i,j+1)) )then
              ii=ind(i,j+1)
              ivalue=idat(i,j+1)-idat(i,j)
              if( ivalue.gt.shearp )then
                  if( ivalue.gt.shearp3 )then
                      numpo3(ii)=numpo3(ii)+1
                  else
                      numpo(ii)=numpo(ii)+1
                  endif
              else
                  if( ivalue.lt.shearn )then
                      if( ivalue.lt.shearn3 )then
                          numne3(ii)=numne3(ii)+1
                      else
                          numne(ii)=numne(ii)+1
                      endif
                  else
                      if( abs(ivalue).lt.shear )numsh(ii)=numsh(ii)+1
                  endif
              endif
          endif
          if( (ibd(i,j-1).eq.1).and.(ind(i,j).ne.ind(i,j-1)) )then
              ii=ind(i,j-1)
              ivalue=idat(i,j-1)-idat(i,j)
              if( ivalue.gt.shearp )then
                  if( ivalue.gt.shearp3 )then
                      numpo3(ii)=numpo3(ii)+1
                  else
                      numpo(ii)=numpo(ii)+1
                  endif
              else
                  if( ivalue.lt.shearn )then
                      if( ivalue.lt.shearn3 )then
                          numne3(ii)=numne3(ii)+1
                      else
                          numne(ii)=numne(ii)+1
                      endif
                  else
                      if( abs(ivalue).lt.shear )numsh(ii)=numsh(ii)+1
                  endif
              endif
          endif
          if( (ibd(i-1,j-1).eq.1).and.(ind(i,j).ne.ind(i-1,j-1)) )then
              ii=ind(i-1,j-1)
              ivalue=idat(i-1,j-1)-idat(i,j)
              if( ivalue.gt.shearp )then
                  if( ivalue.gt.shearp3 )then
                      numpo3(ii)=numpo3(ii)+1
                  else
                      numpo(ii)=numpo(ii)+1
                  endif
              else
                  if( ivalue.lt.shearn )then
                      if( ivalue.lt.shearn3 )then
                          numne3(ii)=numne3(ii)+1
                      else
                          numne(ii)=numne(ii)+1
                      endif
                  else
                      if( abs(ivalue).lt.shear )numsh(ii)=numsh(ii)+1
                  endif
              endif
          endif
          if( (ibd(i+1,j-1).eq.1).and.(ind(i,j).ne.ind(i+1,j-1)) )then
              ii=ind(i+1,j-1)
              ivalue=idat(i+1,j-1)-idat(i,j)
              if( ivalue.gt.shearp )then
                  if( ivalue.gt.shearp3 )then
                      numpo3(ii)=numpo3(ii)+1
                  else
                      numpo(ii)=numpo(ii)+1
                  endif
              else
                  if( ivalue.lt.shearn )then
                      if( ivalue.lt.shearn3 )then
                          numne3(ii)=numne3(ii)+1
                      else
                          numne(ii)=numne(ii)+1
                      endif
                  else
                      if( abs(ivalue).lt.shear )numsh(ii)=numsh(ii)+1
                  endif
              endif
          endif
          if( (ibd(i+1,j+1).eq.1).and.(ind(i,j).ne.ind(i+1,j+1)) )then
              ii=ind(i+1,j+1)
              ivalue=idat(i+1,j+1)-idat(i,j)
              if( ivalue.gt.shearp )then
                  if( ivalue.gt.shearp3 )then
                      numpo3(ii)=numpo3(ii)+1
                  else
                      numpo(ii)=numpo(ii)+1
                  endif
              else
                  if( ivalue.lt.shearn )then
                      if( ivalue.lt.shearn3 )then
                          numne3(ii)=numne3(ii)+1
                      else 
                          numne(ii)=numne(ii)+1
                      endif
                  else
                      if( abs(ivalue).lt.shear )numsh(ii)=numsh(ii)+1
                  endif
              endif
          endif
          if( (ibd(i-1,j+1).eq.1).and.(ind(i,j).ne.ind(i-1,j+1)) )then
              ii=ind(i-1,j+1)
              ivalue=idat(i-1,j+1)-idat(i,j)
              if( ivalue.gt.shearp )then
                  if( ivalue.gt.shearp3 )then
                      numpo3(ii)=numpo3(ii)+1
                  else
                      numpo(ii)=numpo(ii)+1
                  endif 
              else
                  if( ivalue.lt.shearn )then
                      if( ivalue.lt.shearn3 )then
                          numne3(ii)=numne3(ii)+1
                      else
                          numne(ii)=numne(ii)+1
                      endif
                  else
                      if( abs(ivalue).lt.shear )numsh(ii)=numsh(ii)+1
                  endif
              endif
          endif
 10    continue
       if( ivad.eq.0 )then
           do i=1,ng
              if( (ind(i,3).eq.idst).and.(ind(i,nb-2).gt.0).and.
     1            (ind(i,3).ne.ind(i,nb-2)) )then
                  ii=ind(i,nb-2)
                  ivalue=idat(i,nb-2)-idat(i,3)
                  if( ivalue.gt.shearp )then
                      if( ivalue.gt.shearp3 )then 
                          numpo3(ii)=numpo3(ii)+1
                      else
                          numpo(ii)=numpo(ii)+1
                      endif
                  else
                      if( ivalue.lt.shearn )then
                          if( ivalue.lt.shearn3 )then
                              numne3(ii)=numne3(ii)+1
                          else
                              numne(ii)=numne(ii)+1
                          endif
                      else
                      if( abs(ivalue).lt.shear )numsh(ii)=numsh(ii)+1
                      endif
                  endif
              endif
              if( (ind(i,nb-2).eq.idst).and.(ind(i,3).gt.0).and.
     1            (ind(i,3).ne.ind(i,nb-2)) )then
                  ii=ind(i,3)
                  ivalue=idat(i,3)-idat(i,nb-2)
                  if( ivalue.gt.shearp )then
                      if( ivalue.gt.shearp3 )then 
                          numpo3(ii)=numpo3(ii)+1
                      else
                          numpo(ii)=numpo(ii)+1
                      endif
                  else
                      if( ivalue.lt.shearn )then
                          if( ivalue.lt.shearn3 )then
                              numne3(ii)=numne3(ii)+1
                          else
                              numne(ii)=numne(ii)+1
                          endif
                      else
                      if( abs(ivalue).lt.shear )numsh(ii)=numsh(ii)+1
                      endif
                  endif
              endif
           enddo 
       endif

       n=0
       do 20 k=1,num
          if( k.eq.kst )go to 20
          ii=id(k)
          if( (numpo3(ii).gt.8).and.(numpo3(ii).gt.numpo(ii)) )then
              do j=1,nb
              do i=1,ng
                 if( ind(i,j).eq.id(k) )then
                     ind(i,j)=idst 
                     idat(i,j)=idat(i,j)-4*nquist
                     n=n+1
                 endif
              enddo
              enddo
              ngp(kst)=ngp(kst)+ngp(k)
              ngp(k)=0
              go to 20
          endif
          if( (numpo(ii).gt.8).and.(numpo(ii).gt.numsh(ii)) )then
              do j=1,nb
              do i=1,ng
                 if( ind(i,j).eq.id(k) )then
                     ind(i,j)=idst 
                     idat(i,j)=idat(i,j)-2*nquist
                     n=n+1
                 endif
              enddo
              enddo
              ngp(kst)=ngp(kst)+ngp(k)
              ngp(k)=0
              go to 20
          endif
          if( (numne3(ii).gt.8).and.(numne3(ii).gt.numne(ii)) )then
              do j=1,nb
              do i=1,ng
                 if( ind(i,j).eq.id(k) )then
                     ind(i,j)=idst
                     idat(i,j)=idat(i,j)+4*nquist
                     n=n+1
                 endif
              enddo
              enddo
              ngp(kst)=ngp(kst)+ngp(k)
              ngp(k)=0
              go to 20
          endif
          if( (numne(ii).gt.8).and.(numne(ii).gt.numsh(ii)) )then
              do j=1,nb
              do i=1,ng
                 if( ind(i,j).eq.id(k) )then
                     ind(i,j)=idst
                     idat(i,j)=idat(i,j)+2*nquist
                     n=n+1
                 endif
              enddo
              enddo
              ngp(kst)=ngp(kst)+ngp(k)
              ngp(k)=0
              go to 20
          endif
          if( numsh(ii).gt.10 )then
              do j=1,nb
              do i=1,ng
                 if( ind(i,j).eq.id(k) )then
                     ind(i,j)=idst
                     n=n+1
                 endif
              enddo
              enddo
              ngp(kst)=ngp(kst)+ngp(k)
              ngp(k)=0
          endif
 20    continue
       if( n.eq.0 )return
       istatus=1

       return
       end

       subroutine mergegp(ng,nb,idat,ibd,ind,shear,num,id,ngp)
c***********************************************************************
c description      : To merge group index.
c                                                                       
c I/O parameters   :                                                    
c  I/O/W   name,            type,      description
c    I     ng               integer    first dimension of input data array.
c    I     nb               integer    second dimension of input data array.
c    I     idat(ng,ng)      int array  input data array
c    I     ibd(ng,nb)       int array  boundary index.
c   I/O    ind(ng,ng)       int array  group index.
c    I     shear            integer    wind shear.
c    I     num              integer    input dimension.
c    I     id(num)          int array  group index.
c    I     ngp(num)         int array  size of group index.
c
c Called Function: none
c***********************************************************************
       implicit none
       integer ng,nb
       integer idat(ng,nb),ibd(ng,nb),ind(ng,nb)
       integer num,shear
       integer id(num),ngp(num)

       integer nmax
       parameter( nmax=10000 )
       integer numbd(nmax,nmax)
       integer i,j,ii,jj,k,kk,n
       integer nummin

       if( num.le.1 )return 
       if( num.gt.nmax )then
           print*,'Please modify parameter: nmax  in mergegp.'
           print*,'num: ',num,' nmax: ',nmax
           return
       endif
       do j=1,num
       do i=1,num
          numbd(i,j)=0
       enddo
       enddo

       do 10 j=2,nb-1,1
       do 10 i=2,ng-1,1
          if( ibd(i,j).eq.1 )then
              jj=ind(i,j) 
              if( (ibd(i+1,j).eq.1).and.(ind(i,j).ne.ind(i+1,j)).and.
     1            (abs(idat(i+1,j)-idat(i,j)).lt.shear) )then
                  ii=ind(i+1,j)
                  numbd(ii,jj)=numbd(ii,jj)+1
                  go to 10
              endif
              if( (ibd(i-1,j).eq.1).and.(ind(i,j).ne.ind(i-1,j)).and.
     1            (abs(idat(i-1,j)-idat(i,j)).lt.shear) )then
                  ii=ind(i-1,j)
                  numbd(ii,jj)=numbd(ii,jj)+1
                  go to 10
              endif
              if( (ibd(i,j+1).eq.1).and.(ind(i,j).ne.ind(i,j+1)).and.
     1            (abs(idat(i,j+1)-idat(i,j)).lt.shear) )then
                  ii=ind(i,j+1)
                  numbd(ii,jj)=numbd(ii,jj)+1
                  go to 10
              endif
              if( (ibd(i,j-1).eq.1).and.(ind(i,j).ne.ind(i,j-1)).and.
     1            (abs(idat(i,j-1)-idat(i,j)).lt.shear) )then
                  ii=ind(i,j-1)
                  numbd(ii,jj)=numbd(ii,jj)+1
                  go to 10
              endif
              if( (ibd(i-1,j-1).eq.1).and.(ind(i,j).ne.ind(i-1,j-1))
     1            .and.(abs(idat(i-1,j-1)-idat(i,j)).lt.shear) )then
                  ii=ind(i-1,j-1)
                  numbd(ii,jj)=numbd(ii,jj)+1
                  go to 10
              endif
              if( (ibd(i+1,j-1).eq.1).and.(ind(i,j).ne.ind(i+1,j-1))
     1            .and.(abs(idat(i+1,j-1)-idat(i,j)).lt.shear) )then
                  ii=ind(i+1,j-1)
                  numbd(ii,jj)=numbd(ii,jj)+1
                  go to 10
              endif
              if( (ibd(i+1,j+1).eq.1).and.(ind(i,j).ne.ind(i+1,j+1))
     1            .and.(abs(idat(i+1,j+1)-idat(i,j)).lt.shear) )then
                  ii=ind(i+1,j+1)
                  numbd(ii,jj)=numbd(ii,jj)+1
                  go to 10
              endif
              if( (ibd(i-1,j+1).eq.1).and.(ind(i,j).ne.ind(i-1,j+1))
     1            .and.(abs(idat(i-1,j+1)-idat(i,j)).lt.shear) )then
                  ii=ind(i-1,j+1)
                  numbd(ii,jj)=numbd(ii,jj)+1
                  go to 10
              endif
          endif
 10    continue

       do 100 k=num,2,-1
          do 50 kk=1,k-1,1
             ii=id(kk)
             jj=id(k)
             n=numbd(ii,jj)
             if( n.lt.8 )go to 50
             do j=1,nb
             do i=1,ng
                if( ind(i,j).eq.id(k) )then
                    ind(i,j)=id(kk)
                endif
             enddo
             enddo
             ngp(kk)=ngp(kk)+ngp(k)
             ngp(k)=0
             go to 60
 50       continue
 60       continue
 100   continue

       return
       end

       subroutine sortgp(num,id,ngp)
c***********************************************************************
c description      : To sort the size of group index.
c                                                                       
c I/O parameters   :                                                    
c  I/O/W   name,            type,      description
c    I     num              integer    input dimension.
c   I/O    id(num)          int array  group index.
c   I/O    ngp(num)         int array  size of group index.
c
c Called Function: none
c***********************************************************************
       implicit none
       integer num
       integer id(num),ngp(num)
       integer i,j,idtemp,ngptemp

       if( num.le.1 )return

       do i=1,num-1
          do j=i+1,num
             if( ngp(i).lt.ngp(j) )then
                 ngptemp=ngp(i)
                 idtemp=id(i)
                 ngp(i)=ngp(j)
                 id(i)=id(j)
                 ngp(j)=ngptemp
                 id(j)=idtemp
             endif
          enddo 
       enddo

       return
       end

       subroutine get_boundary(ng,nb,ind,ibd)
c***********************************************************************
c description      : To get boundary index.
c                                                                       
c I/O parameters   :                                                    
c  I/O/W   name,            type,      description
c    I     ng               integer    first dimension of input data array.
c    I     nb               integer    second dimension of input data array.
c    I     ind(ng,nb)       int array  group index.
c    O     ibd(ng,nb)       int array  boundary index.
c
c Called Function: none
c***********************************************************************
       implicit none
       integer ng,nb
       integer ind(ng,nb),ibd(ng,nb)
       integer i,j,ii,jj

       do j=1,nb
       do i=1,ng
          ibd(i,j)=0
       enddo
       enddo 
       do 10 j=2,nb-1,1
       do 10 i=2,ng-1,1
          if( ind(i,j).gt.0 )then
              do jj=j-1,j+1,1
              do ii=i-1,i+1,1
                 if( (ind(ii,jj).gt.0).and.
     1               (ind(i,j).ne.ind(ii,jj)) )then
                     ibd(i,j)=1
                     go to 10 
                 endif
              enddo
              enddo
          endif
 10    continue

       return
       end

       subroutine merge_bound(ng,nb,gst,ind,istatus)
c***********************************************************************
c description      : To merge group index for boundary area.
c                                                                       
c I/O parameters   :                                                    
c  I/O/W   name,            type,      description
c    I     ng               integer    first dimension of input data array.
c    I     nb               integer    second dimension of input data array.
c    I     gst              integer    number of group index.
c   I/O    ind(ng,nb)       int array  group index.
c    O     istatus          integer    =0, is done. 
c
c Called Function: none
c***********************************************************************
       implicit none
       integer ng,nb,gst,istatus
       integer ind(ng,nb)
       integer nmax
       parameter( nmax=10000 )
       integer numbd(nmax,nmax)
       integer nlimit
       parameter( nlimit=120 )

       integer nb1,nb2,nb3,i,j,ii,jj

       if( gst.le.1 )return
       if( gst.gt.nmax )then
           print*,'Error in subroutine: merge_bound'
           print*,'Please modify parameter: nmax'
           print*,'nmax: ',nmax,' gst: ',gst
           return
       endif

       istatus=0
       nb1=nb-1
       nb2=nb-2
       nb3=nb-3
       do jj=1,gst
       do ii=1,gst
          numbd(ii,jj)=0
       enddo
       enddo
       do i=1,ng
          if( (ind(i,1).gt.0).and.(ind(i,nb3).gt.0).and.
     1        (ind(i,1).ne.ind(i,nb3)) )then
              ii=ind(i,1)
              jj=ind(i,nb3)
              numbd(ii,jj)=numbd(ii,jj)+1
          endif
          if( (ind(i,2).gt.0).and.(ind(i,nb2).gt.0).and.
     1        (ind(i,2).ne.ind(i,nb2)) )then
              ii=ind(i,2)
              jj=ind(i,nb2)
              numbd(ii,jj)=numbd(ii,jj)+1
          endif
          if( (ind(i,3).gt.0).and.(ind(i,nb1).gt.0).and.
     1        (ind(i,3).ne.ind(i,nb1)) )then
              ii=ind(i,3)
              jj=ind(i,nb1)
              numbd(ii,jj)=numbd(ii,jj)+1
          endif
          if( (ind(i,4).gt.0).and.(ind(i,nb).gt.0).and.
     1        (ind(i,4).ne.ind(i,nb)) )then
              ii=ind(i,4)
              jj=ind(i,nb)
              numbd(ii,jj)=numbd(ii,jj)+1
          endif
       enddo
       do jj=1,gst
       do ii=1,gst
          if( (ii.ne.jj).and.(numbd(ii,jj).gt.nlimit) )then
              istatus=1
              do j=1,nb
              do i=1,ng
                 if( ind(i,j).eq.jj )then
                     ind(i,j)=ii
                 endif
              enddo
              enddo
          endif
       enddo
       enddo

       return
       end

       subroutine subgroup(ng,nb,ishr,ind,ist,jst,gst,ier)
c***********************************************************************
c description      : To get sub-group index.
c                                                                       
c I/O parameters   :                                                    
c  I/O/W   name,            type,      description
c    I     ng               integer    first dimension of input data array.
c    I     nb               integer    second dimension of input data array.
c    I     ishr(8,ng,nb)    int array  shear index.  
c   I/O    ind(ng,nb)       int array  group index.
c    I     ist              integer    i index starting point.
c    I     jst              integer    j index starting point.
c    I     gst              integer    group index.
c    O     ier              integer    error message.
c
c Called Function: none
c***********************************************************************
       implicit none
       integer ng,nb,ist,jst,gst,ier
       integer ishr(8,ng,nb)
       integer ind(ng,nb)  !=0, no group data; =-1, cannot group data;
                             !=-2, missing data.
       integer nog,nob
       parameter( nog=940,nob=430 )
       integer indt(nog,nob)
       integer i,j,k,num,ng1,nb1
       integer nlimit
       parameter( nlimit=5 )

       ier=1
       ng1=ng-1
       nb1=nb-1
       do j=1,nb
       do i=1,ng
          indt(i,j)=ind(i,j)
       enddo
       enddo 

       num=1
       indt(ist,jst)=gst
       if( (indt(ist-1,jst+1).eq.0).and.(ishr(1,ist,jst).eq.1) )then
           indt(ist-1,jst+1)=gst
           num=num+1
       endif
       if( (indt(ist,jst+1).eq.0).and.(ishr(2,ist,jst).eq.1) )then
           indt(ist,jst+1)=gst
           num=num+1
       endif
       if( (indt(ist+1,jst+1).eq.0).and.(ishr(3,ist,jst).eq.1) )then
           indt(ist+1,jst+1)=gst
           num=num+1
       endif
       if( (indt(ist-1,jst).eq.0).and.(ishr(4,ist,jst).eq.1) )then
           indt(ist-1,jst)=gst
           num=num+1
       endif
       if( (indt(ist+1,jst).eq.0).and.(ishr(5,ist,jst).eq.1) )then
           indt(ist+1,jst)=gst
           num=num+1
       endif
       if( (indt(ist-1,jst-1).eq.0).and.(ishr(6,ist,jst).eq.1) )then
           indt(ist-1,jst-1)=gst
           num=num+1
       endif
       if( (indt(ist,jst-1).eq.0).and.(ishr(7,ist,jst).eq.1) )then
           indt(ist,jst-1)=gst
           num=num+1
       endif
       if( (indt(ist+1,jst-1).eq.0).and.(ishr(8,ist,jst).eq.1) )then
           indt(ist+1,jst-1)=gst
           num=num+1
       endif
       if( num.le.1 )return

       do 12 i=ist+2,ng,1
          if( (indt(i-1,jst).ne.gst).and.(indt(i-1,jst+1).ne.gst).and.
     1        (indt(i-1,jst-1).ne.gst) )go to 15
          if( indt(i,jst).eq.0 )then
              if( (indt(i-1,jst).eq.gst).and.
     1            (ishr(4,i,jst).eq.1) )then
                  indt(i,jst)=gst
                  num=num+1
                  go to 10
              endif
              if( (indt(i-1,jst+1).eq.gst).and.
     1            (ishr(1,i,jst).eq.1) )then
                  indt(i,jst)=gst
                  num=num+1
                  go to 10
              endif
              if( (indt(i-1,jst-1).eq.gst).and.
     1            (ishr(6,i,jst).eq.1) )then
                  indt(i,jst)=gst
                  num=num+1
                  go to 10
              endif
          endif
 10       continue
          if( indt(i,jst+1).eq.0 )then
              if( (indt(i,jst).eq.gst).and.
     1            (ishr(7,i,jst+1).eq.1) )then
                  indt(i,jst+1)=gst
                  num=num+1
                  go to 11 
              endif
              if( (indt(i-1,jst+1).eq.gst).and.
     1            (ishr(4,i,jst+1).eq.1) )then
                  indt(i,jst+1)=gst
                  num=num+1
                  go to 11 
              endif
              if( (indt(i-1,jst).eq.gst).and.
     1            (ishr(6,i,jst+1).eq.1) )then
                  indt(i,jst+1)=gst
                  num=num+1
                  go to 11
              endif
          endif
 11       continue
          if( indt(i,jst-1).eq.0 )then
              if( (indt(i,jst).eq.gst).and.
     1            (ishr(2,i,jst-1).eq.1) )then
                  indt(i,jst-1)=gst
                  num=num+1
                  go to 12
              endif 
              if( (indt(i-1,jst-1).eq.gst).and.
     1            (ishr(4,i,jst-1).eq.1) )then
                  indt(i,jst-1)=gst
                  num=num+1
                  go to 12 
              endif
              if( (indt(i-1,jst).eq.gst).and.
     1            (ishr(1,i,jst-1).eq.1) )then
                  indt(i,jst-1)=gst
                  num=num+1
                  go to 12
              endif
          endif
 12    continue 
 15    continue

       do 24 i=ist-2,1,-1
          if( (indt(i+1,jst).ne.gst).and.(indt(i+1,jst-1).ne.gst).and.
     1        (indt(i+1,jst+1).ne.gst) )go to 25
          if( indt(i,jst).eq.0 )then
              if( (indt(i+1,jst).eq.gst).and.
     1            (ishr(5,i,jst).eq.1) )then
                  indt(i,jst)=gst
                  num=num+1
                  go to 20
              endif
              if( (indt(i+1,jst+1).eq.gst).and.
     1            (ishr(3,i,jst).eq.1) )then
                  indt(i,jst)=gst
                  num=num+1
                  go to 20
              endif
              if( (indt(i+1,jst-1).eq.gst).and.
     1            (ishr(8,i,jst).eq.1) )then
                  indt(i,jst)=gst
                  num=num+1
                  go to 20
              endif
          endif
 20       continue
          if( indt(i,jst+1).eq.0 )then
              if( (indt(i,jst).eq.gst).and.
     1            (ishr(7,i,jst+1).eq.1) )then
                  indt(i,jst+1)=gst
                  num=num+1
                  go to 21
              endif 
              if( (indt(i+1,jst+1).eq.gst).and.
     1            (ishr(5,i,jst+1).eq.1) )then
                  indt(i,jst+1)=gst
                  num=num+1
                  go to 21
              endif
              if( (indt(i+1,jst).eq.gst).and.
     1            (ishr(8,i,jst+1).eq.1) )then
                  indt(i,jst+1)=gst
                  num=num+1
                  go to 21
              endif
          endif
 21       continue
          if( indt(i,jst-1).eq.0 )then
              if( (indt(i,jst).eq.gst).and.
     1            (ishr(2,i,jst-1).eq.1) )then
                  indt(i,jst-1)=gst
                  num=num+1
                  go to 24
              endif
              if( (indt(i+1,jst-1).eq.gst).and.
     1            (ishr(5,i,jst-1).eq.1) )then
                  indt(i,jst-1)=gst
                  num=num+1
                  go to 24
              endif
              if( (indt(i+1,jst).eq.gst).and.
     1            (ishr(3,i,jst-1).eq.1) )then
                  indt(i,jst-1)=gst
                  num=num+1
                  go to 24
              endif
          endif
 24    continue
 25    continue

       do 50 i=2,ng1,1
          do 34 j=jst+1,nb,1
             if( (indt(i-1,j-1).ne.gst).and.(indt(i,j-1).ne.gst).and.
     1           (indt(i+1,j-1).ne.gst) )go to 35
             if( indt(i,j).eq.0 )then
                 if( (indt(i,j-1).eq.gst).and.
     1               (ishr(7,i,j).eq.1) )then
                     indt(i,j)=gst
                     num=num+1
                     go to 33   
                 endif
                 if( (indt(i-1,j-1).eq.gst).and.
     1               (ishr(6,i,j).eq.1) )then
                     indt(i,j)=gst
                     num=num+1
                     go to 33   
                 endif
                 if( (indt(i+1,j-1).eq.gst).and.
     1               (ishr(8,i,j).eq.1) )then
                     indt(i,j)=gst
                     num=num+1
                     go to 33   
                 endif
             endif
 33          continue
             if( indt(i+1,j).eq.0 )then
                 if( (indt(i+1,j-1).eq.gst).and.
     1               (ishr(7,i+1,j).eq.1) )then
                     indt(i+1,j)=gst
                     num=num+1
                     go to 34   
                 endif
                 if( (indt(i,j).eq.gst).and.
     1               (ishr(4,i+1,j).eq.1) )then
                     indt(i+1,j)=gst
                     num=num+1
                     go to 34   
                 endif
                 if( (indt(i,j-1).eq.gst).and.
     1               (ishr(6,i+1,j).eq.1) )then
                     indt(i+1,j)=gst
                     num=num+1
                     go to 34   
                 endif
             endif
 34       continue
 35       continue
          do 37 j=jst-2,1,-1
             if( (indt(i-1,j+1).ne.gst).and.(indt(i,j+1).ne.gst).and.
     1           (indt(i+1,j+1).ne.gst) )go to 50
             if( indt(i,j).eq.0 )then
                 if( (indt(i,j+1).eq.gst).and.
     1               (ishr(2,i,j).eq.1) )then
                     indt(i,j)=gst
                     num=num+1
                     go to 36
                 endif
                 if( (indt(i-1,j+1).eq.gst).and.
     1               (ishr(1,i,j).eq.1) )then
                     indt(i,j)=gst
                     num=num+1
                     go to 36
                 endif
                 if( (indt(i+1,j+1).eq.gst).and.
     1               (ishr(3,i,j).eq.1) )then
                     indt(i,j)=gst
                     num=num+1
                     go to 36
                 endif
             endif
 36          continue
             if( indt(i+1,j).eq.0 )then
                 if( (indt(i+1,j+1).eq.gst).and.
     1               (ishr(2,i+1,j).eq.1) )then
                     indt(i+1,j)=gst
                     num=num+1
                     go to 37   
                 endif
                 if( (indt(i,j).eq.gst).and.
     1               (ishr(4,i+1,j).eq.1) )then
                     indt(i+1,j)=gst
                     num=num+1
                     go to 37   
                 endif
                 if( (indt(i,j+1).eq.gst).and.
     1               (ishr(1,i+1,j).eq.1) )then
                     indt(i+1,j)=gst
                     num=num+1
                     go to 37   
                 endif
             endif
 37       continue
 50    continue
       if( num.lt.nlimit )return
       ier=0

       do 80 j=2,nb1,1
       do 80 i=2,ng1,1
          if( indt(i,j).le.0 )go to 80
          if( (indt(i-1,j+1).eq.0).and.(ishr(1,i,j).eq.1) )then
              indt(i-1,j+1)=indt(i,j)
          endif
          if( (indt(i,j+1).eq.0).and.(ishr(2,i,j).eq.1) )then
              indt(i,j+1)=indt(i,j)
          endif
          if( (indt(i+1,j+1).eq.0).and.(ishr(3,i,j).eq.1) )then
              indt(i+1,j+1)=indt(i,j)
          endif
          if( (indt(i-1,j).eq.0).and.(ishr(4,i,j).eq.1) )then
              indt(i-1,j)=indt(i,j)
          endif
          if( (indt(i+1,j).eq.0).and.(ishr(5,i,j).eq.1) )then
              indt(i+1,j)=indt(i,j)
          endif
          if( (indt(i-1,j-1).eq.0).and.(ishr(6,i,j).eq.1) )then
              indt(i-1,j-1)=indt(i,j)
          endif
          if( (indt(i,j-1).eq.0).and.(ishr(7,i,j).eq.1) )then
              indt(i,j-1)=indt(i,j)
          endif
          if( (indt(i+1,j-1).eq.0).and.(ishr(8,i,j).eq.1) )then
              indt(i+1,j-1)=indt(i,j)
          endif
 80    continue
       do j=1,nb
       do i=1,ng
          ind(i,j)=indt(i,j)
       enddo
       enddo
      
       return
       end

       subroutine del_isolated(nogate,nobeam,jend,dat,velny,bmiss,ivad)
c***********************************************************************
c description      : To delete the isolated Vr data.
c                                                                       
c I/O parameters   :                                                    
c  I/O/W   name,            type,      description
c    I     nogate           integer    gate number in one ray.
c    I     nobeam           integer    total ray number of PPI sweep scan.
c    I     jend             integer    the last index along azimthual direction.
c   I/O  dat(nogate,nobeam) real array radial velocity in unit of m/s.
c    I     velny            real       Nquist velocity in unit of m/s.
c    I     bmiss            real       missing value for radial velocity.
c    I     ivad             integer    ppi scanning message.
c                                      =0, is 360 degrees scan;  =1, isnot
c
c Called Function: none
c***********************************************************************
       implicit none
       integer nogate,nobeam,jend,ivad
       real dat(nogate,nobeam),velny,bmiss
       integer i,j,j1,jend1,jend2,miss,shear,n,ij,ivalue
       integer nog,nob
       parameter( nog=940,nob=430 )
       integer idat(nog,nob)

       shear=50*velny
       miss=-32767
       if( ivad.ne.0 )go to 100

       jend1=jend+1
       jend2=jend+2
     
       do j=2,jend1,1
       do i=1,nogate
          j1=j-1
          if( abs(dat(i,j1)).lt.99. )then
              idat(i,j)=100*dat(i,j1)
          else
              idat(i,j)=miss
          endif
       enddo
       enddo
       do i=1,nogate
          idat(i,1)=idat(i,jend1)
          idat(i,jend2)=idat(i,2)
       enddo

       n=0
       do j=2,jend1,1
       do i=2,nogate-1,1
          if( (idat(i,j).ne.miss).and.(abs(idat(i,j)).gt.shear) )then
              ij=0
              ivalue=0
              if( (idat(i+1,j).ne.miss).and.
     1            (abs(idat(i+1,j)).lt.shear) )then
                  ij=ij+1
                  ivalue=ivalue+idat(i+1,j)
              endif
              if( (idat(i-1,j).ne.miss).and.
     1            (abs(idat(i-1,j)).lt.shear) )then
                  ij=ij+1
                  ivalue=ivalue+idat(i-1,j)
              endif
              if( (idat(i,j+1).ne.miss).and.
     1            (abs(idat(i,j+1)).lt.shear) )then
                  ij=ij+1
                  ivalue=ivalue+idat(i,j+1)
              endif
              if( (idat(i,j-1).ne.miss).and.
     1            (abs(idat(i,j-1)).lt.shear) )then
                  ij=ij+1
                  ivalue=ivalue+idat(i,j-1)
              endif
              if( ij.gt.2 )then
                  ivalue=ivalue/ij
                  if( abs(idat(i,j)-ivalue).gt.shear )then
                      n=n+1
                      idat(i,j)=miss
                  endif
              endif 
          endif
       enddo
       enddo
       print*,'number of deleting isolated data: ',n

       do j=2,jend1,1
       do i=2,nogate-1,1
          j1=j-1
          if( idat(i,j).ne.miss )then
              dat(i,j1)=0.01*idat(i,j)
          endif
       enddo
       enddo
       return

 100   continue
       do j=1,nobeam
       do i=1,nogate
          if( abs(dat(i,j)).lt.99. )then
              idat(i,j)=100*dat(i,j)
          else
              idat(i,j)=miss
          endif
       enddo
       enddo
       n=0
       do j=2,nobeam-1,1
       do i=2,nogate-1,1
          if( (idat(i,j).ne.miss).and.(abs(idat(i,j)).gt.shear) )then
              ij=0
              ivalue=0
              if( (idat(i+1,j).ne.miss).and.
     1            (abs(idat(i+1,j)).lt.shear) )then
                  ij=ij+1
                  ivalue=ivalue+idat(i+1,j)
              endif
              if( (idat(i-1,j).ne.miss).and.
     1            (abs(idat(i-1,j)).lt.shear) )then
                  ij=ij+1
                  ivalue=ivalue+idat(i-1,j)
              endif
              if( (idat(i,j+1).ne.miss).and.
     1            (abs(idat(i,j+1)).lt.shear) )then
                  ij=ij+1
                  ivalue=ivalue+idat(i,j+1)
              endif
              if( (idat(i,j-1).ne.miss).and.
     1            (abs(idat(i,j-1)).lt.shear) )then
                  ij=ij+1
                  ivalue=ivalue+idat(i,j-1)
              endif
              if( ij.gt.2 )then
                  ivalue=ivalue/ij
                  if( abs(idat(i,j)-ivalue).gt.shear )then
                      n=n+1
                      idat(i,j)=miss
                  endif
              endif 
          endif
       enddo
       enddo
       print*,'number of deleting isolated data: ',n

       do j=2,nobeam-1,1
       do i=2,nogate-1,1
          if( idat(i,j).ne.miss )then
              dat(i,j)=0.01*idat(i,j)
          endif
       enddo
       enddo

       return
       end
 
       subroutine simple_unf(nogate,nobeam,jend,dat,velny,bmiss,ivad)
c***********************************************************************
c description      : Simple unfold by 3 x 3 cell.
c                                                                       
c I/O parameters   :                                                    
c  I/O/W   name,            type,      description
c    I     nogate           integer    gate number in one ray.
c    I     nobeam           integer    total ray number of PPI sweep scan.
c    I     jend             integer    the last index along azimthual direction.
c   I/O  dat(nogate,nobeam) real array radial velocity in unit of m/s.
c    I     velny            real       Nquist velocity in unit of m/s.
c    I     bmiss            real       missing value for radial velocity.
c    I     ivad             integer    ppi scanning message.
c                                      =0, is 360 degrees scan;  =1, isnot
c
c Called Function: none
c***********************************************************************
       implicit none
       integer nogate,nobeam,jend,ivad
       real dat(nogate,nobeam),velny,bmiss

       integer nog,nob
       parameter( nog=940,nob=430 )
       integer idat(nog,nob)
       integer miss,nx,ny,ivelny2,ivelny5,nvelny5
       integer i,j,i1,i2,j1,j2,ii,jj
       integer np,nn,nzp,nzn,ivalue,n

       miss=-32767
       ivelny2=200*velny
       ivelny5=50*velny
       nvelny5=-ivelny5
       nx=nogate/3
       if( ivad.eq.0 )then
           ny=jend/3+1
       else
           ny=nobeam/3
       endif

       do j=1,jend
       do i=1,nogate
          if( abs(dat(i,j)).lt.99.0 )then
              idat(i,j)=100*dat(i,j)
          else
              idat(i,j)=miss
          endif
       enddo
       enddo
       do i=1,nogate
          idat(i,jend+1)=idat(i,1)
          idat(i,jend+2)=idat(i,2)
          idat(i,jend+3)=idat(i,3)
       enddo

       n=0
       do j=1,ny
          j1=3*(j-1)+1
          j2=j1+2
          do 50 i=2,nx
             i1=3*(i-1)+1
             i2=i1+2 
             np=0
             nn=0
             nzp=0
             nzn=0
             do jj=j1,j2
             do ii=i1,i2
                if( idat(ii,jj).ne.miss )then
                    if( idat(ii,jj).ge.0 )then
                        nzp=nzp+1
                    else
                        nzn=nzn+1
                    endif
                    if( abs(idat(ii,jj)).gt.ivelny5 )then
                        if( idat(ii,jj).gt.0 )then
                            np=np+1
                        else
                            nn=nn+1
                        endif
                    endif
                endif
             enddo 
             enddo
             if( (np.eq.0).or.(nn.eq.0) )go to 50
             if( nzp.gt.nzn )then
                 do jj=j1,j2
                 do ii=i1,i2
                    if( (idat(ii,jj).ne.miss).and.
     1                  (idat(ii,jj).lt.nvelny5) )then
                        idat(ii,jj)=idat(ii,jj)+ivelny2
                        n=n+1
                    endif
                 enddo
                 enddo
             else
                 do jj=j1,j2
                 do ii=i1,i2
                    if( (idat(ii,jj).ne.miss).and.
     1                  (idat(ii,jj).gt.ivelny5) )then
                        idat(ii,jj)=idat(ii,jj)-ivelny2
                        n=n+1
                    endif
                 enddo
                 enddo
             endif
 50       continue
       enddo
       do j=1,jend
       do i=1,nogate
          if( idat(i,j).ne.miss )then
              dat(i,j)=0.01*idat(i,j)
          endif
       enddo
       enddo
       print*,'number of simple unfold: ',n

       return
       end

       subroutine get_azimuindex(nobeam,beam,jend,ivad)
c***********************************************************************
c description      : To get the last index along azimthual direction.
c                                                                       
c I/O parameters   :                                                    
c  I/O/W   name,            type,      description
c    I     nobeam           integer    total ray number of PPI sweep scan.
c    I     beam(420)        real array azimthual angle. (degrees)
c    O     jend             integer    the last index along azimthual direction.
c    O     ivad             integer    ppi scanning message.
c                                      =0, is 360 degrees scan;  =1, isnot
c Called Function: none
c***********************************************************************
       implicit none
       integer nobeam,jend,ivad
       real beam(nobeam)
       real ang_tor,ang_tor1,angmax
       parameter( ang_tor=5.0,ang_tor1=360.-ang_tor )
       parameter( angmax=20. )
       real ang,ang1,ang2
       integer i,nob2

       ivad=0
       jend=nobeam
       nob2=nobeam/2
       do i=2,nobeam
          ang=abs(beam(i)-beam(i-1))
          if( ang.gt.300. )ang=360.-ang
          if( ang.gt.ang_tor )then
              ivad=1
              return
          endif
       enddo
       ang=abs(beam(1)-beam(nobeam))
       if( ang.gt.300. )ang=360.-ang
       if( ang.gt.angmax )then
           ivad=1
           return
       endif 

       ang=beam(2)-beam(1)
       if( ang.gt.0. )then
           if( abs(ang).gt.300. )then
               go to 100
           else
               go to 50
           endif
       else
           if( abs(ang).gt.300. )then
               go to 50
           else
               go to 100 
           endif
       endif

 50    continue
       ang=beam(1)-beam(nobeam)
       if( ang.ge.0. )then
           if( ang.le.ang_tor )then
               jend=nobeam
           else
               if( abs(ang).gt.300. )go to 60
               ivad=1
           endif
           return
       else
           if( abs(ang).ge.ang_tor1 )then
               jend=nobeam
               return
           endif
       endif

 60    continue
       do i=nobeam-1,nob2,-1
          ang=beam(i+1)-beam(i)
          ang1=beam(i+1)-beam(1)
          ang2=beam(i)-beam(1)
          if( abs(ang).gt.300. )then
              if( abs(ang1).gt.300. )then
                  ang1=360.-ang1
              endif
              if( abs(ang2).gt.300. )then
                  ang2=ang2-360.
              endif
          endif
          if( ang1*ang2.le.0. )then
              jend=i
              return
          endif
       enddo
       ivad=1
       return

 100   continue
       ang=beam(nobeam)-beam(1)
       if( ang.ge.0. )then
           if( ang.le.ang_tor )then
               jend=nobeam
           else
               if( abs(ang).gt.300. )go to 110
               ivad=1
           endif
           return
       else
           if( abs(ang).ge.ang_tor1 )then
               jend=nobeam
               return
           endif
       endif
 110   continue
       do i=nobeam-1,nob2,-1
          ang=beam(i+1)-beam(i)
          ang1=beam(i+1)-beam(1)
          ang2=beam(i)-beam(1)
          if( abs(ang).gt.300. )then
              if( abs(ang1).gt.300. )then
                  ang1=ang1-360.
              endif
              if( abs(ang2).gt.300. )then
                  ang2=360.-ang2
              endif
          endif
          if( ang1*ang2.le.0. )then
              jend=i
              return
          endif
       enddo
       ivad=1

       return
       end

       subroutine get_ishr(ng,nb,miss,shear,idat,ishr)
c***********************************************************************
c description      : To get shear index.
c                                                                       
c I/O parameters   :                                                    
c  I/O/W   name,            type,      description
c    I     ng               integer    first dimension of input data array.
c    I     nb               integer    second dimension of input data array.
c    I     miss             integer    missing value.
c    I     shear            integer    shear velocity.
c    I     idat(ng,nb)      int array  input data array
c    O     ishr(8,ng,nb)    int array  shear index. =1 for small shear.
c
c Called Function: none
c***********************************************************************
       implicit none
       integer ng,nb,miss,shear
       integer idat(ng,nb),ishr(8,ng,nb)
       integer i,j,k
       integer ng1,nb1

       do j=1,nb
       do i=1,ng 
       do k=1,8
          ishr(k,i,j)=0
       enddo
       enddo
       enddo
       ng1=ng-1
       nb1=nb-1

       do j=2,nb1,1
       do i=2,ng1,1
          if( idat(i,j).ne.miss )then
              if( (idat(i-1,j+1).ne.miss).and.
     1            (abs(idat(i-1,j+1)-idat(i,j)).lt.shear) )then
                  ishr(1,i,j)=1
                  ishr(8,i-1,j+1)=1
              endif
              if( (idat(i,j+1).ne.miss).and.
     1            (abs(idat(i,j+1)-idat(i,j)).lt.shear) )then
                  ishr(2,i,j)=1
                  ishr(7,i,j+1)=1
              endif
              if( (idat(i+1,j+1).ne.miss).and.
     1            (abs(idat(i+1,j+1)-idat(i,j)).lt.shear) )then
                  ishr(3,i,j)=1
                  ishr(6,i+1,j+1)=1
              endif
              if( (idat(i-1,j).ne.miss).and.
     1            (abs(idat(i-1,j)-idat(i,j)).lt.shear) )then
                  ishr(4,i,j)=1
                  ishr(5,i-1,j)=1
              endif
          endif
       enddo
       enddo

       do j=1,nb1
          if( (idat(ng,j).ne.miss).and.(idat(ng1,j+1).ne.miss).and. 
     1        (abs(idat(ng,j)-idat(ng1,j+1)).lt.shear) )then
              ishr(1,ng,j)=1
              ishr(8,ng1,j+1)=1
          endif
       enddo
       do i=2,ng1
          if( (idat(i,1).ne.miss).and.(idat(i-1,2).ne.miss).and.
     1        (abs(idat(i,1)-idat(i-1,2)).lt.shear) )then
              ishr(1,i,1)=1
              ishr(8,i-1,2)=1
          endif
       enddo

       do j=1,nb1
          if( (idat(1,j).ne.miss).and.(idat(1,j+1).ne.miss).and.
     1        (abs(idat(1,j)-idat(1,j+1)).lt.shear) )then
              ishr(2,1,j)=1
              ishr(7,1,j+1)=1
          endif
          if( (idat(ng,j).ne.miss).and.(idat(ng,j+1).ne.miss).and.
     1        (abs(idat(ng,j)-idat(ng,j+1)).lt.shear) )then
              ishr(2,ng,j)=1
              ishr(7,ng,j+1)=1
          endif
       enddo
       do i=2,ng1
          if( (idat(i,1).ne.miss).and.(idat(i,2).ne.miss).and.
     1        (abs(idat(i,1)-idat(i,2)).lt.shear) )then
              ishr(2,i,1)=1
              ishr(7,i,2)=1
          endif
       enddo

       do j=1,nb1
          if( (idat(1,j).ne.miss).and.(idat(2,j+1).ne.miss).and.
     1        (abs(idat(1,j)-idat(2,j+1)).lt.shear) )then
              ishr(3,1,j)=1
              ishr(6,2,j+1)=1
          endif
       enddo
       do i=2,ng1
          if( (idat(i,1).ne.miss).and.(idat(i+1,2).ne.miss).and.
     1        (abs(idat(i,1)-idat(i+1,2)).lt.shear) )then
              ishr(3,i,1)=1
              ishr(6,i+1,2)=1
          endif
       enddo

       do j=1,nb
          if( (idat(ng,j).ne.miss).and.(idat(ng1,j).ne.miss).and.
     1        (abs(idat(ng,j)-idat(ng1,j)).lt.shear) )then
              ishr(4,ng,j)=1
              ishr(5,ng1,j)=1
          endif
       enddo
       do i=2,ng1
          if( (idat(i,1).ne.miss).and.(idat(i-1,1).ne.miss).and.
     1        (abs(idat(i,1)-idat(i-1,1)).lt.shear) )then
              ishr(4,i,1)=1
              ishr(5,i-1,1)=1
          endif
          if( (idat(i,nb).ne.miss).and.(idat(i-1,nb).ne.miss).and.
     1        (abs(idat(i,nb)-idat(i-1,nb)).lt.shear) )then
              ishr(4,i,nb)=1
              ishr(5,i-1,nb)=1
          endif 
       enddo

       return
       end
