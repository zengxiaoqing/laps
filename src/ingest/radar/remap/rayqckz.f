       subroutine rayqckz(ng,nb,iscale,miss,idat,azimu)
c***********************************************************************
c Subroutine/Function : rayqckz
c
c Usage :
c    call rayqckz(ng,nb,iscale,miss,idat,azimu)
c
c Description      : To QC for ray contining data.
c
c Arguments :
c  I/O/W   name,      type,       description
c    I     ng         integer     the number of gate. (the first dimension) 
c    I     nb         integer     the number of azimuthal angle. (the second dimension)
c    I     iscale     integer     scale value.
c    I     miss       integer     the missing data.
c   I/O   idat(ng,nb) int array   the reflectivity data.
c    I     azimu(nb)  real array  the azimuthal angles.
c
c Called Function : none
c***********************************************************************

       dimension idat(ng,nb),azimu(nb)
       dimension ick(920),jck(450),nck(450)
       dimension jsta(450),nsta(450),jend(450),nend(450)

       nnn=0
       num=0
       ir1=5.*iscale+0.5
       ir2=10.*iscale+0.5
       nmax=50
       izbd1=-10.*iscale
       izbd2=25.*iscale

       do i=1,ng
          ick(i)=ir1+(ir2-ir1)*(i-1)/float(ng-1)+0.5
       enddo
       do i=1,450
          jck(i)=0
          nck(i)=0
       enddo
     
       jb1=nb
       azmin=abs(azimu(nb)-azimu(1))
       if( azmin.gt.330. )azmin=abs(azmin-360.)
       do j=nb-1,nb-10,-1
          az=abs(azimu(j)-azimu(1))
          if( az.gt.330.)az=abs(az-360.)
          if( az.lt.azmin )then
              jb1=j
              azmin=az
          endif
       enddo
       jb2=1
       azmin=abs(azimu(1)-azimu(nb))
       if( azmin.gt.330. )azmin=abs(azmin-360.)
       do j=2,10,1
          az=abs(azimu(j)-azimu(nb))
          if( az.gt.330.)az=abs(az-360.)
          if( az.lt.azmin )then
              jb2=j
              azmin=az
          endif
       enddo

       nn=0
       do 30 j=1,nb

          ii=0
          do 10 i=2,ng
             if( (idat(i,j).eq.miss).or.(idat(i-1,j).eq.miss) )then
                 ii=0
                 go to 10
             endif
             if( (idat(i,j).le.izbd1).or.(idat(i,j).gt.izbd2) )then
                 ii=0
                 go to 10
             endif
             if( abs(idat(i,j)-idat(i-1,j)).lt.ir1 )then
                 ii=ii+1
             else
                 ii=0
             endif
             if( ii.gt.nmax )then
                 go to 15
             endif
 10       continue
          if( ii.le.nmax )go to 30
 15       continue

          jm1=j-1
          jp1=j+1
          if( j.eq.1 )jm1=jb1
          if( j.eq.nb )jp1=jb2

          n1=0
          n2=0
          do 20 i=1,ng
             
             if( idat(i,j).eq.miss )then
                 n1=0
                 n2=0
                 go to 20
             endif

             if( idat(i,jm1).ne.miss )then
                 if( abs(idat(i,j)-idat(i,jm1)).ge.ick(i) )then
                     n1=n1+1
                     if( n1.eq.1 )ndx1=i
                 else
                     n1=0
                 endif
             else
               if( (idat(i,j).gt.izbd1).and.(idat(i,j).le.izbd2) )then
                     n1=n1+1
                     if( n1.eq.1 )ndx1=i
               else
                     n1=0
               endif
             endif
             if( n1.gt.nmax )then
                 nn=nn+1
                 jck(nn)=j
                 jsta(nn)=ndx1
                 jend(nn)=i
                 do ii=i+1,ng,1
                    if( idat(ii,jm1).ne.miss )then
                      if( abs(idat(ii,j)-idat(ii,jm1)).lt.ick(ii) )then
                            jend(nn)=ii
                            go to 30
                        endif
                    else
                        if( idat(ii,j).ge.izbd2 )then
                            jend(nn)=ii-1
                            go to 30
                        endif
                        if( idat(ii,j).eq.miss )then
                            jend(nn)=ii-1
                            go to 30
                        endif
                    endif
                 enddo
                 jend(nn)=ng
                 go to 30
             endif

             if( idat(i,jp1).ne.miss )then
                 if( abs(idat(i,j)-idat(i,jp1)).ge.ick(i) )then
                     n2=n2+1
                     if( n2.eq.1 )ndx2=i
                 else
                     n2=0
                 endif
             else
               if( (idat(i,j).gt.izbd1).and.(idat(i,j).le.izbd2) )then
                     n2=n2+1
                     if( n2.eq.1 )ndx2=i
               else
                     n2=0
               endif
             endif
             if( n2.gt.nmax )then
                 nn=nn+1
                 jck(nn)=j
                 jsta(nn)=ndx2
                 jend(nn)=i
                 do ii=i+1,ng,1
                    if( idat(ii,jp1).ne.miss )then
                      if( abs(idat(ii,j)-idat(ii,jp1)).lt.ick(ii) )then
                            jend(nn)=ii
                            go to 30
                        endif
                    else
                        if( idat(ii,j).gt.izbd2 )then
                            jend(nn)=ii-1
                            go to 30
                        endif
                        if( idat(ii,j).eq.miss )then
                            jend(nn)=ii-1
                            go to 30
                        endif
                    endif
                 enddo
                 jend(nn)=ng
                 go to 30
             endif
 20       continue
          
 30    continue

       if( nn.eq.0 )return
       nnn=nn
       do k=1,nn
          nck(k)=jck(k)
          nsta(k)=jsta(k)
          nend(k)=jend(k)
       enddo

       if( nn.eq.1 )then
           if( (nck(1).lt.3).or.(nck(1).gt.(nb-2)) )go to 34
           jj=nck(1)
           jm2=nck(1)-2
           jm1=nck(1)-1
           jp1=nck(1)+1
           jp2=nck(1)+2
           n1=0
           n2=0
           istart=nsta(1)
           iend=nend(1)
           do i=istart,iend
              if( (idat(i,jj).ne.miss).and.(idat(i,jm1).ne.miss)
     1            .and.(idat(i,jm2).eq.miss) )then
                   n1=n1+1
                   if( n1.eq.1 )ndx1=i
              endif
              if( (idat(i,jj).ne.miss).and.(idat(i,jp1).ne.miss)
     1            .and.(idat(i,jp2).eq.miss) )then
                   n2=n2+1
                   if( n2.eq.1 )ndx2=i
              endif
           enddo
           if( n1.gt.nmax )then
               nnn=nnn+1
               nck(nnn)=jm1
               nsta(nnn)=ndx1
               nend(nnn)=iend
           endif
           if( n2.gt.nmax )then
               nnn=nnn+1
               nck(nnn)=jp1
               nsta(nnn)=ndx2
               nend(nnn)=iend
           endif
           go to 34
       endif

       do  33 k=1,nn
           if( (nck(k).lt.3).or.(nck(k).gt.(nb-2)) )go to 33
           jj=nck(k)
           jm1=nck(k)-1
           jp1=nck(k)+1
           do kk=1,nn
              if( k.ne.kk )then
                  if( nck(kk).eq.jm1 )go to 33
                  if( nck(kk).eq.jp1 )go to 33
              endif
           enddo
           jm2=nck(k)-2
           jp2=nck(k)+2
           n1=0
           n2=0
           istart=nsta(k)
           iend=nend(k)
           do i=istart,iend
              if( (idat(i,jj).ne.miss).and.(idat(i,jm1).ne.miss)
     1            .and.(idat(i,jm2).eq.miss) )then
                   n1=n1+1
                   if( n1.eq.1 )ndx1=i
              endif
              if( (idat(i,jj).ne.miss).and.(idat(i,jp1).ne.miss)
     1            .and.(idat(i,jp2).eq.miss) )then
                   n2=n2+1
                   if( n2.eq.1 )ndx2=i
              endif
           enddo
           if( n1.gt.nmax )then
               nnn=nnn+1
               nck(nnn)=jm1
               nsta(nnn)=ndx1
               nend(nnn)=iend
           endif
           if( n2.gt.nmax )then
               nnn=nnn+1
               nck(nnn)=jp1
               nsta(nnn)=ndx2
               nend(nnn)=iend
           endif
 33    continue

 34    continue

       if( (nnn.gt.0).and.(nnn.le.50) )then
           do k=1,nnn
              j=nck(k)
              i1=nsta(k)
              i2=nend(k)
              if( (i1.ge.1).and.(i2.ge.i1).and.(i2.le.ng) )then
                  do i=i1,i2,1
                     idat(i,j)=miss
                  enddo
              endif
           enddo
       endif

       return
       end
