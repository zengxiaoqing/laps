c
c
      subroutine mvmult(a,b,c,imax,jmax,nmax,m)
c
c*********************************************************************
c     Routine to compute the inner product of A B = C
c
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       21 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c
c*********************************************************************
c
      real a(m,m),b(m,m),c(m,m),sum
      integer imax,jmax,nmax,j,n,i,m
c
      do i=1,imax
      do n=1,nmax
         sum = 0.
         do j=1,jmax
            sum = sum + (a(i,j) * b(j,n))
         enddo !j
         c(i,n) = sum
      enddo !n 
      enddo !i
c
      return
      end
c
c
      subroutine addmv(a,b,c,imax,jmax,m)
c
c*********************************************************************
c
c     Routine adds two vectors or matricies.  Results are output in c.    
c
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       21 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c       13 Oct 1998  Peter Stamus, NOAA/FSL
c          Removed subract option...put in separate routine.
c
c*********************************************************************
c
      real a(m,m),b(m,m),c(m,m)
      integer imax,jmax,i,j,m
c
      do j=1,jmax
      do i=1,imax
         c(i,j) = a(i,j) + b(i,j)
      enddo !i
      enddo !j
c
      return
      end
c
c
      subroutine submv(a,b,c,imax,jmax,m)
c
c*********************************************************************
c
c     Routine subrtacts two vectors or matricies.  Results are output 
c     in c.    
c
c     Original: Peter Stamus, NOAA/FSL  13 Oct 1998 (from 'addmv')
c     Changes:
c
c*********************************************************************
c
      real a(m,m),b(m,m),c(m,m)
      integer imax,jmax,i,j,m
c
      do j=1,jmax
      do i=1,imax
         c(i,j) = a(i,j) - b(i,j)
      enddo !i
      enddo !j
c
      return
      end
c
c
      subroutine trans(a,at,imax,jmax,m)
c
c*********************************************************************
c     Routine creates the transposed matrix of a and places it in 
c     the at array.
c
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       21 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c
c*********************************************************************
c
      real a(m,m),at(m,m)
      integer imax,jmax,i,j,m
c
      do j=1,jmax
      do i=1,imax
         at(i,j) = a(j,i)
      enddo !i
      enddo !j
c
      return
      end
c
c
      subroutine vomult(a,b,c,imax,m)
c
c*********************************************************************
c     Routine provides the outer product of two vectors a bt = c 
c     matrix.
c
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       21 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c
c*********************************************************************
c
      real a(m),b(m),c(m,m)
      integer imax,j,i,m
c
      do j=1,imax
      do i=1,imax
         c(i,j) = a(i) * b(j)
      enddo !i
      enddo !j
c
      return
      end
c
c
      subroutine invert (a,n,np,y,maxsta)  
c
c*********************************************************************
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       21 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c       14 Dec 1999  Peter Stamus, NOAA/FSL
c          Pass in max number of stations for 'ludcmp'. 
c
c*********************************************************************
c
      real a(np,np),y(np,np),d
      integer np,indx(np),n,j,i
c
      do i=1,n
         do j=1,n
            y(i,j)=0.
         enddo !j
         y(i,i)=1.
      enddo !i
      call ludcmp(a,n,np,indx,d,maxsta)  
      do j=1,n
         call lubksb(a,n,np,indx,y(1,j))  
      enddo !j
c
      return
      end
c
      subroutine fastinv(A,imax,jmax,m)
      real A(m,m)      
c     This inverts a matrix that is diagonal only
      if(imax.ne.jmax) print*,'matrix in fast invert is not square'
      print*,'Matrix to be inverted is diagonal...using fast invert'
      do j=1,jmax
        if(A(j,j).eq.0) then
         print*,' zero element in fast matrix inversion at ',j,j
         A(j,j)=1.e-30
        endif
        A(j,j)=1/A(j,j)
      enddo
      return
      end
 
c
      subroutine avgerr(wr,wit,B,c,W,dwt,imax,iav,m,it)
c
c*********************************************************************
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       21 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c       14 Dec 1999  John McGinley and Peter Stamus, NOAA/FSL
c          Add Bias error calculation.  Also some minor housekeeping.
c
c*********************************************************************
c
      real B(m),c(m,m),W(m,m),wit(m,m),wr(m),dwt(m,m)
c
      do l=m,2,-1
         do k=1,imax
            wit(l,k)=wit(l-1,k)
         enddo !k
      enddo !l
      do k=1,m   
         wit(1,k)=wr(k)
      enddo !k
      do i=1,m
         do j=1,m
            W(i,j)=0. 
         enddo !j
      enddo !i
      if(iav.gt.it) then
         ia=it
      else
         ia=iav
      endif
c
      sum = 0.
      cnt = 0.
      do itt=1,ia
         do i=1,m
            B(i)=wit(itt,i)   
            sum = sum + B(i)
            cnt = cnt + 1.
         enddo !i
         do j=1,imax
            do i=1,imax
               C(i,j)=B(i)*B(j)*dwt(i,j)
            enddo !i
         enddo !j
        do j=1,imax
           do i=1,imax
              W(i,j)=W(i,j)+c(i,j)
           enddo !i
        enddo !j
      enddo !itt
      do j=1,imax
         do i=1,imax
            if(i.eq.j) then
               W(i,j)=W(i,j)/float(ia)
            else
               W(i,j)=W(i,j)/float(iav)
            endif
         enddo !i
      enddo !j
c
      if(cnt .gt. 0.) then
         write(6,*) ' Bias Error for W and VV: ', sum / cnt
      else
         write(6,*) ' In routine AVGERR: Cannot calculate bias ',
     &              ' error for W and VV -- cnt = 0'
      endif
c
      return
      end
c
c
      subroutine avgdiagerr(wit,B,c,W,imax,iav,m,it)
c
c*********************************************************************
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       21 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c       14 Dec 1999  John McGinley and Peter Stamus, NOAA/FSL
c          Add Bias error calculation.  Also some minor housekeeping.
c
c*********************************************************************
c
      real B(m),c(m,m),W(m,m),wit(m,m)
c
      call zero(W,m,m)
      ia=iav
      if(ia.gt.it) ia=it
c
      sum = 0.
      cnt = 0.
      do itt=1,ia
         do i=1,imax
            B(i)=wit(i,itt)   
            sum = sum + B(i)
            cnt = cnt + 1.
         enddo !i
         do i=1,imax    
            W(i,i)=W(i,i)+b(i)*b(i)
         enddo !i
      enddo !itt
      do i=1,imax
         W(i,i)=W(i,i)/float(ia)
      enddo !i
c
      if(cnt .gt. 0.) then
         write(6,*) ' Bias Error for W and VV: ', sum / cnt,' for ', ia
     & ,' cycles'
      else
         write(6,*) ' In routine AVGDIAGERR: Cannot calculate ',
     &              ' bias error for W and VV -- cnt = 0'
      endif
c
      return
      end
c
c
      subroutine matrixanal(a,imax,jmax,m,idiag,char)
c
c*********************************************************************
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       21 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c*********************************************************************
c
      real a(m,m)     
      character char*10
c
      dimx=-1.e20
      dimn=1.e20
      elmx=-1.e20
      elmn=1.e20
c
c  diagonal element
c
      sumd=0.
      cnt=0.
      do j=1,jmax
         if(a(j,j).gt.dimx)dimx=a(j,j)
         if(a(j,j).lt.dimn)dimn=a(j,j)
         sumd=sumd+a(j,j)
         cnt=cnt+1.
      enddo !j
      sumd=sumd/cnt
c
c  remander of elements
c
      sumea=0.
      sume=0.
      cnt=0.
      do j=1,jmax
         do i=1,imax
            if (i.eq.j) go to 1
            if(a(i,j).gt.elmx) elmx=a(i,j)
            if(a(i,j).lt.elmn) elmn=a(i,j)
            sume=sume+a(i,j)
            sumea=sumea+abs(a(i,j))
            cnt=cnt+1.
 1       enddo !i
      enddo !j
      sume=sume/cnt
      sumea=sumea/cnt
c
c output results
c
      write(6,1000)  char
 1000 format(1x,'Diagnosis of matrix ',a10)
      write (6,1001) dimx,dimn
 1001 format(1x,'Max diagonal ',f8.3,' Min diagonal ',f8.3)
      write (6,1002) sumd
 1002 format(1x,'Avg diagonal ',f8.3)
      write (6,1003) elmx,elmn
 1003 format(1x,'Max non-diagonal ',f8.3,' Min non-diagonal ',f8.3)
      write (6,1004) sume
 1004 format(1x,'Avg non-diagonal ',f8.3)
      write (6,1005) sumea
 1005 format(1x,'Avg abs non-diagonal ',f8.3)
c
c flag matrix if it is purely diagonal
      if(elmx.eq.0..and.elmn.eq.0.) idiag=1
      return
      end
c
c
      Subroutine svdcmp(a,m,n,mp,np,w,v,nmax)
c
c*********************************************************************
c
c     Given a matrix a(m,n) with physical dimensions 'mp x np',
c     routine calculates a=u.w.vt. u repaces a on output
c     diagonal matrix of singular values w is output as a vector
c     v is output, not vtranspose.
c
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       21 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c       14 Dec 1999  Peter Stamus, NOAA/FSL
c          Pass in max number of stations as nmax.
c
c*********************************************************************
c
       integer m,mp,n,np,nmax
       integer i,its,j,jj,k,l,nm
       real a(mp,np),v(np,np),w(np)
       real anorm,c,f,g,h,s,scale,x,y,z,rv1(nmax),pythag
c
       g=0.
       scale=0.
       anorm=0.
       do i=1,n
          l=i+1
          rv1(i)=scale*g
          g=0.
          s=0.
          scale=0.
          if(i.le.m) then
             do k=i,m
                scale=scale+abs(a(k,i))
             enddo !on k
             if(scale.ne.0.) then
                do k=i,m
                   a(k,i)=a(k,i)/scale
                   s=s+a(k,i)*a(k,i)
                enddo !on k
                f=a(i,i)
                g=-sign(sqrt(s),f)
                h=f*g-s
                a(i,i)=f-g
                do j=l,n
                   s=0.
                   do k=i,m
                      s=s+a(k,i)*a(k,j)
                   enddo !on k
                   f=s/h
                   do k=i,m
                      a(k,j)=a(k,j)+f*a(k,i)
                   enddo !on k
                enddo !on j
                do k=i,m
                   a(k,i)=scale*a(k,i)
                enddo !on k
             endif
          endif
          w(i)=scale*g
          g=0.
          s=0.
          scale=0.
          if((i.le.m).and.(i.ne.n))then
             do k=l,n
                scale=scale+abs(a(i,k))
             enddo !on k
             if(scale.ne.0.)then
                do k=l,n
                   a(i,k)=a(i,k)/scale
                   s=s+a(i,k)*a(i,k)
                enddo !k
                f=a(i,l)
                g=-sign(sqrt(s),f)
                h=f*g-s
                a(i,l)=f-g
                do k=l,n
                   rv1(k)=a(i,k)/h    
                enddo !on k
                do j=l,m
                   s=0.
                   do k=l,n
                      s=s+a(j,k)*a(i,k)
                   enddo !on k
                   do k=l,n
                      a(j,k)=a(j,k)+s*rv1(k)
                   enddo !on k
                enddo !on j
                do k=l,n
                   a(i,k)=scale*a(i,k)
                enddo !on k
             endif
          endif
          anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
       enddo !on i
       do i=n,1,-1
          if(i.lt.n)then
             if(g.ne.0.)then
                do j=l,n
                   v(j,i)=(a(i,j)/a(i,l))/g
                enddo !on j
                do j=l,n
                   s=0.
                   do k=l,n
                      s=s+a(i,k)*v(k,j)
                   enddo !on k
                   do k=l,n
                      v(k,j)=v(k,j)+s*v(k,i)
                   enddo !on k
                enddo !on j
             endif
             do j=l,n
                v(i,j)=0.
                v(j,i)=0.
             enddo !j
          endif
          v(i,i)=1.
          g=rv1(i)
          l=i
       enddo !on i
       do i=min(m,n),1,-1
          l=i+1
          g=w(i)
          do j=l,n
             a(i,j)=0.
          enddo ! on j
          if(g.ne.0.) then
             g=1./g
             do j=l,n
                s=0.
                do k=l,m
                   s=s+a(k,i)*a(k,j)
                enddo !on k
                f=(s/a(i,i))*g
                do k=i,m
                   a(k,j)=a(k,j)+f*a(k,i)
                enddo !on k
             enddo !on j
             do j=i,m
                a(j,i)=a(j,i)*g
             enddo !j
          else
             do j=i,m
                a(j,i)=0.       
             enddo !on j
          endif
          a(i,i)=a(i,i)+1.
       enddo ! on i
       do k=n,1,-1
          do its=1,30
             do l=k,1,-1
                nm=l-1
                if((abs(rv1(l))+anorm).eq.anorm) go to 2
                if((abs(w(nm)+anorm)).eq.anorm) go to 1
             enddo !on l
 1           c=0.
             s=1.
             do i=l,k
                f=s*rv1(i)
                rv1(i)=c*rv1(i)
                if((abs(f)+anorm).eq.anorm) go to 2
                g=w(i)
                h=pythag(f,g)
                w(i)=h
                h=1./h
                c=(g*h)
                s=-(f*h)
                do j=1,m
                   y=a(j,nm) 
                   z=a(j,i)       
                   a(j,nm)=(y*c)+(z*s)                 
                   a(j,i)=-(y*s)+(z*c)
                enddo !on j
             enddo !on i 
 2           z=w(k)  
             if(l.eq.k) then
                if(z.lt.0.) then
                   w(k)=-z
                   do j=1,n
                      v(j,k)=-v(j,k)
                   enddo ! on j
                endif
                go to 3
             endif
             if(its.eq.30) then
                print *,
     &       ' No convergence in sudcmp: inverted matrix may be suspect'
             endif
             x=w(l)
             nm=k-1
             y=w(nm)
             g=rv1(nm) 
             h=rv1(k)
             f=((y-z)*(y+z)+(g+h)*(g-h))/(2.*h*y)
             g=pythag(f,1.)
             f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
             c=1.
             s=1.
             do j=l,nm
                i=j+1
                g=rv1(i)
                y=w(i)
                h=s*g
                g=c*g
                z=pythag(f,h)
                rv1(j)=z
                c=f/z
                s=h/z
                f=(x*c)+(g*s)
                g=-(x*s)+(g*c)
                h=y*s
                y=y*c
                do jj=1,n
                   x=v(jj,j)
                   z=v(jj,i)
                   v(jj,j)=(x*c)+(z*s)
                   v(jj,i)=-(x*s)+(z*c)
                enddo !on jj
                z=pythag(f,h)
                w(j)=z
                if(z.ne.0.) then
                   z=1./z
                   c=f*z
                   s=h*z
                endif
                f=(c*g)+(s*y)
                x=-(s*g)+(c*y)
                do jj=1,m
                   y=a(jj,j)
                   z=a(jj,i)
                   a(jj,j)=(y*c)+(z*s)
                   a(jj,i)=-(y*s)+(z*c)
                enddo !on jj
             enddo !on j
             rv1(l)=0.
             rv1(k)=f
             w(k)=x
          enddo !on its
 3        continue
       enddo ! on k
c
       return
       end
c
c
        Subroutine svbksb(u,w,v,m,n,mp,np,b,x,nmax)
c
c*********************************************************************
c
c     This subroutine solves A.X=B for X where A is specified by
c     arrays u,w,v from svdcmp.
c
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       21 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c       14 Dec 1999  Peter Stamus, NOAA/FSL
c          Pass in max number of stations as nmax.
c
c*********************************************************************
c
        integer i,j,jj
        integer m,mp,n,np,nmax
        real s,tmp(nmax)
        real b(mp),u(mp,np),v(np,np),w(np),x(np)
c
        do j=1,n
           s=0.
           if(w(j).ne.0.) then
              do i=1,m
                 s=s+u(i,j)*b(i)
              enddo !i
              s=s/w(j)
           endif
           tmp(j)=s
        enddo !j
        do j=1,n
           s=0.
           do jj=1,n
              s=s+v(j,jj)*tmp(jj)
           enddo !jj
           x(j)=s
        enddo !on j
c
        return
        end
c
c     
      subroutine ludcmp(a,n,np,indx,d,nmax)  
c
c*********************************************************************
c
c     This subroutine performs lower/upper decomposition of a
c     matrix "a".
c
c     Original: John McGinley, NOAA/FSL  December 1999
c     Changes:
c       14 Dec 1999  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c
c*********************************************************************
c
      real d,a(np,np),tiny
      real aamax,dum,sum,vv(nmax)
      integer i,imax,j,k
      integer n,np,indx(n),nmax
      parameter (tiny=1.e-20)
c
      d=1.
      do i=1,n
         aamax=0.
         do j=1,n
            if(abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
         enddo !j
         if(aamax.eq.0.) then
            print *, 'All 0, singular matrix in ludcmp'
         endif
         vv(i)=1./aamax
      enddo !i
c
      do j=1,n
         do i=1,j-1
            sum=a(i,j)
            do k=1,i-1
               sum=sum-a(i,k)*a(k,j)
            enddo !k
            a(i,j)=sum
         enddo !i
         aamax=0.
         do i=j,n
            sum=a(i,j)
            do k=1,j-1
               sum=sum-a(i,k)*a(k,j)
            enddo !k
            a(i,j)=sum
            dum=vv(i)*abs(sum)
            if(dum.ge.aamax) then
               imax=i
               aamax=dum
            endif
         enddo !i
         if(j.ne.imax) then
            do k=1,n
               dum=a(imax,k)
               a(imax,k)=a(j,k)
               a(j,k)=dum
            enddo !k
            d=-d
            vv(imax)=vv(j)
         endif
         indx(j)=imax
         if(a(j,j).eq.0.)a(j,j)=tiny
         if(j.ne.n)then
            dum=1./a(j,j)
            do i=j+1,n
               a(i,j)=a(i,j)*dum
            enddo !i
         endif
      enddo !j
c
      return
      end
c
c
      subroutine lubksb(a,n,np,indx,b)   
c
c*********************************************************************
c
c     This subroutine solves a set of 'n' linear equationsns (a.x=b).
c
c     Original: John McGinley, NOAA/FSL  December 1999
c     Changes:
c       14 Dec 1999  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c
c*********************************************************************
c
      integer i,ii,j,ll
      integer n,np,indx(n)
      real a(np,np),b(n)  
      real sum
c
      ii=0
      do i=1,n
         ll=indx(i)
         sum=b(ll)
         b(ll)=b(i)
         if(ii.ne.0) then
            do j=ii,i-1
               sum=sum-a(i,j)*b(j)
            enddo !j
         elseif(sum.ne.0.) then
            ii=i
         endif
         b(i)=sum
      enddo !i
c
      do i=n,1,-1
         sum=b(i)
         do j=i+1,n
            sum=sum-a(i,j)*b(j)
         enddo !j
         b(i)=sum/a(i,i)
      enddo !i
c
      return
      end

