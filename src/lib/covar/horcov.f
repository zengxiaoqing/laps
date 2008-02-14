c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine hor_cov(i4time,v_inno
     +                         ,num_of_ens
     +                         ,ix,iy,iz,isoLevel
     +                         ,dx,dy
     +                         ,ilength,jlength
     +                         ,covhor)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c c History:
c Date         Name          Action
c --------     -----------   -------------------------------------------
c 03/05/2007   Ok-Yeon Kim   Created. 
c ----------------------------------------------------------------------

      implicit none

      integer       ix,iy,iz
      integer       num_of_ens,iens
      integer       i,j,k,l,m,n,ii,jj
      integer       ilength,jlength
      integer       ib,ie,jb,je,inc
      integer     i4time
      integer     isoLevel(iz)
      real          dx,dy
      real        v_inno(ix,iy,iz,20)
      real        amean(ix,iy,iz)
      real        dev(ix,iy,iz,20)
      real        var(ix,iy,iz)
      real        covhor(9,9,ix,iy,iz)

c--------------- End of Diclaration -----------------------------------


c Average of innovation at each point  !!
      amean=0.
      do n=1,num_of_ens
        do k=1,iz
          do j=1,iy
            do i=1,ix
              amean(i,j,k)=amean(i,j,k)+v_inno(i,j,k,n)
            enddo 
          enddo
        enddo
      enddo

      do k=1,iz
        do j=1,iy
          do i=1,ix
            amean(i,j,k)=amean(i,j,k)/float(num_of_ens)
          enddo
        enddo
      enddo


c Deviation, variance at each point  !!
      var=0.
      do k=1,iz
        do j=1,iy  
          do i=1,ix
            do n=1,num_of_ens
               dev(i,j,k,n)=v_inno(i,j,k,n)-amean(i,j,k)  ! deviation
               var(i,j,k)=var(i,j,k)+(dev(i,j,k,n)**2.)   ! variance
            enddo
          enddo
        enddo
      enddo


c Covariance in the horizontal length  !!

      iens=num_of_ens-1            
      do k=1,iz
        do j=1,iy
          do i=1,ix

            if((i.le.ilength).and.(j.le.jlength))then    ! horizontal length (radius)
                 ib=i-(i-1)
                 ie=ib+2*ilength
                 jb=j-(j-1)
                 je=jb+2*jlength
            elseif((i.le.ilength).and.(j.gt.iy-jlength))then
                 ib=i-(i-1)
                 ie=ib+2*ilength
                 jb=iy
                 je=iy-2*jlength 
            elseif((i.gt.ix-ilength).and.(j.le.jlength))then
                 ib=ix-2*ilength
                 ie=ix
                 jb=j-(j-1)
                 je=jb+2*jlength
            elseif((i.gt.ix-ilength).and.(j.gt.iy-jlength))then
                 ib=ix-2*ilength
                 ie=ix
                 jb=iy
                 je=iy-2*jlength
            elseif((i.le.ilength).and.
     +             (j.gt.jlength.and.j.le.iy-jlength))then
                 ib=i-(i-1)
                 ie=ib+2*ilength
                 jb=j-jlength
                 je=j+jlength
            elseif((i.ge.ix-ilength).and.
     +             (j.gt.jlength.and.j.le.iy-jlength))then
                 ib=ix-2*ilength
                 ie=ix
                 jb=j-jlength
                 je=j+jlength
            elseif((i.gt.ilength.and.i.le.ix-ilength).and.
     +             (j.le.jlength))then
                 ib=i-ilength
                 ie=i+ilength
                 jb=j-(j-1)
                 je=jb+2*jlength
            elseif((i.gt.ilength.and.i.le.ix-ilength).and. 
     +             (j.gt.iy-jlength))then
                 ib=i-ilength
                 ie=i+ilength
                 jb=iy
                 je=iy-2*jlength
            elseif((i.gt.ilength.and.i.le.ix-ilength).and.
     +             (j.gt.jlength.and.j.le.iy-jlength))then
                 ib=i-ilength
                 ie=i+ilength
                 jb=j-jlength
                 je=j+jlength
            endif


      do l=jb,je 
        do m=ib,ie
           do n=1,num_of_ens
              covhor(m,l,i,j,k)=covhor(m,l,i,j,k)
     +                         +(dev(i,j,k,n)*dev(m,l,k,n))
           enddo     ! n
              covhor(m,l,i,j,k)=covhor(m,l,i,j,k)/float(iens)
        enddo        ! m
      enddo          ! l

           enddo     ! i
        enddo        ! j
      enddo          ! k


      return
      end
