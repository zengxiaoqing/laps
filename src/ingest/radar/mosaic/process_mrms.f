
      subroutine process_mrms(ni,nj,nk)

      real ref_3d(ni,nj,nk)
      real ref3d_column(nk+2,ni*nj)

!     ref3d_column=-999.0
!     numref=0
!     DO j=2,nlat-1
!     DO i=2,nlon-1
!     numlvl=0
!     DO k=1,maxlvl
!       if(abs(ref0(i,j,k)) < 888.0 ) numlvl=numlvl+1
!     ENDDO ! k
!     if(numlvl > 0 ) then
!       numref=numref+1
!       ref3d_column(1,numref)=float(i)
!       ref3d_column(2,numref)=float(j)
!       DO k=1,maxlvl
!          ref3d_column(2+k,numref)=ref0(i,j,k)
!       ENDDO
!     endif
!     ENDDO ! i
!     ENDDO ! j

!     Sample file location
!     /scratch2/portfolios/BMC/rtrr/RAPdev3/cycle/2015032715/obsprd/RefInGSI.dat

      maxlvl = nk
      OPEN(10,file='./'//'RefInGSI.dat',form='unformatted')
        read(10) maxlvl,nj,ni,numref,i1,i2
        write(*,*) 'Dump out results', numref, 'out of', nj,ni
        read(10) ((ref3d_column(k,i),k=1,maxlvl+2),i=1,numref)
      close(10)

      return
      end
