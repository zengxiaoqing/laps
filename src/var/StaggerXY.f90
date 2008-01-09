
SUBROUTINE StaggerXY_3D(vin,nx,ny,nz,onx,ony,onz,vout)

!==========================================================
!  This routine computes a stagger grid 
!  using a given horizontal uniform grid 
!
!  Input:
!	vin:    Data prepare to stagger 	
!	nx:	X grid point number before stagger
!	ny:	Y grid point number before stagger
!	nz:	Z grid point number before stagger
!       onx:    X grid point number after stagger
!       ony:    Y grid point number after stagger
!       onz:    Z grid point number after stagger  
!
!  Output:
!	vout:   Data after stagger	
!
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,nz,onx,ony,onz
  REAL, INTENT(IN) ::    vin(nx,ny,nz)
  REAL, INTENT(OUT) ::   vout(onx,ony,onz)

  ! Local variables:
    
  REAL :: vinx(nx-1,ny,nz),viny(nx,ny-1,nz),vinyy(nx-1,ny-1,nz)

  ! Linear interpolation:
  
    if(nx /= onx .and. ny /= ony .and. onz /= 1) then 
 
        vinx(1:nx-1,1:ny,1:nz) = 0.5*( &
             vin(1:nx-1,1:ny,1:nz)+ &
             vin(2:nx  ,1:ny,1:nz))
        vinyy(1:nx-1,1:ny-1,1:nz) = 0.5*( &
             vinx(1:nx-1,1:ny-1,1:nz)+ &
             vinx(1:nx-1,2:ny  ,1:nz))

        vout(1:onx,1:ony,1:nz) = vinyy(1:nx-1,1:ny-1,1:nz) 

    elseif(nx /= onx .and. ny /= ony .and. onz == 1) then

        vinx(1:nx-1,1:ny,1) = 0.5*( &
             vin(1:nx-1,1:ny,1)+ &
             vin(2:nx  ,1:ny,1))
        vinyy(1:nx-1,1:ny-1,1) = 0.5*( &
             vinx(1:nx-1,1:ny-1,1)+ &
             vinx(1:nx-1,2:ny  ,1))

        vout(1:onx,1:ony,1) = vinyy(1:nx-1,1:ny-1,1)
      
    elseif(nx /= onx .and. ny == ony) then

        vinx(1:nx-1,1:ny,1:nz) = 0.5*( &
             vin(1:nx-1,1:ny,1:nz)+ &
             vin(2:nx  ,1:ny,1:nz))

        vout(1:onx,1:ony,1:nz) = vinx(1:nx-1,1:ny,1:nz)

    elseif(nx == onx .and. ny /= ony) then

        viny(1:nx,1:ny-1,1:nz) = 0.5*( &
             vin(1:nx,1:ny-1,1:nz)+ &
             vin(1:nx,2:ny  ,1:nz))

        vout(1:onx,1:ony,1:nz) = viny(1:nx,1:ny-1,1:nz)

    endif
 

END SUBROUTINE StaggerXY_3D 


SUBROUTINE StaggerXY_2D(vin,nx,ny,onx,ony,vout)

!==========================================================
!  This routine computes a stagger grid 
!  using a given horizontal uniform grid 
!
!  Input:
!	vin:    Data prepare to stagger 	
!	nx:	X grid point number before stagger
!	ny:	Y grid point number before stagger
!       onx:    X grid point number after stagger
!       ony:    Y grid point number after stagger  
!
!  Output:
!	vout:   Data after stagger	
!
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,onx,ony
  REAL, INTENT(IN) ::    vin(nx,ny)
  REAL, INTENT(OUT) ::   vout(onx,ony)

  ! Local variables:
    
  REAL :: vinx(nx-1,ny),viny(nx,ny-1),vinyy(nx-1,ny-1)

  ! Linear interpolation:
  
    if(nx /= onx .and. ny /= ony) then 
 
        vinx(1:nx-1,1:ny) = 0.5*( &
             vin(1:nx-1,1:ny)+ &
             vin(2:nx  ,1:ny))
        vinyy(1:nx-1,1:ny-1) = 0.5*( &
             vinx(1:nx-1,1:ny-1)+ &
             vinx(1:nx-1,2:ny))

        vout(1:onx,1:ony) = vinyy(1:nx-1,1:ny-1) 
      
    elseif(nx /= onx .and. ny == ony) then

        vinx(1:nx-1,1:ny) = 0.5*( &
             vin(1:nx-1,1:ny)+ &
             vin(2:nx  ,1:ny))

        vout(1:onx,1:ony) = vinx(1:nx-1,1:ny)

    elseif(nx == onx .and. ny /= ony) then

        viny(1:nx,1:ny-1) = 0.5*( &
             vin(1:nx,1:ny-1)+ &
             vin(1:nx,2:ny))

        vout(1:onx,1:ony) = viny(1:nx,1:ny-1)

    endif
 

END SUBROUTINE StaggerXY_2D 


SUBROUTINE UntaggerXY_3D(vin,nx,ny,nz,onx,ony,vout)

!==========================================================
!  This routine computes a unstagger grid 
!  using a given horizontal stagger grid 
!
!  Input:
!	vin:    stagger Data prepare to unstagger 	
!	nx:	number of stagger grid in X 
!	ny:     number of stagger grid in Y	
!	nz:     Z grid point number after unstagger 	
!       onx:    X grid point number after unstagger
!       ony:    Y grid point number after unstagger  
!
!  Output:
!	vout:   Data after unstagger	
!
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,nz,onx,ony
  REAL, INTENT(IN) ::    vin(nx,ny,nz)
  REAL, INTENT(OUT) ::   vout(onx,ony,nz)

  ! Local variables:
    
  INTEGER :: i,j,k
  REAL :: vinx(onx,ony,nz),viny(onx,ony,nz),vinyy(nx,ony,nz)

  ! Linear interpolation:
  
    if(nx /= onx .and. ny /= ony) then 

       do k = 1,nz 
          do i=1,nx
             do j=2,ny
                   vinyy(i,j,k) = &
                   0.5*(vin(i,j-1,k)+ vin(i,j,k))
             enddo
             vinyy(i,1,k) = &
                   1.5*vin(i,1,k)-&
                   0.5*vin(i,2,k)
             vinyy(i,ony,k) = &
                   1.5*vin(i,ny,k)-&
                   0.5*vin(i,ny-1,k)
          enddo
       enddo

       do k = 1,nz
          do j=1,ony
             do i=2,nx
                vinx(i,j,k) = &
                0.5*(vinyy(i-1,j,k)+ vinyy(i,j,k))
             enddo
             vinx(1,j,k) = &
                1.5*vinyy(1,j,k)-&
                0.5*vinyy(2,j,k)
             vinx(onx,j,k) = &
                1.5*vinyy(nx,j,k)-&
                0.5*vinyy(nx-1,j,k)
          enddo
       enddo

       vout(1:onx,1:ony,1:nz) = vinx(1:onx,1:ony,1:nz) 
      
    elseif(nx /= onx .and. ny == ony) then

       do k = 1,nz
          do j=1,ony
             do i=2,nx
                vinx(i,j,k) = &
                0.5*(vin(i-1,j,k)+ vin(i,j,k))
             enddo
             vinx(1,j,k) = &
                1.5*vin(1,j,k)-&
                0.5*vin(2,j,k)
             vinx(onx,j,k) = &
                1.5*vin(nx,j,k)-&
                0.5*vin(nx-1,j,k)
          enddo
       enddo 

       vout(1:onx,1:ony,1:nz) = vinx(1:onx,1:ony,1:nz)

    elseif(nx == onx .and. ny /= ony) then

       do k = 1,nz
          do i=1,onx
             do j=2,ny
                   viny(i,j,k) = &
                   0.5*(vin(i,j-1,k)+ vin(i,j,k))
             enddo
             viny(i,1,k) = &
                   1.5*vin(i,1,k)-&
                   0.5*vin(i,2,k)
             viny(i,ony,k) = &
                   1.5*vin(i,ny,k)-&
                   0.5*vin(i,ny-1,k)
          enddo
       enddo 

       vout(1:onx,1:ony,1:nz) = viny(1:onx,1:ony,1:nz)

    endif
 
END SUBROUTINE UntaggerXY_3D 


SUBROUTINE UntaggerXY_2D(vin,nx,ny,onx,ony,vout)

!==========================================================
!  This routine computes a unstagger grid 
!  using a given horizontal stagger grid 
!
!  Input:
!	vin:    stagger Data prepare to unstagger 	
!	nx:	number of stagger grid in X 
!	ny:     number of stagger grid in Y	
!       onx:    X grid point number after unstagger
!       ony:    Y grid point number after unstagger  
!
!  Output:
!	vout:   Data after unstagger	
!
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,onx,ony
  REAL, INTENT(IN) ::    vin(nx,ny)
  REAL, INTENT(OUT) ::   vout(onx,ony)

  ! Local variables:
    
    INTEGER :: i,j
    REAL :: vinx(onx,ony),viny(onx,ony),vinyy(nx,ony)

  ! Linear interpolation:
  
    if(nx /= onx .and. ny /= ony) then 

          do i=1,nx
             do j=2,ny
                   vinyy(i,j) = &
                   0.5*(vin(i,j-1)+ vin(i,j))
             enddo
             vinyy(i,1) = &
                   1.5*vin(i,1)-&
                   0.5*vin(i,2)
             vinyy(i,ony) = &
                   1.5*vin(i,ny)-&
                   0.5*vin(i,ny-1)
          enddo

          do j=1,ony
             do i=2,nx
                vinx(i,j) = &
                0.5*(vinyy(i-1,j)+ vinyy(i,j))
             enddo
             vinx(1,j) = &
                1.5*vinyy(1,j)-&
                0.5*vinyy(2,j)
             vinx(onx,j) = &
                1.5*vinyy(nx,j)-&
                0.5*vinyy(nx-1,j)
          enddo

        vout(1:onx,1:ony) = vinx(1:onx,1:ony) 
      
    elseif(nx /= onx .and. ny == ony) then

          do j=1,ony
             do i=2,nx
                vinx(i,j) = &
                0.5*(vin(i-1,j)+ vin(i,j))
             enddo
             vinx(1,j) = &
                1.5*vin(1,j)-&
                0.5*vin(2,j)
             vinx(onx,j) = &
                1.5*vin(nx,j)-&
                0.5*vin(nx-1,j)
          enddo

        vout(1:onx,1:ony) = vinx(1:onx,1:ony)

    elseif(nx == onx .and. ny /= ony) then

          do i=1,onx
             do j=2,ny
                   viny(i,j) = &
                   0.5*(vin(i,j-1)+ vin(i,j))
             enddo
             viny(i,1) = &
                   1.5*vin(i,1)-&
                   0.5*vin(i,2)
             viny(i,ony) = &
                   1.5*vin(i,ny)-&
                   0.5*vin(i,ny-1)
          enddo
 
        vout(1:onx,1:ony) = viny(1:onx,1:ony)

    endif
 
END SUBROUTINE UntaggerXY_2D 


