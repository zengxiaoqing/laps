SUBROUTINE WriteAnalysis(a,n)

!==========================================================
!  This routine is to write the analysis out in LAPS format.
!
!  HISTORY: MAY. 2004 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  REAL, PARAMETER :: temp0 = 273.16,capr = 287.0, cp=1004.0
  INTEGER, INTENT(IN) :: n(4)
  REAL,    INTENT(IN) :: a(n(1),n(2),n(3),n(4))
  REAL :: te,ept	! Theta_e
  REAL :: dxy(nx,ny,2),flu(n(1)),flv(n(2))
  REAL :: mjohn,r_missing_data

  CHARACTER :: ext*3,varnames(nvlaps+10)*3,vunits(nvlaps+10)*3
  CHARACTER :: lvl_coord(nvlaps+10)*3,comment(nvlaps+10)*125
  INTEGER   :: lvl(nvlaps+10),istatus,it,i,j
  REAL      :: data(n(1)-2*nfic,n(2)-2*nfic,16)
  REAL	    :: pp(n(1),n(2))

  CALL get_directory('lsx', dir_s, len)
  ext = 'lsx'

  PRINT*,'Data directory: ',dir_s(1:len)

  !PRINT*,'MAX U: ',MAXVAL(a(1:n(1),1:n(2),n(3),2))
  !PRINT*,'MIN U: ',MINVAL(a(1:n(1),1:n(2),n(3),2))

  ! Variable names and units:
  varnames(1) = 'T  '	! Temperature
  vunits(1) = 'K  '
  varnames(2) = 'U  '	! U
  vunits(2) = 'M/S'
  varnames(3) = 'V  '	! V
  vunits(3) = 'M/S'
  varnames(4) = 'VIS' 	! Visibility (m)
  ! comment(4) = ''
  vunits(4) = 'M'
  varnames(5) = 'TD '	! Dew point
  vunits(5) = 'K  '
  varnames(6) = 'P  '	! Reduced Pressure
  vunits(6) = 'PA  '
  comment(6)(1:21) = '0  M REDUCED PRESSURE'
  varnames(7) = 'THE'	! Theta_e
  vunits(7) = 'K  '
  varnames(8) = 'DIV'	! Divergence
  vunits(8) = '/S '
  varnames(9) = 'TH '	! Potential temperature
  vunits(9) = 'K  '
  varnames(10) = 'MRC '	! Moisture convergence
  vunits(10) = '/S  '
  varnames(11) = 'PP  '	! Reduced Pressure change
  vunits(11) = 'PA  '
  comment(11)(1:28) = '0  M REDUCED PRESSURE CHANGE'

  DO it=1,nvlaps+9
     lvl(it) = 0
     lvl_coord(it) = 'AGL'
  ENDDO

  PRINT*,'N V n(4) = ',n(4)

  dxy = 5000.0

  DO it=n(3)-1,n(3)

     ! gridpoints/meter:
     mjohn = 1.0/SQRT(grid_spacingx)

     ! Divergence:
     data(1:nx,1:ny,n(4)+2) = &
	(a(nfic+2:n(1)-nfic+1,nfic+1:n(2)-nfic,it,2)- &
         a(nfic  :n(1)-nfic-1,nfic+1:n(2)-nfic,it,2))/grid_spacingx*0.5+ &
        (a(nfic+1:n(1)-nfic,nfic+2:n(2)-nfic+1,it,3)- &
         a(nfic+1:n(1)-nfic,nfic  :n(2)-nfic-1,it,3))/grid_spacingy*0.5

     ! Convert to real*4:
     data(1:nx,1:ny,1:n(4)) = &
         a(1+nfic:n(1)-nfic,1+nfic:n(2)-nfic,it,1:n(4))

     ! Theta_e: equivalent potential temperature
     DO j=1,ny
	DO i=1,nx
	   data(i,j,n(4)+1) = ept(data(i,j,1)-temp0, &
                                  data(i,j,5)-temp0, &
                                  data(i,j,4)/100.0)+temp0
	ENDDO
     ENDDO

     ! Theta:
     data(1:nx,1:ny,n(4)+3) = data(1:nx,1:ny,1)* &
	                (100000.0/data(1:nx,1:ny,4))**(capr/cp)

     ! Flux convergence:
     data(1:nx,1:ny,16) = data(1:nx,1:ny,6)/100.0
     CALL meso_anl(data(1,1,2),data(1,1,3),data(1,1,16),data(1,1,1), &
	           data(1,1,5),data(1,1,n(4)+3),dxy(1,1,1),dxy(1,1,2), &
	           data(1,1,n(4)+5),data(1,1,n(4)+4),data(1,1,n(4)+6), &
                   data(1,1,n(4)+7),data(1,1,n(4)+8),nx,ny)

     ! Reduced Pressure Change:
     ! data(1:nx,1:ny,11) = 0.0
     CALL GridBarnes(a(1,1,it,6),n,n,pp)
     data(1:nx,1:ny,11) = a(1+nfic:n(1)-nfic,1+nfic:n(2)-nfic,it,6)-&
			  pp(1+nfic:n(1)-nfic,1+nfic:n(2)-nfic)

     CALL write_laps_data(istarttime+(it-1)*laps_cycle_time, &
                          dir_s,ext,nx,ny,n(4)+5,n(4)+5,varnames, &
                          lvl,lvl_coord,vunits,comment,data,istatus)
     write(6,*)' LSX file write completed, istatus = ',istatus,nx,ny

  ENDDO
  
END SUBROUTINE
