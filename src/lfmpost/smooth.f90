!dis   
!dis    Open Source License/Disclaimer, Forecast Systems Laboratory
!dis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
!dis    
!dis    This software is distributed under the Open Source Definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis    
!dis    In particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis    
!dis    - Redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis    
!dis    - Redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis    
!dis    - All modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis    
!dis    - If significant modifications or enhancements are made to this
!dis    software, the FSL Software Policy Manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis    
!dis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
!dis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
!dis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
!dis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
!dis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
!dis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
!dis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
!dis   
!dis 

SUBROUTINE smooth(field, ix, iy, iz, smth)

!
! *** Subprogram:  smooth - Smooth a meteorological field.
!     Author:  Stan Benjamin
!     Date  :  90-06-15
!
! *** Abstract:  Shapiro smoother.
!
! *** Program history log:      
!        85-12-09  S. Benjamin - Original version
!        96-06-16  J. Snook    - Modified to do 3d RAMS fields
!                              - hold array is dynamically allocated
!        01-03-07  B. Shaw     - Adapted to free-form F90 as part o
!                                MM5POST. Includes use of array syntax.
!
! *** Usage:  call smooth(field,ix,iy,iz,smth)
!
! *** Input argument list:
!        field    - real array  field(ix,iy,iz)
!                               Meteorological field
!        ix       - integer     x coordinates of field
!        iy       - integer     y coordinates of field
!        iz       - integer     z coordinates of field
!        smth     - real
!
! *** Output argument list:
!        field    - real array  field(ix,iy,iz)
!                               Smoothed meteorological field
!
! *** Remarks:  Reference:  Shapiro, 1970: "Smoothing, filtering, and
!        boundary effects", Rev. Geophys. Sp. Phys., 359-387.
!
!     This filter is of the type
!        z(i) = (1-s)z(i) + s(z(i+1)+z(i-1))/2
!     for a filter which is supposed to damp 2dx waves completely
!     but leave 4dx and longer with little damping,
!     it should be run with 2 passes using smth (or s) of 0.5
!     and -0.5 
!----------------------------------------------------------------------------
  IMPLICIT NONE

  ! Argument list

  INTEGER, INTENT(IN)         :: ix, iy, iz
  REAL, INTENT(INOUT)         :: field(ix,iy,iz)
  REAL, INTENT(IN)            :: smth

  ! Local variables
  INTEGER                     :: i, j, k, i1, i2, it
  REAL                        :: hold(ix,2)
  REAL                        :: smth1, smth2, smth3, smth4, &
                                 smth5, sum1, sum2
!----------------------------------------------------------------------------
  smth1=0.25*smth*smth
  smth2=0.50*smth*(1.-smth)
  smth3=(1.-smth)*(1.-smth)    
  smth4=(1.-smth)
  smth5=0.5*smth 

  DO k = 1, iz
   
    hold(:,1:2) = 0.
    i1 = 2
    i2 = 1
 
    DO j=2,iy-1
      it=i1
      i1=i2
      i2=it
      DO i=2,ix-1
        sum1=field(i-1,j+1,k)+field(i-1,j-1,k) &
            +field(i+1,j+1,k)+field(i+1,j-1,k)
        sum2=field(i  ,j+1,k)+field(i+1,j  ,k) &
            +field(i  ,j-1,k)+field(i-1,j  ,k)
        hold(i,i1)=smth1*sum1+smth2*sum2+smth3*field(i,j,k)
      ENDDO
      IF (j .ne. 2) field(2:ix-1,j-1,k)=hold(2:ix-1,i2)
    ENDDO

    field(2:ix-1,iy-1,k)=hold(2:ix-1,i1)

    field(2:ix-1,1,k)=smth4*field(2:ix-1,1,k) + &
                      smth5*(field(1:ix-2,1,k)+field(3:ix,1,k)) 
    field(2:ix-1,iy,k)=smth4*field(2:ix-1,iy,k) + &
                       smth5*(field(1:ix-2,iy,k)+field(3:ix,iy,k))

    field(1,2:iy-1,k) =  smth4*field(1,2:iy-1,k) + &
                         smth5*(field(1,1:iy-2,k)+field(1,3:iy,k)) 
  
  ENDDO
  RETURN
END SUBROUTINE SMOOTH              
