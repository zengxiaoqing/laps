
      subroutine compute_latlon(nnxp,nnyp,n_staggers
     + ,deltax,xtn,ytn,lats,lons,istatus)


      implicit none

      integer nnxp,nnyp
      integer n_staggers
      integer i,j,k,nc

      real    deltax,deltay
      real    mdlat,mdlon
      real    stagger_ew,stagger_ns
      real    erad

c these returned to gridgen
      real    xtn(nnxp,n_staggers)
      real    ytn(nnyp,n_staggers)

c these are the anchor points (SW corner) of the grid
      real    xmn1,ymn1

c     real    lat(nnxp,nnyp)
c     real    lon(nnxp,nnyp)

c A-c staggers contained within these arrays.
      real    lats(nnxp,nnyp,n_staggers)
      real    lons(nnxp,nnyp,n_staggers)

      character c_dataroot*255
      character c10_grid_fname*10

      integer istatus

      print*,'calculate lat/lon at stagger grid points.'

      deltay=deltax


      call get_grid_center(mdlat,mdlon,istatus)
      if(istatus .ne. 1)then
         write(6,*)' Error returned: get_grid_center'
         return
      endif

      call get_earth_radius(erad,istatus)
      if(istatus .ne. 1)then
         write(6,*)' Error calling get_earth_radius'
         return
      endif

      do k=1,n_staggers

         if(k.eq.1)then

c Get X/Y for lower left corner
            CALL POLAR_GP(mdlat,mdlon,XMN1,YMN1,DELTAX,DELTAY,
     +  NNXP,NNYP)
            stagger_ns=0
            stagger_ew=0
            call get_xytn(deltax,deltay,nnxp,nnyp,stagger_ew
     +,stagger_ns,xmn1,ymn1,xtn(1,k),ytn(1,k))
 
         elseif(k.eq.2)then    !this is A-stagger (.5 E-W stagger) 

            stagger_ew=0.5*deltax
            stagger_ns=0
            call get_xytn(deltax,deltay,nnxp,nnyp,stagger_ew
     +,stagger_ns,xmn1,ymn1,xtn(1,k),ytn(1,k))

         elseif(k.eq.3)then !this is the B-stagger (.5 N-S stagger)

            call get_grid_center(mdlat,mdlon,istatus)
            call get_xytn(deltax,deltay,nnxp,nnyp,0,0,xmn1,ymn1
     +,xtn(1,k),ytn(1,k))
            stagger_ew=0
            stagger_ns=0.5*deltay 
            call get_xytn(deltax,deltay,nnxp,nnyp,stagger_ew
     +,stagger_ns,xmn1,ymn1,xtn(1,k),ytn(1,k)) 

         elseif(k.eq.4)then !this is the C-stagger (.5 both N-S and E-W)

            call get_grid_center(mdlat,mdlon,istatus)
            call get_xytn(deltax,deltay,nnxp,nnyp,0,0,xmn1,ymn1
     +,xtn(1,k),ytn(1,k))
            stagger_ew=0.5*deltax
            stagger_ns=0.5*deltay
            call get_xytn(deltax,deltay,nnxp,nnyp,stagger_ew
     +,stagger_ns,xmn1,ymn1,xtn(1,k),ytn(1,k))

         endif
C
C*****************************************************************
C*  Convert it to lat/lon using the library routines.            *

         Do J = 1,nnyp
         Do I = 1,nnxp

            call xy_to_latlon(xtn(i,k),ytn(j,k),erad ! ,90.,std_lon
     1                                     ,lats(I,J,k),lons(I,J,k))

c             print *,'i,j,xtn,ytn,pla,lplo=',i,j,xtn,ytn,pla,plo
         enddo
         enddo

      enddo

      istatus = 1

      return
      end

c ----------------------------------------------

      subroutine get_xytn(deltax,deltay,nnxp,nnyp,stagger_ew
     +,stagger_ns,xmn1,ymn1,xtn,ytn)

      implicit none

      integer  nnxp,nnyp
      integer  i,j

      real    deltax,deltay
      real    stagger_ew,stagger_ns
      real    xmn1,ymn1
      real    xmn(nnxp)
      real    ymn(nnyp)
      real    xtn(nnxp)
      real    ytn(nnyp)

      xmn(1)=xmn1+stagger_ew
      ymn(1)=ymn1+stagger_ns

      DO 600 I=2,NNXP
         XMN(I)=XMN(I-1)+DELTAX
 600  CONTINUE
      XMN(NNXP)=2*XMN(NNXP-1)-XMN(NNXP-2)
 
      DO 610 J=2,NNYP
         YMN(J)=YMN(J-1)+DELTAY
 610  CONTINUE
      YMN(NNYP)=2*YMN(NNYP-1)-YMN(NNYP-2)
 
      DO 650 I=2,NNXP
         XTN(I)=.5*(XMN(I)+XMN(I-1))
 650  CONTINUE
      XTN(1)=1.5*XMN(1)-.5*XMN(2)

      DO 660 J=2,NNYP
         YTN(J)=.5*(YMN(J)+YMN(J-1))
 660  CONTINUE
      YTN(1)=1.5*YMN(1)-.5*YMN(2)

      return
      end
c
c -------------------------------------------------------------
c
      subroutine get_map_factor_grid(nx,ny,n_staggers
     +,rlats,rlons ,rmap_factors,istatus)

      implicit none

      integer nx,ny
      integer n_staggers
      integer i,j,k
      integer istatus

      real    rlats(nx,ny,n_staggers)
      real    rlons(nx,ny,n_staggers)
      real    rmap_factors(nx,ny,n_staggers)
      real    sigma

      do k=1,n_staggers
      do j=1,ny
      do i=1,nx
         call get_sigma(rlats(i,j,k),rlons(i,j,k)
     +,sigma,istatus)
         if(istatus.ne.1)goto 999
         rmap_factors(i,j,k)=sigma
      enddo
      enddo
      enddo

      return

999   print*,'Error returned: get_sigma'

      return
      end
c
c ------------------------------------------------------------
c
      subroutine get_coriolis_components(nx,ny,lat,coriolis_parms)
c
      include 'trigd.inc'
      implicit none

      integer nx,ny
      integer i,j

      real    lat(nx,ny)
      real    coriolis_parms(nx,ny,2)

      real    omega_ear
      data    omega_ear/7.292e-5/

      do j=1,ny
      do i=1,nx
         coriolis_parms(i,j,1)=2*omega_ear*sind(lat(i,j))
         coriolis_parms(i,j,2)=2*omega_ear*cosd(lat(i,j))
      enddo
      enddo

      return
      end
c
c ----------------------------------------------------------
c
      subroutine get_projrot_grid(nx,ny,lat,lon,projrot_grid
     +,istatus)

      include 'trigd.inc'
      implicit none

      integer nx,ny
      integer istatus
      integer i,j
      real    r_missing_data
      real    projrot_deg
      real    projrot_latlon
      real    projrot_grid(nx,ny,2)
      real    lat(nx,ny),lon(nx,ny)

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus.ne.1)then
         print*,'Error returned: get_r_missing_data'
         return
      endif

      do j=1,ny
      do i=1,nx
         projrot_grid(i,j,1)=r_missing_data
         projrot_grid(i,j,2)=r_missing_data
         projrot_deg=projrot_latlon(lat(i,j),lon(i,j),istatus)
         if(istatus.eq.1)then
            projrot_grid(i,j,1)=sind(projrot_deg)
            projrot_grid(i,j,2)=cosd(projrot_deg)
         endif
      enddo
      enddo

      return
      end
c
c--------------------------------------------------------
c
      subroutine get_static_albedo(nx,ny,lat,lon,landfrac
     +,static_albedo,istatus)

      implicit none

      integer nx,ny
      real    lat(nx,ny),lon(nx,ny)
      real    landfrac(nx,ny)
      real    static_albedo(nx,ny)
      real    water_albedo_cmn
      real    r_missing_data
      integer istatus
      integer i,j
      integer nwater

      data water_albedo_cmn/0.08/

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus.ne.1)then
         print*,'Error returned: get_r_missing_data'
         return
      endif

      nwater = 0
      do j=1,ny
      do i=1,nx
         if(landfrac(i,j).le.0.01)then
            nwater = nwater+1
            static_albedo(i,j)=water_albedo_cmn
         else
            static_albedo(i,j)=r_missing_data
         endif
      enddo
      enddo
      if(nwater .gt. 0)then
         print*,'static_albedo ',water_albedo_cmn,' used ',
     +'at ',nwater,' grid points '
      else
         print*,'No water grid points for water_albedo',
     +' in this domain'
      endif

      return
      end
c
c --------------------------------------------------------
c
      subroutine bilinear_interp(i,j,imax,jmax,array_2d,result)
c
c this used only for getting topography on the c-staggered grid

      implicit none
      integer i,j,imax,jmax
      real*4 result
      real*4 array_2d(imax,jmax)
      real*4 Z1,Z2,Z3,Z4
      real*4 fraci,fracj

      fraci = 0.5
      fracj = 0.5

      Z1=array_2d(i  , j  )
      Z2=array_2d(i-1, j  )
      Z3=array_2d(i-1, j-1)
      Z4=array_2d(i  , j-1)

      result= Z1+(Z2-Z1)*fraci+(Z4-Z1)*fracj
     1     - (Z2+Z4-Z3-Z1)*fraci*fracj

      return
      end

c
c--------------------------------------------------------
c
      subroutine get_gridgen_var(nf,ngrids,var,comment)

      implicit none

      integer        nf,ngrids,i,j
      character*(*)  var(nf)
      character*(*)  comment(nf)
      character*2    cat

      if(ngrids.eq.26)then

         var(1)    = 'LAT'
         var(2)    = 'LON'
         var(3)    = 'AVG'
         var(4)    = 'LDF'
         var(5)    = 'USE'
         var(6)    = 'ALB'
         var(7)    = 'STD'
         var(8)    = 'SLN'
         var(9)    = 'SLT'
         var(10)   = 'STL'
         var(11)   = 'SBL'
         var(12)   = 'LND'  ! Land-Water Mask based on USGS landuse
         i=12
         do j=1,12
            write(cat,'(i2.2)')j
            if(cat(1:1).eq.' ')cat(1:1)='0'
            if(cat(2:2).eq.' ')cat(2:2)='0'
            var(i+j)= 'G'//cat   ! vegetation greenness fraction
         enddo

         var(25)='TMP'
         var(ngrids)   = 'ZIN'
 
         comment(1) = 'Lat: From MODEL by J. Snook/ S. Albers 1-95\0'
         comment(2) = 'Lon: From MODEL by J. Snook/ S. Albers 1-95\0'
         comment(3) = 'Average terrain elevation (m) \0'
         comment(4) = 'Land Fraction \0'
         comment(5) = 'Land Use dominant category (USGS 24 Category) \0'
         comment(6) = 'Clear Sky Albedo - fixed at .04 over water\0'
         comment(7) = 'Standard Deviation of Elevation data (m)\0'
         comment(8) = 'Mean longitudinal terrain slope (m/m)\0'
         comment(9) = 'Mean latitudinal terrain slope (m/m)\0'
         comment(10)= 'Top layer (0-30cm) dominant category soiltype\0'
         comment(11)= 'Bot layer (30-90cm) dominant category soiltype\0'
         comment(12)= 'Land-Water Mask (0=water; 1 otherwise) \0'

         i=12
         do j=1,12
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          comment(i+j)= 'vegetation greenness fraction: mon = '//cat
         enddo

         comment(25)='Mean Annual Soil Temp (deg K)'
         comment(26)='\0'

      elseif(ngrids.eq.97)then

         var(1)    = 'LAT'  ! non-staggered (Analysis-grid) lats
         var(2)    = 'LON'  ! non-staggered (Analysis-grid) lons
         var(3)    = 'LAA'  ! a-stagger (.5*deltax (e-w)) lats
         var(4)    = 'LOA'  ! a-stagger (.5*deltax (e-w)) lons
         var(5)    = 'LAB'  ! b-stagger (.5*deltay (n-s)) lats
         var(6)    = 'LOB'  ! b-stagger (.5*deltay (n-s)) lons
         var(7)    = 'LAC'  ! c-stagger (.5*deltax and .5*deltay) lats
         var(8)    = 'LOC'  ! c-stagger (.5*deltay and .5*deltay) lons

         var(9)    = 'AVG'  ! Topo (m) on Analysis-grid

         var(10)   = 'LDF'  ! Land Fraction
         var(11)   = 'USE'  ! landuse dominant category
         var(12)   = 'LND'  ! Land-Water Mask based on USGS landuse
         var(13)   = 'STL'  ! Soiltype top layer dominant category
         var(14)   = 'SBL'  ! Soiltype bot layer dominant category

         i=14
         do j=1,24
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          var(i+j)= 'U'//cat
         enddo

         i=39
         var(i)     = 'SPR'  ! Sin(projection rotation) from true
         var(i+1)   = 'CPR'  ! Cos(projection rotation) from true
         var(i+2)   = 'MFL'  ! Map factor Analysis grid
         var(i+3)   = 'MFA'  ! Map factor a-stagger grid
         var(i+4)   = 'MFB'  ! Map factor b-stagger grid
         var(i+5)   = 'MFC'  ! Map factor c-stagger grid
         var(i+6)   = 'CPH'  ! Horizontal component of coriolis parameter
         var(i+7)   = 'CPV'  ! Vertical component of coriolis parameter

         var(i+8)   = 'ALB'  ! Static (climatological) albedo
         var(i+9)   = 'STD'  ! Standard Deviation of Elevation Data (m)
         var(i+10)  = 'SLN'  ! Terrain Slope; Longitudinal Component (m/m)
         var(i+11)  = 'SLT'  ! Terrain Slope; Latitudinal Component (m/m)
         var(i+12)  = 'AVC'  ! Topo (m) on c-stagger grid

         i=51
         do j=1,16
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          var(i+j)= 'T'//cat   !top layer (0-30cm) soiltype (% dist)
         enddo
         i=67
         do j=1,16
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          var(i+j)= 'B'//cat   ! bot layer (30-90cm) soiltype (% dist)
         enddo
         i=83
         do j=1,12
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          var(i+j)= 'G'//cat   ! vegetation greenness fraction
         enddo

         var(96)='TMP'

         var(ngrids)   = 'ZIN'

         comment(1) = 'Made from MODEL by J. Snook/ S. Albers 1-95\0'
         comment(2) = 'Made from MODEL by J. Snook/ S. Albers 1-95\0'
         comment(3) = 'a-stagger grid latitudes \0'
         comment(4) = 'a-stagger grid longitudes for WRF_SI \0'
         comment(5) = 'b-stagger grid latitudes for WRF_SI \0'
         comment(6) = 'b-stagger grid longitudes for WRF_SI \0'
         comment(7) = 'c-stagger grid latitudes for WRF_SI \0'
         comment(8) = 'c-stagger grid longitudes for WRF_SI \0'
         comment(9) = 'Average terrain elevation (m) \0'
         comment(10)= 'Land Fraction A-grid \0'
         comment(11)= 'Land Use Dominant category \0'
         comment(12)= 'Land-Water Mask (0=water; 1 otherwise) \0'
         comment(13)= 'Soiltype Top Layer dominant category \0'
         comment(14)= 'Soiltype Bot Layer dominant category \0'

         i=14
         do j=1,24
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          var(i+j)= 'U'//cat
          comment(i+j)= '% Dist Land Use Category '//cat//' \0'
         enddo

         i=39
         comment(i)= 'Sin of projection rotation (rad) \0'
         comment(i+1)= 'Cosine of projection rotation (rad) \0'
         comment(i+2)= 'Map Factor Analysis grid \0'
         comment(i+3)= 'Map Factor a-stagger grid \0'
         comment(i+4)='Map Factor b-stagger grid \0'
         comment(i+5)= 'Map Factor c-stagger grid \0'
         comment(i+6)= 'Horizontal component coriolis parameter \0'
         comment(i+7)= 'Vertical component coriolis parameter \0'

         comment(i+8)= 'Static Albedo (%) valid only over water atm \0'
         comment(i+9)= 'Standard Deviation of Elevation data (m)\0'
         comment(i+10)= 'Mean longitudinal terrain slope (m/m)\0'
         comment(i+11)= 'Mean latitudinal terrain slope (m/m)\0'
         comment(i+12)= 'Average terrain elevation (c-stagger) (m)\0'

         i=51
         do j=1,16
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          comment(i+j)='% Dist Top Layer Soiltype Category '//cat
         enddo

         i=67
         do j=1,16
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          comment(i+j)= '% Dist Bot Layer Soiltype Category '//cat
         enddo

         i=83
         do j=1,12
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          comment(i+j)= 'vegetation greenness fraction: mon = '//cat
         enddo

         comment(96)='1 degree mean annual soiltemp (deg K)'

         comment(ngrids)= '\0'

      endif
      return
      end

c ********************************************************************

	subroutine read_dem(unit_no,unit_name,nn1,nn2,i1,i2
     &,data)
	implicit none

	integer countx,county,unit_no,nn1,nn2

        character  cdata(nn1,nn2)*2
        integer    idata(nn1,nn2)
	real        data(nn1,nn2)
	integer len, lend, i1, i2, ia
        real multiplier
c       logical l1,l2
	character*(*) unit_name
        character*1   ctiletype

C	open(unit_no,file=unit_name,status='old',access='direct',
C	. recl=nn2*nn1*2)
C	inquire(unit_no,exist=l1,opened=l2)
C	read(unit_no,rec=1) idata

	call s_len(unit_name,len)
        call get_directory_length(unit_name,lend)
        ctiletype=unit_name(lend+1:lend+1)

        multiplier=1.0
	if(ctiletype.eq.'T'.or.ctiletype.eq.'U')then
           call read_binary_field(cdata,i1,i2,nn1*nn2,unit_name,len)
           do county=1,nn2
           do countx=1,nn1
              idata(countx,county) = ia (cdata(countx,county),2,0)
           enddo
           enddo
           if(ctiletype.eq.'T')multiplier=.01  !(for T data these are temps * 100)
        else
           call read_binary_field(idata,i1,i2,nn1*nn2,unit_name,len)
        endif

	do county=1,nn2
	do countx=1,nn1
	 if(idata(countx,county).eq.-9999) idata(countx,county)=0
	  data(countx,county)=float(idata(countx,nn2-county+1))*multiplier
c SG97 initial data (DEM format) starts in the lower-left corner;
c SG97 this format is wrapped around to have upper-left corner as its start.
c
c JS00 some machines do not account for signed integers
	   if(data(countx,county).ge.15535.0)
     &data(countx,county)=data(countx,county)-65535

	enddo
	enddo

ccc	 close(unit_no)
	return
	end
c ********************************************************************

        subroutine read_dem_g(unit_no,unit_name,nn1,nn2,nn3,nn4
     &,nofr,i1,i2,data)
        implicit none
        integer countx,county,countz
        integer unit_no,nn1,nn2,nn3,nn4,nofr
        real data(nn1,nn2,nn3,nn4)
        integer idata(nn4,nn1,nn2), len, i1, i2

c       logical l1,l2

        character*(*) unit_name

C       open(unit_no,file=unit_name,status='old',access='direct',
C       . recl=nn2*nn1*2)
C       inquire(unit_no,exist=l1,opened=l2)
C       read(unit_no,rec=1) idata

        call s_len(unit_name,len)

        call read_binary_field(idata,i1,i2,nn1*nn2*nn4,unit_name,len)

        do county=1,nn2
        do countx=1,nn1
        do countz=1,nn4
          if (idata(countz,countx,county).eq.-9999)
     &idata(countz,countx,county)=0

           data(countx,county,nofr,countz)=
     &float(idata(countz,countx,nn2-county+1))

c SG97 initial data (DEM format) starts in the lower-left corner;
c SG97 this format is wrapped around to have upper-left corner as its start.
c
c JS00 some machines do not account for signed integers - shouldnt matter here.

c          if(data(countx,county,nn3,countz).ge.15535.0)
c    &data(countx,county,countz)=data(countx,county,countz)-65535

        enddo
        enddo
        enddo

ccc      close(unit_no)
        return
        end

C +------------------------------------------------------------------+
	SUBROUTINE JCL
	CHARACTER*(*) FILENM,FORMT

C	-------------------------------------------------------
	ENTRY JCLGET(IUNIT,FILENM,FORMT,IPRNT,istatus)
C
C	  This routine access an existing file with the file name of
C	    FILENM and assigns it unit number IUNIT.
C
		IF(IPRNT.EQ.1) THEN
	PRINT*,' Opening input unit ',IUNIT,' file name ',FILENM
	PRINT*,'		format  ',FORMT
	ENDIF

	OPEN(IUNIT,STATUS='OLD',FILE=FILENM,FORM=FORMT,ERR=1)

	istatus=1
	RETURN

 1    istatus = 0
	return

	END
!
!----------------------------------------------------------------------
        FUNCTION IA(CHR,N,ISPVAL)                         
!                                                              
!  PURPOSE: TO CONVERT A N-BYTES CHARACTER (CHR) TO INTEGER IA. 
!        ** THE INTEGER DATA FILE IS SAVED AS A N-BYTE CHARACTER
!           DATA FILE. THIS FUNCTION IS USED TO RECOVER THE    
!           CHARACTER DATA TO THE INTEGER DATA.               
!                                                            
!  N      --- THE NUMBER OF BYTES IN CHR                    
!  ISPVAL --- DEFAULT VALUE FOR THE NEGATIVE INTEGER.      
!                                                       
        CHARACTER*(*) :: CHR                                
        integer  N, II1, II2, JJ, ISN, M, NBIT, MSHFT, IA2, ispval
        INTEGER  BIT_1, BIT_2                            
!                                                    
        BIT_1 = '200'O     ! BINARY '10000000'        
        BIT_2 = '377'O     ! BINARY '11111111'       
        IA    = 0                                   
!                                                
        II1 = ICHAR(CHR(1:1))                     
! .. GET THE SIGN -- ISN=0 POSITIVE, ISN=1 NEGATIVE:
        JJ  = IAND(II1,BIT_1)                        
        ISN = ISHFT(JJ,-7)                          
!                                                
! .. FOR NEGATIVE NUMBER:
!    BECAUSE THE NEGATIVE INTEGERS ARE REPRESENTED BY THE SUPPLEMENTARY
!    BINARY CODE INSIDE MACHINE.
!                              
        IF (ISN.EQ.1) THEN    
          DO M = N+1,4   
           NBIT = (M-1)*8   
           JJ = ISHFT(BIT_2,NBIT)
           IA = IEOR(JJ,IA)     
          END DO                
        ENDIF                   
!                              
!   .. GET THE BYTE FROM CHR: 
        DO M = 1,N          
         II2 = ICHAR(CHR(M:M)) 
         MSHFT = (N-M)*8      
         IA2   = ISHFT(II2,MSHFT)
!   .. THE ABS(INTEGER):          
         IA = IEOR(IA,IA2)     
        END DO                 
!                              
        IF (IA.LT.0) IA = ISPVAL
!                            
        RETURN                
        END
