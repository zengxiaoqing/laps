

		subroutine start_up_mix(ni,nj,niveau)
		
                use lfmgrid, ONLY: hzsig,htsig,hpsig,hmrsig,husig,hvsig,hwsig,llat,llon

                w_to_omegaf(w,pres,scaleht) = - (w * pressure_pa) / scale_height 

		!!-------------------------------------------------
		!! Defined structure needed by read_config_file()
		!!-------------------------------------------------
		TYPE myStruct
      INTEGER :: ip2;                 	 	!! ip2
      INTEGER :: npas; 							!! npas 
      INTEGER :: yyymmddhh;					!! yyyymmddhh
      INTEGER :: hh_init;						!! hh_init
      INTEGER :: hh_fcst;						!! hh_prevu
      INTEGER :: nbr_level_2keep;			!! nbr_niveau_voulu
      INTEGER :: nbr_level_avail;			!! nbr_niveau_disponible
      REAL :: resolution;						!! resolution
      REAL :: thres_meso;						!! seuil_meso
      REAL :: thres_tornade;					!! seuil_tornade
      REAL :: thres_severe_updraft(0:1) 	!! seuil_severe_updraft
      REAL :: thres_couplet_tourbillon(0:1,0:1) !! seuil_couplet_tourbillon
      CHARACTER*10  :: region					!! region
      CHARACTER*200 :: ext_file_out;		!! ext_fichier_out
      CHARACTER*200 :: ext_file_in;			!! ext_fichier_in
      CHARACTER*200 :: path_file_in;		!! path_fichier_in
      CHARACTER*200 :: path_file_out;		!! path_fichier_out
      CHARACTER*100 :: varname3D;			!! nom_var3D
      CHARACTER*50  :: varname2D;			!! nom_var2D
		END TYPE myStruct

                real dx(ni,nj),dy(ni,nj),vort_2d(ni,nj),div_2d(ni,nj)

		!!----------------------------------------------------
		!! Defined const needed throughout the program 
		!!----------------------------------------------------
		INTEGER, PARAMETER :: DEBUG = 1  		!! DEBUG == 0 to not print any of the steps / DEBUG ==1 to print every steps 
		
		!!-------------------------------------------------------
		!! Define 3D variable position in the 3D array called
		!! 3D_fields. The 3D_fields array has 3 dimensions 
		!! (3D_fields[Field_Position][Vertical_level][Horizontal_Grid_Point])
		!!-------------------------------------------------------
		INTEGER, PARAMETER :: GZ   = 0  !! Geopotential height (DAM)
		INTEGER, PARAMETER :: TT   = 1  !! Air Temperatutre (C)
		INTEGER, PARAMETER :: HU   = 2  !! Specific Humidity (kg/kg)
		INTEGER, PARAMETER :: UV   = 3  !! Wind speed (knot)
		INTEGER, PARAMETER :: WW   = 4  !! Vertical velocity (Pa/s)
		INTEGER, PARAMETER :: WD   = 5  !! Wind direction (meteorological degree)
		INTEGER, PARAMETER :: QR   = 6  !! Relative Vorticity (s^-1)
		INTEGER, PARAMETER :: ZET  = 7  !! Equivalent Total Reflectivity (dBZ)
		INTEGER, PARAMETER :: QJT1 = 8  !! Graupel Mixing Ratio (kg/kg)
		INTEGER, PARAMETER :: QIT1 = 9  !! Ice Mixing Ratio (kg/kg)
		INTEGER, PARAMETER :: QHT1 = 10 !! Hail Mixing Ratio (kg/kg)
		INTEGER, PARAMETER :: SLW  = 11 !! Supercooled Liquid Water (kg/m^3)
		INTEGER, PARAMETER :: DMH  = 12 !! Mean Hail Diameter (m)
		
		INTEGER, PARAMETER :: nbr_var_3D = 13
		
		!!-------------------------------------------------------
		!! Define 2D variable position in the 2D array called
		!! 2D_fields. The 3D_fields array has 2 dimensions 
		!! (2D_fields[Field_Position][Horizontal_Grid_Point])
		!!-------------------------------------------------------
		INTEGER, PARAMETER :: SWEAT   = 0 !! SWEAT index
		INTEGER, PARAMETER :: LWC     = 1 !! Integrated Liquid Water Content from the Surface upward
		INTEGER, PARAMETER :: LWC_700 = 2 !! Integrated Liquid Water Content from 700mb upward
		
		INTEGER, PARAMETER :: nbr_var_2D = 3
		
		!!------------------------------------------------------------
		!! Define the position of the severe weather variablesin the
		!! array called TV_Var. The TV_Var array has 2 dimensions 
		!! (TV_Var[SW_Variable_Position][Horizontal_Grid_Point])
		!!-----------------------------------------------------------
		INTEGER, PARAMETER :: Hauteur_MaxR_aboveT0         = 0  !! Height of the Maximun Recflectivity above the freezing level 
		INTEGER, PARAMETER :: MaxR_aboveT0                 = 1  !! Maximun Recflectivity above the freezing level 
		INTEGER, PARAMETER :: Couplet_Tourbillon           = 2  !! Vorticity couplets
		INTEGER, PARAMETER :: Cisaillement_3km             = 3  !! 0-3km wind shear
		INTEGER, PARAMETER :: Cisaillement_6km             = 4  !! 0-6km wind shear
		INTEGER, PARAMETER :: Couplet_Positif              = 5  !! Positive vorticity couplet
		INTEGER, PARAMETER :: Couplet_Negatif              = 6  !! Negative vorticity couplet
		INTEGER, PARAMETER :: zoneSurplomb                 = 7  !! Weak echo region
		INTEGER, PARAMETER :: Severe_Updraft               = 8  !! Severe updraft
		INTEGER, PARAMETER :: Vertical_Moisture_Flux       = 9  !! Vertical Moisture Flux
		INTEGER, PARAMETER :: Mesocyclone                  = 10 !! Mesocyclone
		INTEGER, PARAMETER :: Meso_base                    = 11 !! Height of the mesocyclone's base
		INTEGER, PARAMETER :: VEF                          = 12 !! Bounded Weak Echo Region
		INTEGER, PARAMETER :: Tornade                      = 13 !! Tornado
		INTEGER, PARAMETER :: Grosseur_Grele_Max           = 14 !! Maximum hail size in a column
		INTEGER, PARAMETER :: Grele_sfc                    = 15 !! Hail size reaching the surface 
		INTEGER, PARAMETER :: Supercooled_Liquid_Water     = 16 !! Integrated supercooled liquid water
		INTEGER, PARAMETER :: Vertical_Ice_Flux            = 17 !! Vertical ice flux 
		INTEGER, PARAMETER :: Vertical_Moisture_Flux_NEW   = 18 !! vertical moisture flux (new)
		
		INTEGER, PARAMETER :: nbr_variableTV               = 19
		
		!!-------------------------------------------
		!! Defined variables needed in the program
		!!-------------------------------------------
!       	INTEGER, PARAMETER :: nbr_niveau  = 10
!       	INTEGER, PARAMETER :: ni          = 25
!	        INTEGER, PARAMETER :: nj          = 25
		
		INTEGER :: i, j, k, m, ierr
		
		REAL  fields3d(0:(ni*nj)-1,0:nbr_niveau-1,0:nbr_var_3D-1)
		REAL  fields3d_C(0:nbr_var_3D-1,0:nbr_niveau-1,0:(ni*nj)-1)
		REAL  fields2d(0:(ni*nj)-1,0:nbr_var_2D-1)
		REAL  fields2d_C(0:nbr_var_2D-1,0:(ni*nj)-1)
		REAL  TV_Var(0:nbr_variableTV-1,0:(ni*nj)-1)
		REAL  ProbTV(0:(ni*nj)-1)
		
		character*6 FILE_name
		
		TYPE(myStruct):: ParamConfig;                       
		
		FILE_name="TV.cfg"

		!!--------------------------------
		!! Initialize/Populate the array
		!!--------------------------------
                scaleht = 8000.
                call get_grid_spacing_array(lat,lon,ni,nj,dx,dy)

                do k = 1,niveau

                  call get_grid_spacing_array(llat,llon,ni,nj,dx,dy)
                  call vortdiv(husig(:,:,k),hvsig(:,:,k),vort_2d,div_2d,ni,nj,dx,dy) ! QR

                  do j = 1,nj
                  do i = 1,ni
                    ij = (j-1)*ni + (i-1)
		    fields3d(ij,0:k-1,GZ) = hzsig(i,j,k)
		    fields3d(ij,0:k-1,TT) = htsig(i,j,k)
		    fields3d(ij,0:k-1,HU) = hmrsig(i,j,k) ! mixing ratio is used for now

		    speed = sqrt(husig(i,j,k)**2+hvsig(i,j,k)**2)
		    fields3d(ij,0:k-1,UV) = speed
                    fields3d(ij,0:k-1,WW) = w_to_omegaf(hwsig(i,j,k),hwsig(i,j,k),scaleht)
                    if(speed .gt. 0.)then
                      fields3d(ij,0:k-1,WD) = atan3d(-husig(i,j,k),-hvsig(i,j,k))
                    else
                      fields3d(ij,0:k-1,WD) = 0.
                    endif

                    fields3d(ij,0:k-1,QR) = vort_2d(i,j)                                      

                  enddo ! i
                  enddo ! j
                  

                enddo ! k
                fields2d=0

		!!-------------------------------------------
		!! Read the configuration file
		!!-------------------------------------------
		ierr = Read_Config_File(%VAL(DEBUG),FILE_name,ParamConfig)
		
		!!-----------------------------------------
		!! Re-arrange the multi-dimentional fortran
		!! to have the C array format
		!!-----------------------------------------
		fields3d_C = reshape(fields3d, shape(fields3d_C), ORDER = (/3, 2, 1/))
		
		fields2d_C = transpose(fields2d)
		
     	!!-------------------------------------------
		!! Find the severe weather structural
                !! elements and the environmental ingredients
		!!-------------------------------------------
		ierr = potentielle_temps_violent(%VAL(DEBUG), %VAL(nbr_niveau), &
       %VAL(ni), %VAL(nj), ParamConfig, %VAL(nbr_var_3D), &
       fields3d_C, %VAL(nbr_variableTV), TV_Var)

		!!-------------------------------------------
		!! Compute the severe thunderstorm intensity
                !! for each detected cell
		!!-------------------------------------------
		ierr = CalculProbTempsViolent( %VAL(DEBUG), &
       %VAL(ni), %VAL(nj), ParamConfig, %VAL(nbr_var_2D), &
       %VAL(nbr_variableTV), fields2d_C, TV_Var, ProbTV)
		
                return
                end
		
		
