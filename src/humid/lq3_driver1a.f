
cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis
cdis
cdis 
cdis
cdis
cdis
cdis
cdis
cdis
cdis
      subroutine lq3_driver1a (i4time,ii,jj,kk,mdf,lct,jstatus)

      USE module_sfc_structure

      implicit none

      include 'grid_fname.cmn'

c     parameter variables

      integer ii,jj,kk
      real mdf
      integer lct
      type (lbsi), dimension(ii,jj) :: sfc_data
      real :: pi, d2r, tempz(3,ii,jj), alt, sol_sub_lat, sol_sub_lon
      real :: zenith_angle      !added for testing calculation
      real zenith ! function

     
      
      integer*4
     1     jstatus(3)
      
      
      integer
     1     istatus,
     1     t_istatus,
     1     i4time,
     1     i4timep,
     1     c_istatus,
     1     save_i4time,
     1     ramsi4time
      
      real
     1     ssh2                 !function
c     1        make_ssh !function type

      real sat(ii,jj,kk)        !saturation ssh at each location
      
      
      real*4
     1     data(ii,jj,kk),
     1     tpw(ii,jj), tpw1(ii,jj), tpw2(ii,jj),
     1     tempsh
      
      integer kstart (ii,jj)
      real qs(ii,jj)
      real ps(ii,jj)
      
      real bias_one

      real*4 mask (ii,jj),cg(ii,jj,kk)
      
      real*4 lt1dat(ii,jj,kk)
      
      character
     1     dirlt1*250,dir*250,rhdir*250,dirpw*250,dir3*250,
     1     extlt1*31,ext*50,rhext*50,extpw*50,ext3*50,
     1     varlt1(kk)*3,
     1     lvl_coordlt1(kk)*4,
     1     unitslt1(kk)*10,
     1     commentlt1(kk)*125
      
c     lat lon variables
   
      character*256  directory
c      character*256 grid_fnam_common
      real lat(ii, jj), lon(ii, jj)
      real rspacing_dum
      character*125 comment_2d
      character*10  units_2d
      character*3 var_2d
      integer len_dir
      character*200 fname
      real factor
      
c     rams stuff--------
      character*9
     1     filename,savefilename,ramsfile
      character ramsvar(kk)*3, ramslvlcoord(kk)*4,
     1     ramsunits(kk)*10, ramscomments(kk)*125
      character rams_dir*250, rams_ext*31
      
c     ------------------
      
      real data_in(ii,jj,kk), delta_moisture(kk), avg_moisture(kk)
      real data_start(ii,jj,kk)
      real data_pre_bound (ii,jj,kk)
      real diff_data(ii*jj)
      real ave,adev,sdev,var,skew,curt
      
      
      character*125 commentline
      
      
      integer*4
     1     i,j,k
      
      integer*4 counter
      
      integer*4 lvllm(kk)
      
      real*4 maps_rh(ii,jj,kk)

      character*3 desired_field
   
      real*4 plevel(kk), p_3d(ii,jj,kk)
      integer*4 mlevel(kk)

c     CLOUD variables

      real :: qadjust (ii,jj,kk)

c     SND variables

      real q_snd (ii,jj,kk)
      real weight_snd(ii,jj,kk)

c     gps variables

      real gps_data (ii,jj)
      real gps_w (ii,jj)
      integer istatus_gps
      
c     
c     gvap variables
c     
      real gvap_data(ii,jj)
      real gvap_w (ii,jj)
      real gw1(ii,jj),gww1(ii,jj)
      real gw2(ii,jj),gww2(ii,jj)
      real gw3(ii,jj),gww3(ii,jj)
      real gvap_p (ii,jj)
      integer istatus_gvap
      
      real pressure_of_level    !function call
      
c     namelist data

      integer  raob_switch
      integer  raob_lookback
      real raob_radius
      integer endian  ! 1 = big, 0 = little , big default
      integer goes_switch
      integer cloud_switch
      integer cloud_d
      integer tiros_switch
      integer sounder_switch
      integer sat_skip
      integer gvap_switch
      integer IHOP_flag
      integer time_diff         !time allowed for latency (sec)
      integer sfc_mix
      integer mod_4dda_1
      real    mod_4dda_factor
      real    t_ref
      integer gps_switch
      character*256 path_to_gvap12,path_to_gvap10,path_to_gps
      namelist /moisture_switch_nl/ raob_switch,
     1     raob_lookback, endian,
     1     raob_radius, goes_switch, cloud_switch, cloud_d
     1     ,tiros_switch, sounder_switch, sat_skip
     1     ,gvap_switch, IHOP_flag, time_diff, gps_switch
     1     ,sfc_mix, mod_4dda_1,mod_4dda_factor,
     1     t_ref,path_to_gvap12,path_to_gvap10,path_to_gps
      
      integer len
      character cdomain(ii)
      
      data extpw/'lh1'/
      data ext3/'lh2'/
      data extlt1/'lt1'/
      data ext /'lq3'/
      data rhext /'lh3'/


c     optran 90 coefficient data information

      character*256 sndr_coeff, sndr_trans
      integer sndr_coeff_len, sndr_trans_len
      common /optn_coef/ sndr_coeff, sndr_trans, sndr_coeff_len, 
     1     sndr_trans_len


c----------------------code   ------------------
c     initialization
c

c     define PI
      pi = acos(-1.0)
      d2r = pi/180.
      write (6,*) 'Starting run for I4TIME', i4time
c     initialize IR emissivity (1:n) and reflectance (1:n)

      do i = 1,ii
         do j = 1,jj
            sfc_data(i,j)%sfc_emiss = 0.98
            sfc_data(i,j)%sfc_refl
     1           = (1.0-sfc_data(i,j)%sfc_emiss)/pi ! isotropic
         enddo
      enddo
c     initialize laps field
      
c     call get_laps congif to fill common block used in pressure assignment
c     routine
      
      write (6,*) ' '
      write (6,*) 'Release 3.10.8 successfully incorporates'
      write (6,*) '1) Fix to func_o.f in use of radiance data'
      write (6,*) '2) OPTRAN90 for BOTH imager and sounder data'
      write (6,*) '3) BOTH  Big and little endian formats'
      write (6,*) '4) Initiated removal of old optran code'
      write (6,*) '5) Removal of common block from .f90 modules'
      write (6,*) '   For HP applications'

      call get_directory(extpw,dirpw,len)
      call get_directory(ext3,dir3,len)
      call get_directory(extlt1,dirlt1,len)
      call get_directory(ext,dir,len)
      call get_directory(rhext,rhdir,len)

      call get_laps_config(grid_fnam_common,istatus)

      if(istatus .ne. 1)then
         write(6,*)' error in get_laps_config'
         return
      endif

c     
c     set namelist parameters to defaults 
      cloud_switch = 1
      cloud_d = 1
      raob_switch = 0
      raob_lookback = 0
      endian = 1 ! big endian is default
      raob_radius = 45000.0  ! meters (45km)
      goes_switch = 0
      sounder_switch = 0
      tiros_switch = 0
      sat_skip = 0
      gvap_switch = 1
      IHOP_flag = 0 
      time_diff = 0
      gps_switch = 1
      sfc_mix = 0
      mod_4dda_1 = 0
      mod_4dda_factor = 0.02
      t_ref = -132.0
      path_to_gvap12 = ' '
      path_to_gvap10 = ' '
      path_to_gps = ' '
      
      call get_directory('static',fname,len)
      open (23, file=fname(1:len)//'moisture_switch.nl',
     1     status = 'old', err = 24)
      
      read(23,moisture_switch_nl,end=24)
      
      
      close (23)
      
      
      if (cloud_switch.eq.0) then
         write(6,*) 'Cloud switch off, ignore clouds'
         write(6,*) 'If available, clouds will be used in GOES adjust'
      else
         write (6,*) 'Clouds will be used in the analysis'
      endif

      if (cloud_d.eq.0) then
         write(6,*) 'analysis will be produced even if clouds are'
         write(6,*) 'not available'
      else
         write(6,*) 'analysis depends on cloud presence'
         write(6,*) 'no cloud field.... no moisture field'
      endif
      
      
      if (raob_switch.eq.0) then
         write(6,*) 'raob switch off, ignoring raobs (.snd files)'
      else
         write (6,*) 'Raob switch on... will use raobs if present'
         write (6,*) 'RAOB radius of influ set to..',raob_radius,' m'
      endif
      
      write(6,*) 'RAOB look back set to ', raob_lookback, 'seconds'
      
      if (goes_switch.eq.0) then
         write(6,*) 'GOES switch off, ignoring goes data'
      else
         write(6,*) 'GOES switch on, attempting GOES ', goes_switch
      endif

      if (sounder_switch.eq.0) then
         write(6,*) 'Sounder switch off'
         write(6,*) 'Using IMAGER data only'
         
         if (goes_switch.eq.12 .or. goes_switch.eq.0) then
            
            if(endian.eq.1) then
               
               sndr_coeff = 'imgr_g12_spectral_coefficients_big'
               sndr_trans = 'imgr_g12_transmittance_coefficients_big'
               
            elseif (endian .eq. 0) then

               sndr_coeff = 'imgr_g12_spectral_coefficients_little'
               sndr_trans = 'imgr_g12_transmittance_coefficients_little' 

            endif
         
            
         elseif (goes_switch.eq.11) then

            if(endian.eq.1) then
               
               sndr_coeff = 'imgr_g11_spectral_coefficients_big'
               sndr_trans = 'imgr_g11_transmittance_coefficients_big'

            elseif (endian .eq. 0) then

               sndr_coeff = 'imgr_g11_spectral_coefficients_little'
               sndr_trans = 'imgr_g11_transmittance_coefficients_little'

            endif
            
         elseif (goes_switch.eq.10) then

            if(endian.eq.1) then
               
               sndr_coeff = 'imgr_g10_spectral_coefficients_big'
               sndr_trans = 'imgr_g10_transmittance_coefficients_big'

            elseif (endian .eq. 0) then


               sndr_coeff = 'imgr_g10_spectral_coefficients_little'
               sndr_trans = 'imgr_g10_transmittance_coefficients_little'

            endif
            
         elseif (goes_switch.eq.9)  then
            if(endian.eq.1) then
               
               sndr_coeff = 'imgr_g09_spectral_coefficients_big'
               sndr_trans = 'imgr_g09_transmittance_coefficients_big'

            elseif (endian .eq. 0) then

               sndr_coeff = 'imgr_g09_spectral_coefficients_little'
               sndr_trans = 'imgr_g09_transmittance_coefficients_little'

            endif
            
         else
            write(6,*) 'invalid sat number for this version'
            write(6,*) 'turning goes function off in var. method'

            goes_switch = 0 ! turning satellite off
         endif

         sndr_coeff_len = len_trim (sndr_coeff)
         sndr_trans_len = len_trim (sndr_trans)

      else
         write(6,*) 'Sounder ON using Sounder data'
         write (6,*) 'setting sounder to goes:', goes_switch

         if (goes_switch.eq.12 .or. goes_switch.eq.0) then

            if(endian.eq.1) then
               
               sndr_coeff = 'sndr_g12_spectral_coefficients_big'
               sndr_trans = 'sndr_g12_transmittance_coefficients_big'
               
            elseif (endian .eq. 0) then

               sndr_coeff = 'sndr_g12_spectral_coefficients_little'
               sndr_trans = 'sndr_g12_transmittance_coefficients_little' 

            endif
         
            
         elseif (goes_switch.eq.11) then

            if(endian.eq.1) then
               
               sndr_coeff = 'sndr_g11_spectral_coefficients_big'
               sndr_trans = 'sndr_g11_transmittance_coefficients_big'

            elseif (endian .eq. 0) then

               sndr_coeff = 'sndr_g11_spectral_coefficients_little'
               sndr_trans = 'sndr_g11_transmittance_coefficients_little'

            endif
            
         elseif (goes_switch.eq.10) then

            if(endian.eq.1) then
               
               sndr_coeff = 'sndr_g10_spectral_coefficients_big'
               sndr_trans = 'sndr_g10_transmittance_coefficients_big'

            elseif (endian .eq. 0) then


               sndr_coeff = 'sndr_g10_spectral_coefficients_little'
               sndr_trans = 'sndr_g10_transmittance_coefficients_little'

            endif
            
         elseif (goes_switch.eq.9)  then
            if(endian.eq.1) then
               
               sndr_coeff = 'sndr_g09_spectral_coefficients_big'
               sndr_trans = 'sndr_g09_transmittance_coefficients_big'

            elseif (endian .eq. 0) then

               sndr_coeff = 'sndr_g09_spectral_coefficients_little'
               sndr_trans = 'sndr_g09_transmittance_coefficients_little'

            endif
            
         else
            write(6,*) 'invalid sat number for this version'
            write(6,*) 'turning goes function off in var. method'

            goes_switch = 0 ! turning satellite off
         endif

         sndr_coeff_len = len_trim (sndr_coeff)
         sndr_trans_len = len_trim (sndr_trans)

      endif
      
      if (tiros_switch.eq.0) then
         write(6,*) 'tiros switch off'
      else 
         write (6,*) 'Attempting to run tiros data from NOAA-',  
     1        tiros_switch
      endif
      
      if (tiros_switch.ne.0 .and. goes_switch.ne.0) then
         write(6,*) 'USING BOTH TIROS AND GOES DATA IN THIS RUN'
      endif
      
      if (sat_skip .eq. 1) then
         write(6,*) 'Use full resolution satellite'
      else
         write(6,*) 'Using partial satellite resolution ',sat_skip
      endif
      
      if (gvap_switch .eq. 1) then
         write(6,*) 'Using goes derived pw, assume data connection'
      else
         write(6,*) 'GVAP not used... nominal state'
      endif

      if (IHOP_flag .eq. 1 ) then
         write (6,*) 'IHOP data used over normal NESDIS GVAP data'
      else
         write(6,*) 'Normal NESDIS GVAP data used'
      endif

      if (gps_switch .eq. 1) then
         write(6,*) 'Using GPS IPW data'
      else
         write(6,*) 'GPS data not used'
      endif
      
      if (time_diff .ne. 0) then
         write(6,*) 'GVAP latency assigned to ',time_diff, 'seconds'
      else
         write(6,*) 'NO latency assinged to GVAP data'
      endif
      
      if (sfc_mix .eq. 1) then
         write(6,*) 'Mixing moisture from sfc'
      else
         write(6,*) 'Sfc moisture field ignored'
      endif
      
      if (mod_4dda_1 .eq.1) then
         write(6,*) 'Mod 4dda active, modifying moisture on output'
         write(6,*) 'Mod 4dda factor is set to, ',mod_4dda_factor
      else
         write(6,*) 'Mod 4dda turned off ... nominal state'
      endif
      
      write(6,*) 'T_ref is set to: ',t_ref
      
      if (path_to_gvap12 .eq. ' '.and. path_to_gvap10 .eq. ' ')then
         write(6,*) 'Path to gvap not assigned, assigning gvap switch 0'
         gvap_switch = 0
      else
         write(6,*) 'Gvap switch assigned, using assigned switch'
         write(6,*) 'Path is ', path_to_gvap12, ' ',path_to_gvap10
         write(6,*) 'GVAP switch is set to ',gvap_switch
      endif

      write(6,*) 'GPS path is ', path_to_gps
     

c     initialize field to lq3 internal missing data flag.
c     initialize total pw to laps missing data flag

      data = -1.e+30
      tpw = mdf
      gvap_w = 0.0
      gvap_data = 0.0
      gps_data = 0.0
      gps_w = 0.0
      mask = 0
      
      jstatus = 0               ! %loc(rtsys_abort_prod)

      call get_pres_3d (i4time,ii,jj,kk,p_3d,istatus)

      if(istatus.ne.1) then
         write (6,*) 'error in getting 3d pressures'
         write (6,*) 'aborting'
         istatus = 0
         return
      endif

c     dependence here now remains to one dimension

      p_3d = p_3d * 0.01 ! convert 3d array to hPa

      lvllm (1:kk)  = int( p_3d (1,1,1:kk))

      plevel (1:kk) = float (lvllm(1:kk))
      mlevel = plevel

c     mark the maps gridpoints -- archaic code ???

      do j = 1,jj,(jj-1)/(jj-1)
         do  i = 1,ii,(ii-1)/(ii-1)

            mask (i,j) = 1

         enddo
      enddo

c     translate the i4time into filename to be used for this run
c     and store

      call make_fnam_lp (i4time,filename,istatus)

      savefilename = filename

      write(6,*) 'FILENAME = ',filename

c     preserve the i4time

      save_i4time = i4time

c     Get background field
      
      call get_modelfg_3d(i4time,'sh ',ii,jj,kk
     1     ,data,istatus)
      
      if (istatus.ne.1) then 
         write (6,*) 'getting background field failed... abort'
         return
      endif

      call check_nan3 (data,ii,jj,kk,istatus)
      if (istatus.ne.1) then
         write(6,*) 'NaN detected from RUC/MAPS...abort'
         return
      endif
      
      i4time = save_i4time
      filename = savefilename

c     check for negative input, assign zero value

      where (data < 0.0) 
         data = 0.0
      end where
      
c     **** obtain lat lons for domain
      
      
c      grid_fnam_common = 'nest7grid' ! used in get_directory to modify
                                ! extension based on the grid domain
      ext = 'nest7grid'
      
c     get the location of the static grid directory
      call get_directory(ext,directory,len_dir)
      
      var_2d='lat'
      call rd_laps_static (directory,ext,ii,jj,1,var_2d,
     1     units_2d,comment_2d,
     1     lat,rspacing_dum,istatus)
      if(istatus .ne. 1)then
         write(6,*)' error reading laps static-lat'
         return
      endif

      call check_nan2 (lat,ii,jj,istatus)
      if (istatus.ne.1) then
         write(6,*) 'NaNs in lat file  abort'
         return
      endif
      
      var_2d='lon'
      call rd_laps_static (directory,ext,ii,jj,1,var_2d,
     1     units_2d,comment_2d,
     1     lon,rspacing_dum,istatus)
      if(istatus .ne. 1)then
         write(6,*)' error reading laps static-lon'
         return
      endif
      
      call check_nan2 (lon,ii,jj,istatus)
      if (istatus.ne.1) then
         write(6,*) 'NaNs in lon file  abort'
         return
      endif      
      
      
c     open file for laps temp data
      varlt1 = 't3 '

      call read_laps (save_i4time,save_i4time,
     1     dirlt1,
     1     extlt1,
     1     ii,jj,kk,kk,
     1     varlt1,
     1     lvllm,
     1     lvl_coordlt1,
     1     unitslt1,
     1     commentlt1,
     1     lt1dat,
     1     t_istatus)
      if (t_istatus.ne.1) then
         print*, 'no lt1 quality control performed...'
         print*, 'missing 3-d temp data'
         write(6,*) 'ABORTING MOISTURE RUN...!!!'
         istatus = 0            ! failure
         return
      endif
      
c     fill new data structure surface data

      sfc_data%lat = lat
      sfc_data%lon = lon

      do i = 1, ii
         do j = 1,jj

c     compute the secant of zenith angle for each goes satellite
c     (1) is goes east
c     (2) is goes west
c     (3) is goes 9
c     (x) add additional satellites as needed
            tempz(1,i,j) =  zenith(lat(i,j)*d2r,lon(i,j)*d2r,
     1           0.0,-75.*d2r)
            tempz(2,i,j) = zenith(lat(i,j)*d2r,lon(i,j)*d2r,
     1           0.0,-135.*d2r)
            tempz(3,i,j) = zenith(lat(i,j)*d2r,lon(i,j)*d2r,
     1           0.0,+155.*d2r)

            do k = 1,3
               tempz(k,i,j) = 1./cos(tempz(k,i,j)*d2r)
               sfc_data(i,j)%secza(k) = tempz(k,i,j)
            enddo               !k

         enddo
      enddo

      call solar_position (lat(1,1), lon(1,1), i4time, alt, 
     1     sol_sub_lat, sol_sub_lon)

c     note that the solar hour angle is now used to compute the solar subpoint
c     longitude

      sol_sub_lon = lon(1,1) - sol_sub_lon

c     now the solar subpoint latitude is given by the declination
c     the solar subpoint longitude is given by adjusting the longitude with
c     the solar hour angle offset.

c     the solar zenith angle is now computable from the following.

      do j = 1,jj
         do i = 1,ii

            zenith_angle = zenith ( 
     1           lat(i,j)*d2r,lon(i,j)*d2r, sol_sub_lat*d2r, 
     1           sol_sub_lon*d2r  )
            if(abs (zenith_angle) >= 85.0) then
               sfc_data(i,j)%secsola = 11.47 ! largest allowed value
c     sun is assumed below horizon
            else
               sfc_data(i,j)%secsola = 1./ cos ( d2r*zenith_angle)
            endif

         enddo
      enddo

      call check_nan3 (lt1dat,ii,jj,kk,istatus)
      if (istatus.ne.1) then
         write(6,*) 'NaN detected from lt1...ABORT'
         return
      endif
      
      
c     generate saturation array for the domain based on current temp
c     ans pressure

      do k = 1,kk
         do j = 1,jj
            do i = 1,ii
               sat(i,j,k)=ssh2( p_3d(i,j,k) ,lt1dat(i,j,k)-273.15,
     1              lt1dat(i,j,k)-273.15, t_ref )/1000.
            enddo               ! i
         enddo                  ! j
      enddo                     ! k
      
      
c     perform initialquality control check for supersaturation after ingest

      write(6,*)  'perform qc for supersaturation'
      counter = 0
      do k = 1,kk
         write (6,*)
         write (6,*) 'Level ', k, '   ', p_3d(1,1,k)
         write (6,*)
         
         do j = jj,1,-1
            do i = 1,ii
               if ( data(i,j,k)/sat(i,j,k) .ge. 1.0) then
                  cdomain(i) = 'x'
                  if(data(i,j,k)/sat(i,j,k) .gt. 1.01) then
                     cdomain(i) = 's'
                  endif
                  counter = counter + 1
                  data(i,j,k) = sat(i,j,k)
               else
                  write (cdomain(i),34) int(data(i,j,k)/sat(i,j,k)*10.)
 34               format (i1)
               endif
            enddo
            write(6,*) (cdomain(i),i=1,ii)
         enddo
      enddo
      
      if(counter.gt.0) then
         
         write(6,*) ' '
         write(6,*) 'Questionable INPUT DATA DETECTED'
         write(6,*)  counter,' times.'
         write(6,*) ' '
      endif

c     initial check for computed nans
      call check_nan3(data,ii,jj,kk,istatus)
      if(istatus.ne.1) then
         write(6,*) 'initial data corrupt after checking for'
         write(6,*) 'supsaturation... hosed before beginning'
         write(6,*) 'var:data  routine:lq3driver1a.f'
         return
      endif

      
      
c     record total moisture

      data_in = data
      data_start = data
      
c     ****  execute raob step if switch is on
      
      if(raob_switch.eq.1) then
         write (6,*) 'begin raob insertion'

         call snd_step (i4time,p_3d,raob_lookback, lat,lon,
     1        lt1dat, ii,jj,kk, q_snd,weight_snd, 
     1        raob_radius, raob_switch)
             
      else
         write(6,*) 'the raob switch is off... SND skipped'
      endif
      
c     ****   laps cloud data. used for cloud, bl, goes
      
      call mak_cld_grid (i4time,i4timep,cg,ii,jj,kk,
     1     lct,c_istatus)

      if (cloud_d.eq.1) then
         if(c_istatus .ne. 1) then
            write(6,*) 'cloud data not available'
            write(6,*) 'terminating'
            istatus = 0
            return
         endif
      endif

      c_istatus = 0
      if (i4time.eq.i4timep) c_istatus = 1

      if(cloud_d.eq.1 .and. c_istatus.eq.0) then
         write(6,*) 'cloud data not available for exact time'
         write(6,*) 'cloud_dependence switch is on'
         write(6,*) 'aborting'
         istatus = 0
         return
      endif

      call check_nan3 (cg,ii,jj,kk,istatus)
      if (istatus.ne.1) then
         write(6,*) 'NaN detected from Cloud Grid...ABORT'
         return
      endif

c     ***   insert bl moisture

      data_pre_bound = data

      print*, 'calling lsin'
c     insert boundary layer data
      call lsin (i4time,p_3d,sfc_data,lt1dat,data,cg,tpw,bias_one,
     1     kstart,qs,ps,mdf,ii,jj,kk,istatus)

c     check for supersaturation

      where (data >= sat .and. data /= mdf )
         data = 0.99*sat
      end where

c     check fields after lsin call
      call check_nan3(data,ii,jj,kk,istatus)
      if(istatus.ne.1) then
         write (6,*) 'Nan generated in lsin'
         write (6,*) 'var:data  routine:lq3driver1a.f'
         return
      endif
      call check_nan2(tpw,ii,jj,istatus)
      if(istatus.ne.1) then 
         write(6,*) 'Nan generated in lsin'
         write(6,*) 'var:tpw   routine:lq3driver1a.f'
         return
      endif
      
      write(6,*) 'finished with routine lsin'
      if( sfc_mix.eq.1)then
         write(6,*) 'Lsin allowed to modify data field'
         
         write(6,*) 'Reporting incremental boundary layer effects'
         
         call report_change (data_in, data, p_3d,mdf,ii,jj,kk)
         
         write(6,*) 'Reporting net change'
         call report_change (data_start, data, p_3d, mdf, ii,jj,kk)

         data_in = data

c     end report moisture change block
         
      else
         write(6,*) 'Lsin and sfc mixing step skipped'

         where (data >= 0.0 )
            data = data_pre_bound
         end where
         
      endif

c     call to get cloud adjust parameter

      do i = 1,ii
         do j = 1,jj
c     evaluate the probability of cloud
            qadjust(i,j,kk) = 0.0
            do k = 1,kk-1
               qadjust(i,j,k) = 0.0
               if (ps(i,j) > p_3d(i,j,k)) then
                                !evaluate the probabilty of cloud 
                  call test_cloud (lt1dat(i,j,k),lt1dat(i,j,k+1),
     1              data(i,j,k),data(i,j,k+1),
     1              p_3d(i,j,k), sat(i,j,k),
     1                 qadjust(i,j,k))
               endif
            enddo
         enddo
      enddo
c     gps data inserstion step (bias correction to gvap only)

      istatus_gps = 0

      if (gps_switch .eq. 1) then

         write(6,*) 'Initiate bias correction of gps data'

         call process_gps (ii,jj,gps_data,gps_w,
     1        tpw,lat,lon,time_diff,
     1        path_to_gps,filename,istatus_gps)
         
c     gvap data insertion step
         
      endif

c     gvap data acquisition

      istatus_gvap = 0
      
      if (gvap_switch.eq.1) then
         
         write(6,*) 
         write(6,*) 'Begin GVAP insertion setep'
         write(6,*) 
         
         call process_gvap(ii,jj,gvap_data,gvap_w,
     1        gw1,gw2,gw3,gww1,gww2,gww3,gvap_p,mdf,
     1        lat,lon,time_diff,IHOP_flag,
     1        path_to_gvap12,path_to_gvap10,filename,istatus_gvap)
         

         if(istatus_gvap.eq.1 .and. istatus_gps.eq.1) then ! correct gvap
            continue            ! placeholder for correction call
            write(6,*) 'gvap/gps correction currently disabled'
         endif  
         
      else
         write(6,*) 'Gvap off, not using gvap or attempting' 
         write(6,*) 'any adjustment'
      endif
     
         
c     CHECKING PROCESS OUTPUT

      call check_nan2(gvap_w,ii,jj,istatus)
      if(istatus.ne.1) then
         write(6,*) 'Failed in nan after adjust'
         write(6,*) 'var:gvap_w  routine lq3driver'
         return
      endif
         
      call check_nan2(gps_w,ii,jj,istatus)
      if(istatus.ne.1) then
         write(6,*) 'Failed in nan after adjust'
         write(6,*) 'var:gps_w  routine lq3driver'
         return
      endif
         
      call check_nan2(gps_data,ii,jj,istatus)
      if(istatus.ne.1) then
         write(6,*) 'Failed in nan after adjust'
         write(6,*) 'var:gps_data  routine lq3driver'
         return
      endif
         
      call check_nan2(gvap_data,ii,jj,istatus)
      if(istatus.ne.1) then
         write(6,*) 'Failed in nan after adjust'
         write(6,*) 'var:gsp_data  routine lq3driver'
         return
      endif
         
      call check_nan3(data,ii,jj,kk,istatus)
      if(istatus.ne.1)then
         write(6,*) 'Nan report from TPW processing'
         write(6,*) 'var:data  routine lq3driver1a.f'
         return
      endif
      
      write(6,*) 'Reporting changes from TPW data types'
      
      call report_change (data_in, data, p_3d,mdf,ii,jj,kk)
      
      write(6,*) 'Reporting net change'
      call report_change (data_start, data, p_3d, mdf, ii,jj,kk)

      data_in = data

      write(6,*) 

c     make call to TIROS moisture insertion


      if (tiros_switch .ne. 0) then

         if(c_istatus.eq.1 .and. t_istatus.eq.1) then

        
            write (6,*) 'begin TIROS insertion step'
      
            call tiros (
     1           data,          ! specific humidity g/g
     1           lat,lon,       ! lat and longitude (deg)
     1           i4time,        ! i4time of run (seconds)
     1           p_3d,          ! pressure hpa (laps vert grid)
     1           cg,            ! cloud array
     1           lt1dat,        ! lt1 (laps 3d temps)
     1           14,            ! satellite number
     1           ii,jj,kk       ! grid dimensions
     1           )

         else
            
            write(6,*)
            write(6,*)
            write(6,*)
            write(6,*) 'tiros moisture insertion step skipped'
            write(6,*) 'cloud or lt1 data not current'
            write(6,*) 'cannot assume clear conditions or'
            write(6,*) 'use alternate lt1.. this will create'
            write(6,*) 'forward model problems....'
            write(6,*) 'tiros moisture insertion step skipped'
            write(6,*)
            write(6,*)
            write(6,*)

         endif

      else

         write(6,*) 'tiros switch is off... tiros step skipped...'
         
      endif


      call int_tpw (data,kstart,qs,ps,p_3d,tpw1,mdf,ii,jj,kk)


c     make call to variational solution, radiance no longer a prerequisite


      if(c_istatus.eq.1 .and. t_istatus.eq.1) then
         
         write (6,*) 'Begin variational assimilation'
         call variational (
     1        data,             ! 3-d specific humidity g/g
     1        sfc_data,         ! struct surface data (type lbsi)
     1        i4time,           ! i4time of run
     1        p_3d,             ! pressure mb
     1        cg,               ! 3-e cloud field 0-1 (1=cloudy)
     1        c_istatus,        ! cloud istatus
     1        sat,              ! saturated field
     1        qadjust,          ! q increase needed for cloud formation
     1        lt1dat,           ! laps lt1 (3-d temps)
     1        mdf,
     1        qs,kstart,
     1        goes_switch,      ! goes switch and satellite number
     1        sounder_switch,   ! sounder switch, 0=imager,1=sndr
     1        sat_skip,         ! normally 1 for full resolution
     1        gw1,gw2,gw3,
     1        gww1,gww2,gww3,
     1        gvap_p,istatus_gvap,
     1        gps_data,
     1        gps_w,
     1        istatus_gps,
     1        q_snd,
     1        weight_snd,
     1        raob_switch,
     1        ii,jj,kk
     1        )
         
         write (6,*) 'GOES step complete, effects logged.'
         
         call report_change (data_in, data, p_3d,mdf,ii,jj,kk)
         
         write(6,*) 'Reporting net change'
         call report_change (data_start, data, p_3d, mdf, ii,jj,kk)
         
         data_in = data
         
c     end report moisture change block
         
         
      else
         
         write(6,*)
         write(6,*)
         write(6,*)
         write(6,*) 'variational step skipped'
         write(6,*) 'cloud or lt1 data not current'
         write(6,*) 'cannot assume clear conditions or'
         write(6,*) 'use alternate lt1.. this will create'
         write(6,*) 'forward model problems....'
         write(6,*) 'thus basically returning first guess'
         write(6,*)
         write(6,*)
         write(6,*)
         
      endif
      
c     assess impact of system compared to GPS sites
      
      call int_tpw (data,kstart,qs,ps,p_3d,tpw2,mdf,ii,jj,kk)
      
      call impact_assess (data_start, data_in, tpw1,tpw2,
     1     ii,jj,kk,
     1     gps_data, gps_w, p_3d, mdf)
      
      
c     *** insert cloud moisture, this section now controled by a switch

      if(cloud_switch.eq.0) then
         write(6,*) ' '
         write(6,*) 'Cloud switch  '
         write(6,*) 'Skipping cloud moistening here'
         write(6,*) ' '
         write(6,*) ' '
         c_istatus = 0          !force this to skip here
      endif
      
      if(c_istatus.ne.1 .or. i4time.ne.i4timep)then
         c_istatus = 0
         write(6,*) 'Cloud field not available for exact time'
         write(6,*) 'assume ALL CLEAR values used'
      else                      ! increase moisture based on cloud amount
         
         write(6,*) 'Saturate in cloudy areas'

         call check_nan3(cg,ii,jj,kk,istatus)
         if(istatus.ne.1)then
            write(6,*) 'Nan in cg data prior to saturation'
            return
         endif
         call check_nan3(lt1dat,ii,jj,kk,istatus)
         if(istatus.ne.1) then
            write(6,*) 'NaN in lt1dat prior to saturation'
            return
         endif

         call check_nan3(data,ii,jj,kk,istatus)
         if(istatus.ne.1)then
            write(6,*) 'Nan in data prior to saturation'
            return
         endif

c     saturate in cloudy areas.  

         do k = 1,kk
            write(6,*) lvllm(k),'checking lvllm prior to sat'
            do j = 1,jj
               do i = 1,ii

                  call cloud_sat (cg(i,j,k),qadjust(i,j,k),data(i,j,k))

               enddo
            enddo
         enddo

         call check_nan3(data,ii,jj,kk,istatus)
         if(istatus.ne.1)then
            write(6,*) 'Nan in data AFTER call to saturation'
            return
         endif

         
         write (6,*) 'Reporting cloud effects on analysis'
         call report_change (data_in, data, p_3d,mdf,ii,jj,kk)
         
         write(6,*) 'Reporting net change'
         call report_change (data_start, data, p_3d, mdf, ii,jj,kk)
         
         data_in = data
         
      endif
      
  
c     mod_4dda_1 to decrease overall water in 4dda mode running at AFWA
      
      if(mod_4dda_1 .eq. 1) then ! act to decrease overall water
         
         do k=1,kk
            factor = 1.-(float(k)*mod_4dda_factor)
            do j=1,jj
               do i=1,ii
                  data(i,j,k) = data(i,j,k)*factor
               enddo
            enddo
         enddo
         
         write(6,*) ' mod_4dda loop complete'
         
         call report_change (data_in, data, p_3d,mdf,ii,jj,kk)
         
      endif
      
c     repeat quality control check for supersaturation after pre-analysis
      write (6,*)  'perform qc for supersaturation'
      counter = 0
      do k = 1,kk
         write(6,*) 'Level ',k, '    ', p_3d(1,1,k)
         do j = jj,1,-1
            do i = 1,ii
               
c               tempsh = ssh2( float(lvllm(k)) ,lt1dat(i,j,k)-273.15,
c     1              lt1dat(i,j,k)-273.15,t_ref )/1000.
               
               if ( data(i,j,k)/sat(i,j,k) .ge. 1.0) then
                  cdomain(i) = 'x'
                  if(data(i,j,k)/sat(i,j,k).gt.1.01) cdomain(i:i) = 's'
                  counter = counter + 1
                  data(i,j,k) = sat(i,j,k)
               elseif (data(i,j,k) .lt. 0.0) then
                  cdomain(i) = 'M'
                  
               else
                  write (cdomain(i),35) int(data(i,j,k)/sat(i,j,k)*10.)
 35               format (i1)
                  
               endif

            enddo
            write(6,*) (cdomain(i),i=1,ii)
         enddo
      enddo
      
      if (counter.ne.0) then
         write (6,*) 'supersaturation has been corrected, 
     1        ',counter,' times.'
      endif
     
c     recompute tpw including clouds and supersat corrections
      
      call int_tpw(data,kstart,qs,ps,p_3d,tpw,mdf,ii,jj,kk)
      
c     place the accepted missing data flag in output field
c     sum over the entire grid for a total water sum value for 
c     QC study.
      
c     conform to laps missing data flag (switch from mine)

      where (data < 0.0)
         data = mdf
      endwhere

c     sum total water for simple print diagnostic

      tempsh = 0.0
     
      do i = 1,ii
         do j = 1,jj
            do k = 1,kk
               
               if(data(i,j,k) /= mdf) then
                  tempsh = tempsh + data(i,j,k) ! sum if good data
               endif
               
            enddo
         enddo
      enddo

c     log the amount of water vapor
      
      write (6,*) ' '
      write (6,*) ' '
      write (6,*) '***************************** '
      write (6,*) 'Average water in volume (g/g)*10000'
      write (6,*) tempsh/float(ii)/float(jj)/float(kk)*10000.
      write (6,*) '***************************** '
      write (6,*) ' '
      write (6,*) ' '
      
c     check for NaN values and Abort if found
      
      call check_nan3(data,ii,jj,kk,istatus)
      if(istatus.ne.1) then
         write(6,*) 'NaN values detected (sh array)... aborting'
         return
      endif
      
      call check_nan2(tpw,ii,jj,istatus)
      if(istatus.ne.1) then
         write(6,*) 'NaN values detected (tpw array)... aborting'
         return
      endif
      
      
c     write final 3-d sh field to disk
      commentline = 'maps with clouds and surface effects only'
      call writefile (save_i4time,commentline,mlevel,data,
     1     ii,jj,kk,istatus)
      if(istatus.eq.1)        jstatus(1) = 1
      
c     write total precipitable water field
      call write_lh4 (save_i4time,tpw,bias_one,ii,jj,istatus)
      if(istatus.eq.1) jstatus(3) = 1
      
      
c     generate lh3 file (RH true, RH liquid)
      if (t_istatus.eq.1) then
         call lh3_compress(data,lt1dat,save_i4time,lvllm,t_ref,
     1        ii,jj,kk,istatus)
         if(istatus.eq.1)        jstatus(2) = 1
      else
         print*, 'no lh3 or rh data produced...'
         print*, 'no laps 3-d temp data avail'
         jstatus(2) = 0
      endif
      
      write (6,*) 'Reporting overall changes to moisture'
      
      call report_change (data_start, data, p_3d,mdf,ii,jj,kk)
      
      return
      
 24   write(6,*) 'error finding moisture switch file'
      write(6,*) 'check to see it is under'
      write(6,*) fname(1:len)//'moisture_switch.nl'
      write(6,*) 'aborting'
      
      return
      
      end
