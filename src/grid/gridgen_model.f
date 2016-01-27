cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 

	Program gridgen_model
C**********************************************************************
c Portions of the following code were taken from the RAMS software and 
c   were used by permission of the *ASTER Division/Mission Research 
c   Corporation. 
C**********************************************************************
C*	This program will be used to map model grids based on the RAMS*
C*  version 2b which uses polar stereographic projections. Other      *
C*  projections have since been added.                                *
c*                                                                    *
C*********************************************************************
c Beginning July 2000, this software was updated to accomodate the    *
c WRFSI (Weather Research and Forecasting Standard Initialization)    *
c grid requirements.  This included adding stagger calculations, and  *
c map projection factors (map factor and components of coriolis parm, *
c landuse, soil type, vegetation greenness fraction, mean annual soil *
c temperature, albedo, and sea surface temperature.  Significant work *
c has occurred to reduce memory requirements.  Many sections of code  *
c previously in this file now reside in file gridgen_utils.f          *
c For further information contact                                     *
c Mr. John Smart NOAA/FSL, smart@fsl.noaa.gov (303-497-6590)          *
c*********************************************************************

        use horiz_interp

        integer istag
        integer n_staggers
        integer lf
        character c6_maproj*6
        character c8_maproj*8
        character c10_grid_f*10 ! Type of domain (nest7grid, wrfsi)/
        character c_dataroot*200
        character cstaticdir*200
        logical   localize
        real,     allocatable  ::  mother_lat(:,:)
        real,     allocatable  ::  mother_lon(:,:)
        real,     allocatable  ::  dum2d(:,:)
        real      La1,Lo1,La2,Lo2
        real      LoV,Latin1,Latin2
        real      Dx,Dy

        integer lli_orig
        integer llj_orig
        integer uri_orig
        integer urj_orig

        include 'grid_fname.cmn'
 
        write(6,*)
        write(6,*)' gridgen_model: start'

        call get_config(istatus)
        if(istatus.ne.1)then
           print*,'Error returned: get_config'
           print*,'Abort'
           stop
        endif

        call find_domain_name(c_dataroot,c10_grid_f,istatus)
        if(istatus .ne. 1)then
            write(6,*) 'Error: find domain name - stop'
            stop
        endif
        call s_len(c10_grid_f,lf)

        if(c10_grid_f(1:lf).eq.'wrfsi')then
           call get_num_domains(num_domains,istatus)
           if(istatus.ne.1)then
              print *,'Error returned: get_num_domains'
              print *,'Abort'
              stop
           endif
        else
           num_domains=1
        endif
c -----------------------------------------------------------------
c variable "nest" is common in src/include/grid_fname.cmn
c and is used to communicate the domain # throughout the code.
c -----------------------------------------------------------------
        do nest=1,num_domains
c
c determine if this domain requires localizing or re-localizing
c
           call get_directory('static',cstaticdir,len_dir)

           call rd_static_attr(cstaticdir,nest
     .,c10_grid_f,nx_dom,ny_dom,Dx,Dy,La1,Lo1,Latin1
     .,Latin2,LoV,c8_maproj,istatus)

c -------------------------------------------------------------------------
c routine below accesses environment variable $FORCE_LOCALIZATION which, if
c set, will return variable "localize" = .true.
c -------------------------------------------------------------------------
 
           call eval_localization(cstaticdir,nest,localize
     .,c10_grid_f,La1,Lo1,istatus)

           call get_parent_id(iparent_id,istatus)

           if(localize)then

c get_mother_dims,  uses parent_id as index to determine mother domain specs
           if(nest.gt.1)then
              call get_mother_dims(nx_dom,ny_dom,istatus)
              allocate (mother_lat(nx_dom,ny_dom),
     +                  mother_lon(nx_dom,ny_dom),
     +                  dum2d(nx_dom,ny_dom),stat=istat)
              if(istat.ne.0)then
                 print*,'cannot deallocate mother arrays: ',istat
                 stop
              endif

              call get_mother_domain(nx_dom,ny_dom,mother_lat
     +,mother_lon,istatus)
              if(istatus.lt.1)then
                 print*
                 print *,'Error reading MOTHER lat, lon data.'
                 print *,'Abort'
                 print*
                 stop
              endif
           endif   

           call get_grid_dim_xy(nx_dom,ny_dom,istatus)
           call get_ratio_to_parent(iratio_to_parent,istatus)
           call get_domain_origin(lli_orig,llj_orig,uri_orig,urj_orig
     +,istatus)
c
c use nx_dom and ny_dom generically regardless of a nest or not.
c
           call get_grid_spacing(Dx,istatus)
           call get_c6_maproj(c6_maproj,istatus)
           call downcase(c6_maproj,c6_maproj)

           if(c6_maproj.ne.'rotlat')Dy=Dx

           print*
           print*,'---------------------------------------------'
           print*,'Domain #         = ',nest
           print*,'Parent of nest   = ',iparent_id
           print*,'Ratio to parent  = ',iratio_to_parent
           print*,'Grid_spacing(m)  = ',Dx
           print*,'Grid dims (nx/ny)= ',nx_dom,ny_dom
           print*,'Lower left  i/j  = ',lli_orig,llj_orig
           print*,'Upper right i/j  = ',uri_orig,urj_orig

           if(nest.gt.1)then
              print*,'Calculate center lat/lon '
              La1=mother_lat(lli_orig,llj_orig)
              Lo1=mother_lon(lli_orig,llj_orig)
              La2=mother_lat(uri_orig,urj_orig)
              Lo2=mother_lon(uri_orig,urj_orig)
              print*,'*** SW lat/lon = ',La1,Lo1,'***'
              print*,'*** NE lat/lon = ',La2,Lo2,'***'
              call get_earth_radius(erad,istatus)
c Get X/Y for nested grid center
              call latlon_to_xy(La1,Lo1,ERAD,Xsw,Ysw) 
              call latlon_to_xy(La2,Lo2,ERAD,Xne,Yne) 
              Xcen=(Xne+Xsw)/2.
              Ycen=(Yne+Ysw)/2.
              call xy_to_latlon(Xcen,Ycen,erad,rlatcen,rloncen)
              print*,'*** Center lat/lon= ',rlatcen,rloncen,' ***'
           endif

           call get_n_staggers(n_staggers,istatus)
           call get_stagger_index(istag,istatus)

           if(istatus.ne.1)then
              print*,'Error getting namelist info for',
     +': ',c10_grid_f(1:lf)
              stop
           endif

           if(allocated(mother_lat))then
              deallocate (mother_lat, mother_lon, dum2d)
           endif

           print*
c
c -------------------------------------------------------------------
c
           call Gridmap_sub(nx_dom,ny_dom,n_staggers,istag,nest
     +,rlatcen,rloncen,lli_orig,llj_orig,uri_orig,urj_orig
     +,iparent_id,iratio_to_parent,istatus)
c
c -------------------------------------------------------------------
c
           print*,' Gridmap_sub finished '
           print*
           if(istatus.eq.1)then

              print*,' Gridmap_sub finish successfully'

              if( allocated(mother_lat) )then
               deallocate (mother_lat,mother_lon,dum2d,stat=istat)
               if(istat.ne.0)then
                  print*
                  print*,'Cannot deallocate mother arrays: ',istat
                  print*,'Abort'
                  print*
                  stop
               endif
              endif

           else
              print*,'Problem in Gridmap_sub. Terminating'
              stop
           endif

           else  !localize or not!

           print*
           print*,' No need to re-localize this domain '
           print*,' Localization namelist parms havent changed'
           print*

           if( allocated(mother_lat) )then
              deallocate (mother_lat,mother_lon,dum2d,stat=istat)
              if(istat.ne.0)then
                 print*
                 print*,'Cannot deallocate mother arrays: ',istat
                 print*,'Abort'
                 print*
                 stop
              endif
           endif

           endif !localize case

        enddo

 999	print*,' gridgen_model finish: istatus = ',istatus
        print*

        stop
        end
c
c --------------------------------------------------------------------
c
        subroutine Gridmap_sub(nnxp,nnyp,n_staggers,istag,nest
     +,mdlat,mdlon,lli,llj,iuri,iurj,iparent_id,iratio_to_parent
     +,istatus)

        include 'trigd.inc'

        use mem_namelist, ONLY: laps_cycle_time, l_fsf_gridgen

        logical exist,new_DEM
        logical lforce_ter_2_zero
!mp
	logical categorical, useland
!mp
        logical l_topo_wps, l_parse

        integer nnxp,nnyp,mode
        integer ngrids
        integer n_staggers
        integer nest
        integer iter
        integer istag

        Real  mdlat,mdlon
        Real  xmn(nnxp),ymn(nnyp)
        Real  xtn(nnxp,n_staggers)
        Real  ytn(nnyp,n_staggers)

        real, allocatable ::  topt_10   (:,:)
        real, allocatable ::  topt_10_s (:,:)
        real, allocatable ::  topt_10_ln(:,:)
        real, allocatable ::  topt_10_lt(:,:)
        real, allocatable ::  topt_30   (:,:)
        real, allocatable ::  topt_30_s (:,:)
        real, allocatable ::  topt_30_ln(:,:)
        real, allocatable ::  topt_30_lt(:,:)

        real  topt_out   (nnxp,nnyp)
        real  topt_out_s (nnxp,nnyp)
        real  topt_out_ln(nnxp,nnyp)
        real  topt_out_lt(nnxp,nnyp)

        integer maxdatacat
        parameter (maxdatacat=24)

        real, allocatable ::  adum2d  (:,:)
        real, allocatable ::  adum3d  (:,:,:)

        real, allocatable ::  GEODAT2D(:,:)
        real, allocatable ::  GEODAT3D(:,:,:)

        real lats(nnxp,nnyp,n_staggers)
        real lons(nnxp,nnyp,n_staggers)

	real, allocatable, dimension(:,:):: hlat,hlon,vlat,vlon


        character (len=3),   allocatable :: var(:)
        character (len=125), allocatable :: comment(:)

        character*2   cnest
        character*30  ctopo,clatlon,corners,ctopog
        character*30  ctopos,clatlons,cornerss,ctopogs
        character*30  clatlon2d,clatlons2d
        character*131 model

        character*200 path_to_topt30s
        character*200 path_to_topt10m
        character*200 path_to_pctl10m
        character*200 path_to_soiltype_top_30s
        character*200 path_to_soiltype_bot_30s
        character*200 path_to_green_frac    !no reference to res since there is both
                                            !10m and 8.64m
        character*200 path_to_soiltemp_1deg
        character*200 path_to_luse_30s
        character*200 path_to_albedo        !albedo = 0.144 deg res or 8.64m 
        character*200 path_to_maxsnoalb
        character*200 path_to_islope

        character*255 filename
        character*255 filename_wps
        character*200 c_dataroot
        character*200 cdl_dir
        character*180 static_dir 
        character*10  c10_grid_f           ! Type of domain (nest7grid, wrfsi)
        character*10  c10_grid_fname       ! Actual filename, cdl (nest7grid for now)
        character*6   c6_maproj
        character*1   cdatatype
        integer len,lf,lfn,ns,avgelem,zinelem
        integer ishow_timer,init_timer
        integer itstatus

        real,  allocatable ::  rmapdata(:,:,:)
        real,  allocatable ::  data(:,:,:)     !primary output array storage unit.
 
        interface

          subroutine adjust_geog(nnxp,nnyp,ncat,ctype
     &,istat_dat,lat,topt_out,landmask,geog_data,istatus)

          integer nnxp,nnyp
          integer ncat
          integer istat_dat
          integer istatus

          character*(*) ctype
          real    lat(nnxp,nnyp)
          real    landmask(nnxp,nnyp)
          real    topt_out(nnxp,nnyp)
          real    geog_data(nnxp,nnyp,ncat)
          end subroutine

          subroutine proc_geodat(nx_dom,ny_dom,ncat
     1,path_to_tile_data,dom_lats_in,dom_lons_in,lmask_out
     1,geodat,istatus)

          integer nx_dom
          integer ny_dom
          integer ncat
          character*(*) path_to_tile_data
          real    dom_lats_in(nx_dom,ny_dom)
          real    dom_lons_in(nx_dom,ny_dom)
          real    lmask_out(nx_dom,ny_dom)
          real    geodat(nx_dom,ny_dom,ncat)
          integer istatus
          end subroutine

        end interface

C*********************************************************************

        itstatus = init_timer()

        call find_domain_name(c_dataroot,c10_grid_f,istatus)
        if(istatus .ne. 1)then
            write(6,*) 'Error: returned fro find_domain_name '
            return
        endif

        call s_len(c10_grid_f,lf)

!mp  - need to know map projection information from the start
        call get_c6_maproj(c6_maproj,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error calling laps routine'
            return
        endif
        write(6,*)' c6_maproj = ',c6_maproj
!mp

c add 12 for albedo month 1 - 12.
c add 1 more for islope JS: 2-25-03
        if(c10_grid_f(1:lf).eq.'wrfsi')then
!mp
            IF (c6_maproj .eq. 'rotlat') THEN
               write(6,*) 'ngrids set to 103'
               ngrids=103 !!!
            ELSE !normal
               ngrids=112     !97+12   ! 2d grids (including %dist for landuse
                                       ! and two soiltype categories, green frac, albedo and others).
	    ENDIF

        else
           ngrids=38      !26+12   doesn't include terrain slope index and other stuff.
        endif

        allocate (data(nnxp,nnyp,ngrids),stat=istat)
        
        if(istat.ne.0)then
           print*,'Error allocating data object: array data ',istat
           print*,'Terminating: Maybe not enough memory'
           stop
        endif

        model = 'MODEL 4 delta x smoothed filter\0'

        icount_10 = 0
        icount_30 = 0
        icount_ramp = 0

        call get_directory('static',static_dir,lens)

c iplttopo is 1 if you want to plot the topography
c   the 30s topo covers the continental US
cc        itoptfn_30=static_dir(1:len)//'model/topo_30s/U'
c   the 10m topo covers the world
cc        itoptfn_10=static_dir(1:len)//'model/topo_10m/H'

        call get_path_to_topo_10m(path_to_topt10m,istatus)
        if(istatus .ne. 1)then
           write(6,*) 'Error getting path_to_topt10m'
           return
        endif

        call get_path_to_topo_30s(path_to_topt30s,istatus)
        if(istatus .ne. 1)then
           write(6,*) 'Error getting path_to_topt30s'
           return
        endif

!       l_topo_wps = l_parse(path_to_topt30s,'wps')

        call get_path_to_soiltype_top(path_to_soiltype_top_30s
     +,istatus)
        if(istatus .ne. 1)then
           write(6,*)'Error getting path_to_soiltype_top_30s'
           return
        endif

        call get_path_to_soiltype_bot(path_to_soiltype_bot_30s
     +,istatus)
        if(istatus .ne. 1)then
           write(6,*)'Error getting path_to_soiltype_bot_30s'
           return
        endif

        call get_path_to_landuse_30s(path_to_luse_30s,istatus)
        if(istatus .ne. 1)then
           write(6,*)'Error getting path_to_luse_30s'
           return
        endif

        call get_path_to_green_frac(path_to_green_frac,istatus)
        if(istatus .ne. 1)then
           print*,'Error getting path_to_green_frac'
           return
        endif

        call get_path_to_soiltemp_1deg(path_to_soiltemp_1deg
     &,istatus)
        if(istatus .ne. 1)then
           print*, 'Error getting path_to_soiltemp_1deg'
           return
        endif

        call get_path_to_albedo(path_to_albedo,istatus)
        if(istatus .ne. 1)then
           print*, 'Error getting path_to_albedo'
           return
        endif

        call get_path_to_maxsnoalb(path_to_maxsnoalb,istatus)
        if(istatus .ne. 1)then
           print*, 'Error getting path_to_maxsnoalb'
           return
        endif

        call get_path_to_islope(path_to_islope,istatus)
        if(istatus .ne. 1)then
           print*, 'Error getting path_to_islope'
           return
        endif

        call s_len(path_to_topt10m,len)
        if(len.gt.0)then
           print*,'path to toptl0m:        ',path_to_topt10m(1:len)
        endif
        path_to_topt10m(len+1:len+2)='/H'

        call s_len(path_to_soiltype_top_30s,len)
        print*,'path to soiltype_top:   '
     .,path_to_soiltype_top_30s(1:len)
        path_to_soiltype_top_30s(len+1:len+2)='/O'
        print*,'path to soiltype_bot:   '
     .,path_to_soiltype_bot_30s(1:len)
        path_to_soiltype_bot_30s(len+1:len+2)='/O'

        call s_len(path_to_luse_30s,len)
        print*,'path to landuse 30s: ',path_to_luse_30s(1:len)
        path_to_luse_30s(len+1:len+2)='/V'

        call s_len(path_to_green_frac,len)
        print*,'path to green frac:     ',path_to_green_frac(1:len)
        path_to_green_frac(len+1:len+2)='/G'

        call s_len(path_to_soiltemp_1deg,len)
        print*,'path to soiltemp 1deg:  ',path_to_soiltemp_1deg(1:len)
        path_to_soiltemp_1deg(len+1:len+2)='/T'

        call s_len(path_to_albedo,len)
        print*,'path to albedo:         ',path_to_albedo(1:len)
        path_to_albedo(len+1:len+2)='/A'

        call s_len(path_to_maxsnoalb,len)
        print*,'path to max snow albedo: ',path_to_maxsnoalb(1:len)
        path_to_maxsnoalb(len+1:len+2)='/M'

        call s_len(path_to_islope,len)
        print*,'path to islope categorical: ',path_to_islope(1:len)
        path_to_islope(len+1:len+2)='/I'

        call get_topo_parms(silavwt_parm,toptwvl_parm,istatus)
	if (istatus .ne. 1) then
           write (6,*) 'Error getting terrain smoothing parms'
	   return
	endif

        call get_r_missing_data(r_missing_data,istatus)
	if (istatus .ne. 1) then
           write (6,*) 'Error getting r_missing_data'
	   return
	endif

!       Silhouette weighting parameter
        silavwt=silavwt_parm

!       Terrain wavelength for filtering
        toptwvl=toptwvl_parm

        iplttopo=1

C*********************************************************************
        call get_gridnl(mode) 

c calculate delta x and delta y using grid and map projection parameters
        call get_standard_latitudes(std_lat,std_lat2,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error calling laps routine'
            return
        endif
        write(6,*)' standard_lats = ',std_lat,std_lat2

	IF (c6_maproj .ne. 'rotlat') THEN

        call get_grid_spacing(grid_spacing_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error calling laps routine'
            return
        endif

	ELSE
	
        write(6,*)
     1  'ERROR: get_nmm_grid_spacing not currently supported for rotlat'       
        stop
!       call get_nmm_grid_spacing(dlmd,dphd,grid_spacing_m)

	ENDIF

        write(6,*)' grid_spacing = ',grid_spacing_m

        if(nest.eq.1)then
           print*,'get grid center from namelist'
           call get_grid_center(mdlat,mdlon,istatus)
           if(istatus .ne. 1)then
              write(6,*)' Error calling laps routine'
              return
           endif
        endif
          
        print*,' Gridmap_sub: grid_center = ',mdlat,mdlon

        call get_c6_maproj(c6_maproj,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error calling laps routine'
            return
        endif
        write(6,*)' c6_maproj = ',c6_maproj

        call get_standard_longitude(std_lon,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error calling laps routine'
            return
        endif
        write(6,*)' std_lon = ',std_lon

        call get_earth_radius(erad,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error calling laps routine'
            return
        endif
        write(6,*)' Earth Radius = ',erad

        if(c6_maproj .eq. 'plrstr')then
            call get_ps_parms(std_lat,std_lat2,grid_spacing_m,phi0
     1                       ,grid_spacing_proj_m)

            deltax = grid_spacing_proj_m

            if(phi0 .eq. 90.)then
                write(6,*)' Projection is tangent to earths surface'
            else
                write(6,*)' Projection is secant to earths surface'
            endif

            if(std_lat2 .eq. +90.)then
                write(6,*)' Note, grid spacing will equal '
     1                    ,grid_spacing_m,' at a latitude of ',std_lat
!               write(6,*)' deltax, deltay ',deltax,deltay
!    1                   ,' at the north pole'

            elseif(std_lat2 .eq. -90.)then
                write(6,*)' Note, grid spacing will equal '
     1                    ,grid_spacing_m,' at a latitude of ',std_lat       
!               write(6,*)' deltax, deltay ',deltax,deltay
!    1                   ,' at the south pole'

            else ! abs(std_lat2) .ne. 90. (local stereographic)
                write(6,*)' The standard latitude ',std_lat,' is'
     1                   ,' relative to where the pole'
                write(6,*)' of the map projection is: lat/lon '
     1                   ,std_lat2,std_lon
!               write(6,*)' deltax, deltay ',deltax,deltay
!    1                   ,' at the projection pole'
                if(std_lat .ne. +90.)then
                    write(6,*)' Note: std_lat should usually be set'
     1                       ,' to +90. for local stereographic'
                endif

            endif

        else ! c6_maproj .ne. 'plrstr'
            deltax = grid_spacing_m

        endif ! c6_maproj .eq. 'plrstr'

        deltay = deltax
        write(6,*)' deltax, deltay (in projection plane) '
     1             ,deltax,deltay
 
c*********************************************************************
c in arrays lats/lons, the first stagger is actually not staggered
c but is the usual lat/lon values at non-staggered grid points. These
c are used for the analysis (LAPS-nest7grid). For WRFSI, we use the
c mass c-stagger (#4); staggered in both x and y. 

!mp


	IF (c6_maproj .eq. 'rotlat') THEN
          write(6,*) 'dlmd, dphd, mdlat, mdlon: ',dlmd,dphd,mdlat,mdlon       
                                                                                
          allocate(hlat(nnxp,nnyp),hlon(nnxp,nnyp))
          allocate(vlat(nnxp,nnyp),vlon(nnxp,nnyp))
             
          call etall(nnxp,nnyp,mdlat,mdlon,dlmd,dphd,hlat,hlon
     &                  ,vlat,vlon)
             
	  write(6,*) 'EGRID CORNERS'
          write(6,*) 'hlat(1,1),hlon(1,1): ', hlat(1,1),hlon(1,1)
          write(6,*) 'hlat(1,jm),hlon(1,jm): ',hlat(1,nnyp),hlon(1,nnyp)       
          write(6,*) 'hlat(im,1),hlon(im,1): ', hlat(nnxp,1),
     &			hlon(nnxp,1)
          write(6,*) 'hlat(im,jm),hlon(im,jm): ', hlat(nnxp,nnyp),
     &			hlon(nnxp,nnyp)
             
             
          DO J=1,NNYP
            DO I=1,NNXP
                data(I,J,1)=hlat(I,J)
                data(I,J,2)=hlon(I,J)
                data(I,J,3)=vlat(I,J)
                data(I,J,4)=vlon(I,J)
            ENDDO
          ENDDO

	ELSEIF(c6_maproj .eq. 'icshdr')then ! Read Icosahedral grid from a file
          write(6,*)' Read Icosahedral grid from a file'

          filename=static_dir(1:lens)//'ltln5C.txt'
          open(7,file=filename)
          read(7,*) lats, lons
          close(7) 

          rpd = 3.1415926535897932 / 180.

          DO J=1,NNYP
          DO I=1,NNXP
          DO IS = 1,n_staggers
              lats(i,j,is) = lats(i,j,is) / rpd
              lons(i,j,is) = lons(i,j,is) / rpd
              if(lons(i,j,is) .lt. -180.)then
                  lons(i,j,is) = lons(i,j,is) + 360.
              endif
              if(lons(i,j,is) .gt. +180.)then
                  lons(i,j,is) = lons(i,j,is) - 360.
              endif
          ENDDO
          ENDDO
          ENDDO

        ELSE ! Conformal grid ('plrstr','merctr',lambrt')
          call compute_latlon(nnxp,nnyp,n_staggers,mdlat,mdlon
     +                         ,deltax,xtn,ytn,lats,lons,istatus)
          if(istatus.ne.1)then
             print*,'Error returned: compute_latlon'
             return
          endif
             
        ENDIF

        IF (c6_maproj .ne. 'rotlat' .and. c6_maproj .ne. 'icshdr') THEN       

        ns = istag

        print*,'ns [staggers array index] = ',ns
        if(ns.ne.1)
     &print*,' Using C-stagger grid for LSM fields and other output'

C*****************************************************************
        write(6,*)
        write(6,*)'Corner points.'

        call latlon_to_xy(lats(1,1,ns),lons(1,1,ns)      ,erad,xco,yco)
        write(6,701,err=702)1,1,   lats(1,1,ns)   ,lons(1,1,ns),xco,yco 

        call latlon_to_xy(lats(nnxp,1,ns),lons(nnxp,1,ns),erad,xco,yco)
        write(6,701,err=702)nnxp,1,lats(nnxp,1,ns),lons(nnxp,1,ns)
     +                                                         ,xco,yco   

        call latlon_to_xy(lats(1,nnyp,ns),lons(1,nnyp,ns),erad,xco,yco)
        write(6,701,err=702)1,nnyp,lats(1,nnyp,ns),lons(1,nnyp,ns)
     +                                                         ,xco,yco   

        call latlon_to_xy(lats(nnxp,nnyp,ns),lons(nnxp,nnyp,1)
     +                                                   ,erad,xco,yco)
        write(6,701,err=702)nnxp,nnyp,lats(nnxp,nnyp,ns)
     +                                       ,lons(nnxp,nnyp,1),xco,yco   

 701    format(' lat/lon/x/y at ',i5,',',i5,' =',2f12.5,2f12.0)
 702    continue

        call check_domain(lats(1,1,ns),lons(1,1,ns),nnxp,nnyp
     +,grid_spacing_m,1,istat_chk)  
        if(istat_chk.ne.1)then
           print*,'Error returned from check_domain'
           istatus = istat_chk
           return
        endif

	ENDIF ! NOTE: totally avoided check_domain for rotlat & icshdr


       if(.true.)then       
           print*,'get perimeter of grid'
           call get_domain_perimeter_grid(nnxp,nnyp,c10_grid_f
     1                  ,lats(1,1,ns),lons(1,1,ns)
     1                  ,0.0,rnorth,south,east,west,istatus)
           print*,'static dir = ',static_dir(1:lens)
           open(10,file=static_dir(1:lens)//'/llbounds.dat'
     +         ,status='unknown')

           print*,'write llbounds.dat'
           write(10,*)rnorth                 
           write(10,*)south          
           write(10,*)east
           write(10,*)west
           close(10)
       endif

       write(6,*)'deltax = ',deltax
       write(6,*)'check_domain:status = ',istat_chk
c
C*****************************************************************
c calculate surface static fields
c
c ------------------------------------------------------
c ----------------- Terrain USGS 30 sec ----------------
c type = U
c
       itstatus=ishow_timer()
       print*
       print*,' Processing 30s topo data, l_fsf_gridgen = '
     1                                   ,l_fsf_gridgen      

       allocate (topt_30(nnxp,nnyp),
     +           topt_30_s(nnxp,nnyp),
     +           topt_30_ln(nnxp,nnyp),
     +           topt_30_lt(nnxp,nnyp))

       l_topo_wps = l_fsf_gridgen

 600   if(.not. l_topo_wps)then

        call s_len(path_to_topt30s,len)
        print*,'path to topt30s:        ',path_to_topt30s(1:len)
        path_to_topt30s(len+1:len+2)='/U'

        IF (c6_maproj .eq. 'rotlat') THEN
	 categorical=.false.

	 NCAT=1
	 allocate(adum3d(NNXP,NNYP,NCAT))
         allocate(adum2d(nnxp,nnyp))

         CALL alt_10by10_all(nnxp,nnyp,1200,grid_spacing_m/1000.,
     &   hlat,hlon,PATH_TO_TOPT30S,topt_30,adum3d,NCAT,topt_30_lt,
     &   topt_30_ln,topt_30_s,categorical,2,0)


	 write(6,*) 'topo before smooth'
         do J=nnyp,1,-nnyp/30
           write(6,633) (topt_30(I,J),I=1,nnxp,nnxp/15)
         enddo


         deallocate(adum3d)
                                                                                
  633    format(20(f5.0,1x))
                                                                                
!mp	 bogus status value to keep things smooth below
         istatus_30s=1
!mp

        ELSEIF (c6_maproj .eq. 'icshdr') THEN
         allocate  (GEODAT3D(nnxp,nnyp,1))
         call proc_geodat(nnxp,nnyp,1,path_to_topt30s
     +       ,lats(1,1,1),lons(1,1,1),data(1,1,1)
     +       ,GEODAT3D,istatus_30s)

        ELSE

         CALL GEODAT(nnxp,nnyp,erad,90.,std_lon,xtn(1,1)
     +   ,ytn(1,1),deltax,deltay,TOPT_30,TOPT_30_S,TOPT_30_LN
     +   ,TOPT_30_LT,PATH_TO_TOPT30S,TOPTWVL,SILAVWT,new_DEM,1
     +   ,istatus_30s)

        ENDIF

        if(istatus_30s .ne. 1)then

          print*,'WARNING: File(s) missing for 30s terrain data'
          print*,' >>>> Try Processing 10m topo data <<<<'
          print*

          allocate (topt_10(nnxp,nnyp),
     +              topt_10_s(nnxp,nnyp),
     +              topt_10_ln(nnxp,nnyp),
     +              topt_10_lt(nnxp,nnyp))

          print*
          print*,' Calling GEODAT:  10m terrain data ....'

          CALL GEODAT(nnxp,nnyp,erad,90.,std_lon,xtn(1,1)
     +,ytn(1,1),deltax,deltay,TOPT_10,TOPT_10_S,TOPT_10_LN
     +,TOPT_10_LT,PATH_TO_TOPT10M,TOPTWVL,SILAVWT,new_DEM,1
     +,istatus_10m)

          if(istatus_10m .ne. 1)then

             print*
             print*,'Error: File(s) missing for both 10m and
     +30s terrain data'
             print*,'Error: Aborting --- No static file created'
             return

          else

             print *,'topt_10    =',topt_10(1,1),topt_10(nnxp,nnyp)
             print*
             print*,'try blending 10min and 30s terrain data'
             call blend_topo(nnxp,nnyp,lats(1,1,1),lons(1,1,1)
     1,topt_10,topt_10_s,topt_10_ln,topt_10_lt
     1,topt_30,topt_30_s,topt_30_ln,topt_30_lt
     1,topt_out,topt_out_s,topt_out_ln,topt_out_lt)

             deallocate (topt_10,
     1                   topt_10_s,
     1                   topt_10_ln,
     1                   topt_10_lt)

          endif

        else ! go with 30s topo data since all tiles available

          print*,'Topo SW and NE corner values'
          print *,'topt_30    =',topt_30(1,1),topt_30(nnxp,nnyp)

          topt_out=topt_30
          topt_out_s=topt_30_s
          topt_out_ln=topt_30_ln
          topt_out_lt=topt_30_lt
          icount_30=nnyp*nnxp
          
        endif !istatus_30s ... data processed ok 

       elseif(.false.)then ! read topo data from wps output netCDF file
          write(6,*)' calling read_wrfstatic for wps topo'
          call s_len(path_to_topt30s,lenp)
          filename_wps = path_to_topt30s(1:lenp)//'/geo_em.d01.nc'
          call read_wrfstatic(nnxp,nnyp,lat,lon,filename_wps,topt_out
     1                       ,istatus)
          if(istatus .ne. 1)then
              write(6,*)' Error - no wrf static data: returning'
              return
          endif

       else                ! read topo data from LAPS/WRF FSF file
          call get_systime_i4(i4time_sys,istatus)
          if(istatus .ne. 1)then
              write(6,*)' Error - no systime: returning'
              return
          endif

          istatus = 0
          ntrys = 40
          itry = 0
          do while (itry .le. ntrys .AND. istatus .eq. 0)
              i4time_try = i4time_sys - (itry * laps_cycle_time)
              write(6,*)' Trying for FSF data at',itry,i4time_try
              call get_modelfg_2d(i4time_try,'TER',nnxp,nnyp,topt_out
     1                                                      ,istatus)       
              itry = itry + 1
          enddo

          if(itry .eq. ntrys)then
              write(6,*)' Reached maximum number of tries'
              write(6,*)' Use geog topo data instead'
              l_topo_wps = .false.
              goto 600
          endif

       endif

       deallocate (topt_30,
     1             topt_30_s,
     1             topt_30_ln,
     1             topt_30_lt)

       print*
c
c ---------------------------------------------------------
c ------------------- 30s Landuse -------------------------
c type = V
c
       allocate  (GEODAT2D(nnxp,nnyp))
       allocate  (GEODAT3D(nnxp,nnyp,24))

        IF (c6_maproj .ne. 'rotlat') THEN

       if(.not.allocated(adum2d))allocate (adum2d(nnxp,nnyp))

       itstatus=ishow_timer()
       print*
       print*,' Calling GEODAT: Processing 30s landuse data.'
       print*,' Re-allocate GEODAT3D ',nnxp,nnyp,' 24'
       print*

       CALL GEODAT(nnxp,nnyp,erad,90.,std_lon,xtn(1,ns)
     +,ytn(1,ns),deltax,deltay,GEODAT2D,GEODAT3D
     +,adum2d,adum2d,PATH_TO_LUSE_30S,2.0,0.0,new_DEM,24
     +,istatus)

       if(istatus.ne.1)then
         print*,'********* ERROR ***********'
         print*,'ERROR:  File(s) missing for landuse data'
         print*,'ERROR:  Static file not created'
         print*
         return
       endif
	
	ELSE

	categorical=.true.
        call alt_10by10_all(nnxp,nnyp,1200,grid_spacing_m/1000.,
     &		hlat,hlon,PATH_TO_LUSE_30S,GEODAT2D,
     &          GEODAT3D,24,adum2d,adum2d,adum2d,categorical,1,9999)
 
        do J=nnyp,1,-nnyp/25
        write(6,733) (INT(GEODAT2D(I,J)),I=1,nnxp,nnxp/15)
        enddo
        
  733   format(50(I2,1x))
 
	ENDIF

       ilndmsk=12  !land mask is 12 for both wrfsi and laps
       if(c10_grid_f(1:lf).eq.'wrfsi')then

	 IF (c6_maproj .ne. 'rotlat') THEN

c dominant category landuse
          ilnduse=11
          data(:,:,ilnduse)=GEODAT2D
c landmask for wrfsi
          data(:,:,ilndmsk) = 1.
          where(data(:,:,11) .eq. 16.)data(:,:,ilndmsk)=0.

c grids 15 thru 38 are percent distributions
          i=14
          do j=1,24
             data(:,:,i+j)=GEODAT3D(:,:,j)
          enddo
	ELSE ! rotlat case
                                                                                
        ilndmsk=7
        ilnduse=6
                                                                                
c landmask for wrfsi
           data(:,:,ilnduse)=GEODAT2D
           data(:,:,ilndmsk) = 1.
           where(data(:,:,ilnduse) .eq. 16.)data(:,:,ilndmsk)=0.


!!!! require a large % of water to make a water point
!!	
	fraclk=0.7  !! require that 70% of points be water for output to be h2o

	DO JJ=1,nnyp
	 DO II=1,nnxp
	  IF (GEODAT3D(II,JJ,16) .gt. fraclk) THEN
	    data(II,JJ,ilndmsk)=0.
          ELSE	
	    IF ( data(II,JJ,ilndmsk) .eq. 0.) THEN
!	      write(6,*) 'WATER WATER, BECOMING LAND:: ', II,JJ
	    ENDIF
	    data(II,JJ,ilndmsk)=1.
          ENDIF	
         ENDDO	
	ENDDO


!!	change dominant landuse  (otherwise will have land points with
!!	water as the dominant class)

	DO JJ=2,NNYP-1
        DO II=2,NNXP-1 

	if (data(II,JJ,ilndmsk).eq.1.and.data(II,JJ,6).eq.16.) then
	onepos=amax1(data(II+1,JJ,6),data(II-1,JJ,6),
     +   data(II,JJ+1,6),data(II,JJ-1,6))
	twopos=amin1(data(II+1,JJ,6),data(II-1,JJ,6),
     +   data(II,JJ+1,6),data(II,JJ-1,6))

	if (twopos .ne. 16) then
	 	data(II,JJ,6)=twopos
	elseif (onepos .ne. 16) then
		data(II,JJ,6)=onepos
	else
!	leave data alone, revert the land mask to water
	write(6,*) 'switching back to water at II,JJ: ', II,JJ
		data(II,JJ,ilndmsk)=0.
	endif

	endif

	ENDDO
	ENDDO

                                                                                
! putting land-use percentages in 19-42 for NMM
        I=18
           do j=1,24
              data(:,:,i+j)=GEODAT3D(:,:,j)
           enddo
                                                                                
        ENDIF

       else
c dominant category landuse
          ilnduse=5
          data(:,:,ilnduse)=GEODAT2D   !landuse for laps
c landmask for laps
          data(:,:,ilndmsk)=1.
          where(data(:,:,5) .eq. 16.)data(:,:,ilndmsk)=0.
       endif
c
c **********************************************
c fix water areas to have 0.0m terrain elevation
c currently no fix going on!
c **********************************************
       lforce_ter_2_zero=.false.
       if(lforce_ter_2_zero)then
          where(data(:,:,ilndmsk).eq.0)topt_out=0.0
       endif


c for WRFSI we'll bilinearly interpolate to get the topo
c from the non-staggered topo

       if(c10_grid_f(1:lf).eq.'wrfsi')then

	IF (c6_maproj .ne. 'rotlat') THEN

          data(:,:,9)=topt_out
          call move(topt_out_s,data(1,1,48),nnxp,nnyp)
          call move(topt_out_ln,data(1,1,49),nnxp,nnyp)
          call move(topt_out_lt,data(1,1,50),nnxp,nnyp)

          indxter=51
          data(:,:,indxter)=r_missing_data
          do j=2,nnyp
          do i=2,nnxp
             call bilinear_interp(i,j,nnxp,nnyp,topt_out,result)
             if(data(i,j,51).gt.30000. .and.
     & data(i,j,51).lt.r_missing_data)then
                print*,'Here: Large topo: i/j/topo ',i,j,data(i,j,51)
             endif
             data(i-1,j-1,indxter)=result
          enddo
          enddo

          where((data(:,:,indxter).lt. 0.01).and.
     1          (data(:,:,indxter).gt.-0.01))data(:,:,indxter)=0.0

	ELSE

!       for rotlat, topt_out is already on mass point stagger.
!       No bilinear interp needed


!!	boundary smooth the topo data

	write(6,*) 'call smdhld'
	write(6,*) 'NNXP, NNYP: ', NNXP,NNYP
	call smdhld(nnxp,nnyp,topt_out,1.-data(:,:,ilndmsk),12,12)

c       write(6,*) 'topo after smoothing'
c       do J=nnyp,1,-nnyp/30
c       write(6,633) (topt_out(I,J),I=1,nnxp,nnxp/15)
c       enddo


!-------------4-point averaging of mountains along inner boundary-------

          do i=1,NNXP-1
      topt_out(i,2)=0.25*(topt_out(i,1)+topt_out(i+1,1)+ 
     +                    topt_out(i,3)+topt_out(i+1,3))
          enddo

          do i=1,NNXP-1
      topt_out(i,NNYP-1)=0.25*(topt_out(i,NNYP-2)+topt_out(i+1,NNYP-2)+ 
     +                         topt_out(i,NNYP)+topt_out(i+1,NNYP))
          enddo

          do j=4,NNYP-3,2
      topt_out(1,j)=0.25*(topt_out(1,j-1)+topt_out(2,j-1)+
     +                    topt_out(1,j+1)+topt_out(2,j+1))
          enddo

          do j=4,NNYP-3,2
      topt_out(NNXP,j)=0.25*(topt_out(NNXP-1,j-1)+topt_out(NNXP,j-1)+
     +                         topt_out(NNXP-1,j+1)+topt_out(NNXP,j+1))
          enddo


	indxter=18
        do J=1,nnyp
        do I=1,nnxp
          data(i,j,indxter)=topt_out(i,j)
        enddo
        enddo

          where((data(:,:,indxter).lt. 0.01).and.
     1          (data(:,:,indxter).gt.-0.01))data(:,:,indxter)=0.0

	ENDIF


       else  !must be LAPS

          indxter=3
          call move(topt_out,data(1,1,indxter),nnxp,nnyp)         ! KWD
          call move(topt_out_s,data(1,1,7),nnxp,nnyp)             ! JS
          call move(topt_out_ln,data(1,1,8),nnxp,nnyp)            ! JS
          call move(topt_out_lt,data(1,1,9),nnxp,nnyp)            ! JS

       endif 

       print*
       print *,'terrain  = ',data(1,1,indxter),data(nnxp,nnyp,indxter)

c      print *,'# of grid pts using 30 sec terrain data =  ',icount_30
c      print *,'# of grid pts using blended terrain data = ',icount_ramp
       print*
C
C Output Section for topo, lat-lons and corners
C -----------------------------------------------
       if(c10_grid_f(1:lf).eq.'wrfsi')then
         if(c6_maproj.ne.'rotlat')then
          write(cnest,'(i2.2)')nest
          clatlon='latlon.d'//cnest//'.dat'            !binary non-staggered lat-lon file
          clatlons='latlon-mass.d'//cnest//'.dat'      !binary staggered lat-lon file
          ctopo='topo.d'//cnest//'.dat'                !binary non-staggered topo file
          ctopos='topo-mass.d'//cnest//'.dat'          !binary staggered topo file
          ctopog='topography.d'//cnest//'.dat'         !ascii non-staggered topo file
          ctopogs='topography-mass.d'//cnest//'.dat'   !ascii staggered topo file
          clatlon2d='latlon2d.d'//cnest//'.dat'        !ascii lat-lon non-staggered file
          clatlons2d='latlon2d-mass.d'//cnest//'.dat'  !ascii lat-lon staggered file
          corners='corners.d'//cnest//'.dat'           !non-staggered corners
          cornerss='corners-mass.d'//cnest//'.dat'     !staggered corners
         else
          clatlon='latlon-rotlat.dat'        !binary lat/lon file
          clatlon2d='latlon2d-rotlat.dat'    !ascii lat/lon file
          ctopo='topo-rotlat.dat'            !binary topo file
          ctopog='topography-rotlat.dat'     !ascii topo file
          corners='corners-rotlat.dat'       !domain corners file
         endif
       else
          clatlon='latlon.dat'        !binary lat/lon file
          clatlon2d='latlon2d.dat'    !ascii lat/lon file
          ctopo='topo.dat'            !binary topo file
          ctopog='topography.dat'     !ascii topo file
          corners='corners.dat'       !domain corners file
       endif

       call s_len(static_dir,len)
       if(c10_grid_f(1:lf).eq.'wrfsi')then

         if(c6_maproj.ne.'rotlat')then
           in=0
c write both non-stagger and stagger data for wrfsi, non-rotlat
           do k=1,2
            if(k.eq.1)then
             is=1
             idxt=9
             filename=static_dir(1:len)//clatlon
             open(10,file=filename,status='unknown',form='unformatted')
             filename=static_dir(1:len)//ctopo
             open(11,file=filename,status='unknown',form='unformatted')
             filename=static_dir(1:len)//corners
             open(15,file=filename,status='unknown')
             filename=static_dir(1:len)//ctopog
             open(666,file=filename)
             filename=static_dir(1:len)//clatlon2d
             open(667,file=filename)
            else
             is = 4
             in = 1
             idxt=51
             filename=static_dir(1:len)//clatlons
             open(10,file=filename,status='unknown',form='unformatted')
             filename=static_dir(1:len)//ctopos
             open(11,file=filename,status='unknown',form='unformatted')
             filename=static_dir(1:len)//cornerss
             open(15,file=filename,status='unknown')
             filename=static_dir(1:len)//ctopogs
             open(666,file=filename)
             filename=static_dir(1:len)//clatlons2d
             open(667,file=filename)
            endif

            write(10)lats(1:nnxp-in,1:nnyp-in,is)
     +,lons(1:nnxp-in,1:nnyp-in,is)
            close(10)
            write(11)data(1:nnxp-in,1:nnyp-in,idxt)
            close(11)
            write(15,*)lats(1,1,is),lons(1,1,is)
            write(15,*)lats(1,nnyp-in,is),lons(1,nnyp-in,is)
            write(15,*)lats(nnxp-in,1,is),lons(nnxp-in,1,is)
            write(15,*)lats(nnxp-in,nnyp-in,is),lons(nnxp-in,nnyp-in,is)
            close(15)

            do j=1,nnyp-in
            do i=1,nnxp-in
               write(666,*) data(i,j,idxt)
               write(667,*) lats(i,j,is),lons(i,j,is)
            enddo
            write(666,'()')
            write(667,'()')
            enddo
            close(666)
            close(667)

           enddo

         else   !if(c6_maproj.eq.'rotlat')then

           is=1
           idxt=18
           filename=static_dir(1:len)//clatlon
           open(10,file=filename,status='unknown',form='unformatted')
           filename=static_dir(1:len)//ctopo
           open(11,file=filename,status='unknown',form='unformatted')
           filename=static_dir(1:len)//corners
           open(15,file=filename,status='unknown')
           filename=static_dir(1:len)//ctopog
           open(666,file=filename)
           filename=static_dir(1:len)//clatlon2d
           open(667,file=filename)

           write(10)lats(1:nnxp,1:nnyp,is),lons(1:nnxp,1:nnyp,is)
           close(10)
           write(11)data(:,:,idxt)
           close(11)
           write(15,*)lats(1,1,is),lons(1,1,is)
           write(15,*)lats(1,nnyp,is),lons(1,nnyp,is)
           write(15,*)lats(nnxp,1,is),lons(nnxp,1,is)
           write(15,*)lats(nnxp,nnyp,is),lons(nnxp,nnyp,is)
           close(15)

           do j=1,nnyp
           do i=1,nnxp
              write(666,*) data(i,j,idxt)
              write(667,*) lats(i,j,is),lons(i,j,is)
           enddo
           write(666,'()')
           write(667,'()')
           enddo
           close(666)
           close(667)

         endif

       else !laps case

         is=1
         idxt=3
!        filename=static_dir(1:len)//clatlon
!        open(10,file=filename,status='unknown',form='unformatted')
!        filename=static_dir(1:len)//ctopo
!        open(11,file=filename,status='unknown',form='unformatted')
         filename=static_dir(1:len)//corners
         open(15,file=filename,status='unknown')
         filename=static_dir(1:len)//ctopog
         open(666,file=filename)
         filename=static_dir(1:len)//clatlon2d
         open(667,file=filename)

!        write(10)lats(1:nnxp,1:nnyp,is),lons(1:nnxp,1:nnyp,is)
!        close(10)
!        write(11)data(:,:,idxt)
!        close(11)
         write(15,*)lats(1,1,is),lons(1,1,is)
         write(15,*)lats(1,nnyp,is),lons(1,nnyp,is)
         write(15,*)lats(nnxp,1,is),lons(nnxp,1,is)
         write(15,*)lats(nnxp,nnyp,is),lons(nnxp,nnyp,is)
         close(15)

         do j=1,nnyp
         do i=1,nnxp
            write(666,*) data(i,j,idxt)
            write(667,*) lats(i,j,is),lons(i,j,is)
         enddo
         write(666,'()')
         write(667,'()')
         enddo
         close(666)
         close(667)

       endif
c ----------------------------------------------------------
c
c ------ Land Fraction from 30s Land Use -------------------
c type = L (no longer used ... no land fraction data base)
c
        print*
        print*,' Create from 30s land use fractional dist'
        print*
 
        GEODAT2D(:,:)=1.- GEODAT3D(:,:,16)
c
c this function allows variable amounts of smoothing of
c the land fraction data depending on grid spacing
c
!       iter=min(10,int(8000./deltax))
        iter = 0

	IF (c6_maproj .ne. 'rotlat') THEN
        
          print*,'Filter land fraction with 2dx ',iter
          do k=1,iter
             call filter_2dx(GEODAT2D,nnxp,nnyp,1, 0.5)
             call filter_2dx(GEODAT2D,nnxp,nnyp,1,-0.5)
             GEODAT2D = max(GEODAT2D,0.) 
          enddo

	ENDIF

        if(c10_grid_f(1:lf).eq.'wrfsi')then
	
	    IF (c6_maproj .ne. 'rotlat') THEN
		idx=10
	    ELSE
		idx=5
	    ENDIF

	else ! not wrfsi
		idx=4
	endif

        data(:,:,idx)=GEODAT2D
        print *,'pctlandfrac=',GEODAT2D(1,1),GEODAT2D(nnxp,nnyp)       
c
c -------------------------------------------------------------------
c potential fix of fictitous islands for certain resolution domains.
c Story here is that west coast terrain can be "funky" possibly due
c to steepness and method to compute avg terrain.

	write(6,*) 'idx, indxter: ', idx, indxter
        where(data(:,:,idx).le.0.1 .and. data(:,:,indxter).lt.5.0)
     &data(:,:,indxter)=0.0

c
c ------------------------------------------------------------
c -------------- Top Layer Soil Type -------------------------
c type = O
c
        itstatus=ishow_timer()
        print*
        print*,' Processing 30s soil type top layer data....'

        deallocate(GEODAT3D)
        allocate (GEODAT3D(nnxp,nnyp,16))

        IF (c6_maproj .ne. 'rotlat') THEN

        CALL GEODAT(nnxp,nnyp,erad,90.,std_lon,xtn(1,ns),ytn(1,ns)
     +,deltax,deltay,GEODAT2D,GEODAT3D,adum2d,adum2d
     +,PATH_TO_SOILTYPE_TOP_30S,2.0,0.0,new_DEM,16   !maxdatacat
     +,istatus_soil)

        if(istatus_soil.ne.1)then
           print*
           print*,' Soil type data not processed completely'
           if(c10_grid_f(1:lf).eq.'wrfsi')then
             print*,' File(s) missing for soil type top layer data'
             print*,' Error:  Static file not created'
             print*
             return
           else
             print*,' File(s) missing for soil type top layer data'
             print*,'           *** WARNING ***'
             print*,' Soil Type Top data not added to static file'
             print*
           endif
           GEODAT2D=r_missing_data
           GEODAT3D=r_missing_data
        endif

	ELSE

	categorical=.true.
        call alt_10by10_all(nnxp,nnyp,1200,grid_spacing_m/1000.,
     &		hlat,hlon,PATH_TO_SOILTYPE_TOP_30S,GEODAT2D,
     &          GEODAT3D,16,adum2d,adum2d,adum2d,categorical,1,9999)
                                                                                
        write(6,*) 'soiltype_top'
        do J=nnyp,1,-nnyp/25
        write(6,733) (INT(GEODAT2D(I,J)),I=1,nnxp,nnxp/15)
        enddo
                                                                                
        ENDIF

        if(c10_grid_f(1:lf).eq.'wrfsi')then

         IF (c6_maproj .ne. 'rotlat') THEN
c
c make water points = 0 for adjust_geog. later we'll put it back to original
c
           where(GEODAT2D.eq.14)GEODAT2D=0.0

           call adjust_geog(nnxp-1,nnyp-1,1,'soiltype',istatus_soil
     &,lats(1:nnxp-1,1:nnyp-1,ns),data(1:nnxp-1,1:nnyp-1,indxter)
     &,data(1:nnxp-1,1:nnyp-1,ilndmsk),GEODAT2D(1:nnxp-1,1:nnyp-1)
     &,istatus)

           where(GEODAT2D.eq.0.0)GEODAT2D=14.
           idxstl=13
           data(:,:,idxstl)=GEODAT2D
           i=51
           do j=1,16
              data(:,:,i+j)=GEODAT3D(:,:,j)
           enddo
	 ELSE
           idxstl=8
           data(:,:,8)=GEODAT2D
           I=42
           do j=1,16
              data(:,:,i+j)=GEODAT3D(:,:,j)
           enddo
         ENDIF

        else
c
c make water points = 0 for adjust_geog. After put it back to original
c
           where(GEODAT2D.eq.14)GEODAT2D=0.0

           call adjust_geog(nnxp,nnyp,1,'soiltype',istatus_soil
     &,lats(1,1,ns),data(1,1,indxter),data(1,1,ilndmsk),GEODAT2D
     &,istatus)

           where(GEODAT2D.eq.0.0)GEODAT2D=14.
           idxstl=10
           data(:,:,idxstl)=GEODAT2D

        endif
c
c --------------------------------------------------------------
c -------------- Bottom Layer Soil Type ------------------------
c type = O
c
        itstatus=ishow_timer()
        print*
        print*,' Processing 30s soil type bottom layer data'

        IF (c6_maproj .ne. 'rotlat') THEN

        CALL GEODAT(nnxp,nnyp,erad,90.,std_lon,xtn(1,ns)
     +,ytn(1,ns),deltax,deltay,GEODAT2D,GEODAT3D,adum2d,adum2d
     +,PATH_TO_SOILTYPE_BOT_30S,2.0,0.0,new_DEM,16   !maxdatacat
     +,istatus_soil)

        if(istatus_soil.ne.1)then
           print*
           print*,' Bottom layer soil data not processed completely'
           if(c10_grid_f(1:lf).eq.'wrfsi')then
             print*,' File(s) missing for soil type bot layer data'
             print*,' Error:  Static file not created'
             print*
             return
           else
             print*,' File(s) missing for soil type bot layer data'
             print*,'           *** WARNING ***'
             print*,' Soil Type Bot data not added to static file'
             print*
           endif
           GEODAT2D=r_missing_data
           GEODAT3D=r_missing_data
        endif

	ELSE

	categorical=.true.
        call alt_10by10_all(nnxp,nnyp,1200,grid_spacing_m/1000.,
     &		hlat,hlon,PATH_TO_SOILTYPE_BOT_30S,GEODAT2D,
     &          GEODAT3D,16,adum2d,adum2d,adum2d,categorical,1,9999)

        ENDIF
c
c soiltype bottom layer ... 68 thru 83
c
        if(c10_grid_f(1:lf).eq.'wrfsi')then

          IF (c6_maproj .ne. 'rotlat') THEN
c
c make water points = 0 for adjust_geog. After, put it back to original.
c
            where(GEODAT2D.eq.14)GEODAT2D=0.0

            call adjust_geog(nnxp-1,nnyp-1,1,'soiltype',istatus_soil
     &,lats(1:nnxp-1,1:nnyp-1,ns),data(1:nnxp-1,1:nnyp-1,indxter)
     &,data(1:nnxp-1,1:nnyp-1,ilndmsk),GEODAT2D(1:nnxp-1,1:nnyp-1)
     &,istatus)

            where(GEODAT2D.eq.0.)GEODAT2D=14.0
            idxsbl=14
            data(:,:,idxsbl)=GEODAT2D
            i=67
            do j=1,16
               data(:,:,i+j)=GEODAT3D(:,:,j)
            enddo
	  ELSE
            idxsbl=9
	    data(:,:,9)=GEODAT2D
            i=58
            do j=1,16
               data(:,:,i+j)=GEODAT3D(:,:,j)
            enddo
          ENDIF

        else
c
c make water points = 0 for adjust_geog. After, put it back to original.
c
            where(GEODAT2D.eq.14)GEODAT2D=0.0

            call adjust_geog(nnxp,nnyp,1,'soiltype',istatus_soil
     &,lats(1,1,ns),data(1,1,indxter),data(1,1,ilndmsk),GEODAT2D
     &,istatus)
            where(GEODAT2D.eq.0.)GEODAT2D=14.0
            idxsbl=11
            data(:,:,idxsbl)=GEODAT2D
        endif
c
c ---------------------------------------------------------
c ----------------- Greenness Fraction --------------------
c type = G
c
        itstatus=ishow_timer()
        print*
        print*,' Calling PROC_GEODAT: Process green frac data.'
        print*,' Re-allocate GEODAT3D ',nnxp,nnyp,' 12'
        print*

        deallocate(GEODAT3D)
        allocate  (GEODAT3D(nnxp,nnyp,12),stat=istat)
        if(istat.ne.0)then
           print*,'Error allocating data object GEODAT3D ',istat
           if(c10_grid_f(1:lf).eq.'wrfsi')then
              print*,'Error: Aborting process'
              print*,'Error: Maybe not enough memory'
              print*,'Error: Cannot allocate GEODAT3D'
              stop
           else
              print*,'Warning: Continue without green fraction'         
           endif
        endif

	IF (c6_maproj .ne. 'rotlat') THEN

        print*,'Calling proc_geodat'
        call proc_geodat(nnxp,nnyp,12,path_to_green_frac
     +,lats(1,1,ns),lons(1,1,ns),data(1,1,ilndmsk)
     +,GEODAT3D,istatus_grn)

	ELSE
                                                                                
	categorical=.FALSE.
	useland=.TRUE.
        call alt_hemi_all(nnxp,nnyp,hlat,hlon,
     &                       PATH_TO_GREEN_FRAC,
     &                       type,12,grid_spacing_m/1000.,GEODAT3D,
     &                       GEODAT2D,categorical,data(1,1,ilndmsk),
     &			     useland)


!	convert from 0-100% to 0-1 fraction

	do N=1,12
	do J=1,NNYP
	do I=1,NNXP
	    GEODAT3D(I,J,N)=GEODAT3D(I,J,N)/100.
	enddo
	enddo
	enddo
                                                                                
        istatus_grn=1 !bogus value
                                                                                
        ENDIF

        if(istatus_grn.ne.1)then
         print*
         print*,'greenness fraction data not processed completely'
         if(c10_grid_f(1:lf).eq.'wrfsi')then
            print*,' Error: File(s) missing for green frac data'
            print*,' Error: Static file not created'
            print*
            istatus=0
            return
         else
            print*,' Warning: File(s) missing for green frac data'
            print*,' Warning: green frac not added to static file'
            print*
         endif
         GEODAT3D=r_missing_data
        endif

        if(c10_grid_f(1:lf).eq.'wrfsi')then

           IF (c6_maproj .ne. 'rotlat') THEN

              igrn=83

c             call adjust_geog(nnxp-1,nnyp-1,12,'greenfrac',istatus_grn
c    &,lats(1:nnxp-1,1:nnyp-1,ns),data(1:nnxp-1,1:nnyp-1,indxter)
c    &,data(1:nnxp-1,1:nnyp-1,ilndmsk),GEODAT3D(1:nnxp-1,1:nnyp-1,:)
c    &,istatus)

           ELSE
               igrn=74
           ENDIF

        else

           igrn=12
          
           call adjust_geog(nnxp,nnyp,igrn,'greenfrac',istatus_grn
     &,lats(:,:,ns),data(:,:,indxter),data(:,:,ilndmsk),GEODAT3D
     &,istatus)
        endif
        do j=1,12
           data(:,:,igrn+j)=GEODAT3D(:,:,j)
        enddo
c
c annual max/min greenness fraction in domain
c --------------------------------------------
        if(c10_grid_f(1:lf).eq.'wrfsi')then
           dommaxgf=0.0
           dommingf=99.0
           print*,' compute max/min greenness frac at grid points'
           IF (c6_maproj .ne. 'rotlat') THEN
               mxgiwrf=110
               mngiwrf=111
           ELSE
               mxgiwrf=101
               mngiwrf=102
           ENDIF
           do j=1,nnyp
           do i=1,nnxp
              data(i,j,mxgiwrf)=MAXVAL(GEODAT3D(i,j,:))
              data(i,j,mngiwrf)=MINVAL(GEODAT3D(i,j,:))
              if(data(i,j,mxgiwrf).gt.dommaxgf)
     &dommaxgf=data(i,j,mxgiwrf)
              if(data(i,j,mngiwrf).lt.dommingf)
     &dommingf=data(i,j,mngiwrf)
           enddo
           enddo
           write(6,*)'Domain Annual Max/Min Green Fraction',
     &dommaxgf,dommingf
        endif
c
c ---------------------------------------------------------
c --------------- Deep Soil Temperature -------------------
c type = T
c
	IF (c6_maproj .ne. 'rotlat') THEN
        itstatus=ishow_timer()
        print*
        print*,' Call PROC_GEODAT: Process 1 deg soiltemp data.'

        GEODAT2D=0.0

        call proc_geodat(nnxp,nnyp,1
     1,path_to_soiltemp_1deg,lats(1,1,ns),lons(1,1,ns)
     1,data(1,1,ilndmsk),GEODAT2D,istatus_tmp)

	ELSE

	write(6,*) 'call alt_10by10 for soiltemp'
        GEODAT2D=0.0
	categorical=.false.

	NCAT=1
        call alt_10by10_all(nnxp,nnyp,10,grid_spacing_m/1000.,
     &		hlat,hlon,PATH_TO_SOILTEMP_1DEG,GEODAT2D,
     &          ADUM3D,NCAT,adum2d,adum2d,adum2d,categorical,2,0)

!       convert to degrees
                                                                                
	write(6,*) 'convert to degrees ', nnxp,nnyp
        DO J=1,NNYP
        DO I=1,NNXP
                GEODAT2D(I,J)=GEODAT2D(I,J)/100.
        ENDDO
        ENDDO
                                                                                
        istatus_tmp=1 !bogus
                                                                                
        ENDIF
	
!       do J=NNYP,1,-NNYP/25
!         write(6,833) (GEODAT2D(I,J),I=1,NNXP,NNXP/15)
!       enddo

        if(istatus_tmp.ne.1)then
         print* 
         print*,'soiltemp data not processed completely' 
         if(c10_grid_f(1:lf).eq.'wrfsi')then 
            print*,' Error:  File(s) missing for soiltemp data' 
            print*,' Error:  Static file not created'
            print*
            istatus=0
            return
         else
            print*,' File(s) missing for soiltemp data'
            print*,'           *** WARNING ***'
            print*,' soiltemp data not added to static file'
            print*
         endif
         GEODAT2D=r_missing_data
        endif

        idst=25
        if(c10_grid_f(1:lf).eq.'wrfsi')then 

            IF (c6_maproj .ne. 'rotlat') THEN

                idst=96

                call adjust_geog(nnxp-1,nnyp-1,1,'soiltemp',istatus_tmp
     &,lats(1:nnxp-1,1:nnyp-1,ns),data(1:nnxp-1,1:nnyp-1,indxter)
     &,data(1:nnxp-1,1:nnyp-1,ilndmsk),GEODAT2D(1:nnxp-1,1:nnyp-1)
     &,istatus)

            ELSE
                idst=87
            ENDIF


        else

            call adjust_geog(nnxp,nnyp,1,'soiltemp',istatus_grn
     &,lats(:,:,ns),data(:,:,indxter),data(:,:,ilndmsk),GEODAT2D
     &,istatus)

        endif 
	write(6,*) 'idst= ', idst
        data(:,:,idst)=GEODAT2D
	write(6,*) 'idst written'
c
c -------------- Terrain slope index categories -------------------------
c type = I
c
        itstatus=ishow_timer()
        print*
        print*,' Processing 1 deg islope data'

        allocate (adum3d(nnxp,nnyp,9))  !temporarily holds %dist of islope cat

	IF (c6_maproj .ne. 'rotlat') THEN

        CALL GEODAT(nnxp,nnyp,erad,90.,std_lon,xtn(1,ns)
     +,ytn(1,ns),deltax,deltay,GEODAT2D,adum3d,adum2d,adum2d
     +,PATH_TO_ISLOPE,2.0,0.0,new_DEM,9,istatus_slp)

	ELSE

	categorical=.false.
!!	should be done as n.n. interp by alt_interp code
	useland=.FALSE.
        call alt_hemi_all(nnxp,nnyp,hlat,hlon,
     &                       PATH_TO_ISLOPE,
     &                       type,1,grid_spacing_m/1000.,ADUM3D,
     &                       GEODAT2D,categorical,data(:,:,ilndmsk),
     &			     useland)

	istatus_slp=1 ! bogus

	ENDIF

        deallocate (adum3d)

        if(istatus_slp.ne.1)then
           print*
           print*,' ter slope category data not processed completely'
           if(c10_grid_f(1:lf).eq.'wrfsi')then
             print*,' Error: File(s) missing for islope data'
             print*,' Error: Static file not created'
             print*
             return
           else
             print*,' Warning: File(s) missing for islope data'
             print*,' Warning: islope not added to static file'
             print*
           endif
           GEODAT2D=r_missing_data
           GEODAT3D=r_missing_data
        endif

c put the categories back to the original raw data. if it is a land
c point but islope indicates water, force islope = 1.
        where(GEODAT2D .eq. 8)GEODAT2D=13.0
        where(GEODAT2D .eq. 9)GEODAT2D=0.0

        if(c10_grid_f(1:lf).eq.'wrfsi')then
	
	   IF (c6_maproj .ne. 'rotlat') THEN

               islp=109
               call adjust_geog(nnxp-1,nnyp-1,1,'islope',istatus_slp
     &,lats(1:nnxp-1,1:nnyp-1,ns),data(1:nnxp-1,1:nnyp-1,indxter)
     &,data(1:nnxp-1,1:nnyp-1,ilndmsk),GEODAT2D(1:nnxp-1,1:nnyp-1)
     &,istatus)

           ELSE
               islp=100 ! will need more adjustments down the road for this
           ENDIF

        else

           islp=38
           call adjust_geog(nnxp,nnyp,1,'islope',istatus_slp
     &,lats(1,1,ns),data(1,1,indxter),data(1,1,ilndmsk),GEODAT2D
     &,istatus)

        endif

        if(istatus.ne.1)then
           print*,' Warning: Processing incomplete: adjust_geog_data'
           if(c10_grid_f(1:lf).eq.'wrfsi')then
              print*,'Error: wrfsi static file not generated'
              return
           endif
        endif

        data(:,:,islp)=GEODAT2D

c force land points to have the correct (default) islope value.
        where(data(:,:,islp)   .eq. 0.0 .and.
     .        data(:,:,ilndmsk).eq. 1.0)data(:,:,islp)=1.0
c force water points to have the correct category for islope
        where(data(:,:,ilndmsk).eq. 0.0)data(:,:,islp)=0.0

!write(6,*) 'ISLOPE: '
!      do J=nnyp,1,-nnyp/30
!       write(6,733) (int(data(I,J,islp)),I=1,nnxp,nnxp/20)
!       enddo
c
c    subroutine adjust_geog_data now named adjust_geog and called
c    after each data type is processed (as opposed to all geog data
c    getting adjusted with one call).
c
c!JS: even though Matt Pyle suggests that adjust geog is fine for nmm
c     it currently isn't used.  If it should be used, then it would be
c     similar to the LAPS case which is not a staggered grid.
c
!mp	believe adjust_geog_data should be fine for nmm case
c
c -------------------------------------------------------------------
c ---------------Monthly Albedo Climatology--------------------------
c type = A
c
        itstatus=ishow_timer()
        print*
        print*,' Processing albedo climo data'

	IF (c6_maproj .ne. 'rotlat') THEN

        call proc_geodat(nnxp,nnyp,12,path_to_albedo
     +,lats(1,1,ns),lons(1,1,ns),data(1,1,ilndmsk)
     +,GEODAT3D,istatus_alb)

	ELSE

	categorical=.FALSE.
	useland=.TRUE.
        call alt_hemi_all(nnxp,nnyp,hlat,hlon,
     &                       PATH_TO_ALBEDO,
     &                       type,12,grid_spacing_m/1000.,GEODAT3D,
     &                       ADUM2D,categorical,data(1,1,ilndmsk),
     &			     useland)

        istatus_alb=1 !bogus

	do N=1,12
	 do J=1,NNYP
	  do I=1,NNXP
	   GEODAT3D(I,J,N)=GEODAT3D(I,J,N)*.01
	  enddo
	 enddo
	enddo

        ENDIF

        print*,'Done in proc_geodat: albedo'

        if(istatus_alb.ne.1)then
         print*
         print*,' Error: Albedo climo data not processed completely'
         print*,' Error: File(s) missing for albedo data'
         print*,' Error: Static file not created: albedo missing'
         print*
         istatus=0
         return
        endif

c force water points to 0.08 albedo

        do k=1,12
           where(data(:,:,ilndmsk).eq.0.0)GEODAT3D(:,:,k)=0.08
        enddo

        if(c10_grid_f(1:lf).eq.'wrfsi')then
        IF (c6_maproj .ne. 'rotlat') THEN
           ialb=96
	ELSE
	   ialb=87
	ENDIF
        else
           ialb=25
        endif
        do j=1,12
           data(:,:,ialb+j)=GEODAT3D(:,:,j)
        enddo
c
c ---------------- Max Snow Albedo ------------------
c type = M
c
	IF (c6_maproj .ne. 'rotlat') THEN

        call proc_geodat(nnxp,nnyp,1,path_to_maxsnoalb
     +,lats(1,1,ns),lons(1,1,ns),data(1,1,ilndmsk)
     +,GEODAT2D,istatus_alb)

	ELSE

	categorical=.false.

	useland=.TRUE.
	allocate(adum3d(nnxp,nnyp,1))
        call alt_hemi_all(nnxp,nnyp,hlat,hlon,
     &                       PATH_TO_MAXSNOALB,
     &                       type,1,grid_spacing_m/1000.,ADUM3D,
     &                       GEODAT2D,categorical,data(1,1,ilndmsk),
     &			     useland)

	do J=1,NNYP
	do I=1,NNXP
	GEODAT2D(I,J)=GEODAT2D(I,J)*.01
	enddo
	enddo

	ENDIF

        print*,'Done in proc_geodat: Max Snow Albedo'

        if(istatus_alb.ne.1)then
         print*
         if(c10_grid_f(1:lf).eq.'wrfsi')then
            print*,'--------------- WRFSI ------------------'
            print*,'Error: Max Snow Albedo not processed completely'
            print*,'Error: Static file not created '
            print*
            istatus=0
            return
         else
            print*,'Warning: Max Snow Albedo not processed completely'
         endif
        endif
c
c force max albedo = 0.08 over water. force max albedo = 0.7 over ice
c
        where(data(:,:,ilndmsk).eq.0.0)GEODAT2D=0.08
        where(data(:,:,ilnduse).eq.24.0)GEODAT2D=0.7

        if(c10_grid_f(1:lf).eq.'wrfsi')then

	IF (c6_maproj .ne. 'rotlat') THEN
           mxsnalb=47
           data(:,:,mxsnalb)=GEODAT2D
	ELSE
           mxsnalb=14
	   data(:,:,mxsnalb)=GEODAT2D	
	ENDIF

        else
           mxsnalb=6
           data(:,:,mxsnalb)=GEODAT2D
        endif

        deallocate (GEODAT2D)
        deallocate (GEODAT3D)
c
c ---------------------------------------------------------------------------------
c Let's compare the grids to landmask to assess their consistency (or lack thereof)
c ---------------------------------------------------------------------------------

       if(c10_grid_f(1:lf).eq.'wrfsi')then
        print*,'Total number of gridpoints (nx*ny) = ',nnxp*nnyp
        print*
        do i=1,10
           if(i == 1)ii=indxter  !terrain
           if(i == 2)ii=idxstl   !soil top layer
           if(i == 3)ii=idxsbl   !soil bot layer
           if(i == 4)ii=igrn+6   !mxgiwrf  !max greenness
           if(i == 5)ii=mngiwrf  !min greenness
           if(i == 6)ii=idst     !deep soil temp
           if(i == 7)ii=islp     !terrain slope index
           if(i == 8)ii=ialb+6   !albedo; arbitrarily at month 6
           if(i == 9)ii=mxsnalb  !max snow albedo
           if(i == 10)ii=ilnduse !landuse dominant category
           call gridcompare(nnxp,nnyp,i,data(1,1,ii),data(1,1,ilndmsk)
     &,istatus)
        enddo
       endif
c
c ---------------------------------------------------------------------------------
c --------- This is were we prepare for writing the netCDF static file ------------
c ---------------------------------------------------------------------------------

        if(.not.allocated(var))allocate (var(ngrids)
     +,comment(ngrids))

        if(c10_grid_f(1:lf).eq.'wrfsi')then

	IF (c6_maproj .ne. 'rotlat') THEN

           in1=0
           in2=0
           do j=1,ns
              in1=in2+1
              in2=in1+1
              data(:,:,in1)=lats(:,:,j)
              data(:,:,in2)=lons(:,:,j)
           enddo


c compute map/grid information. Reallocate GEODAT3D array
           allocate (GEODAT3D(nnxp,nnyp,4))
           call get_projrot_grid(nnxp,nnyp,lats(1,1,ns)
     +,lons(1,1,ns),GEODAT3D,istatus)

           i=39

           data(:,:,i)  =GEODAT3D(:,:,1)
           data(:,:,i+1)=GEODAT3D(:,:,2)
 
           call get_map_factor_grid(nnxp,nnyp,n_staggers
     +,lats,lons ,GEODAT3D,istatus)
           if(istatus.ne.1)then
              print*,'Error returned: get_maps_factor_grid'
              return
           endif

           do j=1,ns
              data(:,:,i+j+1)=GEODAT3D(:,:,j)
           enddo
c           
           call get_coriolis_components(nnxp,nnyp,lats(1,1,ns)
     +,GEODAT3D)
           data(:,:,i+6)=GEODAT3D(:,:,1)
           data(:,:,i+7)=GEODAT3D(:,:,2)

           deallocate(GEODAT3D)

c          call move(static_albedo,data(1,1,i+8),nnxp,nnyp)

           call move(topt_out_s,data(1,1,i+9),nnxp,nnyp)
           call move(topt_out_ln,data(1,1,i+10),nnxp,nnyp)
           call move(topt_out_lt,data(1,1,i+11),nnxp,nnyp)
c          call move(topt_stag_out,data(1,1,i+12),nnxp,nnyp) !51

           call get_gridgen_var(ngrids,ngrids,var,comment)

	ELSE ! rotlat case
                                                                                

!!!!!	these doloops possibly redundant.  Defined above????

!        DO J=1,NNYP
!        DO I=1,NNXP
!              data(I,J,1)=hlat(I,J)
!              data(I,J,2)=hlon(I,J)
!              data(I,J,3)=vlat(I,J)
!              data(I,J,4)=vlon(I,J)
!        ENDDO
!        ENDDO
                                                                                
           data(:,:,18)=topt_out
                                                                                
c compute map/grid information. Reallocate GEODAT3D array
           allocate (GEODAT3D(nnxp,nnyp,4))
                                                                                
                                                                                
           call vecrot_rotlat(nnxp,nnyp,mdlat,mdlon
     +,vlat,vlon,GEODAT3D(:,:,1),GEODAT3D(:,:,2))

!        write(6,*) 'cos term'
        do J=nnyp,1,-nnyp/40
!        write(6,833) (GEODAT3D(I,J,1),I=1,nnxp,nnxp/15)
        enddo
!        write(6,*) 'sin term'
        do J=nnyp,1,-nnyp/40
!        write(6,833) (GEODAT3D(I,J,2),I=1,nnxp,nnxp/15)
        enddo

  833	format(30(f5.1,1x))

           i=10
                                                                                
!       sin, then cosine?
           data(:,:,i)  =GEODAT3D(:,:,2)
           data(:,:,i+1)=GEODAT3D(:,:,1)
c
                                                                                
!!!!    coriolis just dependent on lat (hopefully)
           call get_coriolis_components(nnxp,nnyp,vlat
     +,GEODAT3D)
           data(:,:,i+2)=GEODAT3D(:,:,1)
           data(:,:,i+3)=GEODAT3D(:,:,2)

           deallocate(GEODAT3D)

c          call move(static_albedo,data(1,1,i+4),nnxp,nnyp)
                                                                                
           call move(topt_out_s,data(1,1,i+5),nnxp,nnyp)
           call move(topt_out_ln,data(1,1,i+6),nnxp,nnyp)
           call move(topt_out_lt,data(1,1,i+7),nnxp,nnyp)
c          call move(topt_stag_out,data(1,1,i+8),nnxp,nnyp) !18
                                                                                
!030513  Made modifications within the get_gridgen_var code
!	to handle NMM variable set
                                                                                
        write(6,*) 'calling get_gridgen_var with ngrids= ', ngrids
           call get_gridgen_var(ngrids,ngrids,var,comment)

        ENDIF

        else ! non-wrfsi

           call move(lats(1,1,1),data(1,1,1),nnxp,nnyp)            ! KWD
           call move(lons(1,1,1),data(1,1,2),nnxp,nnyp)            ! KWD
           call move(topt_out,data(1,1,3),nnxp,nnyp)               ! KWD
           call move(topt_out_s,data(1,1,7),nnxp,nnyp)             ! JS
           call move(topt_out_ln,data(1,1,8),nnxp,nnyp)            ! JS
           call move(topt_out_lt,data(1,1,9),nnxp,nnyp)            ! JS 

           call get_gridgen_var(ngrids,ngrids,var,comment)
 
        endif

!mp
        if (c6_maproj .ne. 'rotlat') then
          if(c10_grid_f(1:lf).eq.'wrfsi')then
             filename=c10_grid_f(1:lf)//'.d'//cnest//'.cdl'
          else
             filename = c10_grid_f(1:lf)//'.cdl'
          endif
        else
          filename = c10_grid_f(1:lf)//'.rotlat.cdl'
        endif

	write(6,*) 'using filename: ', filename
        call s_len(filename,lfn)
        call get_directory('cdl',cdl_dir,lcdl)

	write(6,*) 'file= ', cdl_dir(1:lcdl)//filename(1:lfn)
        INQUIRE(FILE=cdl_dir(1:lcdl)//filename(1:lfn),EXIST=exist)

        if(.not.exist)then
           print*,'Error: Could not find file '
     +           ,cdl_dir(1:lcdl)//filename(1:lfn)
           print*,'c10_grid_f: ',c10_grid_f(1:lf)
           istatus = 0
           return
        endif
	
        IF (c6_maproj .ne. 'rotlat') THEN

        call check_domain(lats(1,1,ns),lons(1,1,ns)
     +,nnxp,nnyp,grid_spacing_m,1,istat_chk)

        write(6,*)'deltax = ',deltax

        if(istat_chk .eq. 1)then
            write(6,*)'check_domain: status = ',istat_chk
        else
            write(6,*)'ERROR in check_domain: status = ',istat_chk       
        endif

	ENDIF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	write(6,*) 'call put_laps_static ', ngrids

	write(6,*) 'model= ', model

        c10_grid_fname = 'nest7grid'
        call s_len(c10_grid_fname,len_fname)

!       Do zin calc (note this is last [kmax] element in data array)
        if(c10_grid_fname(1:len_fname).eq.'nest7grid')then
          avgelem=3
          zinelem=ngrids
          write(6,*)' Calculating zin field in array element ',zinelem
          do i = 1,nnxp
          do j = 1,nnyp
              psa = ztopsa(data(i,j,avgelem)) ! This is the AVG data (3rd element)
              data(i,j,zinelem) = (20.0 - ((psa - 100.0) * 0.02))
          enddo ! j
          enddo ! i
        endif

        call put_laps_static(grid_spacing_m,model,comment,var
     1,data,nnxp,nnyp,ngrids,ngrids,std_lat,std_lat2,std_lon
     1,c6_maproj,deltax,deltay,mdlat,mdlon,lli,llj,iuri,iurj
     1,iparent_id,iratio_to_parent,nest,c10_grid_fname)

        istatus = istat_chk

        itstatus=ishow_timer()
        print*,' -------------------------------'
        print*,' Total Elapsed time: ', itstatus
        print*,' -------------------------------'

        deallocate (data)

	return
	End
c
c ================================================================
c
      SUBROUTINE GEODAT(n2,n3,erad,rlat,wlon1,xt,yt,deltax,deltay
     1 ,DATR,DATS,DATLN,DATLT,OFN,WVLN,SILWT,which_data,maxdatacat
     1 ,istat_files)

      include 'trigd.inc'

      implicit none


      integer n2,n3
      integer lbf
      integer lb,mof,np,niq,njq,nx,ny,isbego,iwbego,
     1  iblksizo,no,iodim,istat_files,maxdatacat

      real erad,rlat,wlon1,deltax,deltay,wvln,silwt

c wvln = TOPTWVL_PARM_WRF fror wrfsi.nl section &hgridspec
c siwt = SILAVWT_PARM_WRF      "                

      real DATR(N2,N3)
      real DATS(N2,N3,maxdatacat)
      real DATLN(N2,N3)
      real DATLT(N2,N3)

c this parameter can be increased if we ever read data with
c resolution finer than 30 sec, or if the tilesize for 30s
c data becomes greater than 10x10 deg.

      PARAMETER(IODIM=5800000)

      real xt(N2),YT(N3)
      real deltallo,deltaxq,deltayq,
     1  deltaxp,deltayp

      real RSOFF,RWOFF

      real std_lon
      integer istatus
      integer lcat

      CHARACTER*(*) OFN
      character*180 TITLE
      logical which_data
C
c *********************
      nx=n2-1
      ny=n3-1
c ****************************

      lcat = 1
      LBF=INDEX(OFN,' ')-1
      TITLE=OFN(1:LBF)//'HEADER'


      LB=INDEX(TITLE,' ')-1

      CALL JCLGET(29,TITLE(1:LB),'FORMATTED',1,istat_files)
      if(istat_files .ne. 1)then
         write(6,*)' Warning in gridgen_model opening HEADER: '
     1,' geog path = ', title(1:lb)
         return
      endif

      READ(29,*)IBLKSIZO,NO,ISBEGO,IWBEGO,RSOFF,RWOFF
c 2    FORMAT(4I5,2(F10.8))
      print *,'title=',title
      print *,'rsoff,rwoff = ',rsoff,rwoff
      print *,'isbego,iwbego =',isbego,iwbego
      print *,'iblksizo,no =',iblksizo,no
      CLOSE(29)

      if(NO .gt. 1)then
         DELTALLO=FLOAT(IBLKSIZO)/FLOAT(NO-1)
      elseif(NO .eq. 1)then
         DELTALLO=FLOAT(IBLKSIZO)/FLOAT(NO)
      else
         print*,'HEADER value NO = 0'
         return
         istat_files = 0
      endif

      MOF=IODIM/(NO*NO)

      if(ofn(lbf:lbf).eq.'G'.or.ofn(lbf:lbf).eq.'A')then
         lcat=12
         if(no.eq.1250)then
            MOF=1
         endif
      endif
c SG97 MOF determines the number of files held in buffer while reading
c SG97 DEM data; it saves some time when buffer data can be used instead
c SG97 of reading DEM file again. Originally MOF was 4.
      if (MOF.gt.10) MOF=5
      
      DELTAXQ=0.5*WVLN*DELTAX
      DELTAYQ=0.5*WVLN*DELTAY
      print *,'deltaxq,deltayq=',deltaxq,deltayq
      NP=MIN(10,MAX(1,INT(DELTAXQ/(DELTALLO*111000.))))
      print *,' np=',np
      DELTAXP=DELTAXQ/FLOAT(NP)
      DELTAYP=DELTAYQ/FLOAT(NP)
      NIQ=INT(FLOAT(NX)*DELTAX/DELTAXQ)+4
      NJQ=INT(FLOAT(NY)*DELTAY/DELTAYQ)+4
C
      call get_standard_longitude(std_lon,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Error calling laps routine'
          stop 
      endif

      CALL SFCOPQR(NO,MOF,NP,NIQ,NJQ,N2,N3,lcat
     +,XT,YT,90.,std_lon,ERAD,RWOFF,RSOFF
     +,DELTALLO,DELTAXP,DELTAYP,DELTAXQ,DELTAYQ
     +,IBLKSIZO,ISBEGO,IWBEGO,DATR,DATS,DATLN,DATLT
     +,OFN,WVLN,SILWT,which_data,maxdatacat,istat_files)       
      RETURN
      END


C
C     ******************************************************************
C
      SUBROUTINE SFCOPQR(NO,MOF,NP,NIQ,NJQ,N2,N3,lcat
     +          ,XT,YT,RLAT,WLON1,ERAD,RWOFF,RSOFF
     +          ,DELTALLO,DELTAXP,DELTAYP,DELTAXQ,DELTAYQ
     +          ,IBLKSIZO,ISBEGO,IWBEGO,DATR,DATS,DATLN,DATLT
     +          ,OFN,WVLN,SILWT,dem_data,maxdatacat,istat_files)

c JS: removed dato array from subroutine argument list
c JS: added RWOFF/RSOFF - West and South offset of tile data

      real,  allocatable ::  dato(:,:,:,:)    !dato(no,no,mof,lcat)

      real,  allocatable ::  DATP(:,:,:)
      real,  allocatable ::  DATQ(:,:)
      real,  allocatable ::  DATQS(:,:,:)
      real,  allocatable ::  DATSM(:,:)
      real,  allocatable ::  DATSMX(:,:)
      real,  allocatable ::  DATSLN(:,:)
      real,  allocatable ::  DATSLT(:,:)

      real,  intent(inout) ::  DATLN(N2,N3)
      real,  intent(out)   ::  DATLT(N2,N3)
      real,  intent(out)   ::  DATR(N2,N3)
      real,  intent(out)   ::  DATS(N2,N3,maxdatacat)

      real ISO(MOF),IWO(MOF),XT(N2),YT(N3),rlat,wlon1,
     +     erad,deltallo,deltaxp,deltayp,deltaxq,deltayq,
     +     wvln,silwt,xq,yq,xp,yp,xcentr,ycentr,glatp,               ! pla,plo,
     +     glonp,rio,rjo,wio2,wio1,wjo2,wjo1,xq1,yq1,
     +     rwoff,rsoff

      real r_missing_data
      real xr,yr,rval,sh,sha,rh,rha,rhn,rht,shn,sht
      real shln,shlt,rhln,rhlt
      real delta_ln(np,np),delta_lt(np,np)

c     real xpmn,xpmx,ypmn,ypmx
      real xp1,xp2,yp1,yp2
      real xpcentr,ypcentr

      real  pctcat(maxdatacat)

      integer lp
      integer ixr,iyr
      integer lent

      CHARACTER*180 OFN
      CHARACTER*180 TITLE3,TITLE3_last_read,TITLE3_last_inquire
      CHARACTER*3   TITLE1
      CHARACTER*4   TITLE2
      CHARACTER*10  cdatatype

      LOGICAL L1,L2,dem_data,l_string_contains

      data icnt/0/
      save icnt
C
      print *,'no,mof,np,niq,njq=',no,mof,np,niq,njq

      istat_files = 1

      NONO=NO*NO
      XCENTR=0.5*(XT(1)+XT(N2))
      YCENTR=0.5*(YT(1)+YT(N3))
      print *,xt(1),xt(n2),xcentr
      print *,'deltaxp=',deltaxp
      NOFR=0
      DO 11 IOF=1,MOF
         ISO(IOF)=0
         IWO(IOF)=0
  11  continue

      TITLE3_last_read    = '/dev/null'
      TITLE3_last_inquire = '/dev/null'

      lcat=1
      len=index(ofn,' ')
      if(ofn(len-1:len-1).eq.'V'.or.
     &   ofn(len-1:len-1).eq.'I')then
         icnt = 0
         cdatatype='landuse'
         if(ofn(len-1:len-1).eq.'I')cdatatype='islope'
      elseif(ofn(len-1:len-1).eq.'O')then
         icnt = 0
         cdatatype='soiltype'
      elseif(ofn(len-1:len-1).eq.'U' .or.
     &       ofn(len-1:len-1).eq.'H')then
         cdatatype='topography'
      endif

      print*,'SFCOPQR: cdatatype = ',cdatatype

      call s_len(cdatatype,lent)

      allocate(dato(no,no,mof,lcat))

      allocate (DATP(NP,NP,lcat),
     &          DATQ(NIQ,NJQ),
     &          DATSM(NIQ,NJQ),
     &          DATSMX(NIQ,NJQ),
     &          DATSLN(NIQ,NJQ),
     &          DATSLT(NIQ,NJQ),
     &          DATQS(NIQ,NJQ,maxdatacat))

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus.ne.1)then
         print*,'failed to get r_missing_data'
         return
      endif

      DO 15 JQ=1,NJQ
         print *,'jq,njq,niq,nofr=',jq,njq,niq,nofr
         DO 16 IQ=1,NIQ

            XQ=(FLOAT(IQ)-0.5*FLOAT(NIQ+1))*DELTAXQ+XCENTR
            YQ=(FLOAT(JQ)-0.5*FLOAT(NJQ+1))*DELTAYQ+YCENTR

            xpmn=1.0e30
            ypmn=1.0e30
c           xpmx=-1.0e30
c           ypmx=-1.0e30

            DO 17 JP=1,NP
               DO 18 IP=1,NP

                  XP=XQ+(FLOAT(IP)-0.5*FLOAT(NP+1))*DELTAXP
                  YP=YQ+(FLOAT(JP)-0.5*FLOAT(NP+1))*DELTAYP

                  call xy_to_latlon(XP,YP,erad,GLATP,GLONP) 

                  glatp = max(-89.9999,min(89.9999,glatp - rsoff))
                  glonp = glonp - rwoff
                  if(glonp.ge.180.) glonp = glonp - 360.
                  if(glonp.le.-180.) glonp = glonp + 360.

c                 print *,'rlat,wlon1=',rlat,wlon1

                  ISOC=(INT((GLATP-FLOAT(ISBEGO))/FLOAT(IBLKSIZO)
     &          +200.)-200)*IBLKSIZO+ISBEGO
            IWOC=(INT((GLONP-FLOAT(IWBEGO))/FLOAT(IBLKSIZO)
     &          +400.)-400)*IBLKSIZO+IWBEGO

                  DO 19 IOFR=1,NOFR
                     JOFR=IOFR
                     IF(ISO(IOFR).EQ.ISOC.AND.IWO(IOFR).EQ.IWOC)GO TO 10
 19                 continue
                  ISOCPT=ABS(ISOC)/10
                  ISOCPO=ABS(ISOC)-ISOCPT*10
                  IWOCPH=ABS(IWOC)/100
                  IWOCPT=(ABS(IWOC)-IWOCPH*100)/10
                  IWOCPO=ABS(IWOC)-IWOCPH*100-IWOCPT*10
                  IF(ISOC.GE.0)THEN
                     WRITE(TITLE1,'(2I1,A1)')ISOCPT,ISOCPO,'N'
                  ELSE
                     WRITE(TITLE1,'(2I1,A1)')ISOCPT,ISOCPO,'S'
                  ENDIF

                  IF(IWOC.GE.0 
     1               .and. IWOC .ne. 180                    ! 1998 Steve Albers
     1                                      )THEN
                     WRITE(TITLE2,'(3I1,A1)')IWOCPH,IWOCPT,IWOCPO,'E'
                  ELSE
                     WRITE(TITLE2,'(3I1,A1)')IWOCPH,IWOCPT,IWOCPO,'W'
                  ENDIF

                  LB=INDEX(OFN,' ')-1
                  TITLE3=OFN(1:LB)//TITLE1//TITLE2
                  LB=INDEX(TITLE3,' ')-1

                  if(TITLE3 .ne. TITLE3_last_inquire)then
                     INQUIRE(FILE=TITLE3(1:LB),EXIST=L1,OPENED=L2)
                     TITLE3_last_inquire = TITLE3
                  endif

                  IF(.NOT.L1)THEN
                     iwrite = 0

                     if(icnt .le. 100)then ! Reduce the output
                         iwrite=1

                     elseif(icnt .le. 1000)then
                         if(icnt .eq. (icnt/100)*100)iwrite=1

                     elseif(icnt .le. 10000)then
                         if(icnt .eq. (icnt/1000)*1000)iwrite=1

                     elseif(icnt .le. 100000)then
                         if(icnt .eq. (icnt/10000)*10000)iwrite=1

                     else
                         if(icnt .eq. (icnt/100000)*100000)iwrite=1

                     endif

                     if(iwrite .eq. 1)then
                        if(l_string_contains(TITLE3(1:LB),
     1                                       'world_topo_30s',
     1                                       istatus)             )then       
                           PRINT*, ' ERROR: ',TITLE3(1:LB)
     1                            ,' DOES NOT EXIST ',icnt

                        elseif(l_string_contains(TITLE3(1:LB),
     1                                           'topo_30s',
     1                                       istatus)             )then       
                           PRINT*, ' topo_30s file ',TITLE3(1:LB)
     1                            ,' does not exist, using topo_10m '
     1                            ,icnt

                        else ! Generic warning message
                           PRINT*, ' WARNING: ',TITLE3(1:LB)
     1                            ,' DOES NOT EXIST ',icnt

                        endif

                     endif ! iwrite

                     icnt = icnt + 1

c initialize these arrays as they may have some garbage in them
c if we don't actually read in any data.
c
                     DATP(IP,JP,:) = 0.
                     DELTA_LN(IP,JP) = 0.
                     DELTA_LT(IP,JP) = 0.
                     istat_files = 0
                     GO TO 20

                  ENDIF

                  IF(NOFR.GE.MOF)THEN
                     DO 21 IOF=1,MOF
                        ISO(IOF)=0
                        IWO(IOF)=0
21                    continue
                     NOFR=0
                  ENDIF
                  NOFR=NOFR+1
                  JOFR=NOFR

!                 Read the tile
                  if(TITLE3 .ne. TITLE3_last_read)then
                    if( (ofn(len-1:len).eq.'U').and.(no.eq.1200).or.
     .                   no.eq.1201 )then
                         if(no.eq.1201)then
                            print*,'Reading ', title3(1:lb)
                            CALL READ_DEM(29,TITLE3(1:LB),no,no,4,4, ! topo_3s experimental
     .                              DATO(1,1,NOFR,1),istat)
                         else
                            CALL READ_DEM(29,TITLE3(1:LB),no,no,2,2, ! world topo_30s
     .                              DATO(1,1,NOFR,1),istat)
                         endif
                      dem_data=.true.
                    elseif( (ofn(len-1:len).eq.'O') )then      ! soiltype top and bot layer
                      CALL READ_DEM(29,TITLE3(1:LB),no,no,1,4,
     .                              DATO(1,1,NOFR,1),istat)
                      dem_data=.true.
                    elseif( (ofn(len-1:len).eq.'V') )then      ! world USGS 30s landuse
                      CALL READ_DEM(29,TITLE3(1:LB),no,no,1,4,
     .                              DATO(1,1,NOFR,1),istat)
                      dem_data=.true.
                    elseif((ofn(len-1:len).eq.'I'))then      ! only islope in this code section
                      CALL READ_DEM_G(29,TITLE3(1:LB),no,no,1,lcat
     .                     ,nofr, 1,4, DATO,istat)
                      dem_data=.true.
                    elseif( (ofn(len-1:len).eq.'T') )then      ! soiltemp - obsolete in this code
                      CALL READ_DEM(29,TITLE3(1:LB),no,no,2,2,
     .                              DATO(1,1,NOFR,1),istat)
                      dem_data=.true.
                    else                                       ! other
                      CALL JCLGET(29,TITLE3(1:LB),'FORMATTED',0,istatus)      
                      CALL VFIREC(29,DATO(1,1,NOFR,1),NONO,'LIN')
                      if ((ofn(len-1:len).eq.'U').and.(no.eq.121)) then
                        dem_data=.false.                       ! topo_30s
                      endif
                    endif

                    if(istat.ne.0)then
                       print*,'Error returned: SFCOPQR: READ_DEM'
                       return
                    endif


                    TITLE3_last_read = TITLE3

c                   print *,'nofr,dato=',nofr,dato(1,1,nofr)
                    CLOSE(29)

                  else
                    write(6,*)' We have made the code more efficient'

                  endif ! Is this a new file we haven't read yet?

                  ISO(NOFR)=ISOC
                  IWO(NOFR)=IWOC
10		  continue

                  RIO=(GLONP-FLOAT(IWOC))/DELTALLO+1.
                  RJO=(GLATP-FLOAT(ISOC))/DELTALLO+1.

!                 Prevent Bounds Error (Steve Albers)
                  if(RIO .lt. 1.0)then
                      if(RIO .gt. 0.98)then
                          write(6,*)' Reset RIO for Machine Epsilon'      
                          RIO = 1.0
                      elseif(RIO .lt. 0.5)then
                          write(6,*)' ERROR: RIO out of bounds',RIO
                          stop
                      endif
                  endif

                  if(RJO .lt. 1.0)then
                      if(RJO .gt. 0.98)then
                          write(6,*)' Reset RJO for Machine Epsilon'      
                          write(6,*)JQ,IQ,
     1                          IP,JP,IO1,JO1,JOFR,RIO,RJO,GLATP,ISOC
                          RJO = 1.0
                      elseif(RJO .lt. 0.5)then
                          write(6,*)' ERROR: RJO out of bounds',RJO
                          write(6,*)JQ,IQ,
     1                          IP,JP,IO1,JO1,JOFR,RIO,RJO,GLATP,ISOC
                          stop
                      endif
                  endif

C Interp OK for continuous data such as topography

                  if(cdatatype.eq.'topography')then

                   IO1=INT(RIO)
                   JO1=INT(RJO)
                   IO2=IO1+1
                   JO2=JO1+1
                   WIO2=RIO-FLOAT(IO1)
                   WJO2=RJO-FLOAT(JO1)
                   WIO1=1.0-WIO2
                   WJO1=1.0-WJO2

                   do LP = 1,lcat

                   DATP(IP,JP,LP)=WIO1*(WJO1*DATO(IO1,JO1,JOFR,LP)
     +                                 +WJO2*DATO(IO1,JO2,JOFR,LP))
     +                           +WIO2*(WJO1*DATO(IO2,JO1,JOFR,LP)
     +                                 +WJO2*DATO(IO2,JO2,JOFR,LP))

!S & W-facing slopes > 0.
                   DELTA_LN(IP,JP)=
     .           ((DATO(IO2,JO1,JOFR,LP)-DATO(IO1,JO1,JOFR,LP))+
     .            (DATO(IO2,JO2,JOFR,LP)-DATO(IO1,JO2,JOFR,LP)))*.5

                   DELTA_LT(IP,JP)=
     .           ((DATO(IO1,JO2,JOFR,LP)-DATO(IO1,JO1,JOFR,LP))+
     .            (DATO(IO2,JO2,JOFR,LP)-DATO(IO2,JO1,JOFR,LP)))*.5

                   enddo !LP = 1,lcat

                  else

C Nearest grid point for landuse and soiltype

                   IO1=NINT(RIO)
                   JO1=NINT(RJO)
                   do LP = 1,lcat
                    DATP(IP,JP,LP)= DATO(IO1,JO1,JOFR,LP)
                   enddo

                  endif ! cdatatype eq topography.
                   
20               CONTINUE
18             continue ! IP
17           continue ! JP


!           print*,'xpmx/xpmn//ypmx/ypmn/ ',xpmx,xpmn,ypmx,ypmn

! Calculate average and silhouette terrain, then apply SILWT weight

            if(cdatatype(1:lent).eq.'topography')then

             SHA=0.
             RHA=0.
             RHLN=0.
             RHLT=0.
             shmax=0.

             DO 22 JP=1,NP
               SH=0.
               RH=0.
               RHN=0.
               RHT=0.
               DO 23 IP=1,NP
!                 Test for missing - then go to 16?
                  SH=max(SH,DATP(IP,JP,1)) 
                  RH=RH+DATP(IP,JP,1)
                  RHN=RHN+DELTA_LN(IP,JP)
                  RHT=RHT+DELTA_LT(IP,JP)
23             continue ! IP
               SHA=SHA+SH/(2.*FLOAT(NP))
               RHA=RHA+RH
               RHLN=RHLN+RHN
               RHLT=RHLT+RHT
               SHMAX=max(SHMAX,SH)
22           continue ! JP
 
             RHA=RHA/FLOAT(NP*NP)
             RMS=0.0
             DO 24 IP=1,NP 
               SH=0.
               DO 25 JP=1,NP
                  SH=max(SH,DATP(IP,JP,1))
                  RMS=RMS+((DATP(IP,JP,1)-RHA)*(DATP(IP,JP,1)-RHA))
25             continue ! JP
               SHA=SHA+SH/(2.*FLOAT(NP))
24           continue ! IP

             DATQS(IQ,JQ,1)=SQRT(RMS/FLOAT(NP*NP))
             DATQ(IQ,JQ)=SHA*SILWT+RHA*(1.-SILWT)
             DATSM(IQ,JQ)=RHA                           !mean value of points used for IQ,JQ
             DATSMX(IQ,JQ)=SHMAX                        !max value from points used for IQ,JQ
             DATSLN(IQ,JQ)=RHLN/FLOAT(NP*NP)/DELTAXP
             DATSLT(IQ,JQ)=RHLT/FLOAT(NP*NP)/DELTAYP

c            print *,'datq=',datq(iq,jq)

            elseif(cdatatype(1:lent).eq.'islope'    .or.
     &             cdatatype(1:lent).eq.'landuse'   .or.
     &             cdatatype(1:lent).eq.'soiltype'    )then

             call compute_categories(cdatatype,np*np,DATP(1,1,1)
     &               ,maxdatacat,domcat,pctcat)
             datq(iq,jq)=domcat 
             datqs(iq,jq,:)=pctcat(:)

            elseif(cdatatype(1:lent).eq.'greenfrac'   )then

c dominant greenness fraction for each month

             do lp=1,lcat
              call compute_categories(cdatatype,np*np,DATP(1,1,lp)
     &                               ,1,domcat,pctcat)
              datqs(iq,jq,lp)=domcat
             enddo

            endif

16       continue ! IQ
15    continue ! JQ

      print *,'after 15'
 
      XQ1=(1.-0.5*FLOAT(NIQ+1))*DELTAXQ+XCENTR
      YQ1=(1.-0.5*FLOAT(NJQ+1))*DELTAYQ+YCENTR

      if(cdatatype(1:lent).eq.'topography')then

        print*
        print*,'Before GDTOST2'
        print*,'--------------'
        print*,'datq(1,1)/(niq,njq)= ',datq(1,1),datq(niq,njq)
        print*,'datqs(1,1,1)/(niq,njq)= ',datqs(1,1,1),datqs(niq,njq,1)
        print*,'datsln(1,1)/(niq,njq)= ',datsln(1,1),datsln(niq,njq)
        print*,'datslt(1,1)/(niq,njq)= ',datslt(1,1),datslt(niq,njq)
        print*,'Mean/Max topo at IQ,JQ (1,1)/(niq,njq): '
     +,datsm(1,1),datsmx(1,1),datsm(niq,njq),datsmx(niq,njq)

        DO 28 JR=1,N3
         DO 29 IR=1,N2
           XR=(XT(IR)-XQ1)/DELTAXQ+1.
           YR=(YT(JR)-YQ1)/DELTAYQ+1.

           CALL GDTOST2(DATQ,NIQ,NJQ,XR,YR,RVAL)
           DATR(IR,JR)=max(0.,RVAL)
           if( DATR(IR,JR).gt.30000. )then
               print*,'Warning: value out of bounds'
           endif    

           CALL GDTOST2(DATQS,NIQ,NJQ,XR,YR,RVAL)
           DATS(IR,JR,1)=max(0.,RVAL)
           CALL GDTOST2(DATSLN,NIQ,NJQ,XR,YR,RVAL)
           DATLN(IR,JR)=RVAL
           CALL GDTOST2(DATSLT,NIQ,NJQ,XR,YR,RVAL)
           DATLT(IR,JR)=RVAL

 29      CONTINUE
 28     CONTINUE

        print*,'After GDTOST2'
        print*,'-------------'
        print*,'datr(1,1)/(n2,n3)= ',datr(1,1),datr(N2,N3)
        print*,'dats(1,1,1)/(n2,n3)= ',dats(1,1,1),dats(n2,n3,1)
        print*,'datln(1,1)/(n2,n3)= ',datln(1,1),datln(n2,n3)
        print*,'datlt(1,1)/(n2,n3)= ',datlt(1,1),datlt(n2,n3)
 
      elseif(cdatatype(1:lent).eq.'landuse'.or.
     +       cdatatype(1:lent).eq.'islope' .or.
     +       cdatatype(1:lent).eq.'soiltype')then

        DO 38 JR=1,N3
         DO 39 IR=1,N2
            IXR=NINT((XT(IR)-XQ1)/DELTAXQ)+1.
            IYR=NINT((YT(JR)-YQ1)/DELTAYQ)+1.
            if(ixr.lt.1)ixr=1
            if(iyr.lt.1)iyr=1
            if(ixr.gt.n2)ixr=niq
            if(iyr.gt.n3)iyr=njq

            datr(ir,jr)=  datq(ixr,iyr)     !dominant category
            dats(ir,jr,:)=datqs(ixr,iyr,:)  !percent dist for ea category

 39      CONTINUE
 38     CONTINUE

      endif

      deallocate(dato)
      deallocate(DATP,
     &           DATQ,
     &           DATQS,
     &           DATSM, 
     &           DATSMX,
     &           DATSLN,
     &           DATSLT)

      RETURN
      END

      subroutine vfirec(iunit,a,n,type)
      character*1 vc
      character*(*) type
      common/vform/vc(0:63)
      character line*80, cs*1
      dimension a(*)

      if(vc(0).ne.'0') call vfinit

      ich0=ichar('0')
      ich9=ichar('9')
      ichcz=ichar('Z')
      ichlz=ichar('z')
      ichca=ichar('A')
      ichla=ichar('a')
      
      read(iunit,10)nn,nbits,bias,fact
 10   format(2i8,2e20.10)
      if(nn.ne.n) then
         print*,' Word count mismatch on vfirec record '
         print*,' Words on record - ',nn
         print*,' Words expected  - ',n
         stop 'vfirec'
      endif

      nvalline=(78*6)/nbits
      nchs=nbits/6
      do 20 i=1,n,nvalline
         read(iunit,'(a78)') line
         ic=0
         do 30 ii=i,i+nvalline-1
            isval=0
            if(ii.gt.n) go to 20
            do 40 iii=1,nchs
               ic=ic+1
               cs=line(ic:ic)
               ics=ichar(cs)
               if(ics.le.ich9)then
                  nc=ics-ich0
               elseif(ics.le.ichcz) then
                  nc=ics-ichca+10
               else
                  nc=ics-ichla+36
               endif
               isval=ior(ishft(nc,6*(nchs-iii)),isval)
 40         continue
            a(ii)=isval
 30      continue
 20   continue

      facti=1./fact
      if(type.eq.'LIN') then
         do 48 i=1,n
            a(i)=a(i)*facti-bias
 48      continue
      elseif(type.eq.'LOG') then
         scfct=2.**(nbits-1)
         do 55 i=1,n
            a(i)=sign(1.,a(i)-scfct)
     +           *(10.**(abs(20.*(a(i)/scfct-1.))-10.))
 55      continue
      endif

      return
      end
c
cc ------------------------------------------------------------------
c
      subroutine vfinit                                                  
      character*1vc,vcscr(0:63)                                         
      common/vform/vc(0:63)                                             
      data vcscr/'0','1','2','3','4','5','6','7','8','9'                 
     +,'A','B','C','D','E','F','G','H','I','J'                          
     +,'K','L','M','N','O','P','Q','R','S','T'                          
     +,'U','V','W','X','Y','Z','a','b','c','d'                          
     +,'e','f','g','h','i','j','k','l','m','n'                          
     +,'o','p','q','r','s','t','u','v','w','x'                          
     +,'y','z','{','|'/                                                 
                                                                        
      do10n=0,63                                                        
      vc(n)=vcscr(n)                                                    
  10  continue                                                          
                                                                        
      return                                                            
      end
C +------------------------------------------------------------------+
      FUNCTION INTLSHFT(IWORD,NSHFT)
C
C       This function shifts IWORD to the left NSHFT bits in a
C         circular manner.
C
      INTLSHFT=ISHFT(IWORD,NSHFT)
      RETURN
      END
C +------------------------------------------------------------------+
      FUNCTION INTOR(IWORD1,IWORD2)
C
C       This function performs a bit-by-bit OR between IWORD1 and
C         IWORD2.
C
      INTOR=IOR(IWORD1,IWORD2)
      RETURN
      END


      SUBROUTINE BINOM2(X1,X2,X3,X4,Y1,Y2,Y3,Y4,XXX,YYY)
      implicit none
      real x1,x2,x3,x4,y1,y2,y3,y4,xxx,yyy,
     +   wt1,wt2,yz22,yz23,yz24,yz11,yz12,yz13,yoo
      integer istend
c      COMMON/BIN/ITYPP,I0X,I1X,I2X,YOO
       YYY=1E30
       IF(X2.GT.1.E19.OR.X3.GT.1.E19.OR.
     +   Y2.GT.1.E19.OR.Y3.GT.1.E19)RETURN
      WT1=(XXX-X3)/(X2-X3)
      WT2=1.0-WT1
      ISTEND=0
      IF(Y4.LT.1.E19.AND.X4.LT.1.E19) GO TO 410
      YZ22=WT1
      YZ23=WT2
      YZ24=0.0
      ISTEND= 1
410   IF(Y1.LT.1.E19.AND.X1.LT.1.E19) GO TO 430
      YZ11=0.0
      YZ12=WT1
      YZ13=WT2
      IF(ISTEND.EQ.1)GO TO 480
      GO TO 450
430   YZ11=(XXX-X2)*(XXX-X3)/((X1-X2)*(X1-X3))
      YZ12=(XXX-X1)*(XXX-X3)/((X2-X1)*(X2-X3))
      YZ13=(XXX-X1)*(XXX-X2)/((X3-X1)*(X3-X2))
      IF(ISTEND.EQ.  1    ) GO TO 470
450   YZ22=(XXX-X3)*(XXX-X4)/((X2-X3)*(X2-X4))
      YZ23=(XXX-X2)*(XXX-X4)/((X3-X2)*(X3-X4))
      YZ24=(XXX-X2)*(XXX-X3)/((X4-X2)*(X4-X3))
470   YYY=WT1*(YZ11*Y1+YZ12*Y2+YZ13*Y3)+WT2*(YZ22*Y2+YZ23*Y3+YZ24*Y4)
       GO TO 490
480      YYY=WT1*Y2+WT2*Y3
490   YOO=YYY
      RETURN
      END
c
c determine dominant category 1-05-01 JS
c

      subroutine compute_categories(ctype,nnp,data,nlcat,domcat
     +,pctcat)

      implicit none

      character*(*) ctype

      integer nnp
      integer igc
      integer i,j,k
      integer nlcat
      integer maxcat
      integer lcat(nlcat)
      real    pctcat(nlcat)
      real    data(nnp)
      real    domcat
      real    sum_g

c categorical data types
      if(ctype.eq.'landuse'   .or.
     +   ctype.eq.'soiltype'  .or.
     +   ctype.eq.'islope'  )then

         do k=1,nlcat
            lcat(k)=0
            pctcat(k)=0.
         enddo

         do i=1,nnp
         do k=1,nlcat
            if(nint(data(i)).eq.k)then
               lcat(k)=lcat(k)+1
            endif
         enddo
         enddo

         maxcat=-1
         do k=1,nlcat
            pctcat(k)=lcat(k)/float(nnp)
            if(lcat(k).gt.maxcat)then
               maxcat=lcat(k)
               domcat=float(k)
            endif
         enddo
         if(ctype.eq.'landuse')then
          if(pctcat(16).ge.0.5) then                            !!JBresch
	    domcat = 16                                         !!JBresch
          else if(domcat.eq.16.and.pctcat(16).lt.0.5)then       !!JBresch
            maxcat=-1
            do k=1,nlcat
               if(k.ne.16)then
                  if(lcat(k).gt.maxcat)then
                     maxcat=lcat(k)
                     domcat=float(k)
                  endif
               endif
            enddo
          endif
         endif

c quantitative data types
      elseif(ctype.eq.'greenfrac'.or.ctype.eq.'soiltemp')then

         sum_g=0.
         igc=0
         do i=1,nnp
          if(data(i).gt.0.0)then
             sum_g=sum_g+data(i)
             igc = igc+1
          endif
         enddo
         if(igc.gt.0)then
            domcat=sum_g/float(igc)
         else
            domcat=0.0
         endif

      endif

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine zero2d(adum2d,nnxp,nnyp)

	integer I,J,NNXP,NNYP
	real adum2d(nnxp,nnyp)
	
	do J=1,nnyp
	do I=1,nnxp
	adum2d(I,J)=0.
	enddo
	enddo

	end subroutine zero2d

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	subroutine dumtodata(adum2d,nnxp,nnyp,IVAL,data)

	INTEGER NNXP,NNYP,IVAL
	real ADUM2D(:,:),data(:,:,:)

	do J=1,nnyp
	do I=1,nnxp
	data(I,J,IVAL)=adum2d(I,J)
	enddo
	enddo

	end subroutine dumtodata

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine smdhld(ime,jme,h,s,lines,nsmud)
!      parameter(ime=im+2,jme=jm+4)
      dimension ihw(jme),ihe(jme)
      dimension h(ime,jme),s(ime,jme)
     &         ,hbms(ime,jme),hne(ime,jme),hse(ime,jme)
!-----------------------------------------------------------------------
          do j=1,jme
      ihw(j)=-mod(j,2)
      ihe(j)=ihw(j)+1
          enddo
!-----------------------------------------------------------------------

              do j=1,jme
          do i=1,ime
      hbms(i,j)=1.-s(i,j)
          enddo
              enddo
!
      jmelin=jme-lines+1
      ibas=lines/2
      m2l=mod(lines,2)
!
              do j=lines,jmelin
          ihl=ibas+mod(j,2)+m2l*mod(j+1,2)
          ihh=ime-ibas-m2l*mod(j+1,2)
!
          do i=ihl,ihh
      hbms(i,j)=0.
          enddo
              enddo
!-----------------------------------------------------------------------
                  do ks=1,nsmud
!-----------------------------------------------------------------------
              do j=1,jme-1
          do i=1,ime-1
      hne(i,j)=h(i+ihe(j),j+1)-h(i,j)
          enddo
              enddo
              do j=2,jme
          do i=1,ime-1
      hse(i,j)=h(i+ihe(j),j-1)-h(i,j)
          enddo
              enddo
!
              do j=2,jme-1
          do i=1+mod(j,2),ime-1
      h(i,j)=(hne(i,j)-hne(i+ihw(j),j-1)
     &       +hse(i,j)-hse(i+ihw(j),j+1))*hbms(i,j)*0.125+h(i,j)
          enddo
              enddo
!-----------------------------------------------------------------------
                      enddo
!-----------------------------------------------------------------------
      return
      end
