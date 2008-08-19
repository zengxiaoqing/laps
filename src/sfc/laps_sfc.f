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
c
c
	program laps_sfc

!       We are using just a subset of the global parameters in this driver
        use mem_namelist, ONLY: read_namelist_laps
        use mem_namelist, ONLY: NX_L,NY_L,nk_laps,maxstns,
     &                          laps_cycle_time,grid_spacing_m

        use mem_grid, ONLY: lat,lon,topo,ldf

        use mem_sfcanl

        character*150 static_dir,filename
        character*9 a9time
c
!       Read global parameters into module memory structure
        call get_directory('static',static_dir,len_dir)
        filename = static_dir(1:len_dir)//'/nest7grid.parms'
        call read_namelist_laps('lapsparms',filename)

!       Read surface parameters into module memory structure
        filename = static_dir(1:len_dir)//'/surface_analysis.nl'
        call read_namelist_laps('sfc_anal',filename)

!       Allocate static arrays (lat, lon, topo, ldf)
        allocate( lat(NX_L,NY_L), STAT=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' ERROR: Could not allocate lat'
            stop
        endif

        allocate( lon(NX_L,NY_L), STAT=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' ERROR: Could not allocate lon'
            stop
        endif

        allocate( topo(NX_L,NY_L), STAT=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' ERROR: Could not allocate topo'
            stop
        endif

        allocate( ldf(NX_L,NY_L), STAT=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' ERROR: Could not allocate ldf'
            stop
        endif

!       Read static arrays
        call read_static_grid(NX_L,NY_L,'LAT',lat,istatus)
        if(istatus .ne. 1)then
            stop
        endif

        call read_static_grid(NX_L,NY_L,'LON',lon,istatus)
        if(istatus .ne. 1)then
            stop
        endif

        call read_static_grid(NX_L,NY_L,'AVG',topo,istatus)
        if(istatus .ne. 1)then
            stop
        endif

        call read_static_grid(NX_L,NY_L,'LDF',ldf,istatus)
        if(istatus .ne. 1)then
            stop
        endif

!       Get System Analysis Time
        call get_systime(i4time,a9time,istatus)

!       Note these are being allocated later on just before they are filled
!       call alloc_sfc_arrays(NX_L,NY_L)
!       call point_sfcanl_arrays(NX_L,NY_L)

	call laps_sfc_sub(i4time)

        call deallocate_sfcanl_arrays()
c
	end
c
