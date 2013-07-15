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

        include 'trigd.inc'

        use mem_namelist, ONLY: read_namelist_laps
        use mem_namelist, ONLY: NX_L,NY_L,nk_laps
 
        use mem_grid, ONLY: lat,lon,topo

        use mem_temp

        character*9 a9time        
        character*150 static_dir,filename

        real testlat_a(4),testlon_a(4)
c
!       Read global parameters into module memory structure
        call get_directory('static',static_dir,len_dir)
        filename = static_dir(1:len_dir)//'/nest7grid.parms'
        call read_namelist_laps('lapsparms',filename)

!       Read temp parameters into module memory structure
        filename = static_dir(1:len_dir)//'/temp.nl'
        call read_namelist_laps('temp_anal',filename)

!       Allocate static arrays (lat, lon, topo)
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

!       Read static arrays (lat, lon, topo)
        call read_static_grid(NX_L,NY_L,'LAT',lat,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error getting LAPS LAT'
            stop
        endif

        call read_static_grid(NX_L,NY_L,'LON',lon,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error getting LAPS LON'
            stop
        endif

        call read_static_grid(NX_L,NY_L,'AVG',topo,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error getting LAPS topo'
            stop
        endif

        testlat_a(1) = 43.5519
        testlon_a(1) = -136.1232

        testlat_a(2) = 47.308
        testlon_a(2) = -114.4814

        testlat_a(3) = 28.89283
        testlon_a(3) = -129.7618

        testlat_a(4) = 31.74828
        testlon_a(4) = -112.5754

        testlat = 41.0
        testlon = -98.0

        do il = 1,4

          testlat = testlat_a(il)
          testlon = testlon_a(il)

          rclo = 999999.

          do i = 1,NX_L
          do j = 1,NY_L
            deltalon = (lon(i,j) - testlon) * cosd(testlat)
            dist = sqrt( (lat(i,j)-testlat)**2 + deltalon**2) 

!           write(6,*)lat(i,j),lon(i,j),dist
            if(dist .lt. rclo)then
                rclo = dist
!               write(6,*)' new close point ',i,j,lat(i,j),lon(i,j),rclo
                iclo = i
                jclo = j
            endif

          enddo ! j
          enddo ! i

          write(6,*)' corner/i/j/rclo(km) ',il,iclo,jclo,rclo*110.

        enddo ! il

999     continue

        end

