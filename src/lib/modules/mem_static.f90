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



Module mem_static

  type static_type
  integer                     :: nx             ! X grid points
  integer                     :: ny             ! Y grid points
  integer                     :: nz             ! Z grid points
  real                        :: swlat          ! SW Corner Lat
  real                        :: swlon          ! SE Corner Lon
  real                        :: nelat          ! NE Corner Lat
  real                        :: nelon          ! NE Corner Lon
  real                        :: dx             ! X-direction grid spacing
  real                        :: dy             ! Y-direction grid spacing
  real                        :: cenlon         ! center Longitude
  real                        :: cenlat         ! center Latitude
  real                        :: stdlon         ! Orientation Longitude (polelon)
  real                        :: stdlat1        ! Standard Lat 1
  real                        :: stdlat2        ! Standard Lat 2 (polelat)
  real, pointer               :: topo(:,:)      ! LAPS Topographic height
  real, pointer               :: glat(:,:)      ! LAPS Latitude Array
  real, pointer               :: glon(:,:)      ! LAPS Longitude Array
  real, pointer               :: ldf(:,:)       ! Land fraction
  real, pointer               :: pbl(:,:)       ! PBL height
  real, pointer               :: akk(:,:)       ! Drag coefficient
  character (len=132)         :: grid_type      ! Map projection type
  end type

  type(static_type) :: stgrid

  real, pointer, dimension(:,:) :: lat, lon, topo, ldf

Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine alloc_static_arrays (nx, ny)

implicit none
integer :: nx, ny

    allocate ( stgrid%glat (nx,ny) )
    allocate ( stgrid%glon (nx,ny) )
    allocate ( stgrid%topo (nx,ny) )
    allocate ( stgrid%ldf  (nx,ny) )
    allocate ( stgrid%pbl (nx,ny) )
    allocate ( stgrid%akk  (nx,ny) )

return
end subroutine alloc_static_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine point_static_arrays ()

implicit none

lat =>  stgrid%glat
lon =>  stgrid%glon
topo  =>  stgrid%topo
ldf  =>  stgrid%ldf

return
end subroutine point_static_arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dealloc_static_arrays ()

implicit none

    deallocate (stgrid%glat)
    deallocate (stgrid%glon)
    deallocate (stgrid%topo)
    deallocate (stgrid%ldf )
    deallocate (stgrid%pbl)
    deallocate (stgrid%akk )

return
end subroutine dealloc_static_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fill_static      (imax, jmax                           &
                            ,n_flds, dx, dy, lov, latin1         &
                            ,latin2, origin, var, comment         &
                            ,data, model, grid_spacing,map_proj   &
                            ,la1,lo1,la2,lo2                      &
                            ,center_lat, center_lon,lli,llj,uri   &
                            ,urj,parent_id,ratio_2_parent         &
                            ,status)

implicit none

integer     ::  imax, jmax, n_flds,status  &
               ,lli,llj,uri,urj,parent_id,ratio_2_parent
real        ::  dx, dy, lov, latin1, latin2
character(len=*) :: origin,var(n_flds),comment(n_flds),model,map_proj
real        :: data(imax,jmax,n_flds)
real        :: grid_spacing,la1,lo1,la2,lo2,center_lat,center_lon

stgrid%nx               = imax
stgrid%ny               = jmax
stgrid%swlat            = la1
stgrid%swlon            = lo1
stgrid%nelat            = la2
stgrid%nelon            = lo2
stgrid%dx               = grid_spacing
stgrid%dy               = grid_spacing
stgrid%cenlon           = center_lon
stgrid%cenlat           = center_lat
stgrid%stdlon           = lov
stgrid%stdlat1          = latin1
stgrid%stdlat2          = latin2
stgrid%topo(:,:)        = data(:,:,3)
stgrid%glat(:,:)        = data(:,:,1)
stgrid%glon(:,:)        = data(:,:,2)
stgrid%ldf(:,:)         = data(:,:,4)
stgrid%grid_type        = map_proj

!print*,'-----------:',imax,jmax,la1,lo1,la2,lo2,grid_spacing,grid_spacing  &
!,center_lon,center_lat,lov,latin1,latin2,data(1,1,3),data(1,1,1),data(1,1,2) &
!,map_proj


!         var(1)    = 'LAT'
!         var(2)    = 'LON'
!         var(3)    = 'AVG'
!         var(4)    = 'LDF'
!         var(5)    = 'USE'
!         var(6)    = 'ALB'  !now used for max snow alb 2-20-03 JS.
!         var(7)    = 'STD'
!         var(8)    = 'SLN'
!         var(9)    = 'SLT'
!         var(10)   = 'STL'
!         var(11)   = 'SBL'
!         var(12)   = 'LND'  ! Land-Water Mask based on USGS landuse
!         i=12
!         do j=1,12
!            write(cat,'(i2.2)')j
!            if(cat(1:1).eq.' ')cat(1:1)='0'
!            if(cat(2:2).eq.' ')cat(2:2)='0'
!            var(i+j)= 'G'//cat   ! vegetation greenness fraction
!         enddo
!
!         var(25)='TMP'
!         i=25
!         do j=1,12
!            write(cat,'(i2.2)')j
!            if(cat(1:1).eq.' ')cat(1:1)='0'
!            if(cat(2:2).eq.' ')cat(2:2)='0'
!            var(i+j)= 'A'//cat   ! monthly albedo
!         enddo
!
!         var(ngrids)   = 'ZIN'
!
!         comment(1) = 'Lat: From MODEL by J. Smart/ S. Albers 2-03\0'
!         comment(2) = 'Lon: From MODEL by J. Smart/ S. Albers 2-03\0'
!         comment(3) = 'Average terrain elevation (m) \0'
!         comment(4) = 'Land Fraction: derived from USGS land use \0'
!         comment(5) = 'Land Use dominant category (USGS 24 Category) \0'
!         comment(6) = 'Maximum Snow Albedo; defined over land only \0'
!         comment(7) = 'Standard Deviation of Elevation data (m)\0'
!         comment(8) = 'Mean longitudinal terrain slope (m/m)\0'
!         comment(9) = 'Mean latitudinal terrain slope (m/m)\0'
!         comment(10)= 'Top layer (0-30cm) dominant category soiltype\0'
!         comment(11)= 'Bot layer (30-90cm) dominant category soiltype\0'
!         comment(12)= 'Land-Water Mask (0=water; 1 otherwise) \0'

return
end subroutine fill_static

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module mem_static
