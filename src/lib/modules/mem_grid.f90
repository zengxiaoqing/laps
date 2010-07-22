Module mem_grid

real   , allocatable, dimension(:,:) :: lat, lon, topo, ldf
integer, allocatable, dimension(:,:) :: IstartIend
integer, allocatable, dimension(:)   :: recvcounts,displs
integer                              :: nPEs,rank

end Module
