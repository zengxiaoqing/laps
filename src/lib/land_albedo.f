
      subroutine land_albedo(lu,ni,nj,albedo)

!     This routine takes land use categories and broadband albedo relationship
!     to derive a 3 color albedo

      real lu(ni,nj)          ! Land Use (USGS 24 category)
      real albedo(3,ni,nj)    ! Albedo (Red, Green, Blue)

!     See http://www.mmm.ucar.edu/mm5/mm5v2/landuse-usgs-tbl.html

      do i = 1,ni
      do j = 1,nj
          if(lu(i,j) .eq. 1.)then      ! Urban and Built-up Land   (.18 Gray)
              albedo(1,i,j) = 0.18
              albedo(2,i,j) = 0.18
              albedo(3,i,j) = 0.18
          elseif(lu(i,j) .eq. 2.)then  ! Dryland Cropland/Pasture  (.17 Brown)
              albedo(1,i,j) = 0.17
              albedo(2,i,j) = 0.17
              albedo(3,i,j) = 0.01
          elseif(lu(i,j) .eq. 3.)then  ! Irrigtd Cropland/Pasture  (.18 Green)
              albedo(1,i,j) = 0.13
              albedo(2,i,j) = 0.18
              albedo(3,i,j) = 0.09
          elseif(lu(i,j) .eq. 4.)then  ! Mixed                     (.18 Green)
              albedo(1,i,j) = 0.13
              albedo(2,i,j) = 0.18
              albedo(3,i,j) = 0.09
          elseif(lu(i,j) .eq. 5.)then  ! Cropland Grassland Mosaic (.18 Green)
              albedo(1,i,j) = 0.13
              albedo(2,i,j) = 0.18
              albedo(3,i,j) = 0.09
          elseif(lu(i,j) .eq. 6.)then  ! Cropland Woodland Mosaic  (.16 Green)
              albedo(1,i,j) = 0.12
              albedo(2,i,j) = 0.16
              albedo(3,i,j) = 0.08
          elseif(lu(i,j) .eq. 7.)then  ! Grassland                 (.19 Brown)
              albedo(1,i,j) = 0.19
              albedo(2,i,j) = 0.19
              albedo(3,i,j) = 0.01
          elseif(lu(i,j) .eq. 8.)then  ! Shrubland                 (.22 Brown)
              albedo(1,i,j) = 0.22
              albedo(2,i,j) = 0.22
              albedo(3,i,j) = 0.01
          elseif(lu(i,j) .eq. 9.)then  ! Mixed Shrubland/Grassland (.20 Brown)
              albedo(1,i,j) = 0.20
              albedo(2,i,j) = 0.20
              albedo(3,i,j) = 0.01
          elseif(lu(i,j) .eq. 10.)then ! Savanna                   (.20 Brown)
              albedo(1,i,j) = 0.20
              albedo(2,i,j) = 0.20
              albedo(3,i,j) = 0.01
          elseif(lu(i,j) .eq. 11.)then ! Deciduous Broadleaf       (.16 Green)
              albedo(1,i,j) = 0.12
              albedo(2,i,j) = 0.16
              albedo(3,i,j) = 0.08
          elseif(lu(i,j) .eq. 12.)then ! Deciduous Needleleaf      (.14 Green)
              albedo(1,i,j) = 0.10
              albedo(2,i,j) = 0.14
              albedo(3,i,j) = 0.07
          elseif(lu(i,j) .eq. 13.)then ! Evergreen Broadleaf       (.12 Green)
              albedo(1,i,j) = 0.09
              albedo(2,i,j) = 0.12
              albedo(3,i,j) = 0.06
          elseif(lu(i,j) .eq. 14.)then ! Evergreen Needleleaf      (.12 Green)
              albedo(1,i,j) = 0.09
              albedo(2,i,j) = 0.12
              albedo(3,i,j) = 0.06
          elseif(lu(i,j) .eq. 15.)then ! Mixed Forest              (.13 Green)
              albedo(1,i,j) = 0.10
              albedo(2,i,j) = 0.13
              albedo(3,i,j) = 0.07
          elseif(lu(i,j) .eq. 16.)then ! water (.08)
              albedo(1,i,j) = 0.04
              albedo(2,i,j) = 0.04
              albedo(3,i,j) = 0.12
          elseif(lu(i,j) .eq. 19.)then ! Barren or Sparsely Vegetated (.25 Brown)
              albedo(1,i,j) = 0.25
              albedo(2,i,j) = 0.25
              albedo(3,i,j) = 0.10
          else ! default
              albedo(1,i,j) = 0.19
              albedo(2,i,j) = 0.19
              albedo(3,i,j) = 0.01
          endif                        
      enddo ! j
      enddo ! i

      return
      end
