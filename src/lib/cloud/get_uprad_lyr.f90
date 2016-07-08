
        subroutine get_uprad_lyr(nc,ni,nj,gnd_radc,ht,uprad_3d)

        use mem_namelist, ONLY: r_missing_data,earth_radius,grid_spacing_m 
        parameter (rpd = 3.14159/180.)

        real gnd_radc(nc,ni,nj) ! spectral exitance (from sfc lights)
        real sumrad(nc,ni,nj)
        real uprad_3d(ni,nj,nc) ! spectral upward irradiance for layer
        real, allocatable :: drad(:,:)

        radius = 60000.
        iradius = nint(radius / grid_spacing_m)

        write(6,*)' ht/radius/iradius = ',ht,radius,iradius

        allocate(drad(-iradius:+iradius,-iradius:+iradius))

!       Determine radiation weighting function array
        do ii = -iradius,+iradius
        do jj = -iradius,+iradius
            disti = float(ii) * grid_spacing_m
            distj = float(jj) * grid_spacing_m
            distr = sqrt(disti**2+distj**2+ht**2)
            sin_theta_r = ht/distr
            stearadians = sin_theta_r * (grid_spacing_m / distr)**2
            drad(ii,jj) = stearadians
!           if(ii .eq. 0)then
!               write(6,*)'jj|disti|distj|ht|distr|stearadians',jj,disti,distj,ht,distr,stearadians
!           endif
        enddo ! ii
        enddo ! jj

        write(6,*)' range of gnd_radc (red: wm2nm) = ',minval(gnd_radc(1,:,:)),maxval(gnd_radc(1,:,:))
        write(6,*)' range of drad = ',minval(drad),maxval(drad)
        write(6,*)' sum of drad = ',sum(drad)
        sumrad = 0. ! initialize
        uprad_3d(:,:,:) = 0.

        if(iradius .ge. 100)then
          iskip = 10
        elseif(iradius .ge. 50)then
          iskip = 3
        elseif(iradius .ge. 20)then
          iskip = 2
        else
          iskip = 1
        endif

        write(6,*)' iskip = ',iskip

!       Loop through each gridpoint on the layer array
        do i = 1,ni,iskip
        do j = 1,nj,iskip

!         Index limits of layer array
          imin = max(i-iradius,1)
          imax = min(i+iradius,ni)
          jmin = max(j-iradius,1)
          jmax = min(j+iradius,nj)

!         Index limits on the weighting array
          iimin = imin-i
          iimax = imax-i
          jjmin = jmin-j
          jjmax = jmax-j
          do ic = 1,nc
            uprad_3d(i,j,ic) = sum(drad(iimin:iimax,jjmin:jjmax) &
                                 * gnd_radc(ic,imin:imax,jmin:jmax))
          enddo ! ic

        enddo ! j
        enddo ! i

!       Interpolate to fill in the horizontal
!       Note that the integrated radiance over the hemisphere gives
!       irradiance. We can multiply that by the extinction
!       coefficient to obtain the emission density "S", equivalent to
!       radiant flux per unit volume.
        do ic = 1,nc
           write(6,*)' fill horizontal layer to yield uprad (wm2nm) for color',ic
           write(6,*)' range of uprad ',minval(uprad_3d(:,:,ic)),maxval(uprad_3d(:,:,ic))
           call bilinear_fill(uprad_3d(:,:,ic),ni,nj,iskip,r_missing_data)
           write(6,*)' range of uprad ',minval(uprad_3d(:,:,ic)),maxval(uprad_3d(:,:,ic))
        enddo ! ic

        deallocate(drad)

        return
        end
