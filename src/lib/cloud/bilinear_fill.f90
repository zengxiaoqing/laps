
      subroutine bilinear_fill(t,imax,jmax,nskip,r_missing_data)

      integer lowi_lut(imax)
      integer lowj_lut(jmax)

      dimension t(imax,jmax)

!     Bilinearly interpolate to fill in rest of domain
!     Fills in final analysis value and weights from obs alone
!     We may have to extrapolate at the N and E edges
      do i = 1,imax
          lowi_lut(i) = ((i-1)/nskip)*nskip + 1
          il = lowi_lut(i)
          ih = il + nskip
          if(ih .gt. imax)lowi_lut(i) = lowi_lut(i) - nskip
!         write(6,*)' i,il,ih',i,il,ih
      enddo ! i
      do j = 1,jmax
          lowj_lut(j) = ((j-1)/nskip)*nskip + 1
          jl = lowj_lut(j)
          jh = jl + nskip
          if(jh .gt. jmax)lowj_lut(j) = lowj_lut(j) - nskip
      enddo ! i

      do j=1,jmax
          jl = lowj_lut(j)
          jh = jl + nskip
          fracj = dble(j-jl)/dble(nskip)

          do i=1,imax

              il = lowi_lut(i)
              ih = il + nskip
              fraci = dble(i-il)/dble(nskip)

!             write(6,*)' i,il,ih',i,il,ih
!             write(6,*)i,j,il,ih,jl,jh

!             Calculate interpolated cloud cover
              Z1=t(il,jl)
              Z2=t(ih,jl)
              Z3=t(ih,jh)
              Z4=t(il,jh)

              t(i,j) = Z1+(Z2-Z1)*fraci+(Z4-Z1)*fracj &
                     - (Z2+Z4-Z3-Z1)*fraci*fracj

          enddo ! i
      enddo ! j

      return
      end

