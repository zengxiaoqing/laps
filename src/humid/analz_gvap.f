      subroutine analz_gvap (lat,lon,wt,nn,glat,glon,
     1                       data_out,data_weights,ii,jj,istatus)

      integer ii,jj,nn,istatus
      real lat(nn),lon(nn),wt(nn)
      real data_out(ii,jj)
      real data_weights(ii,jj)
      real glat(ii,jj)
      real glon(ii,jj)

c     volitile array points
      real points (nn,3)
      integer mask (ii,jj)


      integer i,j,n,ncount
      real sum

c
      istatus =  0

c     initialaize data_out on first call  (-1 used as missing)
c     since zero water is valid number

      do i = 1,ii
         do j = 1,jj
            data_out(i,j) = -1.
            mask (i,j) = 0
            data_weights(i,j) = 0.0
         enddo
      enddo



c     foreach n element of wt, determine its location in ii,jj space

      ncount = 0

      do n = 1,nn

         call  latlon_to_rlapsgrid(lat(n),lon(n),glat,glon,ii,jj,
     1        ri,rj, istatus)

         if (istatus.eq.1) then
            ncount = ncount + 1
            i = nint (ri)
            j = nint (rj)
            data_out (i,j) = wt(n)
            mask(i,j) = 1
            data_weights (i,j) = 1.
            points(ncount,1) = wt(n)
            points(ncount,2) = i
            points(ncount,3) = j
         else
            continue
         endif

      enddo

c     now that data_out is as full as it is going to get with the water
c     process the remainder.

c     compute the fraction of data_out that is empty

      sum = 0.

      do i = 1,ii
         do j = 1,jj
            if (data_out (i,j) .eq. -1.) then
               sum = sum + 1
               data_out(i,j) = points(1,1) ! helps converge
            endif
         enddo
      enddo

      write (6,*) 'gvap field fraction empty = ', sum/ii/jj

c      if (sum/ii/jj .gt. 0.75) return ! return if not enough data

      
c     now have fairly full data array.  now analyze

      call prep_grid (ii,jj,data_out,nn,points,ncount)
      call slv_laplc (data_out,mask,ii,jj)

c     prep the weighting array for the above analyzed sheet

      call slv_laplc (data_weights, mask, ii,jj)


      istatus = 1

      return

      end


