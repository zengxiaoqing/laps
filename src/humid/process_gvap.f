      subroutine  process_gvap (ii,jj,data_out,
     1     tpw,glat,glon,
     1     filetime,istatus)

      implicit none

c     input variables

      character*9 filename,filetime
      integer ii,jj,istatus
      real data_out(ii,jj),tpw(ii,jj)
      real glat(ii,jj), glon(ii,jj)

      integer nstations,nn
      parameter (nstations = 9000)
      real lat(nstations)
      real lon(nstations)
      real wt(nstations)
      real w1(nstations)
      real w2(nstations)
      real w3(nstations)

c     volitile arrays
      real data_weights(ii,jj)


      integer i,j

      filename = '990531220' 


      call read_gvap (filename, nstations, lat,lon, wt,w1,w2,w3, nn,
     1 istatus)

      write(6,*) nn, ' number of stations read in file'
      write(6,*) w3(nn)

      call analz_gvap (lat,lon,wt,nn,glat,glon,data_out,
     1     data_weights,ii,jj,istatus)

      if(istatus.eq.1) then ! data_out can be used to normalize field
         do i   = 1,ii
            do j  = 1,jj
               data_out(i,j) = data_out(i,j)*0.1/tpw(i,j)
            enddo
         enddo

      endif

c     data_out is now a fractional adjustment
c     data_weights is how much of that fraction should be applied
c     convert data_out to incremental weighted adjustment

      do j = 1,jj
         do i = 1,ii
            data_out(i,j) = (data_out(i,j)-1.0) * data_weights(i,j)
         enddo
      enddo

      write(6,*) 'routine stopping for gvap testing'

      istatus = 1

      end
