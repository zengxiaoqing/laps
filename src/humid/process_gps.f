cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps
cdis
cdis    this software and its documentation are in the public domain and
cdis    are furnished "as is."  the united states government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  they assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  if significant modifications or enhancements
cdis    are made to this software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis
cdis
cdis
cdis
cdis
cdis
cdis
      subroutine  process_gps (ii,jj,data_out,data_weights,
     1     tpw,glat,glon,time_diff,
     1     path,filetime,istatus)

      implicit none

c     input variables


      character*9 filetime
      integer ii,jj,istatus
      integer time_diff         !time allowed for latency (sec)
      real data_out(ii,jj),tpw(ii,jj)
      real glat(ii,jj), glon(ii,jj)
      integer i4time
      character*256 path
      real data_weights(ii,jj)
      
      integer gps_n, gps_num
      parameter (gps_n = 1000)
      real gps_tpw(gps_n)
      real gps_error(gps_n)
      real gps_lat(gps_n)
      real gps_lon(gps_n)

      integer i,j


      call read_gps (path, filetime, time_diff,
     1     gps_tpw, gps_error, gps_lat,
     1     gps_lon, gps_num, gps_n,
     1     istatus)

      if (
     1     istatus .ne. 1
     1     .or.
     1    gps_num .eq. 0
     1     ) then               ! failure

         write(6,*) 'failure to acquire gvap data'
         istatus = 0
         return                 !istatus = fail

      else


         write(6,*) gps_num, ' number of stations read in file'


      endif

      call analz_gvap (gps_lat,gps_lon,gps_tpw,gps_num,glat,
     1     glon,data_out,
     1     data_weights,ii,jj,istatus)

      if(istatus.ne.1) then ! failure to get data
         return
      endif

      if(istatus.eq.1) then ! data_out can be used to normalize field
c     note that the 0.1 factor is to convert mm (gvap) to cm (tpw).
         do i   = 1,ii
            do j  = 1,jj
               data_out(i,j) = data_out(i,j)/tpw(i,j)
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

      istatus = 1

      return

      end
