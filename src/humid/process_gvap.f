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
      subroutine  process_gvap (ii,jj,data_out,data_weights,
     1     tpw,glat,glon,time_diff,
     1     path_to_gvap8,path_to_gvap10,filetime,istatus)

      implicit none

c     input variables

      character*9 filename,filetime
      integer ii,jj,istatus
      integer time_diff         !time allowed for latency (sec)
      real data_out(ii,jj),tpw(ii,jj)
      real glat(ii,jj), glon(ii,jj)
      integer i4time
      character*256 path_to_gvap8,path_to_gvap10
      real data_weights(ii,jj)

      integer nstations,nn
      parameter (nstations = 11000)
      real lat(nstations)
      real lon(nstations)
      real wt(nstations)
      real w1(nstations)
      real w2(nstations)
      real w3(nstations)



      integer i,j

      filename = filetime

      call read_gvap (filename, nstations, path_to_gvap8,path_to_gvap10,
     1     time_diff, lat,lon, wt,w1,w2,w3, nn,
     1     istatus)

      if (
     1     istatus .ne. 1
     1     .or.
     1     nn .eq. 0
     1     ) then               ! failure

         write(6,*) 'failure to acquire gvap data'
         istatus = 0
         return                 !istatus = fail

      else


         write(6,*) nn, ' number of stations read in file'
         if (nn.gt.nstations) then ! exceeded dimension
            istatus  = 0
            write (6,*) 'nstations exceeded (parameter dimension)'
            write (6,*) 'readjust and recompile code'
            write (6,*) 'warning only, not fatal'
            write (6,*) 'gvap data not used'
            return
         endif
c     correct longitute to negative for west
         do i = 1,nn
            lon(i) = lon(i) * (-1.0)
         enddo
         write(6,*) w3(nn)

      endif

      call analz_gvap (lat,lon,wt,nn,glat,glon,data_out,
     1     data_weights,ii,jj,istatus)

      if(istatus.ne.1) then ! failure to get data
         return
      endif


      if(istatus.eq.1) then ! data_out can be used to normalize field
c     note that the 0.1 factor is to convert mm (gvap) to cm (tpw).
         do i   = 1,ii
            do j  = 1,jj
               data_out(i,j) = data_out(i,j)*0.1/tpw(i,j)
            enddo
         enddo

      endif

c     data_out is now a fractional adjustment (weighted)
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
