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
cdis cdis
cdis
cdis
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

      call check_nan2(data_out,ii,jj,istatus)
      if(istatus.ne.1) then
         write(6,*) 'data_out corrupted by tpw divide??'
         write(6,*) 'Nan detected var:data_out  routine:process_gvap.f'
         return
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
