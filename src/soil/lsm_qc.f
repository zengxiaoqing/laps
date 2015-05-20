cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
       subroutine lsm_qc_check(nx_l,ny_l,nz,
     &                         data,r_low,r_high,istatus)
c
c routine determines if grid point has "bad" value. If so then the surrounding
c points are used to derive a good value.  The routine does not consider soil
c type boundaries at this time.
c
       implicit none

       integer*4 nx_l,ny_l,nz
       real*4 data(nx_l,ny_l,nz)
       real*4 temp(9)
       real*4 r_low,r_high
       real*4 r_missing_data
       real*4 ave,adev,sdev,var,skew,curt
       integer*4 istatus
       integer*4 i,j,k,n
       integer*4 icntm,icnts
       integer*4 ii, jj, istat, istat_nan
c
       icnts = 0
       icntm = 0
c
       call get_r_missing_data(r_missing_data,istatus)
       if(istatus .ne. 1)then
          write(6,*)'error getting missing data'
          goto 1000
       endif

       do k=1,nz   !so as not to confuse with laps vertical grid
       do j=2,ny_l-1
       do i=2,nx_l-1

          call check_nan(data(i,j,k),istat_nan)

          if(data(i,j,k).lt.r_low .or. data(i,j,k).gt.r_high .or.
     1       istat_nan .eq. 0)then

             n=0
             do jj = j-1,j+1
             do ii = i-1,i+1

                if(data(ii,jj,k).ge.r_low.and.
     &data(ii,jj,k).le.r_high)then
                   n=n+1
                   temp(n) = data(ii,jj,k)
                endif

             enddo
             enddo

             if(n.ge.2)then
                call moment(temp,n,ave,adev,sdev,var,skew,curt,
     &                      istat)
                data(i,j,k)=ave
                icnts=icnts-1
             else
                data(i,j,k)=r_missing_data
                icntm=icntm-1
             endif

          endif

       enddo !i
       enddo !j
c
c take of boundaries
c
       do i=1,nx_l
          call check_nan(data(i,ny_l,k),istat_nan)
          if(data(i,ny_l,k).lt.r_low.or.
     &       data(i,ny_l,k).gt.r_high.or.istat_nan.eq.0)then
             data(i,ny_l,k)=r_missing_data
             icntm=icntm-1
          endif
          call check_nan(data(i,1,k),istat_nan)
          if(data(i,1,k).lt.r_low.or.
     &       data(i,1,k).gt.r_high.or.istat_nan.eq.0)then
             data(i,1,k)=r_missing_data
             icntm=icntm-1
          endif
       enddo
c
       do j=1,ny_l
          call check_nan(data(nx_l,j,k),istat_nan)
          if(data(nx_l,j,k).lt.r_low.or.
     &       data(nx_l,j,k).gt.r_high.or.istat_nan.eq.0)then
             data(ny_l,j,k)=r_missing_data
             icntm=icntm-1
          endif
          call check_nan(data(1,j,k),istat_nan)
          if(data(1,j,k).lt.r_low.or.
     &       data(1,j,k).gt.r_high.or.istat_nan.eq.0)then
             data(1,j,k)=r_missing_data
             icntm=icntm-1
          endif
       enddo
             
       enddo !k
c
      istatus = 1
      if(icntm .lt. 0) then
        write(6,*)'Set to r_missing_data: ',abs(icntm)
        istatus = icntm
      endif
      if(icnts .lt. 0)then
        write(6,*)'Set to average :',abs(icnts)
        istatus = istatus+icnts
      endif
c
1000  return
      end
