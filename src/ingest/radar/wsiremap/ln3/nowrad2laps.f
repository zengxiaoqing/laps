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
       subroutine NOWRAD_to_LAPS(filename,
     &                    nlines,nelements,nlev,
     &                    imax,jmax,
     &                    lat,lon,
     &                    validTime,
     &                    remapped_prod,
     &                    istatus)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Routine reads netCDF nowrad (WSI) data using subroutine read_wsi_cdf.
c Nowrad data is then remapped to LAPS domain given lat/lon of domain.
c Routine automatically moves the boundary for any domain with the nowrad
c confines.
c
       integer   nlev

       real*4  lat(imax,jmax)
       real*4  lon(imax,jmax)
       real*4  remapped_prod(imax,jmax)
       integer sum(imax,jmax)
       integer count(imax,jmax)
       integer i_max_value(imax,jmax)
       integer validTime
       integer istatus
       integer jstatus
       integer lstatus
       integer data_levels(nlev)
       integer il(nlev)
       integer ul
       integer ndvalue
       integer num_levels
       integer msng_radar

       integer istart,jstart
       integer iend,jend
       integer ickint,itotwait,iageth

       character filename*200
       character c_data_type*2
       character level_prefix(nlev)*50
       character lprefix(nlev)*2

       real*4 dgtord
       real*4 rlat1
       real*4 rlon1
       real*4 rdlat
       real*4 rdlon

       real*4 r_missing_data

       integer image(nelements,nlines)
       integer imdata
c
c ***************************************************************************
c ************************** GET NOWRAD DATA ********************************
c
      istatus=1

      n=index(filename,' ')
      c_data_type=filename(n-2:n-1)

      rdtodg=180.0/3.1415926
      dgtord=0.017453292
c
c eventually this should be put in nest7grid.parms and accessed through wsi_ln3
c i/o routine.
c
      call get_ref_base(ref_base,iostatus)
      call get_ref_base_useable(ref_base_useable,iostatus)
      call get_r_missing_data(r_missing_data,iostatus)

      call get_ln3_parameters(msng_radar,ickint,itotwait,iageth,
     &             istart,jstart,iend,jend,lstatus)
      if(lstatus.ne.0)then
         print*,'Error getting ln3 parameters'
         return
      endif

      write(6,*)'Reading: ',filename(1:n)

      call rd_wsi_3dradar_cdf(filename,nlines,nelements,
     &rdlat,rdlon,rlat1,rlon1,validTime,
     &nlev,data_levels,num_levels,level_prefix,
     &image,jstatus)
      if(jstatus.eq.-1)then
         write(6,*)'Error reading nowrad data'
         goto 19
      end if
      write(6,*)'Found 3d WSI data. Valid Time: ',validTime

      print*
      write(6,*)' rdlat ',rdlat,' rdlon ',rdlon
      write(6,*)' rlat1 ',rlat1,' rlon1 ',rlon1
      write(6,*)' istart/iend ',istart,iend
      write(6,*)' jstart/jend ',jstart,jend
      write(6,*)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     ul=0
c     do i=1,num_levels
c       if(level_prefix(i)(1:2).eq.'ND')then
c          ndvalue=i
c       elseif(level_prefix(i)(1:1).eq.' '.or.
c    &         level_prefix(i)(1:2).eq.'ND')then
c          ul=ul+1
c          il(ul)=data_levels(i)
c       endif
c     enddo

      print*,'useable levels = ',num_levels
      write(6,*)(data_levels(i),i=1,num_levels)

      rdlat=rdlat*rdtodg
      rdlon=rdlon*rdtodg
c
c rlon1 and rlat1 already in degrees
c     rlon1=rlon1*rdtodg
c     rlat1=rlat1*rdtodg
c
c  Remap NOWrad data onto LAPS grid.  We now have the 
c  section of the NOWrad grid in units of dBZ which 
c  is in and near the LAPS domain.  This area is defined by
c  jstart, jend, istart, iend.
c
       do j=1,jmax
       do i=1,imax
          count(i,j)=0
          sum(i,j)=0
          i_max_value(i,j)=0
       enddo
       enddo

       iobcnd=0
       iobcm=0
       iobcn=0
       iobcu=0
       itot =0
       igood=0
       incnt=0
       do j=jstart,jend
          wsi_lat=rlat1-rdlat*(j)
          do i=istart,iend
             wsi_lon=rlon1+rdlon*(i)
c
c correct for signed integer
c
             if(image(i,j).lt.0)image(i,j)=256+image(i,j)

             call latlon_to_rlapsgrid(wsi_lat,wsi_lon,
     &                                lat,lon,imax,jmax,ri,rj,
     &                                jstatus)

             if(jstatus.eq.1)then

               ii=nint(ri)
               jj=nint(rj)

               itot=itot+1
               imdata=image(i,j)+1
               if(imdata.le.num_levels.and.imdata.ne.msng_radar+1)then

                   igood=igood+1
                   sum(ii,jj)=sum(ii,jj)+data_levels(imdata)
                   count(ii,jj)=count(ii,jj)+1
                   i_max_value(ii,jj)=max(i_max_value(ii,jj),
     &                                    data_levels(imdata))
               elseif(imdata.eq.msng_radar+1)then
                 iobcm=iobcm+1
               elseif(imdata.gt.num_levels.and.
     &                imdata.lt.msng_radar+1)then
                 iobcn=iobcn+1
                 incnt=imdata+incnt
               else
                 iobcu=iobcu+1 
               end if

             end if
          end do       !all i within window
       end do       !all j within window.
c
c check for echo tops or vil product
c
       if(c_data_type.ne.'et'.and.c_data_type.ne.'vi')then
          do j = 1,jmax
          do i = 1,imax
             if(count(i,j) .gt. 0)then
c               if(sum(i,j) .gt. 0)then
c               remapped_prod(i,j)=float(sum(i,j))/float(count(i,j))
                remapped_prod(i,j)=float(i_max_value(i,j))
c               else
c                  remapped_prod(i,j)=ref_base_useable
c               endif
             else
                remapped_prod(i,j)=r_missing_data
             end if
          end do
          end do
       else
          do j = 1,jmax
          do i = 1,imax
             if(i_max_value(i,j).gt.0)then
                remapped_prod(i,j)=float(i_max_value(i,j))
             else
                remapped_prod(i,j)=r_missing_data
             endif
          end do
          end do
       endif
       write(6,*)
       write(6,*)'       WSI data attributes:'
       write(6,*)' ------------------------------------'
c      write(6,*)'  Pixels No Data (ND) : ',iobcnd
c      write(6,*)'    (%): ',float(iobcnd)/float(itot)
       write(6,*)'  Pixels = Missing: ',iobcm
       write(6,*)'    (%): ',float(iobcm)/float(itot)
       write(6,*)'  Pixels > Nlevs : ',iobcn
       write(6,*)'    (%): ',float(iobcn)/float(itot)
       if(iobcn.gt.0)then
          write(6,*)'      Avg Pix Value: ',
     &              float(incnt)/float(iobcn) 
       endif
       write(6,*)'  Pixels Unaccounted for: ',iobcu
       write(6,*)'    (%): ',float(iobcu)/float(itot)
       write(6,*)'  Pixels Useable : ',igood
       write(6,*)'    (%): ',float(igood)/float(itot)
       write(6,*)' ------------------------------------'

       goto 16
19     write(6,*)'Error in nowrad_2_laps, terminating'
       goto 16
14     write(6,*)'NOWRAD data not found for given time'

16     write(6,*)'Finished in nowrad_2_laps'
       return
       end
