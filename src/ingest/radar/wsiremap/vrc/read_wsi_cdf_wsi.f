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
	subroutine read_wsi_cdf_wsi(input_name,nlines,nelems,dlat,
     1        dlon,lat2,lon1,validTime,Dx,Dy,image,status)

        
c
c       This routines reads the WSI NOWRAD netcdf files and
c        converts the image value to values of 0-15.
c
c       Mark E. Jackson            17-oct-1994
c       Linda Wharton              05-apr-1996
c         modified to use C version of scan_remap
c
c
        include 'netcdf.inc'
        character*128 dimname                   ! Must match NETCDF.INC's

        common/ncopts/ncopts                    ! NetCDF error handling flag.

        character*200 input_name
        character c_atvalue*80

        real*4 Dx,Dy
       
        integer varid,dimid,count(2),start(2),icstart(2),
     1      icend(2)
     
        integer*4 nelems,nlines
	integer*1 image(nelems,nlines)
  
        integer*2 i_value
        integer*1 b_value(2) 

        INTEGER*1 bad_data_byte
        data bad_data_byte /-1/
        equivalence (i_value,b_value(1))

	integer validTime,lines,elems
	real centerLon,topLat,diffLon,radsPerLine,radsPerElem
        integer status,cdfid
        integer  bad_data_flag
 
        real lat2,lon1,dlon,dlat
 
        status = 0
        call ncpopt(0)

        cdfid =  NCOPN(input_name,NC_NOWRITE,status)
        if (status .ne. 0)then
	   write(*,*)' Error in opening netcdf file'
	   write(*,*)' ...file not found.' 
	   goto 9999 
        else
	   write(*,*)' netcdf nldn file ',input_name,' open'
        endif

        dimid = NCDID(cdfid,'lines',status)
        CALL NCDINQ(cdfid,dimid,dimname,lines,status)
        dimid = NCDID(cdfid,'elems',status)
        CALL NCDINQ(cdfid,dimid,dimname,elems,status)

        start(1) = 1
        start(2) = 1
        count(2) = lines
        count(1) = elems

        icstart(1) = 1
        icend(1) = 50
        icstart(2) = 1
        icend(2) = 1 

        write(*,*)' Getting wsi netcdf data.. '
c        write(*,*)'lines,elems',lines,elems

        varid = NCVID(cdfid,'image',status) 
        CALL NCVGT(cdfid,varid,start,count,image,status)

          varid = NCVID(cdfid,'validTime',status)
          CALL NCVGT1(cdfid,varid,1,validTime,status)

        varid = NCVID(cdfid,'centerLon',status)
        CALL NCVGT1(cdfid,varid,1,centerLon,status)
  
          varid = NCVID(cdfid,'topLat',status)
          CALL NCVGT1(cdfid,varid,1,topLat,status)

        varid = NCVID(cdfid,'diffLon',status)
        CALL NCVGT1(cdfid,varid,1,diffLon,status)
         
          varid = NCVID(cdfid,'radsPerLine',status)
          CALL NCVGT1(cdfid,varid,1,radsPerLine,status)

        varid = NCVID(cdfid,'radsPerElem',status)
        CALL NCVGT1(cdfid,varid,1,radsPerElem,status)
       
        varid = NCVID(cdfid,'Dx',status)
        if(status.ne.0)then
           write(6,*)'Error getting varid - Dx'
           return
        endif

        CALL NCVGT1(cdfid,varid,1,Dx,status)
        if(status.ne.0)then
           write(6,*)'Error reading variable - Dx'
           return
        endif

        CALL NCAINQ(cdfid,varid,'units',itype,attlen,status)

        call NCAGTC(cdfid,varid,'units',c_atvalue,attlen,status)
        if(status.ne.0)then
           write(6,*)'Error getting attribute - Dx'
        endif

  
        varid = NCVID(cdfid,'Dy',status)
        if(status.ne.0)then
           write(6,*)'Error getting varid - Dy'
           return
        endif
        CALL NCVGT1(cdfid,varid,1,Dy,status)
        if(status.ne.0)then
           write(6,*)'Error reading variable - Dy'
           return
        endif

        CALL NCAINQ(cdfid,varid,'units',itype,attlen,status)

        call NCAGTC(cdfid,varid,'units',c_atvalue,attlen,status)
        if(status.ne.0)then
           write(6,*)'Error getting attribute - Dy'
        endif

        if(c_atvalue(1:4).eq.'degr')then
           Dx=Dx*111.1*1000.
           Dy=Dy*111.1*1000.
        endif

          Write(*,*)' closing netcdf file'
           CALL NCCLOS(cdfid,status) 

        dlat = radsPerLine
        lat2 = topLat
        dlon = radsPerElem
        lon1 = -abs(centerLon) - dlon*(elems-1)/2

c       do j = 1,lines
c       do i=1,elems
c          if(image(i,j) .gt. 15)then
c           scanin = image(i,j)
c           call scan_adjust(scanin,scanout)
c           image(i,j) = scanout
c          endif
c       enddo     !all elems/longitudes along a line/latitude in nowrad grid.
c       enddo     !all lines/latitudes in the nowrad grid.

        bad_data_flag=255
        imax_image_value = 0
        imin_image_value = 255
        icount_bad=0
        icount_out=0
        b_value(1)=0
        do j=1,lines
        do i=1,elems
           b_value(2) = image(i,j)
           if(i_value .gt. imax_image_value) imax_image_value = i_value
           if(i_value .lt. imin_image_value) imin_image_value = i_value
           if(i_value .ne. bad_data_flag)then
              image(i,j)=b_value(2)
              if ((i_value .gt. 64) .or. (i_value .lt. 0)) then
                 icount_out = icount_out +1
c                write(6,*) i, j, i_value
              endif

           else
              icount_bad=icount_bad+1
              image(i,j) = bad_data_byte
           endif
        enddo
        enddo
c
        write(6,*)'Number of bad data points (> ',bad_data_flag,' )'
        write(6,*)'prior to calling c_scan_adjust: ',icount_bad
        write(6,*)'Max value found in image array: ',imax_image_value
        write(6,*)'Min value found in image array: ',imin_image_value
        write(6,*)
        write(6,*)'Data found out-of-bounds (icount_out) ',icount_out
        if(icount_out.gt.0)then
           status = -1
           return
        endif


C       new C version of scan_adjust implemented 05-apr-96
	call c_scan_adjust(image,lines,elems,bad_data_flag)

9999    return
        end
