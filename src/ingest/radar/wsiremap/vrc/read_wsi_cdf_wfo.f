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
       subroutine read_wsi_cdf_wfo(input_name,lines,elems,
     1        Dx,Dy,valtime,image,istatus)
c
c       This routines reads the WSI NOWRAD netcdf files and
c        converts the image value to values of 0-15.
c
c       Mark E. Jackson            17-oct-1994
c       Linda Wharton              05-apr-1996
c         modified to use C version of scan_remap
c       J Smart                       sep-1996
c         modified to read conus-c nowrad netCDF data
c
c
        include 'netcdf.inc'
        include 'vrc.inc'
        character*128 dimname                   ! Must match NETCDF.INC's
        integer ncopts
        common/ncopts/ncopts                    ! NetCDF error handling flag.

        character*200 input_name
       
        integer varid,dimid,count(3),start(3)
        real*8 valtime
        integer nrecs
        parameter (nrecs=1)
     
        integer lines,elems
	integer*1 image(nelems,nlines,nrecs)

        integer*4 bad_data_flag
        integer*4 imax_image_value
        integer*4 imin_image_value

	integer ilines,ielems
        integer istatus,cdfid
        integer attlen

        character c_atvalue*80

        integer*2 i_value
        integer*1      b_value(2), bad_data_byte
        data bad_data_byte/-1/
        equivalence (i_value,b_value(1))

        real*4 Dx,Dy
 
        istatus = 0
        bad_data_flag=255
        call ncpopt(0)

        cdfid =  NCOPN(input_name,NC_NOWRITE,istatus)
        if (istatus .ne. 0)then
            write(*,*)' Error in opening netcdf file'
            write(*,*)' ...file not found.' 
            istatus = -1
            goto 9999 
        else
                write(*,*)' netcdf nldn file open'
        endif

        dimid = NCDID(cdfid,'x',istatus)
        CALL NCDINQ(cdfid,dimid,dimname,ielems,istatus)
        dimid = NCDID(cdfid,'y',status)
        CALL NCDINQ(cdfid,dimid,dimname,ilines,istatus)
c
c A little qc check to be sure we are reading the proper data file
c
        if(ilines.lt.nlines .or. ielems.lt.nelems)then
           write(6,*)'WARNING! '
           write(6,*)'n lines from file: ',ilines
           write(6,*)'n elems from file: ',ielems
           write(6,*)'n lines expected : ',nlines
           write(6,*)'n elems expected : ',nelems
        elseif(ilines.gt.nlines .or. ielems.gt.nelems)then
           write(6,*)'TERMINAL ERROR! '
           write(6,*)'n lines from file: ',ilines
           write(6,*)'n elems from file: ',ielems
           write(6,*)'n lines expected : ',nlines
           write(6,*)'n elems expected : ',nelems
           istatus = -1
           return
        else
           write(*,*)' Getting wsi netcdf data.. '
           write(6,*)'lines/elems from netCDF: ',ilines,ielems
        endif
 
        start(1) = 1
        start(2) = 1
        count(2) = ilines
        count(1) = ielems
        start(3) = 1
        count(3) = 1

        varid = NCVID(cdfid,'image',istatus) 
        if(istatus.ne.0)then
           write(6,*)'Error getting varid - image'
           return
        endif
        CALL NCVGT(cdfid,varid,start,count,image,istatus)
        if(istatus.ne.0)then
           write(6,*)'Error reading variable - image'
           return
        endif

        varid = NCVID(cdfid,'valtime',istatus)
        if(istatus.ne.0)then
           varid = NCVID(cdfid,'validTime',istatus)
           if(istatus.ne.0)then
              write(6,*)'Error getting varid - valtime'
              return
           endif
        endif
        CALL NCVGT1(cdfid,varid,1,valtime,istatus)
        if(istatus.ne.0)then
           write(6,*)'Error reading variable - valtime'
           return
        endif

        varid = NCVID(cdfid,'Dx',istatus)
        if(istatus.ne.0)then
           write(6,*)'Error getting varid - Dx'
           return
        endif

        CALL NCVGT1(cdfid,varid,1,Dx,istatus)
        if(istatus.ne.0)then
           write(6,*)'Error reading variable - Dx'
           return
        endif

        CALL NCAINQ(cdfid,varid,'units',itype,attlen,istatus)

        call NCAGTC(cdfid,varid,'units',c_atvalue,attlen,istatus)
        if(istatus.ne.0)then
           write(6,*)'Error getting attribute - Dx'
        endif

  
        varid = NCVID(cdfid,'Dy',istatus)
        if(istatus.ne.0)then
           write(6,*)'Error getting varid - Dy'
           return
        endif
        CALL NCVGT1(cdfid,varid,1,Dy,istatus)
        if(istatus.ne.0)then
           write(6,*)'Error reading variable - Dy'
           return
        endif

        CALL NCAINQ(cdfid,varid,'units',itype,attlen,istatus)

        call NCAGTC(cdfid,varid,'units',c_atvalue,attlen,istatus)
        if(istatus.ne.0)then
           write(6,*)'Error getting attribute - Dy'
        endif

        if(c_atvalue(1:4).eq.'kilo')then
           Dx=Dx*1000.
           Dy=Dy*1000.
        endif

          Write(*,*)' closing netcdf file'
           CALL NCCLOS(cdfid,istatus) 
c
c for wfo data we must first rescale the values back to the original
c wsi form to properly convert to dbz.
c
        imax_image_value = 0
        imin_image_value = 255
        icount_bad=0
        b_value(1)=0
        do j=1,ilines
        do i=1,ielems
           b_value(2) = image(i,j,1)
           if(i_value .gt. imax_image_value) imax_image_value = i_value
           if(i_value .lt. imin_image_value) imin_image_value = i_value
           if(i_value .ne. bad_data_flag)then 
              i_value = i_value/16
              image(i,j,1)=b_value(2)
              if ((i_value .ge. 16) .or. (i_value .lt. 0)) then
                 write(6,*) i, j, i_value
              endif

           else
              icount_bad=icount_bad+1
              image(i,j,1)=bad_data_byte
           endif
        enddo
        enddo
c
        write(6,*)'Number of bad data points (> ',bad_data_flag,' )'
        write(6,*)'prior to calling c_scan_adjust: ',icount_bad
        write(6,*)'Max value found in image array: ',imax_image_value
        write(6,*)'Min value found in image array: ',imin_image_value

9999    return
        end
