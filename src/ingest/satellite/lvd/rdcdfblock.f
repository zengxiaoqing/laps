      Subroutine rdblock_line_elem(csat_id,csat_type,chtype,
     &ncid,varid,n_elems,n_lines,nch,ndsize_x,ndsize_y,
     &data,istatus)
c
c routine designed to read 3-d satellite sounding and image netCDF data with variable
c line and element dimensions.
c if nch = 1 then it is currently assumed to be image data.
c
      Implicit None

      include 'netcdf.inc'
      include 'satellite_dims_lvd.inc'

      integer      RCODE
c
      Integer n_elems,n_lines,nch
      real    data(n_elems, n_lines)
      integer ndsize_x,ndsize_y
c     integer image(ndsize_x,ndsize_y,1)

      Integer n,nn
      Integer varid,ncid
      Integer istatus
      Integer NDSIZE_CH
      Integer dim_id_x
      Integer dim_id_y
      Integer dim_id_k
      integer istart
      integer iend
      integer jstart
      integer jend

      Integer START(10)
      Integer COUNT(10)
      Character*31 DUMMY
      Character*40 c_data_type
      Character*6  csat_id
      Character*3  csat_type
      Character*3  chtype
      integer i,j,ii,jj 
c
C **************************************************************************
c
      istatus = -1
c
c get starting and ending line/pixel values for the satellite data in question
c
      nn=index(chtype,' ')-1
      if(nn.le.0)nn=3

      call getsat_attributes(csat_id,csat_type,chtype,
     &istart,iend,jstart,jend,istatus)
      if(istatus.ne.1)goto 1000
      istatus = -1
c
c Code to get dimension size and read individual element of sounding array
c get dimensions for sounding array (x,y,lambda) [lambda is # of wavelengths]
c
c This is the number of lines
c
c     dim_id_y = NCDID(ncid, 'y', rcode)
c     if(rcode.ne.0)then
c        write(6,*)'Error getting y id code - returning'
c        return
c     endif
c     CALL NCDINQ(NCID, dim_id_y,dummy,NDSIZE_Y,RCODE)
c     if(rcode.ne.0)then
c        write(6,*)'Error getting y dimension size'
c        return
c     endif
c
c get 3rd dimension if it exists. For sounding data this is the number of wavelengths
c set the data variable character type
c
c     if(chtype .eq.'snd')then

c        c_data_type = 'sounding'
c        n=index(c_data_type,' ')-1

c        dim_id_k = NCDID(ncid, 'wavelength', rcode)
c        if(rcode.ne.0)then
c           write(6,*)'Error getting channel id code - returning'
c           return
c        endif

c        CALL NCDINQ(NCID, dim_id_k,dummy,NDSIZE_CH,RCODE)
c        if(rcode.ne.0)then
c           write(6,*)'Error getting channel dimension size'
c           return
c        endif

c        write(6,*)'Reading sounder channel ',nch
c        START(3)=nch
c        COUNT(3)=1

c     else   !assume it is image data. 2D data.
        
c        c_data_type = 'image'
c        n=index(c_data_type,' ')-1
c     endif
c
c Now ready to read the data
c
c get variable id. This is somewhat redundant as the calling program
c does the same. Consider this a safety check.
c
c     rcode=NF_INQ_VARID(ncid,c_data_type(1:n),varid)
c     if(rcode.ne.0)then
c        write(6,*)'Error getting varid ',c_data_type(1:n)
c        return
c     endif
c
c get x dimension id
c
c     dim_id_x = NCDID(ncid, 'x', rcode)
c     if(rcode.ne.0)then
c        write(6,*)'Error getting x id code - returning'
c        return
c     endif
c     call NCDINQ(NCID,dim_id_x,dummy,NDSIZE_X,RCODE)
c     if(rcode.ne.0)then
c        write(6,*)'Error getting x dimension - NDSIZE_X'
c     endif
c
c qc step
c -------
      if(jend.gt.NDSIZE_Y)then
         write(6,*)'jend > NDSIZE_Y - rdcdfblock.f '
         write(6,*)'Returing: istatus = ',istatus
         return
      endif

      if(iend.gt.NDSIZE_X)then
         write(6,*)'iend > NDSIZE_X '
         write(6,*)'Returning to main, istatus = ',istatus
         return
      endif

      START(1)=istart
      COUNT(1)=iend-istart+1
      START(2)=jstart
      COUNT(2)=jend-jstart+1
      if(csat_type.eq.'cdf'.or.csat_type.eq.'wfo')then
         start(3)=1
         count(3)=1
      endif
c
c read line. Switch here discriminates 1-byte versus 2-byte data.
c

      rcode=NF_GET_VARA_REAL(NCID,varid,START,COUNT,data)

      do j=1,n_lines
      do i=1,n_elems
         if(data(i,j).lt.0)data(i,j)=data(i,j)+256.
      enddo
      enddo

c     rcode=nf_get_vara_int(ncid,varid,start,count,image)
c     call main_sub(ncid, 1, ndsize_x, ndsize_y, image, rcode)
c     if(rcode.ne.0)then 
c        write(6,*)'Error reading database'
c        write(6,*)'Returning to main, istatus = ',istatus
c     endif
c
c cdf and wfo data types are unsigned bytes.
c
c     i=0
c     j=0
c     if(csat_type.eq.'cdf'.or.csat_type.eq.'wfo') then
c        do jj=jstart,jend
c           j=j+1
c           do ii=istart,iend
c              i=i+1
c              if(image(ii,jj,1).lt.0)
c    &            image(ii,jj,1)=256+image(ii,jj,1)
c              data(i,j)=float(image(ii,jj,1))
c           enddo
c        enddo
c     endif
c
c -------------------------------------
C
      istatus = 1
1000  Return
      End
