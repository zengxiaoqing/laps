      Subroutine rdblock_line_elem(csat_id,csat_type,chtype,
     &ncid,n_elems,n_lines,nch,ndsize_x,ndsize_y,ndsize_ch,
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
      Integer   n_elems,n_lines,nch
      real   data(n_elems, n_lines)

      Integer n,nn
      Integer varid,ncid
      Integer istatus
      Integer NDSIZE_X
      Integer NDSIZE_Y
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
      integer i,j 
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
      dim_id_y = NCDID(ncid, 'y', rcode)
      if(rcode.ne.0)then
         write(6,*)'Error getting y id code - returning'
         return
      endif
      CALL NCDINQ(NCID, dim_id_y,dummy,NDSIZE_Y,RCODE)
      if(rcode.ne.0)then
         write(6,*)'Error getting y dimension size'
         return
      endif
c
c qc step
c -------
      if(jend.gt.NDSIZE_Y)then
         write(6,*)'jend > NDSIZE_Y - rdcdfblock.f '
         write(6,*)'Returing: istatus = ',istatus
         return
      endif
c
c get 3rd dimension if it exists. For sounding data this is the number of wavelengths
c set the data variable character type
c
      if(chtype .eq.'snd')then

         c_data_type = 'sounding'
         n=index(c_data_type,' ')-1

         dim_id_k = NCDID(ncid, 'wavelength', rcode)
         if(rcode.ne.0)then
            write(6,*)'Error getting channel id code - returning'
            return
         endif

         CALL NCDINQ(NCID, dim_id_k,dummy,NDSIZE_CH,RCODE)
         if(rcode.ne.0)then
            write(6,*)'Error getting channel dimension size'
            return
         endif
c qc step
c -------
         if(nch.gt.NDSIZE_CH)then
            write(6,*)'nch > NDSIZE_CH - rdcdfblock!'
            write(6,*)'Returing to main: istatus = ',istatus
            return
         endif

         write(6,*)'Reading sounder channel ',nch
         START(3)=nch
         COUNT(3)=1

      else   !assume it is image data. 2D data.
        
         c_data_type = 'image'
         n=index(c_data_type,' ')-1

      endif
c
c Now ready to read the data
c
c get variable id. This is somewhat redundant as the calling program
c does the same. Consider this a safety check.
c
      rcode=NF_INQ_VARID(ncid,c_data_type(1:n),varid)
      if(rcode.ne.0)then
         write(6,*)'Error getting varid ',c_data_type(1:n)
         return
      endif
c
c get x dimension id
c
      dim_id_x = NCDID(ncid, 'x', rcode)
      if(rcode.ne.0)then
         write(6,*)'Error getting x id code - returning'
         return
      endif

      call NCDINQ(NCID,dim_id_x,dummy,NDSIZE_X,RCODE)
      if(rcode.ne.0)then
         write(6,*)'Error getting x dimension - NDSIZE_X'
      endif
      if(iend.gt.NDSIZE_X.or.NDSIZE_X.le.0)then
         write(6,*)'iend > NDSIZE_X or iend <= 0'
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

      if(rcode.ne.0)then 
         write(6,*)'Error reading database'
         write(6,*)'Returning to main, istatus = ',istatus
      endif
c
c cdf and wfo data types are unsigned bytes.
c
      if(csat_type.eq.'cdf'.or.csat_type.eq.'wfo') then
         do j=1,n_lines
            do i=1,n_elems
               if(data(i,j).lt.0) data(i,j)=256+data(i,j)
            enddo
         enddo
      endif

c
c -------------------------------------
C
      istatus = 1
1000  Return
      End
