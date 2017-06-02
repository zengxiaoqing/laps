      Subroutine rdblock_line_elem(csat_id,csat_type,chtype,
     &ncid,varid,n_elems,n_lines,data,istatus)
c
c routine designed to read 3-d satellite sounding and image netCDF data with variable
c line and element dimensions.
c
      Implicit None

      include 'netcdf.inc'
      include 'satellite_dims_lvd.inc'
      include 'satellite_common_lvd.inc'

      integer      RCODE
c
      Integer n_elems,n_lines,nch
      real    data(n_elems, n_lines)
      Integer*2 data_int2(n_elems, n_lines)

      Integer n,nn
      Integer varid,ncid
      Integer istatus
      Integer NDSIZE_CH
      Integer NDSIZE_X
      Integer NDSIZE_Y
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
     &istart,iend,jstart,jend,ndsize_x,ndsize_y,istatus)
      if(istatus.ne.1)goto 1000
      istatus = -1
c
c qc step
c -------
      if(jend.gt.NDSIZE_Y)then
         write(6,*)'jend > NDSIZE_Y - rdcdfblock.f ',jend,NDSIZE_Y
         write(6,*)'Returing: istatus = ',istatus
         return
      endif

      if(iend.gt.NDSIZE_X)then
         write(6,*)'iend > NDSIZE_X ',iend,NDSIZE_X
         write(6,*)'Returning to main, istatus = ',istatus
         return
      endif

      START(1)=istart
      COUNT(1)=iend-istart+1
      START(2)=jstart
      COUNT(2)=jend-jstart+1
      if(csat_type.eq.'cdf'.or.
     +   csat_type.eq.'wfo')then
         start(3)=1
         count(3)=1
      elseif(csat_type.eq.'ncp')then
         start(3)=1
         count(3)=1
         start(4)=1
         count(4)=1
      endif
c
c read line. Switch here discriminates 1-byte versus 2-byte data.
c
      if(csat_type.eq.'rll' .or. csat_type.eq.'gnp')then
          write(6,*)'ncid/varid = ',ncid,varid
          write(6,*)'start=',start(1:2),' count=',count(1:2)
          rcode=NF_GET_VARA_INT2(NCID,varid,START,COUNT,data_int2)
          write(6,*)'center pixel i2: ',data_int2(n_elems/2,n_lines/2)
          write(6,*)'rdblock_line_elem i2 data range: '
     1              ,minval(data_int2),maxval(data_int2)
          data(:,:) = data_int2(:,:)
      else
          rcode=NF_GET_VARA_REAL(NCID,varid,START,COUNT,data)
      endif

      write(6,*)'rdblock_line_elem data range: '
     1          ,minval(data),maxval(data)

      if(csat_type.ne.'rll')then
          do j=1,n_lines
          do i=1,n_elems
              if(data(i,j).lt.0)data(i,j)=data(i,j)+256.
          enddo
          enddo
      endif

      write(6,*)'center pixel r4: ',data(n_elems/2,n_lines/2)
C
      istatus = 1
1000  Return
      End
