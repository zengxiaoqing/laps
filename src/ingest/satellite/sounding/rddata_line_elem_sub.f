      Subroutine rddata_line_elem(ncid,imax,jmax,nch,
     &ndsize_x,ndsize_y,ndsize_ch,sounding,istatus)
c
c routine designed to read 3-d satellite sounding netCDF data with variable
c line and element dimensions. Should work for 2d data with variable dimensions
c if nch = 1
c
      Implicit None
      include 'netcdf.inc'
      Integer   RCODE
c
      Integer   imax,jmax,nch
      Integer   sounding(imax, jmax, nch)
      Integer   sounding_rec(imax)

      Integer i,j,k
      Integer varid,ncid
      Integer istatus
*4*4=NF_INQ_VARID,VARID
      Integer NDSIZE_X(jmax)
      Integer NDSIZE_Y,NDSIZE_CH
      Integer dim_id_x
      Integer dim_id_y
      Integer dim_id_k
      Integer i4dum

      Integer START(10)
      Integer COUNT(10)
      Character*31 DUMMY
c     integer     i4val
c     integer     i2val(2)
c     equivalence   (i4val,i2val(1))
c
C **************************************************************************
c
      istatus = 1
c
c Code to get dimension size and read individual element of sounding array
c get dimensions for sounding array (x,y,lambda) [lambda is # of wavelengths]
c
c This is the number of lines
c
      dim_id_y = NCDID(ncid, 'y', rcode)
      if(rcode.ne.0)then
         write(6,*)'Error getting y id code - returning'
         istatus = -1
         return
      endif
      CALL NCDINQ(NCID, dim_id_y,dummy,NDSIZE_Y,RCODE)
      if(rcode.ne.0)then
         write(6,*)'Error getting y dimension size'
         istatus = -1
         return
      endif
c
c qc step
c -------
      if(NDSIZE_Y.gt.jmax)then
         write(6,*)'y dim size > n_lines!'
         istatus = -1
         write(6,*)'Returing to mail: istatus = ',istatus
         return
      endif
c
c get 3rd dimension if it exists. For sounding data this is the number of wavelengths
c
      if(nch .gt. 1)then

         dim_id_k = NCDID(ncid, 'wavelength', rcode)
         if(rcode.ne.0)then
            write(6,*)'Error getting channel id code - returning'
            istatus = -1
            return
         endif

         CALL NCDINQ(NCID, dim_id_k,dummy,NDSIZE_CH,RCODE)
         if(rcode.ne.0)then
            write(6,*)'Error getting channel dimension size'
            istatus = -1
            return
         endif
c qc step
c -------
         if(NDSIZE_CH.gt.nch)then
            write(6,*)'Channel dim size > nchannels!'
            istatus = -1
            write(6,*)'Returing to mail: istatus = ',istatus
            return
         endif

      endif
c
c Now ready to read the data
c
      START(1)=1
      COUNT(2)=1
      COUNT(3)=1
c
c get sounding variable id
c
      rcode=NF_INQ_VARID(ncid,'sounding',varid)
      if(rcode.ne.0)then
         write(6,*)'Error getting sounding varid'
         istatus = -1
         return
      endif
c
c get x dimension id
c
      dim_id_x = NCDID(ncid, 'x', rcode)
      if(rcode.ne.0)then
         write(6,*)'Error getting x id code - returning'
         istatus = -1
         return
      endif

      do k = 1,NDSIZE_CH  !# of channels in sounding database (= 19).

         write(6,*)'Reading channel ',k

         START(3)=k

         do j = 1,NDSIZE_Y

            call NCDINQ(NCID,dim_id_x,dummy,NDSIZE_X(j),RCODE)
            if(rcode.ne.0)then
               write(6,*)'Error getting x dimension - NDSIZE_X'
            endif
            if(NDSIZE_X(j).gt.imax.or.NDSIZE_X(j).le.0)then
               write(6,*)'x dimension size > nlines_max or <= 0'
               istatus = -1
               write(6,*)'Returning to main, istatus = ',istatus
               return
            endif

            COUNT(1)=NDSIZE_X(j)
            START(2)=j

c read line
      rcode=NF_GET_VARA_INT(NCID,varid,START,COUNT,sounding_rec)
            if(rcode.ne.0)then 
               write(6,*)'Error reading sounding database'
               istatus = -1
               write(6,*)'Returning to main, istatus = ',istatus
            endif
c
c load line into the 3d output array
c No more need to use the equivalence to get sounder values. (1-14-98. JRS).
c           i4val=0
            do i = 1,NDSIZE_X(j)
c              i2val(1)=0
c              i2val(2)=sounding_rec(i)
c              if(i4val .lt. 0)then
c                 i4dum = i4val
c              endif
               sounding(i,j,k) = sounding_rec(i)
               if(sounding_rec(i).ge.0)then
c                 sounding(i,j,k) = i4val
               endif
            enddo

         enddo
      enddo
c -------------------------------------
C
      write(6,*)'Sucessful reading sounding - rddata_line_elem '

      Return
      End
