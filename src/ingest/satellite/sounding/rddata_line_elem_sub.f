      Subroutine rddata_line_elem(ncid,imax,jmax,nch,
     &ndsize_x,ndsize_y,ndsize_ch,sounding,istatus)
c
c routine designed to read 3-d satellite sounding netCDF data with variable
c line and element dimensions. Should work for 2d data with variable dimensions
c if nch = 1
c
      Implicit None

      Integer*4   RCODE
c
      Integer*4   imax,jmax,nch
      Integer*4   sounding(imax, jmax, nch)
      Integer*2   sounding_rec(imax)

      Integer*4 i,j,k
      Integer*4 varid,ncid,ncdid
      Integer*4 istatus
      Integer*4 NCVID,VARID
      Integer*4 NDSIZE_X(jmax)
      Integer*4 NDSIZE_Y,NDSIZE_CH
      Integer*4 dim_id_x
      Integer*4 dim_id_y
      Integer*4 dim_id_k
      Integer*4 i4dum

      Integer*4 START(10)
      Integer*4 COUNT(10)
      Character*31 DUMMY
      integer*4     i4val
      integer*2     i2val(2)
      equivalence   (i4val,i2val(1))
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
      dim_id_y = ncdid(ncid, 'y', rcode)
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

         dim_id_k = ncdid(ncid, 'wavelength', rcode)
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
      varid = ncvid(ncid,'sounding',rcode)
      if(rcode.ne.0)then
         write(6,*)'Error getting sounding varid'
         istatus = -1
         return
      endif
c
c get x dimension id
c
      dim_id_x = ncdid(ncid, 'x', rcode)
      if(rcode.ne.0)then
         write(6,*)'Error getting x id code - returning'
         istatus = -1
         return
      endif

      do k = 1,NDSIZE_CH  !# of channels in sounding database (= 19).

         write(6,*)'Reading channel ',k

         START(3)=k

         do j = 1,NDSIZE_Y

            call ncdinq(NCID,dim_id_x,dummy,NDSIZE_X(j),RCODE)
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
            call NCVGT(NCID,varid,START,COUNT,sounding_rec,rcode)
            if(rcode.ne.0)then 
               write(6,*)'Error reading sounding database'
               istatus = -1
               write(6,*)'Returning to main, istatus = ',istatus
            endif
c
c load line into the 3d output array
c
            i4val=0
            do i = 1,NDSIZE_X(j)
               i2val(1)=0
               i2val(2)=sounding_rec(i)
               if(i4val .lt. 0)then
                  i4dum = i4val
               endif
               sounding(i,j,k) = i4val
            enddo

         enddo
      enddo
c -------------------------------------
C
      write(6,*)'Sucessful reading sounding - rddata_line_elem '

      Return
      End
