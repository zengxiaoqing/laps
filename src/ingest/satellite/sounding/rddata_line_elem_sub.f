      Subroutine rddata_line_elem(ncid,imax,jmax,nch
     &,sounding,istatus)
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
      Integer dim_id_x
      Integer dim_id_y
      Integer dim_id_k
      Integer i4dum

      Integer START(10)
      Integer COUNT(10)
      Character*31 DUMMY
c
C **************************************************************************
c
      istatus = 1
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

      do k = 1,nch  !# of channels in sounding database (= 19).

         write(6,*)'Reading channel ',k

         START(3)=k

         do j = 1,jmax

            COUNT(1)=imax
            START(2)=j

c read line
      rcode=NF_GET_VARA_INT(NCID,varid,START,COUNT,sounding_rec)
            if(rcode.ne.0)then 
               print*,'Error reading sounding database ',rcode
               istatus = -1
               print*,'Returning to main, istatus = ',istatus
               return
            endif
c
c load line into the 3d output array
c
            sounding(:,j,k) = sounding_rec(:)

         enddo
      enddo
c -------------------------------------
C
      write(6,*)'Sucessful reading sounding - rddata_line_elem '

      Return
      End
