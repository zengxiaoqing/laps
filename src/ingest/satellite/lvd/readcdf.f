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
      subroutine readcdf(csat_id,csat_type,chtype,
     1record,n_elem,n_lines,r4_image,scale_img,
     1latitude,longitude,
     1La1,Lo1,Dx,Dy,Latin,
     1LoV, ivalidTime , ncid, istatus)

c     written by Dan Birkenheuer February 1994
c     J Smart   5/95    Included code to deal with gvar /public files
c     J Smart   6/96    Included code to deal with WFO satellite netcdf files
c     J Smart   3/97    Included subroutine rdblock_line_elem to read only the
c			appropriate sector of visible satellite data. Argument
c                       list lets subroutine know whether logical or I*2 data type.
c     J Smart   3/97    Converted _vis routine to be the only netCDF reader. Works
c			for vis, ir, wv and sounder. Block reading allows this.
c     J Smart   4/97    Adapted code further to work with gvar. NetCDF headers for
c                       raw gvar are different than for fsl-conus (ie., no La1, Lo1,
c                       Lov, or Latin).
c     J Smart   5/07    Added capability to read FMI (Finish Met Inst) data
c=====================================================================================
      implicit none

      Include 'netcdf.inc'

      Integer n_elem,n_lines,record
      INTEGER RCODE
C
      real      r4_image(n_elem,n_lines)
      real      latitude(n_elem,n_lines)
      real      longitude(n_elem,n_lines)
      integer  ib
      integer       bi (2)
      equivalence (ib,bi (1) )
      integer istatus
      Integer   iDx
      Integer   iD
      INTEGER   Nx
      INTEGER   Ny
      INTEGER   Nz
      REAL      La1     
      REAL      Lo1     
      REAL      Dx      
      REAL      Dy      
      REAL      Latin   
      REAL      LoV     
      REAL      poLat
      REAL      La2,Lo2,level
      REAL      scale_img
      real arg1,arg2
      INTEGER START(10)
      INTEGER COUNT(10)
      integer varid,ncid
      integer center_id, process_id,
     +     wmo_sat_id(record),ivalidtime,nav
      Integer imax,jmax,kmax,kdim
      Integer channel_fcinv

      integer lineRes, elemRes
      common /cdfdata/ lineRes, elemRes

      double precision reftime(record), valtime(record)
      character*30  c_valtime
      character*30  c_Lov
      character*30  c_xres
      character*30  c_yres
      Character*3   csat_type
      character*3   chtype
      character*6   csat_id
      character*132 origin_name
      character*132 x_dim
      character*132 y_dim
      character*132 earth_shape
      character*132 wavelength(record)
      character*132 grid_name
      character*132 process_name
      character*132 grid_type
      character*132 channel_comment_
      character*132 asctime
c---------------------------------------------------------
c   code

      print*,'Subroutine readcdf...'
      print*,'n_elem/n_lines = ',n_elem,n_lines

      istatus = 0 ! bad istatus

      call NCPOPT(0)

      if(csat_type.eq.'rll')then ! read lat/lon
         rcode = NF_INQ_VARID(ncid,'latitude',varid)
         if(rcode.ne.NF_NOERR) then
            print *, NF_STRERROR(rcode)
            print *,'in var latitude'
         else
            write(6,*)
            Write(6,*)'Calling rdblock_line_elem - latitude'

            Call rdblock_line_elem(csat_id,csat_type,chtype,
     &          ncid,varid,n_elem,n_lines,latitude,istatus)

            if(istatus .ne. 1)then
               write(6,*)'Error in rdblock_line_elem'
               return
            endif
         endif
         call clean_nan2(latitude,n_elem,n_lines,istatus)
         arg1 = minval(latitude)
         arg2 = maxval(latitude)
         write(6,*)'readcdf latitude range:       ',arg1,arg2
         if(arg1 .eq. 0. .and. arg2 .eq. 0.)then
            write(6,*)' ERROR in latitude range'
            istatus = 0
            return
         endif

         rcode = NF_INQ_VARID(ncid,'longitude',varid)
         if(rcode.ne.NF_NOERR) then
            print *, NF_STRERROR(rcode)
            print *,'in var longitude'
         else
            write(6,*)
            Write(6,*)'Calling rdblock_line_elem - longitude'

            Call rdblock_line_elem(csat_id,csat_type,chtype,
     &          ncid,varid,n_elem,n_lines,longitude,istatus)

            if(istatus .ne. 1)then
               write(6,*)'Error in rdblock_line_elem'
               return
            endif
         endif
         call clean_nan2(longitude,n_elem,n_lines,istatus)
         write(6,*)'readcdf longitude range:      '
     1          ,minval(longitude),maxval(longitude)

      endif

      if(csat_type.eq.'ncp')then
         rcode = NF_INQ_VARID(ncid,'channel',varid)
      else
         rcode = NF_INQ_VARID(ncid,'image',varid)
      endif
      if(rcode.ne.NF_NOERR) then
         print *, NF_STRERROR(rcode)
         print *,'in var image'
      endif
C
C    statements to fill image                          
C
      write(6,*)
      Write(6,*)'Calling rdblock_line_elem - image read sub'

      Call rdblock_line_elem(csat_id,csat_type,chtype,
     &ncid,varid,n_elem,n_lines,r4_image,istatus)

      if(istatus .ne. 1)then
         write(6,*)'Error in rdblock_line_elem'
         return
      endif

      if(csat_type.eq.'rll')then 

!         read scaling attribute
          rcode=NF_GET_ATT_REAL(ncid,varid,'scale_factor'
     1                              ,scale_img)
          if(rcode.ne.NF_NOERR) then
             write(6,*)'Error reading image scaling attribute'
             scale_img = .01
             write(6,*)' Use default value for image scale ',scale_img
          else
             write(6,*)' Successfully read image scale ',scale_img
          endif

          rcode = NF_INQ_VARID(ncid,'lineRes',varid)
          if(rcode.ne.NF_NOERR) then
             print *, NF_STRERROR(rcode)
             print *,'in var lineRes'
          endif
          rcode = NF_GET_VAR_INT(ncid,varid,lineRes)
          if(rcode.ne.NF_NOERR) then
             print *, NF_STRERROR(rcode)
             print *,'in var lineRes'
          endif
          if(rcode.ne.NF_NOERR) then
             write(6,*)'Error reading lineRes'
          else
             write(6,*)' Successfully read lineRes ',lineRes
          endif

          rcode = NF_INQ_VARID(ncid,'elemRes',varid)
          if(rcode.ne.NF_NOERR) then
             print *, NF_STRERROR(rcode)
             print *,'in var elemRes'
          endif
          rcode = NF_GET_VAR_INT(ncid,varid,elemRes)
          if(rcode.ne.NF_NOERR) then
             print *, NF_STRERROR(rcode)
             print *,'in var elemRes'
          endif
          if(rcode.ne.NF_NOERR) then
             write(6,*)'Error reading elemRes'
          else
             write(6,*)' Successfully read elemRes ',elemRes
          endif
C
      endif

      istatus = 0
c
c JSmart: 6-4-97.  WFO satellite netCDF files changed rather dramatically
c                  such that NO header info exists and we must now jump over
c                  the statements to read the header info. 
c
      if(csat_type.eq.'cdf')then
         call read_netcdf_sat_head(ncid,record,
     + Nx, Ny, center_id,process_id, wmo_sat_id,Dx,Dy,
     + La1, Latin, Lo1, Lov, reftime, valtime, earth_shape,
     + grid_name, grid_type, origin_name, process_name,
     + wavelength, x_dim, y_dim, istatus)

         ivalidtime = int(valtime(1))
      endif

      write(6,*)'readcdf image range: '
     1          ,minval(r4_image),maxval(r4_image)

      istatus = 1  ! ok!

      Return
      END
