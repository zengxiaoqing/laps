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
      subroutine readcdf(csat_id,csat_type,chtype,Nx,Ny,
     1record,n_elem,n_lines,r4_image,La1,Lo1,Dx,Dy,Latin,
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
c================================================================
      implicit none

      Include 'netcdf.inc'

      Integer n_elem,n_lines,record
      INTEGER RCODE
C
      real*4      r4_image(n_elem,n_lines)
      integer  ib
      integer       bi (2)
      equivalence (ib,bi (1) )
      character*30 c_valtime
      character*30 c_Lov
      character*30 c_xres
      character*30 c_yres
      Character*3  csat_type
      character*3  chtype
      character*6  csat_id
      integer istatus
      Integer   iDx
      Integer   iD
      INTEGER   Nx
      INTEGER   Ny
      REAL*4    La1     (record)
      REAL*4    Lo1     (record)
      REAL*4    Dx      (record)
      REAL*4    Dy      (record)
      REAL*4    Latin   (record)
      REAL*4    LoV     (record)
      INTEGER START(10)
      INTEGER COUNT(10)
      integer varid,ncid
      integer center_id, process_id,
     +     wmo_sat_id(record),ivalidtime
      double precision reftime(record), valtime(record)
      character*132 origin_name
      character*132 x_dim
      character*132 y_dim
      character*132 earth_shape
      character*132 wavelength(record)
      character*132 grid_name
      character*132 process_name
      character*132 grid_type
c---------------------------------------------------------
c   code

      print*,'In routine readcdf'
      print*,'n_elem/n_lines = ',n_elem,n_lines

      istatus = 0 ! bad istatus

      call NCPOPT(0)

      rcode = NF_INQ_VARID(ncid,'image',varid)
      if(rcode.ne.NF_NOERR) then
         print *, NF_STRERROR(rcode)
         print *,'in var image'
      endif
C
C    statements to fill image                          
C
      Write(6,*)'Calling rdblock_line_elem - image read sub'

      Call rdblock_line_elem(csat_id,csat_type,chtype,
     &ncid,varid,n_elem,n_lines,1,nx,ny,
     &r4_image,istatus)

      if(istatus .ne. 1)then
         write(6,*)'Error in rdblock_line_elem'
         return
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

      istatus = 1  ! ok!

      rcode= NF_CLOSE(ncid)

      Return
      END
