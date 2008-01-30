      Subroutine Process_SDHS_GVAR_sub(filename,chtype,isat,jtype,
     &nelem,nlines,image_out,nlfi,nefi,depth,width,istatus)

C**********************************************************************
C  Purpose:  Convert the SDHS image data from a cellular format to an 
C  image format
C
C  Method:  read the image header information to determine how many
C  64 scanline by 256 pixel cells make up the image.  Read in the image
C  into a singlely dimnesioned array and convert the data into a two
C  dimensional array where the first index indicates the row (scanline
C  number) and the second dimension the column(pixel number).  Some 
C  Fortran 90 constructs are used to allow dynamic array sizing
C  but the majority of the code is in Fortran 77. 
C  No attempt is made to remove the blank parts of cells
C  which may be in the cells send to satisfy an image request.  Also
C  since AIX fortran will not allow an unsigned byte variable the data
C  is read into an INTEGER*1 variable and the bit string picked off and
C  stored into an INTEGER*2 variable.  Without using this method byte 
C  values such as 255 (all bits on) were being interpreted as -128.
C  
C  References:
C  1.  Final Interface Specification for the Satellite Data Handling
C  System Communications Network (SDHS-COMNET), 01 February 1988
C  2.  AFGWC/DONS PV-Wave program auto_convert.pro
C
C**********************************************************************

      IMPLICIT NONE

      CHARACTER*255 FILENAME
      Character*3   chtype

      INTEGER CELL_WIDTH !Number of pixels in width of a cell
      INTEGER CELL_DEPTH !Number of scanline in depth of cell
      INTEGER HEADER_SIZE  !size of GVAR header
      INTEGER PIXELS_PER_CELL !Number of Pixels per cell
      Integer nelem,nlines
      Integer nefi,nlfi
      Integer imagelen
      Integer nf
      Integer lend
      INTEGER WIDTH      !Number of Cells in Width of image
      INTEGER DEPTH      !Number of Cells in Depth of image
      Integer itempintgr

      PARAMETER (HEADER_SIZE = 256) 
      PARAMETER (CELL_WIDTH = 256)
      PARAMETER (CELL_DEPTH = 64)
      PARAMETER (PIXELS_PER_CELL = CELL_WIDTH*CELL_DEPTH)

c     PARAMETER (WIDTH=5,DEPTH=13)  !this for the test images so far.
c     Parameter (imagelen=width*cell_width*depth*cell_depth)

c     Parameter (imagelen_ir=nlines_ir_max*nelem_ir_max)
c     Parameter (imagelen_vis=nlines_vis_max*nelem_vis_max)
c     Parameter (imagelen_wv=nlines_wv_max*nelem_wv_max)

      INTEGER SIZE       !Size of the SDHS image including the header
c     INTEGER CELL_NUM   !The number of the cell
c     INTEGER CELL_ROW   !The row number of a cell
c     INTEGER CELL_COLUMN !The column number of a cell
c     INTEGER XS         !Pixels in width of image
c     INTEGER YS         !Scanlines depth of image

      integer HEADER(HEADER_SIZE)!array used to read over the header
      INTEGER IOSTATUS   !I/O status variable
      Integer istatus    !Return status
      INTEGER I,J        !DO loop indices
      INTEGER RNUM       !the row number (scanline) of the converted
                           !image
      INTEGER CNUM       !the column number(pixel) of the converted 
                           !image
      Integer ii,jj,m
      Integer ispec
      Integer istat
      Integer istart,jstart
      Integer iend,jend
      Integer isat,jtype

C**********************************************************************
C  Fortran 90 used to declare allocatable arrays 
C**********************************************************************

c     INTEGER*1,DIMENSION(:),ALLOCATABLE:: IMAGEI1 
c     INTEGER*2,DIMENSION(:,:),ALLOCATABLE:: IMAGE 
c static case: JSmart 5-12-97
c     Integer*1 IMAGEI1(WIDTH*CELL_WIDTH*DEPTH*CELL_DEPTH)

      character      IMAGEC(nlfi*nefi)*1
      Integer        image_decell_2d(nefi,nlfi)
      Integer        image_temp(nefi,nlfi)
      Integer        image_decell_1d(nefi*nlfi+HEADER_SIZE*4)
      real           image_out(nelem,nlines)
      real           r8to10

      include 'satellite_dims_lvd.inc'
      include 'satellite_common_lvd.inc'

      istatus=0

      call lvd_file_specifier(chtype,ispec,istat)
      if(ispec.eq.1)then
       istart=i_start_vis(jtype,isat)
       jstart=j_start_vis(jtype,isat)
       iend=i_end_vis(jtype,isat)
       jend=j_end_vis(jtype,isat)
      elseif(ispec.eq.3)then
       istart=i_start_wv(jtype,isat)
       jstart=j_start_wv(jtype,isat)
       iend=i_end_wv(jtype,isat)
       jend=j_end_wv(jtype,isat)
      elseif(ispec.eq.2.or.ispec.eq.4.or.ispec.eq.5)then
       istart=i_start_ir(jtype,isat)
       jstart=j_start_ir(jtype,isat)
       iend=i_end_ir(jtype,isat)
       jend=j_end_ir(jtype,isat)
      endif
C
C**********************************************************************
C  Read the header information and pass back the number of 64 scanline
C  by 256 pixel cells in the width and depth of the image.
C**********************************************************************

C     CALL READHDR(FILENAME,WIDTH, DEPTH)

C**********************************************************************
C  compute the number of pixels in the width of the image and the 
C  number of scanline in the depth of the image.  Also compute the size
C  of the image file including the header. 
C**********************************************************************
C Original constructs:
C     XS =WIDTH*CELL_WIDTH
C     YS = DEPTH*CELL_DEPTH
C     SIZE = YS*XS+HEADER_SIZE

      SIZE = nlfi*nefi+HEADER_SIZE*4   !note HEADER_SIZE = 256*4=1024
C
C**********************************************************************
C  Open the file that contains the GVAR header and pixel data.  The
C  record length in the direct access read is set equal to the SIZE
C**********************************************************************
c
      nf=index(filename,' ')-1
      OPEN(UNIT=8, FILE=FILENAME, ERR=100, IOSTAT=IOSTATUS,
     & ACCESS='DIRECT', FORM='UNFORMATTED',RECL= SIZE)

      if(l_cell_afwa)then
         imagelen=(nlfi*nefi)/4
         READ (8,REC=1,ERR=99)HEADER,IMAGEC
      else
         imagelen=nlfi*nefi
         call read_binary_field(image_decell_1d,1,4,imagelen+
     +HEADER_SIZE*4,filename,nf)
      endif

      Close(8)

      if(l_cell_afwa)then

         call decellularize_image(IMAGEC,imagelen,width,depth,
     &pixels_per_cell,cell_width,cell_depth,nefi,nlfi,image_decell_2d)

c back as floating 10-bit info for the sector in domain.

         jj=0
         do j=jstart,jend
            jj=jj+1
            ii=0
         do i=istart,iend
            ii=ii+1
            image_out(ii,jj)=float(image_decell_2d(i,j))*4.0
         enddo
         enddo

      else ! afwa data not cellularized but bits need moving

c GOES data = 10 bit; METEOSAT data = 8 bit.
         r8to10=4.0
         if(isat.eq.2 .and. jtype.eq.4)r8to10=1.0

         m=HEADER_SIZE*4
         do j=1,nlfi
         do i=1,nefi
            m=m+1
            image_decell_2d(i,j)=image_decell_1d(m)
         enddo
         enddo

         if(isat.eq.2 .and. jtype.eq.4)then

c METEOSAT origin is SE corner (wrt ri/rj lut). Make it SW corner
            do j=1,nlfi
            do i=1,nefi
               image_temp(nefi-i+1,j)=image_decell_2d(i,j)
            enddo
            enddo
            do j=1,nlfi
            do i=1,nefi
               image_decell_2d(i,j)=image_temp(i,j)
            enddo
            enddo

         endif
c
c now load that part within the LAPS domain; convert to 10-bit.
c
         jj=0
         do j=jstart,jend
            jj=jj+1
            ii=0
         do i=istart,iend
            ii=ii+1
            itempintgr=IBITS(image_decell_2d(i,j),0,8)
            image_out(ii,jj)=float(itempintgr)*r8to10
         enddo
         enddo

      endif

      istatus = 1
      goto 1000

99    Write(6,*)'Error Reading GWC file'
      Write(6,*)'Filename = ',filename(1:nf)
      Close(8)
      return 

100   IF (IOSTATUS .NE. 0)THEN
        PRINT *, 'ERROR READING ',FILENAME, ' IO status is', 
     &  IOSTATUS
        istatus = -1
      END IF

1000  RETURN
      END
