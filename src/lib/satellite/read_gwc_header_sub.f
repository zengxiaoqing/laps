      SUBROUTINE read_gwc_header(FILENAME,STRPIX,STRSCNL,STPPIX,STPSCNL,
     &REQOBSTM,IMGTYPE,GOLATSBP,GOLONSBP,iwidth,idepth,GOALPHA,istrbdy1,
     &istrbdy2,istpbdy1,istpbdy2,BEPIXFC,BESCNFC,FSCI,idecimat,IOSTATUS)
C
C***********************************************************************
C  Purpose:  Read and decode the header information for GVAR data coming
C  from SDHS.
C
C  Method:  Read in character data, short and long integer data, floating
C  point data and double precesion data.  All but the character data is 
C  in a vax specific format.  The bits and bytes of this data must be 
C  manipulated to obtain the correct representation.  The method used is
C  platform independent.
C  
C  References:
C  1.  Final Interface Specification for the Satellite Data Handling 
C  System Communications Network (SDHS-COMNET), 01 February 1988
C  2.  AFGWC/DONS PV-Wave program auto_convert.pro 
C  3.  VAX Fortran Reference Manual from the SDHS programmer library  
C***********************************************************************

C***********************************************************************
C  Begin variable declaration section.  The variables are declared in
C  functional groups. The header variables are declared in two groups in 
C  the order that they are stored. Temporary/working variables are then
C  declared.
C***********************************************************************

      IMPLICIT NONE    
      CHARACTER*(*) FILENAME
      CHARACTER*3 DATATYP  !Data Type
      CHARACTER*3 DATASBTP !Data Subtype
cc      character*2 SPARE1     !Spare
      INTEGER   HDRSIZE    !Header Size
      INTEGER DATASTR    !Data Start
      INTEGER SHPTPTR    !Shipping Table Ptr
      Integer lun

      INTEGER STRPIX     !Start Pixel
      INTEGER STRSCNL    !Start Scanline 
      INTEGER STPPIX     !Stop Pixel
      INTEGER STPSCNL    !Stop Scanline

      integer n

cc      CHARACTER*6 SATID    !Satellite Identification
cc      CHARACTER*2 ORGDATTP !Data Type
      INTEGER REQOBSTM   !Requested Observation Time

cc      INTEGER DECIMAT    !Decimation Resolution 
      Integer   idecimat   !Returned variable

cc      CHARACTER*3 ORGSNRES !Original Sensor Resolution

cc      INTEGER STRBDY1    !Requested Start Boundary 1
      Integer   istrbdy1   !Returned variable

cc      INTEGER STRBDY2    !Requested Start Boundary 2
      Integer   istrbdy2   !Returned variable

cc      INTEGER STPBDY1    !Requested Stop Boundary 1 
      Integer   istpbdy1   !Returned Variable

cc      INTEGER STPBDY2    !Requested Stop Boundary 2
      Integer   istpbdy2   !Returned Variable

cc      INTEGER GOVISMD    !GOES Vis Mode
      REAL*4 GOLONSBP      !GOES Longitude Subpoint
      REAL*4 GOLATSBP      !GOES Latitude Subpoint
      REAL*4 GOALPHA       !GOES Alpha
cc      REAL*4 GODELTA       !GOES Delta
cc      REAL*4 GOZETA        !GOES Zeta
cc      REAL*4 GOETA         !GOES Eta
cc      REAL*4 GORHO         !GOES Rho
cc      REAL*4 GOSBSCN       !Scanline of Satellite subpoint
cc      REAL*4 GOSBSAMP      !Pixel of Satellite subpoint
cc      REAL*4 GOGEOCAL      !GOES Geoc altitude????     
cc      INTEGER SPARE2(6)  !Spare Bytes not currently used
      INTEGER BEPIXFC    !Begin Pixel First Cell
      INTEGER BESCNFC    !Begin Scaneline First Cell

cc      INTEGER WIDTH      !Tracks in Width of image
      Integer   iwidth     !Returned variable

cc      INTEGER DEPTH      !Tracks in Depth of image
      Integer   idepth     !Returned Variable

cc      INTEGER LSTPIX     !Last Pixel in complete image
cc      INTEGER LSTSCN     !Last Scanline in complete imagec
cc      BYTE REQDC           !Required data complete
cc      CHARACTER*3 ORTYPE   !Original Type
cc      CHARACTER*3 ORSTYPE  !Original Subtype
cc      REAL*8 POLARIN       !Polar Inclination
      CHARACTER*2 IMGTYPE  !Image Type
cc      BYTE ORPROJ          !Original Reference Projection
cc      CHARACTER*1 ORFIELD  !Orientation Field
cc      REAL*8 GSUBLON       !GOES Subpoint Longitude Spare for DMSP
cc      CHARACTER*1 ORSFLAG  !Original request Size Flag
cc      REAL*8 ASCLON        !Ascending Longitude     
cc      REAL*8 ANOMMM        !Anomolous Mean Motion
cc      REAL*8 RELERR        !Relative Earth Rotation Rate
cc      REAL*8 ARGPER        !ARG Perigee
cc      REAL*8 ECCEN         !Eccentricity
cc      REAL*8 ARGLAT        !ARG Latitude
cc      REAL*8 DTTOP         !DT Top
cc      REAL*8 DTBOT         !DT Bottom
cc      REAL*8 MIDTIME       !Midpoint Time
      INTEGER FSCI       !First Scanline of Complete Image
cc      INTEGER STCI       !Start Time of Complete Image
cc      BYTE CHAN            !Channel 1=IR1, 2=IR2, 3=IR3, 4=Water Vap
cc      REAL*4 ROTANG        !Rotation Angle between Aries and Greenwich
cc      REAL*4 ABSMAG        !Absolute Magnitude of craft in inertial 
                           !reference frame
cc      REAL*4 SATERTH(3,3)  !Satellite to Earth Tranformation Matrix
cc      REAL*4 POINTVC       !Pointing Vector Component
cc      REAL*4 SATVEC(3)     !Satellite Position Vector
cc      REAL*4 RMA           !Role Misalignment Angle
cc      REAL*4 PMA           !Pitch Misalignment Angle


      integer byteswp2, byteswp4
      character input(744) !just read it all and sort it later!

c      BYTE JUNK512(512)     !Byte array used to 'read over' unneeded
c                            !data
c      BYTE JUNK95(95)       !Byte array used to 'read over' unneeded
c                            !data
      INTEGER IOSTATUS    !I/O status flag
      INTEGER I        !Array indices
      logical lopen
C***********************************************************************
C  Open the file that contains the GVAR header and pixel data.  The 
C  record length in the direct access read is 1024 bytes (the number of
C  bytes that make up the header).
C***********************************************************************

      IOSTATUS = 1

      write(6,*)'Opening ',FILENAME
      inquire(file=filename,opened=lopen,number=lun)
      if(lopen)write(6,*)'File already open',lun
      OPEN(UNIT=8, FILE=FILENAME, ERR=100, IOSTAT=IOSTATUS,
     &ACCESS='DIRECT', RECL=744,STATUS='OLD')

cc     &ACCESS='DIRECT', FORM='UNFORMATTED',RECL=1024,STATUS='OLD')

C***********************************************************************
C  Read in the first 20 bytes of header information into the variables 
C  listed.  The rest of the first 512 bytes do not contain useful 
C  information for GVAR data.
C***********************************************************************
      READ (8,REC=1) input

c                      3       3        2       4       4      4
cc    READ (8,REC=1) DATATYP,DATASBTP,SPARE1,HDRSIZE,DATASTR,SHPTPTR

C***********************************************************************
C  Because the multi-byte integer data is coming from a VAX machine the 
C  bytes are not in the correct order for a unix machine.  The routine 
C  BYTESWP4 does the funcitonal equivalent of a byte swap.  The resulting
C  integer representation of the bit string is true for positive integers  
C  but may not be true for negative integer (this depends on whether 
C  the negative integers were stored in complement form).  
C***********************************************************************
      do i=1,3
         datatyp(i:i) = input(i)
         datasbtp(i+3:i+3) = input(i+3)
      enddo
      hdrsize = byteswp4(input(9))
      datastr = byteswp4(input(13))
      shptptr = byteswp4(input(17))      

        
      PRINT*, DATATYP, ' ',DATASBTP, ' ',HDRSIZE, ' ',
     &   DATASTR, ' ', SHPTPTR

C***********************************************************************
C  Read in the bulk of the header data.  Read over the first 512 bytes,
C  then read in the values for the header variables.  These values are
C  in VAX binary format for character, short (2 byte) integer, long (
C  4 byte) integer, floating point, and double precision.  To obtain the 
C  proper resprentation on a unix machine the bits and bytes of this infor-
C  mation must be manipulated.  Bytes must be swapped for short and long
C  integers.  The sign bit, exponent bits and fraction bits in double and
C  floating point values must be properly interpreted.  All floating 
C  point and double precisions data are first read into working integer
C  variables and then the bit strings for the different components are 
C  "picked" out and the floating point number reconstructed.  Two 
C  integer variables must be used to read in the bit strings for double 
C  precision data  
C***********************************************************************
c                     1        513    517      521     525     529
c                     512        4      4       4        4      6
c      READ (8,REC=1) JUNK512, STRPIX, STRSCNL, STPPIX,STPSCNL,SATID,
c       535      537       541       543       546      548
c       c*2       4          2        3         2       2
c     & ORGDATTP, REQOBSTM, DECIMAT, ORGSNRES, STRBDY1, STRBDY2, 
c      550      552      554       556       560       564
c       2         2       2         4          4        4
c     & STPBDY1,STPBDY2, GOVISMD, GOLONSBI, GOLATSBI, GOALPHAI,
c       568         572     576     580      584      588
c         4           4     4        4         4       4
c     & GODELTAI, GOZETAI, GOETAI, GORHOI, GOSBSCNI, GOSBSAMPI, 
c        592        596      620     624      628     630    632
c          4        24          4      4       2      2       4
c     & GOGEOCALI, SPARE2, BEPIXFC, BESCNFC, WIDTH, DEPTH, LSTPIX,
c        636     640   641     644      647         651        655
c          4      1      3       3          4        4           2
c     & LSTSCN, REQDC, ORTYPE, ORSTYPE, POLARINI1, POLARINI2, IMGTYPE,
c       657     658       659         663        667     668
c        1       1         4           4         1        4
c     & ORPROJ, ORFIELD, GSUBLONI1, GSUBLONI2, ORSFLAG, ASCLONI1,
c         672      676       680       684       688
c           4       4         4         4         4
c     & ASCLONI2, ANOMMMI1, ANOMMMI2, RELERRI1, RELERRI2,   
c        692      696        700     704       708       712
c         4         4         4       4         4         4
c     & ARGPERI1,ARGPERI2, ECCENI1, ECCENI2, ARGLATI1, ARGLATI2,   
c         716      720    724       728    732     736       740
c         4         4       4        4      4       4          4
c     & DTTOPI1,DTTOPI2,DTBOTI1,DTBOTI2,MIDTIMEI1,MIDTIMEI2, FSCI,
c        745   749   750   845     849     851      887     891
c         4     1    95     4        4      36        4       12
c     & STCI, CHAN,JUNK95,ROTANGI,ABSMAGI,SATERTHI,POINTVCI,SATVECI,
c       903    907
c        4      4
c     & RMAI, PMAI
      CLOSE(8)
C***********************************************************************
C  Perform a 4 byte swap for long integers and a two byte swap for short
C  integers
C***********************************************************************

      strpix  = BYTESWP4(input(513))
      strscnl = BYTESWP4(input(517))
      stppix  = BYTESWP4(input(521))
      stpscnl = BYTESWP4(input(525))
      reqobstm= BYTESWP4(input(537))

      idecimat = BYTESWP2(input(541))
      istrbdy1  = BYTESWP2(input(546))
      istrbdy2  =  BYTESWP2(input(548))
      istpbdy1 =  BYTESWP2(input(550))
      istpbdy2 =  BYTESWP2(input(552))

      CALL CNVTVXFL(input(556),GOLONSBP)
      CALL CNVTVXFL(input(560),GOLATSBP)
      CALL CNVTVXFL(input(564),GOALPHA)

      bepixfc =  BYTESWP4(input(620))
      bescnfc =  BYTESWP4(input(624))
      iwidth   =  BYTESWP2(input(628))
      idepth   =  BYTESWP2(input(630))

      do i=1,2
         imgtype(i:i)=input(654+i)
      enddo

      FSCI   =  BYTESWP4(input(740))

c
c   No further header variables are required - return here.
c


C***********************************************************************
C  If there is an error opening the data file print out a message and
C  end the program
C***********************************************************************

 100  IF (IOSTATUS .NE. 0)THEN
        PRINT *, 'ERROR READING FILE ', FILENAME , 'IO status is', 
     &  IOSTATUS
      END IF


      RETURN 
      END
