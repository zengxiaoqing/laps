c*****************************************************************
c
c The .swp disc file is read by first reading the header and then
c reading and unpacking the dbz data one row at a time.
c Subroutine READHEAD may be used to open and read the radar
c header.  Subroutines PASSXCMP and UNPAKDBZ2D may be used to read
c and unpack the data respectively.  The READHEAD, PASSXCMP and
c UNPAKDBZ2D subroutines are given below(variables beginning w/
c ijklmn are INTEGER*4 and the rest are REAL unless explicitly
c declared).
c
c The size of the maps in bins is denoted by the variables names
c "imax" and "jmax".  The spacing of the data is denoted by the
c variable names "sx" and "sy".
c
c *****************************************************************

      SUBROUTINE READHEAD(luout, name, lunit, ierr, keyword, fltname,
     +  stmname, radarid, trackname, createtime, ramname, imax, jmax,
     +  kmax, nsweeps, nmosmflag, iunfoldflag, intattflag, ieditflag,
     +  iextra2, iextra3, stime, etime, olat, olon, sx, sy, sz, xdis,
     +  ydis, z0, rot, radaralt, calibco1, calibco2, azmcor, elcor,
     +  thresh, dbzrs, pcor, dcor, rcor, starthorelev, htbb, dbb)
C     OPENS RADAR DATA FILE AND READS RADAR HEADER DATA VALUES.
C     PAUL A. LEIGHTON, USDC/NOAA/AOML/HRD, 4 JUN 1991
      CHARACTER      keyword*4, fltname*8, stmname*12, radarid*4,
     +               trackname*32, createtime*32, ramname*28,
     +               name*140, jfile*140

C     CLOSE PREVIOUSLY OPENED FILES
      CLOSE(lunit)

C     OPEN THE FILE AND READ THE HEADER VARIABLES.
      OPEN(lunit, ERR = 991, FILE = name, IOSTAT = ierr,
     +  STATUS = 'old', FORM = 'unformatted') 
   10 READ(lunit, IOSTAT = ierr, ERR = 993) keyword, fltname, stmname,
     +  radarid, trackname, createtime, ramname, imax, jmax, kmax,
     +  nsweeps, nmosmflag, iunfoldflag, intattflag, ieditflag,
     +  iextra2, iextra3, stime, etime, olat, olon, sx, sy, sz,
     +  xdis, ydis, z0, rot, radaralt, calibco1, calibco2,
     +  azmcor, elcor, thresh, dbzrs, pcor, dcor, rcor,
     +  starthorelev, htbb, dbb
   
      RETURN

C     ERROR MESSAGES: 
  991 WRITE(luout, '(''Error '', I4, '' on file '', A40,
     +  ''trying /hrd/dat/'')') ierr, name 
      jfile = '/hrd/dat/' // name 
      OPEN(lunit, ERR = 995, FILE = jfile, IOSTAT = ierr,
     +  STATUS = 'old', FORM = 'unformatted') 
      GOTO 10

  993 WRITE(luout, '(//, "***********************", //,
     + "READ ERROR ON HEADER", //, "**********************")')

  995 WRITE(luout, '(//, "**********************", //, "Error ", I4,
     +  " on file ", A40, "not on /hrd/dat/", //,
     +  "**********************")') ierr, jfile

      RETURN
      END

c *****************************************************************

      SUBROUTINE PASSXCMP(lunit, imax, jmax, max_x, max_y, zarray)
C     READS AND UNPACKS COMPOSITE DATA VALUES STORED TWO PER WORD IN KPAC AND
C     RETURNS THEM IN ZARRAY FOR LATER??? USE.  ASSUMES HEADER IS PREVIOUSLY
C     READ.  WRITTEN BY PAUL A. LEIGHTON, USDC/NOAA/AOML/HRD, 12 AUG 92.
C     LUOUT - PRINTED OUTPUT LU NUMBER FOR ERROR MESSAGE.
C     LUNIT - DISC LU NUMBER OF COMPOSITE FILE OPENED PREVIOUSLY IN READHEAD.
C     IMAX - ROW DIMENSION OF Z ARRAY,
C          - MAXIMUM NUMBER OF Z VALUES PER ROW IN KPAC
C     JMAX - COLUMN DIMENSION OF Z ARRAY,
C          - COLUMN DIMENSION OF STORED COMPOSITE FILE BEING READ INTO KPAC.
C     KPAC - PACKED ARRAY OF DBZ COMPOSITE VALUES 
C     ZARRAY - BUFFER OF UNPACKED DBZ VALUES TO BE RETURNED
      INTEGER*2    kpac(32767)   ! DBZ COMPOSITE VALUES STORED 2 PER WORD
      INTEGER*2    zarray(0:max_x-1,0:max_y-1) !ROW OF DBZ VALUES TO BE 
                                               !RETURNED
      data luout/6/
c check array boundaries:
      max=(imax*jmax)/2
      if (max .gt. 32767) then
         write(luout,'(" The boundaries for imax or jmax are to big.")')
         stop
      endif
c read packed dbz values from disc:
      do jx=1, jmax
         j=jmax+1-jx
         n1=((j-1)*imax+2)/2
         n2=(imax+(j-1)*imax+1)/2
         read(lunit,iostat=ierr,err=991) (kpac(i),i=n1,n2)
      end do
c unpack data:
      do j = 0, jmax-1
         DO i = 0, imax-1
            CALL UNPAKDBZ2D(i+1, j+1, kpac, imax, jmax, z) ! UNPACKS KPAC.
            zarray(i,j) = nint(z)  ! STORE DBZ VALUE IN OUTPUT ARRAY.
         ENDDO
      end do
      GOTO 9999  ! END SUBROUTINE

C     ERROR MESSAGES: 
  991 WRITE(luout, '(''Error on read:'', I4)') ierr
      write(luout, '("Aborting program!!!!")')
      stop

 9999 RETURN  ! TO CALLING SUBROUTINE
      END

c *****************************************************************

      SUBROUTINE UNPAKDBZ2D(i, j, kpac, imax, jmax, z)

      INTEGER*2 kpac(32767)
      integer   i4_kpac

      n = (i + (j-1)*imax + 1) / 2 
      i4_kpac = kpac(n)
      IF (MOD(i, 2).EQ.0)  THEN
          l = IAND(i4_kpac, 255)
      ELSE
          l = ISHFT(i4_kpac, -8)
      ENDIF
      z = (l - 64) / 2.
      IF (l.EQ.0)  z = -999.0 
      RETURN
      END
      
c *****************************************************************
