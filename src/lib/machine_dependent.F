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
        function i4_to_byte(i4_in)

        character i4_to_byte

        integer*4 i4_in,i4

        character barg(4)
        equivalence (barg,i4)

        i4 = i4_in
        i4_to_byte = barg(4)

        return

        end

        function byte_to_i4(b_in)

        integer*4 i4,byte_to_i4
        character b_in

        character barg(4)
        equivalence (barg,i4)

        i4 = 0

        barg(4) = b_in
#ifndef hpux
        if(i4 .ge. 128)i4 = i4 - 256
#endif

        byte_to_i4 = i4

!       write(6,*)i4,byte_to_i4

        return

        end
#ifdef USE_TRIGD
      module trigd
         public
         private pi_180
         interface sind
            module procedure rsind, dsind
         end interface
         interface cosd
            module procedure rcosd, dcosd
         end interface
         interface tand
            module procedure rtand, dtand
         end interface
         interface asind
            module procedure rasind, dasind
         end interface
         interface acosd
            module procedure racosd, dacosd
         end interface
         interface atand
            module procedure ratand, datand
         end interface

         interface atan2d
            module procedure ratan2d, datan2d
         end interface
         contains
            double precision function pi_180()
              pi_180 = 2.d0*acos(0.d0)/180.d0
            return 
            end function pi_180

            real function rsind(val)
            real, intent(in) :: val
              rsind = sin(val*pi_180())
            return
            end function rsind

            real function rasind(val)
            real, intent(in) :: val
              rasind = asin(val)/pi_180()
            return
            end function rasind

            real function rcosd(val)
            real, intent(in) :: val
              rcosd = cos(val*pi_180())
            return
            end function rcosd

            real function racosd(val)
            real, intent(in) :: val
              racosd = acos(val)/pi_180()
            return
            end function racosd


            real function rtand(val)
            real, intent(in) :: val
              rtand = tan(val*pi_180())
            return
            end function rtand

            real function ratand(val)
            real, intent(in) :: val
              ratand = atan(val)/pi_180()
            return
            end function ratand

            real function ratan2d(val1,val2)
            real, intent(in) ::  val1,val2
              ratan2d = atan2(val1,val2)/pi_180()
            return
            end function ratan2d

            double precision function dsind(val)
            double precision, intent(in) :: val
              dsind = sin(val*pi_180())
            end function dsind

            double precision function dasind(val)
            double precision, intent(in) :: val
              dasind = asin(val)/pi_180()
            return
            end function dasind

            double precision function dcosd(val)
            double precision, intent(in) :: val
              dcosd = cos(val*pi_180())
            end function dcosd

            double precision function dacosd(val)
            double precision, intent(in) :: val
              dacosd = acos(val)/pi_180()
            return
            end function dacosd

            double precision function dtand(val)
            double precision, intent(in) :: val
              dtand = tan(val*pi_180())
            end function dtand
            
            double precision function datand(val)
            double precision, intent(in) :: val
              datand = atan(val)/pi_180()
            return
            end function datand


            double precision function datan2d(val1,val2)
            double precision, intent(in) ::  val1,val2
              datan2d = atan2(val1,val2)/pi_180()
            return
            end function datan2d
      end module

#endif
      subroutine open_append(lun,file,status,istat)
      integer istat
      character*(*) file,status

      istat=1

#if defined(hpux) && !defined(F90)
      open(lun,file=file,status=status,ACCESS='APPEND',err=998)
#elif defined(IRIX) && !defined(F90)
      open(lun,file=file,status=status,ACCESS='APPEND',err=998)
#else
      open(lun,file=file,status=status,POSITION='APPEND',err=998)
#endif

      return
 998  istat=0
      return
      end

#if defined(hpux) && !defined(F90)
      real function transfer(in,out)
      integer in
      real out
      print*, 'TRANSFER is not supported by this compiler'
      print*, 'Please report this error to laps_bugs@fsl.noaa.gov'
      stop
      transfer=0
      return
      end
#endif

C Below are some Taiwan FGGE routines

#if defined(hpux) || defined(alpha)
C---------------------------------------------------------------------------
C PURPOSE :
C   DEMOSTRATION PROGRAM TO USE FGGE DECODE LIBRARY ON ULTRIX 
C BUILD :
C   f77 -o unpk_alpha  unpk_alphx.f unpkfgge_a.o
C USAGE :
C   $ unpk_alpha 
C NOTE :
C   If running on DEC-alpha change the INTEGER*4 to INTEGER*8 on 
C   the definition of FGGE_HEAD structure
C HISTORY
C   MAR-03-1994   C.P.DZEN   ORIGIONAL
C---------------------------------------------------------------------------
C-
C---------------------------------------------------------------------------
C-
      SUBROUTINE TEST_2
      PARAMETER ( MAX_BUF=320000,MAX_PTS=160000 )
      CHARACTER FGGE_TEMP*80,FGGE_NAME*40,FGGE_IPT*30
      CHARACTER*1 FGGE_BUF(MAX_BUF)      
      REAL*4    FDATA(MAX_PTS)
      CHARACTER*1 ANS
C
      STRUCTURE /FGGE_HEAD/
#ifdef alpha
	INTEGER*8 Q,S1,F1,T,C1,E1,M,X,S2,F2,N,C2,E2,CD,CM,KS,K,U1
        INTEGER*8 NW,JJ,MM,YY,GG,R,G,J,B,Z,A,U2,SN,MN,RT,DO,U3
#else
	INTEGER*4 Q,S1,F1,T,C1,E1,M,X,S2,F2,N,C2,E2,CD,CM,KS,K,U1
        INTEGER*4 NW,JJ,MM,YY,GG,R,G,J,B,Z,A,U2,SN,MN,RT,DO,U3
#endif
      END STRUCTURE
      RECORD /FGGE_HEAD/HEADER
C
c      FGGE_NAME = 'gf94031800m.;1'//CHAR(0)
c      FGGE_NAME = '$1$DUA7:[NWPIIDAT]E4032012.;1'//CHAR(0)
c      CALL GETARG(1,FGGE_TEMP)
c      PRINT *,'FGGE_NAME=',FGGE_TEMP
c     FGGE_NAME = 'E4032012'//CHAR(0)
       write(*,'(1x,a$)') 'Enter FGGE file name : '
       read(*,'(a)') FGGE_IPT
c      FGGE_NAME = FGGE_IPT//CHAR(0)
c      write(*,'(''>'',a,''<'')') FGGE_IPT 
c      write(*,'(''>'',a,''<'')') FGGE_NAME
C
      CALL OPEN_FGGE_FILE(FGGE_IPT,ISTAT)
      IF( ISTAT .NE. 0 ) STOP 'open_fgge_file error !!!'
C
      DO 10 I=1,100
C
      CALL READ_FGGE_RECORD(FGGE_BUF,ISTAT)
      IF( ISTAT .NE. 0 ) STOP 'read_fgge_record error !!!'
C
      CALL DECODE_FGGE_HEADER(FGGE_BUF,HEADER,ISTAT)
      IF( ISTAT .NE. 0 ) STOP 'decode_fgge_header error !!!'
      CALL PRT_FGGE_HEADER(HEADER)
C
      CALL DECODE_FGGE_DATA(FGGE_BUF,HEADER,FDATA,ISTAT)
      IF( ISTAT .NE. 0 ) STOP 'decode_fgge_data error !!!'
      CALL PRT_FGGE_DATA(HEADER,FDATA,HEADER.J)
C
      WRITE(*,*) ' '
      WRITE(*,*) ' Read next ?'
      READ(*,'(A1)') ANS
      IF( ANS.NE.'Y' .AND. ANS.NE.'y' ) GOTO 20
   10 CONTINUE
C
   20 CONTINUE
      CALL CLOSE_FGGE_FILE(ISTAT)
C
      RETURN
      END
C---------------------------------------------------------------------------
      SUBROUTINE PRT_FGGE_HEADER(HEADER)
      STRUCTURE /FGGE_HEAD/
#ifdef alpha
        INTEGER*8 Q,S1,F1,T,C1,E1,M,X,S2,F2,N,C2,E2,CD,CM,KS,K,U1
        INTEGER*8 NW,JJ,MM,YY,GG,R,G,J,B,Z,A,U2,SN,MN,RT,DO,U3
#else
	INTEGER*4 Q,S1,F1,T,C1,E1,M,X,S2,F2,N,C2,E2,CD,CM,KS,K,U1
        INTEGER*4 NW,JJ,MM,YY,GG,R,G,J,B,Z,A,U2,SN,MN,RT,DO,U3
#endif
      END STRUCTURE
      RECORD /FGGE_HEAD/HEADER

      PRINT *,'HEADER.Q  : ',HEADER.Q		! DATA TYPE
      PRINT *,'HEADER.S1 : ',HEADER.S1		! TYPE OF SURFACE 1
      PRINT *,'HEADER.F1 : ',HEADER.F1		! FORECAST TAU
      PRINT *,'HEADER.T  : ',HEADER.T		! SET TO 0 FOR FORECAST FIELD
      PRINT *,'HEADER.C1 : ',HEADER.C1		! Pressure or Hight value
      PRINT *,'HEADER.E1 : ',HEADER.E1		! exponential position for C1
      PRINT *,'HEADER.M  : ',HEADER.M
      PRINT *,'HEADER.X  : ',HEADER.X
      PRINT *,'HEADER.S2 : ',HEADER.S2
      PRINT *,'HEADER.F2 : ',HEADER.F2
      PRINT *,'HEADER.N  : ',HEADER.N
      PRINT *,'HEADER.C2 : ',HEADER.C2
      PRINT *,'HEADER.E2 : ',HEADER.E2
      PRINT *,'HEADER.CD : ',HEADER.CD
      PRINT *,'HEADER.CM : ',HEADER.CM
      PRINT *,'HEADER.KS : ',HEADER.KS
      PRINT *,'HEADER.K  : ',HEADER.K
      PRINT *,'HEADER.U1 : ',HEADER.U1
      PRINT *,'HEADER.NW : ',HEADER.NW		!
      PRINT *,'HEADER.JJ : ',HEADER.JJ		! YEAR
      PRINT *,'HEADER.MM : ',HEADER.MM		! MONTH
      PRINT *,'HEADER.YY : ',HEADER.YY		! DAY
      PRINT *,'HEADER.GG : ',HEADER.GG		! HOUR
      PRINT *,'HEADER.R  : ',HEADER.R
      PRINT *,'HEADER.G  : ',HEADER.G
      PRINT *,'HEADER.J  : ',HEADER.J		! TOTAL GRID POINTS
      PRINT *,'HEADER.B  : ',HEADER.B		!
      PRINT *,'HEADER.Z  : ',HEADER.Z
      PRINT *,'HEADER.A  : ',HEADER.A		!
      PRINT *,'HEADER.U2 : ',HEADER.U2
      PRINT *,'HEADER.SN : ',HEADER.SN		!
      PRINT *,'HEADER.MN : ',HEADER.MN		! FLAP ID
      PRINT *,'HEADER.RT : ',HEADER.RT		! RECORD TYPE
      PRINT *,'HEADER.DO : ',HEADER.DO		! DOMAIN ID
      PRINT *,'HEADER.U3 : ',HEADER.U3		! Coding method
      RETURN
      END
C---------------------------------------------------------------------------
      SUBROUTINE PRT_FGGE_DATA(HEADER,FDATA,N)
      REAL*4 FDATA(N)
      STRUCTURE /FGGE_HEAD/
#ifdef alpha
        INTEGER*8 Q,S1,F1,T,C1,E1,M,X,S2,F2,N,C2,E2,CD,CM,KS,K,U1
        INTEGER*8 NW,JJ,MM,YY,GG,R,G,J,B,Z,A,U2,SN,MN,RT,DO,U3
#else
	INTEGER*4 Q,S1,F1,T,C1,E1,M,X,S2,F2,N,C2,E2,CD,CM,KS,K,U1
        INTEGER*4 NW,JJ,MM,YY,GG,R,G,J,B,Z,A,U2,SN,MN,RT,DO,U3
#endif
      END STRUCTURE
      RECORD /FGGE_HEAD/HEADER

      DO 10 I=1,100,5
      WRITE(*,'(4X,5F11.4)')
     *  FDATA(I),FDATA(I+1),FDATA(I+2),FDATA(I+3),FDATA(I+4)
   10 CONTINUE
      WRITE(*,'(/)')
      DO 20 I=HEADER.J-99,HEADER.J,5
      WRITE(*,'(4X,5F11.4)')
     *  FDATA(I),FDATA(I+1),FDATA(I+2),FDATA(I+3),FDATA(I+4)
   20 CONTINUE

      RETURN
      END

#else
C
C     Dummy routines for other platforms

      SUBROUTINE TEST_2
      x = 1.0
      write(6,*)' WARNING: routine TEST_2 not supported '
     1         ,'on this platform'
      RETURN
      END

      SUBROUTINE PRT_FGGE_HEADER(HEADER)
      x = 1.0
      write(6,*)' WARNING: routine PRT_FGGE_HEADER not supported '
     1         ,'on this platform'
      RETURN
      END

      SUBROUTINE PRT_FGGE_DATA(HEADER,FDATA,N)
      x = 1.0
      write(6,*)' WARNING: routine PRT_FGGE_DATA not supported '
     1         ,'on this platform'
      RETURN
      END
#endif
