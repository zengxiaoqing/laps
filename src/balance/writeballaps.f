        subroutine write_bal_laps(i4time,phi,u,v,temp,om
     .  ,imax,jmax,kmax,ip,istatus)
C
       IMPLICIT       NONE
C
       INTEGER      I4TIME,
     1              IMAX,JMAX,KMAX,
     1              KMAX3,kmax2,
     1              IP(KMAX),
     1              LVL(kmax*3),
     1              I,J,K,ERROR(2),
     1              ISTATUS
C
       integer      lend

       REAL*4       PHI(IMAX,JMAX,KMAX),
     1              U(IMAX,JMAX,KMAX),
     1              V(IMAX,JMAX,KMAX),
     1              OM(IMAX,JMAX,KMAX),
     1              temp(imax,jmax,kmax),
     1              bal(imax,jmax,kmax*3)
C
       character*150   dir
       character*150   directory
       CHARACTER*31    EXT
       CHARACTER*3     VAR(kmax*3)
       CHARACTER*4     LVL_COORD(kmax*3)
       CHARACTER*9     fname9
       CHARACTER*10    UNITS(kmax*3)
       CHARACTER*125   COMMENT(kmax*3)
C
C-------------------------------------------------------------------------------
C
c first is the wind (lw3)
c

        ext='lw3'
        call get_directory(ext,directory,lend)
        lend=lend-4  ! get_directory insures that a "/" is at the end of the directory string.

        dir=directory(1:lend)//'balance/'//ext(1:3)//'/'
C
        DO K=1,KMAX
                VAR(K)='U3 '
                LVL(K)=IP(K)
                LVL_COORD(K)='MB  '
                UNITS(K)='M/S   '
                COMMENT(K)='Non-linear balanced u-component wind.'
                DO J=1,JMAX
                DO I=1,IMAX
                        BAL(I,J,K)=U(I,J,K)
                ENDDO
                ENDDO
        ENDDO
C
        DO K=KMAX+1,2*KMAX
                VAR(K)='V3 '
                LVL(K)=IP(K-KMAX)
                LVL_COORD(K)='MB  '
                UNITS(K)='M/S   '
                COMMENT(K)='Non-linear balanced v-component wind.'
                DO J=1,JMAX
                DO I=1,IMAX
                        BAL(I,J,K)=V(I,J,K-KMAX)
                ENDDO
                ENDDO
        ENDDO
C
        DO K=2*KMAX+1,3*KMAX
                VAR(K)='OM '
                LVL(K)=IP(K-2*KMAX)
                LVL_COORD(K)='MB  '
                UNITS(K)='PA/S  '
                COMMENT(K)='Non-linear balanced omega.           '
                DO J=1,JMAX
                DO I=1,IMAX
                        BAL(I,J,K)=OM(I,J,K-2*KMAX)
                ENDDO
                ENDDO
        ENDDO

        call make_fnam_lp(i4time,fname9,istatus)
        write(6,*)' Writing grids ',ext(1:3),' ',fname9

        kmax3=kmax*3
        call write_laps_data(i4time,dir,ext,imax,jmax
     +,kmax3,kmax3,var,lvl,lvl_coord,units,comment,bal,istatus)
c
c----------------------
c now lt1

       EXT='lt1'
       dir=directory(1:lend)//'balance/'//ext(1:3)//'/'
       DO K=1,KMAX
              VAR(K)='HT '
              LVL(K)=IP(K)
              LVL_COORD(K)='MB  '
              UNITS(K)='Meters'
              COMMENT(K)='Non-linear balanced height.'
              DO J=1,JMAX
              DO I=1,IMAX
                     BAL(I,J,K)=PHI(I,J,K)
              ENDDO
              ENDDO
       ENDDO
C
       DO K=KMAX+1,2*KMAX
              VAR(K)='T3'
              LVL(K)=IP(K-KMAX)
              LVL_COORD(K)='MB  '
              UNITS(K)='Kelvin'
              COMMENT(K)='Non-linear balanced temp.'
              DO J=1,JMAX
              DO I=1,IMAX
                     BAL(I,J,K)=TEMP(I,J,K-KMAX)
              ENDDO
              ENDDO
       ENDDO
C
       KMAX2=2*KMAX
C
       write(6,*)' Writing grids ',ext(1:3),' ',fname9

       CALL WRITE_LAPS_DATA(I4TIME,DIR,EXT,IMAX,JMAX,KMAX2,KMAX2,VAR,LVL
     1                  ,LVL_COORD,UNITS,COMMENT,BAL,ISTATUS)
C
       IF (ISTATUS.ne.1) THEN
              PRINT*,'Error writing balanced data.'
              ISTATUS=0
              RETURN
       ENDIF
C
       ISTATUS=1
       RETURN
C
       END
