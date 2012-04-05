        subroutine write_bal_laps(i4time,phi,u,v,temp,om,rh
     .  ,sh,imax,jmax,kmax,p,istatus)
C
       IMPLICIT       NONE
C
       INTEGER      I4TIME,
     1              IMAX,JMAX,KMAX,
     1              KMAX3,kmax2,
     1              IP(KMAX),
     1              LVL(kmax*3),
     1              I,J,K,
     1              ISTATUS
C
       integer      lend

       REAL       PHI(IMAX,JMAX,KMAX),
     1              U(IMAX,JMAX,KMAX),
     1              V(IMAX,JMAX,KMAX),
     1              OM(IMAX,JMAX,KMAX),
     1              temp(imax,jmax,kmax),
     1              bal(imax,jmax,kmax*3),
     1              rh(imax,jmax,kmax),
     1              sh(imax,jmax,kmax),
     1              p(kmax)
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
                IP(K)=int(p(k)/100.)
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
c
c  Next, the rh.
c
       EXT='lh3'
       dir=directory(1:lend)//'balance/'//ext(1:3)//'/'
       DO K=1,KMAX
              VAR(K)='RHL'
              LVL(K)=IP(K)
              LVL_COORD(K)='MB  '
              UNITS(K)='Percent'
              COMMENT(K)='Balanced RH'
       ENDDO

       write(6,*)' Writing grids ',ext(1:3),' ',fname9

       call write_laps_data(i4time,dir,ext,imax,jmax
     +,kmax,kmax,var,lvl,lvl_coord,units,comment,rh,istatus)

       IF (ISTATUS.ne.1) THEN
              PRINT*,'Error writing balanced rh data.'
              ISTATUS=0
              RETURN
       ENDIF
c
c  Finally, the specific humidity
c
       EXT='lq3'
       dir=directory(1:lend)//'balance/'//ext(1:3)//'/'
       DO K=1,KMAX
              VAR(K)='SH '
              LVL(K)=IP(K)
              LVL_COORD(K)='MB  ' 
              UNITS(K)='kg/kg'
              COMMENT(K)='Balanced specific humidity.'
       ENDDO

       write(6,*)' Writing grids ',ext(1:3),' ',fname9

       call write_laps_data(i4time,dir,ext,imax,jmax
     +,kmax,kmax,var,lvl,lvl_coord,units,comment,sh,istatus)

       IF (ISTATUS.ne.1) THEN
              PRINT*,'Error writing balanced sh data.'
              ISTATUS=0
              RETURN
       ENDIF                        
C
       ISTATUS=1
       RETURN
C
       END
!c
!============
! Hongli Jiang add to write sigma_ht output. 11/2/2011
!========
        subroutine write_bal_laps_ht(i4time,p3,u,v,temp,w,rh
     .  ,sh,imax,jmax,kmax,p,istatus)
C
       IMPLICIT       NONE
C
       INTEGER      I4TIME,
     1              IMAX,JMAX,KMAX,
     1              KMAX3,kmax2,
     1              IP(KMAX),
     1              LVL(kmax*3),
     1              I,J,K,
     1              ISTATUS
C
       integer      lend

       REAL        P3(IMAX,JMAX,KMAX),
     1              U(IMAX,JMAX,KMAX),
     1              V(IMAX,JMAX,KMAX),
     1              W(IMAX,JMAX,KMAX),
     1              temp(imax,jmax,kmax),
     1              bal(imax,jmax,kmax*3),
     1              rh(imax,jmax,kmax),
     1              sh(imax,jmax,kmax),
     1              p(kmax)
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
!HJ                IP(K)=int(p(k)/100.)
                IP(K)=int(p(k))
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
                VAR(K)='W3 '
                LVL(K)=IP(K-2*KMAX)
                LVL_COORD(K)='MB  '
                UNITS(K)='PA/S  '
                COMMENT(K)='Non-linear balanced w.           '
                DO J=1,JMAX
                DO I=1,IMAX
                        BAL(I,J,K)=W(I,J,K-2*KMAX)
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
              VAR(K)='P3 '
              LVL(K)=IP(K)
              LVL_COORD(K)='MB  '
              UNITS(K)='Meters'
              COMMENT(K)='Non-linear balanced height.'
              DO J=1,JMAX
              DO I=1,IMAX
                     BAL(I,J,K)=P3(I,J,K)
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
c
c  Next, the rh.
c
       EXT='lh3'
       dir=directory(1:lend)//'balance/'//ext(1:3)//'/'
       DO K=1,KMAX
              VAR(K)='RHL'
              LVL(K)=IP(K)
              LVL_COORD(K)='MB  '
              UNITS(K)='Percent'
              COMMENT(K)='Cloud liquid balanced rh.'
       ENDDO

       write(6,*)' Writing grids ',ext(1:3),' ',fname9

       call write_laps_data(i4time,dir,ext,imax,jmax
     +,kmax,kmax,var,lvl,lvl_coord,units,comment,rh,istatus)

       IF (ISTATUS.ne.1) THEN
              PRINT*,'Error writing balanced rh data.'
              ISTATUS=0
              RETURN
       ENDIF
c
c  Finally, the specific humidity
c
       EXT='lq3'
       dir=directory(1:lend)//'balance/'//ext(1:3)//'/'
       DO K=1,KMAX
              VAR(K)='SH '
              LVL(K)=IP(K)
              LVL_COORD(K)='MB  ' 
              UNITS(K)='kg/kg'
              COMMENT(K)='Cloud liquid balanced specific humidity.'
       ENDDO

       write(6,*)' Writing grids ',ext(1:3),' ',fname9

       call write_laps_data(i4time,dir,ext,imax,jmax
     +,kmax,kmax,var,lvl,lvl_coord,units,comment,sh,istatus)

       IF (ISTATUS.ne.1) THEN
              PRINT*,'Error writing balanced sh data.'
              ISTATUS=0
              RETURN
       ENDIF                        
C
       ISTATUS=1
       RETURN
C
       END
