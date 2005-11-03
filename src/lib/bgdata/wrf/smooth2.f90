  SUBROUTINE smooth2 (nx,ny,dist,data)

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: nx,ny,dist
    REAL,    INTENT(INOUT) :: data(nx,ny)
    REAL                   :: tempdata(nx,ny)

    INTEGER i,j,ii,jj,ipt,jpt
    REAL  :: np, datasum

    tempdata(:,:) = 0.
    DO j = 1, ny
      DO i = 1, nx
        
        np = 0.
        datasum = 0.
        DO jj = -dist,dist,1
          DO ii = -dist,dist,1
            ipt = i + ii
            jpt = j + jj
            IF ((ipt .GE. 1) .AND. (ipt .LE. nx) .AND. &
                (jpt .GE. 1) .AND. (jpt .LE. ny) ) THEN
               datasum = datasum + data(ipt,jpt)
               np = np + 1.
            ENDIF
          ENDDO
        ENDDO

        tempdata(i,j) = datasum / np
      ENDDO
    ENDDO
    data = tempdata
    END SUBROUTINE smooth2
      

