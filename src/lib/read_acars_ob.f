
        subroutine read_acars_ob(lun,c_obtype,xlat,xlon,elev,arg1,arg2       
     1                                   ,asc9_tim_pirep,iwrite,l_eof)       

        real*4 elev ! meters
        real*4 dd   ! degrees (99999. is missing)
        real*4 ff   ! meters/sec (99999. is missing)

        character*9 asc9_tim_pirep,asc9_tim_rcvd
        character*80 string
        character*(*)c_obtype

        logical l_eof

        dd = 99999.
        ff = 99999.

        l_eof = .false.

5       read(lun,101,end=900,err=5)string(1:6)
101     format(a6)

        if(string(2:5) .eq. 'Time')then
!           a9time = string(30:39)
            read(lun,151)asc9_tim_pirep,asc9_tim_rcvd
151         format(1x,a9,2x,a9)
            if(iwrite .eq. 1)write(6,151)asc9_tim_pirep,asc9_tim_rcvd       
        endif

        if(string(2:4) .eq. 'Lat')then
            read(lun,201)xlat,xlon,elev
201         format(2(f8.3,2x), f6.0,2i5)
        endif

        if(string(2:5) .eq. 'Wind' .and. c_obtype .eq. 'wind')then
            read(lun,202)idir_deg,ff
 202        format (1x, i3,7x, f6.1)
 220        format (' ', i3, ' deg @ ', f6.1, ' m/s')
            if(iwrite .eq. 1)write(6,220)idir_deg,ff
            dd = idir_deg
            arg1 = dd
            arg2 = ff
            return
        endif

        if(string(2:5) .eq. 'Temp' .and. c_obtype .eq. 'temp')then
            read(lun,302)temp
 302        format (1x,f10.1)
            write(6,302)temp
            arg1 = temp
            return
        endif

!       if(string(2:5) .eq. 'Clou')then
!           do i = 1,3
!               read(lun,203,err=500)cbase_ft,ctop_ft,icover
!203            format (12x,2f8.0,i5)
!           enddo ! i cloud layer
!       endif ! Cloud Report String

500     goto5

900     l_eof = .true.

        return

        end
