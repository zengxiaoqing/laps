
        subroutine avapsread_sub(filename, lun_out
     1                          ,i4time_drpsnd_earliest
     1                          ,i4time_drpsnd_latest,istatus)

        integer maxlvl
        parameter (maxlvl=10000)

        real elapsed_time(maxlvl)
        real p(maxlvl),  p_out(maxlvl)          ! millibars
        real             ht_out(maxlvl)
        real t(maxlvl),  t_out(maxlvl)
        real td(maxlvl), td_out(maxlvl)
        real rh(maxlvl)
        real u(maxlvl),  u_out(maxlvl)
        real v(maxlvl),  v_out(maxlvl)
        real ws(maxlvl), ws_out(maxlvl)
        real wd(maxlvl), wd_out(maxlvl)
        real dz(maxlvl)
        real lon(maxlvl)
        real lat(maxlvl)
        real rng(maxlvl)
        real az(maxlvl)
        real alt(maxlvl)
        real qp(maxlvl)
        real qt(maxlvl)
        real qh(maxlvl)
        real qu(maxlvl)
        real qv(maxlvl)
        real quv(maxlvl)

        real lat_s,lon_s

        real lat_a(maxlvl),lon_a(maxlvl)

        character*172 header_line(15)
        character*130 line

        character*9 a9time_a(maxlvl), a9time_ob
        character*5 c5_staid
        character*(*)  filename

        logical l_parse

        i_snd_out = 0

        lun = 1
        open ( lun, file=filename, status='old', err=1000 )

        call get_r_missing_data(r_missing_data,istatus)
        if ( istatus .ne. 1 ) then
           write (6,*) ' Error getting r_missing_data'
           return
        endif
        
        nh = 15

!       ---- read the top avaps header ----        

40      do i=1,nh
          read (lun,99,end=1000) header_line(i)
          print*, i,header_line(i)(1:120)
99        format (a)
        enddo

!       ---- read the whole data file  ---- 

        j = 1

        write(6,*)' Reading levels from the sounding'

 2      read (lun,99,end=105)line(1:20)
        if(l_parse(line,'Data') .and. l_parse(line,'Type'))then
            write(6,*)' New Sounding Detected - backspace lun'
            backspace lun
            go to 105

        elseif(l_parse(line,'Project') .and. l_parse(line,'ID'))then
            write(6,*)' New Sounding Detected - backspace lun'
            write(6,*)' WARNING: already reading second header line'
            write(6,*)'          an extra backspace will be done'
            backspace lun
            backspace lun
            go to 105

        else ! reread full line
            backspace lun
            read (lun,99,end=105)line

        endif

        write(6,*)line(1:130)
        read (line,100,end=105)elapsed_time(j),P(j),T(j),Td(j),RH(j),
     1                        u(j),v(j),WS(j),WD(j),dZ(j),lon(j),
     2                        lat(j),rng(j),az(j),alt(j),Qp(j),Qt(j),
     3                        Qh(j),Qu(j),Qv(j),Quv(j)
100     format (f6.1,f7.1,3f6.1,2f7.1,3f6.1,f9.3,f8.3,2f6.1,f8.1,
     1          6f5.1)
        if (lon(j).ne.999.0 .and. lat(j).ne.999.0) then
          hlon = -lon(j)
          hlat = lat(j)
        endif

        if(elapsed_time(j) .le. 5.)then
            write(6,*)'Nearing end of sounding'
        endif
c
        j = j + 1
        go to 2
c
105     j = j - 1

        write(6,*)' Completed sounding read, # of levels = ',j

!       Get Launch Location
        read(header_line(4),4)lon_s,lat_s,ralt
        write(6,*)' lat_s,lon_s,ralt',lat_s,lon_s,ralt
 4      format(59x,f8.0,1x,f7.0,1x,f8.0)

!       if(ralt .eq. 9999.)then
            staelev = -999.        
!       else
!           staelev = ralt
!       endif

!       Get Launch Time
        read(header_line(5),5)iy,mon,id,ih,min,is
 5      format(35x,i4,2x,i2,2x,i2,1x,3(1x,i2))
        write(6,*)'iy,mon,id,ih,min,is',iy,mon,id,ih,min,is

        i4time_1960 = I4TIME_INT_LP (iy,mon,id,ih,min,is,istatus)
!       i4time_1970 = i4time_1960 - 315619200 ! Convert from LAPS i4time

        i4time_launch = i4time_1960
        write(6,*)'i4time_launch = ',i4time_launch

        nlevels = j

        do isort = 1,2
          write(6,*)' isort = ',isort
          write(6,*)' lvl  i4time   p   t   td'
          lvl_out = 0
          isort_status = 1

          if(isort .eq. 1)then
              lvlb = nlevels
              lvle = 1
              lvli = -1
          else
              lvlb = 1
              lvle = nlevels
              lvli = +1
          endif

          do lvl = lvlb,lvle,lvli
            if(p(lvl) .ne. 9999.)then
                lvl_out = lvl_out + 1
                i4time = i4time_launch + nint(elapsed_time(lvl))
                call make_fnam_lp (I4TIME, a9time_a(lvl_out), ISTATUS)
                p_out(lvl_out) = p(lvl) ! *100.
                ht_out(lvl_out) = r_missing_data

                if(t(lvl) .eq. 999.)then
                    t_out(lvl_out) = r_missing_data
                else
                    t_out(lvl_out) = t(lvl)
                endif

                if(td(lvl) .eq. 999.)then
                    td_out(lvl_out) = r_missing_data
                else
                    td_out(lvl_out) = td(lvl)
                endif

                if(wd(lvl) .eq. 999. .or. ws(lvl) .eq. 999.)then
                    wd_out(lvl_out) = r_missing_data
                    ws_out(lvl_out) = r_missing_data
                else
                    wd_out(lvl_out) = wd(lvl)
                    ws_out(lvl_out) = ws(lvl)
                endif

                write(6,*)lvl,lvl_out,i4time,a9time_a(lvl_out)
     1                   ,p_out(lvl_out)
     1                   ,t_out(lvl_out),td_out(lvl_out)
     1                   ,wd_out(lvl_out),ws_out(lvl_out)       

                if(lvl_out .gt. 1)then
                    if(p_out(lvl_out) .ge. p_out(lvl_out-1))then
                        if(isort_status .eq. 1)then
                            write(6,*)
     1                    ' NOTE: AVAPS pressures increase with height'       
                            write(6,*)lvl_out,p_out(lvl_out)
     1                                       ,p_out(lvl_out-1)
                        endif

                        isort_status = 0

                        if(isort .eq. 2)then
                            write(6,*)' ERROR: sounding not monotonic'
                            istatus = 0
                            return
                        endif
                    endif ! pressure trend not monotonically decreasing
                endif

            endif ! valid pressure at this level
          enddo ! lvl

          if(isort_status .eq. 1)then
              write(6,*)' Good sort of levels, we will proceed'
              go to 200
          endif

        enddo ! isort

200     continue

        if(i4time_launch .ge. i4time_drpsnd_earliest .and.
     1     i4time_launch .le. i4time_drpsnd_latest        )then

!           Call write_snd for this sounding
            lat_a = lat_s
            lon_a = lon_s

            i_snd_out = i_snd_out + 1
            write(c5_staid,201)i_snd_out
201         format('AVP',i2.2)

            call write_snd(lun_out                         ! I
     1                    ,1,lvl_out,1                     ! I
     1                    ,iwmostanum                      ! I
     1                    ,lat_a,lon_a,staelev             ! I
     1                    ,c5_staid,a9time_a,'DROPSND '    ! I
     1                    ,lvl_out                         ! I
     1                    ,ht_out                          ! I
     1                    ,p_out                           ! I
     1                    ,t_out                           ! I
     1                    ,td_out                          ! I
     1                    ,wd_out                          ! I
     1                    ,ws_out                          ! I
     1                    ,istatus)                        ! O

        else
            write(6,*)' This dropsonde is outside time window'

        endif ! within time window

        write(6,*)' Looping back to look for rest of new sounding'
        go to 40

1000    write(6,*)' End of file reached'
        
        return
        end
