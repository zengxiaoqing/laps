
        subroutine avapsread_sub(lun)

        integer maxlvl
        parameter (maxlvl=10000)

        real*4 elapsed_time(maxlvl)
        real*4 p(maxlvl),  p_out(maxlvl)
        real*4             ht_out(maxlvl)
        real*4 t(maxlvl),  t_out(maxlvl)
        real*4 td(maxlvl), td_out(maxlvl)
        real*4 rh(maxlvl)
        real*4 u(maxlvl),  u_out(maxlvl)
        real*4 v(maxlvl),  v_out(maxlvl)
        real*4 ws(maxlvl), ws_out(maxlvl)
        real*4 wd(maxlvl), wd_out(maxlvl)
        real*4 dz(maxlvl)
        real*4 lon(maxlvl)
        real*4 lat(maxlvl)
        real*4 rng(maxlvl)
        real*4 az(maxlvl)
        real*4 alt(maxlvl)
        real*4 qp(maxlvl)
        real*4 qt(maxlvl)
        real*4 qh(maxlvl)
        real*4 qu(maxlvl)
        real*4 qv(maxlvl)
        real*4 quv(maxlvl)

        real*4 lat_s,lon_s

        real*4 lat_a(maxlvl),lon_a(maxlvl)

        character*172 header_line(15)

        character*9 a9time_a(maxlvl), a9time_ob

        call get_r_missing_data(r_missing_data,istatus)
        
        do i=1,15
          read (lun,99) header_line(i)
99        format (a172)
        enddo
        j = 1
2       read (lun,100,end=105)elapsed_time(j),P(j),T(j),Td(j),RH(j),
     1                        u(j),v(j),WS(j),WD(j),dZ(j),lon(j),
     2                        lat(j),rng(j),az(j),alt(j),Qp(j),Qt(j),
     3                        Qh(j),Qu(j),Qv(j),Quv(j)
100     format (f6.1,f7.1,3f6.1,2f7.1,3f6.1,f9.3,f8.3,2f6.1,f8.1,
     1          6f5.1)
        if (lon(j).ne.999.0 .and. lat(j).ne.999.0) then
          hlon = -lon(j)
          hlat = lat(j)
        endif
c
        j = j + 1
        go to 2
c
105     j = j - 1

        write(6,*)' Completed sounding read, # of levels = ',j

!       Get Launch Location
        read(header_line(4),4)lat_s,lon_s,ialt
        write(6,*)' lat_s,lon_s,ialt',lat_s,lon_s,ialt
 4      format(60x,f11.0,1x,f9.0,1x,i5)

        elev_s = ialt

!       Get Launch Time
        read(header_line(5),5)iy,mon,id,ih,min,is
 5      format(35x,i4,2x,i2,2x,i2,1x,3(1x,i2))
        write(6,*)'iy,mon,id,ih,min,is',iy,mon,id,ih,min,is

        i4time_1960 = I4TIME_INT_LP (iy,mon,id,ih,min,is,istatus)
        i4time_1970 = i4time_1960 - 315619200 ! Convert to LAPS i4time

        i4time_launch = i4time_1970
        write(6,*)'i4time_launch = ',i4time_launch

        nlevels = j

        write(6,*)' lvl  i4time   p   t   td'

        lvl_out = 0

        do lvl = 1,nlevels
            if(p(lvl) .ne. 9999.)then
                lvl_out = lvl_out + 1
                i4time = i4time_launch + nint(elapsed_time(lvl))
                call make_fnam_lp (I4TIME, a9time_a(lvl_out), ISTATUS)
                p_out(lvl_out) = p(lvl)*100.
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

                write(6,*)lvl,lvl_out,i4time_lvl,a9time_a(lvl_out)
     1                   ,p_out(lvl_out)
     1                   ,t_out(lvl_out),td_out(lvl_out)
     1                   ,wd_out(lvl_out),ws_out(lvl_out)       
            endif
        enddo ! lvl

!       Call write_snd for this sounding
        lat_a = lat_s
        lon_a = lon_s

        call write_snd    (11                              ! I
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

        return
        end
