           
      subroutine write_snd(lun_out                         ! I
     1                    ,maxsnd,maxlvl,nsnd              ! I
     1                    ,iwmostanum                      ! I
     1                    ,stalat,stalon,staelev           ! I
     1                    ,c5_staid,a9time_ob,c8_obstype   ! I
     1                    ,nlvl                            ! I
     1                    ,height_m                        ! I
     1                    ,pressure_pa                     ! I
     1                    ,temp_c                          ! I
     1                    ,dewpoint_c                      ! I
     1                    ,dir_deg                         ! I
     1                    ,spd_mps                         ! I
     1                    ,istatus)                        ! O

!     Steve Albers FSL    2001

!     Write routine for 'snd' file

!     For missing data values, 'r_missing_data' should be passed in 

!.............................................................................

      integer iwmostanum(maxsnd),nlvl(maxsnd)
      real stalat(maxsnd,maxlvl),stalon(maxsnd,maxlvl),staelev(maxsnd)       
      character c5_staid(maxsnd)*5,a9time_ob(maxsnd,maxlvl)*9
     1         ,c8_obstype(maxsnd)*8,c_line*150

      real height_m(maxsnd,maxlvl)
      real pressure_pa(maxsnd,maxlvl)
      real temp_c(maxsnd,maxlvl)
      real dewpoint_c(maxsnd,maxlvl)
      real dir_deg(maxsnd,maxlvl)
      real spd_mps(maxsnd,maxlvl)

!............................................................................

      do isnd = 1,nsnd

!       Write Sounding Header

        write(6,511,err=990)
     1             iwmostanum(isnd),nlvl(isnd)
     1            ,stalat(isnd,1),stalon(isnd,1),staelev(isnd)
     1            ,c5_staid(isnd),a9time_ob(isnd,1),c8_obstype(isnd)

        write(lun_out,511,err=990)
     1             iwmostanum(isnd),nlvl(isnd)
     1            ,stalat(isnd,1),stalon(isnd,1),staelev(isnd)
     1            ,c5_staid(isnd),a9time_ob(isnd,1),c8_obstype(isnd)

  511   format(i12,i12,f11.4,f15.4,f15.0,1x,a5,3x,a9,1x,a8)

        do lvl = 1,nlvl(isnd)

!         Write Sounding Level

          write(c_line,*)height_m(isnd,lvl),pressure_pa(isnd,lvl)
     1              ,temp_c(isnd,lvl)
     1              ,dewpoint_c(isnd,lvl)
     1              ,dir_deg(isnd,lvl),spd_mps(isnd,lvl)

          call s_len2(c_line,len_line)
          write(lun_out,521)c_line(1:len_line)
  521     format(a)

          if(isnd .le. 100)then
              write(c_line,*)height_m(isnd,lvl),pressure_pa(isnd,lvl)
     1              ,temp_c(isnd,lvl)
     1              ,dewpoint_c(isnd,lvl)
     1              ,dir_deg(isnd,lvl),spd_mps(isnd,lvl)
     1              ,a9time_ob(isnd,lvl)
     1              ,stalat(isnd,lvl),stalon(isnd,lvl)
              call s_len2(c_line,len_line)
              write(6,521)c_line(1:len_line)
          endif

        enddo ! lvl
      enddo ! isnd

      go to 999

 990  write(6,*)' ERROR in write_snd'
      istatus=0
      return

 999  istatus = 1
      return
      end

