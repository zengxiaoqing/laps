           
      subroutine write_pro(lun_out                         ! I
     1                    ,maxpro,maxlvl,npro              ! I
     1                    ,iwmostanum                      ! I
     1                    ,stalat,stalon,staelev           ! I
     1                    ,c6_staid,a9time_ob,c8_obstype   ! I
     1                    ,nlvl                            ! I
     1                    ,height_m                        ! I
     1                    ,dir_deg                         ! I
     1                    ,spd_mps                         ! I
     1                    ,istatus)                        ! O

!     Steve Albers FSL    2001

!     Write routine for 'pro' file

!     For missing data values, 'r_missing_data' should be passed in 

!.............................................................................

      integer iwmostanum(maxpro),nlvl(maxpro)
      real stalat(maxpro),stalon(maxpro),staelev(maxpro)       
      character c6_staid(maxpro)*5,a9time_ob(maxpro)*9
     1         ,c8_obstype(maxpro)*8,c_line*150

      real height_m(maxpro,maxlvl)
      real dir_deg(maxpro,maxlvl)
      real spd_mps(maxpro,maxlvl)

      logical l_good_level(maxlvl)
!............................................................................

      call get_r_missing_data(r_missing_data,istatus)

      rms = 1.0

      do ipro = 1,npro

!       Count and index the good levels
        n_good_levels = 0
        do ilvl = 1,nlvl(ipro)
            if(dir_deg(ipro,ilvl) .ne. r_missing_data .and.
     1         spd_mps(ipro,ilvl) .ne. r_missing_data       )then
                n_good_levels = n_good_levels + 1
                l_good_level(ilvl) = .true.
            else
                l_good_level(ilvl) = .false.
            endif
        enddo ! ilvl

!       Write Sounding Header

        call filter_string(c6_staid(ipro))

        write(6      ,401)iwmostanum,n_good_levels,stalat(ipro)
     1                   ,stalon(ipro),staelev(ipro),c6_staid(ipro)
     1                   ,a9time_ob(ipro),c8_obstype(ipro)
        write(lun_out,401)iwmostanum,n_good_levels,stalat(ipro)
     1                   ,stalon(ipro),staelev(ipro),c6_staid(ipro)
     1                   ,a9time_ob(ipro),c8_obstype(ipro)
401     format(i12,i12,f11.3,f15.3,f15.0,5x,a6,3x,a9,1x,a8)

        do ilvl = 1,nlvl(ipro)

          if(l_good_level(ilvl))then

!             Write Profile Level
              write(lun_out,301,err=303)height_m(ipro,ilvl)
     1                                 ,dir_deg(ipro,ilvl)
     1                                 ,spd_mps(ipro,ilvl)
     1                                 ,rms       
301           format(1x,f6.0,f6.0,2f6.1,3f7.1)
303           continue 

              if(ipro .le. 100)then
                  write(6,301,err=313)height_m(ipro,ilvl)
     1                               ,dir_deg(ipro,ilvl)
     1                               ,spd_mps(ipro,ilvl)
     1                               ,rms       
313               continue
              endif

          endif ! good level

        enddo ! ilvl

      enddo ! ipro

      go to 999

 990  write(6,*)' ERROR in write_pro'
      istatus=0
      return

 999  istatus = 1
      return
      end

