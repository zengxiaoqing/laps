


      subroutine get_static_info(c_vars_req,c_values_req_in,n_vars_req      
     1                                                      ,istatus)

      include 'lapsparms.cmn'

!     This routine accesses the run time variables from the 'nest7grid.parms'
!     file via the common block or by reading the actual file.

!     May 1996 Steve Albers - FSL

! ----------------- ARGUMENT LIST ---------------------------------------------

!     Array of 1 or more requested runtime parameter names
      character*40 c_vars_req                             ! Input

!     Array of runtime parameter values - stuffed into character variables.
!     These need to be decoded into integer, real, etc, by the caller.
      character*(*) c_values_req_in                       ! Output
      character*150 c_values_req                          ! Local

!     Number of variables requested
      integer n_vars_req                                  ! Input

      integer istatus                                     ! Output

! -----------------------------------------------------------------------------

      integer max_vars_list
      parameter (max_vars_list = 25)

      character*40 c_vars_list(max_vars_list)

! -----------------------------------------------------------------------------

!     Let's make sure we have the run time variables in the common block

      call get_config(istatus)
      if(istatus .ne. 1)then
          write(6,*)' Error detected in calling get_config'
          istatus = 0
          return
      endif

      if(n_vars_req .ne. 1)then
          write(6,*)' Error detected in get_static_info, n_vars_req = '       
     1              ,n_vars_req
          istatus = 0
          return
      endif

!     List the variable names
      c_vars_list(1) = 'laps_cycle_time'
      c_vars_list(2) = 'i_delta_sat_t_sec'
      c_vars_list(3) = 'r_msng_sat_flag_cdf'

      c_vars_list(4) = 'radarext_3d_accum'
      c_vars_list(5) = 'radarext_3d'

      c_vars_list(6) = 'path_to_raw_pirep'
      c_vars_list(7) = 'path_to_raw_rass'
      c_vars_list(8) = 'path_to_raw_profiler'
      c_vars_list(9) = 'path_to_raw_blprass'
      c_vars_list(10) = 'path_to_raw_blpprofiler'
      c_vars_list(11) = 'path_to_raw_satellite_cdf'
      c_vars_list(12) = 'path_to_ruc'
      c_vars_list(13) = 'path_to_ngm'
      c_vars_list(14) = 'path_to_wsi_2d_radar'
      c_vars_list(15) = 'path_to_wsi_3d_radar'
      c_vars_list(16) = 'path_to_raw_satellite_gvr'
      c_vars_list(17) = 'path_to_qc_acars'
      c_vars_list(18) = 'path_to_raw_raob'
      c_vars_list(19) = 'r_msng_sat_flag_gvr'
      c_vars_list(20) = 'r_msng_sat_flag_asc'
      c_vars_list(21) = 'path_to_raw_sat_wfo_vis'
      c_vars_list(22) = 'path_to_raw_sat_wfo_i39'
      c_vars_list(23) = 'path_to_raw_sat_wfo_iwv'
      c_vars_list(24) = 'path_to_raw_sat_wfo_i11'
      c_vars_list(25) = 'path_to_raw_sat_wfo_i12'

      do i_req = 1,n_vars_req
          call s_len(c_vars_req,len_req)

          do i_list = 1,max_vars_list
              call s_len(c_vars_list(i_list),len_list)

              if(c_vars_req(1:len_req) .eq.
     1           c_vars_list(i_list)(1:len_list))then ! we have a match

               !  Write this variable into the character array

                  if(i_list .eq. 1)goto1
                  if(i_list .eq. 2)goto2
                  if(i_list .eq. 3)goto3
                  if(i_list .eq. 4)goto4
                  if(i_list .eq. 5)goto5
                  if(i_list .eq. 6)goto6
                  if(i_list .eq. 7)goto7
                  if(i_list .eq. 8)goto8
                  if(i_list .eq. 9)goto9
                  if(i_list .eq. 10)goto10
                  if(i_list .eq. 11)goto11
                  if(i_list .eq. 12)goto12
                  if(i_list .eq. 13)goto13
                  if(i_list .eq. 14)goto14
                  if(i_list .eq. 15)goto15
                  if(i_list .eq. 16)goto16
                  if(i_list .eq. 17)goto17
                  if(i_list .eq. 18)goto18
                  if(i_list .eq. 19)goto19
                  if(i_list .eq. 20)goto20
                  if(i_list .eq. 21)goto21
                  if(i_list .eq. 22)goto22
                  if(i_list .eq. 23)goto23
                  if(i_list .eq. 24)goto24
                  if(i_list .eq. 25)goto25

1                 write(c_values_req,'(i10,90x)')
     1                  laps_cycle_time
                  go to 900

2                 write(c_values_req,'(i10,90x)')
     1                  i_delta_sat_t_sec
                  go to 900

3                 write(c_values_req,'(f10.3,90x)')
     1                  r_msng_sat_flag_cdf
                  go to 900

4                 write(c_values_req,'(a8)')
     1                  radarext_3d_accum
                  go to 900

5                 write(c_values_req,'(a8)')
     1                  radarext_3d
                  go to 900

6                 write(c_values_req,'(a150)')
     1                  path_to_raw_pirep
                  go to 900

7                 write(c_values_req,'(a150)')
     1                  path_to_raw_rass
                  go to 900

8                 write(c_values_req,'(a150)')
     1                  path_to_raw_profiler
                  go to 900

9                 write(c_values_req,'(a150)')
     1                  path_to_raw_blprass
                  go to 900

10                write(c_values_req,'(a150)')
     1                  path_to_raw_blpprofiler
                  go to 900

11                write(c_values_req,'(a150)')
     1                  path_to_raw_satellite_cdf
                  go to 900

12                write(c_values_req,'(a150)')
     1                  path_to_ruc
                  go to 900

13                write(c_values_req,'(a150)')
     1                  path_to_ngm
                  go to 900

14                write(c_values_req,'(a150)')
     1                  path_to_wsi_2d_radar
                  go to 900

15                write(c_values_req,'(a150)')
     1                  path_to_wsi_3d_radar
                  go to 900

16                write(c_values_req,'(a150)')
     1                  path_to_raw_satellite_gvr
                  go to 900

17                write(c_values_req,'(a150)')
     1                  path_to_qc_acars
                  go to 900

18                write(c_values_req,'(a150)')
     1                  path_to_raw_raob
                  go to 900

19                write(c_values_req,'(f10.3,90x)')
     1                  r_msng_sat_flag_gvr
                  go to 900

20                write(c_values_req,'(f10.3,90x)')
     1                  r_msng_sat_flag_asc
                  go to 900

21                write(c_values_req,'(a150)')
     1                  path_to_raw_sat_wfo_vis
                  go to 900

22                write(c_values_req,'(a150)')
     1                  path_to_raw_sat_wfo_i39
                  go to 900

23                write(c_values_req,'(a150)')
     1                  path_to_raw_sat_wfo_iwv
                  go to 900

24                write(c_values_req,'(a150)')
     1                  path_to_raw_sat_wfo_i11
                  go to 900

25                write(c_values_req,'(a150)')
     1                  path_to_raw_sat_wfo_i12
                  go to 900

              endif

          enddo ! i_list

        ! No match
          write(6,*)' get_static_info: no match for variable '
     1             ,c_vars_req
          write(6,*)' get_static_info: returning on error '
          istatus = 0
          return

 900      c_values_req_in = c_values_req
      
      enddo ! i_req

      istatus = 1
      return
      end
