

      function get_altitude()
      integer get_altitude          ! Site altitude (meters * 100000)

      return
      end


      function get_latitude()
      integer get_latitude          ! Site latitude (meters * 100000)

      return
      end


      function get_longitude()
      integer get_longitude         ! Site longitude (meters * 100000)

      return
      end


      function get_field_num(ext)
      integer get_field_num
      character*3 ext

      if(ext .eq. 'DBZ')get_field_num = 1
      if(ext .eq. 'VEL')get_field_num = 2

      return
      end


      function read_radial()

      return
      end


      function get_status()

      return
      end


      function get_fixed_angle()
      integer get_fixed_angle     ! Beam tilt angle (degrees * 100)

      return
      end


      function get_scan()
      integer get_scan            ! Scan #

      return
      end


      function get_tilt()
      integer get_tilt            ! Tilt #

      return
      end


      function get_year()

      return
      end


      function get_month()

      return
      end

      function get_day()

      return
      end


      function get_hour()

      return
      end


      function get_min()

      return
      end


      function get_sec()

      return
      end


      function get_vcp()
      integer get_vcp

      get_vcp = 0

      return
      end


      function get_azi()

      return
      end


      function get_nyquist()
      integer get_nyquist        ! Nyquist velocity of the radial (M/S*100)

      return
      end


      function get_number_of_gates(index)
      integer get_number_of_gates

      return
      end


      function get_first_gate()

      return
      end


      function get_data_field(index, n_gates)

      return
      end


      function cvt_fname_data()

      return
      end


