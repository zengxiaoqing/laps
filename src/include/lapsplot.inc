

	type namelist_type
          sequence

          real continent_line_width
          real time_zone

          integer latlon_int

          character*3  c3_time_zone

	end type namelist_type

        type (namelist_type) namelist_parms
