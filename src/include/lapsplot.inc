

	type namelist_type
          sequence

          logical l_discrete

          real continent_line_width
          real time_zone

          integer latlon_int

          character*3 c3_time_zone
          character*9 c_institution
          character*6 c_vnt_units
          character*7 c_units_type
          character*7 c_pbl_depth_units

	end type namelist_type

        type (namelist_type) namelist_parms
