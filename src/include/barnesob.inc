

        integer*4 max_var
        parameter (max_var = 2)       ! max vars (accomodates U & V)

	type barnesob
          sequence
	  integer i              
	  integer j              
	  integer k              
          real value(max_var)         ! Values of incremental obs
          real weight
          real elev                   ! Used for sfc obs
          real ldf                    ! Land fraction for surface obs
          integer mask_sea            ! Used for sfc analyses
	end type barnesob

	type barnesob_qc
          sequence
	  integer i              
	  integer j              
	  integer k              
          real value(max_var)         ! Values of incremental obs
          real weight
          real elev                   ! Used for sfc obs
          real ldf                    ! Land fraction for surface obs
          integer mask_sea            ! Used for sfc analyses
          logical qc
	end type barnesob_qc

