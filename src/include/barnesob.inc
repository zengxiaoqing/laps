

        integer*4 max_var
        parameter (max_var = 2)       ! max vars (accomodates U & V)

	type barnesob
          sequence
	  integer i              
	  integer j              
	  integer k              
          real value(max_var)         ! Values of incremental obs
          real weight
	end type barnesob

	type barnesob_qc
          sequence
	  integer i              
	  integer j              
	  integer k              
          real value(max_var)         ! Values of incremental obs
          real weight
          logical qc
	end type barnesob_qc

