
!       This is a poor man's data structure for raw temperature observations

!       Note that the integer grid values are primarily used in the analysis.
!       The floating point grid values are useful for verification, and just
!       in case we need them for increased accuracy in the analysis.

        integer n_fields
        parameter (n_fields=12)

        real temp_obs(max_obs,n_fields)

        integer i_ri,i_rj,i_rk,i_i,i_j,i_k,i_ob_raw,i_ob_grid
        integer i_bias,i_wt,i_inst_err
        parameter (i_ri = 1)
        parameter (i_rj = 2)
        parameter (i_rk = 3)
        parameter (i_i = 4)
        parameter (i_j = 5)
        parameter (i_k = 6)
        parameter (i_ob_raw = 7)
        parameter (i_ob_grid = 8)
        parameter (i_bias = 9)
        parameter (i_wt = 10)
        parameter (i_inst_err = 11)
        parameter (i_obstype = 12)
