

C       PROF_CDF_COMMON.FOR
C
C Include file for PROF CDF COMMON's
C
        logical open_list(4)
        integer cdfid_list(4),wmoid_list(max_profilers,4),nrec_list(4)
        integer error_flag,res_len(4),recdim_list(4)
        character*6 staname_list(max_profilers,4)
        common/prof_cdf_common/open_list,cdfid_list,staname_list,
     $                         wmoid_list,nrec_list,error_flag,
     $                         res_len,recdim_list
C
