
        integer mxopt, itermx
	parameter (mxopt = 10)
        parameter (itermx = 10000)

	integer cv,cv_ref,modeop,iswp,g2,iterop,itervar,init_optmiz
     1	       ,g6,iterswp,nv,swps_per_extrpl
        integer iv,lswp,kswp,k_ref,qq,num_stepsize_inc,u7

	logical l_aitken,l_aitken_thistime
     1	       ,l_depth,k_,l_swp_extrpl,l_swp_extpl_end

	real*8  j9,l6,l7,m9,increment_ratio,mlow,mhigh,l_m1,l_m2
        integer idir_opt,k3,k4,k5,k8

	common /optmiz/ alow,amid,ahigh,anew,aalow,aamid,aahigh
     1                 ,a8,a9,bottom_last
     1	               ,l_aitken,l_aitken_thistime
     1                 ,l_depth,cfive,c7,dstep_swp
     1	               ,f_ref,f9,g_ref,gone,g2,g5,g6,g8
     1                 ,conv_rate,conv_rate_m1,conv_rate_m2
     1                 ,h9,k_,l_swp_extrpl
     1	               ,iterswp,iterop,itervar,idir_opt,j9,kswp
     1	               ,swps_per_extrpl,k3,k4,k5,k8
     1                 ,k_ref
     1                 ,l6,l7,m9,iv,nv,qq,iswp,lswp
     1	               ,u7,vdepth(mxopt),num_stepsize_inc
     1                 ,increment_ratio,v5,step_factor
     1	               ,f_m1,f_m2,l_m1,l_m2,l_swp_extpl_end
	common /optmiz/ a_ref(mxopt),dstep(mxopt),cv(mxopt)
     1	               ,cv_ref(mxopt)
	common /optmiz/ c1(mxopt),d1(mxopt),e1(mxopt),g1(mxopt),crv(mxopt)
	common /optmiz/ h1(mxopt)
	common /optmiz/ mlow(mxopt),mhigh(mxopt),za(mxopt)
	common /optmiz/ modeop,aref
        common /optmiz/ init_optmiz
