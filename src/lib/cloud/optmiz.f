
	subroutine optimize(a,f_merit,depth_this_run)

	implicit real*8 (a-z)
	include 'optinc.inc'
!       include 'global.inc'
        integer j
        logical l_extra_output,l_extra1 /.true./
	
	logical l_c        

        character*9 c9_outstring
	real*8 a(mxopt)

	data swp_ratio/1.5/

!	Output Logic
	if(iswp .le. 1)then
            l_extra_output = l_extra1 
        else
            l_extra_output = .false.
        endif
        

        if(modeop .ge. 10 .OR. 
     1	     l_extra_output .eqv. .true. 
     1	        .and.  ( iterop .lt. 1000
     1          .and.     kswp .le. 1 )
     1	        .and. iterop .ge. 1                       )then
            write(6,9835)iterop,modeop,iswp,iv,iterswp
     1                  ,cv(iv),dstep(iv),d0,f_merit
9835	    format(1x,'ip,mp,sw,iv,io,cv,dn,da,f'
     1     	  ,i6,5i3,5x,2d10.3,d14.6)
            write(6,9841)(a(j),j=1,nv)
9841        format(1x,'a(n)s ',7d12.5)
        endif

	if(init_optmiz .eq. 0)then
            modeop=1
            iswp=1
            g2=1
            iterop=0
            itervar=0
            g6=0
            iterswp=0
            do j = 1,nv
                mlow(j) = -1d10
                mhigh(j) = 1d10
                crv(j) = 0.
            enddo
            l_aitken = .true.
            l_depth = .true.
            k_ = .true.
            l_swp_extrpl = .false.
            l_swp_extpl_end = .false.
            k4 = 8
            swps_per_extrpl = 7
            init_optmiz = 1
        endif

	iterop = iterop+1
	itervar = itervar+1
	iterswp = iterswp+1
	if(itervar .gt. 250)goto 26509
	goto(21000,21600,22300,22600,23800,25200,25600,26508,27745
     1                                                ,27800),modeop

!****** Modeop = 1 *******************************************************
21000	step_factor = 10**depth_this_run*increment_ratio**num_stepsize_inc
	cfive = 1d-14*step_factor**2
	c7 = 1d+5/step_factor**2
	do j = 1,nv
	    c1(j) = 1
	    d1(j) = 0.
	    e1(j) = 0.
	    g1(j) = 0
	    cv(j) = 2
 	    cv_ref(j) = 2
	    dstep(j) = dstep(j)*step_factor
	enddo ! j
	v5 = num_stepsize_inc
	k8 = 0
	k3 = -1
	conv_rate = 0.
	conv_rate_m1 = 0.
	kswp = 1
	lswp = 1
	iv = 1
	f_ref = f_merit ! New Statement
	g_ref = f_merit
	gone = f_merit
	k5 = swps_per_extrpl
	swps_per_extrpl = k4
	u7 = 0
	g8 = 0

	if(abs(k5 - swps_per_extrpl) .ne. 1)then
	    h9 = -1
	else
	    h9 = k5 - swps_per_extrpl
	    k3 = 0
	endif

!****** Modeop = 2 ***********************************************************
21600	if(g6 .ne. 0)then
            g_ref = f_merit
            g6 = 0
        endif

        if(modeop .ne. 2 .or. itervar .le. 1)then
            a_ref(iv) = a(iv)
            f_ref = f_merit
        endif

        if(.not.(itervar .gt. 2 .or. kswp .lt. 3 
     1                          .or. c1(iv) .eq. 0))then
!           Multidimensional Convergence
            if(itervar .eq. 1 .and. kswp .ge. 3)then
     
                if(kswp .eq. 3)then
            	    a(iv) = a(iv) + d1(iv)
                else if(kswp .gt. 3)then
                    a(iv) = a(iv)+d1(iv)*2.-e1(iv)
                endif

                if(l_extra_output .eqv. .true. 
     1			        .and. iterop .lt. 1000)then
                    write(6,111)iv,a(iv)
111                 format(12x,'Accelerating a(',i2,') to',d15.6)
                endif

	        if( .not. ( abs( d1(iv)/dstep(iv) ) .gt. .4 
     1	                             .OR. 
     1               abs(za(iv)-a_ref(iv)) .gt. .4 * dstep(iv) ) )then       

                    if(l_extra1 .eqv. .true. .and. iterop .lt. 1000)
     1	                               write(6,*)'          aborted'
                    a(iv) = a_ref(iv)
                    goto21760

                endif

!               if(l_extra1 .eqv. .true. 
!	1			.and. iterop .lt. 1000)write(6,*)
                modeop = 2
                return

            endif ! itervar = 1 & kswp >= 3

            if(f_merit .gt. f_ref)then
                if(l_extra1 .eqv. .true. .and. iterop .lt. 2500)then
                    write(6,112)iswp,iv,itervar,f_merit-f_ref
112                 format(' Multidimensional Convergence halted',
     1                 ' swp ',i4,' ivar',i3,' iq ',i4,' f-fref',d10.3)       
                endif
                c1(iv) = 0
                a(iv) = a_ref(iv)
                f_merit = f_ref
            endif

            itervar = 1
        endif

21760   if(f_merit .gt. f_ref)goto26511
        a9 = 0
        anew = 1d10
21900	amid = f_merit
        aref = a(iv)
        a(iv) = a(iv) + dstep(iv)
        modeop = 3
        return

!****** Modeop = 3 ************************************************************
22300	ahigh = f_merit
        a(iv) = aref-dstep(iv)
        modeop = 4
        return

!****** Modeop = 4 ************************************************************
22600	alow = f_merit
	a(iv) = aref
	f_merit = amid
	dhigh = ahigh - amid
	dlow = alow - amid
	if(l_extra_output .eqv. .true.)write(6,102)alow,amid,ahigh
102     format(1x,' alow,amid,ahigh',3d15.5)
!       If step size is too small to affect F, increase it
	if(dhigh .eq. 0. .or. dlow .eq. 0.)then
	    dstep(iv) = dstep(iv)*5.
	    write(6,*)' dn(',iv,') is too small - increased to',dstep(iv)
	    goto21900
        endif

!       Converge parabolically if possible or else linearly
        if(abs(dhigh + dlow) .gt. 1d-6*abs(dhigh-dlow))then
            anew = a(iv)+dstep(iv)*(.5-dhigh/(dhigh+dlow))
        else
            anew = a(iv) - f_merit/((dhigh-dlow)/(2.*dstep(iv)))
	    write(6,*)' Linear Convergence cv(iv) = ',cv(iv),' modeop=10'
        endif

        if(l_extra_output .eqv. .true.)
     1      write(6,103)dhigh,dhigh+dlow,anew
103     format('  dh,crv,anew',3d15.5)

        if((anew - a(iv))*(dhigh-dlow) .gt. 0.)then
            write(6,*)' Parabolic convergence not well behaved'
            write(6,*)' a(iv),anew,alow,amid,ahigh,dh,dl',
 	1	        a(iv),anew,alow,amid,ahigh,dhigh,dlow
            if(cv(iv) .ge. 2)then
                if(abs(anew - a(iv))/dstep(iv) .gt. .01)then
		    write(6,*)' going to mode 8 (26511) anew/a'
     1                       ,anew,a(iv)  
                    goto26511
                else
                    write(6,*)' mode 8 bypassed'
                    goto26508
                endif
            endif
        endif

22900   if(dhigh .lt. dlow)then
            idir_opt = 1
        else
            idir_opt = -1
        endif
 
        crv(iv) = dhigh+dlow
        if(abs(a(iv)-mhigh(iv)) .lt. 1d-16 .and. idir_opt .eq. 1 .or.
     1     abs(a(iv)-mlow(iv))  .lt. 1d-16 .and. idir_opt .eq. 1
     1                                                      )goto26710      
        if(dhigh .gt. 0. .and. dlow .gt. 0.)then
            if(cv(iv) .eq. 0)then
                goto26508
            else
                goto25660
25660           a(iv) = anew
                modeop = 8
                return
            endif
        endif

        if(cv(iv) .ge. 2)then
!           Set variables for convergence speedup
!                      cv(iv) = 2
            a(iv) = anew
            g5 = 1
            goto26538
        endif

        goto23300

23270   g5 = 0
        modeop = 2
        return                                  

23300   if(iswp .eq. 1)then
           dstep_swp = 10.*dstep(iv)*idir_opt
        else
           dstep_swp =  4.*dstep(iv)*idir_opt
        endif

        f_merit = amid
23500   f9 = f_merit
        a(iv) = a(iv) + dstep_swp
        modeop = 5
        return            

!****** Modeop = 5 ***********************************************************
23800	if(f_merit .ge. f9)then
            goto24100
        else
	    g5=2
            a(iv) = a(iv) - dstep_swp
            goto26538
        endif

23900   g5=0
        a(iv) = a(iv) + dstep_swp
        dstep_swp = dstep_swp*1.5
        goto23500

!       cv(iv) = 0 or 1
24100   aahigh = a(iv)
        aamid = aahigh - dstep_swp
        aalow = aamid - dstep_swp
        amid = f9
        goto25200

24600   dstep_swp = dstep_swp * .5
        aahigh = aamid + dstep_swp
        aalow = aamid - dstep_swp
        a(iv) = aahigh
        modeop = 6
        return

 !****** Modeop = 6 ***********************************************************
25200	ahigh = f_merit
        a(iv) = aalow
        modeop = 7
        return

 !****** Modeop = 7 **********************************************************
25600	alow = f_merit
!       Convergence Speedup
        if(cv(iv) .ne. 0)then
            if(a9 .eq. 0.)a8 = anew
            if(ahigh - 2.*amid + alow .ne. 0.)then
                anew = aamid
     1               + dstep_swp*(.5+(amid-ahigh)/(ahigh-2.*amid+alow))       
            else
                anew = 1d10

            endif

            if(l_extra_output .eqv. .true.)then
                write(6,*)' iq,as,da,alow,amid,ahigh,ax',
     1                    itervar,a9,dstep_swp,alow,amid,ahigh,a8
            endif

!           Test for convergence cv(iv) = 1
            if(abs(anew-a8) .lt. dstep(iv))then 
                a9 = a9+1
            else
                a9 = 0
            endif

            if(abs(anew-aamid) .gt. abs(dstep_swp))a9=0

            if(a9 .ge. 2)then
                crv(iv) = (ahigh - 2.*amid + alow) *
     1                    (dstep(iv)/dstep_swp)**2
                a(iv) = anew ! formerly 25660
                modeop = 8
                return
            endif

        endif

        if(alow .lt. amid)then
            m9 = alow
        else
            m9 = amid
        endif

        if(ahigh .lt. m9)m9 = ahigh
        if(amid .ne. m9)then
            if(ahigh .eq. m9)then
                aamid = aahigh
                amid = ahigh
            else
                aamid = aalow
                amid = alow
            endif
        endif

!
!       Test for convergence cv(iv) = 0
        if(abs(dstep_swp) .gt. dstep(iv))goto24600
        crv(iv)=(ahigh - 2.*amid + alow)*(dstep(iv)/dstep_swp)**2
        a(iv) = aamid
        f_merit = amid

!****** Modeop = 8 ***********************************************************
26508	if(f_merit .le. f_ref)then
            goto26538
        else
            goto26511
        endif

26509   if(cv(iv) .le. 1)then
            write(6,*)' Too many iterations on one variable'
            dstep(iv) = dstep(iv) * 5.
            write(6,*)' Increasing step size of variable ',iv,
     1                ' to',dstep(iv)
            goto 26520
        endif

26511   if(abs(a(iv)-a_ref(iv))/dstep(iv) .ge. .01 .or. nv .eq. 1
     1                                                            )then       
            cv(iv) = cv(iv) - 1
            write(6,*)' Switching to convergence mode',cv(iv)
            write(6,*)' Modeop = ',modeop
        else
	    g8 = 1
	    write(6,*)' a is close to a_ref',a(iv),a_ref(iv)
            write(6,*)' Sweep',iswp,'  on variable',iv,
     1                '  aborted cv=',cv(iv)
        endif

26520	write(6,*)' Sweep,f,ft',iswp,f_merit,f_ref
        write(6,*)' f-ft,iv,iq',f_merit-f_ref,iv,itervar
        write(6,*)' a(iv),at(iv),a-at',a(iv),a_ref(iv),a(iv)-a_ref(iv)

        a(iv) = a_ref(iv)
        f_merit = f_ref
        if(g8 .ne. 0)then
            g8 = 0
            goto26538
        else
            modeop = 2
            itervar = 1
            goto21600
        endif

!       Adjust An if out of bounds
26538   if(mhigh(iv) .lt. mlow(iv))then
            modeop = 12
            write(6,*)' Warning, Unreal System, Vanished Constraints'
            return
        endif

        if(a(iv) .gt. mhigh(iv))then
            a(iv) = mhigh(iv)
            modeop = 8
            return
        endif

        if(a(iv) .lt. mlow(iv))then
            a(iv) = mlow(iv)
            modeop = 8
            return
        endif

        if(g5 .eq. 1)then
            goto23270
        else if(g5 .eq. 2)then
            goto23900
        endif

26710   if(l_extra1 .eqv. .true. .and. iterop .le. 2500)then
            write(6,105)iterop,iv,f_merit
105         format(15x,'Iterop = ',i5,'      iv',i4,'     f',d18.9)
!           if(iswp .le. 2)then
!               write(6,101)a(j)
!           endif            
        endif

!       Advance to next variable
        iv = iv + 1
        itervar = 1
        iterswp = 1
        if(iv .le. nv)goto21600

        if(f_merit-g_ref .ge. 1d-14)then
            k_ = .false.
            k8 = 5
            l_aitken_thistime = .false.
            write(6,*)
     1             ' Aitken failure override is cancelled mp,f,g_ref='
     1              ,modeop,f_merit,g_ref
            do j = 1,nv
                a(j) = g1(j)
            enddo
            f_merit = g_ref
            iv = 1
            kswp = k_ref
            goto27720
        endif

27150   iswp = iswp + 1
        iv = 1
        kswp = kswp + 1
        lswp = lswp + 1
        g_ref = f_merit
        if(nv .le. 1)then
            modeop = 10
            write(6,*)
            return
        endif

!       End of Sweep, Reset convergence variables every fourth sweep
        if(kswp * .25 .eq. int(kswp * .25) .or. kswp .eq. 2)then
            do j = 1,nv
                cv(j) = cv_ref(j)
                c1(j) = 1
            enddo
        endif

        modeop = 10
        l_aitken_thistime = .true.
        c9_outstring = ' A(n)s'
        if(l_extra1 .eqv. .true.)write(6,101)c9_outstring,(a(j),j=1,nv)       
101     format(1x,a9,8d14.5/8d14.5)
27508   k_ref = 1
        if(l_swp_extrpl .eqv. .false.)then
            do j = 1,nv
                g1(j) = e1(j)
                e1(j) = d1(j)
                d1(j) = a(j) - a_ref(j)
                if(d1(j) .eq. 0.)c1(j) = 0.
!               write(6,*)' j,e1(j)',j,e1(j)

!               Decide whether to apply Aitken Acceleration
                if(d1(j) .ne. 0. .and. abs(e1(j)) .gt. 1d-14)then
 
 	            if(abs(1.-d1(j)*g1(j)/(e1(j)**2)) .gt. .06
     1          .or.   abs(d1(j)/e1(j))               .gt. .97 )then

				       l_aitken_thistime = .false.
                    endif

                else
			           l_aitken_thistime = .false.
                endif

            enddo
    
!	    Adjust Step Sizes
            if(f_merit .gt. cfive*c7)then
                m9 = f_merit/c7
            else
                m9 = cfive
            endif

            do j = 1,nv
                qq = 0
                if(crv(j) .ne. 0. .and. l_depth)then
                    j9 = .5 * log(m9/crv(j))
                    if(abs(j9) .gt. log(2.) .and. u7 .eq. 0)then
                        j9 = log(2.)*j9/abs(j9)
                        qq = 1
                    endif
                else
                    j9 = 0
                endif

                dstep(j) = dstep(j) * exp(j9)
            enddo
	    u7 = 0
            if(l_extra1 .eqv. .true. .and. 
     1         iterop .le. 1000 .and. l_depth .eqv. .true.)then
                c9_outstring = ' Step(n)s'
                write(6,101)c9_outstring,(d1(j)/dstep(j),j=1,nv)
                c9_outstring = ' d1(n)s'
                write(6,101)c9_outstring,(d1(j),j=1,nv)
                c9_outstring = ' dn(n)s'
                write(6,101)c9_outstring,(dstep(j),j=1,nv)
            endif

!           Test for Convergence
!           Examine maximum change of all variables optimized in this sweep
            zsum = 0.
            do j = 1,nv
                zsum = zsum + ((a(j)-a_ref(j))/dstep(j))**2
            enddo
            if(zsum .gt. 1. .or. qq .eq. 1 .and. v5 .eq. 0)modeop = 2
            write(6,104)iswp-1,iterop,f_merit
104         format(' End of sweep ',i5,'      ip= ',i5,'      f ',d17.7)
            write(6,*)

        else
            goto27725

        endif

!       Apply Aitken acceleration or Sweep extrapolation
27720   if(modeop .gt. 9)return

        if(l_aitken .and. l_aitken_thistime)then
            write(6,*)' Aitken Acceleration'
            goto 27730
        endif

!       write(6,*)' Testing - kswp,swps_per_extrpl'
!	1			,kswp,swps_per_extrpl
27724   if(kswp .lt. swps_per_extrpl)goto21600

!       Initialize Sweep Extrapolation
        l_swp_extrpl = .true.
        l6 = 1
        l7 = 0
        bottom = 0.
        do j = 1,nv
                h1(j) = a(j)
        enddo ! j

!       Loop through Sweep Extrapolation
27725   bottom_last = bottom
	if(l7 .gt. 2.)then
		slope1 = (f_merit - f_m1) / (l7 - l_m1)
		slope2 = (f_m1 - f_m2) / (l_m1 - l_m2)
                crv_swp = (slope1 - slope2) / ((l7 - l_m2) / 2.)
                bottom = (l7 + l_m1)/2. - slope1 / crv_swp
        else
                crv_swp = 0.
                bottom = 0.
        endif


	if(l_extra1 .eqv. .true. .and. iterop .le. 1000)then
!   	        c9_outstring = ' A(n)s'
!               write(6,101)c9_outstring,(a(j),j=1,nv)
                write(6,107)iterop,f_merit,l7,crv_swp,bottom
107             format(' Sweep extpl ip,f,l7,crv,btm='
     1                      ,i5,d15.6,f10.2,d11.3,f10.2)
        endif

27730	k_ref = 0
        kswp = 1
        do j = 1,nv
            g1(j) = a(j)
            a_ref(j) = a(j)
            if(e1(j) .ne. 0)then
                if(l_swp_extrpl .eqv. .true.)then
                    a(j) = a(j)+d1(j)*l6
                else
                    a(j) = a(j)-d1(j)-e1(j) - e1(j)**2/(d1(j)-e1(j))
                endif
            endif
        enddo
        f_ref = f_merit
        l_m2 = l_m1
        f_m2 = f_m1
        l_m1 = l7
        f_m1 = f_merit
        l7 = l7 + l6
        l6 = l6 * swp_ratio
        modeop = 9
        return

!****** Modeop = 9 ***********************************************************
27745	itervar = 1
        iterswp = 1

        if(l_swp_extpl_end .eqv. .true.)then
  	    c9_outstring = ' A(n)s'
            write(6,101)c9_outstring,(a(j),j=1,nv)
            write(6,204)f_merit
204         format('  The bottom value of f is ',d15.6,
     1             ' Setting Modeop to 2'/)
            f_ref = f_merit
            l_swp_extpl_end = .false.
            modeop = 2
            goto21600
        endif

!       Check bounds
        do j = 1,nv
            if(a(j) .gt. mhigh(j) .or. a(j) .lt. mlow(j))then
                if(l_swp_extrpl .eqv. .true. .or. mlow(j) .gt. mhigh(j))     
     1                                                         goto27750       
                if(a(j) .gt. mhigh(j))then
                    a(j) = mhigh(j)
                else
                    a(j) = mlow(j)
                endif
            endif
        enddo

        if(l_swp_extrpl)then
            if(abs(bottom-bottom_last) .lt. 1. .and. bottom .gt. 0.)then
                goto27750
            endif
        endif

        if(f_merit .le. f_ref)then
            if(l_swp_extrpl .eqv. .true.)then
                k_ref = 1
                goto27725
            else
                modeop = 2
                write(6,106)f_merit
106             format(' After Aitken Acceleration, f_merit = ',d17.7)
                u7 = 1
                goto21600
            endif
        endif

27750   if(l_swp_extrpl)then
            write(6,*)' Sweep extrapolation'
        else
            write(6,*)' Aitken Acceleration'
        endif

        write(6,203)f_merit,f_ref,f_merit-f_ref
203     format('  Aborted - ff,ft,ff-ft',3d15.5)

        if(.not.(l_swp_extrpl))then
            k8 = k8-1
            if(k8 .eq. 0)k_ = .true.
            if(k_ .and. mlow(j) .le. mhigh(j))goto27790

            do j = 1,nv
                a(j) = a_ref(j)
            enddo
            f_merit = f_ref
            kswp = k_ref
            modeop = 2
            goto27724
        else
!  	    c9_outstring = ' A(n)s'
!           write(6,101)c9_outstring,(a(j),j=1,nv)
!           write(6,*)


        endif

!       Modify number of sweeps per extrapolation
        if(k3 .ne. 0)then
            conv_rate_m2 = conv_rate_m1
            conv_rate_m1 = conv_rate
        endif

        conv_rate = log(gone/f_ref)/(iterop-g2)

        if(k3 .ne. -1)then
            if(conv_rate_m1 .ne. 0)then
                if(conv_rate_m2 .eq. 0)then
                    if(conv_rate .lt. conv_rate_m1)h9 = -h9
                else
                    if(conv_rate/conv_rate_m1 .lt. 
     1                 conv_rate_m1/conv_rate_m2)h9 = -h9
                endif
            endif
            swps_per_extrpl = swps_per_extrpl + h9
            gone = f_ref
            g2 = iterop
            if(k3 .ge. 5)then
                if(swps_per_extrpl .ge. k4)then
                    k3 = 0
                else
                    k3 = -2
                    k5 = 2
                    swps_per_extrpl = k4
                endif
            endif
        endif

        k3 = k3 + 1
        if(k3 .eq. 0)swps_per_extrpl = k5
        swps_per_extrpl = max(swps_per_extrpl,1)
        write(6,201)conv_rate,swps_per_extrpl,iterop,l7-l6/swp_ratio
201     format('  Conv rate ',d10.3,i4,' sweeps per extrpolation'
     1        ,i7,f9.2)

        if(bottom .eq. 0.)bottom = l_m1
        write(6,205)bottom
205     format('  Resetting a(n)s, Bottom = ',f10.2)
        do j = 1,nv
            a(j) = h1(j)+d1(j)*bottom
        enddo

!       c9_outstring = ' A(n)s'
!       write(6,101)c9_outstring,(a(j),j=1,nv)
!       write(6,*)
!       f_merit = f_ref
        kswp = 1
        l_swp_extrpl = .false.
        l_swp_extpl_end = .true.
        return

27790   write(6,*)' Abort Mode overridden - going for it'
        do j = 1,nv
            a_ref(j) = a(j)
        enddo
        f_ref = f_merit
        modeop = 2
        goto21600             

!****** Modeop = 10 **********************************************************
27800	if(v5 .le. 0)then
            modeop = 11
            return
        endif
        c7 = c7 * increment_ratio**2
        cfive = cfive / increment_ratio**2
        write(6,*)' Incrementing step sizes,c7,crv',c7,cfive
        do j = 1,nv
            dstep(j) = dstep(j)/increment_ratio
        enddo
        v5 = v5-1
        modeop = 2
        itervar = 1
        goto27720

	end

