cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps
cdis
cdis    this software and its documentation are in the public domain and
cdis    are furnished "as is."  the united states government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  they assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  if significant modifications or enhancements
cdis    are made to this software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis
cdis
cdis
cdis
cdis
cdis
cdis
        real function func (x)

c       $log: func.for,v $
c revision 1.1  1996/08/30  20:44:48  birk
c initial revision
c

c       basically a simple cost function for the moisture function giving
c       goes-8 radiances for the moisture channels only  image case,
c       the passed parameter is the variable changed in the iteration
c       other variable (invariant for this problem) are passed via common
c       note: w_cost is same as w, but w has to be in a common for the
c       forwardmodel call reqd by wisconsin code.

c       the function blows up whenever w_cost is negative. therefore we test
c       to see if any w_cost element is negative if it is we assign a high
c       cost and return.

        implicit none
        save
        real x(3)
        integer lsfc
        real psfc
        real emiss
        real gimrad
        integer ngoes,kan(18)
        integer ngoes_cost,isnd_cost,chan_used
        real radiance (18)
        real theta,tau(40),tskin

        integer lsfc_cost
        real radiance_ob(18),p_cost(40),t_cost(40),
     1  ozo_cost(40),tskin_cost,psfc_cost,
     1  theta_cost,w_cost(40)
        real p(40),t(40),w(40),ozo(40)
        common/atmos/p,t,w,ozo
        common/cost_var/radiance_ob, p_cost, t_cost,ozo_cost,
     1  tskin_cost, lsfc_cost,psfc_cost,theta_cost,w_cost,
     1  ngoes_cost,isnd_cost,chan_used

        integer i,j


        ngoes = ngoes_cost

c  set up for sounder if needed instead of imager

        if (isnd_cost.eq.1) then ! sounder radiances used



            kan(1) = 10 ! 7.4 like ch3
            kan(2) = 8  ! 11.02 like ch4
            kan(3) = 7  ! 12.02 like ch5
            kan(4) = 11
            kan(5) = 16
            kan(6) = 6
            kan(7) = 13

          else  ! imager radiances used


            kan(1) = 23 !  7.4  imager ch3
            kan(2) = 24  ! 11.02  imager ch4
            kan(3) = 25  ! 12.02  imager ch5
            

        endif

 
        emiss = .99


c       set up wisconsin variables in common

        do i = 1,40

        if(i.le.40 .and. i.gt. 35) then  ! sfc to 780
        w(i) = abs(x(1)) * w_cost(i)
        elseif (i.le.35 .and. i.gt. 30) then  ! 700 to 500
        w(i) = abs(x(2)) * w_cost(i)
        elseif (i.gt. 20) then   ! between 475 and 100
        w(i) = abs(x(3)) * w_cost(i)
        else
        w(i) =  w_cost(i)
        endif

        tskin = tskin_cost


        t(i) = t_cost(i)
        p(i) = p_cost(i)
        ozo(i) = ozo_cost(i)
        enddo
        lsfc = lsfc_cost
        psfc = psfc_cost
        theta = theta_cost


c       perform forward model computation for radiance


        do j=1,chan_used

        call taugim(t,w,ozo,theta,ngoes,kan(j),tau)
        radiance(j) = gimrad(tau,t,tskin,kan(j),lsfc,psfc,emiss)

        enddo ! j

c       compute cost function

        func = 0.0

        do j = 1,chan_used
        func = func + (radiance_ob(j)-radiance(j))**2/2.
        enddo
        do j = 1,3
        func = func + ((x(j) - 1.)**2 )
        enddo

c       print *, (radiance(i),i=1,3),x


        return
        end

