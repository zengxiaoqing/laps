cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis
cdis
cdis   
cdis

      subroutine barnes_r5(t,imax,jmax,kmax,to,wt_p,cf_modelfg
     1  ,l_perimeter,cld_snd,wt_snd,r_missing_data
     1  ,grid_spacing_m,i_snd,j_snd,n_cld_snd,max_cld_snd
     1  ,max_obs,weight_modelfg,NX_DIM_LUT,NY_DIM_LUT
     1  ,IX_LOW,IX_HIGH,IY_LOW,IY_HIGH,n_fnorm,istatus)

!     1997 Aug 01  K. Dritz  - Added NX_DIM_LUT, NY_DIM_LUT, IX_LOW,
!                              IX_HIGH, IY_LOW, IY_HIGH, and n_fnorm as
!                              dummy arguments.
!     1997 Aug 01  K. Dritz  - Removed PARAMETER statements for the above.
!     1997 Aug 01  K. Dritz  - Changed NX_L_MAX to imax and NY_L_MAX to jmax.
!     1997 Aug 01  K. Dritz  - Added r_missing_data as dummy argument.
!     1997 Aug 01  K. Dritz  - Removed include of lapsparms.for.

      include 'laps_cloud.inc'

      integer NX_DIM_LUT,NY_DIM_LUT,IX_LOW,IX_HIGH,IY_LOW,IY_HIGH

      integer*4 nskip,n_fnorm
!     parameter (nskip = 2)

      integer*4 lowi_lut(imax)
      integer*4 lowj_lut(jmax)

      dimension to(imax,jmax,kmax),t(imax,jmax,kmax),fnorm(n_fnorm)
     1  ,iob(max_obs),job(max_obs),kob(max_obs),nob(max_obs)
     1                          ,cf_modelfg(imax,jmax,kmax)
      dimension wt_p(imax,jmax,kmax)

      logical l_perimeter
      real*4 cld_snd(max_cld_snd,kmax)
      real*4 wt_snd(max_cld_snd,kmax)
      integer*4 i_snd(max_cld_snd)
      integer*4 j_snd(max_cld_snd)

      real*4 sum_a(imax,jmax)
      real*4 sumwt_a(imax,jmax)

      real iiilut(-NX_DIM_LUT:NX_DIM_LUT,-NY_DIM_LUT:NY_DIM_LUT)
      integer nlast(KCLOUD)
      logical l_analyze(KCLOUD)

      write(6,1234)
 1234 format(1x,'barnes_r5 called')

      if(n_cld_snd .gt. max_cld_snd)then
          write(6,*)' barnes_r5: ERROR, too many cloud soundings'
          istatus = 0
          return
      endif 

!     Note that this routine is considerably longer than it might need to
!     be to allow maximum efficiency. 'nskip' controls how many points
!     to skip in the i & j directions during the looping/analysis. The
!     remaining points are filled in later by a faster bilinear interpolation
!     step. Also, l_analyze(k) controls which levels are analyzed. If the
!     obs at a given level are a repeat of the obs at the next lower level,
!     the looping again is skipped and the relevant summations are copied
!     up from the lower level.

!     Obtain and/or iterate for value of nskip
      nskip = nint(20000. / grid_spacing_m)
      nskip = min(max(nskip,1),4)

100   rden_ratioi = ((float(imax)-1.)/float(nskip))
      rden_ratioj = ((float(jmax)-1.)/float(nskip))

      if(abs(rden_ratioi - nint(rden_ratioi))  .gt. .001
     1  .or.(rden_ratioj - nint(rden_ratioj))  .gt. .001)then
          write(6,*)' Bad value of nskip'
          if(nskip .gt. 1)then
              nskip = nskip - 1
              goto 100
          elseif(nskip .eq. 1)then
              write(6,*)' Code error - stop'
              stop
          endif
      else
          write(6,*)' Good value of nskip = ',nskip
      endif

      if(l_perimeter)write(6,*)' Number of cloud soundings = ',n_cld_snd

      istatus = 1
      exponent_distance_wt = 5.0

      iiizero = 0

      do iii = 1,n_fnorm
        bias_iii = 1e0
        fnorm(iii) = 
     1        (100.*bias_iii/float(iii)) ** (exponent_distance_wt / 2.0)       
        if(fnorm(iii) .eq. 0. .and. iiizero .eq. 0)then
            iiizero = iii
            write(6,*)' WARNING: fnorm array reached zero, iii=',iiizero
        endif
      enddo

      write(6,*)' Model background weight = ',weight_modelfg
      write(6,*)' Min possible ob weight = ',fnorm(n_fnorm)
      write(6,*)' Max possible ob weight = ',fnorm(1)

      ncnt=0
      do k=1,kmax

        if(.not. l_perimeter)then
            nlast(k) = ncnt
            do j=1,jmax
            do i=1,imax
              if(to(i,j,k) .eq. r_missing_data) go to 223
              ncnt=ncnt+1
              iob(ncnt)=i
              job(ncnt)=j
              kob(ncnt)=k
              nlast(k) = ncnt
!             if(k .eq. 8)write(6,*)i,j,k,to(i,j,k)
223           continue
            enddo ! i
            enddo ! j

        else ! l_perimeter
            nlast(k) = ncnt
            do n = 1,n_cld_snd
              if(cld_snd(n,k) .eq. r_missing_data) go to 233

!             Test if out of bounds of established perimeter around LAPS domain
              if(i_snd(n) .lt. IX_LOW .or. 
     1           i_snd(n) .gt. IX_HIGH) go to 233
              if(j_snd(n) .lt. IY_LOW .or. 
     1           j_snd(n) .gt. IY_HIGH) go to 233

              ncnt=ncnt+1
              iob(ncnt)=i_snd(n)
              job(ncnt)=j_snd(n)
              kob(ncnt)=k
              nob(ncnt)=n
              nlast(k) = ncnt
233           continue
            enddo ! n

        endif ! l_perimeter

        if(k .gt. 1)then
            km1 = k - 1
            l_analyze(k) = .false.

            if(.not. l_perimeter)then
              do j=1,jmax
              do i=1,imax
                if(to(i,j,k) .ne. to(i,j,km1))then
                    l_analyze(k) = .true.
                    goto250
                endif
              enddo ! i
              enddo ! j

            else
              do n=1,n_cld_snd
                if(cld_snd(n,k) .ne. cld_snd(n,km1))then
                    l_analyze(k) = .true.
                    goto250
                endif
              enddo ! n

            endif ! l_perimeter


250         continue

        else ! k .eq. 1
            l_analyze(k) = .true.
        endif

      enddo ! k

      if(ncnt.eq.0) then
         write(6,1002)
 1002    format(1x,'no data for barnes')
         return
      else
         write(6,*)' Ncnt = ',ncnt
      endif

      spcng = 10.
      radm2=1.533/spcng**2

      rr_max = (NX_DIM_LUT-1)**2 + (NY_DIM_LUT-1)**2
      iii_max = radm2 * 100. * rr_max + 1.

      write(6,*)' radm2*100/iii_max/n_fnorm = '
     1           ,radm2*100.,iii_max,n_fnorm

      if(iii_max .gt. n_fnorm)then
          write(6,*)' iii_max is too large, increase n_fnorm'
     1             ,iii_max,n_fnorm
          istatus = 0
          return
      endif

!     Create a lookup table for fnorm(iii)
      do i = -NX_DIM_LUT,NX_DIM_LUT
      do j = -NY_DIM_LUT,NY_DIM_LUT
              rr=i*i+j*j
              iii=radm2*100.*rr+1.
              if(iii .gt. n_fnorm)iii=n_fnorm
              iiilut(i,j) = fnorm(iii)
      enddo
      enddo

      do k=1,kmax

        if(k .eq. 1)then
          nstart = 1
        else
          nstart = nlast(k-1) + 1
        endif
        nstop = nlast(k)

        nobs = nstop-nstart+1

        if((l_analyze(k) .and. nobs .ge. 1) .or. k .eq. 1)then

          write(6,50)k,nstart,nstop,nobs
50        format(' lvl,nstart,nstop,nobs=',4i6)

!         Analyze every few grid points
          do j=1,jmax,nskip
          do i=1,imax,nskip
            sum=0.
            sumwt=0.
            if(.not. l_perimeter)then
                do n=nstart,nstop
                  ii=iob(n)
                  jj=job(n)
                  weight = iiilut(i-ii,j-jj) * wt_p(ii,jj,k) ! Obs are being weighted
                  sum=weight*to(ii,jj,k)+sum
                  sumwt=sumwt+weight
                enddo ! n
            else
                do n=nstart,nstop
                  ii=iob(n)
                  jj=job(n)
                  nn=nob(n)
                  weight = iiilut(i-ii,j-jj) * wt_snd(nn,k) ! Obs are being weighted
                  sum=weight*cld_snd(nn,k)+sum
                  sumwt=sumwt+weight
                enddo ! n
            endif

!           Save ob weight sums for possible later reuse
            sum_a(i,j) = sum

!           Add in model first guess as an ob
            sum   = sum   + weight_modelfg * cf_modelfg(i,j,k)
            sumwt = sumwt + weight_modelfg

!           Save ob weight sums for possible later reuse
            sumwt_a(i,j) = sumwt

!           Divide weights to get analysis = f(obs + background)
            if (sumwt.eq.0.)then
              t(i,j,k) = r_missing_data
              istatus = 0
            ELSE
              t(i,j,k)=sum/sumwt
            end if

          enddo ! i
          enddo ! j

!         Bilinearly interpolate to fill in rest of domain
!         Fills in final analysis value and weights from obs alone
!         We may have to extrapolate at the N and E edges
          do i = 1,imax
              lowi_lut(i) = (i-1)/nskip*nskip + 1
              il = lowi_lut(i)
              ih = il + nskip
              if(ih .gt. imax)lowi_lut(i) = lowi_lut(i) - nskip
          enddo ! i
          do j = 1,jmax
              lowj_lut(j) = (j-1)/nskip*nskip + 1
              jl = lowj_lut(j)
              jh = jl + nskip
              if(jh .gt. jmax)lowj_lut(j) = lowj_lut(j) - nskip
          enddo ! i

          do j=1,jmax
              jl = lowj_lut(j)
              jh = jl + nskip
              fracj = float(j-jl)/float(nskip)

              do i=1,imax
                  il = lowi_lut(i)
                  ih = il + nskip
                  fraci = float(i-il)/float(nskip)

!                 Calculate interpolated cloud cover
                  Z1=t(il,jl,k)
                  Z2=t(ih,jl,k)
                  Z3=t(ih,jh,k)
                  Z4=t(il,jh,k)

                  t(i,j,k) =  Z1+(Z2-Z1)*fraci+(Z4-Z1)*fracj
     1                - (Z2+Z4-Z3-Z1)*fraci*fracj

!                 Calculate interpolated ob summation
                  Z1=sum_a(il,jl)
                  Z2=sum_a(ih,jl)
                  Z3=sum_a(ih,jh)
                  Z4=sum_a(il,jh)

                  sum_a(i,j) =  Z1+(Z2-Z1)*fraci+(Z4-Z1)*fracj
     1                        - (Z2+Z4-Z3-Z1)*fraci*fracj

!                 Calculate interpolated ob weight summation
                  Z1=sumwt_a(il,jl)
                  Z2=sumwt_a(ih,jl)
                  Z3=sumwt_a(ih,jh)
                  Z4=sumwt_a(il,jh)

                  sumwt_a(i,j) =  Z1+(Z2-Z1)*fraci+(Z4-Z1)*fracj
     1                        - (Z2+Z4-Z3-Z1)*fraci*fracj

              enddo ! i
          enddo ! j

        elseif(nobs .gt. 0)then ! Obs are identical to lvl below; 
                                ! Use analysis weights from last analyzed level
          write(6,51)k,nstart,nstop,nobs
51        format(' lvl,nstart,nstop,nobs=',4i5,
     1           ' Identical Obs; Copy wts from last analyzed lvl')       

          km1 = k - 1

          do j=1,jmax
          do i=1,imax

!             Recover weight summations from last analyzed level below
              sum = sum_a(i,j)
              sumwt = sumwt_a(i,j)

!             Add in model first guess as an ob to 'sum'
!             Note that sumwt does not need to be modified

              sum   = sum   + weight_modelfg * cf_modelfg(i,j,k)

!             Divide weights to get analysis = f(obs + background)
              if (sumwt.eq.0.)then
                t(i,j,k) = r_missing_data
                istatus = 0
              ELSE
                t(i,j,k)=sum/sumwt
              end if

          enddo ! i
          enddo ! j

        else ! No Obs; Set level to model first guess
          write(6,52)k,nstart,nstop,nobs
52        format(' lvl,nstart,nstop,nobs=',4i5,
     1                  ' No Obs; Set level to model fg')

          do j=1,jmax
          do i=1,imax
              if(weight_modelfg .gt. 0.)then
                  t(i,j,k) = cf_modelfg(i,j,k)
              else
                  t(i,j,k) = cf_modelfg(i,j,k)
!                 t(i,j,k) = .01
              endif
          enddo ! i
          enddo ! j

        endif

      enddo ! k

      if(istatus .eq. 0)then
          write(6,*)' WARNING: Distant obs given zero weight in barnes_r
     15'
          write(6,*)' Try increasing bias_iii'
      endif

      return
      end

