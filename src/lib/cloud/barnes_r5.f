cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS
cdis
cdis    This software and its documentation are in the public domain and
cdis    are furnished "as is."  The United States government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  They assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  If significant modifications or enhancements
cdis    are made to this software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis
cdis
cdis
cdis
cdis
cdis
cdis

      subroutine barnes_r5(t,imax,jmax,kmax,to,wt_p,cf_modelfg
     1  ,l_perimeter,cld_snd,wt_snd,r_missing_data
     1  ,grid_spacing_m,i_snd,j_snd,n_cld_snd,max_cld_snd,istatus
     1  ,NX_DIM_LUT,NY_DIM_LUT,IX_LOW,IX_HIGH,IY_LOW,IY_HIGH,n_fnorm)

!     1997 Aug 01  K. Dritz  - Added NX_DIM_LUT, NY_DIM_LUT, IX_LOW,
!                              IX_HIGH, IY_LOW, IY_HIGH, and n_fnorm as
!                              dummy arguments.
!     1997 Aug 01  K. Dritz  - Removed PARAMETER statements for the above.
!     1997 Aug 01  K. Dritz  - Changed NX_L_MAX to imax and NY_L_MAX to jmax.
!     1997 Aug 01  K. Dritz  - Added r_missing_data as dummy argument.
!     1997 Aug 01  K. Dritz  - Removed include of lapsparms.for.

      include 'laps_cloud.inc'

      integer NX_DIM_LUT,NY_DIM_LUT,IX_LOW,IX_HIGH,IY_LOW,IY_HIGH

      integer*4 iskip,n_fnorm
!     parameter (iskip = 2)

      integer*4 lowi_lut(imax)
      integer*4 lowj_lut(jmax)

      integer*4 max_obs
      parameter (max_obs = 9049)

      dimension to(imax,jmax,kmax),t(imax,jmax,kmax),fnorm(n_fnorm)
     1  ,iob(max_obs),job(max_obs),kob(max_obs),nob(max_obs)
     1                          ,cf_modelfg(imax,jmax,kmax)
      dimension wt_p(imax,jmax,kmax)

      logical l_perimeter
      real*4 cld_snd(max_cld_snd,kmax)
      real*4 wt_snd(max_cld_snd,kmax)
      integer*4 i_snd(max_cld_snd)
      integer*4 j_snd(max_cld_snd)


      real iiilut(-NX_DIM_LUT:NX_DIM_LUT,-NY_DIM_LUT:NY_DIM_LUT)
      integer nlast(KCLOUD)
      logical l_analyze(KCLOUD)

      write(6,1234)
 1234 format(1x,'barnes_r5 called')

!     Set weight for using model background clouds beyond a certain effective
!     radius of influence from the sfc obs/pireps
      weight_modelfg = 0.    ! Model wt not active, obs used to infinite radius
!     weight_modelfg = 1.    ! Model used beyond ~100km from nearest obs
!     weight_modelfg = .01   ! Model used beyond ~200km from nearest obs
!     weight_modelfg = .0001 ! Model used beyond ~400km from nearest obs

!     Obtain and/or iterate for value of iskip
      iskip = nint(20000. / grid_spacing_m)
      iskip = max(iskip,1)

100   rden_ratio = ((float(imax)-1.)/float(iskip))
      if(abs(rden_ratio - nint(rden_ratio))  .gt. .001)then
          write(6,*)' Bad value of iskip'
          if(iskip .gt. 1)then
              iskip = iskip - 1
              goto 100
          elseif(iskip .eq. 1)then
              write(6,*)' Code error - stop'
              stop
          endif
      else
          write(6,*)' Good value of iskip = ',iskip
      endif

      if(l_perimeter)write(6,*)' Number of cloud soundings = ',n_cld_snd

      istatus = 1
      exponent_distance_wt = 5.0

      iiizero = 0

      do iii = 1,n_fnorm
        bias_iii = 1e0
        fnorm(iii) = (100.*bias_iii/float(iii)) ** (exponent_distance_wt
     1 / 2.0)
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
              if(i_snd(n) .lt. IX_LOW .or. i_snd(n) .gt. IX_HIGH) go to 
     1233
              if(j_snd(n) .lt. IY_LOW .or. i_snd(n) .gt. IY_HIGH) go to 
     1233

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

      write(6,*)' radm2*100/iii_max/n_fnorm = ',radm2*100.,iii_max,n_fno
     1rm

      if(iii_max .gt. n_fnorm)then
          write(6,*)' iii_max is too large, increase n_fnorm',iii_max,n_
     1fnorm
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
50        format(' lvl,nstart,nstop,nobs=',4i5)

!         height_level = height_of_level(k)

!         Analyze every other grid point
          do j=1,jmax,iskip
          do i=1,imax,iskip
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

            sum   = sum   + weight_modelfg * cf_modelfg(i,j,k)
            sumwt = sumwt + weight_modelfg

            if (sumwt.eq.0.)then
              t(i,j,k) = r_missing_data
              istatus = 0
            ELSE
              t(i,j,k)=sum/sumwt
            end if

          enddo ! i
          enddo ! j

!         Bilinearly interpolate to fill in rest of domain
          do i = 1,imax
              lowi_lut(i) = (i-1)/iskip*iskip + 1
              if(i .eq. imax)lowi_lut(i) = lowi_lut(i) - iskip
          enddo ! i
          do j = 1,jmax
              lowj_lut(j) = (j-1)/iskip*iskip + 1
              if(j .eq. jmax)lowj_lut(j) = lowj_lut(j) - iskip
          enddo ! i

          do j=1,jmax
              jl = lowj_lut(j)
              jh = jl + iskip
              fracj = float(j-jl)/float(iskip)

              do i=1,imax
                  il = lowi_lut(i)
                  ih = il + iskip
                  fraci = float(i-il)/float(iskip)

                  Z1=t(il,jl,k)
                  Z2=t(ih,jl,k)
                  Z3=t(ih,jh,k)
                  Z4=t(il,jh,k)

                  t(i,j,k) =  Z1+(Z2-Z1)*fraci+(Z4-Z1)*fracj
     1                - (Z2+Z4-Z3-Z1)*fraci*fracj

              enddo ! i
          enddo ! j

        elseif(nobs .gt. 0)then ! Obs are identical; Cpy analysis fm 1 lvl below
          write(6,51)k,nstart,nstop,nobs
51        format(' lvl,nstart,nstop,nobs=',4i5,
     1                  ' Identical Obs; Copy from 1 lvl down')

          km1 = k - 1

          do j=1,jmax
          do i=1,imax
              t(i,j,k) = t(i,j,km1)
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

