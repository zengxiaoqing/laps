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

      subroutine barnes_multivariate
     1                     (var                                   ! Output
     1                     ,n_var,max_obs,obs_barnes_in           ! Input
     1                     ,imax,jmax,kmax,grid_spacing_m         ! Inputs
     1                     ,rep_pres_intvl                        ! Input
     1                     ,wt_b                                  ! Input
     1                     ,i4_loop_total                         ! Output
     1                     ,weight_3d,fnorm_dum,n_fnorm_dum       ! Inputs
     1                     ,l_analyze_dum,l_not_struct,rms_thresh ! Input
     1                     ,r0_barnes_max_m,barnes_conv_rate      ! Input
     1                     ,topo,ldf,ni_stat,nj_stat              ! Input
     1                     ,n_obs_lvl,istatus)                    ! Outputs

      use mem_namelist, ONLY: l_superob_barnes

      integer max_obs ! also is the number of input full obs

!     ....................Declare I/O variables...............................

      integer imax,jmax,kmax                                        ! Input
      integer n_var                                                 ! Input

      include 'barnesob.inc'
      type (barnesob) :: obs_barnes_in(max_obs)   ! (input)           Input
      type (barnesob) :: obs_barnes_buff(max_obs) ! (buffer)          Local
      type (barnesob) :: obs_barnes(max_obs)      ! (superobs)        Local

      real    var(imax,jmax,kmax,n_var)                             ! Output
      real, allocatable, dimension(:,:,:,:) :: var_buff             ! Local

      real    grid_spacing_m                                        ! Input

!     Weight for Model Background. Recommended values: 0. to 1e+30.
!     This will make the output values match the background if far from obs.
!     A value of zero means this parameter is not active.
      real    wt_b(imax,jmax,kmax,n_var)                            ! Input?

!     Note that weight_3d is no longer being used
      real    weight_3d(imax,jmax,kmax)                             ! Input?

      integer n_fnorm_dum                                           ! Input
      real    fnorm_dum                                             ! Input
      logical   l_not_struct,l_struct                               ! Input

      integer n_obs_lvl(kmax)                                       ! Phaseout
      integer istatus                                               ! Output

!     Max Barnes radius (meters). Recommended values: 120000. to 300000.
!                                 Best value: 240000.
      real r0_barnes_max_m
      real fac,dfac

      logical ltest_vertical_grid,l_rms_exceeded

!     This static data is used only in the case of the sfc analysis
      real topo(ni_stat,nj_stat), ldf(ni_stat,nj_stat)              ! Input
      ialloc = 0

      l_struct = .not. l_not_struct

      write(6,*)' Subroutine barnes_multivariate..., l_superob_barnes='
     1                                              ,l_superob_barnes

!     vertical radius of influence in gridpoints
      r0_vert_grid = 0.9 * (5000. / rep_pres_intvl)

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Error getting r_missing_data'
          return
      endif

      call get_systime_i4(i4time_sys,istatus)

!     Find minimum value of wt_b array
      weight_bkg_min = r_missing_data
!     fac = 4.0E+28
!     dfac = 2.0E+28/(n_var*kmax*jmax*imax)
      do l = 1,n_var
      do k = 1,kmax
      do j = 1,jmax
      do i = 1,imax
 !        wt_b(i,j,k,l) = fac
 !        fac = fac + dfac
          weight_bkg_min = min(wt_b(i,j,k,l),weight_bkg_min) 
      enddo ! i
      enddo ! j
      enddo ! k
      enddo ! l

!     Tolerance denoting weight sum allowed for observations outside the
!     radius of influence that are omitted, divided by model background weight.
!     This represents the truncation error of the analysis increment.
      rtol = 1e-5 ! 1e-6

!     This is a buffer routine that passes variables into actual
!     barnes routine. This allows for things like dynamic memory allocation

      n_iter = 10

!     if(kmax .gt. 1)then ! 3D temp/wind
!         r0_barnes_max_m = 240000.
!         barnes_conv_rate = 0.8
!     else                ! sfc analysis
!         r0_barnes_max_m = 140000.
!         barnes_conv_rate = 0.5
!     endif

!     This loop does multiple passes at decreasing radii
      do iter = 1,n_iter

          i4_elapsed = ishow_timer()

          r0_barnes_max_in  = r0_barnes_max_m 
     1                      * barnes_conv_rate ** (iter-1)
          r0_vert_grid_in   = r0_vert_grid    
     1                      * barnes_conv_rate ** (iter-1)       

          if(iter .eq. 1)then

!             Determine super obs
              i_ratio_max = 0.2 * (r0_barnes_max_in / grid_spacing_m)
              i_ratio = max(i_ratio_max,1)
              if(l_superob_barnes .and. kmax .gt. 1 
     1                            .and. i_ratio .gt. 1)then
                  call get_barnes_superobs(imax,jmax,kmax,r_missing_data ! I 
     1                        ,n_var,i_var,max_obs,obs_barnes_in         ! I
     1                        ,i_ratio,i4time_sys                        ! I
     1                        ,obs_barnes,ncnt_superobs)                 ! O
                  ncnt_total = ncnt_superobs

              else
                  write(6,*)' Using all obs - skip superobs routine'
                  obs_barnes = obs_barnes_in
                  ncnt_total = max_obs

              endif

              if(.true.)then
!                 ncnt_total = max_obs
                  weight_total = 0.
                  do n = 1,ncnt_total
                      call get_time_wt(i4time_sys,obs_barnes(n)%i4time
     1                                ,time_wt,istatus)
                      weight_total =  weight_total 
     1                             + (obs_barnes(n)%weight * time_wt)
                  enddo ! n

              endif

              i4_elapsed = ishow_timer()
              i4_loop_start = i4_elapsed

              call barnes_multivariate_sub
     1                     (var                              ! Output
     1                     ,wt_b                             ! Input
     1                     ,n_var,max_obs,obs_barnes         ! Input
     1                     ,imax,jmax,kmax,grid_spacing_m    ! Inputs
     1                     ,iter,i4time_sys                  ! Input
     1                     ,weight_bkg_min                   ! Input
     1                     ,rtol,r_missing_data              ! Input
     1                     ,r0_barnes_max_in,r0_vert_grid_in ! Input
     1                     ,ncnt_total,weight_total          ! Input
     1                     ,topo,ldf,ni_stat,nj_stat         ! Input
     1                     ,istatus)                         ! Outputs

              i4_elapsed = ishow_timer()

!             calculate RMS
              i4_loop_end = i4_elapsed
              i4_loop_total = i4_loop_total + 
     1                       (i4_loop_end - i4_loop_start)
              if(i4_elapsed .gt. 0)then
                  pct_loop_total = 
     1                float(i4_loop_total)/float(i4_elapsed) * 100.
              else
                  pct_loop_total = 0.
              endif
              
              write(6,101)            i4_loop_total            
     1                               ,i4_elapsed - i4_loop_total
     1                               ,i4_loop_start
     1                               ,i4_loop_end  
     1                               ,i4_elapsed     
     1                               ,pct_loop_total
101           format(' i4time loop/nonloop/start/end/elapsed/% = ',
     1               5i5,f8.2)

          elseif(iter .gt. 1)then
              if(r0_barnes_max_in .lt. 1.7 * grid_spacing_m)then
                  write(6,*)
     1                ' Radius reached small criterion, exit iter loop '
     1                     ,iter,r0_barnes_max_in,grid_spacing_m
                  goto900
              endif

!             Determine superobs
              i_ratio_max = 0.2 * (r0_barnes_max_in / grid_spacing_m)
              i_ratio = max(i_ratio_max,1)
              if(l_superob_barnes .and. kmax .gt. 1 
     1                            .and. i_ratio .gt. 1)then
                  call get_barnes_superobs(imax,jmax,kmax,r_missing_data ! I 
     1                        ,n_var,i_var,max_obs,obs_barnes_in         ! I
     1                        ,i_ratio,i4time_sys                        ! I
     1                        ,obs_barnes,ncnt_superobs)                 ! O
                  ncnt_total = ncnt_superobs

              else
                  write(6,*)' Using all obs - skip superobs routine'
                  obs_barnes = obs_barnes_in
                  ncnt_total = max_obs

              endif

!             obs_buff = obs-anal            [if iter > 1]

              do n = 1,ncnt_total
                  i = obs_barnes(n)%i
                  j = obs_barnes(n)%j
                  k = obs_barnes(n)%k

!                 Equate the entire structure element initially
                  obs_barnes_buff(n) = obs_barnes(n)

!                 Subtract (obs - analysis) values
                  do l = 1,n_var
                      obs_barnes_buff(n)%value(l) = 
     1                obs_barnes(n)%value(l) - var(i,j,k,l)        
                  enddo ! l

              enddo ! n

              if(ialloc .eq. 0)then
                  allocate(var_buff(imax,jmax,kmax,n_var)
     1                    ,STAT=istat_alloc)
                  if(istat_alloc .ne. 0)then
                      write(6,*)' ERROR: Could not allocate var_buff'
                      stop
                  endif
                  ialloc = 1
              endif

!             anal_buff = analyze(obs_buff)  [if iter > 1]

              i4_elapsed = ishow_timer()
              i4_loop_start = i4_elapsed

              call barnes_multivariate_sub
     1                     (var_buff                          ! Output
     1                     ,wt_b                              ! Input
     1                     ,n_var,max_obs,obs_barnes_buff     ! Input
     1                     ,imax,jmax,kmax,grid_spacing_m     ! Inputs
     1                     ,iter,i4time_sys                   ! Input
     1                     ,weight_bkg_min                    ! Input
     1                     ,rtol,r_missing_data               ! Input
     1                     ,r0_barnes_max_in,r0_vert_grid_in  ! Input
     1                     ,ncnt_total,weight_total           ! Input
     1                     ,topo,ldf,ni_stat,nj_stat          ! Input
     1                     ,istatus)                          ! Output

!             anal = anal + anal_buff        [if iter > 1]
              do l = 1,n_var
              do k = 1,kmax
              do j = 1,jmax
              do i = 1,imax
                  var(i,j,k,l) = var(i,j,k,l) + var_buff(i,j,k,l)
              enddo ! i
              enddo ! j
              enddo ! k
              enddo ! l

              i4_elapsed = ishow_timer()

!             calculate RMS
              i4_loop_end = i4_elapsed
              i4_loop_total = i4_loop_total + 
     1                       (i4_loop_end - i4_loop_start)
              if(i4_elapsed .gt. 0)then
                  pct_loop_total = 
     1                float(i4_loop_total)/float(i4_elapsed) * 100.
              else
                  pct_loop_total = 0.
              endif
              write(6,101)            i4_loop_total             
     1                               ,i4_elapsed - i4_loop_total
     1                               ,i4_loop_start
     1                               ,i4_loop_end  
     1                               ,i4_elapsed     
     1                               ,pct_loop_total

              write(6,*)
              write(6,*)' Calculate RMS of buffer obs'

              do ivar = 1,n_var
                call get_rms_barnesobs(var_buff(1,1,1,ivar)              ! I
     1                ,imax,jmax,kmax,r_missing_data,nobs                ! I/O 
     1                ,n_var,ivar,max_obs,obs_barnes_buff,ncnt_total     ! I
     1                ,rms_obs,rms)                                      ! O

                write(6,1)iter,ivar,nobs,rms_obs,rms

              enddo ! ivar

          endif ! iter > 1

!         calculate RMS
          write(6,*)
          write(6,*)' Calculate RMS of full obs'

          l_rms_exceeded = .false.

          do ivar = 1,n_var
!             Use full (non-superobs) obs set?
              call get_rms_barnesobs(var(1,1,1,ivar)                 ! I
     1                ,imax,jmax,kmax,r_missing_data,nobs            ! I/O 
     1                ,n_var,ivar,max_obs,obs_barnes_in,max_obs      ! I
     1                ,rms_obs,rms)                                  ! O

              write(6,1)iter,ivar,nobs,rms_obs,rms
 1            format(' iter ivar nobs rms_obs rms: ',3i6,2f10.3)

              if(rms .gt. rms_thresh)l_rms_exceeded = .true.

          enddo ! ivar

          i4_elapsed = ishow_timer()

          if(.not. l_rms_exceeded)then
              write(6,*)' RMS criterion reached, exit iter loop '
     1                 ,rms_thresh
              goto900
          endif

      enddo ! iter

 900  continue

      if(ialloc .eq. 1)then
          deallocate(var_buff)
          ialloc = 0
      endif

      write(6,*)' Normal finish of barnes_multivariate.f'


      return
      end



      subroutine barnes_multivariate_sub
     1                     (var                             ! Output
     1                     ,wt_b                            ! Input
     1                     ,n_var,max_obs,obs_barnes        ! Input
     1                     ,imax,jmax,kmax,grid_spacing_m   ! Inputs
     1                     ,iter,i4time_sys                 ! Input
     1                     ,weight_bkg_const                ! Input
     1                     ,rtol,r_missing_data             ! Input
     1                     ,r0_barnes_max_m,r0_vert_grid    ! Input
     1                     ,ncnt_total,weight_total         ! Input
     1                     ,topo,ldf,ni_stat,nj_stat        ! Input
     1                     ,istatus)                        ! Outputs

!     Aug 1995 Steve Albers - This version handles multiple variables
!                             Used within wind, temp, and sfc analyses

!     ....................Declare I/O variables...............................

#ifdef USEMPI
      use mem_grid, ONLY: IstartIend,nPEs,rank,recvcounts,displs
#endif
      implicit none

      integer, INTENT(IN) :: n_var,max_obs                      

      include 'barnesob.inc'
      type (barnesob), INTENT(IN) :: obs_barnes(max_obs)

      integer, INTENT(IN) :: imax,jmax,kmax                         

      real, INTENT(OUT) :: var(imax,jmax,kmax,n_var)        ! Output

      real, INTENT(IN) :: wt_b(imax,jmax,kmax,n_var)  ! Input (for now)

      real, INTENT(IN) :: grid_spacing_m        
      integer, INTENT(IN) :: i4time_sys,iter

!     Weight for Model Background. Recommended values: 0. to 1e+30.
!     This will make the output values match the background if far from obs.
!     A value of zero means this parameter is not active.
      real, INTENT(IN) :: weight_bkg_const 
      real, INTENT(IN) :: rtol,r_missing_data
      real, INTENT(IN) :: r0_barnes_max_m,r0_vert_grid
      integer, INTENT(IN) :: ncnt_total
      real, INTENT(IN) :: weight_total

!     This static data is used only in the case of the sfc analysis
      real, INTENT(IN) :: topo(ni_stat,nj_stat), ldf(ni_stat,nj_stat)
      integer, INTENT(IN)  :: ni_stat,nj_stat
      integer, INTENT(OUT) :: istatus                               

!     ....................Declare Local variables.............................

      integer l,i,ii,j,jj,k,n
      real r0_norm_sq,r0_value,ratio_norm_sq,r_subscript

#ifdef USEMPI
      integer ibuf(4*ncnt_total+10)
      real    rbuf((2+max_var)*ncnt_total+2)
      real buffer(imax*jmax*kmax*n_var)
      real wt_b_local(imax,jmax,kmax,n_var)  
      real, allocatable, dimension(:,:,:,:) :: local_var            ! Local
#endif
      real, allocatable, dimension(:,:,:,:) :: sum_var              ! Local
      real, allocatable, dimension(:,:,:) :: sumwt                  ! Local

      real fnorm_max,weight
      real r0_calc_min
      real varl(2)

!     These have a kmax dimension at least temporarily
      integer n_cross_in

      real    fnorm_lut(-(imax-1):(imax-1),-(jmax-1):(jmax-1))
      real    fnorm_lut_jj(jmax)

      real r0_norm
      real r0_value_min

!     ....................Satisfy Implicit None.............................

      real exp_offset,vertical_weight,exponent,weight_term
      real time_wt,weight_landfrac,weight_varl2,const,relax,diffl
      real weight_varl1,varl2,varl1,weight_i,weight_term_i,expm80
      real fnorm_0,exp_m_offset,r0_vert_ob
      real r0_barnes_max,fnorm_calc,wt_bkg_rel,wt_ob_dist,rdelt_m
      real rdelt_roi

      integer kk_lower_pair, kk_higher_pair,klow,jhigh,khigh,kk,kmin,kkh
      integer kdelt,jdelt,ilow,jlow
      integer ihigh,idelt,istat_alloc
      integer i4_elapsed, ishow_timer
#ifdef USEMPI
      integer flag,ierr,count
      integer :: offset,ioffset,joffset,koffset,noffset,IMAXlocal
      integer :: ip,r
      include 'mpif.h'
#endif

      parameter (exp_offset = 70.)

!     ...................Start Executable Statements..........................

      istatus = 0

      write(6,1234)iter
 1234 format(1x,'barnes_multivariate_sub:',' iter=',i3)

      call get_fnorm_max(imax,jmax,r0_norm,r0_value_min,fnorm_max)
      expm80 = exp(-80.)
      exp_m_offset = exp(-exp_offset)

      r0_norm_sq = r0_norm**2

!     Maximum radius of influence in grid points (for no data coverage)
      r0_barnes_max = r0_barnes_max_m / grid_spacing_m

!     Check for possibility of illegal value in fnorm LUT
      write(6,*)'r0_value_min,fnorm_max'
     1          ,r0_value_min,fnorm_max

      write(6,*)'r0_barnes_max,r0_barnes_max_m'
     1          ,r0_barnes_max,r0_barnes_max_m
      write(6,*)'r0_vert_grid',r0_vert_grid

      fnorm_0 = fnorm_calc(0.,r0_norm_sq,exp_offset,expm80)

      write(6,*)'weight_bkg_const,fnorm_0,ratio'
     1          ,weight_bkg_const,fnorm_0,weight_bkg_const/fnorm_0

!     Create a lookup table for fnorm_lut
      do i = -(imax-1),(imax-1)
      do j = -(jmax-1),(jmax-1)
          if(.true.)then      ! l_3d
!             Nominal radius_sq of influence divided by 
!             desired radius_sq of influence
              ratio_norm_sq = r0_norm_sq / r0_barnes_max**2
              r_subscript = float(i*i+j*j) * ratio_norm_sq

!             Horizontal Weight (Scaled by e^70 ~= 2.5e30)
              fnorm_lut(i,j) = fnorm_calc(r_subscript,r0_norm_sq,
     1                                    exp_offset,expm80)
          endif         
      enddo
      enddo

      write(6,*)'ncnt_total',ncnt_total,'weight_total',weight_total

!     Calculate search radius for gridpoints around each ob based on
!     error tolerance. This will be extra cautious assuming lots of obs
!     associated with 3D weighting.

!     Weight of model background relative to observation at zero distance
      wt_bkg_rel = weight_bkg_const / fnorm_0
      write(6,*)'wt_bkg_rel',wt_bkg_rel

!     Distance weight allowed for individual observations outside the radius
!     of influence that are omitted.
      wt_ob_dist = wt_bkg_rel * (rtol / max(weight_total,1.0)  )

      rdelt_roi = sqrt(-log(wt_ob_dist))
      rdelt_m = r0_barnes_max_m * rdelt_roi
      idelt = nint(rdelt_m / grid_spacing_m)
      jdelt = idelt
      kdelt = int(rdelt_roi * r0_vert_grid) + 1

      write(6,*)'     rtol     wt_ob   grid_spac rdelt_roi  rdelt_m  '
     1         ,' idelt  kdelt'
      write(6,405)rtol,wt_ob_dist,grid_spacing_m,rdelt_roi,rdelt_m,idelt
     1           ,kdelt
 405  format(1x,2e11.3,f10.0,f8.3,f10.0,i8,i5)

      write(6,*)' Summing weights...'
      write(6,*)'   n    i    j    k   kk  vert_wt  time_wt  value'

#ifdef USEMPI
      flag=1
      call mpi_bcast(flag,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if( ierr /= 0) then
        print*,'MPI_BCAST flag error',ierr,count,rank
        stop
      endif
!     Pack and send integer buffer
      ibuf(1)  = n_var
      ibuf(2)  = max_obs
      ibuf(3)  = imax
      ibuf(4)  = jmax
      ibuf(5)  = kmax
      ibuf(6)  = i4time_sys
      ibuf(7)  = ncnt_total
      ibuf(8)  = idelt
      ibuf(9)  = jdelt
      ibuf(10) = kdelt
      call mpi_bcast(ibuf,10,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if( ierr /= 0) then
        print*,'MPI_BCAST ibuf(10) error',ierr,rank
        stop
      endif
      do n=1,ncnt_total
        ibuf(4*(n-1)+1) = obs_barnes(n)%i
        ibuf(4*(n-1)+2) = obs_barnes(n)%j
        ibuf(4*(n-1)+3) = obs_barnes(n)%k
        ibuf(4*(n-1)+4) = obs_barnes(n)%i4time
      enddo
      count = 4*ncnt_total
      call mpi_bcast(ibuf,count,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if( ierr /= 0) then
        print*,'MPI_BCAST ibuf error',ierr,count,rank
        stop
      endif
!     Pack and send real buffer
      rbuf(1) = r_missing_data
      rbuf(2) = r0_vert_grid
      count   = 2
      do n=1,ncnt_total
        rbuf(count+1:count+max_var) = obs_barnes(n)%value(1:max_var)
        rbuf(count+1      +max_var) = obs_barnes(n)%vert_rad_rat 
        rbuf(count+2      +max_var) = obs_barnes(n)%weight
        count                       = count+2+max_var
      enddo
      call mpi_bcast(rbuf,count,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      if( ierr /= 0) then
        print*,'MPI_BCAST rbuf error',ierr,count,rank
        stop
      endif
!     Copy wt_b, with re-ordering, into buffer
      ioffset = 0
      do r=0,nPEs-1
        IMAXlocal = IstartIend(2,r)-IstartIend(1,r)+1
        noffset = 0
        do n=1,n_var
          koffset = 0
          do k=1,kmax
            joffset = 0
            do j=1,jmax
              do i=IstartIend(1,r),IstartIend(2,r)
                ip     = i-ioffset
                offset = joffset+koffset+noffset
                buffer(ip+displs(r)+offset) = wt_b(i,j,k,n)
              enddo
              joffset=joffset+IMAXlocal
            enddo !j
            koffset=koffset+IMAXlocal*jmax
          enddo   !k
          noffset=noffset+IMAXlocal*jmax*kmax
        enddo     !n
        ioffset=ioffset+IMAXlocal
      enddo       !r

      i4_elapsed = ishow_timer()
      write(6,*)' Start mpi_scatterv'

!Send arrays
      call mpi_scatterv(buffer,recvcounts,displs,MPI_REAL,wt_b_local,
     1               recvcounts(rank),MPI_REAL,0,MPI_COMM_WORLD,ierr)
      if( ierr /= 0) then
        print*,'MPI_SCATTERV error in barnes',ierr,rank,recvcounts(rank)
        stop
      endif
      count = (2*(imax-1)+1)*(2*(jmax-1)+1)
      call mpi_bcast(fnorm_lut,count,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      if( ierr /= 0) then
        print*,'MPI_BCAST fnorm_lut error in barnes',ierr,count,rank
        stop
      endif
      allocate(local_var(IstartIend(1,0):IstartIend(2,0),jmax,
     1                         kmax,n_var), STAT=istat_alloc )
      if(istat_alloc .ne. 0)then
        print*,'ERROR: Could not allocate local_var in barnes',
     1         istat_alloc,rank
        print*,IstartIend(1,rank),IstartIend(2,rank),jmax,kmax,n_var
        stop
      endif

      i4_elapsed = ishow_timer()
      write(6,*)' Call analyze_weights_mpi'
      call analyze_weights_mpi(
     1          n_var,max_obs,obs_barnes                                    ! I
     1         ,imax,jmax,kmax                                              ! I
     1         ,wt_b_local                                                  ! I
     1         ,i4time_sys                                                  ! I
     1         ,r_missing_data,r0_vert_grid                                 ! I
     1         ,ncnt_total,idelt,jdelt,kdelt                                ! I
     1         ,fnorm_lut                                                   ! I
     1         ,IstartIend(1,0),IstartIend(2,0)                             ! I
     1         ,local_var)                                                  ! O
      write(6,*)' Returned from analyze_weights_mpi'
      i4_elapsed = ishow_timer()

      call mpi_gatherv(local_var,recvcounts(rank),MPI_REAL,buffer,
     1                 recvcounts,displs,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      if( ierr /= 0) then
        print*,'MPI_GATHERV error in barnes',ierr,recvcounts(rank),rank
        stop
      endif

      write(6,*)' Finished mpi_gatherv'
      i4_elapsed = ishow_timer()

!     Copy buffer, with re-ordering, into var
      ioffset = 0
      do r=0,nPEs-1
        IMAXlocal = IstartIend(2,r)-IstartIend(1,r)+1
        noffset = 0
        do n=1,n_var
          koffset = 0
          do k=1,kmax
            joffset = 0
            do j=1,jmax
              do i=IstartIend(1,r),IstartIend(2,r)
                ip           = i-ioffset
                offset       = joffset+koffset+noffset
                var(i,j,k,n) = buffer(ip+displs(r)+offset)
              enddo
              joffset=joffset+IMAXlocal
            enddo !j
            koffset=koffset+IMAXlocal*jmax
          enddo   !k
          noffset=noffset+IMAXlocal*jmax*kmax
        enddo     !n
        ioffset=ioffset+IMAXlocal
      enddo       !r

      deallocate(local_var)
#else

      allocate(sum_var(max_var,imax,jmax,kmax),STAT=istat_alloc)
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate sum_var'
          stop
      endif

      allocate(sumwt(imax,jmax,kmax),STAT=istat_alloc)
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate sumwt'
          stop
      endif

      sum_var = 0.0         ! Initialize whole array
      sumwt   = 0.0         ! Initialize whole array

!     Sum the weights
      do n=1,ncnt_total
          i = obs_barnes(n)%i
          j = obs_barnes(n)%j
          k = obs_barnes(n)%k

          ilow  = max(i-idelt,1   )
          ihigh = min(i+idelt,imax)
          jlow  = max(j-jdelt,1)
          jhigh = min(j+jdelt,jmax)
          klow  = max(k-kdelt,1)
          khigh = min(k+kdelt,kmax)

!         Calculate J dimension weight ahead of time
          do jj=jlow,jhigh
              fnorm_lut_jj(jj) = fnorm_lut(0,jj-j)
          enddo ! jj

          call get_time_wt(i4time_sys,obs_barnes(n)%i4time,time_wt
     1                    ,istatus)

          do l = 1,n_var ! This is hopefully more efficient
              varl(l) = obs_barnes(n)%value(l)
          enddo ! l

          do kk=klow,khigh
              r0_vert_ob = r0_vert_grid * obs_barnes(n)%vert_rad_rat      
              exponent = -( float(kk-k)/r0_vert_ob )**2
              vertical_weight = exp(exponent)

              if(.true. .and. n .le. 20)then
                  write(6,501)n,i,j,k,kk,vertical_weight,time_wt
     1                       ,obs_barnes(n)%value(1)       
 501              format(5i5,2f8.4,f8.2)
              endif

              weight_term = obs_barnes(n)%weight * vertical_weight 
     1                                           * time_wt      

              if(kmax .gt. 1)then ! 3D field analysis

csms$parallel(grid_dh, <ii>) begin

                  if(n_var .eq. 2)then ! wind analysis
                    if(kk .ge. k)then
                      kk_higher_pair = 0
                    elseif(kk + 2*abs(kk-k) .le. kmax)then
                      kk_higher_pair = 1
                    else
                      kk_higher_pair = 0
                    endif

                    if(kk .le. k)then
                      kk_lower_pair = 0
                    elseif(kk - 2*abs(kk-k) .ge. 1)then
                      kk_lower_pair = 1
                    else
                      kk_lower_pair = 0
                    endif

                    if(n .eq. 1)then
                        write(6,*)'kk,kk_lower_pair,kk_higher_pair',
     1                             kk,kk_lower_pair,kk_higher_pair
                    endif

                    if(kk_lower_pair  .eq. 0 .and. 
     1                 kk_higher_pair .eq. 0)then ! "normal"
                      varl1 = varl(1)
                      varl2 = varl(2)

                      do ii=ilow,ihigh

!                       The exp_m_offset prevents applying the scaling twice
                        weight_i = fnorm_lut(ii-i,0) * exp_m_offset

                        weight_term_i = weight_term * weight_i

                        do jj=jlow,jhigh

                          weight = fnorm_lut_jj(jj) * weight_term_i

                          sum_var(1,ii,jj,kk) = sum_var(1,ii,jj,kk)       
     1                                        + (weight*varl1)       
                          sum_var(2,ii,jj,kk) = sum_var(2,ii,jj,kk)       
     1                                        + (weight*varl2)       

                          sumwt(ii,jj,kk) = sumwt(ii,jj,kk) + weight

                        enddo ! jj
                      enddo ! ii

                    elseif(kk_higher_pair .eq. 1)then ! double up the k s
                      varl1 = varl(1)
                      varl2 = varl(2)

                      kkh = kk + 2*abs(kk-k)

                      do ii=ilow,ihigh

!                       The exp_m_offset prevents applying the scaling twice
                        weight_i = fnorm_lut(ii-i,0) * exp_m_offset

                        weight_term_i = weight_term * weight_i

                        do jj=jlow,jhigh

                          weight = fnorm_lut_jj(jj) * weight_term_i
                          weight_varl1 = weight*varl1
                          weight_varl2 = weight*varl2

                          sum_var(1,ii,jj,kk) = sum_var(1,ii,jj,kk)       
     1                                        + weight_varl1       
                          sum_var(2,ii,jj,kk) = sum_var(2,ii,jj,kk)       
     1                                        + weight_varl2       

                          sumwt(ii,jj,kk) = sumwt(ii,jj,kk) + weight

                          sum_var(1,ii,jj,kkh) = sum_var(1,ii,jj,kkh)       
     1                                        + weight_varl1       
                          sum_var(2,ii,jj,kkh) = sum_var(2,ii,jj,kkh)       
     1                                        + weight_varl2       

                          sumwt(ii,jj,kkh) = sumwt(ii,jj,kkh) + weight

                        enddo ! jj
                      enddo ! ii

                    elseif(kk_lower_pair .eq. 1)then ! skip this k
                      continue

                    else ! this should logically never happen
                      write(6,*)' SEVERE CODE ERROR'
     1                         ,kk_lower_pair,kk_higher_pair
                      stop

                    endif

                  elseif(n_var .eq. 1)then ! temperature analysis
                      varl1 = varl(1)

                      do jj=jlow,jhigh
                      do ii=ilow,ihigh
                          weight = fnorm_lut(ii-i,jj-j) * weight_term

                          sum_var(1,ii,jj,kk) = sum_var(1,ii,jj,kk)       
     1                                        + (weight*varl1)       

                          sumwt(ii,jj,kk) = sumwt(ii,jj,kk) + weight

                      enddo ! ii
                      enddo ! jj

                  else ! generic slightly less efficient method
                      do jj=jlow,jhigh
                      do ii=ilow,ihigh
                          weight = fnorm_lut(ii-i,jj-j) * weight_term

                          do l = 1,n_var
                              sum_var(l,ii,jj,kk) = sum_var(l,ii,jj,kk)       
     1                                            + (weight*varl(l))       
                          enddo ! l
                          sumwt(ii,jj,kk) = sumwt(ii,jj,kk) + weight

                      enddo ! ii
                      enddo ! jj

                  endif
csms$parallel end

              else ! sfc analysis - add in landfrac weight

csms$parallel(grid_dh, <ii>) begin
!                 weight_landfrac = 1.
                  if(n_var .gt. 1)then
                      stop
csms$exit
                  endif

                  do jj=jlow,jhigh
                  do ii=ilow,ihigh

                      if (obs_barnes(n)%mask_sea .ne. 1) then
                          weight_landfrac = 1 ! No landfrac weighting

                      elseif(.true.)then ! threshold method

C
C HACK ALERT: TH Change the way weight_landfrac is calculated.
C
C                         weight_landfrac 
C     1                     = 1.0 - abs( obs_barnes(n)%ldf - ldf(ii,jj) )

                          if (obs_barnes(n)%ldf .gt. 0) then ! land ob
                              if (ldf(ii,jj) .lt. 0.01) then
                                  weight_landfrac = 0.       ! water grid pt
                              else
                                  weight_landfrac = 1.       ! land grid pt
                              endif
                          else                               ! water ob
                              if (ldf(ii,jj) .lt. 0.01) then
                                  weight_landfrac = 1.       ! water grid pt
                              else
                                  weight_landfrac = 0.       ! land grid pt
                              endif
                          endif
C
C END HACK ALERT
C
                      elseif(.false.)then ! hybrid method
                          const = 0.05
                          diffl = ABS(obs_barnes(n)%ldf - ldf(ii,jj))
                          IF (((obs_barnes(n)%ldf .GT. 0.) .AND. 
     1                                (ldf(ii,jj) .EQ. 0.)) .OR.
     1                        ((obs_barnes(n)%ldf .EQ. 0.) .AND. 
     1                                (ldf(ii,jj) .GT. 0.))) THEN
                              weight_landfrac = EXP((-diffl/const)**2)
                          ELSE
                              weight_landfrac = 1
                          ENDIF

                      elseif(.false.)then ! exponential method
                          continue
                      endif

                      weight = fnorm_lut(ii-i,jj-j) * weight_term 
     1                                              * weight_landfrac

                      sum_var(1,ii,jj,1) = sum_var(1,ii,jj,1)
     1                                 + (weight*varl(1))     
                      sumwt(ii,jj,1) = sumwt(ii,jj,1) + weight

                  enddo ! ii
                  enddo ! jj
csms$parallel end

              endif ! 3D or sfc analysis
          enddo ! k
      enddo ! n

      relax = 1.0 ! allows for SUR or SOR with each iteration

      i4_elapsed = ishow_timer()

!     Divide the weights
      do k=1,kmax

csms$parallel(grid_dh, <i>) begin
          do j=1,jmax
          do i=1,imax
              do l = 1,n_var
                  if(sumwt(i,j,k) + wt_b(i,j,k,l) .eq. 0.)then
csms$remove begin
                      var(i,j,k,l) = r_missing_data
csms$remove end
csms$insert           var(i,j,k+(l-1)*kmax) = r_missing_data
                  else
csms$remove begin
                      var(i,j,k,l)=sum_var(l,i,j,k) * relax
csms$remove end
csms$insert           var(i,j,k+(l-1)*kmax)=sum_var(l,i,j,k) * relax
     1                          / (sumwt(i,j,k) + wt_b(i,j,k,l))
                  endif
              enddo ! l
          enddo ! i
          enddo ! j
csms$parallel end

          write(6,*)' Divided level ',k

      enddo ! k

      i4_elapsed = ishow_timer()

      deallocate(sum_var)
      deallocate(sumwt)

#endif

      istatus = 1
      return
      end


      subroutine get_fnorm_max(imax,jmax                         ! I
     1                        ,r0_norm,r0_value_min,fnorm_max)   ! O

      real r0_norm                            ! Output

!     Minimum radius of influence in grid points (for solid data coverage)
      real r0_value_min                       ! Output

      character*8 c8_project

      call get_c8_project(c8_project,istatus)

      if(c8_project .eq. 'CWB')then
          r0_norm = 30.                         ! Grid point lengths
      else
          r0_norm = 10.                         ! Grid point lengths
      endif

      r0_value_min = 1.7

      r0_norm_sq = r0_norm**2
      ratio_fnorm_sq = r0_norm_sq / (r0_value_min**2)
      fnorm_max = (float(imax-1)**2 + float(jmax-1)**2) * ratio_fnorm_sq

      return
      end


        subroutine get_rms_barnesobs(u,ni,nj,nk,r_missing_data,nobs      ! I/O 
     1                        ,n_var,i_var,max_obs,obs_barnes,ncnt_total ! I
     1                        ,rms_obs,rmsu)                             ! O

        include 'barnesob.inc'
        type (barnesob) :: obs_barnes(max_obs)

        dimension u(ni,nj,nk)
	integer nobs
	real residualu
	real sumsq_obs

        nobs = 0
        residualu = 0.
        sumsq_obs = 0.

        write(6,*)' get_rms_barnesobs...'
        write(6,*)'   n    i    j    k          uo         u      diffu'

        if(nk .gt. 1)then
            nprint = 20
        else
            nprint = 100
        endif

        do n = 1,ncnt_total

            i = obs_barnes(n)%i
            j = obs_barnes(n)%j
            k = obs_barnes(n)%k

            ob = obs_barnes(n)%value(i_var)

            if(ob        .ne. r_missing_data .and.
     1          u(i,j,k) .ne. r_missing_data      )then
                nobs = nobs + 1
                diffu = ob - u(i,j,k)
                residualu = residualu + diffu ** 2
                sumsq_obs = sumsq_obs + ob ** 2

                if(nobs .le. nprint)then
                    write(6,1)nobs,i,j,k,ob,u(i,j,k),diffu
 1                  format(4i5,2f12.3,f10.3)
                endif ! nobs
            endif ! Valid Data

        enddo ! n

        if(nobs .gt. 0)then
            rmsu =    sqrt(residualu/float(nobs))
            rms_obs = sqrt(sumsq_obs/float(nobs))
        else
            rmsu = 0.
            rms_obs = 0.
        endif

        return

        end

      
      function fnorm_calc(rrr,r0_norm_sq,exp_offset,expm80)

      real exp_offset

      dist_norm_sq = rrr / r0_norm_sq
      arg = dist_norm_sq - exp_offset

!     Give a floor to the exponential to prevent underflow
      if(arg .le. 80.)then
          fnorm_calc = exp(-arg)
      else
          fnorm_calc = expm80
      endif

      return
      end


      subroutine arrays_to_barnesobs(imax,jmax,kmax                   ! I
     1                              ,r_missing_data                   ! I
     1                              ,varo,weight_3d                   ! I
     1                              ,n_var,max_obs                    ! I
     1                              ,obs_barnes                       ! O
     1                              ,ncnt_total,weight_total          ! O
     1                              ,istatus)                         ! O

!     Map obs from the 3-D arrays into the 'obs_barnes' data structure

      include 'barnesob.inc'
      type (barnesob) :: obs_barnes(max_obs)

      real    varo(imax,jmax,kmax,n_var),weight_3d(imax,jmax,kmax)  

      ncnt_total=0
      weight_total = 0.

csms$insert      call nnt_me(me)
csms$insert      print *, 'got to arrays_to_barnesobs1 processor=',me       

!     Count obs and map obs into data structure
      do k=1,kmax
      do j=1,jmax
      do i=1,imax
          if(varo(i,j,k,1) .ne. r_missing_data)then

              if(ncnt_total .ge. max_obs)then
                  write(6,*)' ERROR: Too many obs for barnes'
     1                     ,ncnt_total,max_obs
                  istatus = 0
                  return
              else

!                 Increment ob counter
                  ncnt_total = ncnt_total + 1
                  weight_total = weight_total + weight_3d(i,j,k)

                  obs_barnes(ncnt_total)%i      = i
                  obs_barnes(ncnt_total)%j      = j
                  obs_barnes(ncnt_total)%k      = k
                  obs_barnes(ncnt_total)%weight = weight_3d(i,j,k)       
                  obs_barnes(ncnt_total)%vert_rad_rat = 1.0

                  do l = 1,n_var
                      obs_barnes(ncnt_total)%value(l) = varo(i,j,k,l)
                  enddo ! l

!                 if(k .eq. 13)write(6,71)ncnt_total,i,j
!       1                           ,uo(i,j,k),vo(i,j,k),weight_3d(i,j,k)
!71               format(1x,i4,2i3,2f11.4,e11.4)

              endif ! Not too many obs

          endif ! We have an ob

      enddo ! i
      enddo ! j
      enddo ! k

csms$insert      call nnt_me(me)
csms$insert      print *, 'got to arrays_to_barnesobs2 processor=',me       

      istatus = 1
      return
      end


      subroutine get_time_wt(i4time_sys,i4time_ob,time_wt,istatus)

      arg = float(i4time_sys-i4time_ob) / 3600.

      if(abs(arg) .gt. 10.)then
          write(6,*)' Error, check times in get_time_wt'
     1             ,i4time_sys,i4time_ob
          time_wt = 1.0
          istatus = 0
          return

      else
          time_wt = exp(-(arg**2))

      endif

!     time_wt = 1.0 ! test 

      istatus = 1
      return
      end

      subroutine analyze_weights_mpi(                           
     1          n_var,max_obs,obs_barnes                                    ! I
     1         ,imax,jmax,kmax                                              ! I
     1         ,wt_b                                                        ! I
     1         ,i4time_sys                                                  ! I
     1         ,r_missing_data,r0_vert_grid                                 ! I
     1         ,ncnt_total,idelt,jdelt,kdelt                                ! I
     1         ,fnorm_lut                                                   ! I
     1         ,Istart,Iend                                                 ! I
     1         ,var)                                                        ! O

      implicit none

!     Argument list declarations

      integer, INTENT(IN) :: n_var,max_obs                      

      include 'barnesob.inc'
      type (barnesob), INTENT(IN) :: obs_barnes(max_obs)

      integer, INTENT(IN) :: imax,jmax,kmax,i4time_sys
      integer, INTENT(IN) :: ncnt_total,idelt,jdelt,kdelt

      real, INTENT(IN) :: r_missing_data,r0_vert_grid

      real, INTENT(IN) :: wt_b(Istart:Iend,jmax,kmax,n_var)         

      real, INTENT(IN)::fnorm_lut(-(imax-1):(imax-1),-(jmax-1):(jmax-1))     
      integer,intent(IN) :: Istart,Iend                         ! Input

      real, INTENT(OUT) :: var(Istart:Iend,jmax,kmax,n_var)     ! Output

!     Local declarations

      real, allocatable, dimension(:,:,:,:) :: sum_var          ! Local
      real, allocatable, dimension(:,:,:) :: sumwt              ! Local
      real varl(2)

      real exponent,vertical_weight,time_wt,weight_term,r0_vert_ob
      real weight

      integer istat_alloc,ilow,ihigh,jlow,jhigh,klow,khigh,ii,jj,kk
      integer i,j,k,l,n,istatus 

      integer i4_elapsed, ishow_timer

!     Start executable statements

      i4_elapsed = ishow_timer()

      allocate(sum_var(max_var,Istart:Iend,jmax,kmax),STAT=istat_alloc)
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate sum_var'
          stop
      endif

      allocate(sumwt(Istart:Iend,jmax,kmax),STAT=istat_alloc)
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate sumwt'
          stop
      endif

      write(6,*)' analyze_weights_mpi: initialize arrays ',Istart,Iend

      sum_var = 0.0         ! Initialize whole array
      sumwt   = 0.0         ! Initialize whole array

      i4_elapsed = ishow_timer()

      write(6,*)' analyze_weights_mpi: start summing ',Istart,Iend

!     Sum the weights
      do n=1,ncnt_total
          i = obs_barnes(n)%i
          j = obs_barnes(n)%j
          k = obs_barnes(n)%k
          ilow  = max(i-idelt,Istart)
          ihigh = min(i+idelt,Iend)
          jlow  = max(j-jdelt,1)
          jhigh = min(j+jdelt,jmax)
          klow  = max(k-kdelt,1)
          khigh = min(k+kdelt,kmax)
          call get_time_wt(i4time_sys,obs_barnes(n)%i4time,time_wt
     1                    ,istatus)

          do l = 1,n_var ! This is hopefully more efficient
              varl(l) = obs_barnes(n)%value(l)
          enddo ! l

          do kk=klow,khigh
              r0_vert_ob = r0_vert_grid * obs_barnes(n)%vert_rad_rat      
              exponent = -( float(kk-k)/r0_vert_ob )**2
              vertical_weight = exp(exponent)

              if(i .ge. Istart .and. i .le. Iend .and. n .le. 20)then
                  write(6,501)n,i,j,k,kk,vertical_weight,time_wt
     1                       ,obs_barnes(n)%value(1)
     1                       ,r0_vert_ob,exponent       
 501              format(5i5,2f8.4,f8.2,2f8.4)
              endif

              weight_term = obs_barnes(n)%weight * vertical_weight 
     1                                           * time_wt      

              do jj=jlow,jhigh
              do ii=ilow,ihigh
                  weight = fnorm_lut(ii-i,jj-j) * weight_term

                  do l = 1,n_var
                      sum_var(l,ii,jj,kk) = sum_var(l,ii,jj,kk)       
     1                                    + (weight*varl(l))       
                  enddo ! l
                  sumwt(ii,jj,kk) = sumwt(ii,jj,kk) + weight

              enddo ! ii
              enddo ! jj
          enddo ! k
      enddo ! n

      write(6,*)' analyze_weights_mpi: finished summing ',Istart,Iend

      i4_elapsed = ishow_timer()

!     Divide the weights
      do k=1,kmax
          do j=1,jmax
          do i=Istart,Iend
              do l = 1,n_var
                  if(sumwt(i,j,k) + wt_b(i,j,k,l) .eq. 0.)then
                      var(i,j,k,l) = r_missing_data
                  else
                      var(i,j,k,l)=sum_var(l,i,j,k) 
     1                          / (sumwt(i,j,k) + wt_b(i,j,k,l))
                  endif
              enddo ! l
          enddo ! i
          enddo ! j

          if(k .eq. 1 .or. k .eq. kmax)then
              write(6,*)' Divided level ',k,Istart,Iend
          endif

      enddo ! k

      i4_elapsed = ishow_timer()

      deallocate(sum_var)
      deallocate(sumwt)

      return
      end
