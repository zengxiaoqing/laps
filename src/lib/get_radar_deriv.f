       subroutine get_radar_deriv(nx,ny,nz,dx,r_miss,                   ! I
     1                           radar_ref_3d,clouds_3d,cld_hts,        ! I
     1                           temp_3d,heights_3d,pres_3d,            ! I
     1                           ibase_array,itop_array,thresh_cvr,     ! I
     1                           vv_to_height_ratio_Cu,                 ! I
     1                           cldpcp_type_3d,                        ! I
     1                           w_3d,istatus)                          ! I/O
!      Adan Teng
!      this rutine calculate the cloud bogus omega in radar echo area
!      First does the convective and stratiform region seperate
!      In convective region using the parabolic profile as in vv.f
!        but using the maximum terminal velocity deduce from reflectivity
!        to be the maximum value in parabolic profile
!      In stratiform region parabolic profile assigns above(behind)
!        the melting level with maximum(minimum) value dependent 
!        on grid space and cloud depth
       Implicit none
       integer nx,ny,nz,istatus                !input
       integer cldpcp_type_3d(nx,ny,nz)        !input
       integer ibase_array(nx,ny)              !input
       integer itop_array(nx,ny)               !input
       real dx,r_miss                        !input
       real radar_ref_3d(nx,ny,nz)           !input (dbZ)
       real clouds_3d(nx,ny,nz)              !input
       real cld_hts(nx,ny,nz)                !input
       real temp_3d(nx,ny,nz)                !input
       real heights_3d(nx,ny,nz)             !input
       real pres_3d(nx,ny,nz)                !input
       real thresh_cvr                       !input
       real w_3d(nx,ny,nz)                   !input and output
       real vv_to_height_ratio_Cu
!      temprary variables
       integer i,j,k
       integer nxx,nyy,ier
       integer str_con_index(nx,ny)
       integer i_strcon
       integer index_random(nx,ny)
       real  radar_2d_max(nx,ny)             ! dbZ (while in this routine)
       logical l_cloud
       logical l_use_random                  ! utilize random indices
       real temp_1d(nz)
       real heights_1d(nz)
       real pressure_mb(nz)
       real pressure_pa(nz)
       integer iarg
       integer cloud_type_1d(nz)
       real  radar_ref_max ! dbZ
       real w_1d(nz)
       real w_to_omega,om_orig
       integer strcon
       integer rand_index
       integer dbz(nx,ny)
       real dx1
       integer nstep 

       l_use_random = .false.
       i_strcon = 2 ! 0 - automatic index, 2 - convective

       istatus = 1
       if (dx .lt. 1500.) then
        dx1 = dx
        nxx = nx
        nyy = ny
       else
        nstep = int((dx+500.)/1000.) 
        dx1 = dx /nstep
        nxx = (nx-1)*nstep+1
        nyy = (ny-1)*nstep+1
       endif

       write(6,*)' calling get_con_str, nxx/nyy = ',nxx,nyy

!      Determine convective and stratiform regions
       call get_con_str(nx,ny,nz,nxx,nyy,radar_ref_3d,           ! I
     1                  pres_3d,temp_3d,                         ! I
     1                  str_con_index,                           ! O
     1                  radar_2d_max,r_miss,ier,index_random,    ! O
     1                  dx1)                                     ! O
       if( ier .eq. 0) then
        write(*,*)'Cannot separate convection and stratiform region'
        istatus = 0
        return
       endif

       do i = 1, nx
         do j = 1, ny
          l_cloud = .false.
          if( ibase_array(i,j) .lt. itop_array(i,j)) then
           l_cloud = .true.
          endif
          do k = 1, nz    ! Initialize
           temp_1d(k) = temp_3d(i,j,k)-273.15
           heights_1d(k) = heights_3d(i,j,k)
           pressure_mb(k) = pres_3d(i,j,k) / 100.
           pressure_pa(k) = pres_3d(i,j,k)
           iarg = cldpcp_type_3d(i,j,k)
           cloud_type_1d(k) = iarg - iarg/16*16
           w_1d(k) = r_miss
          enddo
          radar_ref_max = radar_2d_max(i,j)

          if(i_strcon .eq. 0)then
              strcon = str_con_index(i,j) 
          else
              strcon = i_strcon
          endif

          if(l_use_random)then
              rand_index = index_random(i,j)
          else
              rand_index = 2
          endif

          if ( l_cloud ) then
           if ( strcon .ne. 0 ) then
            if(radar_ref_max .ne. r_miss)then
              call radar_bogus_w
     1           (dx, cloud_type_1d, heights_1d, temp_1d,
     1            radar_ref_max, strcon, nz, w_1d,
     1            vv_to_height_ratio_Cu,                                    ! I
     1            rand_index)
            endif
            do k = 1, nz
             if( w_1d(k) .ne. r_miss ) then
!             om_orig = w_3d(i,j,k)
              w_3d(i,j,k) = w_to_omega(w_1d(k), pressure_pa(k))
!             if(w_3d(i,j,k) .ne. om_orig .or. 
!    1           radar_ref_max .ne. r_miss     )then
!                 write(6,*)' orig/om/w/ref='
!    1                     ,om_orig,w_3d(i,j,k),w_1d(k),radar_ref_max
!    1                     ,i,j,k
!             endif
             endif
            enddo 
           endif   ! strcon
          endif   ! l_cloud
         enddo   !j
       enddo    !i
 
       return
       end

       subroutine radar_bogus_w
     1           (dx, cloud_type, heights, temp,
     1            radar_ref_max, strcon, nz, w,
     1            vv_to_height_ratio_Cu,                                    ! I
     1            rand_index)
       
       Implicit none
       integer nz, cloud_type(nz), strcon
       integer rand_index
       real heights(nz), temp(nz), w(nz)
       real dx, radar_ref_max
       real vv_to_height_ratio
       real vv_to_height_ratio_Cu
 
!      data vv_to_height_ratio /0.5/
       real ratio, vv, parabolic_vv_profile
       real parabolic_vv_profile1
       integer k, k1, kbase, ktop, kmiddle
       real zbase, ztop
       real ratio_radar
       real depth, vvmax
       logical l_always_use_cloud,l_always_use_radar

!   Cloud Type      /'  ','St','Sc','Cu','Ns','Ac','As','Cs','Ci','Cc','Cb'/
!   Integer Value     0     1    2    3    4    5    6    7    8    9   10

       l_always_use_radar = .true. ! original was .false.
       l_always_use_cloud = .true. ! original was .false.

! Put in the vv's for cumuliform clouds with radar reflectivity as determined
! by cloud type at the base
       vv_to_height_ratio = vv_to_height_ratio_Cu
       ratio = vv_to_height_ratio / dx 
    
       if(.not. l_always_use_cloud)then
        Do k = 1, nz
         If (cloud_type(k) .eq. 3  .OR.  cloud_type(k) .eq. 10) then
          kbase = k
          Go to 10
         End if
        End do
       else
        Do k = 1, nz
         If (cloud_type(k) .ne. 0) then
          kbase = k
          Go to 10
         End if
        End do
       endif
        Go to 100

10      Do k = kbase, nz
c        If (cloud_type(k) .eq. 3  .OR.  cloud_type(k) .eq. 10) then
!     change to the cloudtop by Adan
         If (cloud_type(k) .ne. 0) then
          ktop = k
         Else
          Go to 20
         End if
        End do

20      k1 = k          ! save our place in the column
        zbase = heights(kbase)
        ztop  = heights(ktop)
! Add for recomputing ratio by using terminal velocity dervied by radar
!         reflectivity  (Jen-Hsin Teng, Adan)
        if ( strcon .eq. 2 ) then    ! convective
         if (rand_index .eq. 2) then ! higher index 2/3 probability
          ratio_radar=0.
          depth = ztop - zbase
          if (depth.ne.0) then
           if(radar_ref_max .gt. 0.) then
!           vvmax=4.32*radar_ref_max**0.0714286
            vvmax=radar_ref_max * 0.1                
            write(6,*)' radar_ref_max/vvmax = ',radar_ref_max,vvmax
            ratio_radar=vvmax / depth / 1.1
            if(ratio_radar .gt. ratio .or. l_always_use_radar)then
             ratio = ratio_radar
            endif
           endif
          else
           write(6,*) 'depth =',depth, 'No changes to ratio'
          endif

         else          ! for rand_index = 1 (1/3 probability)
          ratio_radar=0.
          depth = ztop - zbase
          if (depth.ne.0) then
           if(radar_ref_max .gt. 0. ) then
!           vvmax=4.32*radar_ref_max**0.0714286
            ratio_radar=vvmax / depth / 1.1
!           if(ratio_radar .gt. ratio) ratio = ratio_radar
           endif
          else
           write(6,*) 'depth =',depth, 'No changes to ratio'
          endif

         endif ! rand_index
 
         Do k = 1, ktop
          vv = Parabolic_vv_profile (zbase, ztop, ratio, heights(k))
          If (vv .gt. 0.) then
           w(k) = vv
          End if
         End do

        elseif (strcon .eq. 1) then    ! stratiform
         if (temp(kbase) .le. 0. .or. 
     1       temp(ktop)  .gt. 0.) then ! doesn't straddle freezing level

                                       ! Are warm layers fully below this level
                                       ! treated the same as cold layers fully
                                       ! above?

          if (rand_index .eq. 2) then  ! higher index 2/3 probability
           ratio_radar=0.
           depth = ztop - zbase
           if (depth.ne.0) then
            vvmax = 1.
            ratio_radar=vvmax / depth / 1.1
            if(ratio_radar .lt. ratio) ratio = ratio_radar
           else
            write(6,*) 'depth =',depth, 'No changes to ratio'
           endif

          else     !  for rand_index = 1 (1/3 probability)
           ratio = ratio * 0.001
           ratio_radar=0.
           depth = ztop - zbase
           if (depth.ne.0) then
            vvmax = 0.001
            ratio_radar=vvmax / depth / 1.1
            if(ratio_radar .lt. ratio) ratio = ratio_radar
           else
            write(6,*) 'depth =',depth, 'No changes to ratio'
           endif

          endif ! rand_index

          Do k = 1, ktop
            vv = Parabolic_vv_profile (zbase, ztop, ratio, heights(k))
            If (vv .gt. 0.) then
             w(k) = vv
            End if
          End do

         else                           ! straddles freezing level
          do k = kbase, ktop
           if (temp(k) .le. 0.) then 
            kmiddle = k
            go to 30
           endif
          enddo
30        continue
          ztop = heights(kmiddle) 

          if( rand_index .eq. 2) then ! 2/3 probability
           ratio_radar=0.
           depth = ztop - zbase
           if (depth.ne.0) then
            vvmax = 1.
            ratio_radar=vvmax / depth /1.1 
            if(ratio_radar .lt. ratio) ratio = ratio_radar
           else
            write(6,*) 'depth =',depth, 'No changes to ratio'
           endif

          else    ! for rand_index = 1 (1/3 probabilty)
           ratio = ratio * 0.001
           ratio_radar=0.
           depth = ztop - zbase
           if (depth.ne.0) then
            vvmax = 0.001
            ratio_radar=vvmax / depth /1.1
            if(ratio_radar .lt. ratio) ratio = ratio_radar
           else
            write(6,*) 'depth =',depth, 'No changes to ratio'
           endif

          endif ! rand_index

          do k = 1, kmiddle
           vv = Parabolic_vv_profile (zbase, ztop, ratio, heights(k))
           if ( vv .gt. 0.) then
            w(k) = -1. * vv
           endif
          enddo

          zbase = ztop
          ztop = heights(ktop)
          if (rand_index .eq. 2) then ! 2/3 probability
           ratio_radar=0.
           depth = ztop - zbase
           if (depth.ne.0) then
            vvmax = 1.
            ratio_radar=vvmax / depth 
            if(ratio_radar .lt. ratio) ratio = ratio_radar
           else
            write(6,*) 'depth =',depth, 'No changes to ratio'
           endif

          else      ! for rand_index = 1 (1/3 probability)
           ratio = ratio * 0.001
           ratio_radar=0.
           depth = ztop - zbase
           if (depth.ne.0) then
            vvmax = 0.001
            ratio_radar=vvmax / depth
            if(ratio_radar .lt. ratio) ratio = ratio_radar
           else
            write(6,*) 'depth =',depth, 'No changes to ratio'
           endif

          endif ! rand_index

!         Above the freezing layer while within the echo apply a parabolic 
!         profile that extends only inside this layer
          do k = kmiddle, ktop
           vv = Parabolic_vv_profile1 (zbase, ztop, ratio, heights(k))
           if ( vv .gt. 0.) then
            w(k) = vv
           endif
          enddo
         endif  ! temp of layer doesn't straddle freezing level?
        endif   ! strcon
100     continue
        return
        end
 

          subroutine get_con_str(nx,ny,nz,nxx,nyy,radar_ref_3d,pres_3d,
     1                  temp_3d,str_con_index,radar_2d_max,r_miss,ier,
     1                  index_random,dx1) 
!
!        Description : Interpolate radar reflectivity to 1 km
!                      resolution, then call dfconstr to define
!                      convective and stratiform region
!                      Also calculate the maximum reflectivity in
!                      vertical column
!        Arguments : 
!      I/O/W    name,          type,        description
!        I      nx             integer      x-dimension
!        I      ny             integer      y-dimension
!        I      nz             integer      z-dimension 
!        I      nxx            integer      x-dimension for 1 km resolution
!        I      nyy            integer      y-dimension for 1 km resolution
!        I      radar_ref_3d   real array   3d reflectivity data (dbZ)
!        O      str_con_index  int array    2d conevction or stratiform index
!                                           =0, missing data area
!                                           =1, stratiform region
!                                           =2, convective region
!        I      pres_3d        real array   3d pressure data (pa)
!        I      temp_3d        real array   3d temperature data (K)
!        O      radar_2d_max   real array   2d reflectivity maximum data
!        I      r_miss         real         missing data
!        O      ier            int          =0, failure
!                                           =1, success 
!        O      index_random   int array    2d random number index
!                                           =0 missing data area
!                                           =1 random number less than 
!                                              devide_random
!                                           =2 random number greater than 
!                                              devide_random
!        I      dx1            real         grid size of small resolution
!
          implicit none
          integer nx, ny, nz, nxx, nyy
          real    radar_ref_3d(nx,ny,nz)
          real    pres_3d(nx,ny,nz)
          real    temp_3d(nx,ny,nz)
          integer str_con_index(nx,ny)
          integer index_random(nx,ny)
          real    radar_2d_max(nx,ny)
          real    r_miss
          real    devide_random
          real    dx1
          integer ier
          integer i,j,k,ii,jj
          real    ref_2d(nxx,nyy)
          integer icon(nxx,nyy)
          integer istep,jstep
          integer ref_max
          integer i1, i2, j1, j2
          integer idx, jdy
          real  a1, a2
          integer ierdf
          integer iconmax, icon0, icon1, icon2, icont
          
!
!         Initialize the data
!
          ier = 1
          ierdf = 1
          devide_random = 0.33
          do i = 1, nx
           do j = 1, ny
            radar_2d_max(i,j) = r_miss
            str_con_index(i,j) = 0
            index_random(i,j) = 0
           enddo
          enddo
 
!
!         calculate the maximum reflectivity in vertical column
!
          do i = 1, nx
           do j = 1, ny 
            ref_max = -10.
            do k = 1, nz
              if( radar_ref_3d(i,j,k) .ne. r_miss .and.
     1            radar_ref_3d(i,j,k) .gt. ref_max ) then
                ref_max = radar_ref_3d(i,j,k)
              endif
            enddo
            if(ref_max .lt. -5.) then
             radar_2d_max(i,j) = r_miss
            else
             radar_2d_max(i,j) = ref_max
            endif
           enddo
          enddo
!
!        define the level of reflectivity for convective and
!               stratiform region seperate
!
          
           
! Change dBZ to Z
          do i = 1, nx
           do j = 1, ny
            if (radar_2d_max(i,j) .ne. r_miss) then
             radar_2d_max(i,j) = 10.**(0.1*radar_2d_max(i,j))
            endif
           enddo
          enddo
!
!         interpolate to small resolution
! 
          
          do i = 1, nxx
           do j = 1, nyy
             ref_2d(i,j) = r_miss
             icon(i,j) = 0
           enddo
          enddo
          istep = (nxx-1)/(nx-1)
          jstep = int(real(nxx-1)/real(nx-1)+0.0001) 
          if (istep .ne. jstep .or. istep .lt. 1) then
           write(*,*) 'Error in subroutine get_con_str'
           ier = 0
           return
          endif
          
          do i = 1, nxx, istep
           ii = 1+(i-1)/istep
           do j = 1, nyy, istep
            jj = 1+(j-1)/istep 
            ref_2d(i,j) = radar_2d_max(ii,jj)
           enddo
          enddo

          if( istep .eq. 1) go to 101

!         
!         interpolate
!
          do i = 1, nxx-istep, istep
           i1 = 1+(i-1)/istep 
           i2 = i1 + 1
           do j = 1, nyy-istep, istep
            j1 = 1+(j-1)/istep
            j2 = j1 + 1
            if( radar_2d_max(i1,j1) .ne. r_miss .and.
     1          radar_2d_max(i1,j2) .ne. r_miss .and.
     1          radar_2d_max(i2,j1) .ne. r_miss .and.
     1          radar_2d_max(i2,j2) .ne. r_miss) then
              do ii = i+1, i+istep-1
               idx = ii - i
               a1 = (real(istep-idx)*radar_2d_max(i1,j1)+
     1              real(idx)*radar_2d_max(i2,j1))/real(istep)
               a2 = (real(istep-idx)*radar_2d_max(i1,j2)+
     1              real(idx)*radar_2d_max(i2,j2))/real(istep)
               do jj = j+1, j+istep-1
                jdy = jj - j
                ref_2d(ii,jj)= (real(istep-jdy)*a1+
     1                          real(jdy)*a2)/real(istep)
               enddo
              enddo
            endif
           enddo
          enddo

 101      continue

! Change Z to dBZ
          do i = 1, nxx
           do j = 1, nyy
            if(ref_2d(i,j) .ne. r_miss) then
             ref_2d(i,j) = 10.* alog10(ref_2d(i,j))
            endif
           enddo
          enddo
          do i = 1, nx
           do j = 1, ny
            if (radar_2d_max(i,j) .ne. r_miss) then
             radar_2d_max(i,j) = 10.*alog10(radar_2d_max(i,j))
            endif
           enddo
          enddo

          call dfconstr(nxx,nyy,r_miss,ref_2d,icon,ierdf,dx1)
          if (ierdf .eq. 0 ) then
           write(*,*)'No radar echo to define convective region'
           ier = 0
           return
          endif

!
!         change convective index to original resolution
!
          do i = 1, nxx, istep
           ii = 1+(i-1)/istep 
           do j = 1, nyy, istep
            jj = 1+(j-1)/istep 
!            if ( ii .eq. 1 .or. ii .eq. nx .or.
!     1           jj .eq. 1 .or. jj .eq. ny .or.
!     1           istep .eq. 1) then
              str_con_index(ii,jj) = icon(i,j)
!            else
!             jstep = (istep+1)/2
!             iconmax = 0
!             icon0 = 0
!             icon1 = 0
!             icon2 = 0
!             icont = 0
!             do i1 = i-jstep, i+jstep
!              do j1 = j-jstep, j+jstep
!                icont = icont + 1
!                if (icon(i1,j1) .eq. 2) then
!                 icon2 = icon2 + 1
!                elseif (icon(i1,j1) .eq. 1) then 
!                 icon1 = icon1 +1
!                else
!                 icon0 = icon0 +1
!                endif
!                if( icon(i1,j1) .gt. iconmax) then
!                 iconmax = icon(i1,j1)
!                endif
!              enddo
!             enddo
!             if( icon(i,j) .eq. 0) then
!              str_con_index(ii,jj) = 0
!             elseif(icon(i,j) .eq. 1) then
!              str_con_index(ii,jj) = 1
!             elseif(icon(i,j) .eq. 2 .and.
!     1        (icon2 .gt. 2*icon1)) then
!              str_con_index(ii,jj) = 2
!             else
!              str_con_index(ii,jj) = 1
!             endif
!            endif
           enddo
          enddo 
          call df_random_index(nx,ny,str_con_index,index_random,
     1                         devide_random)
          return
          end

          subroutine df_random_index(nx,ny,str_con_index,index_random
     1                        ,devide_random)
!        
!        Subroutine/Function : df_random_index
!
!       Usage : call df_random_index(nx,ny,str_con_index,index_random
!    1                        ,devide_random) 
!
!       Description : To define random number index in convective and
!                     stratiform region
!
!       Arguments :
!      I/O/W      name,          type,        description
!        I        nx             integer      x-dimension
!        I        ny             integer      y-dimension
!        I        str_con_index  int array    2d conevction or stratiform index
!                                             =0, missing data area
!                                             =1, stratiform region
!                                             =2, convective region
!        O        index_random   int array    2d random number index
!                                             =0, missing data area
!                                             =1, smaller random number
!                                             =2, larger random number
!        I        devide_random  real number  random number threshhold 
!
         implicit none
         integer nx, ny
         integer str_con_index(nx, ny)
         integer index_random(nx, ny)
         real    devide_random
         integer i, j, k
         integer seed(1)
         character*9 a9time
         integer i4time
         integer istatus
         integer idate
         real    harver
         integer count1
         call get_systime(i4time, a9time, istatus)
         read(a9time, '(i9)') idate
         seed(1) = idate
         count1 = seed(1)
         call random_seed
         call random_seed(size = k)
         do i = 1, count1
          call random_number(harver)
         enddo
!
!        get random number index in stratiform region
!
         do j = 1, ny
          do i = 1, nx
           if(str_con_index(i,j) .eq. 1) then
            call random_number(harver) 
            if(harver .lt. devide_random) then
             index_random(i,j) = 1 ! probability = devide_random
            else
             index_random(i,j) = 2 ! probability = 1. - devide_random
            endif
           endif
          enddo      ! nx
         enddo      ! ny
!
!        get random number index in convective region
!
         do j = 1, ny
          do i = 1, nx
           if(str_con_index(i,j) .eq. 2) then
            call random_number(harver)
            if(harver .lt. devide_random) then
             index_random(i,j) = 1 ! probability = devide_random
            else
             index_random(i,j) = 2 ! probability = 1. - devide_random
            endif
           endif
          enddo      ! nx
         enddo      ! ny

         return
         end

                 
          subroutine dfconstr(nx,ny,r_miss,ref_2d,icon,ier,dx1)
!        
!      Subroutine/Function : dfconstr
!
!      Usage : call dfconstr(nx,ny,r_miss,ref_2d,icon,ier)
!
!      Description : To define convective and stratiform region
!      
!      Arguments :
!      I/O/W      name,        type,        description
!        I        nx           integer      x-dimension
!        I        ny           integer      y-dimension
!        I        r_miss       real         missing value
!        I        ref_2d       real array   2d reflectivity data
!        O        icon         int array    2d conevction or stratiform index
!                                           =0, missing data area
!                                           =1, stratiform region
!                                           =2, convective region
!        O        ier          integer      =0, failure
!                                           =1, success
!        I        dx1          real         grid size
!        Original version: Shiow-Ming Deng
!        Recently version: Jen-Hsin Teng in 15-Aug-2003
!
         implicit none
         integer  nx, ny, ier
         real     r_miss
         real     ref_2d(nx,ny)
         integer  icon(nx,ny) 
         integer  i, j, nn, ii, jj
         real     dzmax, dzmin, sumdz, detz, gdz
         real dx1
         integer  irr
         integer  njump
         
         ier = 0
         njump = int(10000./dx1)
         if( nx .lt. 3*njump .or. ny .lt. 3*njump) then
          write(*,*) 'Error in subroutine dfconstr.'
          write(*,*) 'Your Domain is too small.'
          write(*,*) '   nx:',nx,'   ny:',ny
          return
         endif
         dzmax = 42.
         dzmin = 20.
!     Initial icon data
         do j = 1,ny
          do i = 1,nx
           icon(i,j)=0
          enddo
         enddo

         do 100 j = njump+1, ny-njump, 1
          do 100 i = njump+1, nx-njump, 1
           if(ref_2d(i,j) .eq. r_miss) go to 100
           ier = 1
           if(ref_2d(i,j) .lt. dzmin ) then
            icon(i,j) = 1
            go to 100
           endif
           if(ref_2d(i,j) .ge. dzmax ) then
            icon(i,j) = 2
            go to 100
           endif
           icon(i,j) = 1
           sumdz = 0.
           nn = 0
           do jj = j-njump, j+njump, 1
            do ii = i-njump, i+njump, 1
             irr=(jj-j)**2+(ii-i)**2
             if (irr .le. 100) then
               nn = nn + 1
               if( ref_2d(ii,jj) .ne. r_miss ) then
                sumdz = sumdz + ref_2d(ii,jj)
               else
                sumdz = sumdz - 5.
               endif
             endif
            enddo
           enddo
           if ( nn .gt. 0 ) then
            gdz = sumdz/real(nn)
           else
            go to 100
           endif
           if ( gdz .lt. 0. ) then
            detz = 12 
           else
            detz = (12.-gdz**2./147.)
           endif
           if ( (ref_2d(i,j)-gdz) .gt. detz ) icon(i,j) = 2
 100    continue

        return
        end
