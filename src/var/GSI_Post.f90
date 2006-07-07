      PROGRAM gsi2laps

      IMPLICIT NONE

!--- decalare variables ---------------------------------------

!//common variables
      INTEGER*4 :: imax,jmax,kmax
      INTEGER*4 :: imax1,jmax1,kmax1  
      INTEGER*4 :: i,j,k,ii,jj,kk,istatus
      INTEGER*4 :: i4time,namelen
      CHARACTER :: filename*125

!// for subroutine read_gsi_output
      INTEGER*4,PARAMETER :: ii1=7 
      REAL*4,ALLOCATABLE,DIMENSION(:,:,:) :: uu1,vv1,st,qvar
      REAL*4,ALLOCATABLE,DIMENSION(:,:) :: psfc
      REAL*4,ALLOCATABLE,DIMENSION(:) :: eta
      REAL*4 :: ptop 

!// for subroutine get_systime
      CHARACTER*9 :: a9time

!// for subroutine get_pres_1d
      REAL*4,ALLOCATABLE,DIMENSION(:) :: pres_1d

!// for subroutine get_r_missing_data
      REAL*4 :: r_missing_data 

!// for subroutine get_pres_3d
      REAL*4,ALLOCATABLE,DIMENSION(:,:,:) :: pres_3d 

!// for subroutine get_modelfg_3d
      REAL*4,ALLOCATABLE,DIMENSION(:,:,:) :: heights_3d 

!// for subroutine get_modelfg_3d 
      REAL*4,ALLOCATABLE,DIMENSION(:,:,:) :: t_laps_bkg 

!// for subroutine dryairmass 
      REAL*4,ALLOCATABLE,DIMENSION(:,:) :: dam,pdam 

!// for subroutine stagger
      REAL*4,ALLOCATABLE,DIMENSION(:,:,:) :: uu2,vv2,iqvar 

!// for subroutine mass2p
      REAL*4,ALLOCATABLE,DIMENSION(:,:,:) :: ist,uvar1,vvar1,&
                                             tvar1,qvar1,rh 

!// for subroutine read_static_grid
      REAL*4,ALLOCATABLE,DIMENSION(:,:) :: lat,lon,topo     

!// for subroutine get_grid_spacing_actual
      REAL*4 :: grid_spacing_m
      INTEGER*4 :: icen,jcen

!// for subroutine wind_post_process
      REAL*4,ALLOCATABLE,DIMENSION(:,:) :: rk_terrain,&
                                           uanl_sfcitrp,vanl_sfcitrp
      CHARACTER*3 :: var_3d(3)=(/'U3','V3','OM'/)
      CHARACTER*4 :: units_3d(3)=(/'M/S ','M/S ','PA/S'/) 
      CHARACTER*125 :: comment_3d(3)=(/'3DWIND',&
                                     '3DWIND','3DWIND'/)
      LOGICAL,PARAMETER :: l_grid_north_out = .true.		
      REAL*4 :: height_to_zcoord2,fraclow,frachigh
      integer*4 :: klow,khigh    
 
!// for subroutine make_fnam_lp
      CHARACTER*9 :: asc9_tim

!// for subroutine write_temp_anal
      REAL*4,ALLOCATABLE,DIMENSION(:,:,:,:) ::  output_4d
      CHARACTER*10 :: comment_2d

!// for subroutine writefile
      INTEGER*4,ALLOCATABLE,DIMENSION(:) :: lvl
      CHARACTER*125,PARAMETER :: &
              commentline = 'maps with clouds and surface effects only'

!// for subroutine lh3_compress
      REAL*4 :: t_ref,make_rh
      REAL*4,ALLOCATABLE,DIMENSION(:,:,:) :: rh1

!//Variables for lwm: YUANFU XIE
      REAL*4,ALLOCATABLE,DIMENSION(:,:,:) :: out_sfc_2d
      CHARACTER*3 :: ext_xie,var_a_xie(2),units_a_xie(2)
      CHARACTER*7 :: comment_a_xie(2)

!-----> main program <-----------------------------------------------------
		       
! To read path of wrf_inout.
      CALL get_directory('log',filename,namelen)        ! Yuanfu: use LAPS;
      filename = filename(1:namelen-4)//'tmpdirEjetlaps/wrf_inout'
      !filename = 'wrf_inout.nc'

! To read gsi dimension
! input :: filename
! output:: imax,jmax,kmax
      call read_gsi_dim(imax,jmax,kmax,filename,istatus)
      if( istatus.ne.0 )then
          print*,'Cannot read dimensions of wrf_inout.'
          print*,'filenam: ',filename 
          call exit(1)
      endif
      imax1=imax-1
      jmax1=jmax-1
      kmax1=kmax-1

! To allocate variables  
      allocate(output_4d(imax,jmax,kmax,2))
      allocate(uu1(imax,jmax,kmax))
      allocate(uu2(imax,jmax,kmax))
      allocate(vv1(imax,jmax,kmax))
      allocate(vv2(imax,jmax,kmax))
      allocate(st(imax,jmax,kmax))
      allocate(ist(imax,jmax,kmax))
      allocate(uvar1(imax,jmax,kmax))
      allocate(vvar1(imax,jmax,kmax))
      allocate(tvar1(imax,jmax,kmax))
      allocate(qvar(imax,jmax,kmax))
      allocate(iqvar(imax,jmax,kmax))
      allocate(qvar1(imax,jmax,kmax))
      allocate(rh(imax,jmax,kmax))
      allocate(heights_3d(imax,jmax,kmax))
      allocate(out_sfc_2d(imax,jmax,2))
      allocate(pres_3d(imax,jmax,kmax))
      allocate(t_laps_bkg(imax,jmax,kmax))
      allocate(psfc(imax,jmax))
      allocate(uanl_sfcitrp(imax,jmax))
      allocate(vanl_sfcitrp(imax,jmax))
      allocate(topo(imax,jmax))
      allocate(lat(imax,jmax))
      allocate(lon(imax,jmax))
      allocate(rk_terrain(imax,jmax))
      allocate(dam(imax,jmax))
      allocate(pdam(imax,jmax))
      allocate(eta(kmax))
      allocate(pres_1d(kmax))
      allocate(lvl(kmax))
      allocate(rh1(imax,jmax,kmax))

! To write out variables of wrf_inout(after GSI analysis)
! input :: ii1,imax,jmax,kmax,imax1,jmax1,kmax1,filename
! output:: uu1(u wind),vv1(v wind),st(temperature),
!          psfc(surface pressure),ptop(top pressure),
!          eta,qvar(specific humidity) 
      call read_gsi_output(ii1,imax,jmax,kmax,imax1,jmax1,kmax1, &
                       filename,uu1,vv1,st,psfc,ptop,eta,qvar,istatus)
      if( istatus.ne.0 )then
          print*,'Cannot read in data of wrf_inout.'
          print*,'filenam: ',filename 
          call exit(1)
      endif

! To write out default system time
      call get_systime(i4time,a9time,istatus)     
      if( istatus.ne.1 )then
          print*,'Cannot get system time.'
          call exit(1)
      endif

! To write out 1_D pressure 
      call get_pres_1d(i4time,kmax,pres_1d,istatus)
      if( istatus.ne.1 )then
          print*,'Cannot get vertical coordinates of LAPS.'
          call exit(1) 
      endif 

! To write out background data
      call get_r_missing_data(r_missing_data,istatus)
      
      call get_modelfg_3d(i4time,'HT ',imax,jmax,kmax,heights_3d,istatus)
      
      call get_modelfg_3d(i4time,'T3 ',imax,jmax,kmax,t_laps_bkg,istatus)
    
      call get_pres_3d(i4time,imax,jmax,kmax,pres_3d,istatus)

      call dryairmass(dam,pdam,imax,jmax,kmax,pres_1d,heights_3d,pres_3d,&
                      t_laps_bkg)

! To unstagger output data  of wrf_inout
! input :: imax,jmax,kmax,imax1,jmax1,kmax1,st,uu1,vv1,qvar
! output(unstagger):: ist(temperature),uu2(u wind),
!                     vv2(v wind),iqvar(specific humidity)  
      call stagger(imax,jmax,kmax,imax1,jmax1,kmax1,st,&
                   uu1,vv1,qvar,ist,uu2,vv2,iqvar)

! To transform mass to laps pressure
! input :: pres_1d(laps p coordinate),imax,jmax,kmax,uu2,vv2,ist,
!       dam(background surface p),iqvar,ptop,eta
! output(mass to lasp coordinate):: uvar1(u wind),vvar1(v wind),
!       tavr1(temperature),qvar1(specific humidity),rh(relative humidity)
      call mass2p(pres_1d,imax,jmax,kmax,uu2,vv2,ist,dam,&
           iqvar,ptop,eta,uvar1,vvar1,tvar1,qvar1,rh)


      call read_static_grid(imax,jmax,'LAT',lat,istatus)
      
      call read_static_grid(imax,jmax,'LON',lon,istatus)
      
      call read_static_grid(imax,jmax,'AVG',topo,istatus)
    
      icen = imax/2+1
      jcen = jmax/2+1

      call get_grid_spacing_actual(lat(icen,jcen),lon(icen,jcen),&
                                   grid_spacing_m,istatus) 

      do j = 1,jmax
         do i = 1,imax
            rk_terrain(i,j) = height_to_zcoord2(topo(i,j),&
       	                      heights_3d,imax,jmax,kmax,i,j,istatus)	    
             klow = max(rk_terrain(i,j),1.) 
      	     khigh = klow + 1
      	     fraclow = float(khigh) - rk_terrain(i,j)
      	     frachigh = 1.0 - fraclow
      	     if(uvar1(i,j,klow) .eq. r_missing_data .or. &
      	        vvar1(i,j,klow) .eq. r_missing_data .or. &
      		uvar1(i,j,khigh) .eq. r_missing_data .or. &
      		vvar1(i,j,khigh) .eq. r_missing_data) then
      	        uanl_sfcitrp(i,j) = r_missing_data
      	        vanl_sfcitrp(i,j) = r_missing_data
      	     else
      	        uanl_sfcitrp(i,j) = uvar1(i,j,klow)*fraclow+&
      		                    uvar1(i,j,khigh)*frachigh
      	        vanl_sfcitrp(i,j) = vvar1(i,j,klow)*fraclow+&
      		                    vvar1(i,j,khigh)*frachigh
      	     endif
      	 enddo
      enddo	 
       
! To write wind to lw3
      call wind_post_process(i4time,'lw3',var_3d,units_3d,&
           comment_3d,uvar1,vvar1,imax,jmax,kmax,3,uanl_sfcitrp,&
      	   vanl_sfcitrp,topo,lat,lon,grid_spacing_m,rk_terrain,&
      	   r_missing_data,l_grid_north_out,istatus)

      !----------------------------------------------------------
! Write lwm for display lw3: YUANFU XIE

!       Write out derived winds file (sfc wind)
        ext_xie = 'lwm'

        var_a_xie(1) = 'SU'
        var_a_xie(2) = 'SV'

        do i = 1,2
            units_a_xie(i) = 'm/s'
            comment_a_xie(i) = 'SFCWIND'
        enddo

        call move(uanl_sfcitrp,out_sfc_2d(1,1,1),imax,jmax)
        call move(vanl_sfcitrp,out_sfc_2d(1,1,2),imax,jmax)

        call put_laps_multi_2d(i4time,ext_xie,var_a_xie &
           ,units_a_xie,comment_a_xie,out_sfc_2d,imax,jmax,2,istatus)
      !----------------------------------------------------------

      call make_fnam_lp(i4time,asc9_tim,istatus)
      comment_2d(1:9) = asc9_tim

! To write temperature to lt1
      call write_temp_anal(i4time,imax,jmax,kmax,tvar1,heights_3d,&
                           output_4d,comment_2d,istatus)

      do k = 1,kmax
         lvl(k) =  pres_1d(k)  * .01 
      enddo   

! To write specific humidity to lq3
      call writefile(i4time,commentline,lvl,qvar1,imax,jmax,kmax,istatus)     

     do i = 1,imax
      do j = 1,jmax
       do k = 1,kmax
     rh1(i,j,k) = make_rh(lvl(k),tvar1(i,j,k)-273.15,qvar1(i,j,k)*1000.,0.0) 
       enddo
      enddo
     enddo

     write(*,*) 'rh1=',rh1(1:10,1,1)
     write(*,*) 'rh1=',rh1(1,1,1:10)    

! To write relative humidity to lh3
!      call lh3_compress(rh,tvar1,i4time,lvl,t_ref,imax,jmax,kmax,istatus)

      END PROGRAM



!------------------> subroutine <------------------------------------------



      SUBROUTINE read_gsi_dim(imax,jmax,kmax,filename,istatus)
       
      implicit none
      include 'netcdf.inc'
      character*125,intent(in) :: filename 
      integer*4,intent(out) :: imax,jmax,kmax
      integer*4 :: kk,n(3),ncid
      integer*4 :: istatus,vid,vtype,vn,&
                   vdims(MAXVDIMS),vnatt
      character :: vname*128

      ncid = NCOPN(filename,NCNOWRIT,istatus)
      if( istatus.ne.0 )return

      vid = NCVID(ncid,'U',istatus)
      if( istatus.ne.0 )return
      
      call NCVINQ(ncid,vid,vname,vtype,vn,vdims,vnatt,istatus)
      if( istatus.ne.0 )return

      do kk=1,3
        call NCDINQ(ncid,vdims(kk),vname,n(kk),istatus)
        if( istatus.ne.0 )return
      enddo

      imax = n(1)
      jmax = n(2)+1
      kmax = n(3)+1
       
      call ncclos(ncid,istatus)
       
      RETURN
      END 



      SUBROUTINE read_gsi_output(ii1,imax,jmax,kmax,imax1,jmax1,kmax1, &
                        filename,uu1,vv1,st,psfc,ptop,eta,qvar,istatus)
       
      implicit none
      include 'netcdf.inc'
      integer*4,intent(in) :: ii1,imax,jmax,kmax,imax1,jmax1,kmax1
      character*125,intent(in) :: filename
      real*4,intent(out) :: uu1(imax,jmax1,kmax1),&
                            vv1(imax1,jmax,kmax1),&
                            st(imax1,jmax1,kmax1),&
                            psfc(imax1,jmax1),ptop,&
                            eta(kmax),qvar(imax1,jmax1,kmax1)
      integer*4 :: ncid
      integer*4 :: i,j,k,ii,jj,kk
      integer*4 :: istatus,vid,vtype,vn,&
                   vdims(MAXVDIMS),vnatt
      integer*4,dimension(2) :: start2,count2
      integer*4,dimension(3) :: start1,count1,n
      integer*4,dimension(4) :: start,count
      character :: vname*128
      character*6 :: vargsi(7)=(/'U     ','V     ','T     ','MUB   ',&
                     'P_TOP ','ZNW   ','QVAPOR'/)
     
      ncid = NCOPN(filename,NCNOWRIT,istatus)
      if( istatus.ne.0 )return
 
      GSI_VARS : DO ii = 1,ii1

          vid = NCVID(ncid,TRIM(vargsi(ii)),istatus)
          if( istatus.ne.0 )return

          call NCVINQ(ncid,vid,vname,vtype,vn,vdims,vnatt,istatus)
          if( istatus.ne.0 )return

        SELECT CASE(ii)

        CASE(1)
          do kk=1,3
            call NCDINQ(ncid,vdims(kk),vname,n(kk),istatus)
            if( istatus.ne.0 )return
          enddo
          start = 1
          count = (/n(1),n(2),n(3),1/)
          call NCVGT(ncid,vid,start,count,uu1,istatus)
          if( istatus.ne.0 )return

        CASE(2)
          n(1:3) = (/imax1,jmax,kmax1/)
          start = 1
          count = (/n(1),n(2),n(3),1/)
          call NCVGT(ncid,vid,start,count,vv1,istatus)
          if( istatus.ne.0 )return

        CASE(3)
          do kk=1,3
            call NCDINQ(ncid,vdims(kk),vname,n(kk),istatus)
            if( istatus.ne.0 )return
          enddo
          start = 1
          count = (/n(1),n(2),n(3),1/)
          call NCVGT(ncid,vid,start,count,st,istatus)
          if( istatus.ne.0 )return
      
        CASE(4)
          do kk=1,2
            call NCDINQ(ncid,vdims(kk),vname,n(kk),istatus)
            if( istatus.ne.0 )return
          enddo
          start1 = 1
          count1 = (/n(1),n(2),1/)
          call NCVGT(ncid,vid,start1,count1,psfc,istatus)
          if( istatus.ne.0 )return

        CASE(5)
          call NCDINQ(ncid,vdims,vname,n(1),istatus)
          if( istatus.ne.0 )return
          start2 = 1
          count2 = (/n(1),1/)
          call NCVGT(ncid,vid,start2,count2,ptop,istatus)
          if( istatus.ne.0 )return

        CASE(6)
          start2 = 1
          count2 = (/kmax,1/)
          call NCVGT(ncid,vid,start2,count2,eta,istatus)
          if( istatus.ne.0 )return

        CASE(7)
          do kk=1,3
            call NCDINQ(ncid,vdims(kk),vname,n(kk),istatus)
            if( istatus.ne.0 )return
          enddo
          start = 1
          count = (/n(1),n(2),n(3),1/)
          call NCVGT(ncid,vid,start,count,qvar,istatus)
          if( istatus.ne.0 )return
  
      END SELECT

      ENDDO GSI_VARS
      
      call NCCLOS(ncid,istatus)
      if( istatus.ne.0 )return
        
      RETURN
      END 



      SUBROUTINE stagger(imax,jmax,kmax,imax1,jmax1,kmax1,st,&
                          uu1,vv1,qvar,ist,uu2,vv2,iqvar)
      
       implicit none
       integer :: i,j,k
       integer*4,intent(in) :: imax,jmax,kmax,imax1,jmax1,kmax1
       real*4,intent(in) ::   st(imax1,jmax1,kmax1),&
                             uu1(imax,jmax1,kmax1),&
                             vv1(imax1,jmax,kmax1),&
                             qvar(imax1,jmax1,kmax1)
       real*4,intent(out) :: ist(imax,jmax,kmax),& 
                             uu2(imax,jmax,kmax),&
                             vv2(imax,jmax,kmax),&
                             iqvar(imax,jmax,kmax)  
       real*4 :: sz(imax,jmax,kmax),sy(imax,jmax,kmax),&
                 sx(imax,jmax,kmax),uz(imax,jmax,kmax),&
                 uy(imax,jmax,kmax),vz(imax,jmax,kmax),&
                 vx(imax,jmax,kmax),qz(imax,jmax,kmax),&
                 qy(imax,jmax,kmax),qx(imax,jmax,kmax) 
 
       do i=1,imax1
          do j=1,jmax1
             do k=2,kmax1
                sz(i,j,k) = 0.5*(st(i,j,k-1)+ st(i,j,k))
             enddo
		  sz(i,j,1) = 1.5*st(i,j,1)-&
                  0.5*st(i,j,2)
		  sz(i,j,kmax) = 1.5*st(i,j,kmax1)-&
                  0.5*st(i,j,kmax-2)
          enddo
       enddo 

       do k=1,kmax
          do i=1,imax1
             do j=2,jmax1
                   sy(i,j,k) = &
                   0.5*(sz(i,j-1,k)+ sz(i,j,k))
             enddo 
             sy(i,1,k) = &
                   1.5*sz(i,1,k)-&
                   0.5*sz(i,2,k)
             sy(i,jmax,k) = &
                   1.5*sz(i,jmax1,k)-&
                   0.5*sz(i,jmax-2,k)
          enddo
       enddo
         
       do k=1,kmax
          do j=1,jmax
             do i=2,imax1
                sx(i,j,k) = &
                0.5*(sy(i-1,j,k)+ sy(i,j,k))
             enddo
             sx(1,j,k) = &
                1.5*sy(1,j,k)-&
                0.5*sy(2,j,k)
             sx(imax,j,k) = &
                1.5*sy(imax1,j,k)-&
                0.5*sy(imax-2,j,k)
          enddo
       enddo

       ist = sx
    
       do i=1,imax
          do j=1,jmax1
             do k=2,kmax1
                uz(i,j,k) = 0.5*(uu1(i,j,k-1)+ uu1(i,j,k))
             enddo
             uz(i,j,1) =  1.5*uu1(i,j,1)-&
                  0.5*uu1(i,j,2)
             uz(i,j,kmax) =  1.5*uu1(i,j,kmax1)-&
                  0.5*uu1(i,j,kmax-2)
          enddo
       enddo

       do k=1,kmax
          do i=1,imax
             do j=2,jmax1
                   uy(i,j,k) = 0.5*(uz(i,j-1,k)+ uz(i,j,k))
             enddo
             uy(i,1,k) = &
                   1.5*uz(i,1,k)-&
                   0.5*uz(i,2,k)
             uy(i,jmax,k) = &
                   1.5*uz(i,jmax1,k)-&
                   0.5*uz(i,jmax-2,k)
          enddo
       enddo

       uu2 = uy
        
      do i=1,imax1
          do j=1,jmax

             do k=2,kmax1
                vz(i,j,k) = 0.5*(vv1(i,j,k-1)+ vv1(i,j,k))
             enddo
                  vz(i,j,1) = 1.5*vv1(i,j,1)-&
                  0.5*vv1(i,j,2)
                  vz(i,j,kmax) = 1.5*vv1(i,j,kmax1)-&
                  0.5*vv1(i,j,kmax-2)
          enddo
       enddo
        
       do k=1,kmax
          do j=1,jmax
             do i=2,imax1
                vx(i,j,k) = &
                0.5*(vz(i-1,j,k)+ vz(i,j,k))
             enddo
             vx(1,j,k) = &
                1.5*vz(1,j,k)-&
                0.5*vz(2,j,k)
             vx(imax,j,k) = &
                1.5*vz(imax1,j,k)-&
                0.5*vz(imax-2,j,k)
          enddo
       enddo

       vv2 = vx

       do i=1,imax1
          do j=1,jmax1
             do k=2,kmax1
                qz(i,j,k) = 0.5*(qvar(i,j,k-1)+ qvar(i,j,k))
             enddo
                  qz(i,j,1) = 1.5*qvar(i,j,1)-&
                  0.5*qvar(i,j,2)
                  qz(i,j,kmax) = 1.5*qvar(i,j,kmax1)-&
                  0.5*qvar(i,j,kmax-2)
          enddo
       enddo

       do k=1,kmax
          do i=1,imax1
             do j=2,jmax1
                   qy(i,j,k) = &
                   0.5*(qz(i,j-1,k)+ qz(i,j,k))
             enddo
             qy(i,1,k) = &
                   1.5*qz(i,1,k)-&
                   0.5*qz(i,2,k)
             qy(i,jmax,k) = &
                   1.5*qz(i,jmax1,k)-&
                   0.5*qz(i,jmax-2,k)
          enddo
       enddo
 
       do k=1,kmax
          do j=1,jmax
             do i=2,imax1
                qx(i,j,k) = &
                0.5*(qy(i-1,j,k)+ qy(i,j,k))
             enddo
             qx(1,j,k) = &
                1.5*qy(1,j,k)-&
                0.5*qy(2,j,k)
             qx(imax,j,k) = &
                1.5*qy(imax1,j,k)-&
                0.5*qy(imax-2,j,k)
          enddo
       enddo

       iqvar = ABS(qx)

       RETURN
       END 


       
      SUBROUTINE mass2p(pres_1d,imax,jmax,kmax,uu2,vv2,ist,&
               dam,iqvar,ptop,eta,uvar1,vvar1,tvar1,qvar1,rh)
 
      implicit none
      integer :: i,j,k,k1,k2,kk,jj,ii
      real*4,parameter :: cp=1004.0, rc=287.0,t0=273.15
      real*4 :: e,pres_1dlog(kmax),pvar(imax,jmax,kmax),&
                pvarlog(imax,jmax,kmax),tt2(imax,jmax,kmax)
      integer*4,intent(in) :: imax,jmax,kmax 
      real*4,intent(in) :: pres_1d(kmax),&
                           uu2(imax,jmax,kmax),&
                           vv2(imax,jmax,kmax),&
                           ist(imax,jmax,kmax),&
                           dam(imax,jmax),iqvar(imax,jmax,kmax),&
                           ptop,eta(kmax)
      real*4,intent(out) :: uvar1(imax,jmax,kmax),&
                            vvar1(imax,jmax,kmax),&
                            tvar1(imax,jmax,kmax),&
                            qvar1(imax,jmax,kmax),&
                            rh(imax,jmax,kmax)

      do k = 1,kmax  
         do j = 1,jmax
	    do i = 1,imax
	       pvar(i,j,k) =ptop+eta(k)*(dam(i,j)-ptop)
               pvarlog(i,j,k)=log(pvar(i,j,k))
            enddo
	 enddo
      enddo

      pres_1dlog(1:kmax)=log(pres_1d(1:kmax))

      do i = 1,imax
        do j = 1,jmax
          do k = 1,kmax-1
            if( pres_1d(k).ge.pvar(i,j,1) ) then
                tt2(i,j,k)=( ist(i,j,1)*(pres_1dlog(k)-pvarlog(i,j,2))  &
                            -ist(i,j,2)*(pres_1dlog(k)-pvarlog(i,j,1)) ) &
                          /( pvarlog(i,j,1)-pvarlog(i,j,2) )
                uvar1(i,j,k)=( uu2(i,j,1)*(pres_1dlog(k)-pvarlog(i,j,2))  &
                            -uu2(i,j,2)*(pres_1dlog(k)-pvarlog(i,j,1)) ) &
                          /( pvarlog(i,j,1)-pvarlog(i,j,2) )
                vvar1(i,j,k)=( vv2(i,j,1)*(pres_1dlog(k)-pvarlog(i,j,2))  &
                            -vv2(i,j,2)*(pres_1dlog(k)-pvarlog(i,j,1)) ) &
                          /( pvarlog(i,j,1)-pvarlog(i,j,2) )
                qvar1(i,j,k)=( iqvar(i,j,1)*(pres_1dlog(k)-pvarlog(i,j,2))  &
                            -iqvar(i,j,2)*(pres_1dlog(k)-pvarlog(i,j,1)) ) &
                          /( pvarlog(i,j,1)-pvarlog(i,j,2) )

            else
                do kk=2,kmax
                   if( pres_1d(k).ge.pvar(i,j,kk) )then
                       k1=kk-1
                       k2=kk
                       go to 10
                   endif
                enddo
 10             continue
                tt2(i,j,k)=( ist(i,j,k1)*(pres_1dlog(k)-pvarlog(i,j,k2))  &
                            -ist(i,j,k2)*(pres_1dlog(k)-pvarlog(i,j,k1)) ) &
                          /( pvarlog(i,j,k1)-pvarlog(i,j,k2) )
                uvar1(i,j,k)=( uu2(i,j,k1)*(pres_1dlog(k)-pvarlog(i,j,k2))  &
                            -uu2(i,j,k2)*(pres_1dlog(k)-pvarlog(i,j,k1)) ) &
                          /( pvarlog(i,j,k1)-pvarlog(i,j,k2) )
                vvar1(i,j,k)=( vv2(i,j,k1)*(pres_1dlog(k)-pvarlog(i,j,k2))  &
                            -vv2(i,j,k2)*(pres_1dlog(k)-pvarlog(i,j,k1)) ) &
                          /( pvarlog(i,j,k1)-pvarlog(i,j,k2) ) 
                qvar1(i,j,k)=( iqvar(i,j,k1)*(pres_1dlog(k)-pvarlog(i,j,k2))  &
                            -iqvar(i,j,k2)*(pres_1dlog(k)-pvarlog(i,j,k1)) ) &
                          /( pvarlog(i,j,k1)-pvarlog(i,j,k2) )             
            endif 
         enddo

         tt2(i,j,kmax)=ist(i,j,kmax)
         uvar1(i,j,kmax)=uu2(i,j,kmax)
         vvar1(i,j,kmax)=vv2(i,j,kmax)
         qvar1(i,j,kmax)=iqvar(i,j,kmax)

        enddo
      enddo

      do k=1,kmax
        tvar1(1:imax,1:jmax,k) = (tt2(1:imax,1:jmax,k)+t0)*&
                               (pres_1d(k)/100000.0)**(rc/cp)
      enddo

      do k=1,kmax
         do j=1,jmax
            do i=1,imax
               e = 6.11*10.0**&
                   (7.45*tt2(i,j,k)/tvar1(i,j,k))
               rh(i,j,k) = qvar1(i,j,k)*pres_1d(k)/(6.22*e)           
            enddo
         enddo
      enddo

      RETURN      
      END 

