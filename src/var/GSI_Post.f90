      PROGRAM gsi2laps

!!       HISTORY:
!!       Creation: LungTsung Cheng    8-2006


      IMPLICIT NONE

!--- decalare variables ---------------------------------------

!//common variables

      INTEGER :: imax,jmax,kmax
      INTEGER :: imax1,jmax1,kmax1  
      INTEGER :: i,j,k,ii,jj,kk,istatus
      INTEGER :: i4time,namelen
      CHARACTER :: filename*125

!// for subroutine read_gsi_output

      INTEGER,PARAMETER :: ii1=8 
      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: uu1,vv1,tt1,qq1
      REAL,ALLOCATABLE,DIMENSION(:,:) :: psfc
      REAL,ALLOCATABLE,DIMENSION(:) :: eta,znu
      REAL :: ptop 

!// for subroutine get_systime

      CHARACTER*9 :: a9time

!// for subroutine get_pres_1d

      REAL,ALLOCATABLE,DIMENSION(:) :: pres_1d

!// for subroutine get_r_missing_data

      REAL :: r_missing_data 

!// for subroutine get_pres_3d

      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: pres_3d 

!// for subroutine get_modelfg_3d

      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: heights_3d 

!// for subroutine get_modelfg_3d 

      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: t_laps_bkg 

!// for subroutine dryairmass 

      REAL,ALLOCATABLE,DIMENSION(:,:) :: dam,pdam 

!// for subroutine unstagger

      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: uu2,vv2,qq2 

!// for subroutine mass2laps

      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: tt2,uvar1,vvar1,&
                                        wvar1,tvar1,qvar1,tvar 
      REAL,PARAMETER :: cp=1004.0, rc=287.0,t0=300.0 !273.15

!// for subroutine read_static_grid

      REAL,ALLOCATABLE,DIMENSION(:,:) :: lat,lon,topo     

!// for subroutine get_grid_spacing_actual

      REAL :: grid_spacing_m
      INTEGER :: icen,jcen

!// for subroutine wind_post_process

      REAL,ALLOCATABLE,DIMENSION(:,:) :: rk_terrain,&
                                           uanl_sfcitrp,vanl_sfcitrp
      CHARACTER*3 :: var_3d(3)=(/'U3','V3','OM'/)
      CHARACTER*4 :: units_3d(3)=(/'M/S ','M/S ','PA/S'/) 
      CHARACTER*125 :: comment_3d(3)=(/'3DWIND',&
                                     '3DWIND','3DWIND'/)
      LOGICAL,PARAMETER :: l_grid_north_out = .true.		
      REAL :: height_to_zcoord2,fraclow,frachigh
      integer :: klow,khigh    
 
!// for subroutine make_fnam_lp

      CHARACTER*9 :: asc9_tim

!// for subroutine write_temp_anal

      REAL,ALLOCATABLE,DIMENSION(:,:,:,:) ::  output_4d
      CHARACTER*10 :: comment_2d

!// for subroutine writefile

      INTEGER,ALLOCATABLE,DIMENSION(:) :: lvl
      CHARACTER*125,PARAMETER :: &
              commentline = 'maps with clouds and surface effects only'

!// for subroutine lh3_compress

      REAL :: t_ref

!//Variables for lwm: YUANFU XIE

      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: out_sfc_2d
      CHARACTER*3 :: ext_xie,var_a_xie(2),units_a_xie(2)
      CHARACTER*7 :: comment_a_xie(2)

!-----> main program <-----------------------------------------------------
		       
! To read path of wrf_inout.

      CALL get_directory('log',filename,namelen)        ! Yuanfu: use LAPS;
      filename = filename(1:namelen-4)//'tmpdirEjetlaps/wrf_inout'

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
      allocate(uu1(imax,jmax1,kmax1))
      allocate(uu2(imax,jmax,kmax))
      allocate(vv1(imax1,jmax,kmax1))
      allocate(vv2(imax,jmax,kmax))
      allocate(tt1(imax1,jmax1,kmax1))
      allocate(tt2(imax,jmax,kmax))
      allocate(uvar1(imax,jmax,kmax))
      allocate(vvar1(imax,jmax,kmax))
      allocate(wvar1(imax,jmax,kmax))
      allocate(tvar(imax,jmax,kmax))
      allocate(tvar1(imax,jmax,kmax))
      allocate(qq1(imax1,jmax1,kmax1))
      allocate(qq2(imax,jmax,kmax))
      allocate(qvar1(imax,jmax,kmax))
      allocate(heights_3d(imax,jmax,kmax))
      allocate(out_sfc_2d(imax,jmax,2))
      allocate(pres_3d(imax,jmax,kmax))
      allocate(t_laps_bkg(imax,jmax,kmax))
      allocate(psfc(imax1,jmax1))
      allocate(uanl_sfcitrp(imax,jmax))
      allocate(vanl_sfcitrp(imax,jmax))
      allocate(topo(imax,jmax))
      allocate(lat(imax,jmax))
      allocate(lon(imax,jmax))
      allocate(rk_terrain(imax,jmax))
      allocate(dam(imax,jmax))
      allocate(pdam(imax,jmax))
      allocate(eta(kmax))
      allocate(znu(kmax1))
      allocate(pres_1d(kmax))
      allocate(lvl(kmax))

! To write out variables of wrf_inout(after GSI analysis)
! input :: ii1,imax,jmax,kmax,imax1,jmax1,kmax1,filename
! output:: uu1(u wind),vv1(v wind),tt1(temperature),
!          psfc(surface pressure),ptop(top pressure),
!          eta(on constant p levels),znu(eta on mass levels)),&
!          qq1(specific humidity)
 
      call read_gsi_output(ii1,imax,jmax,kmax,imax1,jmax1,kmax1, &
                   filename,uu1,vv1,tt1,psfc,ptop,eta,znu,qq1,istatus)
      if( istatus.ne.0 )then
          print*,'Cannot read in data of wrf_inout.'
          print*,'filenam: ',filename 
          call exit(1)
      endif

! To write out default system time
!      call get_systime(i4time,a9time,istatus)    
!      if( istatus.ne.1 )then
!          print*,'Cannot get system time.'
!          call exit(1)
!      endif
       CALL get_directory('log',filename,namelen)
       filename = filename(1:namelen)//'i4time.txt'
       ! open(10,file='i4time.txt')
       open(10,file=filename(1:namelen+10))
       read(10,*) i4time
       close(10)

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
! input :: imax,jmax,kmax,imax1,jmax1,kmax1,tt1,uu1,vv1,qq1
!          ptop,psfc,znu
! output(unstagger):: tt2(temperature),uu2(u wind),
!                     vv2(v wind),qq2(specific humidity)  

      call unstagger(imax,jmax,kmax,imax1,jmax1,kmax1,tt1,&
           uu1,vv1,qq1,tt2,uu2,vv2,qq2,ptop,dam,znu)

! XIE: mass to laps transfer:

      CALL Mass2LAPS(tvar,imax,jmax,kmax,pres_1d,dam,eta,4,0,tt2)
      do k=1,kmax
        tvar1(1:imax,1:jmax,k) = (tvar(1:imax,1:jmax,k)+t0)*&
                               (pres_1d(k)/100000.0)**(rc/cp)
      enddo
      CALL Mass2LAPS(uvar1,imax,jmax,kmax,pres_1d,dam,eta,4,0,uu2)
      CALL Mass2LAPS(vvar1,imax,jmax,kmax,pres_1d,dam,eta,4,0,vv2)
      CALL Mass2LAPS(qvar1,imax,jmax,kmax,pres_1d,dam,eta,2,0,qq2)

! for background
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

! old-version-wind_post_process
!      call wind_post_process(i4time,'lw3',var_3d,units_3d,&
!           comment_3d,uvar1,vvar1,imax,jmax,kmax,3,uanl_sfcitrp,&
!      	   vanl_sfcitrp,topo,lat,lon,grid_spacing_m,rk_terrain,&
!      	   r_missing_data,l_grid_north_out,istatus)

! wind_post_process for new version LAPS
! modify : LungTsung Cheng  12-2007 

        call wind_post_process(i4time                         &   ! I 
                             ,uvar1,vvar1                     &   ! I
                             ,wvar1                           &   ! O
                             ,imax,jmax,kmax,3                &   ! I
                             ,heights_3d                      &   ! I
                             ,uanl_sfcitrp,vanl_sfcitrp       &   ! I
                             ,topo,lat,lon,grid_spacing_m     &   ! I
                             ,r_missing_data,l_grid_north_out &   ! I
                             ,istatus)
        
        print*,"wvar1=",wvar1(1,1,1:10)

        call write_wind_output(i4time,'lw3',var_3d            &   ! I
                             ,uvar1,vvar1                     &   ! I
                             ,wvar1                           &   ! I
                             ,uanl_sfcitrp,vanl_sfcitrp       &   ! I
                             ,imax,jmax,kmax,3                &   ! I
                             ,r_missing_data                  &   ! I
                             ,istatus)


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

!
      do k = 1,kmax
         lvl(k) =  pres_1d(k)  * .01 
      enddo   

! To write specific humidity to lq3
    
      call writefile(i4time,commentline,lvl,qvar1,imax,jmax,kmax,istatus)     

!
      t_ref = 0.0   

! To write relative humidity to lh3

      call lh3_compress(qvar1,tvar1,i4time,lvl,t_ref,imax,jmax,kmax,1,istatus)

      END PROGRAM



!------------------> subroutine <------------------------------------------



      SUBROUTINE read_gsi_dim(imax,jmax,kmax,filename,istatus)
       
      implicit none
      include 'netcdf.inc'
      character*125,intent(in) :: filename 
      integer,intent(out) :: imax,jmax,kmax
      integer :: kk,n(3),ncid
      integer :: istatus,vid,vtype,vn,&
                   vdims(MAXVDIMS),vnatt
      character(MAXNCNAM) :: vname

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
                       filename,uu1,vv1,tt1,psfc,ptop,eta,znu,qq1,istatus)
       
      implicit none
      include 'netcdf.inc'
      integer,intent(in) :: ii1,imax,jmax,kmax,imax1,jmax1,kmax1
      character*125,intent(in) :: filename
      real,intent(out) :: uu1(imax,jmax1,kmax1),&
                            vv1(imax1,jmax,kmax1),&
                            tt1(imax1,jmax1,kmax1),&
                            psfc(imax1,jmax1),ptop,&
                            eta(kmax),znu(kmax1),&
                            qq1(imax1,jmax1,kmax1)
      integer :: ncid
      integer :: i,j,k,ii,jj,kk
      integer :: istatus,vid,vtype,vn,&
                   vdims(MAXVDIMS),vnatt
      integer,dimension(2) :: start2,count2
      integer,dimension(3) :: start1,count1,n
      integer,dimension(4) :: start,count
      character(MAXNCNAM) :: vname
      character*6 :: vargsi(8)=(/'U     ','V     ','T     ','MUB   ',&
                     'P_TOP ','ZNW   ','ZNU   ','QVAPOR'/)
     
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
          call NCVGT(ncid,vid,start,count,tt1,istatus)
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
          start2 = 1
          count2 = (/kmax1,1/)
          call NCVGT(ncid,vid,start2,count2,znu,istatus)
          if( istatus.ne.0 )return
        
        CASE(8)
          do kk=1,3
            call NCDINQ(ncid,vdims(kk),vname,n(kk),istatus)
            if( istatus.ne.0 )return
          enddo
          start = 1
          count = (/n(1),n(2),n(3),1/)
          call NCVGT(ncid,vid,start,count,qq1,istatus)
          if( istatus.ne.0 )return
  
      END SELECT

      ENDDO GSI_VARS
      
      call NCCLOS(ncid,istatus)
      if( istatus.ne.0 )return
        
      RETURN
      END 



      SUBROUTINE unstagger(imax,jmax,kmax,imax1,jmax1,kmax1,tt1,&
                 uu1,vv1,qq1,tt2,uu2,vv2,qq2,ptop,psfc,znu)
      
       implicit none
       integer :: i,j,k
       integer,intent(in) :: imax,jmax,kmax,imax1,jmax1,kmax1
       real,intent(in) ::   tt1(imax1,jmax1,kmax1),&
                             uu1(imax,jmax1,kmax1),&
                             vv1(imax1,jmax,kmax1),&
                             qq1(imax1,jmax1,kmax1),&
                             ptop,psfc(imax1,jmax1),znu(kmax1)
       real,intent(out) :: tt2(imax,jmax,kmax),& 
                             uu2(imax,jmax,kmax),&
                             vv2(imax,jmax,kmax),&
                             qq2(imax,jmax,kmax)  
       real :: sz(imax1,jmax1,kmax),uz(imax,jmax1,kmax),&
                 vz(imax1,jmax,kmax),qz(imax1,jmax1,kmax),&
                 qout(imax,jmax,kmax)
        

! get unstagger temperature variable tt2
 
       call UnStaggerZ(tt1,imax1,jmax1,kmax,ptop,psfc,znu,4,0,sz)
       call UntaggerXY_3D(sz,imax1,jmax1,kmax,imax,jmax,tt2)

! get unstagger wind u component uu2 

       call UnStaggerZ(uu1,imax,jmax1,kmax,ptop,psfc,znu,4,0,uz)
       call UntaggerXY_3D(uz,imax,jmax1,kmax,imax,jmax,uu2)

! get unstagger wind v component vv2 
       
       call UnStaggerZ(vv1,imax1,jmax,kmax,ptop,psfc,znu,4,0,vz)
       call UntaggerXY_3D(vz,imax1,jmax,kmax,imax,jmax,vv2)

! get unstagger specific humidity qq2 

       call UnStaggerZ(qq1,imax1,jmax1,kmax,ptop,psfc,znu,4,0,qz)
       call UntaggerXY_3D(qz,imax1,jmax1,kmax,imax,jmax,qout)
       qq2 = ABS(qout)

       RETURN
       END 

