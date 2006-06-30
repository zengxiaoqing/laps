      PROGRAM gsi2laps

      IMPLICIT NONE
      integer*4,parameter :: nxx=300,nyy=300,nzz=30
      INTEGER*4,PARAMETER :: n3d=3,ii1=7
      integer*4 :: imax,jmax,kmax,imax1,jmax1,kmax1
      INTEGER*4 :: i,j,k,ii,jj,kk,ij,ier,ncid
      integer*4 :: ijklow,ijkhigh,icenjcen
      CHARACTER :: filename*125
      REAL*4 :: uu1(nxx*nyy*nzz),&
                uu2(nxx*nyy*nzz),&
                vv1(nxx*nyy*nzz),&
                vv2(nxx*nyy*nzz),&
                st(nxx*nyy*nzz),&
                ist(nxx*nyy*nzz),&
                psfc1(nxx*nyy),&
                ipsfc(nxx*nyy),&
                eta(nzz),lapsp(nzz),&
                uvar1(nxx*nyy*nzz),&
                vvar1(nxx*nyy*nzz),&
                tvar1(nxx*nyy*nzz),&
                qvar(nxx*nyy*nzz),&
                iqvar(nxx*nyy*nzz),qvar1(nxx*nyy*nzz)
      REAL*4 :: uanl_sfcitrp(nxx*nyy),vanl_sfcitrp(nxx*nyy),&
                topo(nxx*nyy),lat(nxx*nyy),lon(nxx*nyy),&    
                rk_terrain(nxx*nyy),heights_3d(nxx*nyy*nzz)
      real*4 :: ptop,grid_spacing_m,r_missing_data,height_to_zcoord2,&
                fraclow,frachigh,output_4d(nxx,nyy,nzz,2)
      INTEGER*4 :: gsi_i4time,i4time_asc_gg,i4time,namelen,istatus,&
                   klow,khigh,icen,jcen,istatus_ht,i4time_nearest
      CHARACTER :: times*19,atime*24,a9time*9
      CHARACTER*3 :: var_3d(n3d)=(/'U3','V3','OM'/)
      character*4 :: units_3d(n3d)=(/'M/S ','M/S ','PA/S'/) 
      character*9 :: asc9_tim
      CHARACTER*10 :: units_2d
      CHARACTER*125 :: comment_3d(n3d)=(/'3DWIND',&
                       '3DWIND','3DWIND'/),comment_2d
      LOGICAL,PARAMETER :: l_grid_north_out = .true.		
      CHARACTER*3,PARAMETER :: EXT = 'lw3'

      ! Variables for lwm: YUANFU XIE
      CHARACTER*3 :: ext_xie,var_a_xie(2),units_a_xie(2)
      CHARACTER*7 :: comment_a_xie(2)
      REAL*4 :: out_sfc_2d(nxx,nyy,2),uma(nxx,nyy,nzz),pi

      real*4 :: dam(nxx*nyy),pdam(nxx*nyy),vdf1(nxx*nyy*nzz),&
                pres_3d(nxx*nyy*nzz),t_laps_bkg(nxx*nyy*nzz),&
                massp(153,141,21),dam1(153,141)
		       
! To read path of wrf_inout.

      CALL get_directory('log',filename,namelen)        ! Yuanfu: use LAPS;
      filename = filename(1:namelen-4)//'tmpdirEjetlaps/wrf_inout'
      !filename = 'wrf_inout.nc'

! To read gsi variables
      call read_gsi_dim(imax,jmax,kmax,ncid,filename,ier)
      if( ier.ne.0 )then
          print*,'Cannot read dimensions of wrf_inout.'
          print*,'filenam: ',filename 
          call exit(1)
      endif
      imax1=imax-1
      jmax1=jmax-1
      kmax1=kmax-1
      call read_gsi_output(ii1,imax,jmax,kmax,imax1,jmax1,kmax1, &
                       ncid,uu1,vv1,st,psfc1,ptop,eta,qvar,ier)
      if( ier.ne.0 )then
          print*,'Cannot read in data of wrf_inout.'
          print*,'filenam: ',filename 
          call exit(1)
      endif

      call get_systime(i4time,a9time,istatus)     
      if( istatus.ne.1 )then
          print*,'Cannot get system time.'
          call exit(1)
      endif

      call get_pres_1d(i4time,kmax,lapsp,istatus)
      if( istatus.ne.1 )then
          print*,'Cannot get vertical coordinates of LAPS.'
          call exit(1) 
      endif 

      call get_r_missing_data(r_missing_data,istatus)
      
      call get_modelfg_3d(i4time,'HT ',imax,jmax,kmax,heights_3d,istatus)
    
      call get_pres_3d(i4time,imax,jmax,kmax,pres_3d,istatus)

      call get_modelfg_3d(i4time,'T3 ',imax,jmax,kmax,t_laps_bkg,istatus)

      call dryairmass(dam,pdam,imax,jmax,kmax,lapsp,heights_3d,pres_3d,&
                      t_laps_bkg)

          ij = 0              
          do j=1,jmax
          do i=1,imax
             ij=ij+1
             dam1(i,j)=dam(ij)
          enddo
          enddo

      !call comdif(imax,jmax,kmax,eta,lapsp,dam,uu1,vdf1) //compare before_GSI and after_GSI              
      !uu1 = vdf1

      do k = 1,kmax
        do j = 1,jmax
         do i = 1,imax
      massp(i,j,k) = eta(k)*(dam1(i,j)-lapsp(kmax))+lapsp(kmax)
         enddo
        enddo
      enddo 

      call stagger(imax,jmax,kmax,imax1,jmax1,kmax1,st,ist,psfc1,ipsfc,&
                   uu1,vv1,uu2,vv2,massp,qvar,iqvar)

      ipsfc = dam

!// uma is to substitute uu2
!       pi = 4.0*atan(1.0)  
!	ij = 0
!	do k = 1,kmax
!        do j = 1,jmax
!         do i = 1,imax
!           uma(i,j,k) = sin(2.0*pi*float(i-1)/float(imax-1))*&
!                    sin(2.0*pi*float(j-1)/float(jmax-1))*&
!                    sin(2.0*pi*(massp(i,j,k)-110000.0)&
!                               /(10000.0-110000.0))
!	   ij = ij+1
!           uu2(ij) = uma(i,j,k)
!         enddo
!        enddo
!       enddo

      CALL mass2p(lapsp,imax,jmax,kmax,uu2,vv2,ist,ipsfc,&
           ptop,eta,uvar1,vvar1,tvar1,iqvar,qvar1)

      call read_static_grid(imax,jmax,'LAT',lat,istatus)
      
      call read_static_grid(imax,jmax,'LON',lon,istatus)
      
      call read_static_grid(imax,jmax,'AVG',topo,istatus)
    
      icen = imax/2+1
      jcen = jmax/2+1
      icenjcen=(jcen-1)*imax+icen
      call get_grid_spacing_actual(lat(icenjcen),lon(icenjcen),&
                                   grid_spacing_m,istatus) 

      ij=0
      do j = 1,jmax
         do i = 1,imax
            ij=ij+1
            rk_terrain(ij) = height_to_zcoord2(topo(ij),&
       	                      heights_3d,imax,jmax,kmax,i,j,istatus)	    
            
             klow = max(rk_terrain(ij),1.) 
      	     khigh = klow + 1
      	     fraclow = float(khigh) - rk_terrain(ij)
      	     frachigh = 1.0 - fraclow

             ijklow=(klow-1)*imax*jmax+ij
             ijkhigh=(khigh-1)*imax*jmax+ij
	     
      	     if(uvar1(ijklow) .eq. r_missing_data .or. &
      	        vvar1(ijklow) .eq. r_missing_data .or. &
      		uvar1(ijkhigh) .eq. r_missing_data .or. &
      		vvar1(ijkhigh) .eq. r_missing_data) then
      	        uanl_sfcitrp(ij) = r_missing_data
      	        vanl_sfcitrp(ij) = r_missing_data
      	     else
      	        uanl_sfcitrp(ij) = uvar1(ijklow)*fraclow+&
      		                    uvar1(ijkhigh)*frachigh
      	        vanl_sfcitrp(ij) = vvar1(ijklow)*fraclow+&
      		                    vvar1(ijkhigh)*frachigh
      	     endif
	    
      	 enddo
      enddo	 
       
      call wind_post_process(i4time,EXT,var_3d,units_3d,&
           comment_3d,uvar1,vvar1,imax,jmax,kmax,n3d,uanl_sfcitrp,&
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

      call write_temp_anal(i4time,imax,jmax,kmax,tvar1,heights_3d,&
                           output_4d,comment_2d,istatus)
      
      END PROGRAM


      subroutine read_gsi_dim(imax,jmax,kmax,ncid,filename,ier)
       
      implicit none
      include 'netcdf.inc'
      integer*4,intent(out) :: imax,jmax,kmax,ncid
      character*125,intent(in) :: filename 
      integer*4 :: kk,n(3)
      integer*4 :: ier,vid,vtype,vn,&
                   vdims(MAXVDIMS),vnatt
      character :: vname*10

      ncid = NCOPN(filename,NCNOWRIT,ier)
      if( ier.ne.0 )return

      vid = NCVID(ncid,'U',ier)
      if( ier.ne.0 )return
      CALL NCVINQ(ncid,vid,vname,vtype,vn,vdims,vnatt,ier)
      if( ier.ne.0 )return

      DO kk=1,3
        CALL NCDINQ(ncid,vdims(kk),vname,n(kk),ier)
        if( ier.ne.0 )return
      ENDDO
      imax = n(1)
      jmax = n(2)+1
      kmax = n(3)+1
       
  !   call ncclos(ncid,rcode)
       
       return
       end 


      subroutine read_gsi_output(ii1,imax,jmax,kmax,imax1,jmax1,kmax1, &
                          ncid,uu1,vv1,st,psfc1,ptop,eta,qvar,ier)
       
      implicit none
      include 'netcdf.inc'
      integer*4,intent(in) :: ii1,imax,jmax,kmax,imax1,jmax1,kmax1,ncid
      integer*4 :: i,j,k,ii,jj,kk
      integer*4 :: ier,vid,vtype,vn,&
                   vdims(MAXVDIMS),vnatt
      integer*4,dimension(2) :: start2,count2
      integer*4,dimension(3) :: start1,count1,n
      integer*4,dimension(4) :: start,count
      integer*4 :: i4time_asc_gg,istatus
      character :: vname*10
      CHARACTER*6 :: vargsi(7)=(/'U     ','V     ','T     ','MUB   ',&
                     'P_TOP ','ZNW   ','QVAPOR'/)
      REAL*4,intent(out) :: uu1(imax,jmax1,kmax1),&
                vv1(imax1,jmax,kmax1),&
                st(imax1,jmax1,kmax1),&
                psfc1(imax1,jmax1),&
                ptop,eta(kmax),qvar(imax1,jmax1,kmax1)
     
  !   ncid = NCOPN(filename,NCNOWRIT,ier)
 
      GSI_VARS : DO ii = 1,ii1

      vid = NCVID(ncid,TRIM(vargsi(ii)),ier)
      if( ier.ne.0 )return

      CALL NCVINQ(ncid,vid,vname,vtype,vn,vdims,vnatt,ier)
      if( ier.ne.0 )return

      SELECT CASE(ii)

      CASE(1)
      DO kk=1,3
        CALL NCDINQ(ncid,vdims(kk),vname,n(kk),ier)
        if( ier.ne.0 )return
      ENDDO
      start = 1
      count = (/n(1),n(2),n(3),1/)
      CALL NCVGT(ncid,vid,start,count,uu1,ier)
      if( ier.ne.0 )return

      CASE(2)
      n(1:3) = (/imax1,jmax,kmax1/)
      start = 1
      count = (/n(1),n(2),n(3),1/)
      CALL NCVGT(ncid,vid,start,count,vv1,ier)
      if( ier.ne.0 )return

      CASE(3)
      DO kk=1,3
        CALL NCDINQ(ncid,vdims(kk),vname,n(kk),ier)
        if( ier.ne.0 )return
      ENDDO
      start = 1
      count = (/n(1),n(2),n(3),1/)
      CALL NCVGT(ncid,vid,start,count,st,ier)
      if( ier.ne.0 )return
      
      CASE(4)
      DO kk=1,2
         CALL NCDINQ(ncid,vdims(kk),vname,n(kk),ier)
         if( ier.ne.0 )return
      ENDDO
      start1 = 1
      count1 = (/n(1),n(2),1/)
      CALL NCVGT(ncid,vid,start1,count1,psfc1,ier)
      if( ier.ne.0 )return

      CASE(5)
      CALL NCDINQ(ncid,vdims,vname,n(1),ier)
      if( ier.ne.0 )return
      start2 = 1
      count2 = (/n(1),1/)
      CALL NCVGT(ncid,vid,start2,count2,ptop,ier)
      if( ier.ne.0 )return

      CASE(6)
      start2 = 1
      count2 = (/kmax,1/)
      CALL NCVGT(ncid,vid,start2,count2,eta,ier)
      if( ier.ne.0 )return

      CASE(7)
      DO kk=1,3
        CALL NCDINQ(ncid,vdims(kk),vname,n(kk),ier)
        if( ier.ne.0 )return
      ENDDO
      start = 1
      count = (/n(1),n(2),n(3),1/)
      CALL NCVGT(ncid,vid,start,count,qvar,ier)
      if( ier.ne.0 )return
     
      do k = 1,20
         write(*,*) 'qvar=',qvar(4,4,k)
      enddo 
  
      END SELECT

      ENDDO GSI_VARS
      
      call ncclos(ncid,ier)
      if( ier.ne.0 )return
        
      return
      end 

       subroutine stagger(imax,jmax,kmax,imax1,jmax1,kmax1,st,ist,psfc,&
                          ipsfc1,uu1,vv1,uu2,vv2,massp,qvar,iqvar)
      
       implicit none
       integer :: i,j,k
       integer*4,intent(in) :: imax,jmax,kmax,imax1,jmax1,kmax1
       real*4 :: sz(imax,jmax,kmax),sy(imax,jmax,kmax),&
                 sx(imax,jmax,kmax),psx(imax,jmax),&
		 psy(imax,jmax),uz(imax,jmax,kmax),&
                 uy(imax,jmax,kmax),vz(imax,jmax,kmax),&
                 vx(imax,jmax,kmax),qz(imax,jmax,kmax),&
                 qy(imax,jmax,kmax),qx(imax,jmax,kmax) 

       real*4,intent(in) :: st(imax1,jmax1,kmax1),&
			    psfc(imax1,jmax1),&
                            uu1(imax,jmax1,kmax1),&
                            vv1(imax1,jmax,kmax1),&
                            massp(imax,jmax,kmax),&
                            qvar(imax1,jmax1,kmax1)
       real*4,intent(out) :: ist(imax,jmax,kmax),& 
                             ipsfc1(imax,jmax),&
                             uu2(imax,jmax,kmax),&
                             vv2(imax,jmax,kmax),&
                             iqvar(imax,jmax,kmax)  
       real*4 :: unstag_diffu(imax,jmax,kmax),unstagu(imax,jmax,kmax),&
                 err,ftn,pi,maxdata
       integer*4 :: ii,jj,kk
 
       do i = 1,imax1
          do j = 2,jmax1
	     psy(i,j) = &
	     0.5*(psfc(i,j-1)+psfc(i,j))
	  enddo
	  psy(i,1) = &
	     1.5*psfc(i,1)-&
	     0.5*psfc(i,2)
	  psy(i,jmax) = &
	     1.5*psfc(i,jmax1)-&
	     0.5*psfc(i,jmax-2)
       enddo    	  

       do j = 1,jmax
          do i = 2,imax1
	     psx(i,j) = &
	     0.5*(psy(i-1,j)+psy(i,j))
	  enddo
          psx(1,j) = &
	     1.5*psy(1,j)-&
	     0.5*psy(2,j)
          psx(imax,j) = &
	     1.5*psy(imax1,j)-&
	     0.5*psy(imax-2,j)
       enddo	  

       ipsfc1 = psx

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

!// compare unstagger_data and real_function_data
!      write(*,*) 'unstag_mass_pressure==',massp(4,4,1:21) 
!      pi = 4.0*ATAN(1.0)
!
!       do kk = 1,kmax
!        do jj = 1,jmax
!         do ii = 1,imax
!           unstagu(ii,jj,kk) = sin(2.0*pi*float(ii-1)/float(imax-1))*&
!                    sin(2.0*pi*float(jj-1)/float(jmax-1))*&
!                    sin(2.0*pi*(massp(ii,jj,kk)-110000.0)&
!                               /(massp(ii,jj,kmax)-110000.0))
!         enddo
!        enddo
!       enddo
!
!      unstag_diffu(1:imax,1:jmax,1:kmax) = uu2(1:imax,1:jmax,1:kmax)&
!                                       -unstagu(1:imax,1:jmax,1:kmax)
!
!  do kk = 1,kmax
!  write(*,*) 'k=',kk
!  write(*,*) 'stag2unstag_mass_u_diff(4,4,k)(%)=',unstag_diffu(4,4,kk)*100.0
!  write(*,*) 'stag2unstag_mass_u=',uu2(4,4,kk)
!  write(*,*) 'real_stag2unstag_mass_unstag_u=',unstagu(4,4,kk)
!  write(*,*) '                                                         '
!  enddo
!       err = 0.0
!       ftn = 0.0
!       maxdata = 0.0
!       DO kk=1,kmax
!         DO jj=1,jmax
!           DO ii=1,imax
!             err = err+unstag_diffu(ii,jj,kk)**2
!             ftn = ftn+unstagu(ii,jj,kk)**2
!             maxdata = max(maxdata,ABS(unstag_diffu(ii,jj,kk)))
!           ENDDO
!         ENDDO
!       ENDDO
!        PRINT*,'Relative error of unstagger conversion: ',sqrt(err/ftn)
!
!       DO kk=1,kmax
!         DO jj=1,jmax
!           DO ii=1,imax
!              if(ABS(unstag_diffu(ii,jj,kk)) == maxdata) then
!                 print*,'i=',ii,'j=',jj,'k=',kk,'maxdata=',maxdata
!              endif
!           ENDDO
!         ENDDO
!       ENDDO
        
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

       iqvar = qx

 
       return
       end 

       
      SUBROUTINE mass2p(vlp,imax,jmax,kmax,uu1,vv1,tt1,&
                 psfc1,ptop1,eta,uu2,vv2,tt3,iqvar,qvar1)
      
      IMPLICIT NONE
      
      INTEGER :: i,j,k,k1,k2,kk,jj,ii
      real*4,parameter :: cp=1004.0, rc=287.0,t0=273.15
      INTEGER*4,INTENT(IN) :: imax,jmax,kmax 
      REAL*4 :: a,b,vlplog(kmax),pvar(imax,jmax,kmax),&
                pvarlog(imax,jmax,kmax),tt2(imax,jmax,kmax)
      REAL*4,INTENT(OUT) :: uu2(imax,jmax,kmax),&
                            vv2(imax,jmax,kmax),&
                            tt3(imax,jmax,kmax),&
                            qvar1(imax,jmax,kmax)
      REAL*4,INTENT(IN) :: uu1(imax,jmax,kmax),&
                vv1(imax,jmax,kmax),&
                tt1(imax,jmax,kmax),&
                psfc1(imax,jmax),iqvar(imax,jmax,kmax)
      REAL*4,INTENT(IN) :: ptop1,&
                eta(kmax),vlp(kmax)
      real*4 :: m2p(imax,jmax,kmax),m2p_diffu(imax,jmax,kmax),&
                pi,err,ftn,maxdata

      DO k = 1,kmax  
         DO j = 1,jmax
	    DO i = 1,imax
	       pvar(i,j,k) =ptop1+eta(k)*(psfc1(i,j)-ptop1)
               pvarlog(i,j,k)=log(pvar(i,j,k))
            ENDDO
	 ENDDO
      ENDDO
      do k = 1,kmax
         vlplog(k)=log(vlp(k))
      enddo

      do i = 1,imax
      do j = 1,jmax
         DO k = 1,kmax-1
            if( vlp(k).ge.pvar(i,j,1) ) then
                tt2(i,j,k)=( tt1(i,j,1)*(vlplog(k)-pvarlog(i,j,2))  &
                            -tt1(i,j,2)*(vlplog(k)-pvarlog(i,j,1)) ) &
                          /( pvarlog(i,j,1)-pvarlog(i,j,2) )
                uu2(i,j,k)=( uu1(i,j,1)*(vlplog(k)-pvarlog(i,j,2))  &
                            -uu1(i,j,2)*(vlplog(k)-pvarlog(i,j,1)) ) &
                          /( pvarlog(i,j,1)-pvarlog(i,j,2) )
                vv2(i,j,k)=( vv1(i,j,1)*(vlplog(k)-pvarlog(i,j,2))  &
                            -vv1(i,j,2)*(vlplog(k)-pvarlog(i,j,1)) ) &
                          /( pvarlog(i,j,1)-pvarlog(i,j,2) )
                qvar1(i,j,k)=( iqvar(i,j,1)*(vlplog(k)-pvarlog(i,j,2))  &
                            -iqvar(i,j,2)*(vlplog(k)-pvarlog(i,j,1)) ) &
                          /( pvarlog(i,j,1)-pvarlog(i,j,2) )

            else
                do kk=2,kmax
                   if( vlp(k).ge.pvar(i,j,kk) )then
                       k1=kk-1
                       k2=kk
                       go to 10
                   endif
                enddo
 10             continue
                tt2(i,j,k)=( tt1(i,j,k1)*(vlplog(k)-pvarlog(i,j,k2))  &
                            -tt1(i,j,k2)*(vlplog(k)-pvarlog(i,j,k1)) ) &
                          /( pvarlog(i,j,k1)-pvarlog(i,j,k2) )
                uu2(i,j,k)=( uu1(i,j,k1)*(vlplog(k)-pvarlog(i,j,k2))  &
                            -uu1(i,j,k2)*(vlplog(k)-pvarlog(i,j,k1)) ) &
                          /( pvarlog(i,j,k1)-pvarlog(i,j,k2) )
                vv2(i,j,k)=( vv1(i,j,k1)*(vlplog(k)-pvarlog(i,j,k2))  &
                            -vv1(i,j,k2)*(vlplog(k)-pvarlog(i,j,k1)) ) &
                          /( pvarlog(i,j,k1)-pvarlog(i,j,k2) ) 
                qvar1(i,j,k)=( iqvar(i,j,k1)*(vlplog(k)-pvarlog(i,j,k2))  &
                            -iqvar(i,j,k2)*(vlplog(k)-pvarlog(i,j,k1)) ) &
                          /( pvarlog(i,j,k1)-pvarlog(i,j,k2) )             
            endif 
         ENDDO

         tt2(i,j,kmax)=tt1(i,j,kmax)
         uu2(i,j,kmax)=uu1(i,j,kmax)
         vv2(i,j,kmax)=vv1(i,j,kmax)
         qvar1(i,j,kmax)=iqvar(i,j,kmax)

      enddo
      enddo

      do k=1,kmax
        tt3(1:imax,1:jmax,k) = (tt2(1:imax,1:jmax,k)+t0)*&
       (vlp(k)/100000.0)**(rc/cp)
      enddo

!      pi = 4.0*ATAN(1.0)
!
!      do k = 1,kmax
!        do j = 1,jmax
!         do i = 1,imax
!           m2p(i,j,k) = sin(2.0*pi*float(i-1)/float(imax-1))*&
!                    sin(2.0*pi*float(j-1)/float(jmax-1))*&
!                    sin(2.0*pi*float(k-1)/float(kmax-1))
!         enddo
!        enddo
!       enddo
!
!       m2p_diffu(1:imax,1:jmax,1:kmax) = uu2(1:imax,1:jmax,1:kmax)-&
!                                        m2p(1:imax,1:jmax,1:kmax)
!       do k = 1,kmax
!       if (m2p(4,4,k) .ne. 0.0) then
!       write(*,*) ' i=4 ',' j=4 ',' k= ',k,'laps_pressure=',vlp(k)
!       write(*,*) 'abs_laps_u_diff(%)=',ABS(m2p_diffu(4,4,k)/m2p(4,4,k))*100.0
!       write(*,*) 'mass2laps_u=',uu2(4,4,k)
!       write(*,*) 'real_laps_u=',m2p(4,4,k)
!       write(*,*) '                                                '
!       else
!       write(*,*) ' i=4 ',' j=4 ',' k= ',k,'laps_pressure=',vlp(k)
!       write(*,*) 'laps_u_diff(%)=',m2p_diffu(4,4,k)*100.0
!       write(*,*) 'mass2laps_u=',uu2(4,4,k)
!       write(*,*) 'real_laps_u=',m2p(4,4,k)
!       write(*,*) '                                                '
!       endif
!       enddo
!
!       err = 0.0
!       ftn = 0.0
!       maxdata = 0.0
!       DO kk=1,kmax
!         DO jj=1,jmax
!           DO ii=1,imax
!             err = err+m2p_diffu(ii,jj,kk)**2
!             ftn = ftn+m2p(ii,jj,kk)**2
!             maxdata = max(maxdata,ABS(m2p_diffu(ii,jj,kk)))
!           ENDDO
!         ENDDO
!       ENDDO
!        PRINT*,'Relative error of mass to pressure conversion: ',sqrt(err/ftn)
!       
!      DO kk=10,kmax
!         DO jj=1,jmax
!           DO ii=1,imax
!              if(ABS(m2p_diffu(ii,jj,kk)) == maxdata) then
!                 print*,'i=',ii,'j=',jj,'k=',kk,'maxdata=',maxdata
!              endif
!           ENDDO
!         ENDDO
!       ENDDO

      RETURN      
      END 

subroutine comdif(imax,jmax,kmax,znw,pres_1d,dam,uu1,vdf1)

!==========================================================
!  This program plots the increment of GSI analyzed field
!  and LAPS analysis.
!
!  HISTORY:
! Creation: MAR-2006 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  ! NetCDF:
  INCLUDE 'netcdf.inc'

  integer,intent(in) :: imax,jmax,kmax
  real*4,intent(in) :: znw(kmax),pres_1d(kmax),dam(imax,jmax),&
                       uu1(imax,jmax-1,kmax-1)
  real*4,intent(out) :: vdf1(imax,jmax-1,kmax-1)
  CHARACTER :: filename*60,varsname*10,variable*30
  CHARACTER :: filenamegsi*60,variablegsi*30
  INTEGER :: i,j,k,error,n(3),m(3)
  INTEGER :: ncid,vid,type,numd,dims(MAXVDIMS)
  INTEGER :: ncidgsi,vidgsi,typegsi,numdgsi,dimsgsi(MAXVDIMS)
  INTEGER :: nvatts,start(4),count(4)
  REAL,ALLOCATABLE,DIMENSION(:,:,:) :: v,vgsi,vdf
  REAL,ALLOCATABLE,DIMENSION(:,:) :: lat,lon

  REAL :: x,y,z,pi,tlat1,tlat2,emlon,bllat,bllon,trlat,trlon
  REAL :: fl,fr,fb,ft,ul,ur,ub,ut

  ! LAPS analysis:
  INTEGER :: n_var,ierr
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:,:) :: uvlaps,uvgsi
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:) :: parray

  ! Fort 10 and 12 are analysis files from wind.exe and GSIprep.exe:
  ALLOCATE(uvlaps(imax,jmax,kmax,2), &
           uvgsi (imax,jmax,kmax,2), &
           parray(imax,jmax,kmax))

  ! Convert to GSI coordinate:
  call laps2mass(uvlaps(1,1,1,1),imax,jmax,kmax, &
                 uvgsi(1,1,1,1),znw,pres_1d,dam)
  call laps2mass(uvlaps(1,1,1,2),imax,jmax,kmax, &
                 uvgsi(1,1,1,2),znw,pres_1d,dam)

  pi = 4.0*ATAN(1.0)

  fl = 0.0
  fr = 1.0
  fb = 0.0
  ft = 1.0
  ul = 0.01
  ur = 0.99
  ub = 0.01
  ut = 0.99

  ! Filename:
  filename = 'wrf_inout.before'         ! Before GSI
  filenamegsi = 'wrf_inout.after'         ! GSI analysis: use the differece

  ! Open the netcdf file:
  ncid = NCOPN(filename,NCWRITE,error)
  ncidgsi = NCOPN(filenamegsi,NCWRITE,error)

  ! Get variable ids:
  varsname = 'U'
  if (varsname(1:1) .EQ. 'U') then
    parray(1:imax,1:jmax,1:kmax) = uvgsi(1:imax,1:jmax,1:kmax,1)
  else
    parray(1:imax,1:jmax,1:kmax) = uvgsi(1:imax,1:jmax,1:kmax,2)
  endif

  vid = NCVID(ncid,varsname,error)
  vidgsi = NCVID(ncidgsi,varsname,error)

  ! Get variable info:
  CALL NCVINQ(ncid,vid,variable,type,numd,dims,nvatts,error)
  CALL NCVINQ(ncidgsi,vidgsi,variablegsi,typegsi, &
              numdgsi,dimsgsi,nvatts,error)

  ! Get the dimensions:
  DO i=1,numd
    CALL NCDINQ(ncid,dims(i),variable,n(i),error)
  ENDDO
  ! Allocate memory for u and v
  ALLOCATE(v(n(1),n(2),n(3)),vgsi(n(1),n(2),n(3)), &
           vdf(n(1),n(2),n(3)), STAT=error)

  ! Read the full grid:
  start = 1
  count(4) = 1
  count(1:3) = n

  ! Get the field:
  CALL NCVGT(ncid,vid,start,count,v(1,1,1),error)
  CALL NCVGT(ncidgsi,vidgsi,start,count,vgsi(1,1,1),error)
  vdf = vgsi-v
  vdf1 = uu1-v

  ! Get the true latitude:
  tlat1 = 38.9
  tlat2 = 90.0
  emlon = -121.05

  ! Lat/lon for CWB 9km domain:
  tlat1 = 10.0
  tlat2 = 40.0
  emlon = 120.9100

  vid = NCVID(ncid,'XLAT',error)
  ! Get variable info:
  CALL NCVINQ(ncid,vid,variable,type,numd,dims,nvatts,error)
  ! Get the dimensions:
  DO i=1,numd
    CALL NCDINQ(ncid,dims(i),variable,m(i),error)
  ENDDO
  ! Allocate memory for lat/lon:
  ALLOCATE(lat(m(1),m(2)),lon(m(1),m(2)),  STAT=error)
  start = 1
  count = 1
  count(1:numd) = m(1:numd)
  CALL NCVGT(ncid,vid,start,count,lat,error)

  vid = NCVID(ncid,'XLONG',error)
  ! Get variable info:
  CALL NCVINQ(ncid,vid,variable,type,numd,dims,nvatts,error)
  ! Get the dimensions:
  DO i=1,numd
    CALL NCDINQ(ncid,dims(i),variable,m(i),error)
  ENDDO
  start = 1
  count = 1
  count(1:numd) = m(1:numd)
  CALL NCVGT(ncid,vid,start,count,lon,error)

  bllat = lat(1,1)
  bllon = lon(1,1)
  trlat = lat(m(1),m(2))
  trlon = lon(m(1),m(2))

  PRINT*,'True Lat: ',tlat1,tlat2,emlon
  PRINT*,'Corners: ',bllat,bllon,trlat,trlon


  ! Close the netcdf file:
  CALL NCCLOS(ncid,error)

return

END  
