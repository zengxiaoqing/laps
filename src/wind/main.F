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

      program laps_wind_driver

      use mem_namelist, ONLY: read_namelist_laps
      use mem_namelist, ONLY: NX_L,NY_L,nk_laps,r_missing_data
     1                         ,iwrite_output

      use mem_grid,ONLY:lat,lon,topo
#ifdef USEMPI
     1                 ,nPEs,rank,IstartIend,recvcounts,displs
#endif

!      use mem_wind, wind%uanl->uanl,wind%vanl->vanl,wind%wanl->wanl
!     +          , wind%uanl_sfcitrp->uanl_sfcitrp
!     +          , wind%vanl_sfcitrp->vanl_sfcitrp
      use mem_wind

      implicit none

!     General declarations
      integer N_3D_FIELDS,len_dir,ntmin,ntmax,n_meso,n_sao,n_pirep
      integer istatus,istatus_alloc,istat_alloc,istat_ht,istat_wind
      integer istat_lw3
      parameter (N_3D_FIELDS = 3)

      character*3 var_3D(N_3D_FIELDS) 
      character*9 a9_time
      character*3 var_2d
      character*31 EXT

!     LAPS Grid Dimensions

      real, allocatable, dimension(:,:,:) :: heights_3d
      !real, allocatable, dimension(:,:,:) :: uanl
      !real, allocatable, dimension(:,:,:) :: vanl
      !real, allocatable, dimension(:,:) :: uanl_sfcitrp
      !real, allocatable, dimension(:,:) :: vanl_sfcitrp

!     Wind parameters

      character*150 static_dir,filename
!     logical l_use_raob, l_use_cdw, l_use_radial_vel

      integer i4time_sys

#ifdef USEMPI
!     MPI stuff
      include 'mpif.h'
      integer :: ierr,n,flag,count
      integer :: STATUS(MPI_STATUS_SIZE)
      integer :: PointsPerPE,extra,adjust

      integer :: max_obs,imax,jmax,kmax,ncnt_total
      integer :: idelt,jdelt,kdelt
      real    :: r_missing_data_local,r0_vert_grid

      integer,allocatable :: ibuf(:)
      real   ,allocatable :: rbuf(:)
      real   ,allocatable,dimension(:,:,:,:) :: wt_b_local,local_var
      real   ,allocatable :: fnorm_lut(:,:)
      real :: buffer !argument to mpi_gatherv ignored on non-root

      include 'barnesob.inc'
      include 'windparms.inc'

      type (barnesob),allocatable :: obs_barnes(:)

!     Set up the mpi tasks and spit process 0 from the rest
      call MPI_INIT( ierr )
      if( ierr /= 0) then
        print*,'MPI_INIT error',ierr
        stop
      endif
      call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
      if( ierr /= 0) then
        print*,'MPI_COMM_RANK error',ierr
        stop
      endif
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nPEs, ierr )
      if( ierr /= 0) then
        print*,'MPI_COMM_SIZE error',ierr,rank
        stop
      endif
      print *, 'Process ', Rank, ' of ', nPEs, ' is alive'
      allocate(IstartIend(2,0:nPEs-1), STAT=istat_alloc )
      if(istat_alloc .ne. 0)then
        print*,' ERROR: Could not allocate IstartIend',istat_alloc,rank
        stop
      endif
      allocate(recvcounts(  0:nPEs-1), STAT=istat_alloc )
      if(istat_alloc .ne. 0)then
        print*,' ERROR: Could not allocate recvcounts',istat_alloc,rank
        stop
      endif
      allocate(displs    (  0:nPEs-1), STAT=istat_alloc )
      if(istat_alloc .ne. 0)then
        print*,' ERROR: Could not allocate displs',istat_alloc,rank
        stop
      endif
      if(Rank > 0) then
       call mpi_bcast(IstartIend,2*nPEs,MPI_INTEGER,0,MPI_COMM_WORLD,
     1                                                          ierr)
       if( ierr /= 0) then
         print*,'MPI_BCAST IstartIend error',ierr,rank
         stop
       endif
       print*,'rank,Istart,Iend',rank,IstartIend(1,rank)
     1                               ,IstartIend(2,rank)
       call mpi_bcast(displs,nPEs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       if( ierr /= 0) then
         print*,'MPI_BCAST displs error',ierr,rank
         stop
       endif
       call mpi_bcast(recvcounts,nPEs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       if( ierr /= 0) then
         print*,'MPI_BCAST recvcounts error',ierr,rank
         stop
       endif
       print*,'Process',Rank,' is waiting'
       flag = 1
       do while(flag==1) 
         call MPI_BCAST(flag,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         if( ierr /= 0) then
           print*,'MPI_BCAST flag error',ierr,rank
           stop
         endif
         if(flag==-1) exit
         print*,'Process',Rank,' received go-ahead',flag,ierr
!        Receive and unpack integer buffers
         allocate(ibuf(10), STAT=istat_alloc )
         if(istat_alloc .ne. 0)then
           print*,' ERROR: Could not allocate ibuf(10)',istat_alloc,rank
           stop
         endif
         call mpi_bcast(ibuf,10,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         if( ierr /= 0) then
           print*,'MPI_BCAST ibuf error',ierr,rank
           stop
         endif
         if(n_var /= ibuf(1)) then
           print*,'Error in main.f: n_var /= ibuf(1)',n_var,ibuf(1)
           stop
         endif
         max_obs    = ibuf(2)
         imax       = ibuf(3)
         jmax       = ibuf(4)
         kmax       = ibuf(5)
         i4time_sys = ibuf(6)
         ncnt_total = ibuf(7)
         idelt      = ibuf(8)
         jdelt      = ibuf(9)
         kdelt      = ibuf(10)
         deallocate  (ibuf)
         count      = 4*ncnt_total
         allocate(ibuf(count), STAT=istat_alloc )
         if(istat_alloc .ne. 0)then
           print*,'ERROR Could not allocate ibuf',istat_alloc,rank,count
           stop
         endif
         call mpi_bcast(ibuf,count,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         if( ierr /= 0) then
           print*,'MPI_BCAST ibuf 1 error',ierr,count,rank
           stop
         endif
         allocate(obs_barnes(max_obs), STAT=istat_alloc )
         if(istat_alloc .ne. 0)then
           print*,' ERROR: Could not allocate abs_barnes',istat_alloc,
     1                         rank,max_obs
           stop
         endif
         do n=1,ncnt_total
           obs_barnes(n)%i      = ibuf(4*(n-1)+1)
           obs_barnes(n)%j      = ibuf(4*(n-1)+2)
           obs_barnes(n)%k      = ibuf(4*(n-1)+3)
           obs_barnes(n)%i4time = ibuf(4*(n-1)+4)
         enddo
         deallocate(ibuf)
!        Receive an unpack real buffer
         count = (2+max_var)*ncnt_total+2
         allocate(rbuf(count), STAT=istat_alloc )
         if(istat_alloc .ne. 0)then
           print*,'ERROR Could not allocate rbuf',istat_alloc,rank,count
           stop
         endif
         call mpi_bcast(rbuf,count,MPI_REAL,0,MPI_COMM_WORLD,ierr)
         if( ierr /= 0) then
           print*,'MPI_BCAST rbuf error',ierr,count,rank
           stop
         endif
         r_missing_data_local = rbuf(1)
         r0_vert_grid         = rbuf(2)
         count                = 2
         do n=1,ncnt_total
           obs_barnes(n)%value(1:max_var)=rbuf(count+1:count+max_var)
           obs_barnes(n)%vert_rad_rat    =rbuf(count+1      +max_var)
           obs_barnes(n)%weight          =rbuf(count+2      +max_var)
           count                         =count+2+max_var
         enddo
         deallocate(rbuf)

!        Receive arrays
         allocate(wt_b_local(imax,jmax,kmax,n_var), STAT=istat_alloc )
         if(istat_alloc .ne. 0)then
           print*,' ERROR: Could not allocate wt_b_local in main',
     1              istat_alloc
           print*,rank,imax,jmax,kmax,n_var
           stop
         endif
         allocate(fnorm_lut(-(imax-1):(imax-1),-(jmax-1):(jmax-1)), 
     1                                           STAT=istat_alloc )
         if(istat_alloc .ne. 0)then
           print*,' ERROR: Could not allocate fnorm_lut in main',
     1              istat_alloc,rank,imax,jmax
           stop
         endif
         call mpi_scatterv(buffer,recvcounts,displs,MPI_REAL,wt_b_local,
     1                  recvcounts(rank),MPI_REAL,0,MPI_COMM_WORLD,ierr)
         if( ierr /= 0) then
          print*,'MPI_SCATTERV error in main',ierr,rank,recvcounts(rank)
          stop
         endif
         count = (2*(imax-1)+1)*(2*(jmax-1)+1)
         call mpi_bcast(fnorm_lut,count,MPI_REAL,0,MPI_COMM_WORLD,ierr)
         if( ierr /= 0) then
           print*,'MPI_BCAST fnorm_lut error in main',ierr,count,rank
           stop
         endif
         allocate(local_var(IstartIend(1,rank):IstartIend(2,rank),jmax,
     1                                  kmax,n_var), STAT=istat_alloc )
         if(istat_alloc .ne. 0)then
           print*,'ERROR: Could not allocate local_var in main',
     1            istat_alloc,rank
           print*,IstartIend(1,rank),IstartIend(2,rank),jmax,kmax,n_var
           stop
         endif
         call analyze_weights_mpi(
     1        n_var,max_obs,obs_barnes                                    ! I
     1       ,imax,jmax,kmax                                              ! I
     1       ,wt_b_local                                                  ! I
     1       ,i4time_sys                                                  ! I
     1       ,r_missing_data_local,r0_vert_grid                           ! I
     1       ,ncnt_total,idelt,jdelt,kdelt                                ! I
     1       ,fnorm_lut                                                   ! I
     1       ,IstartIend(1,rank),IstartIend(2,rank)                       ! I
     1       ,local_var)                                                  ! O
         call mpi_gatherv(local_var,recvcounts(rank),MPI_REAL,buffer,
     1               recvcounts,displs,MPI_REAL,0,MPI_COMM_WORLD,ierr)
         if( ierr /= 0) then
           print*,'MPI_GATHERV error',ierr,rank,recvcounts(rank)
           stop
         endif

         deallocate(wt_b_local,fnorm_lut,obs_barnes,local_var)
       enddo !while
       call MPI_FINALIZE(ierr)
       stop
      endif !Rank > 0
      print*,'Process',Rank,' is proceeding'
#endif

!     Obtain global parameters

!     Read global parameters into module memory structure
      call get_directory('static',static_dir,len_dir)
      filename = static_dir(1:len_dir)//'/nest7grid.parms'
      call read_namelist_laps('lapsparms',filename)

!     Read wind parameters into module memory structure
      filename = static_dir(1:len_dir)//'/wind.nl'
      call read_namelist_laps('wind',filename)

csms$create_decomp(grid_dh, <nx_l>, <0>)
csms$serial(<r_missing_data, out>:default=ignore)  begin              
#ifdef USEMPI
!     Calculate the distribution of points in the I direction on each processor
      PointsPerPE = NX_L/nPEs
      if(PointsPerPE < 1) then
        print*,'Error in main.f'
        print*,'rank,PointsPerPe,NX_L,nPEs',rank,PointsPerPe,NX_L,nPEs
        call flush(6)
        stop
      endif
      extra = NX_L - PointsPerPE*nPEs
      if(extra > 0) then
        adjust=1
        extra = extra-1
      else
        adjust = 0
      endif
      IstartIend(1,0) = 1
      IstartIend(2,0) = PointsPerPE+adjust
      do n=1,nPEs-1
        if(extra > 0) then
          adjust=1
          extra = extra-1
        else
          adjust = 0
        endif
        IstartIend(1,n) = IstartIend(2,n-1)+1
        IstartIend(2,n) = IstartIend(1,n)+PointsPerPE-1+adjust
      enddo
      do n=0,nPEs-1
        print*,'Distribution:',n,IstartIend(1,n),IstartIend(2,n),NX_L
      enddo
      call mpi_bcast(IstartIend,2*nPEs,MPI_INTEGER,0,MPI_COMM_WORLD,
     1                                                         ierr)
      if(ierr/=0) then
        print*,'bcast IstartIend error',rank,2*nPEs,ierr
        stop
      endif
      do n=0,nPEs-1
        recvcounts(n)=(IstartIend(2,n)-IstartIend(1,n)+1)*NY_L*nk_laps
     1                                                           *n_var
      enddo
      displs(0) = 0
      do n=1,nPEs-1
        displs(n) = displs(n-1) + recvcounts(n-1)
      enddo
      call mpi_bcast(displs,nPEs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if( ierr /= 0) then
        print*,'MPI_BCAST displs',ierr,nPEs,rank
        stop
      endif
      call mpi_bcast(recvcounts,nPEs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if( ierr /= 0) then
        print*,'MPI_BCAST recvcounts',ierr,nPEs,rank
        stop
      endif
#endif

      NTMIN = -1
      NTMAX = +1

	call get_meso_sao_pirep(N_MESO,N_SAO,N_PIREP,istatus)
	if (istatus .ne. 1) then
	   write (6,*) 'Error obtaining N_MESO, N_SAO, and N_PIREP'
	   stop
	endif

!     Allocate static arrays (lat, lon, topo)
      allocate( lat(NX_L,NY_L), STAT=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate lat'
          stop
      endif

      allocate( lon(NX_L,NY_L), STAT=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate lon'
          stop
      endif

      allocate( topo(NX_L,NY_L), STAT=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate topo'
          stop
      endif

!     Read static arrays (lat, lon, topo)
      call read_static_grid(NX_L,NY_L,'LAT',lat,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Error getting LAPS LAT'
          stop
      endif

      call read_static_grid(NX_L,NY_L,'LON',lon,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Error getting LAPS LON'
          stop
      endif

      call read_static_grid(NX_L,NY_L,'AVG',topo,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Error getting LAPS topo'
          stop
      endif

!     Get system time
      call get_systime(i4time_sys,a9_time,istatus)
      if(istatus .ne. 1)go to 999
      write(6,*)' systime = ',a9_time

!     Allocate global arrays
      allocate( heights_3d(NX_L,NY_L,nk_laps), STAT=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate heights_3d'
          stop
      endif

      allocate( wind%uanl(NX_L,NY_L,nk_laps), STAT=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate uanl'
          stop
      endif

      allocate( wind%vanl(NX_L,NY_L,nk_laps), STAT=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate vanl'
          stop
      endif

      allocate( wind%uanl_sfcitrp(NX_L,NY_L), STAT=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate uanl_sfcitrp'
          stop
      endif

      allocate( wind%vanl_sfcitrp(NX_L,NY_L), STAT=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate vanl_sfcitrp'
          stop
      endif

csms$serial end

      write(6,*)
      write(6,*)' Getting 3_D Height Analysis from MODEL:'

      istat_ht = 0
      var_2d = 'HT'
      call get_modelfg_3d(i4time_sys,var_2d
     1                             ,NX_L,NY_L,nk_laps
     1                             ,heights_3d,istat_ht)

      if(istat_ht .ne. 1)then
          write(6,*)
     1    ' Aborting from LAPS Wind Anal - Error in getting MDL heights'
          goto 999
      endif

      call lapswind_anal  (i4time_sys,
     1            NTMIN,NTMAX,                                     ! I 
     1            N_MESO,N_SAO,                                    ! I
     1            heights_3d,                                      ! I
     1            istat_wind)                                      ! O

      if(istat_wind .eq. 1 .and. iwrite_output .ge. 0)then
          EXT = 'lw3'
          call write_wind_output(i4time_sys,EXT,var_3d                
     1                            ,wind%uanl,wind%vanl                   ! I
     1                            ,wind%wanl                             ! I
     1                            ,wind%uanl_sfcitrp,wind%vanl_sfcitrp   ! I
     1                            ,NX_L,NY_L,nk_laps,N_3D_FIELDS         ! I
     1                            ,r_missing_data                        ! I
     1                            ,istat_lw3)
      endif

!     if (allocated(wind%wanl) ) deallocate(wind%wanl)

999   continue
#ifdef USEMPI
      flag=-1
      call mpi_bcast(flag,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if( ierr /= 0) then
        print*,'MPI_BCAST flag final error',ierr,count,rank
        stop
      endif
      call MPI_FINALIZE(ierr)
      stop
#endif

csms$exit
      end

