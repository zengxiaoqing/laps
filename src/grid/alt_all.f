CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine alt_10by10_all(im,jm,ID2,testrad,
     &          glatr,glonr,path,out2d,out3d,ncat,
     &		slope_lt,slope_ln,stdev,categorical,WRDLEN,
     &		IGNOREVAL)

        parameter(r2d=57.29577951,dtr=0.017453293)

        INTEGER GDS(200),WRDLEN
        character*2 uname
        character*8 filename
        character*(*) path
        character*256 fullname

        REAL out2d(im,jm),out3d(im,jm,NCAT),stdev(im,jm)
	REAL slope_lt(im,jm),slope_ln(im,jm)
	logical pure_nn,categorical,topogen
        real glatr(im,jm),glonr(im,jm)
        integer*2, allocatable:: topoin(:,:)
        integer*2  topotmp(ID2,ID2), IGNOREVAL
	integer*1 topotmp_1(ID2,ID2)
        real dum2d(im,jm)
	real, allocatable:: topolat(:,:),topolon(:,:)
	real, allocatable:: xpts(:,:),ypts(:,:),dum(:,:)



!!	ID2=1200 is 30" topo, landuse, soil
!!	ID2=10 is 1 degree soil temp


	if (ID2 .eq. 1200) then
	IPERDEG=ID2/10
	ISPAN=10
	pure_nn=.false.
	elseif (ID2 .eq. 10) then
	IPERDEG=1
	ISPAN=10
	pure_nn=.true.
	else
	write(6,*) 'not appropriate...quitting'
	STOP
	endif
!
	out2d=-999.
!
        apmx=maxval(glatr)
        apmn=minval(glatr)

	if (glatr(im/2,jm/2) .gt. 0) then
        aloneast=glonr(im,jm)
	else
	aloneast=glonr(im,1)
	endif

        alonwest1=glonr(1,1)
        alonwest2=glonr(1,jm)

	if (apmn .gt. 0) then
!NH
        if ( (alonwest1 .lt. 0. .and. alonwest2 .lt. 0.) .or. 
     &	     (alonwest1 .gt. 0. .and. alonwest2 .gt. 0.) .or. 
     &	     (alonwest1 .gt. 0. .and. alonwest2 .lt. 0.) ) then
          alonwest=amin1(alonwest1,alonwest2)
	elseif (alonwest1 .lt. 0 .and. alonwest2 .gt. 0) then
! NH western boundary straddles DL
	  alonwest=amax1(alonwest1,alonwest2)
	endif

	else
!SH
        if ( (alonwest1 .lt. 0 .and. alonwest2 .lt. 0) .or. 
     &	     (alonwest1 .gt. 0 .and. alonwest2 .gt. 0) .or. 
     &	     (alonwest1 .lt. 0. .and. alonwest2 .gt. 0.) ) then
          alonwest=amin1(alonwest1,alonwest2)
	elseif (alonwest1 .gt. 0 .and. alonwest2 .lt. 0) then
! NH western boundary straddles DL
	  alonwest=amax1(alonwest1,alonwest2)
	endif

	endif



        ISTRAD=0

        if (alonwest .gt. 0 .and. aloneast .lt. 0) then
        write(6,*) 'straddling the DL!!!'
        ISTRAD=1
        endif

        if (alonwest .lt. -180.) then
        write(6,*) 'straddling the DL!!!'
        ISTRAD=1
        endif

        if (alonwest .lt. 0 .and. aloneast .gt. 0) then
        write(6,*) 'straddling the PM!!!'
        ISTRAD=2
        endif

        rnlat=apmx+0.1
        slat=apmn-0.1
        west=almn-0.1
        west=alonwest-0.1
        east=aloneast+0.1

        itemp=int(west/ISPAN)
        if (west .lt. 0) then
        west=float(itemp-1)*ISPAN
        else
        west=float(itemp)*ISPAN
        endif

        itemp=int(east/ISPAN)
        if (east .gt. 0) then
        east=float(itemp+1)*ISPAN
        else
        east=float(itemp)*ISPAN
        endif

        itemp=int(rnlat/ISPAN)
        rnlat=float(itemp+1)*ISPAN
        itemp=int(slat/ISPAN)

	if (slat .gt. 0) then
	slat=float(itemp)*ISPAN
	else
	slat=float(itemp)*ISPAN-ISPAN
	endif


        write(6,*) 'limits on discretized grid'
        write(6,*) 'west= ', west
        write(6,*) 'east= ', east
        write(6,*) 'north= ', rnlat
        write(6,*) 'south= ', slat


        if (ISTRAD .eq. 0 .or. ISTRAD .eq. 2) then
        ILIM=IPERDEG*(east-west)
        JLIM=IPERDEG*(rnlat-slat)
        endif
        if (ISTRAD .eq. 1) then
        ILIM=IPERDEG*((east+360.)-west)
        JLIM=IPERDEG*(rnlat-slat)
        endif

!
!	Create 10 degree longitude strip
!
	ILIM=ID2

        allocate(topoin(ilim,jlim))
        allocate(topolat(ilim,jlim))
        allocate(topolon(ilim,jlim))
        allocate(xpts(ilim,jlim))
        allocate(ypts(ilim,jlim))
        allocate(dum(ilim,jlim))

        do j=1,jlim
        do i=1,ilim
        xpts(i,j)=i
        ypts(i,j)=j
        enddo
        enddo

        if (ISTRAD .eq. 1) then
        easttmp=east+360
        else
        easttmp=east
        endif

         write(6,*) 'loop I from: ', int(west), int(easttmp)-ISPAN

!!!!!!!!
!!!!!!!!
        do IIII=int(west),int(easttmp)-ISPAN,ISPAN
!!!!!!!!
!!!!!!!!

        if (east .ne. easttmp) then
          if (IIII .ge. 180) then
           ITMP=IIII-360
          else
           ITMP=IIII
          endif
        else
          ITMP=IIII
        endif

!!!!!!!!
!!!!!!!!
        do J=int(slat),int(rnlat)-ISPAN,ISPAN
!!!!!!!!
!!!!!!!!

        call gen_filename(J,ITMP,filename)
        write(6,*) 'want to process ', filename


	call s_len(path,len)

        fullname=path(1:len)//filename(1:7)

	IF (WRDLEN .eq. 2) then

!ENDIANISSUE
        open(unit=1,file=trim(fullname),status='old'
     .          ,form='unformatted',access='direct',recl=2*ID2)
        do nrrd=1,ID2
        read(1,rec=nrrd) (topotmp(ii,nrrd),ii=1,ID2)

!!! ADD AN #IFDEF structure, SWAP if needed?

	if (nrrd .eq. ID2/2) then
	write(6,*) 'read vals: ', (topotmp(II,nrrd),II=1,ID2,ID2/10)
	endif

        enddo
        close(1)

	ELSEIF (WRDLEN .eq. 1) then
!ENDIANISSUE
        open(unit=1,file=trim(fullname),status='old'
     .          ,form='unformatted',access='direct',recl=1*ID2)
        do nrrd=1,ID2
        read(1,rec=nrrd) (topotmp_1(ii,nrrd),ii=1,ID2)
        enddo
        close(1)
	
	ENDIF

CC      put this tile into longitudinal strip

        do JJJ=(J-int(slat))*(IPERDEG)+1,(J-int(slat))*(IPERDEG)+ID2
        do III=1,ID2
        JLOC=JJJ-(J-int(slat))*(IPERDEG)
		if (WRDLEN .eq. 2) then
       			 topoin(III,JJJ)=topotmp(III,ID2-JLOC+1)
		elseif (WRDLEN .eq. 1) then
        		 topoin(III,JJJ)=topotmp_1(III,ID2-JLOC+1)
		endif
        enddo
        enddo

        enddo ! enddo for J

        write(6,*) 'ILIM, JLIM: ', ILIM, JLIM


CC      come up with an estimate of how many input 30 s points should
CC      be included.  Consider roughly enough to average over the size of
CC      the gridbox.

       limit=1

CC
CC      30" ~ 0.927 km
CC      10' ~ 18.5 km
CC

CC      limit represents an intentionally over-large search radius to be
CC      used below to specify which input points define the output topography

!        limit=int((testrad/2.)/0.927)+2

        resfac=(1./IPERDEG)*(111.2)
        limit=int((testrad/2.)/resfac)+2

!        write(6,*) 'if possible, will search +/- ', limit

        GDS=0

        GDS(1)=0
        GDS(2)=ILIM
        GDS(3)=JLIM
        GDS(4)=INT(SLAT*1000)
        GDS(5)=INT(IIII*1000)
        GDS(6)=128
        GDS(7)=INT((RNLAT-1./(IPERDEG))*1000)
        GDS(8)=INT((IIII+10-1./(IPERDEG))*1000)
        GDS(9)=INT(1./(IPERDEG)*1000)
        GDS(10)=INT(1./(IPERDEG)*1000)

	write(6,*) 'GDS= ', (GDS(I),I=1,10)

	CALL GDSWIZ(GDS,1,ILIM*JLIM,-9999.,XPTS,YPTS,
     &	topolon,topolat,NRET,0,dum,dum)

	do J=1,JLIM
	 do I=1,ILIM
	  if (TOPOLON(I,J) .gt. 180) then
	   TOPOLON(I,J)=TOPOLON(I,J)-360.
	  endif
         enddo
	enddo


!!!!
	if (categorical) then
	write(6,*) 'calling alt_categories ', PURE_NN
		call alt_categories(IM,JM,GLATR,GLONR,
     &		ILIM,JLIM,TOPOIN,TOPOLAT,TOPOLON,PURE_NN,OUT3D,
     &		NCAT,OUT2D,GDS,limit,testrad)
	else
	write(6,*) 'calling alt_interp ', PURE_NN

	if (path(len:len) .eq. 'U') then
		topogen=.true.
	else
		topogen=.false.
	endif
		call alt_interp(IM,JM,GLATR,GLONR,
     &		ILIM,JLIM,TOPOIN,TOPOLAT,TOPOLON,PURE_NN,OUT2D,SLOPE_LT,
     &		SLOPE_LN,STDEV,GDS,limit,testrad,IGNOREVAL,topogen)
	endif

!	write(6,*) 'after this tile: '

!	do J=JM,1,-JM/20
!	write(6,633) (out2d(I,J),I=1,IM,IM/15)
!	enddo

  633	format(15(f5.0,1x))


!!!!!!!!
!!!!!!!!
	enddo ! enddo for IIII (looping over strips)
!!!!!!!!
!!!!!!!!

CCCCCCCCCCCCC

	do J=1,JM
	do I=1,IM
	if (out2d(I,J) .eq. -999.) then
	write(6,*) 'undefined at ', i,j
	write(6,*) 'glat,glon: ', glatr(i,j),glonr(i,j)

	rmax=-999.
	do JJ=J-1,J+1
	do II=I-1,I+1
		if (out2d(II,JJ) .gt. rmax) then
			rmax=out2d(II,JJ)
		endif
	enddo
	enddo

	write(6,*) 'redefined out2d to be ', rmax
	out2d(I,J)=rmax

	endif
	enddo
	enddo


        deallocate(topoin,topolat,topolon,xpts,ypts,dum)

        END subroutine alt_10by10_all

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	subroutine alt_categories(IM,JM,GLATR,GLONR,
     &          ILIM,JLIM,DATAIN,DATALAT,DATALON,
     &		PURE_NN,CAT3D,
     &          NCAT,CAT2D,GDS,LIMIT,res)

	REAL GLATR(IM,JM),GLONR(IM,JM)
	REAL DATALAT(ILIM,JLIM),DATALON(ILIM,JLIM)
	INTEGER*2 DATAIN(ILIM,JLIM)
	REAL CAT2D(IM,JM),CAT3D(IM,JM,NCAT)

	INTEGER GDS(200), COUNTER(NCAT)
	LOGICAL PURE_NN

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	write(6,*) 'search limit is : ', limit
	write(6,*) 'NCAT= ', ncat
        do J=1,JM
        do I=1,IM
                                                                                
C DETERMINE LAT/LON OF TARGET (output) GRID POINT
C
        GLATD = GLATR(i,j)
        GLOND = GLONR(I,J)
                                                                                
        call ced_ij(glatd,glond,x,y,gds)
                                                                                
!! kludgy fix for dead space "between" strips (like -60.003 longitude)
         
        if (X .le. 1200.5) then
        II=INT(X+0.5)
        else
!       write(6,*) 'kludgy fix applied...'
        II=INT(X)
        endif
         
        JI=INT(Y+0.5)
         
!       II,JI represents nearest neighbor from input array
         
        IF (II .ge. 1 .and. II .le. ILIM .and.
     +      JI .ge. 1 .and. JI .le. JLIM ) THEN
         
	IF (.NOT. PURE_NN) THEN

        ICOUNT=0
        sum=0.
        counter=0
         
         
!       no strictly n.n. now.  Test as many points within the limits
!       as possible.

	
        do JJJ=JI-limit,JI+limit
           do III=II-limit,II+limit
                                                                                
        if (JJJ .lt. 1 .or. JJJ .gt. JLIM .or.
     &      III .lt. 1 .or. III .gt. ILIM) goto 47
                                                                                
        call greatcir(glatd,glond,
     &          datalat(III,JJJ),datalon(III,JJJ),dist)
                                                                                

        if ( dist .le. res/2.) then
        counter(datain(III,JJJ))=counter(datain(III,JJJ))+1
        ICOUNT=ICOUNT+1
        endif
                                                                                
  47    continue
                                                                                
           enddo
        enddo
C
                                                                                
	IF (ICOUNT .gt. 0) THEN

        imaxcount=-9
        do N=1,NCAT
                                                                                
        cat3d(i,j,N)=float(counter(N))/ICOUNT
                                                                                
        if (counter(N) .gt. imaxcount) then
                imaxcount=counter(N)
                IBEST=N
        endif
                                                                                
        enddo
                                                                                
        cat2d(i,j)=IBEST

	if (MOD(I,10) .eq. 0 .and. mod(J,10) .eq. 0) then
!	write(6,*) 'I,J,cat2d: ', i,j,cat2d(i,j)
!	write(6,*) 'ICOUNT,IMAXCOUNT ', ICOUNT,IMAXCOUNT
	endif
	
	ENDIF

	ELSE ! pure N.N. case

!!! dont believe this section will/should ever be used

	cat2d(i,j)=datain(II,JI)
	
	do N=1,NCAT
	if (cat2d(I,J) .eq. N) then
	  cat3d(i,j,N)=1.0
	else
	  cat3d(i,j,N)=0.0
	endif
	enddo
	
	ENDIF
                                                                                
        ELSE
! outside of this strip...ignore
                                                                                
        ENDIF
                                                                                
                                                                                
        enddo ! enddos for im,jm loop
        enddo

	end subroutine alt_categories
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

         subroutine alt_hemi_all(im,jm,glatr,glonr,path,
     &     type,NCAT,res,cat3d,cat2d,categorical,outland,useland)

	integer counter(NCAT)

        parameter(r2d=57.29577951,dtr=0.017453293)

        character*12 fnameout
        INTEGER GDS(200)
        character*2 uname
	character*1 type
        character*8 filename
	character*7 geoloc(2)
        character*(*) path
        character*256 fullname
        character*256 title,unit_name

        real cat3d(im,jm,NCAT)
        real cat2d(im,jm),x(im,jm),y(im,jm)
        real glatr(im,jm),glonr(im,jm)
        real dum2d(im,jm)
	real, allocatable:: datalat(:,:),datalon(:,:)
	real, allocatable:: xpts(:,:),ypts(:,:),dum(:,:)
	integer, allocatable:: idata(:,:,:)

!mp	respecify with proper "kind" type statment????
	integer*2, allocatable:: data2d(:,:)

	integer*2 ignoreval

!	
	real outland(im,jm) ! (land=1, water=0)

	real leftlon, rightlon

	logical have_left_east, have_left_west
	logical have_right_east, have_right_west
	logical read_west, read_east, read_both
	logical categorical, pure_nn, useland
	logical redo(im,jm), topogen

!
	cat3d=-99999.
	cat2d=-99999.
!



	have_left_east=.false.
	have_left_west=.false.
	have_right_east=.false.
	have_right_west=.false.

	do J=1,JM

	if (glonr(1,J) .lt. 0) then
	  have_left_west=.true.
	elseif (glonr(1,J) .gt. 0) then
	  have_left_east=.true.
	endif

	if (glonr(IM,J) .lt. 0) then
	  have_right_west=.true.
	elseif (glonr(IM,J) .gt. 0) then
	  have_right_east=.true.
	endif

	enddo


!	from this information can gather the straddling situation, and
!	devise a way to find west and east longitude


	read_west=.false.
	read_east=.false.
	read_both=.false.

	geoloc(1)='1234567'
	geoloc(2)='1234567'

	if (have_left_west .and. have_right_east) then
	   write(6,*) 'straddles PM'
	   read_both=.true.
	geoloc(1)='90S180W'
	geoloc(2)='90S000E'
	elseif (have_left_east .and. have_right_west) then
	   write(6,*) 'straddles DL'
	   read_both=.true.
	geoloc(1)='90S180W'
	geoloc(2)='90S000E'
	else
	   read_both=.false.
	   write(6,*) 'purely in one hemisphere'
	   if (have_left_west .and. have_right_west) then
	      read_west=.true.
	geoloc(1)='90S180W'
	   elseif (have_left_east .and. have_right_east) then
	      read_east=.true.
	geoloc(1)='90S000E'
	   endif
	endif

        call s_len(path,lenp)
        TITLE=path(1:lenp)//'HEADER'
              lentd=INDEX(TITLE,' ')-1
      CALL JCLGET(29,TITLE(1:lentd),'FORMATTED',1,istatus)
      if(istatus .ne. 1)then
         write(6,*)'Warning: proc_geodat_tiles opening HEADER: check'
     1            ,'geog paths and HEADER file'
         return
      endif
      READ(29,*)IBLKSIZO,NO,ISBEGO,IWBEGO,RWOFF,RSOFF
      print *,'title=',title(1:lentd)
      print *,'RWOFF,RSOFF = ',RWOFF,RSOFF
      print *,'isbego,iwbego=',isbego,iwbego
      print *,'iblksizo,no=',iblksizo,no
      CLOSE(29)

	write(6,*) 'IBLKSIZO= ', IBLKSIZO

	ILIM=NO
	JLIM=NO

        write(6,*) 'ILIM, JLIM: ', ILIM,JLIM
	if (allocated(datalat)) deallocate(datalat)
	if (allocated(datalon)) deallocate(datalon)
	if (allocated(xpts)) deallocate(xpts)
	if (allocated(ypts)) deallocate(ypts)
	if (allocated(dum)) deallocate(dum)

        allocate(datalat(ilim,jlim))
        allocate(datalon(ilim,jlim))
        allocate(xpts(ilim,jlim))
        allocate(ypts(ilim,jlim))
        allocate(dum(ilim,jlim))

        do j=1,jlim
        do i=1,ilim
        xpts(i,j)=i
        ypts(i,j)=j
        enddo
        enddo

	IF (read_both) THEN
	  NTILES=2
        ELSE
	  NTILES=1
	ENDIF	

!!!!
!!!!
	DO NN=1,NTILES
!!!!
!!!!

        nn1=ILIM
        nn2=ILIM
	nn4=NCAT

	if (allocated(idata)) deallocate(idata)
	if (allocated(data2d)) deallocate(data2d)

	write(6,*) 'allocated with dims: ', nn4,nn1,nn2
        allocate (idata(nn4,nn1,nn2))

	write(6,*) 'geoloc(NN): ', geoloc(NN)

	if (geoloc(NN)(4:7) .eq. '180W') then
		 WSTART=-180
	elseif (geoloc(NN)(4:7) .eq. '000E') then
		 WSTART=0
	else
		write(6,*) 'not right'
		STOP
	endif

	write(6,*) 'WSTART= ', WSTART

	unit_name=path(1:LENP)//geoloc(NN)
	write(6,*) 'unit_name= ', unit_name

        call s_len(unit_name,len)

	i1=1
	i2=4

!	write(6,*) 'nn1*nn2*nn4: ', nn1*nn2*nn4

!ENDIANISSUE???
        call read_binary_field(idata,i1,i2,nn1*nn2*nn4,
     &          unit_name,len)

!        N=1
!        write(6,*) 'NDATA: ', N
!        do J=NN2,1,-NN2/35
!        write(6,633) (idata(N,I,J),I=1,NN1,NN1/20)
!        enddo


  633	format(40I3)

	
CC      put this tile into longitudinal strip


!        write(6,*) 'ILIM, JLIM: ', ILIM, JLIM


CC      come up with an estimate of how many input 30 s points should
CC      be included.  Consider roughly enough to average over the size of
CC      the gridbox.

       limit=1

CC	limit represents an intentionally over-large search radius to be
CC	used below to specify which input points define the output point

!!	are these defs proper??

	DELTLAT=float(IBLKSIZO)/NO
	DELTLON=float(IBLKSIZO)/NO

	write(6,*) 'NO,IBLKSIZO, DELTLAT: ', NO,IBLKSIZO, DELTLAT

        limit=int((res/2.)/(DELTLAT*111.2))+2

        GDS=0
	
	SLAT=ISBEGO+RSOFF
	WLAT=WSTART+RWOFF

	RNLAT=90-RSOFF
	ELAT=WLAT+180

        GDS(1)=0
        GDS(2)=ILIM
        GDS(3)=JLIM
        GDS(4)=INT(SLAT*1000)
        GDS(5)=INT(WLAT*1000)
        GDS(6)=128
        GDS(7)=INT((RNLAT)*1000)
        GDS(8)=INT(ELAT*1000)
	GDS(9)=INT(DELTLAT*1000)
	GDS(10)=INT(DELTLON*1000)

	write(6,*) 'deltalat, deltalon: ', deltlat, deltlon

!	write(6,*) 'GDS= ', GDS
!	write(6,*) 'ILIM*JLIM= ', ILIM*JLIM

	CALL GDSWIZ(GDS,1,ILIM*JLIM,-9999.,XPTS,YPTS,
     &	datalon,datalat,NRET,0,dum,dum)

	do J=1,JLIM
	 do I=1,ILIM
	  if (DATALON(I,J) .gt. 180) then
	   DATALON(I,J)=DATALON(I,J)-360.
	  endif
         enddo
	enddo


	if (DELTLAT .lt. 1) then
	PURE_NN=.false.
	else
	PURE_NN=.true.
	endif


	allocate(data2d(ILIM,JLIM))

        if (categorical) then

        write(6,*) 'calling alt_categories!! ', PURE_NN
	write(6,*) 'dont think this is proper!!!'

                call alt_categories(IM,JM,GLATR,GLONR,
     &          ILIM,JLIM,IDATA,DATALAT,DATALON,PURE_NN,CAT3D,
     &          NCAT,CAT2D,GDS,limit,res)

        else

	do N=1,NCAT

	do J=1,JLIM
	do I=1,ILIM
	data2d(I,J)=idata(N,I,J)
	enddo
	enddo

!	always want to ignore zeros for albedo, green frac, slope

	ignoreval=0
	
	topogen=.false.

                call alt_interp(IM,JM,GLATR,GLONR,
     &          ILIM,JLIM,DATA2D,DATALAT,
     &		DATALON,PURE_NN,CAT2D,DUM2D,
     &          DUM2D,DUM2D,GDS,limit,res,ignoreval,topogen)

	
	IF (USELAND) THEN

! check if any output land points have a zero value for field.
! If so, assign an average of nearest non-zero value points within
! search radius

	ISRCH=3

	DO J=1,JM
	DO I=1,IM

	redo(I,J)=.false.

	if (outland(I,J) .eq. 1 .and. cat2d(I,J) .eq. 0) then

!! trouble point

	ipts=0
	sum=0.

	do JJ=J-ISRCH,J+ISRCH
	do II=I-ISRCH,I+ISRCH

	if (II .ge. 1 .and. II .le. IM .and. 
     &		JJ .ge. 1 .and. JJ .le. JM) then

	if (cat2d(II,JJ) .gt. 0) then
	sum=sum+cat2d(II,JJ)
	ipts=ipts+1
	endif

	endif

	enddo
	enddo ! end II, JJ loops

	if (ipts .gt. 0) then
!	  write(6,*) 'based on  ', ipts, 'valid points'
!	  write(6,*) 'derived new value of ', sum/ipts
	  cat2d(I,J)=sum/ipts
	else
	  write(6,*) 'found no valid points at ' , i,j
	  redo(i,j)=.true.
	endif

	endif

	ENDDO
	ENDDO


!new

	

	do J=1,JM
	 do I=1,IM
	ISRCH=3
	  do while (redo(I,J) .and. ISRCH .lt. IM/2)

! spiral outward as needed to find a valid value
	ISRCH=ISRCH+1

	do JJ=J-ISRCH,J+ISRCH
	do II=I-ISRCH,I+ISRCH

	IF (II .ge. 1 .and. II .le. IM .and. 
     &		JJ .ge. 1 .and. JJ .le. JM) then

	if (cat2d(II,JJ) .gt. 0) then
	  cat2d(I,J)=cat2d(II,JJ)
	  write(6,*) 'defined ', i,j, 'to be: ', cat2d(i,j)
	  write(6,*) 'found value at search : ', ISRCH
	  redo(i,j)=.false.
	  goto 99
	endif

	ENDIF

	enddo
	enddo


  99	continue

	enddo

	if (redo(I,J)) then
	write(6,*) 'giving up on point...'
	cat2d(I,J)=0.25 ! a reasonable default for albedo/gfrac ???
	write(6,*) 'setting to default val ', cat2d(I,J)
	endif

	enddo
	enddo
	ENDIF ! endif useland block

!endnew

	do J=1,JM
	do I=1,IM
	CAT3D(I,J,N)=CAT2D(I,J)
	enddo
	enddo

	enddo ! ncat loop

  467	format(25(f3.0,1x))

!	each month will be interpolated, then dumped into the 3D
!	array.  Will it work for slope/maxsnowalb?
	
        endif ! categorical/interp branch

	ENDDO ! for tiles

CCCCCCCCCCCCC

        deallocate(datalat,datalon,xpts,ypts,dum)

        END subroutine alt_hemi_all

C +------------------------------------------------------------------+
	subroutine alt_interp(IM,JM,GLATR,GLONR,
     &          ILIM,JLIM,DATAIN,DATALAT,DATALON,
     &		PURE_NN,OUT2D,SLOPE_LT,
     &          SLOPE_LN,STDEV,GDS,LIMIT,TESTRAD,
     &		IGNOREVAL,topogen)


	REAL GLATR(IM,JM),GLONR(IM,JM)
	REAL DATALAT(ILIM,JLIM),DATALON(ILIM,JLIM)
	INTEGER*2 DATAIN(ILIM,JLIM),IGNOREVAL
	REAL OUT2D(IM,JM),SLOPE_LT(IM,JM),SLOPE_LN(IM,JM)
	REAL STDEV(IM,JM)

	INTEGER GDS(200)

	LOGICAL PURE_NN, topogen

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!	write(6,*) 'testrad = ', testrad

!	write(6,*) 'datain for alt_interp, lims ', ILIM, JLIM
!	do J=JLIM,1,-JLIM/20
!	write(6,647) (DATAIN(I,J),I=1,ILIM,ILIM/15)
!	enddo
  647	format(25I3)
	

        do J=1,JM
        do I=1,IM
                                                                                
C DETERMINE LAT/LON OF TARGET (output) GRID POINT

        GLATD = GLATR(i,j)
        GLOND = GLONR(I,J)
                                                                                
        call ced_ij(glatd,glond,x,y,gds)
                                                                                
!! kludgy fix for dead space "between" strips (like -60.003 longitude)

        if (X .le. ID2+0.5) then
        II=INT(X+0.5)
        else
!       write(6,*) 'kludgy fix applied...'
        II=INT(X)
        endif
                                                                                
        JI=INT(Y+0.5)
                                                                                
!       II,JI represents nearest neighbor from input data array
                                                                                
        IF (II .ge. 1 .and. II .le. ILIM .and.
     +      JI .ge. 1 .and. JI .le. JLIM ) THEN
                                                                                
                IF (.NOT. PURE_NN) THEN
                                                                                
        ICOUNT=0
	IGNORCNT=0
        sum=0.
        icnt_lt=0
        icnt_ln=0
        sum_lt=0.
        sum_ln=0.
        do JJJ=JI-limit,JI+limit
           do III=II-limit,II+limit
                                                                                
        if (JJJ .lt. 1 .or. JJJ .gt. JLIM .or.
     &      III .lt. 1 .or. III .gt. ILIM) goto 47
                                                                                
        call greatcir(glatd,glond,
     &          datalat(III,JJJ),datalon(III,JJJ),dist)

        IF ( dist .le. testrad/2.) THEN
	  IF (datain(III,JJJ) .ne. IGNOREVAL) then
              sum=sum+real(datain(III,JJJ))
              ICOUNT=ICOUNT+1
	  ELSE
              IGNORCNT=IGNORCNT+1
	  ENDIF

!	
                                                                                
        if (III+1 .le. ILIM .and. JJJ+1 .le. JLIM ) then
        icnt_ln=icnt_ln+1
        icnt_lt=icnt_lt+1
        sum_ln = sum_ln + 0.5*(datain(III+1,JJJ)-datain(III,JJJ)+
     &                         datain(III+1,JJJ+1)-datain(III,JJJ+1))
        sum_lt = sum_lt + 0.5*(datain(III,JJJ+1)-datain(III,JJJ)+
     &                         datain(III+1,JJJ+1)-datain(III+1,JJJ))
        endif
                                                                                
        ENDIF
                                                                                
  47    continue
           enddo
        enddo
C
                                                                                
        slope_ln(i,j)=(sum_ln/float(icnt_ln))/(1000.*testrad)
        slope_lt(i,j)=(sum_lt/float(icnt_lt))/(1000.*testrad)
        out2d(i,j)=sum/float(icount)

	if (icount .eq. 0) then ! force nearest neighbor
	out2d(i,j)=datain(ii,ji)

!!	in this case, can feel comfortable that is not topo
!!	so dont worry about stdev and slope computations

	endif


!!	if predominantly zero values, dominant land/soil type will be water
!!	so set topo to zero here.  Looking for better topo/landmask 
!!	agreement

	if (topogen .and. IGNORCNT .gt. ICOUNT  .and. 
     &		out2d(i,J) .gt. 0) then
	write(6,*) 'set topo to zero from value ', out2d(i,j)
	out2d(I,J)=0.
	endif

                                                                                
!!      now essentially redo to get standard deviation info

!! 1002 - NOT USING BELOW FOR ANYTHING...COMMENT OUT?
                                                                                
        ICOUNT=0
        sum=0.
                                                                                
        do JJJ=JI-limit,JI+limit
           do III=II-limit,II+limit
                                                                                
        if (JJJ .ge. 1 .and. JJJ .le. JLIM .and.
     &      III .ge. 1 .and. III .le. ILIM) then
                                                                                
        call greatcir(glatd,glond,
     &          datalat(III,JJJ),datalon(III,JJJ),dist)
                                                                                
        IF ( dist .le. testrad/2.) THEN
              sum=sum+(out2d(I,J)-real(datain(III,JJJ)))**2.
              ICOUNT=ICOUNT+1
        ENDIF
                                                                                
        endif
                                                                                
                                                                                
        enddo
        enddo
                                                                                
        stdev(I,J)=(sum/float(icount))**(0.5)
                                                                                
!! end standard deviation
                ELSE ! pure N.N. case
                                                                                
                stdev(i,j)=0.
                                                                                
                out2d(i,j)=datain(II,JI)
                slope_ln(i,j)=0.
                slope_lt(i,j)=0.
                                                                                
                ENDIF
                                                                                
        ELSE
! outside of this strip...ignore
                                                                                
        ENDIF
                                                                                
        enddo ! enddos for IM,JM loop
        enddo

	end subroutine alt_interp


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        subroutine ced_ij(RLAT,RLON,XPTS,YPTS,KGDS)

        integer kgds(200)

        IM=KGDS(2)
        JM=KGDS(3)
        RLAT1=KGDS(4)*1.E-3
        RLON1=KGDS(5)*1.E-3
        RLAT2=KGDS(7)*1.E-3
        RLON2=KGDS(8)*1.E-3
        ISCAN=MOD(KGDS(11)/128,2)
        JSCAN=MOD(KGDS(11)/64,2)
        NSCAN=MOD(KGDS(11)/32,2)
        HI=(-1.)**ISCAN
        HJ=(-1.)**(1-JSCAN)
        DLON=HI*(MOD(HI*(RLON2-RLON1)-1+3600,360.)+1)/(IM-1)
        DLAT=(RLAT2-RLAT1)/(JM-1)
        XMIN=0
        XMAX=IM+1
        IF(IM.EQ.NINT(360/ABS(DLON))) XMAX=IM+2
        YMIN=0
        YMAX=JM+1
        NRET=0
        LROT=0
        FILL=-999

            IF(ABS(RLON).LE.360.AND.ABS(RLAT).LE.90) THEN
              XPTS=1+HI*MOD(HI*(RLON-RLON1)+3600,360.)/DLON
              YPTS=1+(RLAT-RLAT1)/DLAT

              IF(XPTS.GE.XMIN.AND.XPTS.LE.XMAX.AND.
     &           YPTS.GE.YMIN.AND.YPTS.LE.YMAX) THEN
                NRET=NRET+1
                IF(LROT.EQ.1) THEN
                  CROT=1
                  SROT=0
                ENDIF
              ELSE
                XPTS=FILL
                YPTS=FILL
              ENDIF
            ELSE
              XPTS=FILL
              YPTS=FILL
            ENDIF


        return
        end

C --------------------

        subroutine gen_filename(lat,lon,name)

        integer lat,lon
        character*7 name
        character*1 latswitch,lonswitch,type
        character*2 latval
        character*3 lonval

        if (lat .ge. 0) latswitch='N'
        if (lat .lt. 0) latswitch='S'

        if (lon .ge. 0) lonswitch='E'
        if (lon .lt. 0) lonswitch='W'


         write(latval,'(i2.2)') abs(lat)
         write(lonval,'(i3.3)') abs(lon)

        name=latval//latswitch//lonval//lonswitch

        end subroutine gen_filename

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


	subroutine greatcir(rlat1,rlon1,rlat2,rlon2,dist)

        parameter (pie=3.141592654,rconv=3.141592654/180.)
        parameter (erad=6.3712E+6)

C	formulation from Steers (1965)

        costerm=cos(rlat2*rconv)*cos(rlat1*rconv)
        absdelong=abs(rlon2-rlon1)*rconv
        havdelong=(1.-cos(absdelong))/2.
        colatf=rconv*(90-rlat1)
        colatt=rconv*(90-rlat2)
        factor=amax1(colatf,colatt)-amin1(colatf,colatt)
        havfactor=(1.-cos(factor))/2.
        havft=havdelong*costerm+havfactor
        angdist=acos(1-2*havft)/rconv
        distm=angdist*(2*pie*erad)/360.
        dist=distm/1000.

	end subroutine greatcir
