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

	Program Gridmap
C**********************************************************************
c Portions of the following code were taken from the RAMS software and 
c   were used by permission of the *ASTER Division/Mission Research 
c   Corporation. 
C**********************************************************************
C*	This program will be used to map model grids based on the RAMS*
C*  version 2b which uses polar stereographic projections.            *
c
C*********************************************************************
        integer NX_L,NY_L

	call get_grid_dim_xy(NX_L,NY_L,istatus)
	if (istatus .ne. 1) then
           write (6,*) 'Error getting horizontal domain dimensions'
	   stop
	endif

        call Gridmap_sub(NX_L,NY_L)

        end
       
        subroutine Gridmap_sub(nnxp,nnyp)
c       include 'lapsparms.for'
        logical exist,new_DEM
        integer nnxp,nnyp
c       Parameter(nnxp=NX_L,nnyp=NY_L)
	Real mdlat,mdlon
	Real xmn(nnxp),ymn(nnyp)
	Real xtn(nnxp),ytn(nnyp)
	Real lat(nnxp,nnyp),lon(nnxp,nnyp)
        Real sw(2),nw(2),ne(2),se(2),pla,plo
        real  nboundary
        real  topt_30(nnxp,nnyp)
        real  topt_10(nnxp,nnyp)
        real  topt_out(nnxp,nnyp)
        real  topt_pctlfn(nnxp,nnyp)
        character*80 itoptfn_10,itoptfn_30,ipctlfn
        character*3 swt,twt
        character*6 c6_maproj

C*********************************************************************
c set nnxp,nnyp,mdlat,deltax,deltay, and ngrids
c       Data mdlat/39.32573/,mdlon/-104.24382/

c********************************************************************
c       Declarations for wrt_laps_static
        integer*4    ni,nj,nf
c       parameter (ni = NX_L)
c       parameter (nj = NY_L)
        parameter (nf = 4)
        
        character*3 var(nf) 
        data var /'LAT','LON','AVG','LDF'/

        character*125 comment(nf)
        character*131 model

        character*180 static_dir 
        integer len
        real*4 data(nnxp,nnyp,nf)
        real*4 zin_dum(nnxp,nnyp)
c       equivalence(data(1,1,1),lat)
c       equivalence(data(1,1,2),lon)
c       equivalence(data(1,1,3),topt_out)
c       equivalence(data(1,1,4),topt_pctlfn)

C*********************************************************************

        call get_laps_config('nest7grid',istatus)
        if(istatus .ne. 1)then
            write(6,*)' Bad status from get_laps_config'
            stop
        endif

        model = 'MODEL 4 delta x smoothed filter'

        comment(1) = 'Made from MODEL by J. Snook/ S. Albers 1-95'
        comment(2) = 'Made from MODEL by J. Snook/ S. Albers 1-95'

        icount_10 = 0
        icount_30 = 0
        icount_ramp = 0

        call get_directory('static',static_dir,len)

c ipltgrid is 1 if you want to plot the grid itself
c iplttopo is 1 if you want to plot the topography
c   the 30s topo covers the continental US
        itoptfn_30=static_dir(1:len)//'model/topo_30s/U'
c   the 10m topo covers the world
        itoptfn_10=static_dir(1:len)//'model/topo_10m/H'
c   the 10m pctl covers the world
        ipctlfn=static_dir(1:len)// 'model/land_10m/L'

        call get_topo_parms(silavwt_parm,toptwvl_parm,istatus)
	if (istatus .ne. 1) then
           write (6,*) 'Error getting terrain smoothing parms'
	   stop
	endif

!       Silhouette weighting parameter
        silavwt=silavwt_parm

!       Terrain wavelength for filtering
        toptwvl=toptwvl_parm

        ipltgrid=1
        iplttopo=1

C*********************************************************************

c calculate delta x and delta y using grid and map projection parameters
        call get_standard_latitudes(std_lat,std_lat2,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error calling laps routine'
            stop 
        endif
        write(6,*)' standard_lat = ',std_lat

        call get_grid_spacing(grid_spacing_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error calling laps routine'
            stop 
        endif
        write(6,*)' grid_spacing = ',grid_spacing_m

        call get_grid_center(mdlat,mdlon,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error calling laps routine'
            stop 
        endif
        write(6,*)' grid_center = ',mdlat,mdlon

        call get_c6_maproj(c6_maproj,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error calling laps routine'
            stop 
        endif
        write(6,*)' c6_maproj = ',c6_maproj

        call get_standard_longitude(std_lon,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error calling laps routine'
            stop 
        endif
        write(6,*)' std_lon = ',std_lon

        if(c6_maproj .eq. 'plrstr')then
            deltax = 2.0 / (1. + sind(std_lat)) * grid_spacing_m
            deltay = deltax
            if(std_lat2 .eq. +90.)then
                write(6,*)' Note, grid spacing will equal '
     1                    ,grid_spacing_m,' at a latitude of ',std_lat
                write(6,*)' deltax, deltay ',deltax,deltay
     1                   ,' at the north pole'
            else
                write(6,*)' This latitude is relative to where the pole'
     1                   ,' of the map projection is: lat/lon '
     1                   ,std_lat,std_lon
                write(6,*)' deltax, deltay ',deltax,deltay
     1                   ,' at the projection pole'
            endif
        else
            deltax = grid_spacing_m
            deltay = deltax
            write(6,*)' deltax, deltay ',deltax,deltay
        endif


c*********************************************************************

        erad=6367000.
        CALL POLAR_GP(mdlat,mdlon,XMN,YMN,DELTAX,DELTAY,
     1  NNXP,NNYP)
	DO 600 I=2,NNXP
	   XMN(I)=XMN(I-1)+DELTAX
 600	CONTINUE
        XMN(NNXP)=2*XMN(NNXP-1)-XMN(NNXP-2)
C
	DO 610 J=2,NNYP
	   YMN(J)=YMN(J-1)+DELTAY
 610	CONTINUE
        YMN(NNYP)=2*YMN(NNYP-1)-YMN(NNYP-2)
C
	   DO 650 I=2,NNXP
	      XTN(I)=.5*(XMN(I)+XMN(I-1))
 650	   CONTINUE
	   XTN(1)=1.5*XMN(1)-.5*XMN(2)

	   DO 660 J=2,NNYP
	      YTN(J)=.5*(YMN(J)+YMN(J-1))
 660	   CONTINUE
	   YTN(1)=1.5*YMN(1)-.5*YMN(2)

C*****************************************************************
C*  Convert it to lat/lon using the library routines.            *
        
           Do J = 1,nnyp
	      Do I = 1,nnxp
!	         call xytops(xtn(i),ytn(j),pla,plo,erad)
!                call pstoge(pla,plo,lat(I,J),lon(I,J),90.,std_lon)           

                 call xy_to_latlon(xtn(i),ytn(j),erad ! ,90.,std_lon
     1                                          ,lat(I,J),lon(I,J))

c              print *,'i,j,xtn,ytn,pla,lplo=',i,j,xtn,ytn,pla,plo
    	      enddo
           enddo
       print *,'lat,lon at 1,1 =',lat(1,1),lon(1,1)
       print *,'lat,lon at nnxp,nnyp =',lat(nnxp,nnyp),lon(nnxp,nnyp)
c
C*****************************************************************
c calculate topography
c
       if(iplttopo.eq.1)then

           call get_standard_longitude(std_lon,istatus)
           if(istatus .ne. 1)then
               write(6,*)' Error calling laps routine'
               stop 
           endif

           write(6,*)' Standard lon = ',std_lon

           write(6,*)
           write(6,*)' Processing 30s topo data....'
           CALL GEODAT(nnxp,nnyp,erad,90.,std_lon,xtn,ytn,
     +         deltax,deltay,TOPT_30,ITOPTFN_30,TOPTWVL,SILAVWT,new_DEM)

           if (.not.new_DEM) then
             write(6,*)
             write(6,*)' Processing 10m topo data....'
             CALL GEODAT(nnxp,nnyp,erad,90.,std_lon,xtn,ytn,
     +         deltax,deltay,TOPT_10,ITOPTFN_10,TOPTWVL,SILAVWT,new_DEM)
           endif

           write(6,*)
           write(6,*)' Processing 10m land data....'
           CALL GEODAT(nnxp,nnyp,erad,90.,std_lon,xtn,ytn,
     +        deltax,deltay,TOPT_PCTLFN,IPCTLFN,TOPTWVL,SILAVWT,new_DEM)

         if (.not.new_DEM) then

         do i = 1,nnxp
           do j = 1,nnyp

!              Select 30s or 10m topo data for this grid point (or a blend)

!              Check whether 30s data is missing or zero
               if(topt_30(i,j) .eq. 1e30 .or. topt_30(i,j) .eq. 0.
!    1                              .or.
!                  Are we in the Pittsburgh data hole?
!    1            (lat(i,j) .gt. 39.7 .and. lat(i,j) .lt. 41.3 .and.
!    1             lon(i,j) .gt.-79.3 .and. lon(i,j) .lt.-77.7)
     1                                                          )then 

!                  Use 10 min data
                   topt_out(i,j) = topt_10(i,j)
                   icount_10 = icount_10 + 1

               else ! Use 30s data, except ramp to 10m if near data boundary

!                  Determine the northern boundary of the 30s data at this lon
                   if    (lon(i,j).ge.-129..and.lon(i,j).le.-121.)then       
                       nboundary = 51.
                   elseif(lon(i,j).ge.-121..and.lon(i,j).le.-120.)then     
                       nboundary = 51. - lon(i,j) - (-121.)
                   elseif(lon(i,j).ge.-120..and.lon(i,j).le.-118.)then     
                       nboundary = 50.
                   elseif(lon(i,j).ge.-118..and.lon(i,j).le.-117.)then     
                       nboundary = 50. + lon(i,j) - (-118.)
                   elseif(lon(i,j).ge.-117..and.lon(i,j).le. -89.)then     
                       nboundary = 51.
                   elseif(lon(i,j).ge. -89..and.lon(i,j).le. -85.)then     
                       nboundary = 50.
                   elseif(lon(i,j).ge. -85..and.lon(i,j).le. -83.)then     
                       nboundary = 49.
                   elseif(lon(i,j).ge. -83..and.lon(i,j).le. -81.)then     
                       nboundary = 48.
                   elseif(lon(i,j).ge. -81..and.lon(i,j).le. -73.)then     
                       nboundary = 46.
                   elseif(lon(i,j).ge. -73..and.lon(i,j).le. -67.)then     
                       nboundary = 47.
                   elseif(lon(i,j).ge. -67.                      )then     
                       nboundary = 46.
                   else
                       nboundary = 51.
                   endif

                   alat1n = nboundary - 0.3
                   alat2n = nboundary - 0.1

!                  Determine the southern boundary of the 30s data at this lon
                   if(lon(i,j) .le. -103.)then         !        lon < -103
                       sboundary = 28. 
                   elseif(lon(i,j).ge.-103. .and. lon(i,j).le.-102.)then       
                       sboundary = 25. +  (-102. - lon(i,j)) * 3.
                   elseif(lon(i,j).ge.-102. .and. lon(i,j).le. -99.)then       
                       sboundary = 25.
                   elseif(lon(i,j).ge.-99.  .and. lon(i,j).le. -98.)then       
                       sboundary = 24. +  ( -98. - lon(i,j))
                   elseif(lon(i,j).ge.-98.                         )then       
                       sboundary = 24.
                   endif

                   alat1s = sboundary + 0.3
                   alat2s = sboundary + 0.1

!                  Decide whether to use 30s or 10m data (or a blend)

                   if  (  lat(i,j) .ge. alat2n)then    ! Use 10m data
                       topt_out(i,j) = topt_10(i,j)
                       icount_10 = icount_10 + 1

                   elseif(lat(i,j) .ge. alat1n .and. 
     1                    lat(i,j) .le. alat2n)then

!                      Between alat1n and alat2n,        Use weighted average

                       width = alat2n - alat1n
                       frac10 = (lat(i,j) - alat1n) / width
                       topt_out(i,j) = topt_10(i,j) * frac10 
     1                               + topt_30(i,j) * (1. - frac10)
                       icount_ramp = icount_ramp + 1

                       if(icount_ramp .eq. (icount_ramp/5) * 5 )then       
                           write(6,*)
                           write(6,*)'In blending zone, nboundary = '
     1                                       ,nboundary,alat1n,alat2n       
                           write(6,*)'lat/lon/frac =',lat(i,j),lon(i,j)
     1                                               ,frac10
                           write(6,*)'topt_30      =',topt_30(i,j)
                           write(6,*)'topt_10      =',topt_10(i,j)
                           write(6,*)'topt_out     =',topt_out(i,j)
                       endif

                   elseif(lat(i,j) .ge. alat1s .and. 
     1                    lat(i,j) .le. alat1n)then
                       topt_out(i,j) = topt_30(i,j)
                       icount_30 = icount_30 + 1       ! Use 30s data

                   elseif(lat(i,j) .ge. alat2s .and. 
     1                    lat(i,j) .le. alat1s)then

!                      Between alat1s and alat2s,        Use weighted average

                       width = alat1s - alat2s
                       frac10 = (alat1s - lat(i,j)) / width
                       topt_out(i,j) = topt_10(i,j) * frac10 
     1                               + topt_30(i,j) * (1. - frac10)
                       icount_ramp = icount_ramp + 1

                       if(icount_ramp .eq. (icount_ramp/5) * 5 )then       
                           write(6,*)
                           write(6,*)'In blending zone, sboundary = '
     1                                       ,sboundary,alat1s,alat2s       
                           write(6,*)'lat/lon/frac =',lat(i,j),lon(i,j)
     1                                               ,frac10
                           write(6,*)'topt_30      =',topt_30(i,j)
                           write(6,*)'topt_10      =',topt_10(i,j)
                           write(6,*)'topt_out     =',topt_out(i,j)
                       endif

                   elseif(lat(i,j) .le. alat2s)then    
                       topt_out(i,j) = topt_10(i,j)    ! Use 10m data
                       icount_10 = icount_10 + 1

                   else
                       write(6,*)' Software error in gridgen_model.f'
                       write(6,*)' lat/lon = ',lat(i,j),lon(i,j)
                       stop

                   endif ! Test to see if we blend the data

               endif

           enddo ! j
           enddo ! i
         else
           do j=1,nnyp
             do i=1,nnxp
               topt_out(i,j)=topt_30(i,j)
             enddo
           enddo
           icount_30=nnyp*nnxp
         endif
       endif

       write(6,*)
       print *,'topt_30    =',topt_30(1,1),topt_30(nnxp,nnyp)
       print *,'topt_10    =',topt_10(1,1),topt_10(nnxp,nnyp)
       print *,'topt_out   =',topt_out(1,1),topt_out(nnxp,nnyp)
       print *,'topt_pctlfn=',topt_pctlfn(1,1),topt_pctlfn(nnxp,nnyp)       
       print *,'# of grid pts using 30 sec data =  ',icount_30
       print *,'# of grid pts using 10 min data =  ',icount_10
       print *,'# of grid pts using blended data = ',icount_ramp
C
C*****************************************************************
c now make plot
c rbord is extra border width around grid to make plot look nice
c
	   rbord = 0.0
	   sw(1) = lat(1,1) - rbord
	   sw(2) = lon(1,1) - rbord
	   nw(1) = lat(1,nnyp) + rbord
	   nw(2) = lon(1,nnyp) - rbord
	   ne(1) = lat(nnxp,nnyp) + rbord
	   ne(2) = lon(nnxp,nnyp) + rbord
	   se(1) = lat(nnxp,1) - rbord
	   se(2) = lon(nnxp,1) + rbord

!	   Call opngks
!          call get_standard_longitude(std_lon,istatus)
!          if(istatus .ne. 1)then
!              write(6,*)' Error calling laps routine'
!              stop 
!          endif
!          Call Map(90.,std_lon,sw,nw,ne,se)
!          print *,'ipltgrid=',ipltgrid
!          if(ipltgrid.eq.1)then
!          print *,'ipltgrid=',ipltgrid
!          Call Grid(nnxp,nnyp,lat,lon)
!          endif
!          if(iplttopo.eq.1)then
!            cinc=100.
!            call conrec(topt_out,nnxp,nnxp,nnyp,
!    +          0.,0.,cinc,-1,0,0)
!           endif
c           call getset(a1,a2,a3,a4,b1,b2,b3,b4,c)
c           call set(a1,a2,a3,a4,b1,b2,b3,b4,c)
c           print *,a1,a2,a3,a4
c           print *,b1,b2,b3,b4
c           print *,c
c           call wtstr(0.2,0.2,'XXXXXXX=',20,0,-1)
c            call pwritx(800,800,'SILAVWT=',8,20,0,-1)
C

!          if(.true.)then
!            call frame
!            cinc=5.
!            call conrec(topt_pctlfn,nnxp,nnxp,nnyp,
!    +          0.,0.,cinc,-1,0,0)
!          endif

!          Call frame
!	   Call clsgks

        open(10,file=static_dir(1:len)//'latlon.dat'
     +         ,status='unknown',form='unformatted')
        open(11,file=static_dir(1:len)//'topo.dat'
     +         ,status='unknown',form='unformatted')
        open(15,file=static_dir(1:len)//'corners.dat'
     +         ,status='unknown')
        write(10)lat,lon
        write(11)topt_out
        write(12,*)topt_30
        write(13,*)topt_10
        write(14,*)topt_out
        write(15,*)lat(1,1),lon(1,1)
        write(15,*)lat(1,nnyp),lon(1,nnyp)
        write(15,*)lat(nnxp,1),lon(nnxp,1)
        write(15,*)lat(nnxp,nnyp),lon(nnxp,nnyp)
        write(16,*)topt_pctlfn
        close(10)
        close(11)
        close(15)

c SG97  topography.dat is written to be displayed with gnuplot
c SG97  this will make an elevated view from the SW direction over the domain
c SG97  start gnuplot and type in the following commands:
c SG97  set data style lines
c SG97  set view 30,330
c SG97  splot 'topography.dat'

	open(666,file=static_dir(1:len)//'topography.dat')
        do j=1,nnyp
	  do i=1,nnxp
	    write(666,*) topt_out(i,j)
          enddo
            write(666,'()')
        enddo
        close(666)
        call move(lat,data(1,1,1),nnxp,nnyp)            ! KWD
        call move(lon,data(1,1,2),nnxp,nnyp)            ! KWD
        call move(topt_out,data(1,1,3),nnxp,nnyp)       ! KWD
        call move(topt_pctlfn,data(1,1,4),nnxp,nnyp)    ! KWD

        INQUIRE(FILE=static_dir(1:len)//'nest7grid.cdl',EXIST=exist)
        if(.not.exist) then
	   print*,'Could not find file '
     +           ,static_dir(1:len)//'nest7grid.cdl '
           stop
	endif

        call check_domain(lat,lon,nnxp,nnyp,grid_spacing_m,5,istat_chk)

        write(6,*)'deltax = ',deltax

        write(6,*)'check_domain:status = ',istat_chk

        call put_laps_static(grid_spacing_m,model,comment,data
     1                                           ,nnxp,nnyp)


	Stop
	End

      SUBROUTINE GEODAT(n2,n3,erad,rlat,wlon1,xt,yt,deltax,deltay,
     1  DATR,OFN,WVLN,SILWT,which_data)
      implicit none
      integer n2,n3,n23,lb,mof,np,niq,njq,nx,ny,isbego,iwbego,
     1  iblksizo,no,iodim
      parameter (n23=20000)
      real vt3da(500),vt3db(n23)
      real vctr1(n23),
     1 vctr21(n23),erad,rlat,wlon1,deltax,deltay,wvln,silwt
      real DATR(N2,N3)
c      PARAMETER(IODIM=59000)
c SG97 iodim increased, to be able to read larger blocks of data
      PARAMETER(IODIM=5800000)
      real DATO(IODIM)
      real xt(N2),YT(N3),deltallo,deltaxq,deltayq,
     1  deltaxp,deltayp
      real std_lon
      integer istatus
      CHARACTER*80 OFN,TITLE
      logical which_data
C
c *********************
      nx=n2-1
      ny=n3-1
c ****************************
      LB=INDEX(OFN,' ')-1
      TITLE=OFN(1:LB)//'HEADER'
      LB=INDEX(TITLE,' ')-1
      CALL JCLGET(29,TITLE(1:LB),'FORMATTED',0)
      READ(29,2)IBLKSIZO,NO,ISBEGO,IWBEGO
 2    FORMAT(4I5)
      print *,'title=',title
      print *,'isbego,iwbego=',isbego,iwbego
      print *,'iblksizo,no=',iblksizo,no
      CLOSE(29)
      DELTALLO=FLOAT(IBLKSIZO)/FLOAT(NO-1)
      MOF=IODIM/(NO*NO)
c SG97 MOF determines the number of files held in buffer while reading
c SG97 DEM data; it saves some time when buffer data can be used instead
c SG97 of reading DEM file again. Originally MOF was 4.
      if (MOF.gt.10) MOF=5
      DELTAXQ=0.5*WVLN*DELTAX
      DELTAYQ=0.5*WVLN*DELTAY
      print *,'deltaxq,deltayq=',deltaxq,deltayq
      NP=MIN(10,MAX(1,INT(DELTAXQ/(DELTALLO*111000.))))
      print *,' np=',np
      DELTAXP=DELTAXQ/FLOAT(NP)
      DELTAYP=DELTAYQ/FLOAT(NP)
      NIQ=INT(FLOAT(NX)*DELTAX/DELTAXQ)+4
      NJQ=INT(FLOAT(NY)*DELTAY/DELTAYQ)+4
C
      call get_standard_longitude(std_lon,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Error calling laps routine'
          stop 
      endif
      CALL SFCOPQR(NO,MOF,NP,NIQ,NJQ,N2,N3,XT,YT,90.,std_lon,ERAD
     +            ,DELTALLO,DELTAXP,DELTAYP,DELTAXQ,DELTAYQ,IBLKSIZO
     +            ,ISBEGO,IWBEGO,DATO,VT3DA,VT3DB,DATR
     +            ,VCTR1,VCTR21,OFN,WVLN,SILWT,which_data)
      RETURN
      END


C
C     ******************************************************************
C
      SUBROUTINE SFCOPQR(NO,MOF,NP,NIQ,NJQ,N2,N3,XT,YT,RLAT,WLON1,ERAD
     +          ,DELTALLO,DELTAXP,DELTAYP,DELTAXQ,DELTAYQ,IBLKSIZO
     +          ,ISBEGO,IWBEGO,DATO,DATP,DATQ,DATR,ISO,IWO
     +          ,OFN,WVLN,SILWT,dem_data)
      real dato(no,no,mof)
      real DATP(NP,NP),DATQ(NIQ,NJQ),DATR(N2,N3)
      real ISO(MOF),IWO(MOF),XT(N2),YT(N3),rlat,wlon1,
     +    erad,deltallo,deltaxp,deltayp,deltaxq,deltayq,
     +   wvln,silwt,xq,yq,xp,yp,xcentr,ycentr,pla,plo,glatp,
     + glonp,rio,rjo,wio2,wio1,wjo2,wjo1,xq1,yq1
      real xr,yr,rval,sh,sha,rh,rha
      CHARACTER*80 OFN,TITLE3
      CHARACTER*3 TITLE1
      CHARACTER*4 TITLE2
      LOGICAL L1,L2,dem_data
C
      print *,'no,mof,np,niq,njq=',no,mof,np,niq,njq
c      stop

      NONO=NO*NO
      XCENTR=0.5*(XT(1)+XT(N2))
      YCENTR=0.5*(YT(1)+YT(N3))
      print *,xt(1),xt(n2),xcentr
      print *,'deltaxp=',deltaxp
      NOFR=0
      DO 11 IOF=1,MOF
         ISO(IOF)=0
         IWO(IOF)=0
  11    continue
      DO 15 JQ=1,NJQ
       print *,'jq,njq,niq=',jq,njq,niq
         DO 16 IQ=1,NIQ
            XQ=(FLOAT(IQ)-0.5*FLOAT(NIQ+1))*DELTAXQ+XCENTR
            YQ=(FLOAT(JQ)-0.5*FLOAT(NJQ+1))*DELTAYQ+YCENTR
            DO 17 JP=1,NP
               DO 18 IP=1,NP
                  XP=XQ+(FLOAT(IP)-0.5*FLOAT(NP+1))*DELTAXP
                  YP=YQ+(FLOAT(JP)-0.5*FLOAT(NP+1))*DELTAYP
!                 CALL XYTOPS(XP,YP,PLA,PLO,ERAD)
!                 CALL PSTOGE(PLA,PLO,GLATP,GLONP,rlat,wlon1)

c                 call xy_to_latlon(XP,YP,erad,rlat,wlon1,GLATP,GLONP) 
                  call xy_to_latlon(XP,YP,erad,GLATP,GLONP) 

c         print *,'rlat,wlon1=',rlat,wlon1
                  ISOC=(INT((GLATP-FLOAT(ISBEGO))/FLOAT(IBLKSIZO)+200.)
     +                -200)*IBLKSIZO+ISBEGO
                  IWOC=(INT((GLONP-FLOAT(IWBEGO))/FLOAT(IBLKSIZO)+400.)
     +                -400)*IBLKSIZO+IWBEGO
                  DO 19 IOFR=1,NOFR
                     JOFR=IOFR
                     IF(ISO(IOFR).EQ.ISOC.AND.IWO(IOFR).EQ.IWOC)GO TO 10
 19                 continue
                  ISOCPT=ABS(ISOC)/10
                  ISOCPO=ABS(ISOC)-ISOCPT*10
                  IWOCPH=ABS(IWOC)/100
                  IWOCPT=(ABS(IWOC)-IWOCPH*100)/10
                  IWOCPO=ABS(IWOC)-IWOCPH*100-IWOCPT*10
                  IF(ISOC.GE.0)THEN
                     WRITE(TITLE1,'(2I1,A1)')ISOCPT,ISOCPO,'N'
                  ELSE
                     WRITE(TITLE1,'(2I1,A1)')ISOCPT,ISOCPO,'S'
                  ENDIF
                  IF(IWOC.GE.0)THEN
                     WRITE(TITLE2,'(3I1,A1)')IWOCPH,IWOCPT,IWOCPO,'E'
                  ELSE
                     WRITE(TITLE2,'(3I1,A1)')IWOCPH,IWOCPT,IWOCPO,'W'
                  ENDIF
                  LB=INDEX(OFN,' ')-1
                  TITLE3=OFN(1:LB)//TITLE1//TITLE2
                  LB=INDEX(TITLE3,' ')-1
                  INQUIRE(FILE=TITLE3(1:LB),EXIST=L1,OPENED=L2)
                  IF(.NOT.L1)THEN
                     PRINT*, ' FILE',TITLE3(1:LB),' DOES NOT EXIST'
                     DATP(IP,JP)=0. ! set to missing?
                     GO TO 20
                  ENDIF

                  IF(NOFR.GE.MOF)THEN
                     DO 21 IOF=1,MOF
                        ISO(IOF)=0
                        IWO(IOF)=0
21                    continue
                     NOFR=0
                  ENDIF
                  NOFR=NOFR+1
                  JOFR=NOFR
                  len=index(ofn,' ')
                  if ((ofn(len-1:len).eq.'U').and.(no.eq.1200)) then 
                    CALL READ_DEM(29,TITLE3(1:LB),no,no,
     .              DATO(1,1,NOFR))
                    dem_data=.true.
                  else
                    CALL JCLGET(29,TITLE3(1:LB),'FORMATTED',0)
                    CALL VFIREC(29,DATO(1,1,NOFR),NONO,'LIN')
                    if ((ofn(len-1:len).eq.'U').and.(no.eq.121)) then
                      dem_data=.false.
                    endif
                  endif
c              print *,'nofr,dato=',nofr,dato(1,1,nofr)
                  CLOSE(29)
                  ISO(NOFR)=ISOC
                  IWO(NOFR)=IWOC
10		  continue
                  RIO=(GLONP-FLOAT(IWOC))/DELTALLO+1.
                  RJO=(GLATP-FLOAT(ISOC))/DELTALLO+1.
!                 Prevent Bounds Error (Steve Albers)
                  if(RIO .lt. 1.0)then
                      if(RIO .gt. 0.98)then
                          write(6,*)' Reset RIO for Machine Epsilon'      
                          RIO = 1.0
                      else
                          write(6,*)' ERROR: RIO out of bounds',RIO
                          stop
                      endif
                  endif

                  if(RJO .lt. 1.0)then
                      if(RJO .gt. 0.98)then
                          write(6,*)' Reset RJO for Machine Epsilon'      
                          write(6,*)JQ,IQ,
     1                          IP,JP,IO1,JO1,JOFR,RIO,RJO,GLATP,ISOC
                          RJO = 1.0
                      else
                          write(6,*)' ERROR: RJO out of bounds',RJO
                          write(6,*)JQ,IQ,
     1                          IP,JP,IO1,JO1,JOFR,RIO,RJO,GLATP,ISOC
                          stop
                      endif
                  endif

C
                  IO1=INT(RIO)
                  JO1=INT(RJO)
                  IO2=IO1+1
                  JO2=JO1+1
                  WIO2=RIO-FLOAT(IO1)
                  WJO2=RJO-FLOAT(JO1)
                  WIO1=1.0-WIO2
                  WJO1=1.0-WJO2
                  DATP(IP,JP)=WIO1*(WJO1*DATO(IO1,JO1,JOFR)
     +                             +WJO2*DATO(IO1,JO2,JOFR))
     +                       +WIO2*(WJO1*DATO(IO2,JO1,JOFR)
     +                             +WJO2*DATO(IO2,JO2,JOFR))
20               CONTINUE
18             continue
17           continue
C
            SHA=0.
            RHA=0.
            DO 22 JP=1,NP
               SH=0.
               RH=0.
               DO 23 IP=1,NP
!                 Test for missing - then go to 16?
                  SH=max(SH,DATP(IP,JP)) 
                  RH=RH+DATP(IP,JP)
23            continue
               SHA=SHA+SH/(2.*FLOAT(NP))
               RHA=RHA+RH
22           continue
            RHA=RHA/FLOAT(NP*NP)
            DO 24 IP=1,NP
               SH=0.
               DO 25 JP=1,NP
                  SH=max(SH,DATP(IP,JP))
25              continue
               SHA=SHA+SH/(2.*FLOAT(NP))
24           continue
            DATQ(IQ,JQ)=SHA*SILWT+RHA*(1.-SILWT)
c        print *,'datq=',datq(iq,jq)
16       continue
15     continue
       print *,'after 15'
c       stop
C
      XQ1=(1.-0.5*FLOAT(NIQ+1))*DELTAXQ+XCENTR
      YQ1=(1.-0.5*FLOAT(NJQ+1))*DELTAYQ+YCENTR
      print *,'datq=',datq(1,1),datq(niq,njq)
      DO 26 JR=1,N3
         DO 27 IR=1,N2
            XR=(XT(IR)-XQ1)/DELTAXQ+1.
            YR=(YT(JR)-YQ1)/DELTAYQ+1.
            CALL GDTOST(DATQ,NIQ,NJQ,XR,YR,RVAL)
            DATR(IR,JR)=max(0.,RVAL)
 27        continue
 26     continue
      print *,'datr=',datr(1,1),datr(N2,N3)
      RETURN
      END

      SUBROUTINE BINOM(X1,X2,X3,X4,Y1,Y2,Y3,Y4,XXX,YYY)
      implicit none
      real x1,x2,x3,x4,y1,y2,y3,y4,xxx,yyy,
     +   wt1,wt2,yz22,yz23,yz24,yz11,yz12,yz13,yoo
      integer istend
c      COMMON/BIN/ITYPP,I0X,I1X,I2X,YOO
       YYY=1E30
       IF(X2.GT.1.E19.OR.X3.GT.1.E19.OR.
     +   Y2.GT.1.E19.OR.Y3.GT.1.E19)RETURN
      WT1=(XXX-X3)/(X2-X3)
      WT2=1.0-WT1
      ISTEND=0
      IF(Y4.LT.1.E19.AND.X4.LT.1.E19) GO TO 410
      YZ22=WT1
      YZ23=WT2
      YZ24=0.0
      ISTEND= 1
410   IF(Y1.LT.1.E19.AND.X1.LT.1.E19) GO TO 430
      YZ11=0.0
      YZ12=WT1
      YZ13=WT2
      IF(ISTEND.EQ.1)GO TO 480
      GO TO 450
430   YZ11=(XXX-X2)*(XXX-X3)/((X1-X2)*(X1-X3))
      YZ12=(XXX-X1)*(XXX-X3)/((X2-X1)*(X2-X3))
      YZ13=(XXX-X1)*(XXX-X2)/((X3-X1)*(X3-X2))
      IF(ISTEND.EQ.  1    ) GO TO 470
450   YZ22=(XXX-X3)*(XXX-X4)/((X2-X3)*(X2-X4))
      YZ23=(XXX-X2)*(XXX-X4)/((X3-X2)*(X3-X4))
      YZ24=(XXX-X2)*(XXX-X3)/((X4-X2)*(X4-X3))
470   YYY=WT1*(YZ11*Y1+YZ12*Y2+YZ13*Y3)+WT2*(YZ22*Y2+YZ23*Y3+YZ24*Y4)
       GO TO 490
480      YYY=WT1*Y2+WT2*Y3
490   YOO=YYY
      RETURN
      END

      SUBROUTINE GDTOST(A,IX,IY,STAX,STAY,STAVAL)
*  SUBROUTINE TO RETURN STATIONS BACK-INTERPOLATED VALUES(STAVAL)
*  FROM UNIFORM GRID POINTS USING OVERLAPPING-QUADRATICS.
*  GRIDDED VALUES OF INPUT ARRAY A DIMENSIONED A(IX,IY),WHERE
*  IX=GRID POINTS IN X, IY = GRID POINTS IN Y .  STATION
*  LOCATION GIVEN IN TERMS OF GRID RELATIVE STATION X (STAX)
*  AND STATION COLUMN.
*  VALUES GREATER THAN 1.0E30 INDICATE MISSING DATA.
*
      real A(IX,IY),R(4),SCR(4),stax,stay,staval
     +  ,fixm2,fiym2,yy,xx
      IY1=INT(STAY)-1
      IY2=IY1+3
      IX1=INT(STAX)-1
      IX2=IX1+3
      STAVAL=1E30
      FIYM2=FLOAT(IY1)-1
      FIXM2=FLOAT(IX1)-1
      II=0
      DO 100 I=IX1,IX2
      II=II+1
      IF(I.GE.1.AND.I.LE.IX) GO TO 101
      SCR(II)=1E30
      GO TO 100
101   JJ=0
      DO 111 J=IY1,IY2
      JJ=JJ+1
      IF(J.GE.1.AND.J.LE.IY) GO TO 112
      R(JJ)=1E30
      GO TO 111
112   R(JJ)=A(I,J)
111   CONTINUE
      YY=STAY-FIYM2
      CALL BINOM(1.,2.,3.,4.,R(1),R(2),R(3),R(4),YY,SCR(II))
100   CONTINUE
      XX=STAX-FIXM2
      CALL BINOM(1.,2.,3.,4.,SCR(1),SCR(2),SCR(3),SCR(4),XX,STAVAL)
      RETURN
      END

       SUBROUTINE POLAR_GP(LAT,LON,X,Y,DX,DY,NX,NY)
C
       REAL*4 LAT,LON,X,Y,DX,DY,
     1 ERAD,TLAT,TLON,PLAT,PLON,
     1 XDIF,YDIF
C
       INTEGER*4 NX,NY
C
       RAD=3.141592654/180.
       ERAD=6367000.
       TLAT=90.0
       call get_standard_longitude(std_lon,istatus)
       if(istatus .ne. 1)then
           write(6,*)' Error calling laps routine'
           stop 
       endif
       TLON=std_lon
C
C      CALL GETOPS(PLAT,PLON,LAT,LON,TLAT,TLON)
C      CALL PSTOXY(XDIF,YDIF,PLAT,PLON,ERAD)

C      call latlon_to_xy(LAT,LON,TLAT,TLON,ERAD,XDIF,YDIF)
       call latlon_to_xy(LAT,LON,ERAD,XDIF,YDIF)

C
       X=XDIF+(1.-FLOAT(NX)/2.)*DX
       Y=YDIF+(1.-FLOAT(NY)/2.)*DY
C
       RETURN
C
       END
C +------------------------------------------------------------------+
      SUBROUTINE JCL
      CHARACTER*(*) FILENM,FORMT
      CHARACTER CFNAME*16,TEXTSTR*40
      LOGICAL EXANS

C     -------------------------------------------------------
      ENTRY JCLGET(IUNIT,FILENM,FORMT,IPRNT)
C
C         This routine access an existing file with the file name of
C           FILENM and assigns it unit number IUNIT.
C
      IF(IPRNT.EQ.1) THEN
      PRINT*,' Opening input unit ',IUNIT,' file name ',FILENM
      PRINT*,'         format  ',FORMT
      ENDIF

      OPEN(IUNIT,STATUS='OLD',FILE=FILENM,FORM=FORMT)

C
      RETURN
      END
c--------------------------------------------------------               
      subroutinevfirec(iunit,a,n,type)                                  
      character*1vc                                                     
      character*(*)type                                                 
      common/vform/vc(0:63)                                             
      characterline*80,cs*1                                             
      dimensiona(*)                                                     
                                                                        
      if(vc(0).ne.'0')callvfinit                                        
                                                                        
      ich0=ichar('0')                                                   
      ich9=ichar('9')                                                   
      ichcz=ichar('Z')                                                  
      ichlz=ichar('z')                                                  
      ichca=ichar('A')                                                  
      ichla=ichar('a')                                                  
                                                                        
      read(iunit,10)nn,nbits,bias,fact                                  
 10   format(2i8,2e20.10)                                               
      if(nn.ne.n)then                                                   
      print*,' Word count mismatch on vfirec record '                   
      print*,' Words on record - ',nn                                   
      print*,' Words expected  - ',n                                    
      stop'vfirec'                                                      
      endif                                                             
                                                                        
      nvalline=(78*6)/nbits                                             
      nchs=nbits/6                                                      
      do20i=1,n,nvalline                                                
      read(iunit,'(a78)', end=15)line                                          
      go to 16

 15   write(6,*)' Warning, incomplete terrain file detected'
      
 16   ic=0                                                              
      do30ii=i,i+nvalline-1                                             
      isval=0                                                           
      if(ii.gt.n)goto20                                                 
      do40iii=1,nchs                                                    
      ic=ic+1                                                           
      cs=line(ic:ic)                                                    
      ics=ichar(cs)                                                     
      if(ics.le.ich9)then                                               
      nc=ics-ich0                                                       
      elseif(ics.le.ichcz)then                                          
      nc=ics-ichca+10                                                   
      else                                                              
      nc=ics-ichla+36                                                   
      endif                                                             
      isval=intor(intlshft(nc,6*(nchs-iii)),isval)                      
 40   continue                                                          
      a(ii)=isval                                                       
 30   continue                                                          
 20   continue                                                          
                                                                        
      facti=1./fact                                                     
      if(type.eq.'LIN')then                                             
      do48i=1,n                                                         
      a(i)=a(i)*facti-bias                                              
 48   continue                                                          
      elseif(type.eq.'LOG')then                                         
      scfct=2.**(nbits-1)                                               
      do55i=1,n                                                         
      a(i)=sign(1.,a(i)-scfct)                                          
     +*(10.**(abs(20.*(a(i)/scfct-1.))-10.))                            
 55   continue                                                          
      endif                                                             
      return                                                  
      end                            

      subroutinevfinit                                                  
      character*1vc,vcscr(0:63)                                         
      common/vform/vc(0:63)                                             
      datavcscr/'0','1','2','3','4','5','6','7','8','9'                 
     +,'A','B','C','D','E','F','G','H','I','J'                          
     +,'K','L','M','N','O','P','Q','R','S','T'                          
     +,'U','V','W','X','Y','Z','a','b','c','d'                          
     +,'e','f','g','h','i','j','k','l','m','n'                          
     +,'o','p','q','r','s','t','u','v','w','x'                          
     +,'y','z','{','|'/                                                 
                                                                        
      do10n=0,63                                                        
      vc(n)=vcscr(n)                                                    
  10  continue                                                          
                                                                        
      return                                                            
      end

c ********************************************************************

      subroutine read_dem(unit_no,unit_name,nn1,nn2,data)
      implicit none
      integer countx,county,unit_no,nn1,nn2
      real data(nn1,nn2)
      integer*2 idata(nn1,nn2)
      logical l1,l2
      character*(*) unit_name

      open(unit_no,file=unit_name,status='old',access='direct',
     . recl=nn2*nn1*2)
      inquire(unit_no,exist=l1,opened=l2)
      read(unit_no,rec=1) idata
      do county=1,nn2
        do countx=1,nn1
          if (idata(countx,county).eq.-9999) idata(countx,county)=0
           data(countx,county)=float(idata(countx,nn2-county+1))
c SG97 initial data (DEM format) starts in the lower-left corner;
c SG97 this format is wrapped around to has upper-left corner as its start.
        enddo
      enddo
      close(unit_no)
      return
      end

C +------------------------------------------------------------------+
      FUNCTION INTLSHFT(IWORD,NSHFT)
C
C       This function shifts IWORD to the left NSHFT bits in a
C         circular manner.
C
      INTLSHFT=ISHFT(IWORD,NSHFT)
      RETURN
      END
C +------------------------------------------------------------------+
      FUNCTION INTOR(IWORD1,IWORD2)
C
C       This function performs a bit-by-bit OR between IWORD1 and
C         IWORD2.
C
      INTOR=IOR(IWORD1,IWORD2)
      RETURN
      END




