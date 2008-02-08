
      program pl_wideband
c
c * Read and plot netcdf wsr88d level 2 data
c   19 Nov 2001  Dongsoo Kim  Original 
c
      integer    rad_max, ang_max, nz_max
      parameter (rad_max = 999, ang_max = 999, nz_max = 16)
      character  YYDDDHHMM*9,FILENAME*150,ATYPE*7(16),XXXX*6(12)
      character  path*150
      integer    dim_dim, rad_dim(16), Z_dim, V_dim, VCP
c        rad_dim  !number of radials, slightly larger than 360
c        Z_dim    !number reflectivity bins in each radial
c        V_dim    !number wind bins in each radial
c        VCP      !Volume coverage pattern
      real  Refl(460,rad_max,16)
     &       ,Dwind(920,rad_max,16)
     &       ,SpecW(920,rad_max,16)
      real  A_ang(rad_max,16), E_ang(rad_max,16)
      integer   nsites
     &       ,vcpmode31(16),vcpmode11(16),vcpmode21(16)
      real    site_lat, site_lon, site_alt
      real    radar_lat(99), radar_lon(99), radar_alt(99)
      character radarname*5(99)
     &         ,sitename*132(99)

c * Working
      character CHMM*2, ch_ang*5, ch_alt*5, stname*4
      integer*2 image_Z(460, rad_max)
     &       ,image_V(920, rad_max)
     &       ,image_W(920, rad_max)
      real  Azim_ang(rad_max), Elev_ang(rad_max), resolV
      integer   colia(460,rad_max),ibin,vlvl,nlvl,ELELVL
      real      xloc(460,460), yloc(460,460)
      integer   colia_v(920,rad_max)
      real      uloc(920,920), vloc(920,920)
      real      deg2rad

c * Output:
      data ATYPE/'_elev01','_elev02','_elev03','_elev04','_elev05'
     &          ,'_elev06','_elev07','_elev08','_elev09','_elev10'
     &          ,'_elev11','_elev12','_elev13','_elev14','_elev15'
     &          ,'_elev16'/
      data  XXXX/'/kama/','/kcys/','/kddc/','/kftg/','/kfws/'
     &          ,'/kgld/','/kict/','/kinx/','/klbb/','/klzk/'
     &          ,'/ksgf/','/ktlx/'/
      data vcpmode31/1,3,5,6,8,0,0,0,0,0,0,0,0,0,0,0/
      data vcpmode11/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/
      data vcpmode21/1,3,5,6,7,8,9,10,11,0,0,0,0,0,0,0/
c          VCP = 11 !precipitation/severe weather mode
c * Graphics for reflectivity
      integer   LND2(16)
      character*3  LLB2(17)
      data LND2 /3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18/
      data LLB2 /' 0',' 0',' 5','10','15','20','25','30','35','40',
     &           '45','50','55','60','65','70','75'/
      character title*25, subtit1*25, subtit2*25
c * Graphics for Doppler winds
      integer   LND3(16)
      character*3  LLB3(17)
      data LND3 /3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18/
      data LLB3 /' 0',' 0',' 5','10','15','20','25','30','35','40',
     &           '45','50','55','60','65','70','75'/
c
      deg2rad = 3.1415926/180.d0
      print *,' Choose radar station (a4)'
      write(6,'(12a6)') (XXXX(k),k=1,12)
      read(5,'(a4)') stname
      
      print *,' Enter YYDDDHHMM in /public/data/radar/wsr88d/... '
      read(5,'(a9)') YYDDDHHMM
c
      ncount = 1
      DO NSITE = 1,1
         ncount = ncount + 1
         do 10 itry=1,1 ! 10,19
!           encode(2,'(i2)',CHMM) itry
            write(CHMM,1)itry
1           format(i2)

!           YYDDDHHMM(9:9) = CHMM(2:2)
!           XXXX(NSITE)(2:5) = stname

            if(stname .ne. 'RCWF')then
                FILENAME = '/public/data/radar/wsr88d/wideband/'
     1                 //stname//'/netcdf/'//YYDDDHHMM//ATYPE(1)
            else
                FILENAME = 
     1          '/data/lapb/import/lapsdat/tw/radar/RCWF/lvl2/'
     1                 //YYDDDHHMM//ATYPE(1)
            endif

            call get_wideband_hdr(FILENAME,dim_dim,Z_dim,V_dim
     &                   ,site_lat,site_lon,site_alt,VCP,istat)
            if(istat .eq. -1)  then
               goto 10
            elseif(istat .ne. -1) then
               goto 20
            endif
  10     continue
  20     continue
c
        radar_lat(ncount) = site_lat
        radar_lon(ncount) = site_lon
        radar_alt(ncount) = site_alt
         write(6,*) FILENAME
         print *,' Site name=', stname
         print *,' Site Lat=', radar_lat(ncount)
         print *,' Site Lon=', radar_lon(ncount)
         print *,' Site Alt=', radar_alt(ncount)
         print *,' Scan mode =', VCP

        vlvl = 0
        DO 30 n16=1,nz_max
         if(VCP .eq. 31 .or. VCP .eq. 32) nelev=vcpmode31(n16)
         if(VCP .eq. 11) nelev=vcpmode11(n16)
         if(VCP .eq. 21) nelev=vcpmode21(n16)
         if(nelev.eq. 0) goto 40

!        FILENAME =
c     &     'WIDEBAND'//XXXX(NSITE)//'netcdf/'//YYDDDHHMM//ATYPE(nelev)
!    &        'WIDEBAND/'//YYDDDHHMM//ATYPE(1)

         if(stname .ne. 'RCWF')then
             FILENAME = '/public/data/radar/wsr88d/wideband/'
     1                  //stname//'/netcdf/'//YYDDDHHMM//ATYPE(nelev)
         else
             FILENAME = 
     1          '/data/lapb/import/lapsdat/tw/radar/RCWF/lvl2/'
     1                 //YYDDDHHMM//ATYPE(nelev)
         endif

         call get_wideband_hdr(FILENAME,dim_dim,Z_dim,V_dim
     &                ,site_lat,site_lon,site_alt,VCP,istat)
         call get_wideband_netcdf(FILENAME,dim_dim,Z_dim,V_dim
     &         ,image_Z,image_V,image_W,Azim_ang,Elev_ang,resolV
     &         ,istat)
            if(istat .eq. -1) then
cdk                print *,' No file exist in this elev ', ATYPE(nelev)
                goto 30
            endif
            vlvl = vlvl + 1
            rad_dim(vlvl) = dim_dim
            do j=1,rad_dim(vlvl)
               A_ang(j,vlvl) = Azim_ang(j)
               E_ang(j,vlvl) = Elev_ang(j)
               do i=1,Z_dim
                  Refl(i,j,vlvl) = ((image_Z(i,j)-2.)/2.) - 32.
               enddo
               do i=1,V_dim
                  Dwind(i,j,vlvl) = (image_V(i,j)-129)*resolV
                  SpecW(i,j,vlvl) = image_W(i,j)
               enddo
            enddo
cdk        write(6,*) 'Azim_ang at Elev_ang ', Elev_ang(1) 
cdk        write(6,'(10f8.2)') (Azim_ang(irad),irad=1,rad_dim(vlvl))
cdk        write(6,*) 'image_Z at radial=100 '
        write(6,'(10f6.1)') (Dwind(100,j,vlvl),j=1,rad_dim(vlvl))
   30   CONTINUE
   40   continue
        nlvl = vlvl   !number of vertical levels
c * 
c * Convert (r,theta) |--> (x0,y0) per elevation angle
c * 
      print *,' Number of vertical levels is ', nlvl
cdk      ELELVL = 1 
cdk      print *,' Enter vertical level 1 -',nlvl
cdk      read(5,'(i2)') nlvl

        write(6,*)' Enter range of levels to plot '
        read(5,*)lvl_start,lvl_end

c
        do ELELVL=lvl_start,lvl_end
c
        print *,' plotting elevation level ',E_ang(1,ELELVL),' deg'      
!       encode(5,'(f5.2)',ch_ang) E_ang(1,ELELVL)
        write(ch_ang,2)E_ang(1,ELELVL)
 2      format(f5.2)

!       encode(5,'(i5)',ch_alt) int(radar_alt(ncount))
        write(ch_alt,3)int(radar_alt(ncount))
 3      format(i5)


c * First reflectivity
      write(6,*)' reflectivity'
      do j=1,rad_dim(ELELVL)
      do i=1,Z_dim
        xloc(i,j) = i*sin(A_ang(j,ELELVL)*deg2rad)/1200.
        yloc(i,j) = i*cos(A_ang(j,ELELVL)*deg2rad)/1200.
          if(Refl(i,j,ELELVL) .gt. 0. .and. 
     &       Refl(i,j,ELELVL).lt. 255. ) then
            colia(i,j) = ir_color(Refl(i,j,ELELVL))
          else
            colia(i,j) = 1
          endif
      enddo
      enddo

c * Second Dopple winds
      write(6,*)' velocity'
      do j=1,rad_dim(ELELVL)
      do i=1,V_dim
        uloc(i,j) = i*sin(A_ang(j,ELELVL)*deg2rad)/2400.
        vloc(i,j) = i*cos(A_ang(j,ELELVL)*deg2rad)/2400.
          if(Dwind(i,j,ELELVL) .gt. -10. .and. 
     &       Dwind(i,j,ELELVL).lt. 30. ) then
            colia_v(i,j) = dw_color(Dwind(i,j,ELELVL))
          else
            colia_v(i,j) = 1
          endif
      enddo
      enddo
c
      CALL OPNGKS
      call RADAR_COLOR
      title = 'WSR88D WIDEBAND AT '//stname
      subtit1 = 'ELEV ANGLE '//ch_ang
      subtit2 = 'SITE ELEV(m) '//ch_alt
  
      call set(0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,1)
      do j=1,365,1
      do i=10,Z_dim,2
         xo = xloc(i,j) + 0.5
         yo = yloc(i,j) + 0.5
         scale = 0.004*sqrt(2.35*i/rad_dim(ELELVL))
         call ngdots(xo,yo,1,scale,colia(i,j))
      enddo
      enddo
      call PLCHLQ(0.05,0.95,title,0.02,0.,-1.)
      call PLCHLQ(0.95,0.95,YYDDDHHMM,0.02,0.,1.)
      call PLCHLQ(0.05,0.92,subtit1,0.015,0.,-1.)
      call PLCHLQ(0.05,0.89,subtit2,0.015,0.,-1.)
      call LBLBAR(0,0.05,0.95,0.01,0.06,16,1.,0.4,LND2,0,LLB2,17,1)
      call frame
c
      if(.false.)then

      call set(0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,1)
      do j=1,365,1
      do i=10,V_dim,4
         xo = uloc(i,j) + 0.5
         yo = vloc(i,j) + 0.5
         scale = 0.003*sqrt(2.35*i/rad_dim(ELELVL))
         call ngdots(xo,yo,1,scale,colia_v(i,j))
      enddo
      enddo
      call PLCHLQ(0.05,0.95,title,0.02,0.,-1.)
      call PLCHLQ(0.95,0.95,YYDDDHHMM,0.02,0.,1.)
      call PLCHLQ(0.05,0.92,subtit1,0.015,0.,-1.)
      call PLCHLQ(0.05,0.89,subtit2,0.015,0.,-1.)
      call LBLBAR(0,0.05,0.95,0.01,0.06,16,1.,0.4,LND2,0,LLB2,17,1)
      call frame

      endif ! plot velocity
c
        enddo !ELELVL
c
      ENDDO ! nsite
      CALL CLSGKS

      stop
      end
c
      SUBROUTINE RADAR_COLOR
C
C  The color scheme is customized for NEXRAD image convention
c * See OFCM FM Handbook No.11 Part C A-7
c
      CALL GSCR (1,0,1.,1.,1.) !White
C
C    FOREGROUND COLORS
c
      CALL GSCR (1,1,0.,0.,0.) !Black

      CALL GSCR  (1,2, 0.00,0.00,0.00) !black
      CALL GSCR  (1,3, 0.61,0.61,0.61) !med gray
      CALL GSCR  (1,4, 0.46,0.46,0.46) !dk gray
      CALL GSCR  (1,5, 1.00,0.67,0.67) !lt pink
      CALL GSCR  (1,6, 0.93,0.55,0.55) !med pink
      CALL GSCR  (1,7, 0.79,0.44,0.44) !dk pink
      CALL GSCR  (1,8, 0.00,0.98,0.56) !lt green
      CALL GSCR  (1,9, 0.00,0.98,0.00) !med green
      CALL GSCR  (1,10,1.00,1.00,0.44) !lt yellow
      CALL GSCR  (1,11,0.82,0.82,0.38) !dark yellow
      CALL GSCR  (1,12,1.00,0.38,0.38) !lt red
      CALL GSCR  (1,13,0.85,0.00,0.00) !med red
      CALL GSCR  (1,14,0.70,0.00,0.00) !dk red
      CALL GSCR  (1,15,0.00,0.00,1.00) !blue
      CALL GSCR  (1,16,1.00,1.00,1.00) !white
      CALL GSCR  (1,17,0.91,0.00,1.00) !purple
      CALL GSCR  (1,18,0.91,0.00,1.00) !purple
      RETURN
      END
c
      integer function ir_color(echo)
      real echo

      ir_color = int(echo/5)+2

      return
      end
c
      integer function dw_color(echo)
      real echo

      dw_color = int(echo/2)+2

      return
      end
