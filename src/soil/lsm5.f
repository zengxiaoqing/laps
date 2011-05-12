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
      Program LSM5

      implicit none

      integer    istatus
      integer    nx_l,ny_l
       
      call get_grid_dim_xy(NX_L,NY_L,istatus)
      if(istatus.eq.1)then
        write(6,*)'LAPS Parameters obtained'
      else
        write(6,*)'IStatus = ',IStatus,'Error - Get_LAPS_Config'
        write(6,*)'Terminating LAPS-LSM5. Soil Moisture'
        stop
      end if
      call lsm5_sub(NX_L,NY_L)
      stop
      end

c
c ************************************************************************** 
c
      subroutine LSM5_sub(imax,jmax)

C     Program LAPS SoilMoisture Analysis
C     Changed to enable gridded analysis and real time
C     Chandran Subramaniam
C     2/8/93
      
C     John Smart 12/1/93: Adapting to run in real time on the
C     			  UNIX platform.  Set up LAPS standard I/O.
C                         Enhanced modularity.
C     J Smart    9/22/97  Adapt for dynamic array memory allocation
C
      integer   imax,jmax

      Include 'soilm.inc'
C
C**** MODEL12 is a Soil Moisture Content model developed in June l986, and
C     based upon the Remote Sensing Based Watershed Model developed at
C     the University of Maryland from 1982 to 1985.  This is Version 12.4,
C     as submitted to Water Resources Bulletin in January 1989.
C     Created by Groves.

      REAL 	KSAT,LAMDA,IN
      DATA 	DAY,SUMR,IN,OLDWEA/4*0./
      Integer   IStatus, i, j, ii
      Integer   Istatus_precip, istatus_m
      Integer   Istatus_w,Istatus_n,Istatus_e
      Integer   SoilType(Imax,Jmax)
      REAL      Laps_u(Imax,Jmax)     !u-component
      REAL      Laps_v(Imax,Jmax)     !v-component
      REAL      Laps_T(Imax,Jmax)     !Temperature
      REAL      Laps_TD(Imax,Jmax)    !Dewpoint Temperature
      REAL      Laps_IN(Imax,Jmax)    !Infiltration
      REAL      Laps_Wx(Imax,Jmax)    !Weather data; grid pt wet or dry.
      REAL      Laps_WFZ(Imax,Jmax)   !Wetting Front Depth, z
      REAL      Laps_MWF(Imax,Jmax)   !Wetting Front Moisture Content
      REAL      Laps_MWF_pre(Imax,Jmax)!Previous Wetting Front Moist Content.
      REAL      Laps_Evap(Imax,Jmax)  !Evaporation
      REAL      Laps_Rain(Imax, Jmax) !Radar-estimated 1hr liq. precip
      REAL      Laps_SC(Imax,Jmax)    !Snow Cover; fractional  0.0 to 1.0.
      REAL      Laps_SC_pre(Imax,Jmax)!Previous Hour Snow Cover   " .
c     REAL      Laps_SM(Imax,Jmax)    !Snow Melt, undefined currently
      REAL      Laps_SMC_3D(Imax,Jmax,3)!Three layer Soil Moisture Content
      REAL      soilm_field_cap(Imax,Jmax)!Field cap soil moistr (m**3/m**3)
      REAL      soilm_sat(Imax,Jmax) !Saturated soil moistr (m**3/m**3)
      REAL      data(Imax,Jmax,7)    !Holds LM2 output data - current time.
      REAL      data_s(Imax,Jmax,4)  !Holds LSX input data - current time

      Logical   Griddry,Filefound
      Integer   loop_bound
      Integer   i4time_smcur, i4time_smpre(25),
     &          lvl_s(4),lvl_1(3),lvl_2(7),lvl_l
      Character ftime_smcur*9
c laps precip
      Character ext_l*31, dir_l*150, var_l*3, lvl_coord_l*4,
     &          units_l*10, comment_l*125
c laps surface
      Character ext_s*31, dir_s*150, var_s(4)*3, lvl_coord_s(4)*4,
     &          units_s(4)*10, comment_s(4)*125
c background
      Character var_bkg(4)*3
c laps lm1...3-layer % soil saturation
      Character dir_1*150, ext1*31, var_1(3)*3, lvl_coord1(3)*4,
     &          units1(3)*10, comment1(3)*125
c laps lm2...variety of soil characteristic variables
      Character dir_2*150, ext2*31, var_2(7)*3, lvl_coord2(7)*4, 
     &          units2(7)*10, comment2(7)*125
      Character atime_smpre(25)*24
      Character fname*200

      DATA	ext_s/'lsx'/
      DATA	ext_l/'l1s'/
      DATA	ext1/'lm1'/
      DATA	ext2/'lm2'/
      DATA      var_s/'U  ','V  ','T  ','TD '/
      DATA      var_bkg/'USF','VSF','TSF','DSF'/
      DATA      var_l/'R01'/
      DATA      var_1/'LSM', 'LSM', 'LSM'/
      DATA      var_2/'CIV', 'DWF', 'MWF', 'WX ', 'EVP', 'SC ', 'SM'/
      DATA      units1/'%','%','%'/
      DATA      units2/'m**3','m','m**3/m**3',' ','m/s','%',' '/
      DATA      comment1/
     &          'layer 1 (0-6in [0-0.152m]) % soil saturation',
     &          'layer 2 (6-12in [.152-.305m]) % soil saturation',
     &          'layer 3 (12-36in [.305-0.914m]) % soil saturation'/
      DATA comment2(1)/'Cumulative Infiltration Volume (m)'/
      DATA comment2(2)/'Depth of Wetting Front (m)'/
      DATA comment2(3)/'Moisture content of wetting front (m**3/m**3)'/
      DATA comment2(4)/'Weather indicator - pos = wet, neg = dry'/
      DATA comment2(5)/'Evaporation (m/s)'/
      DATA comment2(6)/'Snow Cover Value (fractional, 0 to 1.0)'/
      DATA comment2(7)/'Snow Melt Value, (gm/s)'/
      DATA	Lvl_1/-1,-2,-3/
c 
c first read in current time
c
      call get_directory('etc',fname,len)
!     open(11, file=fname(1:len)//'systime.dat',status='unknown')
!     read(11,*,err=999)i4time_smcur
!     read(11,22,err=999)ftime_smcur
!22   format(1x,a9)
!     close(11)
      call get_systime(i4time_smcur,ftime_smcur,istatus)
c -------------------------
      print*,'systime: ',ftime_smcur,' i4time: ',i4time_smcur

C**** Read soil description and simulation time step 

      call get_directory('lm1',dir_1,len)

      Call Soil_In5(imax,jmax,SoilType,IStatus)   	! Get soil texture group
      if(IStatus.ne.1)then
         write(6,*)'Soil Textures Obtained'
         write(6,*)'Soil type = ',SoilType(1,1)
      else
         write(6,*)'Soil Textures Not Obtained - Terminating'
         stop
      end if
c
c  Get the current and previous i4time and ascii times. Allow for 4 previous
c  times for the previous soil moisture product.  Mandatory to have 1 previous
c  time so if the 4th previous time is the current time then the 5th is prev.

      loop_bound=5
      icnt=0
      do i=1,loop_bound
         i4time_smpre(i) = i4time_smcur - 3600*icnt
         call cv_i4tim_asc_lp(i4time_smpre(i),
     &                        atime_smpre(i),IStatus)
         icnt = icnt + 1
      end do
      print*,'current time: ',atime_smpre(1)

c  Get current surface data: LAPS LSX; u- v-component, temp and dew point.
       call get_directory('lsx',dir_s,len)
       do i=1,4
          lvl_s(i)=0
       end do

       Filefound = .false.
       ii = 0
       do while (.not.Filefound)
          ii = ii + 1
          Call Read_Laps_Data(i4time_smpre(ii), dir_s, ext_s,
     &              Imax, Jmax, 4, 4, var_s, lvl_s, lvl_coord_s, 
     &              units_s, comment_s, data_s, IStatus)

         if(IStatus.eq.1)then
            Write(6,*)'LAPS Surface data retrieved'
            write(6,*)'Sfc Directory [dir_s] = ',dir_s
            write(6,*)'time retrieved:',atime_smpre(ii)
            Filefound = .true.

            Laps_u=data_s(:,:,1)
            Laps_v=data_s(:,:,2)
            Laps_T=data_s(:,:,3)
            Laps_TD=data_s(:,:,4)

         elseif(ii .ge. loop_bound)then
            Filefound = .true.
            Write(6,*)'LAPS Surface data not available'
            print*,'Lets try model bkg for surface u/v/T/Td'
            call get_modelfg_2d(i4time_smcur,var_bkg(1)
     1,imax,jmax,Laps_u,istatus)
            if(istatus.ne.1)goto 59
            call get_modelfg_2d(i4time_smcur,var_bkg(2)
     1,imax,jmax,Laps_v,istatus)
            if(istatus.ne.1)goto 59
            call get_modelfg_2d(i4time_smcur,var_bkg(3)
     1,imax,jmax,Laps_T,istatus)
            if(istatus.ne.1)goto 59
            call get_modelfg_2d(i4time_smcur,var_bkg(4)
     1,imax,jmax,Laps_Td,istatus)

59          if(istatus.ne.1)then
               print*,'Failed to get background in get_modelfg_2d'
               print*,'Terminating LAPS Soil Moisture'
               return
            endif
         end if
       end do
c
c Get LAPS snow cover
c -------------------
       write(6,*)'Compute snow cover'
       write(6,*)

       Call readcsc(i4time_smcur
     &              ,imax,jmax
     &              ,LAPS_SC)

       write(6,*)'Snow cover computed'
       Write(6,*)'***************************************'
       Write(6,*)'Getting Field Cap and Saturated Soil Moisture'

       Do J = 1, Jmax
          Do I = 1, Imax
             ISOIL = SoilType(I,J)
C   Get default soil hydraulic parameters and initial soil moisture
             CALL SOILS(ISOIL,KSAT,THS,THR,PSIF,PSIAE,LAMDA)
             CALL AMC(Th_field_cap,ISOIL,THS)
             soilm_sat(i,j) = THS
             soilm_field_cap(i,j) = Th_field_cap
          Enddo
       Enddo
c
c Get previous soil moisture data: LAPS LM2; infiltration, depth of wetting
c front (WF), previous weather (ie, wet or dry), and moisture content of WF.
c
       Write(6,*) 'Obtaining Soil Moisture data'

       call get_directory('lm2',dir_2,len)
c       dir_2 = '../lapsprd/lm2/'
       do i=1,7
          lvl_2(i) = 0
       end do
c
       filefound = .false.
       ii=2
       do while(.not.filefound)
c
          call Read_Laps_Data(i4time_smpre(ii), dir_2, ext2,
     &       Imax, Jmax, 7, 7, var_2, lvl_2, lvl_coord2,
     &       units2, comment2, data, istatus)
c
          if(IStatus.eq.1)then
            Write(6,*)'LAPS lm2 data retrieved'
            write(6,*)'lm2 Directory [dir_2] = ',dir_2
            write(6,*)'time retrieved:',i4time_smpre(ii)
             do j = 1,jmax
             do i = 1,imax
                Laps_in(i,j) = data(i,j,1)
                Laps_wfz(i,j) = data(i,j,2)
                Laps_mwf(i,j) = data(i,j,3)
                Laps_wx(i,j) = data(i,j,4)
                Laps_Evap(i,j) = data(i,j,5)
                Laps_SC_pre(i,j) = data(i,j,6)
                Laps_mwf_pre(i,j) = data(i,j,7)
             end do
             end do
             filefound = .true.
c
c             Write(6,*)'Current Soil Moisture Product Retrieved'
c
          else
c
             if(ii.gt.14)then
c
c   Use Initial SM from Field Cap
c
             Write(6,*) '***************************************'
             Write(6,*) 'Using Default (Field Cap) Soil Moisture'
             Do J = 1, Jmax
             Do I = 1, Imax
                Laps_MWF(I,J) = soilm_field_cap(i,j)
                Laps_MWF_pre(I,J) = soilm_field_cap(i,j)
             Enddo
             Enddo
             filefound = .true.
c
             end if
c
             ii = ii + 1
c
          end if
       end do
c
c   Get Radar Estimated Precipitation.  LAPS L1S. Use "data_s" to hold.
c
       write(6,*)'Get LAPS precip and snow total'
       call get_directory('l1s',dir_l,len)
c       dir_l = '../lapsprd/l1s/'
       lvl_l = 0
       Call Read_Laps_Data(i4time_smcur, dir_l, ext_l,
     &           Imax, Jmax, 1, 1, var_l, lvl_l, lvl_coord_l, 
     &           units_l, comment_l, laps_rain, IStatus_precip)
       if(IStatus_precip.eq.1)then
c
          Write(6,*)
          Write(6,*)'Retrieved LAPS Rain data'
          write(6,*)'l1s Directory [dir_l] = ',dir_l
          write(6,*)'time retrieved:',i4time_smcur
          do j = 1,jmax
             do i = 1,imax
                laps_rain(i,j)=laps_rain(i,j)*100.    !cm/hr
             end do
          end do
          GridDry = .FALSE.
c
       else
c
          Write(6,*)'LAPS precip not available - assume dry'
          write(6,*)'l1s Directory [dir_l] = ',dir_l
          write(6,*)'time attempted:',i4time_smcur
          GridDry = .TRUE.

       end if
c
c This concludes getting the initial data for the soil moisture model
c ****************************************************************************
c ********************** Soil Moisture Subroutine ****************************
c
       Call Soil_Moisture(imax,jmax,Laps_u,Laps_v,
     &      Laps_T,Laps_TD,Laps_Rain,Laps_sc,
     &      Laps_IN,Laps_WFZ,Laps_MWF,Laps_MWF_pre,
     &      Laps_Wx,SoilType,GridDry,Laps_Evap,Laps_SMC_3D,
     &                  IStatus)
c
c Arrays Laps_IN, _WFZ, _MWF, _Wx return the new values of IN, WFZ, MWF and Wx
c from the above subroutine, thus they are both input and output.  Laps_SMC_3D
c is returned with current three layer soil moisture content.
c
       If(Istatus.eq.1)then
          Write(6,*)' Soil Moisture Computation complete'
          Write(6,*)' Check results for bad values'
c
          bad_lower_bndry = 0.0
          bad_upper_bndry = 100.0
          kmax=1
c
          call lsm_qc_check(imax,jmax,kmax,
     &                      laps_mwf,
     &                      bad_lower_bndry,
     &                      bad_upper_bndry,
     &                      istatus_m)

          if(istatus_m .lt. 0) write(6,910) istatus_m
c
          call lsm_qc_check(imax,jmax,kmax,
     &                      laps_wfz,
     &                      bad_lower_bndry,
     &                      bad_upper_bndry,
     &                      istatus_w)

          if(istatus_w .lt. 0) write(6,911) istatus_w
c
          call lsm_qc_check(imax,jmax,kmax,
     &                      laps_in,
     &                      bad_lower_bndry,
     &                      bad_upper_bndry,
     &                      istatus_n)

          if(istatus_n .lt. 0) write(6,912) istatus_n
c
          call lsm_qc_check(imax,jmax,kmax,
     &                      laps_evap,
     &                      bad_lower_bndry,
     &                      bad_upper_bndry,
     &                      istatus_e)

          if(istatus_e .lt. 0) write(6,913) istatus_e
c
          kmax=3
          call lsm_qc_check(imax,jmax,kmax,
     &                      laps_smc_3d,
     &                      bad_lower_bndry,
     &                      bad_upper_bndry,
     &                      istatus_3d)
c
          if(istatus_3d .lt. 0) write(6,914) istatus_3d
c
910       format(' +++ WARNING.  MWF:  istatus_m = ',i8,'+++')
911       format(' +++ WARNING.  WFZ:  istatus_w = ',i8,'+++')
912       format(' +++ WARNING.  CIV:  istatus_n = ',i8,'+++')
913       format(' +++ WARNING.  EVAP:  istatus_e = ',i8,'+++')
914       format(' +++ WARNING.  SMC-3D:  istatus_3d = ',i8,'+++')
c
          Write(6,*)' Writing current hour Soil Moisture for grid'
          do j=1,jmax
          do i=1,imax
             data(i,j,1)=Laps_IN(i,j)
             data(i,j,2)=Laps_WFZ(i,j)
             data(i,j,3)=Laps_MWF(i,j)
             data(i,j,4)=Laps_Wx(i,j)
             data(i,j,5)=Laps_Evap(i,j)
             data(i,j,6)=Laps_SC(i,j)
             data(i,j,7)=Laps_mwf_pre(i,j)
          end do
          end do
c
          Call Write_Laps_Data(i4Time_smcur, dir_1, ext1,
     &             Imax, Jmax,3, 3, var_1, lvl_1, lvl_coord1,
     &             units1, comment1, Laps_SMC_3D, IStatus)
          if(IStatus.eq.1)then
             write(6,*)'3D Soil Moisture Successfully Written'
          else
             write(6,*)'Error Writing 3D Soil Moisture - Laps_SMC_3D'
          end if
c
          Call Write_Laps_Data(I4Time_smcur, dir_2, ext2,
     &             Imax, Jmax,7, 7, var_2, lvl_2, lvl_coord2,
     &             units2, comment2, data, IStatus)
c
          if(IStatus.eq.1)then
             write(6,*)'2D Soil Moisture Fields Successfully Written'
          else
             write(6,*)'Error Writing 2D Soil Moisture Fields'
          end if
c
       else
c
          write(6,*)'Error during soil moisture computation'
c
       end if
       goto 100
c
 999   Write(6,*)'Error reading SYSTIME.DAT file'
       go to 100
 998   write(6,*)'error opening log file at:',ftime_smcur
c
100    WRITE(6,*)'End of Soil Moisture Simulation'
       close(6)
c
       STOP
       END
