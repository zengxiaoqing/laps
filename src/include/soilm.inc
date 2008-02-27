C    Include file for Soilmoisture Routine
c Laps soilmoisture parameters
c----------------------------
      Integer DEPTH
      Real   HOR1,HOR2,HOR3
      Real   XDay,DCM
      Real   RainThres,SnowThres
      Real   Delt

c      Parameter(Imax=61,Jmax = 61)     ! Laps grid dimension
c      Parameter(HOR1=2., HOR2=6., HOR3=12.) ! Soil Horizons
c      Parameter(HOR1=2., HOR2=8., HOR3=20.) ! Soil Horizons

      Parameter(HOR1=6., HOR2=12., HOR3=36.) ! Soil Horizons
      Parameter(XDay = 1.0/24.0)

c      Parameter(DEPTH = 12, DCM = 30.48)
c      Parameter(DEPTH = 20, DCM = 50.80)

      Parameter(DEPTH = 36, DCM = 91.4)
      Parameter(RainThres = 0.075)    ! cm/hr ... like 0.03 in/hr
      Parameter(N=1)			!Number of rain periods
      Parameter(Snowthres=0.5)	!Snow cover threshold
c
C    End of include file
