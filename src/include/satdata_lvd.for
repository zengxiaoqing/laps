!     This file is included in subroutine 'config_satellite_lvd' in file
!     'src/lib/read_namelist.f'. Within the above subroutine these data are
!     aligned with the /lvd_namelist_cmn/ common block.      
      
!     C_SAT_ID(1)='goes08'
!     C_SAT_ID(2)='meteos'
!     C_SAT_ID(3)='goes10'
!     C_SAT_ID(4)='gmssat'
!     C_SAT_ID(5)='goes12'
!     C_SAT_ID(6)='goes09'
!     C_SAT_ID(7)='goes11'
!     C_SAT_ID(8)='noaapo'
!     C_SAT_ID(9)='mtsat'
!     C_SAT_ID(10)='fy'
!     C_SAT_ID(11)='goeswe'
!     C_SAT_ID(12)='goesea'

      C_SAT_TYPES = '   '
      C_CHANNEL_TYPES = '   '

c first satellite (goes08)
      C_SAT_TYPES(1,1)='gvr'
      C_SAT_TYPES(2,1)='wfo'
      C_SAT_TYPES(3,1)='cdf'
      C_SAT_TYPES(4,1)='rll'   !new type = raw lat lon: replaced gwc 5-17-07 (JS)

c second satllite (meteosat)
      C_SAT_TYPES(4,2)='rll'

c third satllite (goes10)
      C_SAT_TYPES(1,3)='gvr'
      C_SAT_TYPES(2,3)='wfo'

c forth satellite (gms)
      C_SAT_TYPES(2,4)='hko'
      C_SAT_TYPES(3,4)='twn'
      C_SAT_TYPES(4,4)='rll'
c fifth satellite (goes12)
      C_SAT_TYPES(1,5)='gvr'
      C_SAT_TYPES(2,5)='wfo'
      C_SAT_TYPES(3,5)='cdf'
      C_SAT_TYPES(4,5)='rll'
c sixth satllite (goes09)
      C_SAT_TYPES(1,6)='gvr'
      C_SAT_TYPES(3,6)='cdf'
c seventh satellite (goes11)
      C_SAT_TYPES(1,7)='gvr'
      C_SAT_TYPES(2,7)='wfo'
      C_SAT_TYPES(3,7)='cdf'
      C_SAT_TYPES(4,7)='rll'
c eighth satellite (noaa polar orbiter)
      C_SAT_TYPES(1,8)='ncp'
      C_SAT_TYPES(4,8)='rll'
c ninth satellite (mtsat)
      C_SAT_TYPES(4,9)='rll'
c
c ----
c goes08 (first satellite type)
c ===============================
c     format type 1 (gvr)
c     -------------------
      C_CHANNEL_TYPES(1,1,1)='vis'
      C_CHANNEL_TYPES(2,1,1)='4u '
      C_CHANNEL_TYPES(3,1,1)='wv '
      C_CHANNEL_TYPES(4,1,1)='11u'
      C_CHANNEL_TYPES(5,1,1)='12u'

c     format type 2 (wfo)
c     -------------------
      C_CHANNEL_TYPES(1,2,1)='vis'
      C_CHANNEL_TYPES(2,2,1)='i39'
      C_CHANNEL_TYPES(3,2,1)='iwv'
      C_CHANNEL_TYPES(4,2,1)='i11'
      C_CHANNEL_TYPES(5,2,1)='i12'

c     format type 3 (cdf)
c     -------------------
      C_CHANNEL_TYPES(1,3,1)='vis'
      C_CHANNEL_TYPES(2,3,1)='4u '
      C_CHANNEL_TYPES(3,3,1)='wv '
      C_CHANNEL_TYPES(4,3,1)='11u'
      C_CHANNEL_TYPES(5,3,1)='12u'

c     format type 4 (rll)
c     -------------------
      C_CHANNEL_TYPES(1,4,1)='vis'
      C_CHANNEL_TYPES(2,4,1)='4u '
      C_CHANNEL_TYPES(3,4,1)='wv '
      C_CHANNEL_TYPES(4,4,1)='11u'
      C_CHANNEL_TYPES(5,4,1)='12u'  !<-- end 1st type
c
c meteos    #goes09 --- 5-12-99 J. Smart changed to AFWA METEOSAT
c ===============================
c     format type 4 (rll)
c     -------------------
      C_CHANNEL_TYPES(1,4,2)='vis'
      C_CHANNEL_TYPES(2,4,2)='4u '
      C_CHANNEL_TYPES(3,4,2)='wv '
      C_CHANNEL_TYPES(4,4,2)='11u'
      C_CHANNEL_TYPES(5,4,2)='12u'  !<-- end 2nd type
c
c goes10 (third satellite type)
c ===============================
c     format type 1 (gvr)
c     -------------------
      C_CHANNEL_TYPES(1,1,3)='vis'
      C_CHANNEL_TYPES(2,1,3)='4u '
      C_CHANNEL_TYPES(3,1,3)='wv '
      C_CHANNEL_TYPES(4,1,3)='11u'
      C_CHANNEL_TYPES(5,1,3)='12u'

c     format type 2 (wfo)
c     -------------------
      C_CHANNEL_TYPES(1,2,3)='vis'
      C_CHANNEL_TYPES(2,2,3)='i39'
      C_CHANNEL_TYPES(3,2,3)='iwv'
      C_CHANNEL_TYPES(4,2,3)='i11'
      C_CHANNEL_TYPES(5,2,3)='i12'

c     format type 4 (rll)
c     -------------------
      C_CHANNEL_TYPES(1,4,3)='vis'
      C_CHANNEL_TYPES(2,4,3)='4u '
      C_CHANNEL_TYPES(3,4,3)='wv '
      C_CHANNEL_TYPES(4,4,3)='11u'
      C_CHANNEL_TYPES(5,4,3)='12u'

c gmssat (fourth satellite type)
c ===============================
c     format type 4 (rll)
c     -------------------
      C_CHANNEL_TYPES(1,4,4)='vis'
      C_CHANNEL_TYPES(3,4,4)='wvp'
      C_CHANNEL_TYPES(4,4,4)='11u'
      C_CHANNEL_TYPES(5,4,4)='12u'  
c     format type 2 (hko):  for HKO (JS  and PW Li 3-20-03)
c     -------------------
      C_CHANNEL_TYPES(1,2,4)='vis'
      C_CHANNEL_TYPES(3,2,4)='wv '  !water vapor
      C_CHANNEL_TYPES(4,2,4)='11u'  !ir1 
      C_CHANNEL_TYPES(5,2,4)='12u' 
c     format type 3 (twn):  for taiwan (JS  and BS Wang 6-7-01)
c     -------------------
      C_CHANNEL_TYPES(1,3,4)='vis'
      C_CHANNEL_TYPES(3,3,4)='wv '  !water vapor
      C_CHANNEL_TYPES(4,3,4)='11u'  !ir1 
      C_CHANNEL_TYPES(5,3,4)='12u'  !<-- end 4th type 

c goes12 (fifth satellite type)
c ===============================
c     format type 1 (gvr)
c     -------------------
      C_CHANNEL_TYPES(1,1,5)='vis'
      C_CHANNEL_TYPES(2,1,5)='4u '
      C_CHANNEL_TYPES(3,1,5)='wv '
      C_CHANNEL_TYPES(4,1,5)='11u'
      C_CHANNEL_TYPES(5,1,5)='   '  !doesn't appear to be 12u for goes12
      C_CHANNEL_TYPES(6,1,5)='13u'

c     format type 2 (wfo)
c     -------------------
      C_CHANNEL_TYPES(1,2,5)='vis'
      C_CHANNEL_TYPES(2,2,5)='i39'
      C_CHANNEL_TYPES(3,2,5)='iwv'
      C_CHANNEL_TYPES(4,2,5)='i11'
      C_CHANNEL_TYPES(5,2,5)='i12'

c     format type 3 (cdf)
c     -------------------
      C_CHANNEL_TYPES(1,3,5)='vis'
      C_CHANNEL_TYPES(2,3,5)='4u '
      C_CHANNEL_TYPES(3,3,5)='wv '
      C_CHANNEL_TYPES(4,3,5)='11u'
      C_CHANNEL_TYPES(5,3,5)='12u'

c goes09 (sixth satellite type)
c ===============================
c     format type 1 (gvr)
c     -------------------
      C_CHANNEL_TYPES(1,1,6)='vis'
      C_CHANNEL_TYPES(2,1,6)='4u '
      C_CHANNEL_TYPES(3,1,6)='wv '
      C_CHANNEL_TYPES(4,1,6)='11u'
      C_CHANNEL_TYPES(5,1,6)='12u'
c     format type 3 (cdf)
c     -------------------
      C_CHANNEL_TYPES(1,3,6)='vis'
      C_CHANNEL_TYPES(2,3,6)='4u '
      C_CHANNEL_TYPES(3,3,6)='wv '
      C_CHANNEL_TYPES(4,3,6)='11u'
      C_CHANNEL_TYPES(5,3,6)='12u'

c goes11 (seventh satellite type)
c ===============================
c     format type 1 (gvr)
c     -------------------
      C_CHANNEL_TYPES(1,1,7)='vis'
      C_CHANNEL_TYPES(2,1,7)='4u '
      C_CHANNEL_TYPES(3,1,7)='wv '
      C_CHANNEL_TYPES(4,1,7)='11u'
      C_CHANNEL_TYPES(5,1,7)='   '  !doesn't appear to be 12u for goes12
      C_CHANNEL_TYPES(6,1,7)='13u'

c     format type 2 (wfo)
c     -------------------
      C_CHANNEL_TYPES(1,2,7)='vis'
      C_CHANNEL_TYPES(2,2,7)='i39'
      C_CHANNEL_TYPES(3,2,7)='iwv'
      C_CHANNEL_TYPES(4,2,7)='i11'
      C_CHANNEL_TYPES(5,2,7)='i12'

c     format type 3 (cdf)
c     -------------------
      C_CHANNEL_TYPES(1,3,7)='vis'
      C_CHANNEL_TYPES(2,3,7)='4u '
      C_CHANNEL_TYPES(3,3,7)='wv '
      C_CHANNEL_TYPES(4,3,7)='11u'
      C_CHANNEL_TYPES(5,3,7)='12u'

c     format type 4 (cdf)
c     -------------------
      C_CHANNEL_TYPES(1,4,7)='vis'
      C_CHANNEL_TYPES(2,4,7)='4u '
      C_CHANNEL_TYPES(3,4,7)='wv '
      C_CHANNEL_TYPES(4,4,7)='11u'
      C_CHANNEL_TYPES(5,4,7)='12u'
c
c noaa polar orbiter (eighth satellite type)
c ===============================
c     format type 1 (ncp) -> only three channels are used.
c     -------------------
      C_CHANNEL_TYPES(1,1,6)='vis'
      C_CHANNEL_TYPES(2,1,6)='i39'
      C_CHANNEL_TYPES(4,1,6)='i11'
c     format type 4 (rll)
c     -------------------
      C_CHANNEL_TYPES(1,3,6)='vis'
      C_CHANNEL_TYPES(2,3,6)='4u '
      C_CHANNEL_TYPES(3,3,6)='wv '
      C_CHANNEL_TYPES(4,3,6)='11u'
      C_CHANNEL_TYPES(5,3,6)='12u'

c
c everything from here down should eventually be removed
c and obtained automatically from the input (netCDF) file.
c-----
      R_LATIN=0.0
      R_LAP=0.0
      R_LOV=0.0
      I_EWCYCLES=0
      I_EWINCS=0
      I_NSCYCLES=0
      I_NSINCS=0

c these values come from the netcdf files for type gvr (GVAR)
c and satellites goes08 and goes10.  We'll also need them for
c GOES12 and GOES09

c type wfo and cdf for goes08, goes10, and goes12
      R_LATIN(2,1)=25.00000 !goes08/wfo
      R_LATIN(2,3)=25.00000 !goes10/wfo
      R_LATIN(2,5)=25.00000 !goes12/wfo
      R_LATIN(3,1)=25.00000 !goes08/cdf
      R_LATIN(3,5)=25.00000 !goes12/cdf
      R_LATIN(3,6)=25.00000 !goes09/cdf
      R_LAP(2,1)=25.00000   !goes08/wfo
      R_LAP(2,3)=25.00000   !goes10/wfo
      R_LAP(2,5)=25.00000   !goes12/wfo
      R_LAP(3,1)=25.00000   !goes08/cdf
      R_LAP(3,5)=25.00000   !goes12/cdf
      R_LAP(:,17)=25.00000  !goes16/gnp
      R_LOV(2,1)=-95.00000  !goes08/wfo
      R_LOV(2,3)=-95.00000  !goes10/wfo
      R_LOV(2,5)=-95.00000  !goes12/wfo
      R_LOV(3,1)=-95.00000  !goes08/cdf
      R_LOV(3,6)=+120.00000  !goes08/cdf
      R_LOV(:,17)=-95.00000  !goes16/gnp

c type gvr and rll for goes08 and goes10. (Soon add goes09 and goes12)
c                                          but only for gvr [no rll])
      I_EWCYCLES(1,1)=2
      I_EWCYCLES(1,3)=2
      I_EWCYCLES(1,5)=2

      I_EWCYCLES(4,1)=2
      I_EWCYCLES(4,3)=2
      I_EWCYCLES(4,5)=2

      I_EWINCS(1,1)=3068
      I_EWINCS(1,3)=3025
      I_EWINCS(1,5)=3096
      I_EWINCS(4,1)=3068
      I_EWINCS(4,3)=3068
      I_EWINCS(4,5)=3096


      I_NSCYCLES(1,1)=4
      I_NSCYCLES(1,3)=4
      I_NSCYCLES(1,5)=4

      I_NSCYCLES(4,1)=4
      I_NSCYCLES(4,3)=4
      I_NSCYCLES(4,5)=4

      I_NSINCS(1,1)=3487
      I_NSINCS(1,3)=3299
      I_NSINCS(1,5)=2997
      I_NSINCS(4,1)=3487
      I_NSINCS(4,3)=2795
      I_NSINCS(4,5)=2997

      I_NWLINE_VIS(1,1)=3084
      I_NWLINE_VIS(1,3)=3060
      I_NWLINE_VIS(1,5)=2960
      I_NWLINE_VIS(4,1)=2984
 
      I_NWLINE_IR(1,1)=3080
      I_NWLINE_IR(1,3)=3056
      I_NWLINE_IR(1,5)=2960
      I_NWLINE_IR(4,1)=2984

      I_NWLINE_WV(1,1)=3080
      I_NWLINE_WV(1,3)=3056
      I_NWLINE_WV(1,5)=2960
      I_NWLINE_WV(4,1)=2984
c ------------------------------
      I_NWPIX_VIS(1,1)=9050
      I_NWPIX_VIS(1,3)=10696
      I_NWPIX_VIS(1,5)=9450
      I_NWPIX_VIS(4,1)=9248

      I_NWPIX_IR(1,1)=9050
      I_NWPIX_IR(1,3)=10696
      I_NWPIX_IR(1,5)=9450
      I_NWPIX_IR(4,1)=9248

      I_NWPIX_WV(1,1)=9050
      I_NWPIX_WV(1,3)=13500
      I_NWPIX_WV(1,5)=9450
      I_NWPIX_WV(4,1)=9248
c ------------------------------
c -- normalization values for all satellites except polar orbiter --
c -- Inputs to stretch for albedo calc --
      VIS_CNT_RANGE_IN(1,1)=0.00
      VIS_CNT_RANGE_IN(2,1)=303.57
      VIS_CNT_RANGE_IN(1,2)=0.00
      VIS_CNT_RANGE_IN(2,2)=0.00
      VIS_CNT_RANGE_IN(1,3)=0.00
      VIS_CNT_RANGE_IN(2,3)=305.00
      VIS_CNT_RANGE_IN(1,4)=38.00
      VIS_CNT_RANGE_IN(2,4)=400.00
      VIS_CNT_RANGE_IN(1,5)=0.00
      VIS_CNT_RANGE_IN(2,5)=364.00
      VIS_CNT_RANGE_IN(1,6)=0.00
      VIS_CNT_RANGE_IN(2,6)=305.00
      VIS_CNT_RANGE_IN(1,7)=0.00
      VIS_CNT_RANGE_IN(2,7)=364.00
      VIS_CNT_RANGE_IN(1,8)=0.00
      VIS_CNT_RANGE_IN(2,8)=303.57

      VIS_CNT_RANGE_OUT(1,1)=0.00
      VIS_CNT_RANGE_OUT(2,1)=255.00
      VIS_CNT_RANGE_OUT(1,2)=0.00
      VIS_CNT_RANGE_OUT(2,2)=0.00
      VIS_CNT_RANGE_OUT(1,3)=0.00
      VIS_CNT_RANGE_OUT(2,3)=255.00
      VIS_CNT_RANGE_OUT(1,4)=68.00
      VIS_CNT_RANGE_OUT(2,4)=220.00
      VIS_CNT_RANGE_OUT(1,5)=0.00
      VIS_CNT_RANGE_OUT(2,5)=255.00
      VIS_CNT_RANGE_OUT(1,6)=0.00
      VIS_CNT_RANGE_OUT(2,6)=255.00
      VIS_CNT_RANGE_OUT(1,7)=0.00
      VIS_CNT_RANGE_OUT(2,7)=255.00
      VIS_CNT_RANGE_OUT(1,8)=0.00
      VIS_CNT_RANGE_OUT(2,8)=255.00

c -- geostationary satellite height (m) above earth surface
!     SAT_RANGE_M(:)=42155680.00
      SAT_RANGE_M(:)=35786000.00
c -- sub lat/lon for each geostationary satellite (radians)
      R_SAT_SUB_LAT(1)=0.0
      R_SAT_SUB_LAT(2)=0.0
      R_SAT_SUB_LAT(3)=0.0
      R_SAT_SUB_LAT(4)=0.0
      R_SAT_SUB_LAT(5)=0.0
      R_SAT_SUB_LAT(6)=0.0
      R_SAT_SUB_LAT(7)=0.0
      R_SAT_SUB_LAT(8)=1.051
      R_SAT_SUB_LAT(9)=0.0
      R_SAT_SUB_LAT(17)=0.0
      R_SAT_SUB_LON(1)=-1.30900
      R_SAT_SUB_LON(2)= -2.35619
      R_SAT_SUB_LON(3)= -2.35619
      R_SAT_SUB_LON(4)= 2.44346
      R_SAT_SUB_LON(5)= -1.30900
      R_SAT_SUB_LON(6)= 2.530727
      R_SAT_SUB_LON(7)= -2.35619
      R_SAT_SUB_LON(8)= 0.42935
      R_SAT_SUB_LON(9)= 2.530727
      R_SAT_SUB_LON(17)= -1.5621
      R_SAT_SUB_LON(18)= 2.443
c -- these resolution (grid spacing) settings should be available in
c -- the input satellite data file; but not always so we hardwire here.
      R_RESOLUTION_X_VIS(1,1)=1264.2157
      R_RESOLUTION_X_VIS(2,1)=1015.9000
      R_RESOLUTION_X_VIS(3,1)=5079.2998
      R_RESOLUTION_X_VIS(4,1)=6607.7588
      R_RESOLUTION_X_VIS(1,2)=1075.8813
      R_RESOLUTION_X_VIS(2,2)=981.7122
      R_RESOLUTION_X_VIS(3,2)=0.0
      R_RESOLUTION_X_VIS(4,2)=6607.7588
      R_RESOLUTION_X_VIS(1,3)=1074.9694
      R_RESOLUTION_X_VIS(2,3)=0.0
      R_RESOLUTION_X_VIS(3,2)=0.0
      R_RESOLUTION_X_VIS(4,2)=0.0
      R_RESOLUTION_X_VIS(1,4)=0.0
      R_RESOLUTION_X_VIS(2,4)=0.0
      R_RESOLUTION_X_VIS(3,4)=0.0
      R_RESOLUTION_X_VIS(4,4)=0.0
      R_RESOLUTION_X_VIS(1,5)=1263.8979
      R_RESOLUTION_X_VIS(2,5)=0.0
      R_RESOLUTION_X_VIS(3,5)=0.0
      R_RESOLUTION_X_VIS(4,5)=0.0
      R_RESOLUTION_X_VIS(1,6)=0.0
      R_RESOLUTION_X_VIS(2,6)=0.0
      R_RESOLUTION_X_VIS(3,6)=0.0
      R_RESOLUTION_X_VIS(4,6)=0.0
      R_RESOLUTION_X_VIS(:,17)=507.94
      R_RESOLUTION_X_VIS(:,18)=5000.0
      R_RESOLUTION_Y_VIS(1,1)=1264.2157
      R_RESOLUTION_Y_VIS(2,1)=1015.9000
      R_RESOLUTION_Y_VIS(3,1)=5079.2998
      R_RESOLUTION_Y_VIS(4,1)=6607.7588
      R_RESOLUTION_Y_VIS(1,2)=1075.8813
      R_RESOLUTION_Y_VIS(2,2)=981.7125
      R_RESOLUTION_Y_VIS(3,2)=0.0
      R_RESOLUTION_Y_VIS(3,2)=6607.7588
      R_RESOLUTION_Y_VIS(1,3)=1074.9694
      R_RESOLUTION_Y_VIS(2,3)=0.0
      R_RESOLUTION_Y_VIS(3,3)=0.0
      R_RESOLUTION_Y_VIS(4,3)=0.0
      R_RESOLUTION_Y_VIS(1,4)=0.0
      R_RESOLUTION_Y_VIS(2,4)=0.0
      R_RESOLUTION_Y_VIS(3,4)=0.0
      R_RESOLUTION_Y_VIS(4,4)=10000.0
      R_RESOLUTION_Y_VIS(1,5)=1263.8979
      R_RESOLUTION_Y_VIS(2,5)=0.0
      R_RESOLUTION_Y_VIS(3,5)=0.0
      R_RESOLUTION_Y_VIS(4,5)=0.0
      R_RESOLUTION_Y_VIS(1,6)=0.0
      R_RESOLUTION_Y_VIS(2,6)=0.0
      R_RESOLUTION_Y_VIS(3,6)=0.0
      R_RESOLUTION_Y_VIS(4,6)=0.0
      R_RESOLUTION_Y_VIS(:,17)=507.94
      R_RESOLUTION_Y_VIS(:,18)=5000.0
      R_RESOLUTION_X_IR(1,1)=5053.0898
      R_RESOLUTION_X_IR(2,1)=4063.6001
      R_RESOLUTION_X_IR(3,1)=5079.2998
      R_RESOLUTION_X_IR(4,1)=6607.7588
      R_RESOLUTION_X_IR(1,2)=4302.6284
      R_RESOLUTION_X_IR(2,2)=3930.7061
      R_RESOLUTION_X_IR(3,2)=0.0
      R_RESOLUTION_X_IR(4,2)=6607.7588
      R_RESOLUTION_X_IR(1,3)=4298.6274
      R_RESOLUTION_X_IR(2,3)=0.0
      R_RESOLUTION_X_IR(3,3)=0.0
      R_RESOLUTION_X_IR(4,3)=0.0
      R_RESOLUTION_X_IR(1,4)=0.0
      R_RESOLUTION_X_IR(2,4)=0.0
      R_RESOLUTION_X_IR(3,4)=0.0
      R_RESOLUTION_X_IR(4,4)=10000.0
      R_RESOLUTION_X_IR(1,5)=5052.5850
      R_RESOLUTION_X_IR(2,5)=0.0
      R_RESOLUTION_X_IR(3,5)=0.0
      R_RESOLUTION_X_IR(4,5)=0.0
      R_RESOLUTION_X_IR(1,6)=0.0
      R_RESOLUTION_X_IR(2,6)=0.0
      R_RESOLUTION_X_IR(3,6)=0.0
      R_RESOLUTION_X_IR(4,6)=0.0
      R_RESOLUTION_X_IR(4,9)=8000.0
      R_RESOLUTION_X_IR(:,17)=2031.7625
      R_RESOLUTION_X_IR(:,18)=5000.0
      R_RESOLUTION_Y_IR(1,1)=5053.0898
      R_RESOLUTION_Y_IR(2,1)=4063.6001
      R_RESOLUTION_Y_IR(3,1)=5079.2998
      R_RESOLUTION_Y_IR(4,1)=6607.7588
      R_RESOLUTION_Y_IR(1,2)=4302.6284
      R_RESOLUTION_Y_IR(2,2)=3938.4451
      R_RESOLUTION_Y_IR(3,2)=0.0
      R_RESOLUTION_Y_IR(4,2)=6607.7588
      R_RESOLUTION_Y_IR(1,3)=4298.6274
      R_RESOLUTION_Y_IR(2,3)=0.0
      R_RESOLUTION_Y_IR(3,3)=0.0
      R_RESOLUTION_Y_IR(4,3)=0.0
      R_RESOLUTION_Y_IR(1,4)=0.0
      R_RESOLUTION_Y_IR(2,4)=0.0
      R_RESOLUTION_Y_IR(3,4)=0.0
      R_RESOLUTION_Y_IR(4,4)=10000.0
      R_RESOLUTION_Y_IR(1,5)=5052.5850
      R_RESOLUTION_Y_IR(2,5)=0.0
      R_RESOLUTION_Y_IR(3,5)=0.0
      R_RESOLUTION_Y_IR(4,5)=0.0
      R_RESOLUTION_Y_IR(1,6)=0.0
      R_RESOLUTION_Y_IR(2,6)=0.0
      R_RESOLUTION_Y_IR(3,6)=0.0
      R_RESOLUTION_Y_IR(4,6)=0.0
      R_RESOLUTION_Y_IR(4,9)=8000.0
      R_RESOLUTION_Y_IR(:,17)=2031.7625
      R_RESOLUTION_Y_IR(:,18)=5000.0
      R_RESOLUTION_X_WV(1,1)=8373.9775
      R_RESOLUTION_X_WV(2,1)=8127.2002
      R_RESOLUTION_X_WV(3,1)=5079.2998
      R_RESOLUTION_X_WV(4,1)=10126.4473
      R_RESOLUTION_X_WV(1,1)=7248.4502
      R_RESOLUTION_X_WV(2,2)=7876.8872
      R_RESOLUTION_X_WV(3,2)=0.0
      R_RESOLUTION_X_WV(4,2)=10126.4473
      R_RESOLUTION_X_WV(1,1)=7241.7622
      R_RESOLUTION_X_WV(2,3)=0.0
      R_RESOLUTION_X_WV(3,3)=0.0
      R_RESOLUTION_X_WV(4,3)=0.0
      R_RESOLUTION_X_WV(1,1)=0.0
      R_RESOLUTION_X_WV(1,4)=0.0
      R_RESOLUTION_X_WV(2,4)=0.0
      R_RESOLUTION_X_WV(3,4)=10000
      R_RESOLUTION_X_WV(4,4)=0.0
      R_RESOLUTION_X_WV(1,5)=8373.6602
      R_RESOLUTION_X_WV(2,5)=0.0
      R_RESOLUTION_X_WV(3,5)=0.0
      R_RESOLUTION_X_WV(4,5)=0.0
      R_RESOLUTION_X_WV(1,6)=0.0
      R_RESOLUTION_X_WV(2,6)=0.0
      R_RESOLUTION_X_WV(3,6)=0.0
      R_RESOLUTION_X_WV(4,6)=0.0
      R_RESOLUTION_X_WV(1:4,18)=5000.0
      R_RESOLUTION_Y_WV(1,1)=8373.9775
      R_RESOLUTION_Y_WV(2,1)=8127.2002
      R_RESOLUTION_Y_WV(3,1)=5079.2998
      R_RESOLUTION_Y_WV(4,1)=10126.4473
      R_RESOLUTION_Y_WV(1,2)=7248.4502
      R_RESOLUTION_Y_WV(2,2)=7876.8901
      R_RESOLUTION_Y_WV(3,2)=0.0
      R_RESOLUTION_Y_WV(4,2)=10126.4473
      R_RESOLUTION_Y_WV(1,3)=7241.7622
      R_RESOLUTION_Y_WV(2,3)=0.0
      R_RESOLUTION_Y_WV(3,3)=0.0
      R_RESOLUTION_Y_WV(4,3)=0.0
      R_RESOLUTION_Y_WV(1,4)=0.0
      R_RESOLUTION_Y_WV(2,4)=0.0
      R_RESOLUTION_Y_WV(3,4)=0.0
      R_RESOLUTION_Y_WV(4,4)=10000.0
      R_RESOLUTION_Y_WV(1,5)=8373.6602
      R_RESOLUTION_Y_WV(2,5)=0.0
      R_RESOLUTION_Y_WV(3,5)=0.0
      R_RESOLUTION_Y_WV(4,5)=0.0
      R_RESOLUTION_Y_WV(1,6)=0.0
      R_RESOLUTION_Y_WV(2,6)=0.0
      R_RESOLUTION_Y_WV(3,6)=0.0
      R_RESOLUTION_Y_WV(4,6)=0.0
      R_RESOLUTION_Y_WV(1:4,18)=5000.0
c -- navigation parms 1st lat/lon for each type
      R_LO1(1,1)=0.0
      R_LO1(2,1)=-116.48560
      R_LO1(3,1)=-126.13780
      R_LO1(4,1)=0.0
      R_LO1(1,2)=0.0
      R_LO1(2,2)=-120.47117
      R_LO1(3,2)=0.0
      R_LO1(4,2)=0.0
      R_LO1(1:4,3)=0.0
      R_LO1(1:4,4)=0.0
      R_LO1(1:4,5)=0.0
      R_LO1(1:4,6)=0.0
      R_LA1(1,1)=0.0
      R_LA1(2,1)=29.39183
      R_LA1(3,1)=16.28100
      R_LA1(4,1)=0.0
      R_LA1(1,2)=0.0
      R_LA1(2,2)=47.15661
      R_LA1(3,2)=0.0
      R_LA1(4,2)=0.0
      R_LA1(1:4,3)=0.0
      R_LA1(1:4,4)=0.0
      R_LA1(1:4,5)=0.0
      R_LA1(1:4,6)=0.0
c -- number of lines and pixels (elements) for each type
      N_PIXELS_VIS(1,1)=8001
      N_PIXELS_VIS(2,1)=2042
      N_PIXELS_VIS(3,1)=1201
      N_PIXELS_VIS(4,1)=1280
      N_PIXELS_VIS(1,2)=8013
      N_PIXELS_VIS(2,2)=2039
      N_PIXELS_VIS(3,2)=0
      N_PIXELS_VIS(4,2)=1280
      N_PIXELS_VIS(1,3)=8013
      N_PIXELS_VIS(2,3)=0
      N_PIXELS_VIS(3,3)=0
      N_PIXELS_VIS(4,3)=0
      N_PIXELS_VIS(1,4)=0
      N_PIXELS_VIS(2,4)=0
      N_PIXELS_VIS(3,4)=0
      N_PIXELS_VIS(4,4)=1
      N_PIXELS_VIS(1,5)=8001
      N_PIXELS_VIS(2,5)=0
      N_PIXELS_VIS(3,5)=0
      N_PIXELS_VIS(4,5)=1
      N_PIXELS_VIS(1,6)=0
      N_PIXELS_VIS(2,6)=0
      N_PIXELS_VIS(3,6)=0
      N_PIXELS_VIS(4,6)=1
      N_PIXELS_VIS(4,9)=1375
      N_PIXELS_VIS(1:4,17)=3072
      N_PIXELS_VIS(1:4,18)=1934
      N_LINES_VIS(1,1)=3290
      N_LINES_VIS(1,1)=2042
      N_LINES_VIS(1,1)=897
      N_LINES_VIS(1,1)=832
      N_LINES_VIS(1,2)=3944
      N_LINES_VIS(2,2)=2039
      N_LINES_VIS(3,2)=0
      N_LINES_VIS(4,2)=832
      N_LINES_VIS(1,3)=3944
      N_LINES_VIS(2,3)=0
      N_LINES_VIS(3,3)=0
      N_LINES_VIS(4,3)=0
      N_LINES_VIS(1,4)=0
      N_LINES_VIS(2,4)=0
      N_LINES_VIS(3,4)=0
      N_LINES_VIS(4,4)=1
      N_LINES_VIS(1,5)=3290
      N_LINES_VIS(2,5)=0
      N_LINES_VIS(3,5)=0
      N_LINES_VIS(4,5)=0
      N_LINES_VIS(1:4,6)=0
      N_LINES_VIS(1:4,9)=1375
      N_LINES_VIS(1:4,17)=2048
      N_LINES_VIS(1:4,18)=1544

      N_PIXELS_IR(1,1)=2000
      N_PIXELS_IR(2,1)=511
      N_PIXELS_IR(3,1)=1201
      N_PIXELS_IR(4,1)=1280
      N_PIXELS_IR(1,2)=2003
      N_PIXELS_IR(2,2)=510
      N_PIXELS_IR(3,2)=0
      N_PIXELS_IR(4,2)=1237
      N_PIXELS_IR(1,3)=2003
      N_PIXELS_IR(2:4,3)=0
      N_PIXELS_IR(1:4,4)=0
      N_PIXELS_IR(1,5)=2000
      N_PIXELS_IR(2:4,5)=0
      N_PIXELS_IR(1:4,6)=0
      N_PIXELS_IR(1:4,9)=1375
      N_PIXELS_IR(1:4,10)=1250 ! 1140
      N_PIXELS_IR(1:4,17)=2048
      N_PIXELS_IR(1:4,18)=1934
      N_LINES_IR(1,1)=823
      N_LINES_IR(2,1)=511
      N_LINES_IR(3,1)=897
      N_LINES_IR(4,1)=832
      N_LINES_IR(1,2)=986
      N_LINES_IR(2,2)=509
      N_LINES_IR(3,2)=0
      N_LINES_IR(4,2)=1237
      N_LINES_IR(1,3)=986
      N_LINES_IR(2:4,3)=0
      N_LINES_IR(1:4,4)=0
      N_LINES_IR(1,5)=823
      N_LINES_IR(2:4,5)=0
      N_LINES_IR(1:4,6)=0
      N_LINES_IR(1:4,9)=1375
      N_LINES_IR(1:4,10)=1250 ! 1140
      N_LINES_IR(1:4,17)=1024
      N_LINES_IR(1:4,18)=1544

      N_PIXELS_WV(1,1)=2000
      N_PIXELS_WV(2,1)=256
      N_PIXELS_WV(3,1)=1201
      N_PIXELS_WV(4,1)=1280
      N_PIXELS_WV(1,2)=2003
      N_PIXELS_WV(1,2)=255
      N_PIXELS_WV(1,2)=0
      N_PIXELS_WV(1,2)=1280
      N_PIXELS_WV(1,3)=2003
      N_PIXELS_WV(2:4,3)=0
      N_PIXELS_WV(1:4,4)=0
      N_PIXELS_WV(1,5)=2000
      N_PIXELS_WV(2:4,5)=0
      N_PIXELS_WV(1:4,6)=0
      N_PIXELS_WV(1:4,18)=1934
      N_LINES_WV(1,1)=412
      N_LINES_WV(2,1)=255
      N_LINES_WV(3,1)=897
      N_LINES_WV(4,1)=832
      N_LINES_WV(1,2)=493
      N_LINES_WV(2,2)=255
      N_LINES_WV(3,2)=0
      N_LINES_WV(4,2)=832
      N_LINES_WV(1,3)=493
      N_LINES_WV(2:4,3)=0
      N_LINES_WV(1:4,4)=0
      N_LINES_WV(1,5)=823
      N_LINES_WV(2:4,5)=0
      N_LINES_WV(1:4,6)=0
      N_LINES_WV(1:4,18)=1544
      IMC(1:6)=0
