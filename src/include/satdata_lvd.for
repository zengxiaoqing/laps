      C_SAT_ID(1)='goes08'
      C_SAT_ID(2)='meteos'
      C_SAT_ID(3)='goes10'
      C_SAT_ID(4)='gmssat'
      C_SAT_ID(5)='goes12'
      C_SAT_ID(6)='goes09'

      C_SAT_TYPES = '   '
      C_CHANNEL_TYPES = '   '

c     do k=1,maxsat
c     do j=1,maxtype
c        C_SAT_TYPES(j,k)='   '
c     do i=1,maxchannel
c        C_CHANNEL_TYPES(i,j,k)='   '
c     enddo
c     enddo
c     enddo

c first satellite (goes08)
      C_SAT_TYPES(1,1)='gvr'
      C_SAT_TYPES(2,1)='wfo'
      C_SAT_TYPES(3,1)='cdf'
      C_SAT_TYPES(4,1)='gwc'

c second satllite (meteosat)
      C_SAT_TYPES(4,2)='gwc'

c third satllite (goes10)
      C_SAT_TYPES(1,3)='gvr'
      C_SAT_TYPES(2,3)='wfo'

c forth satellite (gms)
      C_SAT_TYPES(3,4)='twn'
      C_SAT_TYPES(4,4)='gwc'
c fifth satellite (goes12)
      C_SAT_TYPES(1,5)='gvr'
      C_SAT_TYPES(2,5)='wfo'
      C_SAT_TYPES(3,5)='cdf'
      C_SAT_TYPES(4,5)='gwc'
c sixth satllite (goes09)
      C_SAT_TYPES(1,6)='gvr'

c ----
c goes08 (first satellite type)
c     format type 1 (gvr)
      C_CHANNEL_TYPES(1,1,1)='vis'
      C_CHANNEL_TYPES(2,1,1)='4u '
      C_CHANNEL_TYPES(3,1,1)='wv '
      C_CHANNEL_TYPES(4,1,1)='11u'
      C_CHANNEL_TYPES(5,1,1)='12u'

c     format type 2 (wfo)
      C_CHANNEL_TYPES(1,2,1)='vis'
      C_CHANNEL_TYPES(2,2,1)='i39'
      C_CHANNEL_TYPES(3,2,1)='iwv'
      C_CHANNEL_TYPES(4,2,1)='i11'
      C_CHANNEL_TYPES(5,2,1)='i12'

c     format type 3 (cdf)
      C_CHANNEL_TYPES(1,3,1)='vis'
      C_CHANNEL_TYPES(2,3,1)='4u '
      C_CHANNEL_TYPES(3,3,1)='wv '
      C_CHANNEL_TYPES(4,3,1)='11u'
      C_CHANNEL_TYPES(5,3,1)='12u'

c     format type 4 (gwc)
      C_CHANNEL_TYPES(1,4,1)='vis'
      C_CHANNEL_TYPES(2,4,1)='4u '
      C_CHANNEL_TYPES(3,4,1)='wv '
      C_CHANNEL_TYPES(4,4,1)='11u'
      C_CHANNEL_TYPES(5,4,1)='12u'  !<-- end 1st type
c
c meteos    #goes09 --- 5-12-99 J. Smart changed to AFWA METEOSAT
c     format type 4 (gwc)
      C_CHANNEL_TYPES(1,4,2)='vis'
      C_CHANNEL_TYPES(2,4,2)='4u '
      C_CHANNEL_TYPES(3,4,2)='wv '
      C_CHANNEL_TYPES(4,4,2)='11u'
      C_CHANNEL_TYPES(5,4,2)='12u'  !<-- end 2nd type
c
c goes10 (third satellite type)
c     format type 1 (gvr)
      C_CHANNEL_TYPES(1,1,3)='vis'
      C_CHANNEL_TYPES(2,1,3)='4u '
      C_CHANNEL_TYPES(3,1,3)='wv '
      C_CHANNEL_TYPES(4,1,3)='11u'
      C_CHANNEL_TYPES(5,1,3)='12u'

c     format type 2 (wfo)
      C_CHANNEL_TYPES(1,2,3)='vis'
      C_CHANNEL_TYPES(2,2,3)='i39'
      C_CHANNEL_TYPES(3,2,3)='iwv'
      C_CHANNEL_TYPES(4,2,3)='i11'
      C_CHANNEL_TYPES(5,2,3)='i12'

c     format type 4 (gwc)
      C_CHANNEL_TYPES(1,4,3)='vis'
      C_CHANNEL_TYPES(2,4,3)='4u '
      C_CHANNEL_TYPES(3,4,3)='wv '
      C_CHANNEL_TYPES(4,4,3)='11u'
      C_CHANNEL_TYPES(5,4,3)='12u'

c gmssat (fourth satellite type)
c     format type 4 (gwc)
      C_CHANNEL_TYPES(1,4,4)='vis'
      C_CHANNEL_TYPES(3,4,4)='wvp'
      C_CHANNEL_TYPES(4,4,4)='11u'
      C_CHANNEL_TYPES(5,4,4)='12u'  !<-- end 4th type
c     format type 3 (twn):  for taiwan (JS  and BS Wang 6-7-01)
      C_CHANNEL_TYPES(1,3,4)='vis'
      C_CHANNEL_TYPES(3,3,4)='wv '  !water vapor
      C_CHANNEL_TYPES(4,3,4)='11u'  !ir1 
      C_CHANNEL_TYPES(5,3,4)='12u'  !<-- end 4th type ... not known yet

c goes12 (fifth satellite type)
c     format type 1 (gvr)
      C_CHANNEL_TYPES(1,1,5)='vis'
      C_CHANNEL_TYPES(2,1,5)='4u '
      C_CHANNEL_TYPES(3,1,5)='wv '
      C_CHANNEL_TYPES(4,1,5)='11u'
      C_CHANNEL_TYPES(5,1,5)='12u'

c     format type 2 (wfo)
      C_CHANNEL_TYPES(1,2,5)='vis'
      C_CHANNEL_TYPES(2,2,5)='i39'
      C_CHANNEL_TYPES(3,2,5)='iwv'
      C_CHANNEL_TYPES(4,2,5)='i11'
      C_CHANNEL_TYPES(5,2,5)='i12'

c     format type 3 (cdf)
      C_CHANNEL_TYPES(1,3,5)='vis'
      C_CHANNEL_TYPES(2,3,5)='4u '
      C_CHANNEL_TYPES(3,3,5)='wv '
      C_CHANNEL_TYPES(4,3,5)='11u'
      C_CHANNEL_TYPES(5,3,5)='12u'

c goes09 (sixth satellite type)
c     format type 1 (gvr)
      C_CHANNEL_TYPES(1,1,6)='vis'
      C_CHANNEL_TYPES(2,1,6)='4u '
      C_CHANNEL_TYPES(3,1,6)='wv '
      C_CHANNEL_TYPES(4,1,6)='11u'
      C_CHANNEL_TYPES(5,1,6)='12u'

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

c     do j=1,maxsat
c     do i=1,maxtype
c        R_LATIN(i,j)=0.0
c        R_LAP(i,j)=0.0
c        R_LOV(i,j)=0.0
c        I_EWCYCLES(i,j)=0
c        I_EWINCS(i,j)=0
c        I_NSCYCLES(i,j)=0
c        I_NSINCS(i,j)=0
c     enddo
c     enddo

c these values come from the netcdf files for type gvr (GVAR)
c and satellites goes08 and goes10.  We'll also need them for
c GOES12 and GOES09

c type wfo and cdf for goes08, goes10, and goes12
      R_LATIN(2,1)=25.00000 !goes08/wfo
      R_LATIN(2,3)=25.00000 !goes10/wfo
      R_LATIN(2,5)=25.00000 !goes12/wfo
      R_LATIN(3,1)=25.00000 !goes08/cdf
      R_LATIN(3,5)=25.00000 !goes12/cdf
      R_LAP(2,1)=25.00000   !goes08/wfo
      R_LAP(2,3)=25.00000   !goes10/wfo
      R_LAP(2,5)=25.00000   !goes12/wfo
      R_LAP(3,1)=25.00000   !goes08/cdf
      R_LAP(3,5)=25.00000   !goes12/cdf
      R_LOV(2,1)=-95.00000  !goes08/wfo
      R_LOV(2,3)=-95.00000  !goes10/wfo
      R_LOV(2,5)=-95.00000  !goes12/wfo
      R_LOV(3,1)=-95.00000  !goes08/cdf

c type gvr and gwc for goes08 and goes10. (Soon add goes09 and goes12)
c                                          but only for gvr [no gwc])
      I_EWCYCLES(1,1)=2
      I_EWCYCLES(1,3)=2

      I_EWCYCLES(4,1)=2
      I_EWCYCLES(4,3)=2

      I_EWINCS(1,1)=3068
      I_EWINCS(1,3)=3025
      I_EWINCS(4,1)=3068
      I_EWINCS(4,3)=3068

      I_NSCYCLES(1,1)=4
      I_NSCYCLES(1,3)=4
      I_NSCYCLES(4,1)=4
      I_NSCYCLES(4,3)=4

      I_NSINCS(1,1)=3487
      I_NSINCS(1,3)=3299
      I_NSINCS(4,1)=3487
      I_NSINCS(4,3)=2795

      I_NWLINE_VIS(1,1)=3084
      I_NWLINE_VIS(1,3)=3060
      I_NWLINE_VIS(4,1)=2984
 
      I_NWLINE_IR(1,1)=3080
      I_NWLINE_IR(1,3)=3056
      I_NWLINE_IR(4,1)=2984

      I_NWLINE_WV(1,1)=3080
      I_NWLINE_WV(1,3)=3056
      I_NWLINE_WV(4,1)=2984
c ------------------------------
      I_NWPIX_VIS(1,1)=9050
      I_NWPIX_VIS(1,3)=10696
      I_NWPIX_VIS(4,1)=9248

      I_NWPIX_IR(1,1)=9050
      I_NWPIX_IR(1,3)=10696
      I_NWPIX_IR(4,1)=9248

      I_NWPIX_WV(1,1)=9050
      I_NWPIX_WV(1,1)=13500
      I_NWPIX_WV(1,1)=10696
      I_NWPIX_WV(1,1)=9248
