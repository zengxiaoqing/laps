      C_SAT_ID(1)='goes08'
      C_SAT_ID(2)='meteos'
      C_SAT_ID(3)='goes10'
      C_SAT_ID(4)='gmssat'

      do k=1,maxsat
      do j=1,maxtype
         C_SAT_TYPES(j,k)='   '
      do i=1,maxchannel
         C_CHANNEL_TYPES(i,j,k)='   '
      enddo
      enddo
      enddo

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
c ----
c goes08
      C_CHANNEL_TYPES(1,1,1)='vis'
      C_CHANNEL_TYPES(2,1,1)='4u '
      C_CHANNEL_TYPES(3,1,1)='wv '
      C_CHANNEL_TYPES(4,1,1)='11u'
      C_CHANNEL_TYPES(5,1,1)='12u'

      C_CHANNEL_TYPES(1,2,1)='vis'
      C_CHANNEL_TYPES(2,2,1)='i39'
      C_CHANNEL_TYPES(3,2,1)='iwv'
      C_CHANNEL_TYPES(4,2,1)='i11'
      C_CHANNEL_TYPES(5,2,1)='i12'

      C_CHANNEL_TYPES(1,3,1)='vis'
      C_CHANNEL_TYPES(2,3,1)='4u '
      C_CHANNEL_TYPES(3,3,1)='wv '
      C_CHANNEL_TYPES(4,3,1)='11u'
      C_CHANNEL_TYPES(5,3,1)='12u'

      C_CHANNEL_TYPES(1,4,1)='vis'
      C_CHANNEL_TYPES(2,4,1)='4u '
      C_CHANNEL_TYPES(3,4,1)='wv '
      C_CHANNEL_TYPES(4,4,1)='11u'
      C_CHANNEL_TYPES(5,4,1)='12u'  !<-- end 1st type
c
c goes09 --- No. 5-12-99 J. Smart changed to AFWA METEOSAT
      C_CHANNEL_TYPES(1,4,2)='vis'
      C_CHANNEL_TYPES(2,4,2)='4u '
      C_CHANNEL_TYPES(3,4,2)='wv '
      C_CHANNEL_TYPES(4,4,2)='11u'
      C_CHANNEL_TYPES(5,4,2)='12u'  !<-- end 2nd type
c goes10
      C_CHANNEL_TYPES(1,1,3)='vis'
      C_CHANNEL_TYPES(2,1,3)='4u '
      C_CHANNEL_TYPES(3,1,3)='wv '
      C_CHANNEL_TYPES(4,1,3)='11u'
      C_CHANNEL_TYPES(5,1,3)='12u'

      C_CHANNEL_TYPES(1,2,3)='vis'
      C_CHANNEL_TYPES(2,2,3)='i39'
      C_CHANNEL_TYPES(3,2,3)='iwv'
      C_CHANNEL_TYPES(4,2,3)='i11'
      C_CHANNEL_TYPES(5,2,3)='i12'

      C_CHANNEL_TYPES(1,4,3)='vis'
      C_CHANNEL_TYPES(2,4,3)='4u '
      C_CHANNEL_TYPES(3,4,3)='wv '
      C_CHANNEL_TYPES(4,4,3)='11u'
      C_CHANNEL_TYPES(5,4,3)='12u'
c gmssat
      C_CHANNEL_TYPES(1,4,4)='vis'
      C_CHANNEL_TYPES(3,4,4)='wvp'
      C_CHANNEL_TYPES(4,4,4)='11u'
      C_CHANNEL_TYPES(5,4,4)='12u'  !<-- end 4th type
c gmssat for taiwan (JS  and BS Wang 6-7-01)
      C_CHANNEL_TYPES(1,3,4)='vis'
      C_CHANNEL_TYPES(3,3,4)='wv '  !water vapor
      C_CHANNEL_TYPES(4,3,4)='11u'  !ir1 
      C_CHANNEL_TYPES(5,3,4)='12u'  !<-- end 4th type ... not known yet

c-----
      do j=1,maxsat
      do i=1,maxtype
         R_LATIN(i,j)=0.0
         R_LAP(i,j)=0.0
         R_LOV(i,j)=0.0
         I_EWCYCLES(i,j)=0
         I_EWINCS(i,j)=0
         I_NSCYCLES(i,j)=0
         I_NSINCS(i,j)=0
      enddo
      enddo
      R_LATIN(2,1)=25.00000
      R_LATIN(2,3)=25.00000
      R_LATIN(3,1)=25.00000 
      R_LAP(2,1)=25.00000
      R_LAP(2,3)=25.00000
      R_LAP(3,1)=25.00000
      R_LOV(2,1)=-95.00000
      R_LOV(2,3)=-95.00000
      R_LOV(3,1)=-95.00000
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
