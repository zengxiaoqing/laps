      subroutine rewrite_satellite_lvd_nl(istatus)
c
c This routine writes the satellite namelist ( static/satellite_lvd.nl),
c updating it with new values from the data files where necessary.
c
      implicit none

      integer    len_dir
      integer    istatus
      integer    lun
      integer    i,j,k
      integer    nc
      integer    slen

      character  nest7grid*150
      character  c_cell_afwa*7
      character  c_national*7

      include 'satellite_dims_lvd.inc'
      include 'satellite_common_lvd.inc'
      include 'satellite_namelist_lvd.cmn'
c
c start
c -----
      istatus=0
      call get_directory('nest7grid',nest7grid,len_dir)

      nest7grid = nest7grid(1:len_dir)//'/satellite_lvd.nl'

      lun=99
      open(lun,file=nest7grid,status='old',form='formatted',err=900)
      rewind(lun)

      if(l_national)then
         c_national='.TRUE.'
      else
         c_national='.FALSE.'
      endif
      if(l_cell_afwa)then
         c_cell_afwa='.TRUE.'
      else
         c_cell_afwa='.FALSE.'
      endif
      nc=0
      do k=1,maxsat
      do j=1,maxtype
      do i=1,maxchannel
         call s_len(path_to_raw_sat(i,j,k),slen)
         nc=max(slen,nc)
      enddo
      enddo
      enddo
      write(lun,*,err=950)' &satellite_lvd_nl'
      if(nc.eq.200)print*,'WARNING! Max string length = ',nc
      write(lun,1,err=950)(((path_to_raw_sat(i,j,k)(1:nc),
     &i=1,maxchannel),j=1,maxtype),k=1,maxsat)
c     write(lun,64)iflag_lvd_common
      write(lun,65)c_cell_afwa
      write(lun,66)c_national
      write(lun,61)(isats(i),i=1,maxsat)
      write(lun,62)((itypes(i,j),i=1,maxtype),j=1,maxsat)
      write(lun,63)(((ichannels(i,j,k),i=1,maxchannel),j=1,maxtype),
     &k=1,maxsat)
      write(lun,2,err=950)i_delta_sat_t_sec
      write(lun,3,err=950)((i_msng_sat_flag(i,j),i=1,maxtype),
     &j=1,maxsat)
      write(lun,5,err=950)(sat_range_m(i),i=1,maxsat)
      write(lun,7,err=950)n_images

c     write(lun,8,err=950)(c_sat_id(i),i=1,maxsat)
c     write(lun,9,err=950)((c_sat_types(i,j),i=1,maxtype),j=1,maxsat)
c     write(lun,10,err=950)(((c_channel_types(i,j,k),
c    &i=1,maxchannel),j=1,maxtype),k=1,maxsat)

      write(lun,11,err=950)(r_sat_sub_lat(i),i=1,maxsat)
      write(lun,12,err=950)(r_sat_sub_lon(i),i=1,maxsat)

      write(lun,13,err=950)((r_resolution_x_vis(i,j),i=1,maxtype),
     &j=1,maxsat)
      write(lun,14,err=950)((r_resolution_y_vis(i,j),i=1,maxtype),
     &j=1,maxsat)
      write(lun,15,err=950)((r_resolution_x_ir(i,j),i=1,maxtype),
     &j=1,maxsat)
      write(lun,16,err=950)((r_resolution_y_ir(i,j),i=1,maxtype),
     &j=1,maxsat)
      write(lun,17,err=950)((r_resolution_x_wv(i,j),i=1,maxtype),
     &j=1,maxsat)
      write(lun,18,err=950)((r_resolution_y_wv(i,j),i=1,maxtype),
     &j=1,maxsat)
      write(lun,20,err=950)((r_lo1(i,j),i=1,maxtype),j=1,maxsat)
      write(lun,21,err=950)((r_la1(i,j),i=1,maxtype),j=1,maxsat)

c     write(lun,19,err=950)((r_latin(i,j),i=1,maxtype),j=1,maxsat)
c     write(lun,22,err=950)((r_lap(i,j),i=1,maxtype),j=1,maxsat)
c     write(lun,23,err=950)((r_lov(i,j),i=1,maxtype),j=1,maxsat)

      write(lun,30,err=950)((i_start_vis(i,j),i=1,maxtype),
     &j=1,maxsat)
      write(lun,31,err=950)((i_end_vis(i,j),i=1,maxtype),
     &j=1,maxsat)
      write(lun,32,err=950)((j_start_vis(i,j),i=1,maxtype),
     &j=1,maxsat)
      write(lun,33,err=950)((j_end_vis(i,j),i=1,maxtype),
     &j=1,maxsat)
      write(lun,34,err=950)((i_start_ir(i,j),i=1,maxtype),
     &j=1,maxsat)
      write(lun,35,err=950)((i_end_ir(i,j),i=1,maxtype),
     &j=1,maxsat)
      write(lun,36,err=950)((j_start_ir(i,j),i=1,maxtype),
     &j=1,maxsat)
      write(lun,37,err=950)((j_end_ir(i,j),i=1,maxtype),
     &j=1,maxsat)
      write(lun,38,err=950)((i_start_wv(i,j),i=1,maxtype),
     &j=1,maxsat)
      write(lun,39,err=950)((i_end_wv(i,j),i=1,maxtype),
     &j=1,maxsat)
      write(lun,40,err=950)((j_start_wv(i,j),i=1,maxtype),
     &j=1,maxsat)
      write(lun,41,err=950)((j_end_wv(i,j),i=1,maxtype),
     &j=1,maxsat)

      write(lun,42,err=950)((n_pixels_vis(i,j),i=1,maxtype),
     &j=1,maxsat)
      write(lun,43,err=950)((n_lines_vis(i,j),i=1,maxtype),
     &j=1,maxsat)
      write(lun,44,err=950)((n_pixels_ir(i,j),i=1,maxtype),
     &j=1,maxsat)
      write(lun,45,err=950)((n_lines_ir(i,j),i=1,maxtype),
     &j=1,maxsat)
      write(lun,46,err=950)((n_pixels_wv(i,j),i=1,maxtype),
     &j=1,maxsat)
      write(lun,47,err=950)((n_lines_wv(i,j),i=1,maxtype),
     &j=1,maxsat)

c     write(lun,48,err=950)((i_ewCycles(i,j),i=1,maxtype),
c    &j=1,maxsat)
c     write(lun,49,err=950)((i_ewIncs(i,j),i=1,maxtype),
c    &j=1,maxsat)
c     write(lun,50,err=950)((i_nsCycles(i,j),i=1,maxtype),
c    &j=1,maxsat)
c     write(lun,51,err=950)((i_nsIncs(i,j),i=1,maxtype),
c    &j=1,maxsat)
c     write(lun,52,err=950)((i_nwline_vis(i,j),i=1,maxtype),
c    &j=1,maxsat)
c     write(lun,53,err=950)((i_nwline_ir(i,j),i=1,maxtype),
c    &j=1,maxsat)
c     write(lun,54,err=950)((i_nwline_wv(i,j),i=1,maxtype),
c    &j=1,maxsat)
c     write(lun,55,err=950)((i_nwpix_vis(i,j),i=1,maxtype),
c    &j=1,maxsat)
c     write(lun,56,err=950)((i_nwpix_ir(i,j),i=1,maxtype),
c    &j=1,maxsat)
c     write(lun,57,err=950)((i_nwpix_wv(i,j),i=1,maxtype),
c    &j=1,maxsat)

      write(lun,58,err=950)(imc(j),j=1,maxsat)
      write(lun,59,err=950)


1     format(1x,'PATH_TO_RAW_SAT=',5("'",a,"',",/))   !channel x type x sat
61    format(1x,'ISATS=',4(i1,","))
62    format(1x,'ITYPES=',4(i1,","),1x)
63    format(1x,'ICHANNELS=',5(i1,","),1x)
64    format(1x,'IFLAG_LVD_COMMON= ',i1,",")
65    format(1x,'L_CELL_AFWA= ',a,",")
66    format(1x,'L_NATIONAL= ',a,",")
2     format(1x,'I_DELTA_SAT_T_SEC=',i6,",")
3     format(1x,'I_MSNG_SAT_FLAG=',4(i4,","))
5     format(1x,'SAT_RANGE_M=',4(1x,f11.2,','))       !maxsat = 4
7     format(1x,'N_IMAGES=',i2,',')
8     format(1x,'C_SAT_ID=',4("'",a,"',"))            !maxsat = 4
9     format(1x,'C_SAT_TYPES=',4("'",a,"',"))
10    format(1x,'C_CHANNEL_TYPES=',5("'",a,"',"))     !channel x type x sat
11    format(1x,'R_SAT_SUB_LAT=',4(f10.5,","))        !all the rest are maxtype
12    format(1x,'R_SAT_SUB_LON=',4(f10.5,","))
13    format(1x,'R_RESOLUTION_X_VIS=',4(f10.4,','))
14    format(1x,'R_RESOLUTION_Y_VIS=',4(f10.4,','))
15    format(1x,'R_RESOLUTION_X_IR=',4(f10.4,','))
16    format(1x,'R_RESOLUTION_Y_IR=',4(f10.4,','))
17    format(1x,'R_RESOLUTION_X_WV=',4(f10.4,','))
18    format(1x,'R_RESOLUTION_Y_WV=',4(f10.4,','))
19    format(1x,'R_LATIN=',4(f10.5,','))
20    format(1x,'R_LO1=',4(f10.5,','))
21    format(1x,'R_LA1=',4(f10.5,','))
22    format(1x,'R_LAP=',4(f10.5,','))
23    format(1x,'R_LOV=',4(f10.5,','))

30    format(1x,'I_START_VIS=',4(i5,','))
31    format(1x,'I_END_VIS=',4(i5,','))
32    format(1x,'J_START_VIS=',4(i5,','))
33    format(1x,'J_END_VIS=',4(i5,','))
34    format(1x,'I_START_IR=',4(i5,','))
35    format(1x,'I_END_IR=',4(i5,','))
36    format(1x,'J_START_IR=',4(i5,','))
37    format(1x,'J_END_IR=',4(i5,','))
38    format(1x,'I_START_WV=',4(i5,','))
39    format(1x,'I_END_WV=',4(i5,','))
40    format(1x,'J_START_WV=',4(i5,','))
41    format(1x,'J_END_WV=',4(i5,','))
42    format(1x,'N_PIXELS_VIS=',4(i5,','))
43    format(1x,'N_LINES_VIS=',4(i5,','))
44    format(1x,'N_PIXELS_IR=',4(i5,','))
45    format(1x,'N_LINES_IR=',4(i5,','))
46    format(1x,'N_PIXELS_WV=',4(i5,','))
47    format(1x,'N_LINES_WV=',4(i5,','))
48    format(1x,'I_EWCYCLES=',4(i5,','))
49    format(1x,'I_EWINCS=',4(i5,','))
50    format(1x,'I_NSCYCLES=',4(i5,','))
51    format(1x,'I_NSINCS=',4(i5,','))
52    format(1x,'I_NWLINE_VIS=',4(i5,','))
53    format(1x,'I_NWLINE_IR=',4(i5,','))
54    format(1x,'I_NWLINE_WV=',4(i5,','))
55    format(1x,'I_NWPIX_VIS=',4(i5,','))
56    format(1x,'I_NWPIX_IR=',4(i5,','))
57    format(1x,'I_NWPIX_WV=',4(i5,','))
58    format(1x,'IMC=',4(i4,","))
59    format(1x,"/")
      close(lun)
      istatus=1
      return

 900  print*,'error opening file ',nest7grid
      return
 950  print*,'error writing satellite_lvd.nl in ',nest7grid
      return
      end
