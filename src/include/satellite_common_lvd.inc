      character*200 path_to_raw_sat(maxchannel,maxtype,maxsat)

      integer iflag_lvd_common
      integer isats(maxsat)
      integer itypes(maxtype,maxsat)
      integer ichannels(maxchannel,maxtype,maxsat)
      integer i_delta_sat_t_sec
      integer i_msng_sat_flag(maxtype,maxsat)
      integer n_images

      integer image_width_ir(maxtype,maxsat)
      integer image_depth_ir(maxtype,maxsat)
      integer image_width_vis(maxtype,maxsat)
      integer image_depth_vis(maxtype,maxsat)
      integer image_width_wv(maxtype,maxsat)
      integer image_depth_wv(maxtype,maxsat)

      character c_sat_id(maxsat)*6
      character c_sat_types(maxtype,maxsat)*3
      character c_channel_types(maxchannel,maxtype,maxsat)*3

      real*4 r_sat_sub_lat(maxsat)
      real*4 r_sat_sub_lon(maxsat)
      real*4 sat_range_m(maxsat)

      real*4 r_la1(maxtype,maxsat)
      real*4 r_lo1(maxtype,maxsat)
      real*4 r_lap(maxtype,maxsat)
      real*4 r_lov(maxtype,maxsat)
      real*4 r_latin(maxtype,maxsat)

      integer i_start_vis(maxtype,maxsat)
      integer i_end_vis(maxtype,maxsat)
      integer j_start_vis(maxtype,maxsat)
      integer j_end_vis(maxtype,maxsat)
      integer i_start_ir(maxtype,maxsat)
      integer i_end_ir(maxtype,maxsat)
      integer j_start_ir(maxtype,maxsat)
      integer j_end_ir(maxtype,maxsat)
      integer i_start_wv(maxtype,maxsat)
      integer i_end_wv(maxtype,maxsat)
      integer j_start_wv(maxtype,maxsat)
      integer j_end_wv(maxtype,maxsat)

      integer n_pixels_vis(maxtype,maxsat)
      integer n_lines_vis(maxtype,maxsat)
      integer n_pixels_ir(maxtype,maxsat)
      integer n_lines_ir(maxtype,maxsat)
      integer n_pixels_wv(maxtype,maxsat)
      integer n_lines_wv(maxtype,maxsat)
      integer i_ewCycles(maxtype,maxsat)
      integer i_ewIncs(maxtype,maxsat)
      integer i_nsCycles(maxtype,maxsat)
      integer i_nsIncs(maxtype,maxsat)
      integer i_nwline_vis(maxtype,maxsat)
      integer i_nwline_ir(maxtype,maxsat)
      integer i_nwline_wv(maxtype,maxsat)
      integer i_nwpix_vis(maxtype,maxsat)
      integer i_nwpix_ir(maxtype,maxsat)
      integer i_nwpix_wv(maxtype,maxsat)
      integer imc

      real*4  r_resolution_x_vis(maxtype,maxsat)
      real*4  r_resolution_y_vis(maxtype,maxsat)
      real*4  r_resolution_x_ir(maxtype,maxsat)
      real*4  r_resolution_y_ir(maxtype,maxsat)
      real*4  r_resolution_x_wv(maxtype,maxsat)
      real*4  r_resolution_y_wv(maxtype,maxsat)


      common /lvd_namelist_cmn/path_to_raw_sat,
     +iflag_lvd_common,isats,itypes,ichannels,
     +i_delta_sat_t_sec,i_msng_sat_flag,
     +sat_range_m,n_images,
     +image_width_ir,image_depth_ir,
     +image_width_vis,image_depth_vis,
     +image_width_wv,image_depth_wv,
     +c_sat_id, c_sat_types,c_channel_types,
     +r_sat_sub_lat, r_sat_sub_lon,
     +r_la1, r_lo1,
     +r_resolution_x_vis,r_resolution_y_vis,r_resolution_x_ir,
     +r_resolution_y_ir,r_resolution_x_wv,r_resolution_y_wv,
     +r_latin,r_Lov,r_lap,
     +i_start_vis,i_end_vis,
     +j_start_vis,j_end_vis,
     +i_start_ir, i_end_ir,
     +j_start_ir, j_end_ir,
     +i_start_wv, i_end_wv,
     +j_start_wv, j_end_wv,
     +n_pixels_vis,n_lines_vis,
     +n_pixels_ir,n_lines_ir,
     +n_pixels_wv,n_lines_wv,
     +i_ewCycles,i_ewIncs,
     +i_nsCycles,i_nsIncs,
     +i_nwline_vis,i_nwline_ir,i_nwline_wv,
     +i_nwpix_vis,i_nwpix_ir,i_nwpix_wv,
     +imc
