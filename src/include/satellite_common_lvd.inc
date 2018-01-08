      character path_to_raw_sat(maxchannel,maxtype,maxsat)*200

      integer isats(maxsat)
      integer itypes(maxtype,maxsat)
      integer ichannels(maxchannel,maxtype,maxsat)
      integer i_delta_sat_t_sec
      integer i_msng_sat_flag(maxchannel,maxtype,maxsat)
      integer iskip_bilin
      integer nsats,ntypes(maxsat),nchannel(maxtype*maxsat)
      integer n_images

!these are in the common block and used throughout the lvd code
      character c_sat_id(maxsat)*6
      character c_sat_types(maxtype,maxsat)*3
      character c_channel_types(maxchannel,maxtype,maxsat)*3

      real   r_sat_sub_lat(maxsat)
      real   r_sat_sub_lon(maxsat)
      real   sat_range_m(maxsat)

      real   r_la1(maxtype,maxsat)   ! lower left latitude
      real   r_lo1(maxtype,maxsat)   ! lower left longitude
      real   r_lap(maxtype,maxsat)   ! standard latitude
      real   r_lov(maxtype,maxsat)   ! standard longitude
      real   r_latin(maxtype,maxsat) 

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
      integer imc(maxsat)

      logical l_national
      logical l_cell_afwa

      real    r_resolution_x_vis(maxtype,maxsat)
      real    r_resolution_y_vis(maxtype,maxsat)
      real    r_resolution_x_ir(maxtype,maxsat)
      real    r_resolution_y_ir(maxtype,maxsat)
      real    r_resolution_x_wv(maxtype,maxsat)
      real    r_resolution_y_wv(maxtype,maxsat)
      real    vis_cnt_range_in(2,maxsat)
      real    vis_cnt_range_out(2,maxsat)
      real    solalt_thr_vis


      common /lvd_namelist_cmn/
     +l_cell_afwa,l_national,iskip_bilin,
     +isats,itypes,ichannels,
     +nsats,ntypes,nchannel,
     +i_delta_sat_t_sec,i_msng_sat_flag,
     +sat_range_m,n_images,solalt_thr_vis,
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
     +imc,path_to_raw_sat,
     +c_sat_id,c_sat_types,c_channel_types,
     +vis_cnt_range_in,vis_cnt_range_out
