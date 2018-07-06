       subroutine set_vis_flag(i4time,lat,lon,nx_l,ny_l,lvis_flag)
c
c routine sets flag for visible data availability. When the solar
c elevation is less than 14 degrees then the visible data is useless
c even if the file is available. This routine returns lvis_flag = true
c when the solar elevation angle is less than the threshold.
c
       use mem_namelist, ONLY: solalt_thr_vis

       implicit None
c
       integer nx_l,ny_l
       integer i4time
       logical   lvis_flag
       real    lat(nx_l,ny_l)
       real    lon(nx_l,ny_l)
       real    alt(nx_l,ny_l)
       real    azi(nx_l,ny_l)

       integer i,j
       real    solar_alt_max
c
       call get_solaltaz_2d(lat,lon,i4time,nx_l,ny_l,alt,azi)
       solar_alt_max = maxval(alt)

       write(6,*)' max / thresh solar altitude is ',solar_alt_max
     1                                             ,solalt_thr_vis

       if(solar_alt_max.lt.solalt_thr_vis)lvis_flag=.true.

       return
       end
