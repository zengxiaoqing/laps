       subroutine set_vis_flag(i4time,lat,lon,nx_l,ny_l,lvis_flag)
c
c routine sets flag for visible data availability. When the solar
c elevation is less than 14 degrees then the visible data is useless
c even if the file is available. This routine returns lvis_flag = true
c when the solar elevation angle is less than the threshold.
c
       implicit None
c
       integer nx_l,ny_l
       integer i4time
       logical   lvis_flag
       real    lat(nx_l,ny_l)
       real    lon(nx_l,ny_l)

       integer i,j
       real    solar_alt_thresh
       real    solar_alt_max
       real    alt,dec,hrangle
c
c
       solar_alt_thresh=14.0
       solar_alt_max=0.0
       do j=1,ny_l
       do i=1,nx_l

          call solar_position(lat(i,j),lon(i,j),i4time,alt,dec,hrangle)
          if(alt.gt.solar_alt_max)solar_alt_max=alt

       enddo
       enddo

       if(solar_alt_max.lt.solar_alt_thresh)lvis_flag=.true.

       return
       end
