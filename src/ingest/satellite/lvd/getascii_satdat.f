       subroutine getascii_satdat(i4time_current,
     &                        lvis_flag,i_delta_t,
     &                        nirlines, nirelem,
     &                        nvislines,nviselem,
     &                        nwvlines,nwvelem,
     &                        sat_dir_path,
     &                        ntm,c_type,max_channels,
     &                        image_ir,image_vis,
     &                        image_12,image_39,image_67,
     &                        i4time_data,
     &                        grid_spacing_km,
     &                        istatus)
c
c
c
       implicit none

       include 'satsector.inc'

       integer i,j,n
       integer ntm
       integer nirelem
       integer nirlines
       integer nviselem
       integer nvislines
       integer nwvelem
       integer nwvlines
       integer max_channels

       integer i4time_current
       integer i4time_data
       integer i4time_data_ir
       integer i4time_data_vis
       integer i4time_diff
       integer i_delta_t
       integer nelem_ir
       integer nlines_ir
       integer nelem_vis
       integer nlines_vis

       integer istatus
       integer istatus_ir
       integer istatus_vis

       real image_ir  (nirelem,nirlines) 
       real image_12  (nirelem,nirlines)
       real image_39  (nirelem,nirlines)
       real image_67  (nwvelem,nwvlines)
       real image_vis (nviselem,nvislines)

       real rlat_ir(n_ir_elem,n_ir_lines)
       real rlon_ir(n_ir_elem,n_ir_lines)
       real rlat_vis(n_vis_elem,n_vis_lines)
       real rlon_vis(n_vis_elem,n_vis_lines)

       integer imageline_ir(n_ir_elem,n_ir_lines)
       integer imageelem_ir(n_ir_elem,n_ir_lines)
       integer imageline_vis(n_vis_elem,n_vis_lines)
       integer imageelem_vis(n_vis_elem,n_vis_lines)

       integer image_data(n_ir_elem,n_ir_lines)

       real    grid_spacing
       real    grid_spacing_km(max_channels)
       real    r4time_data_ir
       real    r4time_data_vis

       logical   lvis_flag

       character c_type(max_channels)*3
       character c_filetime*9
       character sat_dir_path*200
       character c_filename*255
c ==============================================
c first try for the ir
c
       istatus = 1
       n=index(sat_dir_path,' ')
       c_filename=sat_dir_path(1:n-1)//'LAPSI.ASC'
       n=index(c_filename,' ')
       write(6,*)'Reading: ',c_filename(1:n)
       call read_ascii_satdat(c_filename,
     &                        i4time_current,
     &                        c_filetime,i_delta_t,
     &                        n_ir_lines,
     &                        n_ir_elem,
     &                        nelem_ir,nlines_ir,
     &                        rlat_ir,rlon_ir,
     &                        imageline_ir,imageelem_ir,
     &                        image_data,
     &                        i4time_data_ir,
     &                        grid_spacing,
     &                        istatus_ir)
       if(istatus_ir .eq. 1)then
          call loadarray_i2r(image_data,n_ir_elem,n_ir_lines,
     &                    image_ir)
          ntm=ntm+1
          c_type(ntm)='ir '
          grid_spacing_km(ntm)=grid_spacing
       else
          write(6,*)'Not successful getting IR image'
       endif
c
c now try for the vis
c
       if(.not.lvis_flag)then

          n=index(sat_dir_path,' ')
          c_filename=sat_dir_path(1:n-1)//'LAPSV.ASC'
          n=index(c_filename,' ')
          write(6,*)'Reading: ',c_filename(1:n)
          call read_ascii_satdat(c_filename,
     &                        i4time_current,
     &                        c_filetime,i_delta_t,
     &                        n_vis_lines,
     &                        n_vis_elem,
     &                        nelem_vis,nlines_vis,
     &                        rlat_vis,rlon_vis,
     &                        imageline_vis,imageelem_vis,
     &                        image_data,
     &                        i4time_data_vis,
     &                        grid_spacing,
     &                        istatus_vis)
          if(istatus_vis .eq. 1)then
             call loadarray_i2r(image_data,
     &                          n_vis_elem,n_vis_lines,
     &                          image_vis)
             ntm=ntm+1
             c_type(ntm)='vis'
             grid_spacing_km(ntm)=grid_spacing

          else
             write(6,*)'Not successful getting VIS image'
          endif

       else

          write(6,*)'lvis_flag set! Not reading vis data'
          istatus_vis=-1

       endif
c
c this section needed to prevent either vis or ir which is actually old
c but satisfies the i_delta_t criterion, to be neglected.
c
       if(ntm.gt.1)then
          i4time_diff=i4time_data_ir-i4time_data_vis
          if(i4time_diff.gt.0)then
             ntm=ntm-1
             c_type(ntm)='ir '
             write(6,*)'VIS is actually old'
             istatus_vis = -1
          elseif(i4time_diff.lt.0)then
             ntm=ntm-1
             c_type(ntm)='vis'
             write(6,*)'IR is actually old'
             istatus_ir = -1
          elseif(i4time_diff.eq.0)then
             write(6,*)'IR and VIS time are equal'
             write(6,*)'Using both IR and VIS'
          endif
       endif
       write(6,*)

       if(ntm.gt.0)then
          if(lvis_flag)then
             i4time_data=i4time_data_ir
          else

             i4time_diff=i4time_data_ir-i4time_data_vis

             if(i4time_diff.eq.0)then                   !ir and vis have same (current) time.
                i4time_data=i4time_data_ir
             elseif(ntm.eq.1.and.i4time_diff.gt.0)then  !ir is current data
                i4time_data=i4time_data_ir
                write(6,*)'Using i4time_data_ir'
             elseif(ntm.eq.1.and.i4time_diff.lt.0)then  !vis is current data
                i4time_data=i4time_data_vis
                write(6,*)'Use i4time_data_vis for filename'
             elseif(ntm.eq.2)then                       !ir and vis current but different times
c               r4time_data_ir = i4time_data_ir/10000.
c               r4time_data_vis = i4time_data_vis/10000.
                if(i4time_diff.gt.0)then
                   i4time_data=i4time_data_ir
                else
                   i4time_data=i4time_data_vis
                endif
                write(6,*)'Using most current time between ir and vis'
             else
                write(6,*)'Something is wrong, ntm > 2!'
                istatus=-1
             endif
          endif
       endif

       if(istatus_vis.eq.-1.and.istatus_ir.eq.-1)then
          istatus=-1
       endif

       return
       end
