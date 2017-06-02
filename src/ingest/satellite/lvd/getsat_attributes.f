      subroutine getsat_attributes(csat_id,csat_type,chtype,
     &istart,iend,jstart,jend,nx,ny,istatus)

      implicit none

      integer   i,j,k
      integer   ispec
      integer   istart,iend
      integer   jstart,jend
      integer   nx,ny
      integer   istatus

      character csat_id*6
      character csat_type*3
      character chtype*3

      include 'satellite_dims_lvd.inc'
      include 'satellite_common_lvd.inc'

      istatus = 0

      do k=1,maxsat
       if(c_sat_id(k).eq.csat_id)then

        do j=1,maxtype
         if(c_sat_types(j,k).eq.csat_type)then

          do i=1,maxchannel
           if(c_channel_types(i,j,k).eq.chtype)then

            call lvd_file_specifier(chtype,ispec,istatus)
            if(ispec.eq.1)then
              istart=i_start_vis(j,k)
              iend  =i_end_vis(j,k)
              jstart=j_start_vis(j,k)
              jend  =j_end_vis(j,k)
              nx    =n_pixels_vis(j,k)
              ny    =n_lines_vis(j,k)
            elseif(ispec.eq.2.or.ispec.eq.4.or.ispec.eq.5)then
              istart=i_start_ir(j,k)
              iend  =i_end_ir(j,k)
              jstart=j_start_ir(j,k)
              jend  =j_end_ir(j,k)
              nx    =n_pixels_ir(j,k)
              ny    =n_lines_ir(j,k)
            elseif(ispec.eq.3)then
              istart=i_start_wv(j,k)
              iend  =i_end_wv(j,k)
              jstart=j_start_wv(j,k)
              jend  =j_end_wv(j,k)
              nx    =n_pixels_wv(j,k)
              ny    =n_lines_wv(j,k)
            endif
           endif
          enddo
         endif
        enddo
       endif
      enddo 

      if(min(istart,iend,jstart,jend).le.0)then
         write(6,*)'Error in getsat_attributes'
         write(6,*)'istart, iend, jstart, or jend = 0'
         return
      else
         write(6,1)istart,iend,jstart,jend
1        format(' getsat_attributes i/j start/end',4i6)       
      endif

      istatus = 1

      return
      end
