      subroutine set_missing_flag(csatid,csat_type,nir_elem,nir_lines,
     &                   nvis_elem,nvis_lines,nwv_elem,nwv_lines,
     &                   nft,ntm,c_type,smsng,maxchannels,nimages,
     &                   image_ir,image_39,image_12,image_67,
     &                   image_vis,
     &                   mstatus)
c
c
c
      implicit none

      integer    maxchannels
      integer    nimages

      character*3   c_type(maxchannels,nimages)
      character*3   csat_type
      character*(*) csatid

      integer    nwv_elem,nwv_lines
      integer    nir_elem,nir_lines
      integer    nvis_elem,nvis_lines

      real image_vis (nvis_elem,nvis_lines,nimages) !ispan vis image polar NH
      real image_ir  (nir_elem,nir_lines,nimages) !ispan ir image polar NH
      real image_12  (nir_elem,nir_lines,nimages)
      real image_39  (nir_elem,nir_lines,nimages)
      real image_67  (nwv_elem,nwv_lines,nimages)
      real image_data_ir(nir_elem,nir_lines)
      real image_data_vis(nvis_elem,nvis_lines)
      real image_data_wv(nwv_elem,nwv_lines)

      Real*4      smsng
      Real*4      r_missing_data
      Real*4      tot_ir_pix,tot_vis_pix,tot_wv_pix
      Real*4      mstatus(maxchannels,nimages)

      integer i,j
      integer ispec
      integer nft,ntm(nimages)
      integer istatus
      integer istat
c
c =====================================================
c
c get r_missing_data from nest7grid.parms
c
      call get_r_missing_data(r_missing_data, istatus)
      if(istatus.ne.1)goto 900

      tot_ir_pix   = float(nir_lines*nir_elem)
      tot_vis_pix  = float(nvis_lines*nvis_elem)
      tot_wv_pix   = float(nwv_lines*nwv_elem)

      do i=1,nft
         do j=1,ntm(i)
            call lvd_file_specifier(c_type(j,i),ispec,istat)

            if(ispec.eq.4)then
               call move(image_ir(1,1,i),image_data_ir,
     &                      nir_elem,nir_lines)
            elseif(ispec.eq.5)then
               call move(image_12(1,1,i),image_data_ir,
     &                      nir_elem,nir_lines)
            elseif(ispec.eq.2)then
               call move(image_39(1,1,i),image_data_ir,
     &                      nir_elem,nir_lines)
            elseif(ispec.eq.3)then
               call move(image_67(1,1,i),image_data_wv,
     &                      nwv_elem,nwv_lines)
            elseif(ispec.eq.1)then
               call move(image_vis(1,1,i),image_data_vis,
     &                      nvis_elem,nvis_lines)
            endif

            write(6,*)'Enter set missing sat: ',c_type(j,i)

            if(ispec.eq.2.or.ispec.eq.4.or.ispec.eq.5)then

               call set_missing_sat(csatid,csat_type,c_type(j,i),
     &               image_data_ir,
     &               nir_elem,nir_lines,
     &               smsng,r_missing_data,
     &               istatus)

               mstatus(j,i)=float(abs(istatus))/tot_ir_pix
               write(6,*)'Missing status (%): ',mstatus(j,i)

            elseif(ispec.eq.3)then

               call set_missing_sat(csatid,csat_type,c_type(j,i),
     &               image_data_wv,
     &               nwv_elem,nwv_lines,
     &               smsng,r_missing_data,
     &               istatus)

               mstatus(j,i)=float(abs(istatus))/tot_wv_pix
               write(6,*)'Missing status (%): ',mstatus(j,i)

            else    !must be the vis data

               call set_missing_sat(csatid,csat_type,c_type(j,i),
     &               image_data_vis,
     &               nvis_elem,nvis_lines,
     &               smsng,r_missing_data,
     &               istatus)

               mstatus(j,i)=float(abs(istatus))/tot_vis_pix
               write(6,*)'Missing status (%): ',mstatus(j,i)

            endif

            if(ispec.eq.4)then
               call move(image_data_ir,image_ir(1,1,i),
     &                      nir_elem,nir_lines)
            elseif(ispec.eq.5)then
               call move(image_data_ir,image_12(1,1,i),
     &                      nir_elem,nir_lines)
            elseif(ispec.eq.2)then
               call move(image_data_ir,image_39(1,1,i),
     &                      nir_elem,nir_lines)
            elseif(ispec.eq.3)then
               call move(image_data_wv,image_67(1,1,i),
     &                      nwv_elem,nwv_lines)
            elseif(ispec.eq.1)then
               call move(image_data_vis,image_vis(1,1,i),
     &                      nvis_elem,nvis_lines)
            endif

         enddo
      enddo

      goto 1000

900   write(6,*)'Did not get r_missing_data - set_missing_flag'

1000  return
      end
