      subroutine satdatfill(csat_id,csat_type,nft,ntm,
     &    nir_elem,nir_lines,nvis_elem,nvis_lines,
     &    nwv_elem,nwv_lines,c_type,smsng,maxchannels,
     &    maximage,maxfiles,i4time_data,sat_dir_path,
     &    image_ir,image_39,image_12,image_67,image_vis,
     &    scale_img,mstatus)
c
      implicit none

      integer    maxchannels
      integer    maximage
      integer    maxfiles
      integer    nir_elem,nir_lines
      integer    nvis_elem,nvis_lines
      integer    nwv_elem,nwv_lines

      character*3 c_type(maxchannels,maxfiles)
      character*3 csat_type
      character*6 csat_id
      character*200 sat_dir_path(maxchannels)
      character*200 c_satdir

      real image_vis (nvis_elem,nvis_lines,maximage)
      real image_ir  (nir_elem,nir_lines,maximage)
      real image_12  (nir_elem,nir_lines,maximage)
      real image_39  (nir_elem,nir_lines,maximage)
      real image_67  (nwv_elem,nwv_lines,maximage)

      Real        r_missing_data, scale_img
      Real        rmsng
      Real        smsng(maxchannels)
      Real        mstatus(maxchannels,maxfiles)
      Real        percent_missing

      integer i4time_data(maxfiles)
      integer istatus
      integer jstatus
      integer istat
      integer ispec
      integer i,j
      integer nft,ntm(maxfiles)
c
c-----------------------------------------------------------------------
c------------------------------------------------------------------------
c ===================================================== 72 chars to here>
c
      call get_r_missing_data(r_missing_data, jstatus)
      if(jstatus.ne.1)then
         write(6,*)'Error getting r_missing_data - satdatfill'
         return
      endif
c
c fill missing pixels and pixels determined as bad with r_missing_data
c
      do i=1,nft
         do j=1,ntm(i)

            write(6,*)
            write(6,*)'Satellite quality control - ',c_type(j,i)
            write(6,*)'------------------------------'
            write(6,*)'first set missing sat: '

            call lvd_file_specifier(c_type(j,i),ispec,istat)
            if(ispec.eq.2)then
c -------------------------------------------
            rmsng=smsng(ispec)
            call set_missing_sat(csat_id,csat_type,c_type(j,i),
     &               image_39(1,1,i),nir_elem,nir_lines,
     &               rmsng,r_missing_data,scale_img,istatus)

            write(6,*)'  Missing status : ',abs(istatus)
            write(6,*)'  Enter Satfill1 for: ',c_type(j,i)

            c_satdir=sat_dir_path(ispec)
            call satfill1(csat_id,csat_type,
     &               i4time_data(i),rmsng,
     &               c_type(j,i),
     &               nir_lines,nir_elem,
     &               c_satdir,
     &               r_missing_data,
     &               image_39(1,1,i),
     &               percent_missing,
     &               istatus)
            if(istatus .le. -1)then
               write(6,*)'  Previous IR image bad' 
            endif

            mstatus(j,i)=percent_missing

            elseif(ispec.eq.4)then
c -------------------------------------------
            rmsng=smsng(ispec)
            call set_missing_sat(csat_id,csat_type,c_type(j,i),
     &               image_ir(1,1,i),nir_elem,nir_lines,
     &               rmsng,r_missing_data,scale_img,
     &               istatus)

            write(6,*)'  Missing status : ',abs(istatus)
            write(6,*)'  Enter Satfill1 for: ',c_type(j,i)

            c_satdir=sat_dir_path(ispec)
            call satfill1(csat_id,csat_type,
     &               i4time_data(i),rmsng,
     &               c_type(j,i),
     &               nir_lines,nir_elem,
     &               c_satdir,
     &               r_missing_data,
     &               image_ir(1,1,i),
     &               percent_missing,
     &               istatus)
            if(istatus .le. -1)then
               write(6,*)'  Previous IR image bad'
            endif

            mstatus(j,i)=percent_missing

            elseif(ispec.eq.5)then
c -------------------------------------------
            rmsng=smsng(ispec)
            call set_missing_sat(csat_id,csat_type,c_type(j,i),
     &               image_12(1,1,i),nir_elem,nir_lines,
     &               rmsng,r_missing_data,scale_img,
     &               istatus)

            write(6,*)'  Missing status : ',abs(istatus)
            write(6,*)'  Enter Satfill1 for: ',c_type(j,i)

            c_satdir=sat_dir_path(ispec)
            call satfill1(csat_id,csat_type,
     &               i4time_data(i),rmsng,
     &               c_type(j,i),
     &               nir_lines,nir_elem,
     &               c_satdir,
     &               r_missing_data,
     &               image_12(1,1,i),
     &               percent_missing,
     &               istatus)
            if(istatus .le. -1)then
               write(6,*)'  Previous IR image bad'
            endif

            mstatus(j,i)=percent_missing

            elseif(ispec.eq.3)then
c -------------------------------------------
            rmsng=smsng(ispec)
            call set_missing_sat(csat_id,csat_type,c_type(j,i),
     &               image_67(1,1,i),nwv_elem,nwv_lines,
     &               rmsng,r_missing_data,scale_img,
     &               istatus)

            write(6,*)'  Missing status : ',abs(istatus)
            write(6,*)'  Enter Satfill1 for: ',c_type(j,i)

            c_satdir=sat_dir_path(ispec)
            call satfill1(csat_id,csat_type,
     &               i4time_data(i),rmsng,
     &               c_type(j,i),
     &               nwv_lines, nwv_elem,
     &               c_satdir,
     &               r_missing_data,
     &               image_67(1,1,i),
     &               percent_missing,
     &               istatus)
            if(istatus .le. -1)then
               write(6,*)'  Previous WV image bad'
            endif

            mstatus(j,i)=percent_missing

            elseif(ispec.eq.1)then
c -------------------------------------------
            rmsng=smsng(ispec)
            c_satdir=sat_dir_path(ispec)  !Particularly wfo! The channel types are in order.
            call set_missing_sat(csat_id,csat_type,c_type(j,i),
     &               image_vis(1,1,i),nvis_elem,nvis_lines,
     &               rmsng,r_missing_data,scale_img,
     &               istatus)

            write(6,*)'  Missing status : ',abs(istatus)

            write(6,*)'  Enter Satellite fill '
            call satfill1(csat_id,csat_type,
     &               i4time_data(i),rmsng,
     &               c_type(j,i),
     &               nvis_lines, nvis_elem,
     &               c_satdir,
     &               r_missing_data,
     &               image_vis(1,1,i),
     &               percent_missing,
     &               istatus)
            if(istatus .le. -1)then
               write(6,*)'  Previous VIS image bad'
            endif

            mstatus(j,i)=percent_missing

            endif

         enddo
      enddo

1000  return
      end
