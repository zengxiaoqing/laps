      subroutine read_gms_taiwan(path_to_raw_sat,nxt,nyt
     &,maxchannels,max_files,nchannels,csatid,csattype
     &,chtype,i4time_needed,nelems_ir,nlines_ir,nelems_vis
     &,nlines_vis,nelems_wv,nlines_wv,image_11u,image_vis
     &,image_wv,image_12u,nimages,nft,ntm,c_type,i4time_data
     &,istatus)
c
c routine to read binary satellite data file from
c Taiwan CWB.
c B.S. Wang (CWB)  and J.R. Smart (FSL)   June 2001
c
      implicit none

      integer     nxt,nyt
      integer     maxchannels
      integer     max_files
      integer     nchannels
      integer     nimages
      integer     nelems_ir,nlines_ir
      integer     nelems_vis,nlines_vis
      integer     nelems_wv,nlines_wv
      integer     i4time_needed
      integer     i4time_latest_sat
      integer     i4time_latest_lvd
      integer     ispec

      character*1   c_imagedata(nxt,nyt)
      character*3   c_type(maxchannels,max_files)
      character*3   csattype
      character*3   chtype(maxchannels)
      character*6   csatid
      character*9   fname9_sat
      character*200 fname_sat
      character*200 fname_lvd
      character*200 path_to_raw_sat(maxchannels)
      character*255 filename

      integer     istatus,lf,lp,ls,lt,lv
      integer     i,j,ii,jj,icnt,k
      integer     istart,jstart
      integer     iend,jend
      integer     nft,ntm(maxchannels)
      integer     i4time_data(max_files)
      integer     nx,ny
      real        image_11u(nelems_ir,nlines_ir,nimages)
      real        image_12u(nelems_ir,nlines_ir,nimages)
      real        image_vis(nelems_vis,nlines_vis,nimages)
      real        image_wv(nelems_wv,nlines_wv,nimages)

      istatus = 0

      nft=1

      do k = 1,nchannels

         call getsat_attributes(csatid, csattype, chtype(k)
     &,istart,iend,jstart,jend,nx,ny,istatus)

         call lvd_file_specifier(chtype(k),ispec,istatus)
         if(istatus.ne.0)then
            print*,'Error setting ispec from lvd_file_specifier'
            return
         endif

c add file-time source here. CWB names are yyjjjhhmm_"type"

         call s_len(path_to_raw_sat(ispec),lp)
         fname_sat=path_to_raw_sat(ispec)(1:lp)//'/*_'
         call s_len(fname_sat,ls)
         fname_sat=fname_sat(1:ls)//chtype(k)
         call s_len(fname_sat,ls)

         print*
         print*,'Get latest satellite time'
         call get_latest_file_time(fname_sat,i4time_latest_sat) 
         call make_fnam_lp(i4time_latest_sat,fname9_sat,istatus)

         print*,'Get latest lvd time'
         call get_directory('lvd',fname_lvd,lv)
         fname_lvd=fname_lvd(1:lv)//csatid//'/*.lvd'
         call get_latest_file_time(fname_lvd,i4time_latest_lvd)

c check if this is new data

         if(i4time_latest_sat.le.i4time_latest_lvd)then
            print*,'No new data found for twn/gmssat/',chtype(k)
            goto 20
         endif

         filename=path_to_raw_sat(ispec)(1:lp)//'/'//fname9_sat
         call s_len(filename,lf)
         call s_len(chtype(k),lt)
         filename=filename(1:lf)//'_'//chtype(k)(1:lt)
         call s_len(filename,lf)


         call s_len(filename,lf)

         open(10,file=filename,access='direct',recl=nxt,err =99)

         call read_binary_field(c_imagedata,1,1,nxt*nyt,filename,lf)

c load block for known i-j start/end locations

         if(ispec.eq.1)then

            call load_image_twn(nelems_vis,nlines_vis
     &,istart,iend,jstart,jend,nxt,nyt,c_imagedata,image_vis(1,1,1))

            ntm(nft)=ntm(nft)+1
            c_type(ntm(nft),nft)=chtype(k)
            i4time_data(nft)=i4time_latest_sat

         elseif(ispec.eq.3)then

            call load_image_twn(nelems_wv,nlines_wv
     &,istart,iend,jstart,jend,nxt,nyt,c_imagedata,image_wv(1,1,1))

            ntm(nft)=ntm(nft)+1
            c_type(ntm(nft),nft)=chtype(k)
            i4time_data(nft)=i4time_latest_sat

         elseif(ispec.eq.4)then

            call load_image_twn(nelems_ir,nlines_ir
     &,istart,iend,jstart,jend,nxt,nyt,c_imagedata
     &,image_11u(1,1,1))

            ntm(nft)=ntm(nft)+1
            c_type(ntm(nft),nft)=chtype(k)
            i4time_data(nft)=i4time_latest_sat

         elseif(ispec.eq.5)then

            call load_image_twn(nelems_ir,nlines_ir
     &,istart,iend,jstart,jend,nxt,nyt,c_imagedata
     &,image_12u(1,1,1))

            ntm(nft)=ntm(nft)+1
            c_type(ntm(nft),nft)=chtype(k)
            i4time_data(nft)=i4time_latest_sat

         else

            print*,'read_gms_taiwan: unknown type = ',chtype(k)
            print*,'returning no image data'
            return 

         endif

         close (10)

20    enddo

      istatus = 1
      return

99    print*,'Error opening twn file: ',filename(1:lf)
      nft=0
      return
      end 

      subroutine load_image_twn(nelems,nlines
     &,istart,iend,jstart,jend,nxt,nyt,c_imagedata,r_image)

      implicit none

      integer nxt,nyt
      integer nelems,nlines
      integer istart,iend
      integer jstart,jend
      integer i,j,ii,jj
      integer icnt
      integer istatus

      character*(*) c_imagedata(nxt,nyt)

      real    r_image(nelems,nlines)
      real    r_imagedata(nelems,nlines)
      real    r_imagetmp(nxt,nyt)
      real    rmsngdata
      real    r_missing_data

      rmsngdata = 0.0
      call get_r_missing_data(r_missing_data,istatus)
      if(istatus.ne.1)then
         print*,'Error getting r_missing_data'
         return
      endif

      if(jend-jstart+1.ne.nlines)then
         print*,'jend - jstart + 1 /!=/ nlines: ',jend-jstart+1,nlines
      endif

      if(iend-istart+1.ne.nelems)then
         print*,'iend - istart + 1 /!=/ nelems: ',iend-istart+1,nelems
      endif

c reverse order (CWB data is 1,1 in NW corner)
      do j = 1,nyt
      do i = 1,nxt
         r_imagetmp(i,j)=float(ichar(c_imagedata(i,nyt-j+1)))
      enddo
      enddo

c convert to counts and load block.
 
      icnt=0
      jj=0

      do j=jstart,jend

         jj=jj+1
         ii=0

         do i=istart,iend

            ii=ii+1
            r_image(ii,jj)=r_imagetmp(i,j) 

            if(r_image(ii,jj).le. 0.0    .or.
     &         r_image(ii,jj).gt. 255.0  .or.
     &         r_image(ii,jj).eq. rmsngdata)then

               r_image(ii,jj)=r_missing_data
               icnt=icnt+1

            endif

         enddo
      enddo

      if(icnt.gt.0)then
         print*,'found missing data: twn ',icnt
      endif

      return
      end
