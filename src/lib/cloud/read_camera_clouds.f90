
       subroutine get_camera_clouds(minalt,maxalt,minazi,maxazi,alt_scale,azi_scale, & ! I
                                    i4time,fname_ppm,mode, &                           ! I
                                    mask_cyl,istatus)                                  ! O
 
       use ppm

       parameter (nip = 511)
       parameter (njp = 511)

       integer mask_cyl(minalt:maxalt,minazi:maxazi)  ! O
                                                      ! 0 means unknown
                                                      ! 1 means clear
                                                      ! 2 means cloud
       integer image(minalt:maxalt,minazi:maxazi) ! L
       integer img_polar(3,511,511)
       integer mask_polar(511,511)
       integer, allocatable :: img(:,:,:)
       integer u /12/ 

       character*255 imgfile,img_png,img_ppm,convert_cmd,camera_path
       character*13 a13name, cvt_i4time_wfo_fname13
       character*9 a9time
       character*10 fname_ppm

!      Statement functions
       i2ialt(i) = minalt+i-1
       altf(ialt) = (ialt-minalt) * alt_scale

       j2jazi(i) = minazi+j-1
       azif(jazi) = (jazi-minazi) * azi_scale

!      Executable statements
       camera_path = '.' ! not needed for mode=2
       mask_cyl = 0

       a13name = cvt_i4time_wfo_fname13(i4time)
       call make_fnam_lp(i4time,a9time,istatus)

       if(.false.)then      ! Read polar observed image in PPM format 

         imgfile = '/w3/lapb/allsky/sites/dsrc/output3'//'/'//a13name//'_dsrc.ppm'

         open(u,file=trim(imgfile),status='old',err=999)
         read(u,*)   
         read(u,*)iwidth,iheight
         rewind(u)
         write(6,*)' dynamic dims ',ncol,iwidth,iheight
         allocate(img(ncol,iwidth,iheight))
         call read_ppm(u,img,ncol,iwidth,iheight)
         close(u)
         write(6,*)' polar image has been read in'

       elseif(.false.)then ! Read cyl observed image in PPM format 
                           ! (remapped to alt/az grid)

         imgfile = '/w3/lapb/allsky/sites/dsrc/output3c'//'/'//a13name//'_dsrc.ppm'

         open(u,file=trim(imgfile),status='old',err=999)
         read(u,*)   
         read(u,*)iwidth,iheight
         rewind(u)
         write(6,*)' dynamic dims ',ncol,iwidth,iheight
         allocate(img(ncol,iwidth,iheight))
         call read_ppm(u,img,ncol,iwidth,iheight)
         close(u)
         write(6,*)' cyl image has been read in'

       elseif(mode .eq. 1)then ! Read polar (contingency table) mask
         img_png = trim(camera_path)//'/verif_allsky_mask.dsrc.'//a9time//'.png'
         img_ppm = trim(camera_path)//'/verif_allsky_mask.dsrc.'//a9time//'.ppm'
        
         convert_cmd = 'convert -compress none '//trim(img_png)//' '//trim(img_ppm)
         write(6,*)trim(convert_cmd)
         call system(trim(convert_cmd))

         open(u,file=trim(img_ppm),status='old',err=999)
         read(u,*)   
         read(u,*)iwidth,iheight
         rewind(u)
         write(6,*)' dynamic dims ',ncol,iwidth,iheight
         ncol = 3
         call read_ppm(u,img_polar,ncol,iwidth,iheight)
         close(u)
         write(6,*)' polar mask (mode 1) has been read into img_polar array'

         do ic = 1,ncol
           write(6,*)
           write(6,*)' polar view for color ',ic
           do j = 1,nip,25
             write(6,21)ic,j,img_polar(ic,1:nip:25,j)
21           format(30i4)
           enddo ! j
         enddo ! ic

!        Convert to camera cloud mask image
         do i = 1,nip
         do j = 1,njp
             if(img_polar(3,i,j) .eq. 220)then
                 mask_polar(i,j) = 2 ! cloud
             elseif(img_polar(3,i,j) .eq. 255 .and. img_polar(2,i,j) .eq. 95)then
                 mask_polar(i,j) = 1 ! clear
             else
                 mask_polar(i,j) = 0 ! unknown
             endif
         enddo ! j
         enddo ! i

         do j = 1,njp,10
             write(6,*)
             write(6,*)' polar categorical mask'
             write(6,31)j,mask_polar(1:nip:6,j)
31           format(1x,i3,1x,100i1)
         enddo ! j

         write(6,*)' Sum of mask_polar is ',sum(mask_polar)

!        Convert to cylindrical image
         write(6,*)' Projecting to Cylindrical Mask using polar_to_cyl'
         call polar_to_cyl(mask_polar,mask_cyl,nip,njp,minalt,maxalt,minazi,maxazi,alt_scale,azi_scale)

       elseif(mode .eq. 2)then ! Read cyl mask (if produced by IDL code)
         if(.false.)then
           img_png = trim(camera_path)//'/camera_allsky_mask.dsrc.'//a9time//'.png'
           img_ppm = trim(camera_path)//'/camera_allsky_mask.dsrc.'//a9time//'.ppm'
        
           convert_cmd = 'convert -compress none '//trim(img_png)//' '//trim(img_ppm)
           write(6,*)trim(convert_cmd)
           call system(trim(convert_cmd))
           open(u,file=trim(img_ppm),status='old',err=999)
         else
           write(6,*)' Opening ',trim(fname_ppm)
           open(u,file=trim(fname_ppm),status='old',err=999)
         endif

         read(u,*)   
         read(u,*)iwidth,iheight
         rewind(u)
         ncol = 3
         write(6,*)' dynamic dims ',ncol,iwidth,iheight
         call read_ppm(u,img_polar,ncol,iwidth,iheight)
         close(u)
         write(6,*)' polar mask (mode 2) has been read into img_polar array'

         do ic = 1,ncol
           write(6,*)
           do j = 1,nip,25
             write(6,21)ic,j,img_polar(ic,1:nip:25,j)
           enddo ! j
         enddo ! ic

!        Convert to camera cloud mask image
         do i = 1,nip
         do j = 1,njp
             if(img_polar(3,i,j) .eq. 220)then
                 mask_polar(i,j) = 2 ! cloud
             elseif(img_polar(3,i,j) .eq. 255 .and. img_polar(2,i,j) .eq. 95)then
                 mask_polar(i,j) = 1 ! clear
             else
                 mask_polar(i,j) = 0 ! unknown
             endif
         enddo ! j
         enddo ! i

         do j = 1,njp,10
             write(6,31)j,mask_polar(1:nip:6,j)
         enddo ! j

         write(6,*)' Sum of mask_polar is ',sum(mask_polar)

!        Convert to cylindrical image
         write(6,*)' Projecting to Cylindrical Mask using polar_to_cyl'
         call polar_to_cyl(mask_polar,mask_cyl,nip,njp,minalt,maxalt,minazi,maxazi,alt_scale,azi_scale)

       endif

!      Mask out label at the top or clone nearby data (a small distance)

!      Mask out horizon area

!      Apply cloud algorithm
       do ialt = minalt,maxalt
       do jazi = minazi,maxazi
       enddo ! jazi
       enddo ! ialt

       istatus = 1
       write(6,*)' success in get_camera_clouds'
       return

999    istatus = 0
       write(6,*)' ERROR opening image file in read_camera_clouds'
       return

       end

       subroutine polar_to_cyl(img_p,img_c,nip,njp,minalt,maxalt,minazi,maxazi,alt_scale,azi_scale)

       include 'trigd.inc'

       integer img_p(nip,njp)
       integer img_c(minalt:maxalt,minazi:maxazi)

!      write(6,*)' sum of input polar is ',sum(img_p)

       img_c = 0

!      Loop through cylindrical array
       do ialt = minalt,maxalt
       do jazi = minazi,maxazi

!          Determine altitude and azimuth
           alt_d = float(ialt) * alt_scale
           azi_d = float(jazi) * azi_scale

!          Calculate Polar Image Coordinates (zero azimuth is up)
           ipcen = 256
           jpcen = 256
           radius = 255. * (90. - alt_d)/90.

           rip = float(ipcen) - sind(azi_d) * radius
           rjp = float(jpcen) - cosd(azi_d) * radius

           ip = nint(rip)
           jp = nint(rjp)

!          write(6,*)ialt,jazi,ip,jp,img_p(ip,jp)
           img_c(ialt,jazi) = img_p(ip,jp)

       enddo ! jazi
       enddo ! ialt

       return
       end

       subroutine get_camera_image(minalt,maxalt,minazi,maxazi,nc,alt_scale,azi_scale, & ! I
                                   i4time,fname_ppm,mode,site, &                         ! I
                                   rcam_rgb,istatus)                                     ! O
 
       use ppm

       parameter (nip = 511)
       parameter (njp = 511)

       real    rcam_rgb(nc,minalt:maxalt,minazi:maxazi) ! O
       integer icam_rgb(nc,minazi:maxazi,minalt:maxalt) ! L (flipped dimensions)

       integer img_polar(3,511,511)
       integer mask_polar(511,511)
       integer, allocatable :: img(:,:,:)
       integer u /12/ 

       character*275 imgfile,img_png,img_ppm,convert_cmd,imgdir
       character*13 a13name, cvt_i4time_wfo_fname13
       character*9 a9time
       character*10 fname_ppm
       character*10 site

!      Statement functions
       i2ialt(i) = minalt+i-1
       altf(ialt) = (ialt-minalt) * alt_scale

       j2jazi(i) = minazi+j-1
       azif(jazi) = (jazi-minazi) * azi_scale

!      Executable statements
       iverbose = 0

       a13name = cvt_i4time_wfo_fname13(i4time)
!      a13name = '20171017_1459'

       if(mode .eq. 1)then      ! Read polar observed image in PPM format 

!        imgdir = '/Users/albers/noaa/180309/wwwallsky/cases/17290/output3'
         imgdir = '/data/fab/projects/allsky/sites/dsrc/output3'
         img_png = trim(imgdir)//'/'//a13name//'_'//trim(site)//'.png'
         img_ppm = trim(imgdir)//'/'//a13name//'_'//trim(site)//'.ppm'
         convert_cmd = 'convert -compress none '//trim(img_png)//' '//trim(img_ppm)
         write(6,*)trim(convert_cmd)
         call system(trim(convert_cmd))

         open(u,file=trim(imgfile),status='old',err=999)
         read(u,*)   
         read(u,*)iwidth,iheight
         rewind(u)
         write(6,*)' dynamic dims ',ncol,iwidth,iheight
         allocate(img(ncol,iwidth,iheight))
         call read_ppm(u,img,ncol,iwidth,iheight)
         close(u)
         write(6,*)' polar image has been read in'

       elseif(mode .eq. 2)then ! Read cyl observed image in PPM format 
                               ! (remapped to alt/az grid)

         call get_directory('www',imgdir, imgdir_len)
         imgdir = trim(imgdir)//'als/'//trim(site)//'/observed/cyl'
!        imgdir = '/Users/albers/noaa/180323/wwwallsky/cases/17290/output3c'
!        imgdir = '/data/fab/projects/allsky/sites/dsrc/output3c'

         img_png = trim(imgdir)//'/'//a13name//'_'//trim(site)//'.png'
         img_ppm = trim(imgdir)//'/'//a13name//'_'//trim(site)//'.ppm'
!        convert_cmd = 'convert -crop 21x11+0+0 -compress none '//trim(img_png)//' '//trim(img_ppm)
         convert_cmd = 'convert -resize 50% -compress none '//trim(img_png)//' '//trim(img_ppm)
         write(6,*)trim(convert_cmd)
         call system(trim(convert_cmd))

         open(u,file=trim(img_ppm),status='old',err=999)
         read(u,*)   
         read(u,*)iwidth,iheight
         rewind(u)
         ncol = 3
         write(6,*)' dynamic dims ',ncol,iwidth,iheight
         if(iwidth .ne. maxazi-minazi+1)then
            write(6,*)' ERROR width discrepancy ',iwidth,maxazi-minazi+1
            istatus = 0
            return
         endif
         if(iheight .ne. maxalt-minalt+1)then
            write(6,*)' ERROR height discrepancy ',iheight,maxalt-minalt+1
            istatus = 0
            return
         endif

         call read_ppm(u,icam_rgb,ncol,iwidth,iheight)
         close(u)
         write(6,*)' cyl image has been read in'

         if(.true.)then ! clean ppm
           convert_cmd = 'rm -f '//trim(img_ppm)
           write(6,*)trim(convert_cmd)
           call system(trim(convert_cmd))
         endif

         if(iverbose .gt. 0)then
           do ih = 0,maxalt
             write(6,1)icam_rgb(:,:,ih)
1            format(10000(3i4))
           enddo

           write(6,*)' 0,0',minazi,maxazi
           write(6,*)icam_rgb(:,minazi,minalt)
           write(6,*)' minalt',minalt
           write(6,*)icam_rgb(1,minazi:maxazi:120,minalt)
           write(6,*)' maxalt',maxalt
           write(6,*)icam_rgb(1,minazi:maxazi:120,maxalt)
         endif
         
!        Rotate 90 degrees and flip
         do i = minalt,maxalt
           ialt_flip = minalt+(maxalt-i)
           do j = minazi,maxazi
             rcam_rgb(:,ialt_flip,j) = min(float(icam_rgb(:,j,i)),255.)
             if(iverbose .gt. 0)then
               write(6,15)j,i,min(icam_rgb(:,j,i),255),ialt_flip,j,rcam_rgb(:,ialt_flip,j)        
15             format('j/i -> ialt_flip/j',2i5,2x,3i4,5x,2i5,2x,3f5.0)
             endif
           enddo ! j
         enddo ! i

!      return ! test

         do ic = 1,ncol
           write(6,*)
           write(6,*)' cyl view for color ',ic,minazi,maxazi,minalt,maxalt
           do i = maxalt,minalt,-20
!            write(6,21)ic,i,nint(rcam_rgb(ic,i,minazi:maxazi:60))
!            write(6,21)ic,i,nint(rcam_rgb(ic,i,0:1140:60))
             write(6,21)ic,i,nint(rcam_rgb(ic,i,0:10))
21           format(30i4)
           enddo ! j
         enddo ! ic

         if(iverbose .gt. 0)then
           do i = maxalt,maxalt-5,-1
             write(6,*)' upper row ',i
             do j = minazi,maxazi
               write(6,22)i,j,nint(rcam_rgb(:,i,j))
22             format(i5,i5,2x,3i4)
             enddo ! j
           enddo ! i

           write(6,*)' 1st column'
           do i = maxalt,maxalt-40,-1
             write(6,*)i,nint(rcam_rgb(:,i,0))
           enddo ! i
         endif
       endif ! mode

       istatus = 1
       write(6,*)' success in get_camera_image'
       return

999    istatus = 0
       write(6,*)' ERROR opening image file in get_camera_image'
       return

       end
