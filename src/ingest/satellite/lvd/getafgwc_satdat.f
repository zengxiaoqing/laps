       subroutine getafgwc_satdat(isat,jtype,
     &                        max_channels,nchannels,chtype,
     &                        i4time_current,lvis_flag,
     &                        nirlines, nirelem,
     &                        nvislines,nviselem,
     &                        nwvlines,nwvelem,
     &                        ntm,max_files,c_type,
     &                        image_11,image_vis,
     &                        image_12,image_39,image_67,
     &                        i4time_data,
     &                        istatus)
c
c
c Note: if isat and jtype are both = 4, this indicates the
c       processing of GMS satellite data from AFWA. These settings
c       are defined in /data/static/satellite_lvd.nl.
c
       implicit none

       include 'satellite_dims_lvd.inc'
       include 'satellite_common_lvd.inc'

       integer isat,jtype
       integer ntm
       integer nirelem
       integer nirlines
       integer nwvelem
       integer nwvlines
       integer nviselem
       integer nvislines
       integer nchannels
       integer max_channels
       integer max_files
       integer i4time_data
       integer i4time_current
       integer istatus

       logical lvis_flag

       real image_11  (nirelem,nirlines)
       real image_12  (nirelem,nirlines)
       real image_39  (nirelem,nirlines)
       real image_67  (nwvelem,nwvlines)
       real image_vis (nviselem,nvislines)

       character c_type(max_files)*3
       character chtype(max_channels)*3

       istatus = 0

       if(isat.eq.4 .and. jtype.eq.4)then
          call getgmsdata(isat,jtype,
     &                    max_channels,nchannels,chtype,
     &                    i4time_current,lvis_flag,
     &                    nirlines, nirelem,
     &                    ntm,max_files,c_type,
     &                    image_11,image_vis,
     &                    image_12,           !image_39,
     &                    image_67,
     &                    i4time_data,
     &                    istatus)
          if(istatus.ne.0)then
             print*,'Error returned from getgmsdata'
          endif

       else

          call getgoesdata(isat,jtype,
     &                     max_channels,nchannels,chtype,
     &                     i4time_current,lvis_flag,
     &                     nirlines, nirelem,
     &                     nvislines,nviselem,
     &                     nwvlines,nwvelem,
     &                     ntm,max_files,c_type,
     &                     image_11,image_vis,
     &                     image_12,image_39,image_67,
     &                     i4time_data,
     &                     istatus)
          if(istatus.ne.0)then 
             print*,'Error returned from getgoesdata'
          endif

       endif

       return
       end
c
c =====================================================================
       subroutine getgoesdata(isat,jtype,
     &                        max_channels,nchannels,chtype,
     &                        i4time_current,lvis_flag,
     &                        nirlines, nirelem,
     &                        nvislines,nviselem,
     &                        nwvlines,nwvelem,
     &                        ntm,max_files,c_type,
     &                        image_11,image_vis,
     &                        image_12,image_39,image_67,
     &                        i4time_data,
     &                        istatus)
c
c Note: if isat and jtype are both = 4, then this indicates that
c       we are dealing with GMS satellite data from AFWA. Since
c       /data/static/satellite_lvd.nl will not change (other than
c       to add another satellite), this will remain a hardwire situtation.
c
       implicit none

       include 'satellite_dims_lvd.inc'
       include 'satellite_common_lvd.inc'

       integer isat,jtype
       integer i,j,n,nn
       integer nt
       integer ntm
       integer nirelem
       integer nirlines
       integer nwvelem
       integer nwvlines
       integer nviselem
       integer nvislines
       integer nchannels
       integer max_channels
       integer istart,jstart
       integer iend,jend
       integer ifiles_sat_raw
       integer ispec,istat
       integer lend
       integer jsave
       integer itm
       integer max_files

       integer i4time_current
       integer i4time_data
       integer i4time_data_io
       integer i4time_diff1
       integer i4time_diff2
       integer i4time_diff_min
       integer i4time_data_int(max_channels)

       integer istatus
       integer fstatus
       integer iostatus

       real image_11  (nirelem,nirlines) 
       real image_12  (nirelem,nirlines)
       real image_39  (nirelem,nirlines)
       real image_67  (nwvelem,nwvlines)
       real image_vis (nviselem,nvislines)

       real*4    grid_spacing_ir
       real*4    grid_spacing_wv
       real*4    grid_spacing_vis
       real*4    r4time_data_ir
       real*4    r4time_data_wv
       real*4    r4time_data_vis

       logical   lvis_flag
       logical   lfound1

       character c_type(max_files)*3
       character chtype(max_channels)*3
       character cid4*4
       character cmt*3
       character cname*100
       character csname*6
       character c_filetime*9
       character c_fname_data(max_channels)*9
       character cfd*9
       character path*200
       character cfname*11
       character c_afwa_fname*11
       character cfilename*255
c
      istatus = -1   !bad status return
      ntm=0

c     cid4='go'//c_sat_id(isat)(5:6)

      do i = 1,nchannels

         call lvd_file_specifier(chtype(i),ispec,istat)
         cfname=c_afwa_fname(c_sat_id(isat),chtype(i))
         n=index(path_to_raw_sat(ispec,jtype,isat),' ')-1
         cfilename=path_to_raw_sat(ispec,jtype,isat)(1:n)//cfname

         n=index(cfilename,' ')
         write(6,*)'Reading: ',cfilename(1:n)

         if(ispec.eq.1)then

            call read_afgwc_satdat(cfilename,isat,jtype,
     &l_cell_afwa,chtype(i),i_delta_sat_t_sec,i4time_current,
     &nvislines,nviselem,image_vis,i4time_data_io,iostatus)

         elseif(ispec.eq.2)then

            call read_afgwc_satdat(cfilename,isat,jtype,l_cell_afwa
     &,chtype(i),i_delta_sat_t_sec,i4time_current,nirlines,nirelem
     &,image_39,i4time_data_io,iostatus)

         elseif(ispec.eq.3)then

            call read_afgwc_satdat(cfilename,isat,jtype,l_cell_afwa
     &,chtype(i),i_delta_sat_t_sec,i4time_current,nwvlines,nwvelem
     &,image_67,i4time_data_io,iostatus)

         elseif(ispec.eq.4)then

            call read_afgwc_satdat(cfilename,isat,jtype,l_cell_afwa
     &,chtype(i),i_delta_sat_t_sec,i4time_current,nirlines,nirelem
     &,image_11,i4time_data_io,iostatus)

         elseif(ispec.eq.5)then

            call read_afgwc_satdat(cfilename,isat,jtype,l_cell_afwa
     &,chtype(i),i_delta_sat_t_sec,i4time_current,nirlines,nirelem
     &,image_12,i4time_data_io,iostatus)

         endif

         if(iostatus .eq. 1)then
            ntm=ntm+1
            c_type(ntm)=chtype(i)
            i4time_data_int(ntm)=i4time_data_io
            call make_fnam_lp(i4time_data_io,c_fname_data(ntm)
     &,fstatus)
            write(6,*)'gwc data loaded: ',c_type(ntm)
         else
            write(6,*)'No data loaded - ',chtype(i)
            goto 1000
         endif  

      enddo
c
c In this section we determine the i4time of the data and also check
c whether the i4times satisfy i_delta_sat_t_sec criterion.
c
c sort times in ascending order
c
      if(ntm.gt.2)then

         do j=2,ntm
            itm=i4time_data_int(j)
            cfd=c_fname_data(j)
            cmt=c_type(j)
            do i=j-1,1,-1
               if(i4time_data_int(i).le.itm)goto 100
               i4time_data_int(i+1)=i4time_data_int(i)
               c_fname_data(i+1)=c_fname_data(i)
               c_type(i+1)=c_type(i)
            enddo
            i=0
100         i4time_data_int(i+1)=itm
            c_fname_data(i+1)=cfd
            c_type(i+1)=cmt
         enddo

         lfound1=.false.
         jsave=0
         do j=1,ntm
            i4time_diff1=(i4time_data_int(j)-i4time_current)
            if(i4time_diff1.ge.i_delta_sat_t_sec.and.
     &(.not.lfound1))then
               jsave=j-1
               lfound1=.true.
            endif
         enddo

         i4time_diff_min=100000.
         if(jsave.gt.0)then
          ntm=jsave
          if(ntm.gt.1)then
           do i=2,ntm
            i4time_diff1=abs(i4time_data_int(i-1)-i4time_data_int(i))
            if(i4time_diff1.gt.600)then
               write(6,*)'Warning: Sat data not concurrent'
            endif
            if(i4time_diff1.lt.i4time_diff_min)then
               i4time_diff_min = i4time_diff1
               i4time_data = i4time_data_int(i)
            endif
           enddo
          else
           i4time_data=i4time_data_int(1)
          endif
         else
          i4time_data=i4time_data_int(1)
         endif

      elseif(ntm.eq.2)then

         i4time_diff1=i4time_data_int(1)-i4time_data_int(2)
         if(i4time_diff1.eq.0)then
            i4time_diff1=abs(i4time_data_int(1)-i4time_current)
            if(i4time_diff1.lt.i_delta_sat_t_sec)then
               i4time_data=i4time_data_int(1)
            else
               write(6,*)'Data too old'
               istatus=1
            endif
         else

            i4time_diff1=abs(i4time_current-i4time_data_int(1))
            i4time_diff2=abs(i4time_current-i4time_data_int(2))

            if(i4time_diff1.gt.i_delta_sat_t_sec)then
               if(i4time_diff2.gt.i_delta_sat_t_sec)then
                  istatus=1
                  write(6,*)'N files = ',ntm, 'Data too old'
               else
                  ntm=ntm-1
                  i4time_data=i4time_data_int(2)
                  c_type(ntm)=c_type(2)
               endif
            elseif(i4time_diff2.le.i_delta_sat_t_sec)then
               write(6,*)'Found ',ntm,' current files'
            else
               ntm=ntm-1
               i4time_data=i4time_data_int(1)
               c_type(ntm)=c_type(1)
            endif
         endif
      else   !ntm = 1
         i4time_data=i4time_data_int(1)
      endif

      istatus = 0  !successful return status
1000  return
      end
c
c ============================================================
c
      subroutine getgmsdata(isat,jtype,
     &                      max_channels,nchannels,chtype,
     &                      i4time_current,lvis_flag,
     &                      nlines, nelem,
     &                      ntm,max_files,c_type,
     &                      image_11,image_vis,
     &                      image_12,           !no image_39 data!,
     &                      image_67,
     &                      i4time_data,
     &                      istatus)
c
c
c
      implicit none

      include 'satellite_dims_lvd.inc'
      include 'satellite_common_lvd.inc'

      integer isat,jtype
      integer i,j,k,l,m
      integer nf,nc,np
      integer ierr
      integer ntm
      integer nelem
      integer nlines
      integer nchannels
      integer max_channels
      integer lvdindex,lstatus
      integer lend,lenf
      integer nfiles
      integer npixgms,nlingms
      integer max_files
      integer nxny(2)

      integer i4time_current
      integer i4time_data
      integer i4time_latest_lvd

      integer istatus
      integer fstatus
      integer gstatus

      real image_11  (nelem,nlines)
      real image_12  (nelem,nlines)
      real image_67  (nelem,nlines)
      real image_vis (nelem,nlines)
      integer image_gms (nelem*nlines+8)

      logical   lvis_flag

      character c_type(max_files)*3
      character chtype(max_channels)*3
      character cfname_cur*9
      character c_fname_data(max_channels)*9
      character cfd*9
      character cfiletime*9
      character cjjjhhmm*7
      character cjjjhr*5
      character ct*3
      character cpath*200
      character c_afwa_fname*100
      character cdir_lvd*150
      character cfilenames(max_files)*255
      character cname*100

      print*,'Subroutine getgmsdata'

      istatus=1

c get latest lvd file in lvd directory

      call get_directory('lvd',cdir_lvd,lend)
      call make_fnam_lp(i4time_current,cfname_cur,istatus)
      cdir_lvd=cdir_lvd(1:lend)//c_sat_id(isat)
      call get_file_names(cdir_lvd,nfiles,cfilenames,max_files
     +,gstatus)

      if(nfiles.gt.0)then
         call get_directory_length(cfilenames(nfiles),lenf)
         cname=cfilenames(nfiles)(1:lenf)//cfname_cur(1:5)//'*'
         call s_len(cname,lenf)
         call get_latest_file_time(cname,cfiletime)
         call i4time_fname_lp(cfiletime,i4time_latest_lvd,istatus)
      else
         i4time_latest_lvd=0
      endif

c     cjjjhr=cfname_cur(3:7)
      ntm=0
 
c hardwire for testing
c     cjjjhr='19607'
c     cjjjhr='15614'

      do i=1,nchannels

        if(chtype(i).eq.'vis'.and.lvis_flag)goto 90

        cpath=path_to_raw_sat(i,jtype,isat)
        np=index(cpath,' ')-1
        ct=chtype(i)
c       c_afwa_fname=cpath(1:np)//'GM5_'//cjjjhr//'*'//ct//'.unf'
        c_afwa_fname=cpath(1:np)//'GM5_'//'*'//ct//'.unf'
        np=index(c_afwa_fname,' ')-1
        call get_file_names(c_afwa_fname,nfiles,cfilenames,
     +max_files,istatus)

        if(istatus.ne.1)then
          print*,'Error status returned from get_file_names'
          print*,'Error getting filename for: ',c_afwa_fname(1:np)
          goto 1000
        elseif(nfiles.gt.0)then
          print*,'found ',nfiles,' matching ',c_afwa_fname(1:np)
          call s_len(cfilenames(nfiles),lenf)
          cjjjhr=cfilenames(nfiles)(lenf-19:lenf-13)
        else
          print*,'No files found ',cpath(1:np)
          goto 1000
        endif

        call lvd_file_specifier(ct,lvdindex,lstatus)
        if(lstatus.ne.0)goto 1000
c
c build filename and check to see if this time has already been
c processed. Assuming for the time being that nfiles should be = 1.
c

c       do j=1,nfiles

          call get_directory_length(cfilenames(nfiles),lenf)
          cjjjhhmm=cfilenames(nfiles)(lenf+5:lenf+11)
          cfd=cfname_cur(1:2)//cjjjhhmm
          call i4time_fname_lp(cfd,i4time_data,istatus)
          if(istatus.ne.1)then
            print*,'error returned from i4time_fname_lp'
            return
          endif

          if(i4time_data .gt. i4time_latest_lvd)then
c
            print*,'opening ',cfilenames(nfiles)

            call s_len(cfilenames(nfiles),lenf)
            call read_binary_field(nxny,4,4,2,cfilenames(nfiles),lenf)
            npixgms=nxny(1)
            nlingms=nxny(2)

            if(npixgms.ne.nelem.or.nlingms.ne.nlines)goto 900

            print*,'reading npix/nline gms ',npixgms,nlingms

c =====================================================================
            m=8
            call read_binary_field(image_gms,1,4,npixgms*nlingms+8,
     +cfilenames(nfiles),lenf)

          if(lvdindex.eq.1)then
            do k=1,nlingms
            do l=1,npixgms
               m=m+1
               image_vis(l,k)=float(image_gms(m))
            enddo
            enddo

          elseif(lvdindex.eq.3)then

            do k=1,nlingms
            do l=1,npixgms
               m=m+1
               image_67(l,k)=float(image_gms(m))
            enddo
            enddo

          elseif(lvdindex.eq.4)then

            do k=1,nlingms
            do l=1,npixgms
               m=m+1
               image_11(l,k)=float(image_gms(m))
            enddo
            enddo

          elseif(lvdindex.eq.5)then

            do k=1,nlingms
            do l=1,npixgms
               m=m+1
               image_12(l,k)=float(image_gms(m))
            enddo
            enddo

          endif

          ntm=ntm+1
          c_type(ntm)=ct
          call make_fnam_lp(i4time_data,c_fname_data(ntm),fstatus)

         endif
 
c       enddo

        goto 100

90    print*,'lvis_flag set = true '

100   enddo
c
c might want to put wait for data functionality in here (or near here).
c
      istatus=0
      goto 1000
900   print*,'Error: Mismatch between elem/lines for gms data'
      goto 1000
995   print*,'Error: reading file ',cfilenames(nfiles)(1:nf)
1000  return
      end
