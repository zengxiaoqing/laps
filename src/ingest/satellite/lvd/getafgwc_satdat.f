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
       integer istatus_gfn

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
c first try for the ir data
c
      istatus = -1   !bad status return
      ntm=0

c     cid4='go'//c_sat_id(isat)(5:6)

      do i = 1,nchannels

         call lvd_file_specifier(chtype(i),ispec,istat)
         cfname=c_afwa_fname(c_sat_id(isat),chtype(i))
         n=index(path_to_raw_sat(ispec,jtype,isat),' ')-1
         cfilename=path_to_raw_sat(ispec,jtype,isat)(1:n)//cfname

c        goto(10,20,20,20,20)ispec   !(vis,3.9,wv,11.0,12.0)

c10          if(lvis_flag)goto 30
c           path=path_to_raw_sat(ispec,jtype,isat)
c           csname='u'//cid4//'v1'
c           n=index(path,' ')-1
c           cfname=path(1:n)//csname//'*_'//chtype(i)
c           goto 15

c20          csname='u'//cid4/'i1'
c           path=path_to_raw_sat(ispec,jtype,isat)
c           n=index(path,' ')-1
c           cfname=path(1:n)//csname//'*_'//chtype(i)
c
               n=index(cfilename,' ')
               write(6,*)'Reading: ',cfilename(1:n)

               goto(21,22,23,24,25)ispec

21             call read_afgwc_satdat(cfilename,isat,jtype,
     &chtype(i),i_delta_sat_t_sec,i4time_current,nvislines,nviselem,
     &image_vis,i4time_data_io,iostatus)
               goto 26

22             call read_afgwc_satdat(cfilename,isat,jtype,
     &chtype(i),i_delta_sat_t_sec,i4time_current,nirlines,nirelem,
     &image_39,i4time_data_io,iostatus)
               goto 26

23             call read_afgwc_satdat(cfilename,isat,jtype,
     &chtype(i),i_delta_sat_t_sec,i4time_current,nwvlines,nwvelem,
     &image_67,i4time_data_io,iostatus)
               goto 26

24             call read_afgwc_satdat(cfilename,isat,jtype,
     &chtype(i),i_delta_sat_t_sec,i4time_current,nirlines,nirelem,
     &image_11,i4time_data_io,iostatus)
               goto 26

25             call read_afgwc_satdat(cfilename,isat,jtype,
     &chtype(i),i_delta_sat_t_sec,i4time_current,nirlines,nirelem,
     &image_12,i4time_data_io,iostatus)

26             continue

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

30    enddo
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
               istatus=0
            endif
         else

            i4time_diff1=abs(i4time_current-i4time_data_int(1))
            i4time_diff2=abs(i4time_current-i4time_data_int(2))

            if(i4time_diff1.gt.i_delta_sat_t_sec)then
               if(i4time_diff2.gt.i_delta_sat_t_sec)then
                  istatus=0
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

      if(ntm.gt.0)istatus = 1

1000  return
      end
