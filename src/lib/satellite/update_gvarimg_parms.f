      subroutine update_gvarimg_parms(cd6,
     &                        cstype,l_cell_afwa,chtype,
     &                        cdir_path,
     &                        ewCycles,
     &                        ewIncs,
     &                        nsCycles,
     &                        nsIncs,
     &                        frameStartTime,
     &                        imc,x_res,y_res,
     &                        nw_vis_pix,
     &                        nw_vis_line,
     &                        orbitAttitude,
     &                        SatSubLAT,
     &                        SatSubLON,
     &                        x_step,y_step,
     &                        nx,ny,
     &                        istatus)
c
c Program reads the netcdf gvarimage satellite data file header
c or the AFWA file header and returns the relevant information for updating
c the satellite namelist (/static/satellite.nl).
c
      implicit none

      Integer   max_files
      parameter   (max_files = 500)

c     Integer   i
c     Integer   max_sat
c     Integer   max_channels

      INTEGER   nw_vis_pix
      INTEGER   nw_vis_line
      INTEGER   nw_vis_pix_gwc
      INTEGER   nw_vis_line_gwc
      INTEGER   se_vis_pix
      INTEGER   se_vis_line
      Integer   image_depth
      Integer   image_width
      Integer   strbdy1,strbdy2
      Integer   stpbdy1,stpbdy2
      Integer   strtpix,strtline
      Integer   stoppix,stopline
      Integer   bepixfc,bescnfc,fsci
      Integer   decimat
      Integer   x_step,y_step

      Real*8      frameStartTime 
      Real*8      getftime
      Real*8      SatSubLAT,SatSubLON
      Real*4      golonsbp
      Real*4      golatsbp
      Real*4      goalpha
      CHARACTER*1 imc                            (4)
      CHARACTER*4 c_imc
      INTEGER   ewCycles                 
      INTEGER   ewIncs                   
      INTEGER   nsCycles               
      INTEGER   nsIncs                 
      Integer   x_res,y_res
      Integer   imci4
c     Integer   nch

      logical   l_cell_afwa

      REAL*8      orbitAttitude                  (336)

      Integer   i4time_now_gg
      Integer   i4time_cur
      Integer   i4time_nearest
      Integer   i_obstime
      Integer   ld
      Integer   numoffiles
      Integer   nf,nl,nn
c     Integer   nsat
c     Integer   n,ns,nc
      Integer   n
      Integer   istatus
      Integer   gstatus
      Integer   lstatus
      Integer   nx,ny
c     Integer   lend
c     Integer   nstypes
c     Integer   indx
c     Integer   idum

c     character*150 dir
      character*150 c_filespec
      character*255 filename_cdf
      character*200 cdir_path
      character*255 fname_sat
      character*255 c_filenames(max_files)
c     character*200 cdum
      character     cfname*100
      character     c_afwa_fname*100
c     character*10  cmode
      character*9   c_fname_cur
      character*9   c_fname
      character*3   chtype
      character*3   cstype
      character*6   cd6
      character*2   image_type
c     character*1   ct
c
c ========================================================
c
      istatus = -1  !bad status return

      i4time_cur = i4time_now_gg()
      call make_fnam_lp(i4time_cur,c_fname_cur,lstatus)
      if(lstatus.ne.1)goto 997
c
c find latest gvarimage data.  cstype = 'gvr' is FSL GVAR data
c
      if(cstype.eq.'gvr')then

         nn=index(chtype,' ')-1
         if(nn.le.0)nn=3

         n=index(cdir_path,' ')-1
         fname_sat=cdir_path(1:n)//'*_'//chtype(1:nn)
         call get_file_time(fname_sat,i4time_cur,i4time_nearest)
c
c just in case there is no current data, skip this channel.
c
         if(i4time_nearest.gt.0)then

            call make_fnam_lp(i4time_nearest,c_fname,lstatus)
            if(lstatus.ne.1)goto 997
            fname_sat=cdir_path(1:n)//c_fname//'*_'//chtype(1:nn)

            call get_file_names(fname_sat,numoffiles,c_filenames
     1,max_files,lstatus)
            if(lstatus.ne.1)goto 998
            if(numoffiles.gt.0)then
               if(numoffiles.eq.1)then
                  filename_cdf=c_filenames(numoffiles)
               else
                  write(6,*)'Odd, numoffiles > 1: = ',numoffiles
                  filename_cdf=c_filenames(1)
               endif
            else
               nf=index(filename_cdf,' ')
               write(6,*)'NO files found ',filename_cdf(1:nf)
               goto 899
            endif

         else

            nf=index(cdir_path,' ')
            write(6,*)'No data in ',cdir_path(1:nf)
            istatus = 1  !status flag indicating gvar header not read
            goto 900

         endif
c
c read header of current 11u gvarimage
c
         write(6,*)'Reading gvar cdf header = ',chtype(1:nn)

         call rd_gvarimg_cdf_header(filename_cdf,
     &                        nw_vis_pix,
     &                        nw_vis_line,
     &                        se_vis_pix,
     &                        se_vis_line,
     &                        ewCycles,
     &                        ewIncs,
     &                        nsCycles,
     &                        nsIncs,
     &                        x_res,y_res,
     &                        frameStartTime,
     &                        imc,
     &                        orbitAttitude,
     &                        nx,ny,
     &                        x_step,y_step,
     &                        lstatus)

         if(lstatus .ne. 1)then
            write(6,*)'Error reading gvar cdf header ',chtype(1:nn)
            goto 996
         else
            write(6,*)'Header info obtained = ',chtype(1:nn)
         endif
c
c compute satellite sub latitude and sub longitude.
c
         read(imc(4),100)imci4
100      format(i1)
cc
c the bottom line:
c
         imci4 = 0

         write(6,*)'Compute SatSubLAT/SatSubLON - enter subroutine'
         call sat_sublatlon(ewCycles,ewIncs,nsCycles,nsIncs,
     &frameStartTime,imci4,orbitAttitude,SatSubLAT,SatSubLON,lstatus)
         if(lstatus.ne.1)then
            goto 999
         endif
c
c ====================== gwc switch =======================
c
      elseif(cstype .eq. 'gwc')then

         n=index(cdir_path,' ')-1
         cfname=c_afwa_fname(cd6,chtype)
         call s_len(cfname,nl)
         filename_cdf=cdir_path(1:n)//cfname(1:nl)

         call read_gwc_header(filename_cdf,l_cell_afwa,strtpix,
     &strtline,stoppix,stopline,i_obstime,image_type,golatsbp,golonsbp,
     &image_width,image_depth,goalpha,strbdy1,strbdy2,stpbdy1,
     &stpbdy2,bepixfc,bescnfc,fsci,decimat,gstatus)
         if(gstatus.ne.0)then
            write(6,*)'Error in read_gwc_header'
            istatus=-1 
            goto 900
         endif

         SatSubLAT=golatsbp
         SatSubLON=golonsbp
c
c no water vapor switch atm
c
         if(cd6.eq.'meteos')then
            nw_vis_pix=stoppix
            nw_vis_line=-(fsci+strtline)
         else
            nw_vis_pix = nw_vis_pix_gwc(chtype,decimat,bepixfc,goalpha)
            nw_vis_line= nw_vis_line_gwc(chtype,decimat,bescnfc,fsci)
         endif

c        if(chtype.eq.'vis')then
c           nw_vis_pix= (bepixfc/4+goalpha)*8
c           nw_vis_line=(bescnfc+fsci)
c        else
c           nw_vis_pix=(bepixfc+goalpha)*8
c           nw_vis_line=(bescnfc+fsci)*4 
c        endif

         if(l_cell_afwa)then
            nx = image_width*256
            ny = image_depth*64
         else
            nx = stoppix-strtpix+1
            ny = stopline-strtline+1
         endif
         
         write(6,*)'GWC nw_vis_pix/nw_vis_line: ',nw_vis_pix,nw_vis_line
         write(6,*)

c Currently no O&A data for METEOSAT
         if(cd6.ne.'meteos')then

           filename_cdf=cdir_path(1:n)//'*_OA_01.DAT'
           nf=index(filename_cdf,' ')-1
           call get_file_names(filename_cdf,numoffiles,c_filenames
     1        ,max_files,istatus)

           if(istatus.eq.1.and.numoffiles.gt.0)then

           call get_gwc_oa(c_filenames(1),c_imc,orbitAttitude,336,
     &gstatus)
           if(gstatus.ne.0)then
             write(6,*)'error: get_gwc_oa '
             call get_directory('static',c_filespec,ld)
             c_filespec=c_filespec(1:ld)//'/lvd'
             ld=index(c_filespec, ' ')-1
             print*,'try ',c_filespec(1:ld),'lvd/',cd6,'_orbatt.dat'
             call read_orb_att(c_filespec(1:ld),cd6,336,orbitAttitude,
     &istatus)
             if(istatus.ne.0)then
                write(6,*)'O&A Data not obtained',c_filespec(1:ld)
                goto 900
             endif
           else
c            call make_fnam_lp(i4time_nearest,c_fname,istatus)
c            filename_cdf=cdir_path(1:n)//c_fname//'.oad'
c            ld=index(c_filespec, ' ')-1
c            write(6,*)'Using: ',filename_cdf(1:ld)
             write(6,*)'gwc O&A obtained '
           endif

           else

             write(6,*)'No O&A files exist ',filename_cdf(1:nf)
             call get_directory('static',c_filespec,ld)
             c_filespec=c_filespec(1:ld)//'/lvd'
             ld=index(c_filespec, ' ')-1
             call read_orb_att(c_filespec(1:ld),cd6,336,orbitAttitude,
     &istatus)
             if(istatus.ne.0)then
                write(6,*)'O&A Data not obtained',c_filespec(1:ld)
                goto 900
             endif
           endif
         endif
         frameStartTime=getftime()
         x_step=decimat
         y_step=decimat
c
      endif
      istatus = 0

      goto 900
 
101   write(6,*)'Error getting ri/rj luts'
      goto 900

996   write(6,*)'error in rd_gvarimg_cdf_header: terminating'
      istatus = 1
      goto 900

997   write(6,*)'error in make_fname_lp'
      goto 900

998   write(6,*)'error in get_file_names'
      goto 900

999   write(6,*)'error in satsublatlon'
      istatus = 1
      goto 900

899   write(6,*)'Parmfile not updated'

900   return
      end

c ===============================================================
c
      subroutine get_wfo_nav_parms(cpath,chtype,
     &rlat,rlon,rlatnxny,rlonnxny,rlatdxdy,rlondxdy,
     &dx,dy,nx,ny,istatus)

      character*(*) cpath
      character     cfname_sat*200
      character     fname9_to_wfo_fname13*13
      character     cfname9*9
      character     cfname13*13
      character     chtype*3
      integer       istatus
      integer       i4time_cur
      integer       nx,ny
      real*4        rlat
      real*4        rlon
      real*4        rlatdxdy
      real*4        rlondxdy
      real*4        rlatnxny
      real*4        rlonnxny
      real*4        dx,dy

      istatus = -1  !bad status return

      i4time_cur = i4time_now_gg()
c
c find latest wfo satellite data.
c
      nn=index(chtype,' ')-1
      if(nn.le.0)nn=3
      n=index(cpath,' ')-1
      cfname_sat=cpath(1:n)//'*'

      call get_file_time(cfname_sat,i4time_cur,i4time_nearest)

      if(i4time_nearest.gt.0)then
         call make_fnam_lp(i4time_nearest,cfname9,lstatus)
         cfname13=fname9_to_wfo_fname13(cfname9)
         if(lstatus.ne.1)goto 997
         cfname_sat=cpath(1:n)//cfname13
      else
         write(6,*)'No data in ',cpath(1:n)
         goto 900
      endif
c
c read header of current 11u gvarimage
c
      write(6,*)'calling get_attribute_wfo ',chtype(1:nn)
      call get_attribute_wfo(cfname_sat,rlat,rlon,rlatnxny,rlonnxny,
     &rlatdxdy,rlondxdy,dx,dy,nx,ny,lstatus)
      if(lstatus .lt. 0)then
         write(6,*)'No attributes returned: get_attribute_wfo ',chtype
     &(1:nn)
         goto 900

c     elseif(lstatus.gt.0)then
c        print*,'No attributes returned from get_attribute_wfo'
c        istatus=1
c        goto 900

      endif

      dx=dx*1000.0
      dy=dy*1000.0

      istatus = 0
      goto 1000

900   write(6,*)'Returning without new attributes'
      goto 1000

997   write(6,*)'Error returned from make_fnam_lp'

1000  return
      end
