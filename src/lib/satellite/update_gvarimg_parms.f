      subroutine update_gvarimg_parms(cd6,
     &                        cstype,chtype,
     &                        cdir_path,
     &                        ewCycles,
     &                        ewIncs,
     &                        nsCycles,
     &                        nsIncs,
     &                        frameStartTime,
     &                        imc,
     &                        nw_vis_pix,
     &                        nw_vis_line,
     &                        orbitAttitude,
     &                        SatSubLAT,
     &                        SatSubLON,
     &                        decimat,
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

      Integer   i
      Integer   max_sat
      Integer   max_channels

      INTEGER   nw_vis_pix
      INTEGER   nw_vis_line
      INTEGER   se_vis_pix
      INTEGER   se_vis_line
      Integer   image_depth
      Integer   image_width
      Integer   strbdy1,strbdy2
      Integer   stpbdy1,stpbdy2
      Integer   decimat
      Integer   strtpix,strtline
      Integer   stoppix,stopline
      Integer   bepixfc,bescnfc,fsci
      REAL*8      frameStartTime 
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
      Integer   nch

      REAL*8      orbitAttitude                  (336)

      Integer   i4time_now_gg
      Integer   i4time_cur
      Integer   i4time_nearest
      Integer   i_obstime
      Integer   numoffiles
      Integer   nn,nf
      Integer   nsat
      Integer   n,ns,nc
      Integer   istatus
      Integer   gstatus
      Integer   lstatus
      Integer   nx,ny
      Integer   lend
      Integer   nstypes
      Integer   indx
      Integer   idum

      character*100 dir
      character*255 filename_cdf
      character*200 cdir_path
      character*255 fname_sat
      character*255 c_filenames(max_files)
      character*200 cdum
      character*10  cmode
      character*9   c_fname_cur
      character*9   c_fname
      character*3   chtype
      character*3   cstype
      character*6   cd6
      character*2   image_type
      character*1   ct
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
c === gwc switch ===
c
      elseif(cstype .eq. 'gwc')then

         nn=index(chtype,' ')-1
         if(nn.le.0)nn=3
         call lvd_file_specifier(chtype,indx,lstatus)
         goto(5,6,6,6,6)indx

5        ct='v'      !visible
         goto 8
6        ct='i'      !ir - either 3.9, 6.7, 11.2, or 12.0
8        continue

         fname_sat='u'//cd6(1:2)//cd6(5:6)//ct//'1_'//chtype(1:nn)
         nf=index(fname_sat,' ')
         n=index(cdir_path,' ')-1

         filename_cdf=cdir_path(1:n)//fname_sat(1:nf)

         call read_gwc_header(filename_cdf,strtpix,strtline,
     &stoppix,stopline,i_obstime,image_type,golatsbp,golonsbp,
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
         if(ct.eq.'i')then
            nw_vis_pix=(bepixfc+goalpha)*8
            nw_vis_line=(bescnfc+fsci)*4 
         elseif(ct.eq.'v')then
            nw_vis_pix= (bepixfc/4+goalpha)*8
            nw_vis_line=(bescnfc+fsci)
         endif
         nx = image_width*256
         ny = image_depth*64
         
         write(6,*)'GWC nw_vis_pix/nw_vis_line: ',nw_vis_pix,nw_vis_line
         write(6,*)

         filename_cdf=cdir_path(1:n)//'*.oad'

         call get_file_time(filename_cdf,i4time_cur,i4time_nearest)

         if(i4time_nearest.eq.0)then
            filename_cdf=cdir_path(1:n)//'971841336.oad'
         else
            call make_fnam_lp(i4time_nearest,c_fname,istatus)
            filename_cdf=cdir_path(1:n)//c_fname//'.oad'
         endif

         write(6,*)'Using: ',filename_cdf(1:100)
         call get_gwc_oa(filename_cdf,c_imc,orbitAttitude,336,
     &gstatus)
         if(gstatus.ne.0)then
            write(6,*)'Error in get_gwc_oanda'
            istatus=-1 
            goto 900
         endif
c
         frameStartTime=getftime()
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
     &rlat,rlon,dx,dy,nx,ny,istatus)

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
         istatus = 0
         goto 900
      endif
c
c read header of current 11u gvarimage
c
      write(6,*)'Reading wfo cdf header = ',chtype(1:nn)
      call get_attribute_wfo(cfname_sat,rlat,rlon,dx,dy,nx,ny,lstatus)
      if(lstatus .ne. 0)then
         write(6,*)'Error getting wfo attributes = ',chtype(1:nn)
         goto 900
      endif

      dx=dx*1000.0
      dy=dy*1000.0

      istatus = 1
      goto 1000

900   write(6,*)'Returning without new attributes'
      goto 1000

997   write(6,*)'Error returned from make_fnam_lp'

1000  return
      end
