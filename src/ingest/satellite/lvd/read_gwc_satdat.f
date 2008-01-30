      subroutine read_afgwc_satdat(c_filename,
     &                       isat,jtype,l_cell_afwa,
     &                       chtype,i_delta_t,
     &                       i4time_current,
     &                       nlines,nelems,
     &                       image_data,
     &                       i4time_data,
     &                       istatus)
c
      Implicit None

      Integer     cell_depth
      Integer     cell_width
      parameter  (cell_depth=64,
     &            cell_width=256)

      Integer     nlines,nelems

      Integer     n
      Integer     i4time_current
      Integer     i4time_data
      Integer     i4time_diff
      Integer     istatus
      Integer     i_delta_t
      Integer     isat,jtype
      Integer     nefi,nlfi
      Integer     lun

      INTEGER DECIMAT
      INTEGER STRPIX     !Start Pixel
      INTEGER STRSCNL    !Start Scanline 
      INTEGER STPPIX     !Stop Pixel
      INTEGER STPSCNL    !Stop Scanline
      INTEGER REQOBSTM   !Requested Observation Time
      INTEGER BEPIXFC
      INTEGER BESCNFC
      INTEGER FSCI
      INTEGER STRBDY1    !Requested Start Boundary 1
      INTEGER STRBDY2    !Requested Start Boundary 2
      INTEGER STPBDY1    !Requested Stop Boundary 1 
      INTEGER STPBDY2    !Requested Stop Boundary 2
      REAL   GOLONSBP      !GOES Longitude Subpoint
      REAL   GOLATSBP      !GOES Latitude Subpoint
      REAL   GOALPHA       !GOES Alpha

      REAL   image_data  (nelems,nlines)

      INTEGER WIDTH      !Tracks in Width of image
      INTEGER DEPTH      !Tracks in Depth of image
      CHARACTER*2 IMGTYPE  !Image Type

      Character     c_filetime*9
      Character     cfname*9
      Character     chtype*3
      Character     c_filename*255

      logical       lopen,lext
      logical       l_cell_afwa

      istatus=0
 
      n=index(c_filename,' ')-1

      inquire(file=c_filename,exist=lext,opened=lopen,number=lun)
      if(.not.lext)then
         print*,'File does not exist: ',c_filename(1:n)
         goto 1000
      endif
      if(lopen)then
         print*,'File is already open: ',c_filename(1:n)
         goto 1000
      endif
      call read_gwc_header(c_filename(1:n),l_cell_afwa,STRPIX,STRSCNL,
     +   STPPIX,STPSCNL,REQOBSTM,IMGTYPE,GOLATSBP,GOLONSBP,WIDTH,DEPTH, 
     +   GOALPHA,STRBDY1,STRBDY2,STPBDY1,STPBDY2,BEPIXFC,BESCNFC,
     +   FSCI,DECIMAT,istatus)
      if(istatus.eq.0)then
         write(6,*)'Got gwc header info'
      else
         write(6,*)'Error reading gwc header '
         goto 1000
      endif

c     write(c_filetime(1:7),111)reqobstm
c111   format(i7)
c     do i=1,7
c        if(c_filetime(i:i).eq.' ')then
c           c_filetime(i:i)='0'
c        endif
c     enddo
c     call make_fnam_lp(i4time_current,cfname,istatus)
c     c_yr = cfname(1:2)
c     c_jday = c_filetime(1:3)
c     c_hhmm = c_filetime(4:7)
c     c_filetime=c_yr//c_jday//c_hhmm
c     call cv_asc_i4time(c_filetime,i4time_data)

c
c put check in here for latest file time in lvd subdirectory.
c
      i4time_data=reqobstm

      i4time_diff = 0  !i4time_current-i4time_data    <---- Don't forget to un-comment.
      if(i4time_diff.gt.i_delta_t)then
         write(6,*)'Data is too old!'
         write(6,*)'Data    Time: ',c_filetime
         call make_fnam_lp(i4time_current,c_filetime,istatus)
         write(6,*)'Current Time: ',c_filetime
         goto 902
      elseif(i4time_diff.lt.0)then
         write(6,*)'Data time is more recent than current time!?'
         write(6,*)'This does not make sense!'
      else
         write(6,*)'Found AFWA data'
      endif
c
c read afgwc binary data file
c
      if(l_cell_afwa)then
         nlfi=depth*cell_depth
         nefi=width*cell_width
      else
         nlfi=stpscnl-strscnl+1
         nefi=stppix-strpix+1
      endif

      call Process_SDHS_GVAR_sub(c_filename,chtype,isat,jtype,
     &nelems,nlines,image_data,nlfi,nefi,depth,width,istatus)

      if(istatus.eq.1)then
         write(6,*)'Got AFWA satellite data'
      else
         write(6,*)'Error reading gwc sat data - Process_SDHS_sub'
      endif

      goto 1000

902   write(6,*)'Returning without new data'

1000  return
      end
