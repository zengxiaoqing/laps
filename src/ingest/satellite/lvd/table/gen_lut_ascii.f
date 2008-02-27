
      Subroutine gen_ascii_lut(data_path,csatid,
     &csattype,chtype,nelem,nlines,nx_l,ny_l,lat,lon,
     &istatus)
c
c routine generates the ri/rj arrays for mapping laps grid points
c to satellite grid points. 
c
      implicit none

      integer nx_l,ny_l

      integer mxfiles
      parameter (mxfiles=500)

      integer i,j,k,l,n
      integer nc,nf,lend
      integer nlines
      integer nelem
      integer lvd_index
      integer nfiles
      integer i4time_current
      integer i4time_data
      integer istatus

      real      rlat(nelem,nlines)
      real      rlon(nelem,nlines)
      real      ri(nx_l,ny_l)
      real      rj(nx_l,ny_l)
      real      iline,jline
      real      idiff,jdiff

      real      grid_spacing_km

      real      lat(nx_l,ny_l)
      real      lon(nx_l,ny_l)

      logical   firsttime_for_ir
      data      firsttime_for_ir/.true./

      real   img_line(nelem,nlines)
      real   img_elem(nelem,nlines)
      real   image_data(nelem,nlines)

      integer gfn_status

      integer n_vars_req

      character     cname*14
      character     chtype*3
      character     csattype*3
      character     csatid*2
      character*9   c_filetime
      character*200 table_path
      character*150 cdir
      character*(*) data_path
      character*255 cfnames_asc(mxfiles)
c
c -------------------------------------------------------------------
      call get_file_names(data_path,
     &                    nfiles,
     &                    cfnames_asc,
     &                    mxfiles,
     &                    gfn_status)
      if(gfn_status.eq.1)then
         write(*,*)'Success GFN (lvd)'
      else
         write(6,*)'Error GFN (lvd)'
         istatus=-1
         goto 996
      endif
      nf=index(chtype,' ')-1
      if(nf.le.0)nf=3
      do i=1,nfiles
         j=255
         do while(j.gt.0.and.cfnames_asc(i)(j:j).ne.'_')
               j=j-1
         enddo
         if(j.gt.0)then 
            nc=index(cfnames_asc(i),' ')-1
            if(cfnames_asc(i)(j+1:nc).eq.chtype(1:nf))goto 10
         endif
      enddo
      goto 997
    
c     i4time_current=1142776920      !this set to the i4time corresponding
c                                     to date/time within the file LAPSI.ASC
c10    call lvd_file_specifier(chtype,lvd_index,istatus)
c      goto(1,2,3,2,2)lvd_index

10    if(firsttime_for_ir)then

         call  read_ascii_satdat(cfnames_asc(i),
     &                  i4time_current,
     &                  c_filetime,
     &                  nlines,nelem,
     &                  rlat,rlon,
     &                  img_line,img_elem,
     &                  image_data,
     &                  i4time_data,
     &                  grid_spacing_km,
     &                  istatus)
          if(chtype.eq.'11u'.or.chtype.eq.'12u'.or.
     &chtype.eq.'4u')then
             firsttime_for_ir=.false.
          endif

          call compute_rirj(nx_l,ny_l,nelem,nlines,
     &rlat,rlon,lat,lon,ri,rj,istatus)

      else

         write(6,*)'Already Computed ir lut'
         goto 900

      endif

      do j = 1,ny_l,10
      do i = 1,nx_l,10

         write(6,*)'i,j,ri,rj: ',i,j,ri(i,j),rj(i,j)

      enddo
      enddo

      cname='/lvd/sat-llij-'
      call get_directory('static',cdir,lend)
      table_path = cdir(1:lend)//cname//chtype(1:nc)//'-asc.lut'
      write(6,*)'Write lat/lon to i/j look up table'
      write(6,*)table_path(1:35)

      call write_table (table_path,nx_l,ny_l,lat,lon,ri,rj,istatus)
      if(istatus .ne. 1)then
         write(6,*)'Error writing look-up table'
         goto 900
      endif
c
c     call write_satsector_incfile(chtype,csattype,
c    &1,nlines,1,nelem,istatus)
c     if(istatus.ne.1)then
c        write(6,*)'Error write satsector_incfile'
c        goto 900
c     endif

c     write(6,*)'Sector file written'
c     write(6,*)'Done generating lut', chtype
c     write(6,*)
c
c ------------------------------------------------------------------------------
      goto 900

898   write(6,*)'Error get_static_info - path-to-raw-satellite'
      goto 900

901   write(6,*)'Error openning static/lvd/g8ir.parms file'
      goto 900

996   write(6,*)'Something wrong in get_file_names'
      goto 900

997   write(6,*)'No files of type: ',chtype

900   return
      end
