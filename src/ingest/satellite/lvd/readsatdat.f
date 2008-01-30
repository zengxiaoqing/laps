      subroutine readsatdat(csat_id,
     &                      csat_type, 
     &                      c_dir_path,
     &                      c_fname_data,
     &                      c_type,
     &                      ntf,ntm,
     &                      maxchannels,
     &                      max_files,
     &                      maximages,
     &                      nlineswv,nelemwv,
     &                      nlinesir,nelemir,
     &                      nlinesvis,nelemvis,
     &                      image_67,
     &                      image_ir,image_12,
     &                      image_39,image_vis,
     &                      istatus)
c
c
c
      implicit none

      include 'netcdf.inc'

      integer i,j,n,id
      integer ntf
      integer nlinesir,nelemir
      integer nlinesvis,nelemvis
      integer nlineswv,nelemwv
      integer max_files,maximages,maxchannels
      integer ntm(max_files)
      integer istatus
      integer istatus_ir
      integer istatus_vis
      integer istatus_wv
      integer wstatus
      integer istat
      integer ispec

      real    image_ir(nelemir,nlinesir,maximages)
      real    image_12(nelemir,nlinesir,maximages)
      real    image_39(nelemir,nlinesir,maximages)
      real    image_67(nelemwv,nlineswv,maximages)
      real    image_vis(nelemvis,nlinesvis,maximages)
      real    wv_image(nelemwv,nlineswv)
      real    ir_image(nelemir,nlinesir)
      real    vis_image(nelemvis,nlinesvis)

      character c_fname_data(max_files)*9
      character c_type(maxchannels,max_files)*3
      character c_dir_path(maxchannels)*200
      character c_filename*200
      character c_wfo_fname*13
      character fname9_to_wfo_fname13*13
      character csat_id*6
      character csat_type*3
c
      INTEGER   record
      INTEGER   ncid
      INTEGER   rcode
      INTEGER   ivalidTime
      doubleprecision validtime
      REAL    dummy
      REAL    la1_vis,la1_ir,la1_wv
      REAL    lo1_vis,lo1_ir,lo1_wv
      REAL    dx_vis,dx_ir,dx_wv
      REAL    dy_vis,dy_ir,dy_wv
      REAL    latin_vis,latin_ir,latin_wv
      REAL    lov_vis,lov_ir,lov_wv

      istatus=1

      record=1

      do i=1,ntf
         do j=1,ntm(i)

            call lvd_file_specifier(c_type(j,i),ispec,istat)

            if(csat_type.eq.'wfo'.or.csat_type.eq.'ncp')then
               n=index(c_dir_path(ispec),' ')-1
               c_wfo_fname = fname9_to_wfo_fname13(c_fname_data(i))
               c_filename=c_dir_path(ispec)(1:n)//c_wfo_fname
               n=index(c_filename,' ')
            else
               n=index(c_dir_path(1),' ')-1
               c_filename=c_dir_path(1)(1:n)//c_fname_data(i)//
     &'_'//c_type(j,i)
            endif
            n=index(c_filename,' ')
            print*,'Reading: ',c_filename(1:n)

            rcode=NF_OPEN(c_filename,NF_NOWRITE,NCID)
            if(rcode.ne.nf_noerr) return

            if(ispec.ne.1.and.ispec.ne.3)then    !check for visible and water vapor

               call readcdf(csat_id,
     &                    csat_type,
     &                    c_type(j,i),
     &                    record,
     &                    nelemir,nlinesir,
     &                    ir_image,
     &                    la1_ir,lo1_ir,
     &                    Dx_ir,Dy_ir,
     &                    Latin_ir,Lov_ir,
     &                    ivalidTime,
     &                    ncid,
     &                    istatus_ir)
               if(istatus_ir.eq.1)then
                  write(6,*)'Successful'
               else
                  Write(6,*)'NOT Successful!!'
                  istatus=-1
                  goto 125
               endif


cisido          
		write(6,*)'la1_ir',la1_ir,'lo1_ir',lo1_ir,
     &                     'Dx_ir',Dx_ir,'Dy_ir',Dy_ir,
     &                     'Latin_ir',Latin_ir,'Lov_ir',Lov_ir 
cisid

c
               if(ispec .eq. 2)then 
                  call move(ir_image,image_39(1,1,i),nelemir,nlinesir)
               elseif(ispec .eq. 4)then
                  call move(ir_image,image_ir(1,1,i),nelemir,nlinesir)
               elseif(ispec .eq. 5)then
                  call move(ir_image,image_12(1,1,i),nelemir,nlinesir)
               endif
c
            elseif(ispec.eq.1)then

               call readcdf(csat_id,
     &                    csat_type,
     &                    c_type(j,i),
     &                    record,
     &                    nelemvis,nlinesvis,
     &                    vis_image,
     &                    la1_vis,lo1_vis,
     &                    Dx_vis,Dy_vis,
     &                    Latin_vis,Lov_vis,
     &                    ivalidTime,
     &                    ncid,
     &                    istatus_vis)
               if(istatus_vis.eq.1)then
                  write(6,*)'Successful'
                  call move(vis_image,image_vis(1,1,i),nelemvis,
     &                      nlinesvis)
               else
                  Write(6,*)'NOT Successful!!'
                  istatus=-1
               endif

cisido
               write(6,*)'la1_vis',la1_vis,'lo1_vis',lo1_vis,
     &                   'Dx_vis',Dx_vis,'Dy_vis',Dy_vis,
     &                   'Latin_vis',Latin_vis,'Lov_vis',Lov_vis
cisid

c
c load water vapor attributes
c
            elseif(ispec.eq.3)then

               call readcdf(csat_id,
     &                    csat_type,
     &                    c_type(j,i),
     &                    record,
     &                    nelemwv,nlineswv,
     &                    wv_image,
     &                    la1_wv,lo1_wv,
     &                    Dx_wv,Dy_wv,
     &                    Latin_wv,Lov_wv,
     &                    ivalidTime,
     &                    ncid,
     &                    istatus_wv)
               if(istatus_wv.eq.1)then
                  write(6,*)'Successful'
                  call move(wv_image,image_67(1,1,i),nelemwv,
     &                      nlineswv)
               else
                  Write(6,*)'NOT Successful!!'
                  istatus=-1
               endif

cisido
                write(6,*)'la1_wv',la1_wv,'lo1_wv',lo1_wv,
     &                     'Dx_wv',Dx_wv,'Dy_wv',Dy_wv,
     &                     'Latin_wv',Latin_wv,'Lov_wv',Lov_wv
cisid





c
c load vis attributes
c
            endif

            rcode= NF_CLOSE(ncid)

125      enddo
      enddo
c
      return
      end
