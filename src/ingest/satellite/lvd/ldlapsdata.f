      subroutine loadlapsdata(nx_l,ny_l,maxchannels,n_lvd_fields_max,
     &                        ntm,c_type,r_image_status,csatid,
     &                        image_vis,
     &                        image_39,image_67,
     &                        image_ir,image_12,
     &                        var_lvd,c_lvd,units_lvd,
     &                        nlf,laps_data,istatus)
c
c routine loads image satellite data that has already been
c navigated to the laps domain. Currently this is only from AFWA.
c
      implicit none

      integer   nx_l,ny_l
      integer   maxchannels
      integer   n_lvd_fields_max
      integer   i
      integer   ispec
      integer   ntm
      integer   nlf
      integer   istatus
      integer   lstatus

      real      r_image_status(maxchannels)
      real      laps_data(nx_l,ny_l,n_lvd_fields_max)
      real      image_vis(nx_l,ny_l)
      real      image_39(nx_l,ny_l)
      real      image_67(nx_l,ny_l)
      real      image_ir(nx_l,ny_l)
      real      image_12(nx_l,ny_l)

      character c_lvd(n_lvd_fields_max)*125
      character c_type(maxchannels)*3
      character var_lvd(n_lvd_fields_max)*3
      character units_lvd(n_lvd_fields_max)*10
      character csatid*(*)

      istatus=1
      do i=1,ntm

         call lvd_file_specifier(c_type(i),ispec,lstatus)
         if(lstatus.ne.0)goto 900

         goto(10,11,12,13,14)ispec

10       print*,'Not ready to deal with gms visible data'
         goto 99

c           call some-routine-to-deal-with-gms-visible-data
c           if(r_image_status(i).lt.0.3333)then
c           nlf=nlf+1
c           call move(visraw,laps_data(1,1,nlf),nx_l,ny_l)
c           var_lvd(nlf) = 'SVS'       ! satellite, visible
c           c_lvd(nlf)=csatid//' (VISIBLE) SATELLITE - RAW'
c           units_lvd(nlf) = 'COUNTS'
c and more
c
c goes only has 3.9u
c
11       if(csatid.ne.'gmssat'.or.csatid.ne.'metsat')then
            if(r_image_status(i).le.0.3333)then
               nlf=nlf+1
               call move(image_39,laps_data(1,1,nlf),nx_l,ny_l)
               var_lvd(nlf) = 'S3A'       ! satellite, , averaged
               c_lvd(nlf)=csatid//' (3.9u) IR B-TEMPS'
               nlf=nlf+1
               call move(image_39,laps_data(1,1,nlf),nx_l,ny_l)
               var_lvd(nlf)  = 'S3C'       ! satellite, , filtered
               c_lvd(nlf)=csatid//' (3.9u) IR B-TEMPS'
            else
               write(6,*)'39u image not processed: missing data'
            endif
         endif
         goto 99

12       if(r_image_status(i).le.0.3333)then
            nlf=nlf+1
            call move(image_67,laps_data(1,1,nlf),nx_l,ny_l)
            var_lvd(nlf) = 'S4A'       ! satellite, averaged
            c_lvd(nlf)=csatid//' (6.7u) IR B-TEMPS'
            nlf=nlf+1
            call move(image_67,laps_data(1,1,nlf),nx_l,ny_l)
            var_lvd(nlf) = 'S4C'       ! satellite, filtered
            c_lvd(nlf)=csatid//' (6.7u) IR B-TEMPS'
         else
            write(6,*)'wv image not processed: missing data'
         endif
         goto 99

13       if(r_image_status(i).lt.0.3333)then
            nlf=nlf+1
            call move(image_ir,laps_data(1,1,nlf),nx_l,ny_l)
            var_lvd(nlf)  = 'S8A'       ! satellite, channel-4, averaged
            c_lvd(nlf)=csatid//' (11.0u) IR B-TEMPS'
            nlf=nlf+1
            call move(image_ir,laps_data(1,1,nlf),nx_l,ny_l)
            var_lvd(nlf)='S8W'       ! satellite, channel-4, warm pixel
            c_lvd(nlf)=csatid//' (11.0u) IR B-TEMPS'
            nlf=nlf+1
            call move(image_ir,laps_data(1,1,nlf),nx_l,ny_l)
            var_lvd(nlf)='S8C'       ! satellite, channel-4, warm pixel
            c_lvd(nlf)=csatid//' (11.0u) IR B-TEMPS'
         else
            write(6,*)'IR image not processed: missing ir data'
         endif
         goto 99

14       if(r_image_status(i).lt.0.3333)then
            nlf=nlf+1
            call move(image_12,laps_data(1,1,nlf),nx_l,ny_l)
            var_lvd(nlf) = 'SCA'       ! satellite, averaged
            c_lvd(nlf)=csatid//' (12.0u) IR B-TEMPS'
            nlf=nlf+1
            call move(image_12,laps_data(1,1,nlf),nx_l,ny_l)
            var_lvd(nlf) = 'SCC'       ! satellite, averaged
            c_lvd(nlf)=csatid//' (12.0u) IR B-TEMPS'
         else
            write(6,*)'12u image not processed: missing data'
         endif

99    enddo

      istatus=0

900   print*,'error in lvd_file_specifier'

      return
      end
