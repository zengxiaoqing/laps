      subroutine write_lga(nx_laps,ny_laps,nz_laps,bgtime,bgvalid,
     .cmodel,missingflag,pr,ht,tp,sh,uw,vw,ww,istatus)

      implicit none

      integer nx_laps,ny_laps,nz_laps
      integer i,j,k,ic,lendir
      integer bgtime,bgvalid
      integer istatus
      integer warncnt
      integer ip(nz_laps)

      real    missingflag 
      real    pr(nz_laps)                  !Pressure or Sigma-P grid levels
      real    ht(nx_laps,ny_laps,nz_laps), !Height (m)
     .        tp(nx_laps,ny_laps,nz_laps), !Temperature (K)
     .        sh(nx_laps,ny_laps,nz_laps), !Specific humidity (kg/kg)
     .        uw(nx_laps,ny_laps,nz_laps), !U-wind (m/s)
     .        vw(nx_laps,ny_laps,nz_laps), !V-wind (m/s)
     .        ww(nx_laps,ny_laps,nz_laps)  !W-wind (pa/s)

      character*256 outdir
      character*31  ext
      character*3   var(nz_laps)
      character*4   lvl_coord(nz_laps)
      character*10  units(nz_laps)
      character*125 comment(nz_laps)
      character*5 fcst_hh_mm
      character*9 a9time
      character*256 outfile
      character*(*) cmodel

      character*40   v_g

      logical l_exist

      call get_vertical_grid(v_g,istatus)
      call upcase(v_g, v_g)

      ext='lga'
      call get_directory('lga',outdir,lendir)

c     Test for existence of output file
      call make_fnam_lp(bgtime,a9time,istatus)
      if (istatus .ne. 1) then
        write (6,*)
     1'Error converting i4time to file name...write aborted.'
        istatus=0
        return
      endif

      call make_fcst_time(bgvalid,bgtime,fcst_hh_mm,istatus)

      outfile = trim(outdir)//'/'//a9time//trim(fcst_hh_mm)//'.lga'
      inquire(file=outfile,exist=l_exist)
      if(l_exist .eqv. .false.)then
          write(6,*)' LGA file = ',trim(outfile)
      else
          write(6,*)' LGA file (already exists) = ',trim(outfile)
          write(6,*)
     1         ' Returning from write_lga without writing new LGA file'     
          istatus = 1
          return
      endif
c

      call s_len(cmodel,ic)
      do k=1,nz_laps                 ! Height
         do j=1,ny_laps
         do i=1,nx_laps
            if(ht(i,j,k) .ge. missingflag .and. warncnt.lt.100) 
     +              then
               print*,'Missingflag at ',i,j,k,' in ht'
               warncnt=warncnt+1
            endif
         enddo
         enddo
c only need to do some of these once.
         if (v_g .eq. 'PRESSURE') then     
            ip(k)=int(pr(k)) ! integer mb
            lvl_coord(k)='mb  '
            comment(k)=cmodel(1:ic)//' interpolated to LAPS isobaric.'
         elseif (v_g .eq. 'SIGMA_P') then
            ip(k)=nint(pr(k)*1000.) ! integer sigma * 1000
            lvl_coord(k)='    '
            comment(k)=cmodel(1:ic)//' interpolated to LAPS SIGMA P'
         endif

         var(k)='HT '
         units(k)='Meters'
      enddo
      print*,'HT'
      call write_laps(bgtime,bgvalid,outdir,ext,
     .                nx_laps,ny_laps,nz_laps,nz_laps,var,
     .                ip,lvl_coord,units,comment,ht,istatus)
      if (istatus .ne. 1) then
          print *,'Error writing interpolated data to LAPS database.'
          return
      endif

      do k=1,nz_laps                 ! Temperature
         do j=1,ny_laps
         do i=1,nx_laps
            if(tp(i,j,k) .ge. missingflag .and. warncnt.lt.100) 
     +              then
               print*,'Missingflag at ',i,j,k,' in tp'
               warncnt=warncnt+1
            endif
         enddo
         enddo
         var(k)='T3 '
         units(k)='Kelvin'
      enddo
      print*,'T3'
      call write_laps(bgtime,bgvalid,outdir,ext,
     .                nx_laps,ny_laps,nz_laps,nz_laps,var,
     .                ip,lvl_coord,units,comment,tp,istatus)
      if (istatus .ne. 1) then
          print *,'Error writing interpolated data to LAPS database.'
          return
      endif

      do k=1,nz_laps                 ! Specific humidity
         do j=1,ny_laps
         do i=1,nx_laps
            if(sh(i,j,k) .ge. missingflag .and. warncnt.lt.100) 
     +           then
               print*,'Missingflag at ',i,j,k,' in sh'
               warncnt=warncnt+1
            endif
         enddo
         enddo
         var(k)='SH '
         units(k)='kg/kg'
      enddo
      print*,'SH'
      call write_laps(bgtime,bgvalid,outdir,ext,
     .                nx_laps,ny_laps,nz_laps,nz_laps,var,
     .                ip,lvl_coord,units,comment,sh,istatus)
      if (istatus .ne. 1) then
          print *,'Error writing interpolated data to LAPS database.'
          return
      endif

      do k=1,nz_laps                 ! u-component wind
         do j=1,ny_laps
         do i=1,nx_laps
            if(uw(i,j,k) .ge. missingflag .and. warncnt.lt.100) 
     +              then
               print*,'Missingflag at ',i,j,k,' in uw'
               warncnt=warncnt+1
            endif
         enddo
         enddo
         var(k)='U3 '
         units(k)='m/s'
      enddo
      print*,'U3'
      call write_laps(bgtime,bgvalid,outdir,ext,
     .                nx_laps,ny_laps,nz_laps,nz_laps,var,
     .                ip,lvl_coord,units,comment,uw,istatus)
      if (istatus .ne. 1) then
          print *,'Error writing interpolated data to LAPS database.'
          return
      endif

      do k=1,nz_laps                 ! v-component wind
         do j=1,ny_laps
         do i=1,nx_laps
            if(vw(i,j,k) .ge. missingflag .and. warncnt.lt.100) 
     +              then
               print*,'Missingflag at ',i,j,k,' in vw'
               warncnt=warncnt+1
            endif
         enddo
         enddo
         var(k)='V3 '
      enddo
      print*,'V3'
      call write_laps(bgtime,bgvalid,outdir,ext,
     .                nx_laps,ny_laps,nz_laps,nz_laps,var,
     .                ip,lvl_coord,units,comment,vw,istatus)
      if (istatus .ne. 1) then
          print *,'Error writing interpolated data to LAPS database.'
          return
      endif
c
c if LAPS is background then we'll allow writing missing ww data.
c --------------------------------------------------------------
      if(cmodel.eq.'LAPS')then
         do k=1,nz_laps                 ! w-component wind
            do j=1,ny_laps
            do i=1,nx_laps
             if(ww(i,j,k) .ge. missingflag .and. warncnt.lt.50)
     +              then
               print*,'Missingflag at ',i,j,k,' in ww'
               warncnt=warncnt+1
             endif
            enddo
            enddo
            var(k)='OM '
            units(k)='pa/s'
         enddo


      else

         do k=1,nz_laps                 ! w-component wind
            do j=1,ny_laps
            do i=1,nx_laps
             if(ww(i,j,k) .ge. missingflag .and. warncnt.lt.100)
     +              then
               print*,'Missingflag at ',i,j,k,' in ww'
               warncnt=warncnt+1
             endif
            enddo
            enddo
            var(k)='OM '
            units(k)='pa/s'
         enddo
      endif

      print*,'OM'
      call write_laps(bgtime,bgvalid,outdir,ext,
     .                nx_laps,ny_laps,nz_laps,nz_laps,var,
     .                ip,lvl_coord,units,comment,ww,istatus)
      if (istatus .ne. 1) then
          print *,'Error writing interpolated data to LAPS database.'
          return
      endif

      return
      end
c
c -------------------------------------------------------------------------------
c
      subroutine write_lgap(nx_laps,ny_laps,nz_laps,bgtime,bgvalid,
     .cmodel,missingflag,ht,pr,tp,sh,uw,vw,ww,istatus)

      implicit none

      integer nx_laps,ny_laps,nz_laps
      integer i,j,k,ic,lendir
      integer bgtime,bgvalid
      integer istatus
      integer warncnt
      integer ip(nz_laps)

      real    missingflag
      real    ht(nz_laps)                  !SIGMA_HT height levels
      real    pr(nx_laps,ny_laps,nz_laps), !Pressure (Pa)
     .        tp(nx_laps,ny_laps,nz_laps), !Temperature (K)
     .        sh(nx_laps,ny_laps,nz_laps), !Specific humidity (kg/kg)
     .        uw(nx_laps,ny_laps,nz_laps), !U-wind (m/s)
     .        vw(nx_laps,ny_laps,nz_laps), !V-wind (m/s)
     .        ww(nx_laps,ny_laps,nz_laps)  !W-wind (m/s)

      real scale_height

      character*256 outdir
      character*31  ext
      character*3   var(nz_laps)
      character*4   lvl_coord(nz_laps)
      character*10  units(nz_laps)
      character*125 comment(nz_laps)
      character*(*) cmodel

      ext='lga'
      call get_directory('lga',outdir,lendir)

      call s_len(cmodel,ic)
      do k=1,nz_laps                 ! Pressure
         do j=1,ny_laps
         do i=1,nx_laps
            if(pr(i,j,k) .ge. missingflag .and. warncnt.lt.100) 
     +              then
               print*,'Missingflag at ',i,j,k,' in pr'
               warncnt=warncnt+1
            endif
         enddo
         enddo
c only need to do some of these once.
         ip(k)=nint(ht(k))
         var(k)='P3 '
         lvl_coord(k)='    '
         units(k)='Pascals'
         comment(k)=cmodel(1:ic)//' interpolated to LAPS sigma_ht.'
      enddo
      print*,'P3'
      call write_laps(bgtime,bgvalid,outdir,ext,
     .                nx_laps,ny_laps,nz_laps,nz_laps,var,
     .                ip,lvl_coord,units,comment,pr,istatus)
      if (istatus .ne. 1) then
          print *,'Error writing interpolated data to LAPS database.'
          return
      endif

      do k=1,nz_laps                 ! Temperature
         do j=1,ny_laps
         do i=1,nx_laps
            if(tp(i,j,k) .ge. missingflag .and. warncnt.lt.100) 
     +              then
               print*,'Missingflag at ',i,j,k,' in tp'
               warncnt=warncnt+1
            endif
         enddo
         enddo
         var(k)='T3 '
         units(k)='Kelvin'
      enddo
      print*,'T3'
      call write_laps(bgtime,bgvalid,outdir,ext,
     .                nx_laps,ny_laps,nz_laps,nz_laps,var,
     .                ip,lvl_coord,units,comment,tp,istatus)
      if (istatus .ne. 1) then
          print *,'Error writing interpolated data to LAPS database.'
          return
      endif

      do k=1,nz_laps                 ! Specific humidity
         do j=1,ny_laps
         do i=1,nx_laps
            if(sh(i,j,k) .ge. missingflag .and. warncnt.lt.100) 
     +           then
               print*,'Missingflag at ',i,j,k,' in sh'
               warncnt=warncnt+1
            endif
         enddo
         enddo
         var(k)='SH '
         units(k)='kg/kg'
      enddo
      print*,'SH'
      call write_laps(bgtime,bgvalid,outdir,ext,
     .                nx_laps,ny_laps,nz_laps,nz_laps,var,
     .                ip,lvl_coord,units,comment,sh,istatus)
      if (istatus .ne. 1) then
          print *,'Error writing interpolated data to LAPS database.'
          return
      endif

      do k=1,nz_laps                 ! u-component wind
         do j=1,ny_laps
         do i=1,nx_laps
            if(uw(i,j,k) .ge. missingflag .and. warncnt.lt.100) 
     +              then
               print*,'Missingflag at ',i,j,k,' in uw'
               warncnt=warncnt+1
            endif
         enddo
         enddo
         var(k)='U3 '
         units(k)='m/s'
      enddo
      print*,'U3'
      call write_laps(bgtime,bgvalid,outdir,ext,
     .                nx_laps,ny_laps,nz_laps,nz_laps,var,
     .                ip,lvl_coord,units,comment,uw,istatus)
      if (istatus .ne. 1) then
          print *,'Error writing interpolated data to LAPS database.'
          return
      endif

      do k=1,nz_laps                 ! v-component wind
         do j=1,ny_laps
         do i=1,nx_laps
            if(vw(i,j,k) .ge. missingflag .and. warncnt.lt.100) 
     +              then
               print*,'Missingflag at ',i,j,k,' in vw'
               warncnt=warncnt+1
            endif
         enddo
         enddo
         var(k)='V3 '
      enddo
      print*,'V3'
      call write_laps(bgtime,bgvalid,outdir,ext,
     .                nx_laps,ny_laps,nz_laps,nz_laps,var,
     .                ip,lvl_coord,units,comment,vw,istatus)
      if (istatus .ne. 1) then
          print *,'Error writing interpolated data to LAPS database.'
          return
      endif
c
c if LAPS is background then we'll allow writing missing ww data.
c --------------------------------------------------------------
      if(cmodel.eq.'LAPS')then
         do k=1,nz_laps                 ! w-component wind
            do j=1,ny_laps
            do i=1,nx_laps
             if(ww(i,j,k) .ge. missingflag .and. warncnt.lt.50)
     +              then
               print*,'Missingflag at ',i,j,k,' in ww'
               warncnt=warncnt+1
             endif

!            Convert from Omega to W (approximately - not yet coded)            

            enddo
            enddo
            var(k)='W3 '
            units(k)='m/s'
         enddo


      else

         do k=1,nz_laps                 ! w-component wind
            do j=1,ny_laps
            do i=1,nx_laps
             if(ww(i,j,k) .ge. missingflag)then
               if(warncnt.lt.100)
     +              then
                 print*,'Missingflag at ',i,j,k,' in ww'
                 warncnt=warncnt+1
               endif
             else
!              Convert from Omega to W (approximately)            
               scale_height = 8000.
               ww(i,j,k) = - (ww(i,j,k) / pr(i,j,k)) * scale_height

             endif

            enddo
            enddo
            var(k)='W3 '
            units(k)='m/s'
         enddo
      endif

      print*,'W3'
      call write_laps(bgtime,bgvalid,outdir,ext,
     .                nx_laps,ny_laps,nz_laps,nz_laps,var,
     .                ip,lvl_coord,units,comment,ww,istatus)
      if (istatus .ne. 1) then
          print *,'Error writing interpolated data to LAPS database.'
          return
      endif

      return
      end
c
c -------------------------------------------------------------------------------
c
      subroutine write_lgb(nx_laps,ny_laps,bgtime,bgvalid,cmodel
     .,missflag,uw_sfc,vw_sfc,tp_sfc,t_sfc,qsfc,pr_sfc,mslp,td_sfc
     .,rp_sfc,pcp_sfc,cw_sfc,istatus)

      implicit none

      integer nx_laps,ny_laps
      integer bgtime,bgvalid
      integer istatus,ic,len_dir
      integer i,j
      integer warncnt

      real      missflag
      real      tp_sfc(nx_laps,ny_laps),
     .          t_sfc (nx_laps,ny_laps),
     .          td_sfc(nx_laps,ny_laps),
     .          qsfc(nx_laps,ny_laps),
     .          uw_sfc(nx_laps,ny_laps),
     .          vw_sfc(nx_laps,ny_laps),
     .          pr_sfc(nx_laps,ny_laps),
     .          rp_sfc(nx_laps,ny_laps),
     .          pcp_sfc(nx_laps,ny_laps),
     .          cw_sfc(nx_laps,ny_laps),
     .          mslp(nx_laps,ny_laps)

      character*256 outdir
      character*31  ext
      character*3   var
      character*4   lvl_coord
      character*10  units
      character*125 comment
      character*(*) cmodel

      ext = 'lgb'
      call get_directory(ext,outdir,len_dir)
      call s_len(cmodel,ic)

c        print *,'writing to ',outdir(1:len_dir),fname//af(nf)(3:4)
c    &,'00.',ext
c        call write_laps(bgtime,bgvalid,outdir,ext
c    .           ,nx_laps,ny_laps,nsfc_fields,nsfc_fields
c    .           ,var,ip,lvl_coord
c    .           ,units,comment,sfcgrid,istatus)
c        if(istatus.eq.1)then
c           goto 88
c        else
c           print*,'Error writing interpolated data to LAPS lgb'
c           return
c        endif

      comment=cmodel(1:ic)//' interpolated to LAPS surface.'

c --- U
      warncnt=0
      do i=1,nx_laps
      do j=1,ny_laps
         if(uw_sfc(i,j).ge.missflag.and.warncnt.lt.100)
     +              then
            print*,'Missing data at ',i,j,' in uw_sfc'
            warncnt=warncnt+1
         endif
      enddo
      enddo
      units='m/s'
      var='USF'
      print*,'USF'
      call write_laps(bgtime,bgvalid,outdir,ext
     .           ,nx_laps,ny_laps,1,1,var,0,lvl_coord
     .           ,units,comment,uw_sfc,istatus)
      if (istatus .ne. 1) then
         print*,'Error writing interpolated data to LAPS lgb - U'
         return
      endif
c --- V
      warncnt=0
      do i=1,nx_laps
      do j=1,ny_laps
         if(vw_sfc(i,j).ge.missflag.and.warncnt.lt.100)
     +              then
            print*,'Missing data at ',i,j,' in vw_sfc'
            warncnt=warncnt+1
         endif
      enddo
      enddo
      var='VSF'
      print*,'VSF'
      call write_laps(bgtime,bgvalid,outdir,ext
     .           ,nx_laps,ny_laps,1,1,var,0,lvl_coord
     .           ,units,comment,vw_sfc,istatus)
      if (istatus .ne. 1) then
         print*,'Error writing interpolated data to LAPS lgb - V'
         return
      endif
c --- T
      warncnt=0
      do i=1,nx_laps
      do j=1,ny_laps
         if(tp_sfc(i,j).ge.missflag.and.warncnt.lt.100)
     +              then
            print*,'Missing data at ',i,j,' in tp_sfc'
            warncnt=warncnt+1
         endif
      enddo
      enddo
      units='K'
      var='TSF'
      print*,'TSF'
      comment=cmodel(1:ic)//' 2-m Temp FG'
      call write_laps(bgtime,bgvalid,outdir,ext
     .           ,nx_laps,ny_laps,1,1,var,0,lvl_coord
     .           ,units,comment,tp_sfc,istatus)
      if (istatus .ne. 1) then
         print*,'Error writing interpolated data to LAPS lgb - T'
         return
      endif
c --- T at Sfc
      warncnt=0
      do i=1,nx_laps
      do j=1,ny_laps
         if(t_sfc(i,j).ge.missflag.and.warncnt.lt.100)
     +              then
            print*,'Missing data at ',i,j,' in t_sfc'
            warncnt=warncnt+1
         endif
      enddo
      enddo
      units='K'
      var='TGD'
      print*,'TGD'
      comment=cmodel(1:ic)//' Ground Temp FG'
      call write_laps(bgtime,bgvalid,outdir,ext
     .           ,nx_laps,ny_laps,1,1,var,0,lvl_coord
     .           ,units,comment,t_sfc,istatus)
      if (istatus .ne. 1) then
         print*,'Error writing interpolated data to LAPS lgb - Tsfc'
         return
      endif
c --- q
      warncnt=0
      do i=1,nx_laps
      do j=1,ny_laps
         if(qsfc(i,j).ge.missflag.and.warncnt.lt.100)
     +              then
            print*,'Missing data at ',i,j,' in qsfc'
            warncnt=warncnt+1
         endif
      enddo
      enddo
      units='kg/kg'
      var='RSF'
      print*,'RSF'
      call write_laps(bgtime,bgvalid,outdir,ext
     .           ,nx_laps,ny_laps,1,1,var,0,lvl_coord
     .           ,units,comment,qsfc,istatus)
      if (istatus .ne. 1) then
         print*,'Error writing interpolated data to LAPS lgb - SH'
         return
      endif
c --- P sfc
      warncnt=0
      do i=1,nx_laps
      do j=1,ny_laps
         if(pr_sfc(i,j).ge.missflag.and.warncnt.lt.100)
     +              then
            print*,'Missing data at ',i,j,' in pr_sfc'
            warncnt=warncnt+1
         endif
      enddo
      enddo
      units='PA'
      var='PSF'
      print*,'PSF'
      call write_laps(bgtime,bgvalid,outdir,ext
     .           ,nx_laps,ny_laps,1,1,var,0,lvl_coord
     .           ,units,comment,pr_sfc,istatus)
      if (istatus .ne. 1) then
         print*,'Error writing interpolated data to LAPS lgb - PSF'
         return
      endif
c --- MSLP
      warncnt=0
      do i=1,nx_laps
      do j=1,ny_laps
         if(mslp(i,j).ge.missflag.and.warncnt.lt.100)
     +              then
            print*,'Missing data at ',i,j,' in mslp'
            warncnt=warncnt+1
         endif
      enddo
      enddo
      var='SLP'
      print*,'SLP'
      call write_laps(bgtime,bgvalid,outdir,ext
     .           ,nx_laps,ny_laps,1,1,var,0,lvl_coord
     .           ,units,comment,mslp,istatus)
      if (istatus .ne. 1) then
         print*,'Error writing interpolated data to LAPS lgb - SLP'
         return
      endif
c --- Td
      warncnt=0
      do i=1,nx_laps
      do j=1,ny_laps
         if(td_sfc(i,j).ge.missflag.and.warncnt.lt.100)
     +              then
            print*,'Missing data at ',i,j,' in td_sfc'
            warncnt=warncnt+1
         endif
      enddo
      enddo
      units='K'
      var='DSF'
      print*,'DSF'
      call write_laps(bgtime,bgvalid,outdir,ext
     .           ,nx_laps,ny_laps,1,1,var,0,lvl_coord
     .           ,units,comment,td_sfc,istatus)
      if (istatus .ne. 1) then
         print*,'Error writing interpolated data to LAPS lgb - DSF'
         return
      endif
c --- Reduced Pressure
      warncnt=0
      do i=1,nx_laps
      do j=1,ny_laps
         if(rp_sfc(i,j).ge.missflag.and.warncnt.lt.100)
     +              then
            print*,'Missing data at ',i,j,' in reduced press'
            warncnt=warncnt+1
         endif
      enddo
      enddo
      var='P'
      units='PA'
      print*,'P'
      call write_laps(bgtime,bgvalid,outdir,ext
     .           ,nx_laps,ny_laps,1,1,var,0,lvl_coord
     .           ,units,comment,rp_sfc,istatus)
      if (istatus .ne. 1) then
         print*,'Error writing interpolated data to LAPS lgb - P'
         return
      endif

c --- Precipitation 
      warncnt=0
      do i=1,nx_laps
      do j=1,ny_laps
         if(pcp_sfc(i,j).ge.missflag.and.warncnt.lt.100) ! modified rp_sfc -> pcp_sfc by Wei-Ting (130312)
     +              then
            print*,'Missing data at ',i,j,' in precip - pcp '
            warncnt=warncnt+1
         endif
      enddo
      enddo
      var='R01'
      units='M' 
      print*,'R01'
      call write_laps(bgtime,bgvalid,outdir,ext
     .           ,nx_laps,ny_laps,1,1,var,0,lvl_coord
     .           ,units,comment,pcp_sfc,istatus)
      if (istatus .ne. 1) then
         print*,'Error writing interpolated data to LAPS lgb - R01'
         return
      endif

c --- Cloud Water
      warncnt=0
      do i=1,nx_laps
      do j=1,ny_laps
         if(cw_sfc(i,j).ge.missflag.and.warncnt.lt.100)
     +              then
            print*,'Missing data at ',i,j,' in cw_sfc'
            warncnt=warncnt+1
         endif
      enddo
      enddo
      var='CW'
      units='M' 
      print*,'CW'
      call write_laps(bgtime,bgvalid,outdir,ext
     .           ,nx_laps,ny_laps,1,1,var,0,lvl_coord
     .           ,units,comment,cw_sfc,istatus)
      if (istatus .ne. 1) then
         print*,
     .      'Warning: not writing interpolated data to LAPS lgb - CW'
         istatus = 1 ! keep the program running normally
         return
      endif

      return
      end
