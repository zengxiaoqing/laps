 
      subroutine put_laps_static(grid_spacing,model,comment,var
     1                          ,data,imax,jmax,mkmax,ngrids
     1                          ,std_lat,std_lat2,std_lon      
     1                          ,c6_maproj,deltax,deltay,center_lat
     1                          ,center_lon,lli,llj,uri,urj
     1                          ,parent_id,ratio_2_parent,nest
     1                          ,c10_grid_fname)
 
      include 'trigd.inc' 

      integer   mkmax                   !This should be max as defined in calling routine.
      integer   ngrids

      INTEGER  	IMAX,		        !I4time of data
     1		JMAX,KMAX,	        !# cols, # rows, # fields
     1		KDIM,i,j,k,		!K dimension of DATA array
     1		ISTATUS
C
      REAL  	DATA(imax,jmax,mkmax),	!Raw data to be written
     1		grid_spacing,val,
     1          center_lat,center_lon,
     1          la1,lo1,la2,lo2
C
      CHARACTER*150     DIR_out		!Directory to be written to
      CHARACTER*31	EXT		!File name ext (up to 31 chars)
      CHARACTER*3	VAR(mkmax)	!3 letter ID of each field
      CHARACTER*10	UNITS(mkmax)	!units of each field
      CHARACTER*(*)	COMMENT(mkmax)	!Comments for each field
      character*(*)     model
      character*80      origin          !Run time parameter - c80_description
      character*10  c10_grid_f          ! Type of domain (nest7grid, wrfsi)
      character*10  c10_grid_fname      ! Actual filename, cdl (nest7grid for now)
!mp      character*9       laps_dom_file
      character*16       laps_dom_file
      character*6       c6_maproj       !Map projection
      character*200     cdataroot       !Dummy Variable used in find_domain_name.
      character*2       cnest           !domain number (including MOAD)
      integer len,lf 
      integer avgelem
      integer parent_id,ratio_2_parent
      integer uri,urj,lli,llj
      integer indxt,indxl,indxlat,indxlon
      integer istat

      call get_directory('static',dir_out,len) 
      kmax = ngrids

      call find_domain_name(cdataroot,c10_grid_f,istatus)
      if(istatus.ne.1)then
         print*,'Error returned from find_domain_name'
         return
      endif
      call s_len(c10_grid_f,lf)

      call get_c80_description(origin,istatus)

      call s_len(c10_grid_fname,len_fname)
      if(len_fname.gt.9)then
         print*,'domain filename > 9  characters'
         print*,'shorten this to <= 9 characters',
     +' in either c80_description or simulation_name variables'
         stop
      endif

      laps_dom_file = c10_grid_fname(1:len_fname)

      if(c10_grid_f(1:lf).eq.'wrfsi')then
        IF (c6_maproj .ne. 'rotlat') THEN
         indxt=51  !c-stagger
         indxl=12
         indxlat=7
         indxlon=8
         write(cnest,'(i2.2)')nest
         laps_dom_file=laps_dom_file(1:len_fname)//'.d'//cnest
         avgelem=9
         la1=data(1,1,7)
         la2=data(imax-1,jmax-1,7)
         lo1=data(1,1,8)
         lo2=data(imax-1,jmax-1,8)
        ELSE
         laps_dom_file = c10_grid_fname(1:len_fname)//'.rotlat'
         avgelem=18 ! actually AVC, not AVG.
        ENDIF
      else
        indxt=3  !laps analysis grid
        indxl=12
        indxlat=1
        indxlon=2

c commented JRS 5-10-05
c       la1=data(1,1,1)
c       la2=data(imax,jmax,1)
c       lo1=data(1,1,2)
c       lo2=data(imax,jmax,2)

        la1=-999.0
        la2=-999.0
        lo1=-999.0
        lo2=-999.0
      endif

      print*,'Static domain file name: ',laps_dom_file

!     Calculate deltax_cdf and deltay_cdf
      deltax_cdf = deltax
      deltay_cdf = deltay
      std_lat_cdf = std_lat

      if(c6_maproj .eq. 'plrstr')then

! Per Steve Albers, standardize polar stereo grids to +/-60.0 degrees
! LW this is only for writing Dx and Dy to the static file!
            if(std_lat2 .eq. 90.)then
              call get_grid_spacing_actual(60.0,std_lon,deltax_cdf,
     1                                     istatus)
              std_lat_cdf = +60.
            else
              call get_grid_spacing_actual(-60.0,std_lon,deltax_cdf,
     1                                     istatus)
              std_lat_cdf = -60.
            endif

            if(istatus .ne. 1)then
                write(6,*) ' Error calling get_grid_spacing_actual '
     1          ,'from put_laps_static'
                return
            endif

            deltay_cdf = deltax_cdf
            write(6,*)' Polar Stereographic grid identified'
            write(6,*)' Dx and Dy for writing to static file:'
            write(6,*)' Dx = deltax_cdf, Dy = deltay_cdf ',
     1                 deltax_cdf, deltay_cdf

      endif

      write(6,*) dir_out(1:len),len
	
      write(6,*)' List of vars:'
      do i = 1,kmax
          write(6,11)i,var(i),comment(i)(1:60)
11        format(i4,1x,a3,1x,a60)
      enddo

      write(6,*) 'calling wrt_laps_static, kmax = ',kmax

      call s_len(laps_dom_file,lapslen)

      call wrt_laps_static (dir_out(1:len),laps_dom_file(1:lapslen),
     1                      imax,jmax,
     1                      kmax,deltax_cdf,deltay_cdf,
     1                      std_lon,std_lat_cdf,       
     1                      std_lat2,origin,var,comment,
     1                      data,model,grid_spacing,
     1                      c6_maproj,la1,lo1,la2,lo2,
     1                      center_lat, center_lon,lli,llj,uri,
     1                      urj,parent_id,ratio_2_parent,istatus)
      write(6,*) 'return wrt_laps_static ', istatus
      if(istatus .ne. 1)then
          write (6,*)'ERROR wrt_laps_static: status = ',istatus
      else
          write (6,*)'wrt_laps_static: status = ',istatus


          call rd_laps_static(dir_out(1:len),laps_dom_file(1:lapslen)
     1,imax,jmax,kmax,var,units,comment,data,grid_spacing
     1,istatus)
          if(istatus .eq. 1)then
             write (6,*)'rd_laps_static: status = ',istatus
          else
             write (6,*)'ERROR rd_laps_static: status = ',istatus
          endif

      endif

      if(.false.)then
         call terrain_stats(imax,jmax,data(1,1,indxt),data(1,1,indxl)
     &,data(1,1,indxlat),data(1,1,indxlon),c10_grid_fname(1:lf),istat)
         if(istat.eq.0)then
            print*,'Error: returned from terrain_stats'
            istatus=0
         endif
      endif

      return
      end
c
c --------------------------------------------------------
c
      subroutine terrain_stats(nx,ny,tera,landmask,rlat,rlon,
     &cgname,istat)
      implicit none
      character*(*) cgname
      integer       nx,ny,i,ii,j,jj,ic,icl,icg
      integer       istat,istatus
      integer       imax,jmax
      real          tera(nx,ny)
      real          rlat(nx,ny)
      real          rlon(nx,ny)
      real          landmask(nx,ny)
      real          r_missing_data
      real          maxt,mint,sum,avgt,ave,adev,sdev,var,skew,curt
      real          maxtl,mintl,suml,avgtl
      real          maxgl2,mingl2,sumg2,avggl2
      real          maxluadvt,maxlvadvt,maxluvadvt
      real          minluadvt,minlvadvt,minluvadvt
      real          sumuadvt,sumvadvt,sumuvadvt
      real          avgul,avgvl,avguvl
      real          dxdyp2,dterdx2,dterdy2,gspacing,gspace2
      real          dterdx,dterdy,uxdterdx,vxdterdy,uvXdter_tot
      real          uconst,vconst

      real, allocatable :: dxpdy2(:)
      real, allocatable :: udterdx(:)
      real, allocatable :: vdterdy(:)
      real, allocatable :: uvdter_tot(:)

      istat = 0
      call get_r_missing_data(r_missing_data,istatus)
      if(istatus.eq.0)then
         print*,'Error: get_r_missing_data'
         return
      endif

      imax=nx
      jmax=ny
      if(cgname.eq.'wrfsi')then
         imax=nx-1
         jmax=ny-1
      endif

      maxt=-9999.
      maxtl=-9999.
      mint=+9999.
      mintl=+9999.
      mingl2=+9999.
      maxgl2=-9999.
      maxluadvt=-9999.
      maxlvadvt=-9999.
      maxluvadvt=-9999.
      minluadvt=+9999.
      minlvadvt=+9999.
      minluvadvt=+9999.


      sum=0.0
      suml=0.0
      sumg2=0.0
      sumg2=0.0
      sumuadvt=0.0
      sumvadvt=0.0
      sumuvadvt=0.0

      ic=0
      icl=0
      icg=0
      do j=1,jmax
       do i=1,imax
          if(tera(i,j).lt.r_missing_data)then
             if(tera(i,j).gt.maxt)maxt=tera(i,j)
             if(tera(i,j).lt.mint)mint=tera(i,j)
             sum=sum+tera(i,j)
             ic=ic+1
             if(landmask(i,j) .eq. 1)then
                if(tera(i,j).gt.maxtl)maxtl=tera(i,j)
                if(tera(i,j).lt.mintl)mintl=tera(i,j)
                suml=suml+tera(i,j)
                icl=icl+1
             endif
          endif
       enddo
      enddo
      avgt=sum/float(ic)
      avgtl=suml/float(icl)

      call get_grid_spacing(gspacing,istatus)
      if(istatus .eq. 0)then
         print*,'Error: returned from get_grid_spacing'
         return
      endif
      allocate (dxpdy2(icl))
      allocate (udterdx(icl))
      allocate (vdterdy(icl))
      allocate (uvdter_tot(icl))

      dxpdy2=r_missing_data
      udterdx=r_missing_data
      vdterdy=r_missing_data
      uvdter_tot=r_missing_data

      uconst=20.0  !m/s
      vconst=20.0  !m/s
 
      do j=2,jmax-2    !jmax-2 because wrfsi is staggered
       jj=jj+1
       ii=0
       do i=2,imax-2   !imax-2    ditto 
          ii=ii+1
          call get_grid_spacing_actual(rlat(i,j),rlon(i,j)
     1,gspacing,istat)
          gspace2=(gspacing*gspacing)            !m
          dterdx2=tera(i+1,j)-2*tera(i,j)+tera(i-1,j)
          dterdy2=tera(i,j+1)-2*tera(i,j)+tera(i,j-1)
          dxdyp2=dterdx2/gspace2 + dterdy2/gspace2  !m/m

          dterdx=(tera(i+1,j)-tera(i,j))/gspacing
          dterdy=(tera(i,j+1)-tera(i,j))/gspacing
          uxdterdx=dterdx*uconst
          vxdterdy=dterdy*vconst
          uvXdter_tot=uxdterdx+vxdterdy

c land-only stats
          if(landmask(i,j) .eq. 1)then

c laplacian stats
             if(dxdyp2.gt.maxgl2)maxgl2=dxdyp2
             if(dxdyp2.lt.mingl2)mingl2=dxdyp2
             sumg2=sumg2+abs(dxdyp2)
             icg=icg+1
             if(icg.gt.icl)then
                print*,'icg > icl'
                return
             endif
             dxpdy2(icg)=dxdyp2

c terrain advection stats
             if(uxdterdx.gt.maxluadvt)maxluadvt=uxdterdx
             if(vxdterdy.gt.maxlvadvt)maxlvadvt=vxdterdy
             if(uvXdter_tot.gt.maxluvadvt)maxluvadvt=uvXdter_tot
             if(uxdterdx.lt.minluadvt)minluadvt=uxdterdx
             if(vxdterdy.lt.minlvadvt)minlvadvt=vxdterdy
             if(uvXdter_tot.lt.minluvadvt)minluvadvt=uvXdter_tot

             sumuadvt=sumuadvt+abs(uxdterdx)
             sumvadvt=sumvadvt+abs(vxdterdy)
             sumuvadvt=sumuvadvt+abs(uvXdter_tot)

             udterdx(icg)=uxdterdx
             vdterdy(icg)=vxdterdy
             uvdter_tot(icg)=uvXdter_tot

          endif
       enddo
      enddo

      avggl2=sumg2/float(icg)
      avgul=sumuadvt/float(icg)
      avgvl=sumvadvt/float(icg)
      avguvl=sumuvadvt/float(icg)

      call moment(tera,nx*ny,
     &            ave,adev,sdev,var,skew,curt,
     &            istatus)
      if(istatus.eq.1)then
         print*,'Error: returned from sub moment'
         return
      endif

      print*
      print*,'Terrain Statistics: # = ',ic
      print*,'-------------------'
      print*,'Terrain: Max: ',maxt
      print*,'Terrain: Min: ',mint
      print*,'Terrain: Avg: ',avgt
      print*,'Terrain: Average deviation: ', adev
      print*,'Terrain: Standard deviation: ', sdev
      print*,'Terrain: Variance: ', var
      print*

      print*
      print*,'Terrain Statistics - land only # = ',icl
      print*,'----------------------------------------'
      print*,'Land Only: Max Terrain: ',maxtl
      print*,'Land Only: Min Terrain: ',mintl
      print*,'Land Only: Avg Terrain: ',avgtl
      print*

      call moment(dxpdy2,icg,
     &            ave,adev,sdev,var,skew,curt,
     &            istatus)

      print*
      print*,'Terrain Laplacian (land only) Statistics '
      print*,'Units: m of terrain change per km grid space'
      print*,'---------------------------------------------'
      print*,'Laplacian: Max Terrain: ',maxgl2
      print*,'Laplacian: Min Terrain: ',mingl2
      print*,'Laplacian: Avg Terrain: ',avggl2
      print*,'Laplacian: Avg Deviation: ', adev
      print*,'Laplacian: Std Deviation: ', sdev
      print*,'Laplacian: Variance:  ', var
      print*

      call moment(udterdx,icg,
     &            ave,adev,sdev,var,skew,curt,
     &            istatus)
      print*
      print*,'Terrain Advection (land only) '
      print*,'Units: m**2/s of terrain change '
      print*,'---------------------------------------------'
      print*,'U-comp Terrain Adv: Max Terrain: ',maxluadvt
      print*,'U-comp Terrain Adv: Min Terrain: ',minluadvt
      print*,'U-comp Terrain Adv: Avg Terrain: ',ave
      print*,'U-comp Terrain Adv: Avg Deviation: ', adev
      print*,'U-comp Terrain Adv: Std Deviation: ', sdev
      print*,'U-comp Terrain Adv: Variance:  ', var
      print*

      call moment(vdterdy,icg,
     &            ave,adev,sdev,var,skew,curt,
     &            istatus)
      print*
      print*,'Terrain Advection (land only) '
      print*,'Units: m**2/s of terrain change '
      print*,'---------------------------------------------'
      print*,'V-comp Terrain Adv: Max Terrain: ',maxlvadvt
      print*,'V-comp Terrain Adv: Min Terrain: ',minlvadvt
      print*,'V-comp Terrain Adv: Avg Terrain: ',ave
      print*,'V-comp Terrain Adv: Avg Deviation: ', adev
      print*,'V-comp Terrain Adv: Std Deviation: ', sdev
      print*,'V-comp Terrain Adv: Variance:  ', var
      print*

      call moment(uvdter_tot,icg,
     &            ave,adev,sdev,var,skew,curt,
     &            istatus)
      print*
      print*,'Terrain Advection (land only) '
      print*,'Units: m**2/s of terrain change '
      print*,'---------------------------------------------'
      print*,'Total Terrain Adv: Max Terrain: ',maxluvadvt
      print*,'Total Terrain Adv: Min Terrain: ',minluvadvt
      print*,'Total Terrain Adv: Avg Terrain: ',ave
      print*,'Total Terrain Adv: Avg Deviation: ', adev
      print*,'Total Terrain Adv: Std Deviation: ', sdev
      print*,'Total Terrain Adv: Variance:  ', var
      print*

      deallocate (dxpdy2)
      deallocate (udterdx)
      deallocate (vdterdy)
      deallocate (uvdter_tot)
 
      istat=1

      return
      end
