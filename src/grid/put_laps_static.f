
      subroutine put_laps_static(grid_spacing,model,comment,var
     1                          ,data,imax,jmax,mkmax,ngrids
     1                          ,std_lat,std_lat2,std_lon      
     1                          ,c6_maproj,deltax,deltay)
 
      include 'trigd.inc' 

      integer   mkmax                   !This should be max as defined in calling routine.
      integer   ngrids

      INTEGER  	IMAX,		        !I4time of data
     1		JMAX,KMAX,	        !# cols, # rows, # fields
     1		KDIM,i,j,k,		!K dimension of DATA array
     1		ISTATUS
C
      REAL*4	DATA(imax,jmax,mkmax),	!Raw data to be written
     1		grid_spacing,val
C
      CHARACTER*150     DIR_out		!Directory to be written to
      CHARACTER*31	EXT		!File name ext (up to 31 chars)
      CHARACTER*3	VAR(mkmax)	!3 letter ID of each field
      CHARACTER*10	UNITS(mkmax)	!units of each field
      CHARACTER*(*)	COMMENT(mkmax)	!Comments for each field
      character*(*)     model
      character*80      origin          !Run time parameter - c80_description
      character*10      c10_grid_fname  !The name associated with static files,
!                                        namelists (ie., 'nest7grid'), cdl's, etc.
      character*9       laps_dom_file
      character*6       c6_maproj       !Map projection
      character*200     cdataroot       !Dummy Variable used in find_domain_name.
      integer len,lf 
      integer avgelem
      integer zinelem

      call get_directory('static',dir_out,len) 
      kmax = ngrids
      zinelem = ngrids
      kdim = kmax

      call find_domain_name(cdataroot,c10_grid_fname,istatus)
      if(istatus.ne.1)then
         print*,'Error returned from find_domain_name'
         return
      endif
      call s_len(c10_grid_fname,lf)

      call get_c80_description(origin,istatus)

      call s_len(c10_grid_fname,len_fname)
      if(len_fname.gt.9)then
         print*,'domain filename > 9  characters'
         print*,'shorten this to <= 9 characters',
     +' in either c80_description or simulation_name variables'
         stop
      endif
      laps_dom_file = c10_grid_fname(1:len_fname)

! for LAPS, AVG data is elem 3, otherwise 7.
      avgelem=3
      if(c10_grid_fname(1:lf).eq.'wrfsi')avgelem=7

!     Do zin calc (note this is last [kmax] element in data array)
      do i = 1,imax
      do j = 1,jmax
          psa = ztopsa(data(i,j,avgelem)) 
          data(i,j,zinelem) = (20.0 - ((psa - 100.0) * 0.02))
      enddo ! j
      enddo ! i

!     Calculate deltax_cdf and deltay_cdf
      deltax_cdf = deltax
      deltay_cdf = deltay

      if(.false.)then
!     if(c6_maproj .eq. 'plrstr')then
          call get_ps_parms(std_lat,std_lat2,grid_spacing,phi0
     1                     ,grid_spacing_proj_m)

          if(phi0 .lt. 90.)then
              write(6,*)' Calculate Polar Stereo NetCDF parameters on'
     1                 ,' equivalent projection tangent'
              write(6,*)' to pole. Internal LAPS projection is secant'
              factor = 2. / (1. + sind(phi0))
              deltax_cdf = deltax * factor
              deltay_cdf = deltay * factor
              write(6,*)' deltax_cdf, deltay_cdf',deltax_cdf,deltay_cdf   
          endif

      endif

      write(6,*) dir_out(1:len),len
      call wrt_laps_static (dir_out(1:len),laps_dom_file,imax,jmax,
     1                      kmax,deltax_cdf,deltay_cdf,std_lon,std_lat,       
     1                      std_lat2,origin,var,comment,
     1                      data,model,grid_spacing,
     1                      c6_maproj,istatus)
      if(istatus .eq. 1)then
          write (6,*)'wrt_laps_static: status = ',istatus
      else
          write (6,*)'ERROR wrt_laps_static: status = ',istatus
      endif


      call rd_laps_static(dir_out(1:len),laps_dom_file,imax,jmax,kdim,
     1            var,units,comment,data,grid_spacing,istatus)
      if(istatus .eq. 1)then
          write (6,*)'rd_laps_static: status = ',istatus
      else
          write (6,*)'ERROR rd_laps_static: status = ',istatus
      endif

      return
      end

