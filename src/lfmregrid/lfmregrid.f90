program lfmregrid

use mem_namelist

character*256 fname_in, fullname_in, lfmprd_dir, laps_data_root
character*150 static_dir,filename

character*16 atime,afcst
character*9 a9time
character*30  projname
character*2   gproj
character*1 cgrddef 

integer fcsttime,fcsthr,fcstmn

real La1, Lo1, La1in, La2in, LoV, La2, Lo2
real Lat0, Lat1, Lon0
real sw(2),ne(2)

!
! *** Common block variables for polar stereographic grid.
!
integer nx_ps,ny_ps,nz_ps    !No. of PS domain grid points
real lat0_ps,lon0_ps,rota &  !Pol ste. std lat, lon and rotation
            ,sw_ps(2),ne_ps(2)     !SW lat, lon, NE lat, lon
common /psgrid/nx_ps,ny_ps,nz_ps,lat0_ps,lon0_ps &
                    ,rota,sw_ps,ne_ps

!
! *** Common block variables for lambert-conformal grid.
!
integer   nx_lc,ny_lc,nz_lc
real    lat1_lc,lat2_lc,lon0_lc,sw_lc(2),ne_lc(2)
common /lcgrid/nx_lc,ny_lc,nz_lc,lat1_lc,lat2_lc &
,lon0_lc,sw_lc,ne_lc

call getarg(1,fname_in)       ! Input file name (without the .f?? extension)
call getarg(2,a9time)         ! Ascii 9 character time for initialization
call getarg(3,afcst)          ! HHMM of the forecast                                      
call getarg(4,lfmprd_dir)     ! Output directory
call getarg(5,laps_data_root) 

call cv_asc_i4time(a9time,laps_reftime)
read(afcst,'(i2)',err=900) fcsthr        
read(afcst,'(2x,i2)',err=900) fcstmn        
fcsttime=fcsthr*3600 + fcstmn*60
laps_valtime=laps_reftime+fcsttime        

write(6,*)' a9time/laps_reftime ',a9time,' ',laps_reftime

fullname_in = trim(fname_in)//'.fsf'

! Obtain input FSF file horizontal dimensions, note that z will be missing
write(6,*)' Get dimensions from input file: ',trim(fullname_in)
call getdims_lapsprd(fullname_in,nxbg,nybg,nzbg,istatus)
if(istatus .ne. 1)then
    goto900
endif

write(6,*)' input nxbg,nybg,nzbg = ',nxbg,nybg,nzbg

if(.false.)then
    write(6,*)' Get nav info from input file... '
    call read_lapsprd_attr(fullname_in, &
           dxbg, dybg, La1, Lo1, La1in, La2in, LoV, &
           projname, La2,Lo2, istatus) 
    if(istatus.ne.1)then
        print*,'error returned: read_lapsprd_attr'
        goto900
    endif
else
    write(6,*)' Assume hard wired nav info from input file'
    dxbg = 3000.
    dybg = 3000.
    La1 = 21.13812
    Lo1 = -122.7195
    La1in = 38.500
    La2in = 38.500
    Lov = 262.500
    projname = 'tangential lambert conformal'
    La2 = 47.84364        
    Lo2 = -60.90137      
endif

if(Lo1.gt.180)Lo1=Lo1-360
if(Lo2.gt.180)Lo2=Lo2-360
if(LoV.gt.180)LoV=LoV-360
nzbg_ht=nzbg
nzbg_tp=nzbg
nzbg_sh=nzbg
nzbg_uv=nzbg
nzbg_ww=nzbg
sw(1)=La1
sw(2)=Lo1
ne(1)=La2
ne(2)=Lo2
Lon0=LoV
Lat0=La1in
cenlat=La1in
cenlon=LoV
Lat1=La2in

! specify whether the standard lat is Southern or Northern boundary
cgrddef = 'N'

call s_len2(projname,l)

if(projname(1:l).eq. 'polar'.or. &
   projname(1:l).eq.'polar stereographic')then
    gproj='PS'
elseif(projname(1:l).eq. 'mercator')then
    gproj='MC'
elseif(projname(1:l).eq.'lambert')then
    gproj='LC'
elseif(projname(1:l).eq.'secant lambert conformal')then
    gproj='LC'
elseif(projname(1:l).eq.'tangential lambert conformal')then
    gproj='LC'
else
    print*,'ERROR: unable to determine gproj setting ',projname
endif

write(6,*)' Lat0/Lat1/Lon0 = ',Lat0,Lat1,Lon0
write(6,*)' sw = ',sw
write(6,*)' ne = ',ne

call init_gridconv_cmn(gproj,nxbg,nybg,nzbg_ht &
     ,dlat,dlon,cenlat,cenlon,Lat0,Lat1,Lon0 &
     ,sw(1),sw(2),ne(1),ne(2),cgrddef,istatus)

write(6,*)
if(gproj .eq. 'PS')then
    write(6,*)' psgrid common block variables'
    write(6,*)' nx_ps,ny_ps,nz_ps ',nx_ps,ny_ps,nz_ps
    write(6,*)' lat0_ps,lon0_ps,rota,sw_ps,ne_ps ',lat0_ps,lon0_ps,rota,sw_ps,ne_ps
elseif(gproj .eq. 'LC')then
    write(6,*)' lcgrid common block variables'
    write(6,*)' nx_lc,ny_lc,nz_lc ',nx_lc,ny_lc,nz_lc
    write(6,*)' lat1_lc,lat2_lc,lon0_lc,sw_lc,ne_lc ',lat1_lc,lat2_lc,lon0_lc,sw_lc,ne_lc
endif

! Read global LAPS parameters into module memory structure
call get_directory('static',static_dir,len_dir)
filename = static_dir(1:len_dir)//'/nest7grid.parms'
call read_namelist_laps('lapsparms',filename)

write(6,*)' input NX_L,NY_L = ',NX_L,NY_L

call lfmregrid_sub(nxbg,nybg,nzbg,fname_in,NX_L,NY_L,gproj &
                  ,lfmprd_dir,laps_data_root,laps_reftime,laps_valtime)

write(6,*)' Returned from lfmregrid_sub'

900 continue

end
