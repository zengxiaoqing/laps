
subroutine COMScount2tbNrad_sub(path_to_raw_sat,max_files,n_lines_ir,n_pixels_ir,r_missing_data,image_lat_ir,image_lon_ir,ir1_tb_out,i4time_data,istatus)


!****************************************************************************
!*This Fortran code converts from COMS count value (Binary)                 *
!*                                     to temperature and radiance (ASCII). *
!*                                                                          *
!*INPUT                                                                     *
!*nx and ny : Array size of COMS data                                       *
!* FD :  1375 * 1429                                                        * 
!* ENH : 1934 * 1544                                                        *
!*Input units number 10 to 11 : Temperature and radiance conversion tables  *
!*Input units number 20 to 24 : COMS count data                             *
!*Each unit is assigned chanel in ascending order(IR1, IR2, WV, SWIR, VIS). *
!*Example code uses Jan 1, 2016, 02:45 UTC data.                            *
!*Input units number 40 (optional) : latitude and longitude conversion table*
!*                                                                          *
!*OUTPUT                                                                    *
!*chanel_tb.dat : Temperature of each chanel obtained from count value      *
!*chanel_rad.dat : Radiance of each chanel obtained from count value        *
!*Their unit numbers are 30 to 39.                                          *
!*                                                                          *
!*Optional commands                                                         *
!*There are several lines for additional analysis.                          *
!*Line number 42, 43 and 154 to 162 are commented.                          *
!*If users want to know latitude and longitude of each cell,                *
!*                                     relevant to same row of latticepoint.*
!*Each column of latticepoint denotes column number and row number of cell. *
!*And each column of latlon represents latitude and longitude of cell,      *
!*                                     relevant to same row of latticepoint.*
!*                                                                          *
!*Author : Gyeong-Gyeun Ha                                                  *
!*Affiliation : National Meteorological Satellite Center (of KMA)           *
!*Date: Nov.2016                                                            *
!****************************************************************************

integer :: i, j
integer, parameter :: nx=1934, ny=1544, nbyte=2, numofch=5, sizeofctab=1024
integer*2, dimension(nx,ny) :: ir1, ir2, wv, swir, vis
real, dimension(sizeofctab,numofch) :: tb, rad
real, dimension(nx,ny) :: ir1_tb, ir2_tb, wv_tb, swir_tb, vis_tb, ir1_rad, ir2_rad, wv_rad, swir_rad, vis_rad
integer, dimension(nx*ny,2) :: latticepoint
real, dimension(nx*ny,2) :: latlon
real, dimension(nx,ny) :: image_lat_ir,image_lon_ir
character*12 a12time_coms
character*13 a13time
integer cvt_wfo_fname13_i4time

! Inputs
character*200 path_to_raw_sat

! Outputs
integer i4time_data(max_files)
real, dimension(n_pixels_ir,n_lines_ir) :: ir1_tb_out
real, dimension(n_pixels_ir,n_lines_ir) :: ir2_tb_out
real, dimension(n_pixels_ir,n_lines_ir) :: wv_tb_out
real, dimension(n_pixels_ir,n_lines_ir) :: swir_tb_out
real, dimension(n_pixels_ir,n_lines_ir) :: vis_rad_out

write(6,*)' Subroutine COMScounttbNrad ',nx,ny,n_pixels_ir,n_lines_ir

a12time_coms = '201701010000'
a13time = a12time_coms(1:8)//'_'//a12time_coms(9:12)
i4time_data(1)=cvt_wfo_fname13_i4time(a13time)
write(6,*)' Read binary data for ',a12time_coms,' ',a13time,i4time_data(1)

! Read binary data

open(20,file=trim(path_to_raw_sat)//'/coms_le1b_ir1_ch1_cn_'//a12time_coms//'.bin',& 
        form='unformatted', recl=nx*ny*nbyte, access='direct',&
        status='old', convert='big_endian')
open(21,file=trim(path_to_raw_sat)//'/coms_le1b_ir2_ch2_cn_'//a12time_coms//'.bin',& 
        form='unformatted', recl=nx*ny*nbyte, access='direct',&
        status='old', convert='big_endian')
open(22,file=trim(path_to_raw_sat)//'/coms_le1b_wv_ch3_cn_'//a12time_coms//'.bin',& 
        form='unformatted', recl=nx*ny*nbyte, access='direct',&
        status='old', convert='big_endian')
open(23,file=trim(path_to_raw_sat)//'/coms_le1b_swir_ch4_cn_'//a12time_coms//'.bin',& 
        form='unformatted', recl=nx*ny*nbyte, access='direct',&
        status='old', convert='big_endian')
open(24,file=trim(path_to_raw_sat)//'/coms_le1b_vis_ch5_cn_'//a12time_coms//'.bin',& 
        form='unformatted', recl=nx*ny*nbyte, access='direct',&
        status='old', convert='big_endian')

read(20, rec=1) ((ir1(i,j), i=1,nx), j=1,ny)
read(21, rec=1) ((ir2(i,j), i=1,nx), j=1,ny)
read(22, rec=1) ((wv(i,j), i=1,nx), j=1,ny)
read(23, rec=1) ((swir(i,j), i=1,nx), j=1,ny)
read(24, rec=1) ((vis(i,j), i=1,nx), j=1,ny)
close(20)
close(21)
close(22)
close(23)
close(24)


write(6,*)' Read temperature and radiance table'

! Read temperature and radiance table

open(10,file=trim(path_to_raw_sat)//'/coms_mi_conversion_table_TB.dat')
open(11,file=trim(path_to_raw_sat)//'/coms_mi_conversion_table_radiance.dat')

do i=1, sizeofctab
	read(10,*) tb(i,1), tb(i,2), tb(i,3), tb(i,4), tb(i,5)
	read(11,*) rad(i,1), rad(i,2), rad(i,3), rad(i,4), rad(i,5)
enddo 

close(10)
close(11)

write(6,*)' Convert count value to TB and Albedo'

! Convert count value to TB and Albedo

do i=1,nx
	do j=1,ny
		ir1_tb(i,j)=tb(ir1(i,j)+1,1)      ! Add 1 to assign array index start from 1 in Fortran 
		ir1_rad(i,j)=rad(ir1(i,j)+1,1)

		ir2_tb(i,j)=tb(ir2(i,j)+1,2)
		ir2_rad(i,j)=rad(ir2(i,j)+1,2)

		wv_tb(i,j)=tb(wv(i,j)+1,3)
		wv_rad(i,j)=rad(wv(i,j)+1,3)

		swir_tb(i,j)=tb(swir(i,j)+1,4)
		swir_rad(i,j)=rad(swir(i,j)+1,4)

		vis_tb(i,j)=tb(vis(i,j)+1,5)
		vis_rad(i,j)=rad(vis(i,j)+1,5)
	enddo
enddo

iwrite = 0 

if(iwrite .eq. 1)then

! Save output in ASCII (Example) 

    open(30,file='./ir1_tb.dat')
    open(31,file='./ir1_rad.dat')
    open(32,file='./ir2_tb.dat')
    open(33,file='./ir2_rad.dat')
    open(34,file='./wv_tb.dat')
    open(35,file='./wv_rad.dat')
    open(36,file='./swir_tb.dat')
    open(37,file='./swir_rad.dat')
    open(38,file='./vis_tb.dat')
    open(39,file='./vis_rad.dat')

    do i=1,nx
	write(30, '(1543(F8.4,x),F8.4)') ir1_tb(i,:)    ! To deal with FD or LA, put 'ny-1' value in 1543.
	write(31, '(1543(ES11.5,x),ES11.5)') ir1_rad(i,:)

	write(32, '(1543(F8.4,x),F8.4)') ir2_tb(i,:)
	write(33, '(1543(ES11.5,x),ES11.5)') ir2_rad(i,:)

	write(34, '(1543(F8.4,x),F8.4)') wv_tb(i,:)
	write(35, '(1543(ES11.5,x),ES11.5)') wv_rad(i,:)

	write(36, '(1543(F8.4,x),F8.4)') swir_tb(i,:)
	write(37, '(1543(ES11.5,x),ES11.5)') swir_rad(i,:)

	write(38, '(1543(F8.4,x),F8.4)') vis_tb(i,:)
	write(39, '(1543(ES11.5,x),ES11.5)') vis_rad(i,:)
    enddo

    close(30)
    close(31)
    close(32)
    close(33)
    close(34)
    close(35)
    close(36)
    close(37)
    close(38)
    close(39)

else
    ir1_tb_out(:,:) = ir1_tb(:,:) 
    ir2_tb_out(:,:) = ir2_tb(:,:) 
    wv_tb_out(:,:) = wv_tb(:,:) 
    swir_tb_out(:,:) = swir_tb(:,:) 
    vis_rad_out(:,:) = vis_rad(:,:) 

endif

! Read latitude and longitude (optional)
! To deal with FD or LA, put proper latitude-longitude table in open command.
! Example code uses latitude-longitude table for ENH mode.  

write(6,*)' Read latlon data'
image_lat_ir = r_missing_data
image_lon_ir = r_missing_data

open(40,file=trim(path_to_raw_sat)//'/cn_latlon.txt',err=998)
do i=1, nx*ny 
  read(40,*)lattice1, lattice2, rlat, rlon
  ii = lattice1+1
  jj = lattice2+1
! write(6,*)lattice1,lattice2,ii,jj,rlat,rlon
  image_lat_ir(ii,jj) = rlat
  image_lon_ir(ii,jj) = rlon
enddo 
close(40)

write(6,*)' Center ir value ',ir1_tb_out(nx/2,ny/2)
write(6,*)' Corner ir value ',ir1_tb_out(1,1)
write(6,*)' Center lat/lon ',image_lat_ir(nx/2,ny/2),image_lon_ir(nx/2,ny/2)
write(6,*)' Corner lat/lon ',image_lat_ir(1,1),image_lon_ir(1,1)

istatus = 1
return ! normal return

998 write(6,*)' ERROR reading COMS lat/lon cn_latlon.txt data'
istatus = 0
return ! error return

end subroutine COMScount2tbNrad_sub
