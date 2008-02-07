      subroutine lprep_eta_conusc(nx,ny,nz,pr,tp,rh
     .         ,tp_sfc,pr_sfc,rh_sfc,istatus)

c    .         ,gproj,lon0_lc,lat1_lc,lat2_lc,istatus)
c
      implicit none

c JS (4-01) removed projection info
c
c     integer nx,ny,nz,nvars,recdim,i,j,k,istatus

      integer nx,ny,nz,i,j,k,istatus
c
      real   tp(nx,ny,nz),
     .       rh(nx,ny,nz),   ! in as rh out as sh
     .       pr(nx,ny,nz),
     .      tp_sfc(nx,ny),
     .      rh_sfc(nx,ny),   ! in as rh out as sh
     .      pr_sfc(nx,ny)

c     real make_ssh
c
c
      include 'bgdata.inc'
      real xe, mrsat
      integer it
c
c____________________________________________________________
c  pr is input as a single verticle column output as 3d array
c  No longer true.
      istatus = 1

c
c convert ua rh to sh
c
      do k=1,nz
         do j=1,ny
            do i=1,nx
               it=tp(i,j,k)*100
               it=min(45000,max(15000,it))
               xe=esat(it)
               mrsat=0.00622*xe/(pr(i,j,k)-xe)
               rh(i,j,k)=rh(i,j,k)*mrsat
               rh(i,j,k)=rh(i,j,k)/(1.+rh(i,j,k)) 

c              rh(i,j,k)=make_ssh(pr(i,j,k),tp(i,j,k)-273.15,
c     +              rh(i,j,k)/100.,0.0)*0.001
            enddo
         enddo
      enddo
c
c convert sfc rh to sh
c
      do j=1,ny
         do i=1,nx
               it=tp_sfc(i,j)*100
               it=min(45000,max(15000,it))
               xe=esat(it)
               mrsat=0.00622*xe/(pr_sfc(i,j)*0.01-xe)
               rh_sfc(i,j)=rh_sfc(i,j)*mrsat
               rh_sfc(i,j)=rh_sfc(i,j)/(1.+rh_sfc(i,j)) 

c              rh_sfc(i,j)=make_ssh(pr_sfc(i,j)/100.,tp_sfc(i,j)-273.15,
c     +              rh_sfc(i,j)/100.,0.0)*0.001
         enddo
      enddo
      

c
c *** Fill Lambert-conformal common block variables.
c *** removed JS (4-01). Now in get_bkgd_model_info.f

      istatus = 0
      return
      end
