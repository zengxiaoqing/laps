      subroutine lprep_ruc2_hybrid(nx,ny,nz,ht,pr,sh,uw,vw,th
     +     ,usfc,vsfc,tsfc,prsfc,shsfc,htsfc)

c    +     ,gproj,lon0_lc,lat1_lc,lat2_lc)
c
      implicit none
      include 'bgdata.inc'
      integer nx,ny,nz,i,j,k
      real ht(nx,ny,nz),pr(nx,ny,nz),sh(nx,ny,nz),uw(nx,ny,nz)
     +     ,vw(nx,ny,nz),th(nx,ny,nz)
      real usfc(nx,ny), vsfc(nx,ny), tsfc(nx,ny), prsfc(nx,ny),
     +     shsfc(nx,ny),htsfc(nx,ny)
      real cp,g,r,cpog,kappa
      parameter (cp=1004.686,g=9.80665,r=287.053,cpog=cp/g,kappa=r/cp)
      real tv, psi(nx,ny),psj(nx,ny),lat(nx,ny),lon(nx,ny)
      real missing
      
      
c
c *** Common block variables for Lambert-conformal grid.
c *** removed JS (4-01)

c *** Convert Pascals to mb.
c *** Compute tv from thetav.
c *** Compute height from msf.
c *** Compute tp (returned in th) from tv.
c *** Compute sh from mr.
c
      do k=1,nz
         do j=1,ny
            do i=1,nx
               missing = max(pr(i,j,k),th(i,j,k),sh(i,j,k))
               if(missing.lt.missingflag) then
                  pr(i,j,k)=pr(i,j,k)*0.01
                  tv=th(i,j,k)*(pr(i,j,k)*0.001)**kappa
                  th(i,j,k)=tv/(1.+0.61*sh(i,j,k))
                  sh(i,j,k)=sh(i,j,k)/(1.+sh(i,j,k))
               endif
            enddo
         enddo
      enddo
c
c Copy first level into surface fields (this is 5m AGL)
c
      do j=1,ny
         do i=1,nx
            prsfc(i,j) = pr(i,j,1)*100.
            usfc(i,j) = uw(i,j,1)
            vsfc(i,j) = vw(i,j,1)
            tsfc(i,j) = th(i,j,1)
            shsfc(i,j) = sh(i,j,1)
            htsfc(i,j) = ht(i,j,1)
         enddo
      enddo

c
c *** Fill Lambert-conformal common block variables.
c *** removed JS (4-01)

c **** No Longer Needed *****
c *** Convert ruc winds from grid north to true north.
c
c      do j=1,ny
c         do i=1,nx
c            psi(i,j)=float(i)
c            psj(i,j)=float(j)
c         enddo
c      enddo
c      call psij_2_latlon(nx*ny,psi,psj,lat,lon)
c
c      call uvgrid_to_uvtrue_a(uw,vw,lon,lon0,nx,ny,nz)
c
      return
      end

