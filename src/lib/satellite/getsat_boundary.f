      subroutine get_sat_boundary(nx,ny,max_lines,max_elems,
     &r_llij_lut_ri,r_llij_lut_rj,
     &ilinestart,ilineend,ielemstart,ielemend,
     &rlinestart,rlineend,relemstart,relemend,istatus)

      implicit none

      integer nx,ny
      integer max_lines
      integer max_elems
      integer ilinestart
      integer ilineend
      integer ielemstart
      integer ielemend
c     integer*4 i_grid_spacing_2
      integer istatus
c     integer*4 jstatus

c     real*4 grid_spacing
      real*4 r_llij_lut_ri(nx,ny)
      real*4 r_llij_lut_rj(nx,ny)
      real*4 rlinestart
      real*4 rlineend
      real*4 relemstart
      real*4 relemend
      real*4 r_missing_data

      integer i,j
c
c -------------------------------------
c
      istatus=1

      ilinestart = 1000000
      ilineend   = 0
      ielemstart = 1000000
      ielemend   = 0
      rlinestart = 1000000.
      rlineend   = 0.
      relemstart = 1000000.
      relemend   = 0.

      do j = 1,ny
      do i = 1,nx

         rlinestart=min(r_llij_lut_rj(i,j),rlinestart)
         rlineend  =max(r_llij_lut_rj(i,j),rlineend)
         relemstart=min(r_llij_lut_ri(i,j),relemstart)
         relemend  =max(r_llij_lut_ri(i,j),relemend)

         ilinestart=nint(rlinestart)
         ilineend  =nint(rlineend)
         ielemstart=nint(relemstart)
         ielemend  =nint(relemend)

      enddo
      enddo

      if(ilinestart.lt.0)then
         write(6,*)'WARNING: LAPS exceeds nrthrn Sat bndry!'
         print*,'linestart = ',ilinestart
         rlinestart=0.0
         istatus = 0
      endif
      if(ilineend.gt.max_lines)then
         write(6,*)'WARNING: LAPS exceeds sothrn Sat bndry!'
         print*,'lineend = ',ilineend
         istatus = 0
      endif
      if(ielemstart.lt.0)then
         write(6,*)'WARNING: LAPS exceeds westrn Sat bndry!'
         print*,'elemstart = ',ielemstart
         relemstart=0.0
         istatus = 0
      endif
      if(ielemend.gt.max_elems)then
         write(6,*)'WARNING: LAPS exceeds eastrn Sat bndry!'
         print*,'elemend = ',ielemend
         istatus = 0
      endif

      if(istatus.eq.0)then
         print*,'set look-up values outside bndry to missing'
         call get_r_missing_data(r_missing_data,istatus)
         do j=1,ny
         do i=1,nx
            if(r_llij_lut_rj(i,j).le.0.0.or.
     &         r_llij_lut_rj(i,j).ge.max_lines)then
               r_llij_lut_rj(i,j)=r_missing_data
            endif
            if(r_llij_lut_ri(i,j).le.0.0.or.
     &         r_llij_lut_ri(i,j).ge.max_elems)then
               r_llij_lut_ri(i,j)=r_missing_data
            endif
         enddo
         enddo
      endif
      return
      end
