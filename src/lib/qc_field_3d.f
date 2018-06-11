
      subroutine qc_field_2d(var_2d,field_2d,ni,nj,istatus)
      character*(*) var_2d
      real field_2d(ni,nj)

      call qc_field_3d(var_2d,field_2d,ni,nj,1,istatus)

      return
      end

      subroutine qc_field_3d(var_2d,field_3d,ni,nj,nk,istatus)
      implicit none
      integer ni,nj,nk, istatus, i, j, k
      character*(*) var_2d
      real field_3d(ni,nj,nk), lower_bound, upper_bound, arg
      real r_missing_data, arg_max

      if(var_2d .eq. 'U3' .or. var_2d .eq. 'V3'.or.
     1   var_2d .eq.'USF' .or. var_2d .eq. 'VSF')then
          lower_bound = -200.
          upper_bound = +200.
      elseif(var_2d .eq. 'T3'.or.var_2d.eq.'TSF')then
          lower_bound = +173.
          upper_bound = +400.
      elseif(var_2d .eq. 'DSF')then
          lower_bound = +70.
          upper_bound = +350.
      elseif(var_2d .eq.'PSF'.or.var_2d.eq.'SLP')then
          lower_bound = 30000.  ! including Mount Everest
          upper_bound = 110000.
      else
          lower_bound = -1e10
          upper_bound = +1e10
      endif

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus .ne. 1)return

      call check_nan3(field_3d,ni,nj,nk,istatus) ! relatively efficient
      if(istatus .eq. 0)then ! identify the grid point with the Nan
          do k=1,nk
          do j=1,nj
          do i=1,ni
              call check_nan(field_3d(i,j,k),istatus)
              if(istatus .ne. 1)then
                  write(6,*)' QC Error, Nan detected in ',var_2d,' at '       
     1                     ,i,j,k
                  istatus = 0
                  return
              endif
          enddo ! i
          enddo ! j
          enddo ! k
      endif

      do k=1,nk
      do j=1,nj
      do i=1,ni
          arg = field_3d(i,j,k)
          if(arg .gt. upper_bound)then
            if(arg .ne. r_missing_data)then
              write(6,*)' QC Error detected in ',var_2d,' at ',i,j,k
              write(6,*)' Value exceeded upper bound of '
     1                 ,upper_bound,', value = ',field_3d(i,j,k)       
              istatus = 0
              return
            endif
          endif

          if(arg .lt. lower_bound)then
            if(arg .ne. r_missing_data)then
              write(6,*)' QC Error detected in ',var_2d,' at ',i,j,k
              write(6,*)' Value exceeded lower bound of '
     1                 ,lower_bound,', value = ',field_3d(i,j,k)       
              istatus = 0
              return
            endif
          endif

      enddo ! i
      enddo ! j
      enddo ! k

      do k=1,nk
      do j=1,nj
      do i=1,ni
          if(field_3d(i,j,k) .eq. r_missing_data)then
              write(6,*)' QC Warning detected in ',var_2d,' at ',i,j,k       
              write(6,*)' Value equals r_missing_data or '
     1                 ,r_missing_data      
              istatus = -1
              return
          endif

      enddo ! i
      enddo ! j
      enddo ! k

      if(var_2d.eq.'SLP')then ! check max value
          arg_max = maxval(field_3d)
          if(arg_max .lt. 80000.)then
              write(6,*)' Max SLP is less than 80000.',arg_max
              istatus = 0
              return
          endif
      endif

      istatus = 1
      return

      end 
