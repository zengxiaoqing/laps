
      function l_closest_time_i(wmoid,a9time_ob,nobs
     1                        ,i,i4time_sys,istatus)      

cdoc  Determine if the ob time is the closest time for that station to systime

      logical l_closest_time_i

      character*9 a9time_ob(nobs)
      integer wmoid(nobs)
!     character*(*)wmoid(nobs)    ! Allows arbitrary variable type to compare

      i4_closest = 99999

      do j = 1,nobs
          if(wmoid(j) .eq. wmoid(i))then
!             Calculate time of station j
              call i4time_fname_lp(a9time_ob(j),i4time_j,istatus)
              i4_diff = abs(i4time_j - i4time_sys)
              if(i4_diff .lt. i4_closest)then
                  j_closest = j
                  i4_closest = i4_diff
              endif
          endif
      enddo ! j

      if(i .eq. j_closest)then
          l_closest_time_i = .true.
          write(6,*)' Closest time: ',a9time_ob(i)
     1             ,i,wmoid(i),j_closest,i4_closest
      else
          l_closest_time_i = .false.
      endif

      return
      end


      subroutine convert_array(array_in,array_out,n,string
     1                        ,r_missing_data,istatus)       

cdoc  QC the observation array and convert units if needed
cdoc  If 'string' is 'none', then do just the QC without conversion

      character*(*) string

      real k_to_c

      real array_in(n),array_out(n)
      real array_buf(n)

      do i = 1,n
          if(abs(array_in(i)) .ge. 1e10 .or. 
     1           array_in(i)  .eq. r_missing_data )then
              array_out(i) = r_missing_data
          elseif(string .eq. 'k_to_c')then
              array_out(i) = k_to_c(array_in(i))
          elseif(string .eq. 'pa_to_mb')then
              array_out(i) = array_in(i) / 100.
          elseif(string .eq. 'reverse')then
              array_out(i) = array_in(i)
          elseif(string .eq. 'none')then
              array_out(i) = array_in(i)
          elseif(string .ne. 'reverse')then
              write(6,*)' Unknown operator in convert_array: ',string
              istatus = 0
              return
          endif
      enddo ! i
     
      if(string .eq. 'reverse')then ! Reverse output array 
          do i = 1,n
              array_buf(i) = array_out(n+1-i)
          enddo
          do i = 1,n
              array_out(i) = array_buf(i)
          enddo
      endif ! i

      istatus = 1

      return
      end

      subroutine convert_array_i2r(array_in,array_out,n,string
     1                        ,r_missing_data,istatus)       

cdoc  QC the observation array and convert units if needed
cdoc  If 'string' is 'none', then do just the QC without conversion

      character*(*) string

      real k_to_c

      integer array_in(n)
      real array_out(n)
      real array_buf(n)

      do i = 1,n
          if(abs(array_in(i)) .ge. 1000000 .or. 
     1           abs(array_in(i)) .eq. 9999 )then
              array_out(i) = r_missing_data
          elseif(string .eq. 'k_to_c')then
              array_out(i) = k_to_c(array_in(i))
          elseif(string .eq. 'pa_to_mb')then
              array_out(i) = array_in(i) / 100.
          elseif(string .eq. 'reverse')then
              array_out(i) = array_in(i)
          elseif(string .eq. 'none')then
              array_out(i) = array_in(i)
          elseif(string .ne. 'reverse')then
              write(6,*)' Unknown operator in convert_array: ',string
              istatus = 0
              return
          endif
      enddo ! i

      if(string .eq. 'reverse')then ! Reverse output array 
          do i = 1,n
              array_buf(i) = array_out(n+1-i)
          enddo
          do i = 1,n
              array_out(i) = array_buf(i)
          enddo
      endif ! i
     
      istatus = 1

      return
      end

      subroutine apply_qc_rsa(iflag_rsa,variable,nobs)

      integer iflag_rsa(nobs)
      real variable(nobs)

      call get_r_missing_data(r_missing_data,istatus)

      do iob = 1,nobs
          if(iqc_rsa(iflag_rsa(iob)) .eq. -1)then
              variable(iob) = r_missing_data
          endif
      enddo ! iob

      return
      end
