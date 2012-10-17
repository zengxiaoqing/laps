
        subroutine calc_contable_3d(
     1             o,f,thresh,ni,nj,nk,                    ! I
     1             lmask_rqc_3d,r_missing_data,            ! I
     1             cont_3d)                                ! O

        character*150 cont_dir,filename
        character*31 ext

        real o(ni,nj,nk)
        real f(ni,nj,nk)
        real cont_3d(ni,nj,nk)

        logical lmask_rqc_3d(ni,nj,nk)

!       Calculate 3-D contingency table
        do k = 1,nk
        do i = 1,ni
        do j = 1,nj
  
          if(lmask_rqc_3d(i,j,k))then

            if(o(i,j,k) .ge. thresh)then
                index1 = 0
            else
                index1 = 1
            endif

            if(f(i,j,k) .ge. thresh)then
                index2 = 0
            else
                index2 = 1
            endif

            rcont = index1*2 + index2

            cont_3d(i,j,k) = rcont

          else
            cont_3d(i,j,k) = r_missing_data
 
          endif

        enddo ! j
        enddo ! i
        enddo ! k        


        return
        end
