
        subroutine calc_contable(i4_initial,i4_valid,
     1             o,f,thresh,contable,ni,nj,nk,
     1             cont_3d)  

        character*150 cont_dir,filename
        character*31 ext

        real o(ni,nj,nk)
        real f(ni,nj,nk)
        real cont_3d(ni,nj,nk)

        integer contable(0:1,0:1)

!       Calculate 3-D contingency table
        do k = 1,nk
        do i = 1,nj
        do j = 1,ni

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

        enddo ! j
        enddo ! i
        enddo ! k        


        return
        end
