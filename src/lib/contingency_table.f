
        subroutine contingency_table(o,f,ni,nj,nk             ! I
     1                              ,thresh_o,thresh_f        ! I
     1                              ,lun_out                  ! I
     1                              ,ilow,ihigh,jlow,jhigh    ! I
     1                              ,lmask_rqc_3d             ! I
     1                              ,contable)                ! O

        real o(ni,nj,nk)
        real f(ni,nj,nk)

        logical lmask_rqc_3d(ni,nj,nk)

!       First index is observed, second index is forecast
!       0 is Yes, 1 is No
        integer,parameter :: k12 = selected_int_kind(12)
        integer (kind=k12) :: contable(0:1,0:1)

        contable = 0 ! initialize

        do k = 1,nk
        do i = ilow,ihigh
        do j = jlow,jhigh

          if(lmask_rqc_3d(i,j,k))then

            if(o(i,j,k) .ge. thresh_o)then
                index1 = 0
            else
                index1 = 1
            endif

            if(f(i,j,k) .ge. thresh_f)then
                index2 = 0
            else
                index2 = 1
            endif

            contable(index1,index2) = contable(index1,index2) + 1

          endif

        enddo ! j
        enddo ! i
        enddo ! k

!       Write contingency table
        write(lun_out,*)
        write(lun_out,1)
        write(lun_out,2)
        write(lun_out,3)contable(0,0),contable(1,0) 
        write(lun_out,4)contable(0,1),contable(1,1) 

 1      format('                    Obs')
 2      format('                  Y         N     ')
 3      format('  Fcst  Y',  i11,         i11)
 4      format('        N',  i11,         i11)

        return
        end
