
        subroutine contingency_table(o,f,ni,nj,nk,thresh
     1                              ,lun_out
     1                              ,ilow,ihigh,jlow,jhigh
     1                              ,contable)       

        real o(ni,nj,nk)
        real f(ni,nj,nk)

!       First index is observed, second index is forecast
!       0 is Yes, 1 is No
        integer contable(0:1,0:1)

        contable = 0 ! initialize

        do k = 1,nk
        do i = ilow,ihigh
        do j = jlow,jhigh
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

            contable(index1,index2) = contable(index1,index2) + 1

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
