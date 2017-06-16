
        subroutine region_grow(array,ni,nj,r_missing_data,nbor,npass)

        real array(ni,nj),array_buff(ni,nj)

        do ipass = 1,npass
           array_buff(:,:) = array(:,:)

           nmiss = 0
           nconv = 0

           do i = 1,ni
              iil = max(i-1,1)
              iih = min(i+1,ni)                 
              do j = 1,nj
                 jjl = max(j-1,1)
                 jjh = min(j+1,nj)                 

                 if(array(i,j) .eq. r_missing_data)then ! input is missing

                    nmiss = nmiss + 1

                    icount = 0
                    sum = 0.

                    do ii = iil,iih
                    do jj = jjl,jjh
                       if(array(ii,jj) .ne. r_missing_data)then
                          icount = icount + 1
                          sum = sum + array(ii,jj)
                       endif
                    enddo ! jj
                    enddo ! ii

                    if(icount .ge. nbor)then ! filled in with data
                       array_buff(i,j) = sum / float(icount)
                       nconv = nconv + 1
                    else
                       array_buff(i,j) = r_missing_data
                    endif

                 else
                    continue

                 endif

              enddo ! j
           enddo ! i            

           array(:,:) = array_buff(:,:)

           write(6,*)'ipass/nmiss/nconv',ipass,nmiss,nconv

        enddo ! ipass

        return
        end
