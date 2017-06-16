

       subroutine smooth_box_3d(a,ni,nj,nk,kernsize)

!      3D smoother using a box filter  

       real a(ni,nj,nk) 

       real a_2d(ni,nj)

       idelt = kernsize / 2

       do k = 1,nk
           do i = 1,ni   
           do j = 1,nj   
               il = max(i-idelt,1)
               ih = min(i+idelt,ni)
               jl = max(j-idelt,1)
               jh = min(j+idelt,nj)

               isum = 0
               sum = 0.
                
               do ii = il,ih
               do jj = jl,jh
                   isum = isum + 1
                   sum = sum + a(ii,jj,k)
               enddo ! jj
               enddo ! ii

               a_2d(i,j) = sum / float(isum)

           enddo ! j
           enddo ! i

           a(:,:,k) = a_2d(:,:)

       enddo ! k

       return
       end 

       subroutine smooth_box_2d(a,ni,nj,kernsize)

!      2D smoother using a box filter  

       real a(ni,nj) 

       real a_2d(ni,nj)

       idelt = kernsize / 2

       do i = 1,ni   
       do j = 1,nj   
           il = max(i-idelt,1)
           ih = min(i+idelt,ni)
           jl = max(j-idelt,1)
           jh = min(j+idelt,nj)

           isum = 0
           sum = 0.
                
           do ii = il,ih
           do jj = jl,jh
               isum = isum + 1
               sum = sum + a(ii,jj)
           enddo ! jj
           enddo ! ii

           a_2d(i,j) = sum / float(isum)

       enddo ! j
       enddo ! i

       a(:,:) = a_2d(:,:)

       return
       end 
