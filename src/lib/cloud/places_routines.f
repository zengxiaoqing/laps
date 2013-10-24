
      function v_to_b(v)

!     Convert visual magnitude per square second (v) to brightness in 
!     nanolamberts (b)

      v_to_b = 34.08 * exp(20.7233 - 0.92104 * v)

      return
      end


      function b_to_v(b)

!     Convert brightness in nanolamberts (b) to visual magnitude per square 
!     second (v) 

      b_to_v = (log(b/34.08) - 20.7233) / (- 0.92104) 

      return
      end     

      function b_to_maglim(b)

!     Convert brightness in nanolamberts (b) to visual magnitude limit.
!     This starts with Weavers formula (1947). Also Eq. 20 of Garstang (1986).
!     Then, a Schaaf - Albers correction of +0.23 magnitude is applied.

      if(b .le. 1479.)then
          b_to_maglim = 7.930 - 2.171 * log(1. + .1122   * sqrt(b))       
      else ! b > 1479.
          b_to_maglim = 4.305 - 2.171 * log(1. + .001122 * sqrt(b))
      endif

!     b_to_maglim = b_to_maglim + 0.23

!     write(13,1)b,b_to_maglim                               
!     format('  b/b_to_maglim',2f9.2)

      return
      end     


      function b_to_maglim_hecht(b)

!     Convert brightness in nanolamberts (b) to visual magnitude limit.
!     This starts with Hecht's formula (1947) as quoted in Schaefer using
!     Schaefer's Equations. 


      real*4 k

      blog = log10(b)

      if(blog .ge. 3.17)then
          c = 10. ** (-8.35)
          k = 10. ** (-5.90)
      else ! blog < 3.17
          c = 10. ** (-9.80)
          k = 10. ** (-1.90)
      endif

      rith = c * (1. + sqrt(k * b)) ** 2       ! Eq 34
      rmag_v = -16.57 - 2.5 * alog10(rith)     ! Eq 10
      b_to_maglim_hecht = rmag_v - 0.3         ! Eq 9,48

!     write(22,*)b,c,k,rith,rmag_v,b_to_maglim_hecht

      return
      end     

      function rmaglim_to_b(rmaglim)

      bestval = 0.
      resid_min = 999.
!     do i = 1000,10000,2           
      do i = 1000,10000,5           
          testb = 10. ** (float(i) / 1000.)
          testmag = b_to_maglim(testb)
          resid = testmag - rmaglim
          if(abs(resid) .lt. resid_min)then
              resid_min = resid
              rmaglim_to_b = testb
              ibest = i
          endif
!         write(6,1)rmaglim,ibest,resid_min,rmaglim_to_b        
      enddo ! i

!     write(13,1)rmaglim,ibest,resid_min,rmaglim_to_b        
1     format('  rmaglim/ibest/resid_min/rmaglim_to_b ',f9.2,i6,2f9.2)
          
      return
      end


      subroutine addmaglims(rmaglim1,rmaglim2,rmaglim_out)

      b1 = rmaglim_to_b(rmaglim1)
      b2 = rmaglim_to_b(rmaglim2)

      bsum = b1 + b2

      rmaglim_out = b_to_maglim(bsum)

      return
      end
