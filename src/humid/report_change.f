      subroutine report_change (data_in, data, plevel,mdf,ii,jj,kk)

      implicit none

      integer ii,jj,kk
      real data_in (ii,jj,kk)
      real data (ii,jj,kk)
      real plevel (kk)
      real mdf ! missing data flag



      integer i,j,k
      integer counter
      integer istatus
      real delta_moisture (kk)
      real avg_moisture (kk)
      real diff_data (ii*jj)
      real ave,adev,sdev,var,skew,curt




c     report moisture change
c     this is a generic loop that can be place about anywhere in the
c     module to help track changes in moisture from any stage
c     this block is planned for a future suboutine.
        
        write(6,*)
        write(6,*)
        
        write(6,*) 'Delta moisture stats:'
        write(6,*) 'Avg = Average difference (g/g) Q'
        write(6,*) 'Std Dev = +/- difference (g/g) Q'
        
        do k = 1,kk
           delta_moisture(k) = 0.0
           avg_moisture(k) = 0.0
           counter = 0
           do i = 1,ii
              do j = 1,jj
                 if( data(i,j,k) .gt. 0.0.and.data(i,j,k).ne.mdf) then
                    counter = counter+1
                    diff_data(counter) = (data(i,j,k) - data_in(i,j,k))
                    delta_moisture(k) = 
     1                   diff_data(counter) + delta_moisture(k)
                    avg_moisture(k) = avg_moisture(k) + data_in(i,j,k)
                 endif
              enddo
           enddo
           if(avg_moisture(k).ne.0) then
              delta_moisture(k) = delta_moisture(k)/avg_moisture(k)
              call moment_b (diff_data,counter,ave,adev,sdev,
     1             var,skew,curt,istatus)
              write(6,*) plevel(k), ave, ' +/-', sdev,' g/g Q'  
           endif
           
        enddo
        write(6,*)
        write(6,*)
        
        write(6,*) 'Relative moisture change % each level'
        write(6,*) 'Avg delta moisture/Avg moisture for level*100'
        
        do k = 1,kk
           write(6,*)plevel(k), delta_moisture(k)*100.,'%'
        enddo
        write(6,*)
        write(6,*)
        
        return
        
        end
