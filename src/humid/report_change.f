cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
      subroutine report_change (data_in, data, p_3d,mdf,ii,jj,kk)

      implicit none

      integer ii,jj,kk
      real data_in (ii,jj,kk)
      real data (ii,jj,kk)
      real p_3d (ii,jj,kk)
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

      write(6,*) 'precheck data in report'
      call check_nan3(data_in,ii,jj,kk,istatus)
      if(istatus.ne.1)then 
         write(6,*) 'NaN in var:data_in routine:report_change'
         return
      endif

      call check_nan3(data,ii,jj,kk,istatus)
      if(istatus.ne.1)then 
         write(6,*) 'NaN in var:data routine:report_change'
         return
      endif     

      call check_nan3(p_3d,ii,jj,kk,istatus)
      if(istatus.ne.1)then 
         write(6,*) 'NaN in var:p_3d routine:report_change'
         return
      endif

      call check_nan (mdf, istatus)
      if(istatus.ne.1)then 
         write(6,*) 'NaN in var:mdf routine:report_change'
         return
      endif     
      write(6,*) 'data_in, data, p_3d, mdf check okay'

 



       
        
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
              write(6,*) 'level ',k, 'approx pressure ',
     1             p_3d(1,1,k), ave, ' +/-', sdev,' g/g Q'  
           endif
           
        enddo
        write(6,*)
        write(6,*)
        
        write(6,*) 'Relative moisture change % each level'
        write(6,*) 'Avg delta moisture/Avg moisture for level*100'
        
        do k = 1,kk
           write(6,*) 'level ',k, delta_moisture(k)*100.,'%'
        enddo
        write(6,*)
        write(6,*)
        
        return
        
        end
