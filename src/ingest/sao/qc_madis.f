       
         subroutine madis_qc_r(var,DD,level_qc,badflag)

         real var
         character*1 DD
         real badflag

         if(level_qc .ge. 1)then

             if(DD .eq. 'X')then ! Failed Level 1 QC
                 var = badflag
             endif

             if(level_qc .ge. 2)then
                 if(DD .eq. 'Q')then ! Failed Level 2 or 3 QC
                     var = badflag
                 endif
             endif             

         endif

         if(DD .eq. 'B')then ! Subjective QC (reject list)
             var = badflag
         endif
 
         return
         end
                 

         subroutine madis_qc_b(var,iqc,ibmask,idebug,badflag)

         integer ibmask(8)

         logical l_bit(8), btest

         if(iqc .eq. -99)then
             return
         endif

         do i = 1,8
             l_bit(i) = btest(iqc,i)
         enddo

         do i = 1,8
             if(l_bit(i) .AND. ibmask(i) .eq. 1)then
                 if(idebug .eq. 1)then
                     write(6,*)' Failed test ',iqc,i
                 endif
                 var = badflag
             endif
         enddo

         return
         end          
