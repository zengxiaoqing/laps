
      subroutine build_comment(cref_comment,cvel_comment,
     &cnyq_comment,nxv01,nyv01,nzv01,max_fields,comment_a,
     &out_array_4d,istatus)

      implicit none

      integer i,j,k,l
      integer nxv01,nyv01,nzv01
      integer max_fields
      integer istatus

      real    out_array_4d(nxv01,nyv01,nzv01,max_fields)
      real    r_missing_data

      character cref_comment*126
      character cvel_comment*126
      character cnyq_comment*126
      character comment_a(max_fields)*125

      integer   ire
      integer   ive
      integer   ine
      integer   icnt(max_fields)

      call get_r_missing_data(r_missing_data,istatus)

      do i=1,126
         if(cref_comment(i:i).eq.'K')then
            ire=i-1
         endif
         if(cvel_comment(i:i).eq.'K')then
            ive=i-1
         endif
         if(cnyq_comment(i:i).eq.'K')then
            ine=i-1
         endif
      enddo

      do l=1,max_fields
         icnt(l)=0
         do k=1,nzv01
         do j=1,nyv01
         do i=1,nxv01
            if(out_array_4d(i,j,k,l).ne.r_missing_data)then
               icnt(l)=icnt(l)+1
            endif
         enddo
         enddo
         enddo
      enddo

      write(cref_comment(ire-5:ire),100)icnt(1) 
      write(cvel_comment(ive-5:ive),100)icnt(2)
      write(cnyq_comment(ine-5:ine),100)icnt(3)
100   format(i6)

      comment_a(1)=cref_comment
      comment_a(2)=cvel_comment
      comment_a(3)=cnyq_comment

      return
      end
