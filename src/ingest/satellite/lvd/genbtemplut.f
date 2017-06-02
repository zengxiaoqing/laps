      subroutine genbtemplut(cstype,ctype,chid,rcal,cnt2btemp,
     &istatus)
c
c
c
      implicit none
c
      character*(*) cstype
      character*(*) ctype
      integer chid,i
      integer istatus
      real  a,b,spc,rad,rcal
      real    cnt2btemp(0:4095)
      real offset,scale

      istatus=1
c
c Meteosat coefficients are for METEOSAT-7 IR2 and WV1
c
      if(cstype.eq.'meteos')then

         if(chid.eq.3)then
            a =9.2477
            b=-2233.4882
            spc=6.0
         elseif(chid.eq.4)then
            a =6.9618
            b=-1255.5465
            spc=5.0
         endif
         do i=0,255
            rad = (i-spc)*rcal
            rad = max(0.0,rad)
            cnt2btemp(i)=b/(alog(rad)-a)
         enddo

      elseif(ctype.eq.'cdf'.and.chid.eq.3)then   !this is the fsl-conus /public switch
                                                 !channel 3 has different enhancement curve.
         do i=0,255
            cnt2btemp(i) = (1348.925 - float(i) ) / 5.1417
         enddo

      elseif(ctype.eq.'gnp')then
         write(6,*)' Generating btemp lut for gnp type'
         offset = 173.15
         scale = 0.03931624
         do i = 0,4095
           cnt2btemp(i) = float(i) * scale + offset
         enddo ! i
         
      elseif(chid.ne.2)then                      !WFO switch; channels 3,4, and 5 have same
                                                 !enhancement curve atm (10-27-99).
         do i=0,180
            cnt2btemp(i)=(660.0-float(i))/2.0
         enddo
         do i=181,255
            cnt2btemp(i)=420.0-float(i)
         enddo

      elseif(chid.eq.2)then                      !channel 2 (3.9u) has different enchancement
                                                 !curve than the other ir channels.
         do i=0,183
            cnt2btemp(i)=(660.4-float(i))/2.0
         enddo
         do i=184,216
            cnt2btemp(i)=421.7-float(i)
         enddo
         do i=217,255
            cnt2btemp(i)=0.0
         enddo
c
c commented out the separate wv btemp calc. Now the same as ch 4 and 5.
c 10-21-99.
c        elseif(chid.eq.3)then
Count range 255 --> 0, T = (1349.27 - C)/5.141.   4/26/96. Recommendation from D. Birkenheuer
Channel 3 
Count range 255 --> 0, T = (1354.235 - C)/5.1619  5/14/96.   "
c           do i=0,255
C           cnt2btemp(i) = 249.346 - 0.12945*float(i)
C           cnt2btemp(i) = (1354.235 - float(i))/5.1619
C           cnt2btemp(i) = (1349.27 - float(i))/5.141
C           cnt2btemp(i) = (1344.38 - float(i))/5.12 
c
c new as of 10-3-96.
c              cnt2btemp(i) = (1348.925 - float(i) ) / 5.1417
c           enddo
c
      else

         write(6,*)'Channel # error, it doesnt exist'
         istatus=0

      endif

c
c test output of cnt2btemp table
c
c     write(6,*)'cnt 2 btemp for channel: ',chid
c     write(6,*)'-----------------------------'
c     do i=0,255,8
c        write(6,33)i,i+7,(cnt2btemp(j),j=i,i+7)
c     enddo
33    format(2x,i3,'-',i3,1x,'||',2x,8(f5.1,1x))
c     write(6,*)'-----------------------------'

      return
      end
