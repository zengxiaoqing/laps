       subroutine dualwind_noz
     1          (u,v,x,y,ht,vr1,vr2,x1,y1,x2,y2,ht1,ht2,ier)     
c***********************************************************************
c description      : To derive Dual-Doppler wind from radar radial
c                    velocity at constant height.
c                                                                       
c I/O parameters   :                                                    
c  I/O/W   name,            type,      description
c    O     u,v              real       the horizontal wind component of x,y coordinate
c                                      in unit of m/s. 
c    I     x,y              real       the (x,y) coordinate at the wind position.
c                                      in unit of km.
c    I     ht               real       the height at the wind position 
c                                      in unit of km.
c    I     vr1,vr2          real       the radial velocity observed by radar1 and
c                                      radar2 in unit of m/s.
c    I     x1,y1            real       the x,y coordinate at radar 1 in unit of km.
c    I     x2,y2            real       the x,y coordinate at radar 2 in unit of km.
c    I     ht1,ht2          real       the height of radar1 and radar2 in unit of km.
c    O     ier              integer    =0, success message.
c                                      =1, cannot derive horizontal wind.
c                                      =2, (ht1 or ht2) > ht.
c
c Date :
c   Mar. 9, 2004 (S.-M. Deng)
c***********************************************************************

      det=(x-x1)*(y-y2)-(x-x2)*(y-y1)
      if( abs(det).lt.2. )then
          ier=1
          return
      endif
      ier=0
      d1=sqrt( (x-x1)**2 + (y-y1)**2 )
      d2=sqrt( (x-x2)**2 + (y-y2)**2 )
      if( (ht1.gt.ht).or.(ht2.gt.ht) )then
          ier=2
          return
      endif
      cos1=d1/sqrt( d1*d1 + (ht-ht1)**2 ) 
      cos2=d2/sqrt( d2*d2 + (ht-ht2)**2 ) 
      vh1=vr1/cos1
      vh2=vr2/cos2
      u=(vh1*d1*(y-y2)-vh2*d2*(y-y1))/det
      v=(vh2*d2*(x-x1)-vh1*d1*(x-x2))/det

      return
      end
