      FUNCTION PLANGO(T,K)                                                      
C $ GOES/I-M RADIANCE FROM BRIGHTNESS TEMPERATURE IN CHANNEL 'K' (HMW)          
C * CALL 'PFCGIM' BEFORE USING THIS OR 'BRITGO' OR 'DBDTGO'                     
      COMMON/PLNCGO/WNUM(25),FK1(25),FK2(25),TC(2,25)                           
      TT=TC(1,K)+TC(2,K)*T                                                      
      EXPN=EXP(FK2(K)/TT) - 1.                                                  
      PLANGO=FK1(K)/EXPN                                                        
      RETURN                                                                    
      END                                                                       
