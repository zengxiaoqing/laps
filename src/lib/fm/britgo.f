      FUNCTION BRITGO(R,K)                                                      
C $ GOES/I-M BRIGHTNESS TEMPERATURE FROM RADIANCE IN CHANNEL 'K' (HMW)          
C * CALL 'PFCGIM' BEFORE USING THIS OR 'PLANGO' OR 'DBDTGO'                     
C $   This is experimental software, unsupported by SSEC                        
      COMMON/PLNCGO/WNUM(25),FK1(25),FK2(25),TC(2,25)                           
      EXPN=FK1(K)/R+1.                                                          
      TT=FK2(K)/ALOG(EXPN)                                                      
      BRITGO=(TT-TC(1,K))/TC(2,K)                                               
      RETURN                                                                    
      END                                                                       
