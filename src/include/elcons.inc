C***********************************************************************        
C**                                                                             
C**   INTEGRAL SYSTEMS, INC.                                                    
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT                     
C**   SYSTEM    : EARTH LOCATION USERS GUIDE                                    
C**   NAME      : ELCONS                                                        
C**   TYPE      : DATA AREA                                                     
C**   SOURCE    : I.ELCONS                                                      
C**                                                                             
C**   VER.    DATA    BY                COMMENT                                 
C**   ----  --------  ----------------  --------------------------------        
C**   A     01/09/89  I. LEVINE         INITIAL CREATION                        
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   DESCRIPTION                                                               
C**   MATHEMATICAL AND EARTH-RELATED CONSTANTS                                  
C**                                                                             
C***********************************************************************        
C***********************************************************************        
C                                                                               
      REAL*8 PI                                                                 
           PARAMETER (PI=3.141592653589793D0)                                   
      REAL*8 DEG                                                                
           PARAMETER (DEG=180.D0/PI)                                            
      REAL*8 RAD                                                                
           PARAMETER (RAD=PI/180.D0)                                            
C                    DEGREES TO RADIANS CONVERSION PI/180                       
      REAL*8 NOMORB                                                             
           PARAMETER (NOMORB=42164.365D0)                                       
C                    NOMINAL RADIAL DISTANCE OF SATELLITE (km)                  
      REAL*8 AE                                                                 
           PARAMETER (AE=6378.137D0)                                            
C                    EARTH EQUATORIAL RADIUS (km)                               
      REAL*8 FER                                                                
           PARAMETER (FER=1.D0-(6356.7533D0/AE))                                
C                    EARTH FLATTENING COEFFICIENT = 1-(BE/AE)                   
      REAL*4 AEBE2                                                              
           PARAMETER (AEBE2=1.D0/(1.D0-FER)**2)                                 
      REAL*4 AEBE3                                                              
           PARAMETER (AEBE3=AEBE2-1.)                                           
      REAL*4 AEBE4                                                              
           PARAMETER (AEBE4=(1.D0-FER)**4-1.)

