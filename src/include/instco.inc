C***********************************************************************        
C***********************************************************************        
C**                                                                             
C**   INTEGRAL SYSTEMS, INC.                                                    
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT                     
C**   SYSTEM    : EARTH LOCATION USERS GUIDE                                    
C**   NAME      : INSTCOMM                                                      
C**   TYPE      : DATA AREA                                                     
C**   SOURCE    : I.INSTCOMM                                                    
C**                                                                             
C**   VER.    DATA    BY                COMMENT                                 
C**   ----  --------  ----------------  --------------------------------        
C**   A     02/16/89  I. LEVINE         INITIAL CREATION                        
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   DESCRIPTION                                                               
C**   COMMON AREA FOR INSTRUMENT-RELATED CONTROL PARAMETERS                     
C**                                                                             
C***********************************************************************        
C***********************************************************************        
C                                                                               
C     VARIABLES                                                                 
C        CONSTANTS NEEDED TO PERFORM TRANSFORMATIONS BETWEEN THE                
C        LATITUDE/LONGITUDE, LINE/PIXEL AND INSTRUMENT CYCLES/INCREMENTS        
C        COORDINATES.                                                           
C                                                                               
      INTEGER INCMAX(2)                                                       
C                       NUMBER OF INCREMENTS PER CYCLE                          
      REAL*4 ELVMAX(2)                                                          
C                       BOUNDS IN ELEVATION (RADIANS)                           
      REAL*4 SCNMAX(2)                                                          
C                       BOUNDS IN SCAN ANGLE (RADIANS)                          
      REAL*4 ELVINC(2)                                                          
C                       CHANGE IN ELEVATION ANGLE PER INCREMENT (RAD)           
      REAL*4 SCNINC(2)                                                          
C                       CHANGE IN SCAN ANGLE PER INCREMENT (RADIANS)            
      REAL*4 ELVLN(2)                                                           
C                       ELEVATION ANGLE PER DETECTOR LINE (RADIANS)             
      REAL*4 SCNPX(2)                                                           
C                       SCAN ANGLE PER PIXEL (RADIANS)                          
                                                                                
      COMMON /INSTCO/ INCMAX,ELVMAX,SCNMAX,                                     
     1   ELVINC,SCNINC,ELVLN,SCNPX

