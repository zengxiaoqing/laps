C***********************************************************************        
C***********************************************************************        
C**   INTEGRAL SYSTEMS, INC.                                                    
C***********************************************************************        
C**   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT                     
C**   SYSTEM    : EARTH LOCATION USERS GUIDE                                    
C**   NAME      : ELCOMM                                                        
C**   TYPE      : DATA AREA                                                     
C**   SOURCE    : I.ELCOMM                                                      
C**                                                                             
C**   VER.    DATA    BY                COMMENT                                 
C**   ----  --------  ----------------  --------------------------------        
C**   A     01/09/89  I. LEVINE         INITIAL CREATION                        
C***********************************************************************        
C**   DESCRIPTION                                                               
C**   INSTRUMENT POSITION AND ATTITUDE VARIABLES AND TRANSFORMATION             
C**   MATRIX                                                                    
C***********************************************************************        
C                                                                               
C     COMMON VARIABLES                                                          
C                                                                               
      REAL*8 XS(3)                                                              
C                      NORMALIZED S/C POSITION IN ECEF COORDINATES              
      REAL*8 BT(3,3)                                                            
C                      ECEF TO INSTRUMENT COORDINATES TRANSFORMATION            
      REAL*8  Q3                                                                
C                      USED IN SUBROUTINE LPOINT                                
      REAL*8 PITCH,ROLL,YAW                                                     
C                          PITCH,ROLL,YAW ANGLES OF INSTRUMENT (RAD)            
      REAL*4 PMA,RMA                                                            
C                          PITCH,ROLL MISALIGNMENTS OF INSTRUMENT (RAD)         
         COMMON /ELCOMM/ XS,BT,Q3,PITCH,ROLL,YAW,PMA,RMA

