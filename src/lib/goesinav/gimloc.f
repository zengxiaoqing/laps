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

C***********************************************************************
C**   NOAA/NESDIS/SOCC/SOFTWARE BRANCH AND ISI                                  
C***********************************************************************        
C**                                                                             
C**   PROJECT     OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT                     
C**   SYSTEM      EARTH LOCATION USERS GUIDE                                    
C**   ROUTINE     GIMLOC                                                        
C**   SOURCE      F.GIMLOC                                                      
C**   LOAD NAME   ANY                                                           
C**   PROGRAMMER  THOMAS I. BABICKI                                             
C**                                                                             
C**   VER.    DATA    BY   COMMENT                                              
C**   ----  --------  --   ---------------------------------------------        
C**    A     5/01/89  TB   INITIAL CREATION(SOCC/ISI)                           
C**
C**    B     2/19/93  JH   ADD ACCESS TO BOTH IMAGER AND SOUNDER SETS
C**
C***********************************************************************        

C
C AWS -- I eliminated the PROGRAM code here...
C

C***********************************************************************        
C**                                                                             
C**   NOAA/NESDIS/SOCC/SOFTWARE BRANCH                                          
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT                     
C**   SYSTEM    : EARTH LOCATION USERS GUIDE                                    
C**   ROUTINE   : TIME50                                                        
C**   SOURCE    : F.TIME50                                                      
C**   LOAD NAME : ANY                                                           
C**   PROGRAMMER: THOMAS I. BABICKI                                             
C**                                                                             
C**   VER.    DATA    BY   COMMENT                                              
C**   ----  --------  ---  ---------------------------------------------        
C**   A      2/17/89  TB   INITIAL CREATION                                     
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   TIME50 ACCEPTS TWO WORDS CONTAINING DATE AND TIME                         
C**   AND RETURNS TIME EXPRESSED AS DOUBLE PRECISION MINUTES FROM               
C**   1950 JAN. 1.0                                                             
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   CALLED BY       : ANY                                                     
C**   COMMONS MODIFIED: NONE                                                    
C**   INPUTS          : NONE                                                    
C**   OUTPUTS         : NONE                                                    
C**   ROUTINES CALLED : NONE                                                    
C**                                                                             
C***********************************************************************        
      FUNCTION TIME50(I)                                                        
C                                                                               
C     CALLING PARAMETERS                                                        
C                                                                               
      INTEGER*4 I(2)                                                            
C                                                                               
C     LOCAL VARIABLES                                                           
C                                                                               
      INTEGER*4 NY                                                              
C                    YEAR                                                       
      INTEGER*4 ND                                                              
C                    DAY OF YEAR                                                
      INTEGER*4 NH                                                              
C                    HOUR                                                       
      INTEGER*4 NM                                                              
C                    MINUTE                                                     
      REAL*8    S                                                               
C                    SECONDS                                                    
      INTEGER J                                                                 
C                                                                               
C     INCLUDE FILES - REC                                                       
C                                                                               
C                                                                               
C***********************************************************************        
C                                                                               
C     CONVERT INPUT YEAR, DAY OF YEAR, HOUR AND MINUTE                          
C     TO INTEGER VALUES.                                                        
C                                                                               
      NY=I(1)/10000                                                          
      IAA=I(1)-(NY*10000)                                                    
      ND=(I(1)-(NY*10000))*.1                                                
      IAB=(IAA-(ND*10))*10                                                      
      NBC=I(2)/10000000.                                                     
      IAC=I(2)-(NBC*10000000)                                                
      NH=IAB+NBC                                                                
      DEF=I(2)-IAC                                                           
      NM=IAC*.00001                                                             
      S=(I(2)-(DEF+(NM*100000)))*.001D0
      PRINT 1000,NY,ND,NH,NM,S                                                  
 1000 FORMAT (1H1,'YEAR =',I4,/,1X,'JDAY =',I3,/,1X,                            
     *            'HOUR =',I2,/,1X,'MIN  =',I2,/,1X,                            
     *            'SEC  =',F6.3)                                                
C                                                                               
C***********************************************************************        
C                                                                               
C     HERE WE CONVERT INTEGER YEAR AND DAY OF YEAR TO NUMBER OF                 
C     DAYS FROM 0 HOUR UT, 1950 JAN. 1.0                                        
C     THIS CONVERTION IS BASED ON AN ALGORITHM BY FLIEGEL AND VAN               
C     FLANDERN, COMM. OF ACM, VOL.11, NO. 10, OCT. 1968 (P.657)                 
C                                                                               
      J=ND+1461*(NY+4799)/4-3*((NY+4899)/100)/4-2465022                         
C                                                                               
C     COMPUTE TIME IN MINUTES FROM JANUARY 1.0, 1950                            
C                                                                               
      TIME50=J*1440.D0+NH*60.D0+NM+S/60.D0                                      
      RETURN                                                                    
      END                                                                       
C***********************************************************************        
C**   NOAA/NESDIS/SOCC/SOFTWARE BRANCH                                          
C***********************************************************************        
C**                                                                             
C**   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT                     
C**   SYSTEM    : EARTH LOCATION USERS GUIDE                                    
C**   ROUTINE   : TIMEX                                                         
C**   SOURCE    : F.TIMEX                                                       
C**   LOAD NAME : ANY                                                           
C**   PROGRAMMER: THOMAS I. BABICKI                                             
C**                                                                             
C**   VER.    DATA    BY   COMMENT                                              
C**   ----  --------  ---  ---------------------------------------------        
C**   A      4/25/89  TB   INITIAL CREATION                                     
C***********************************************************************        
      FUNCTION TIMEX(NY,ND,NH,NM,S)                                             
C     
      REAL*8 TIMEX
C                                                                          
C     CALLING PARAMETERS                                                        
C                                                                               
      INTEGER*4 NY                                                              
C                    YEAR                                                       
      INTEGER*4 ND                                                              
C                    DAY OF YEAR                                                
      INTEGER*4 NH                                                              
C                    HOUR                                                       
      INTEGER*4 NM                                                              
C                    MINUTE                                                     
      REAL*8    S
C                    SECONDS                                                    
      INTEGER J                                                                 
C                                                                               
C***********************************************************************        
C                                                                               
      PRINT 1000,NY,ND,NH,NM,S                                                  
 1000 FORMAT (/,1X,'YEAR =',I4,/,1X,'JDAY =',I3,/,1X,                           
     *            'HOUR =',I2,/,1X,'MIN  =',I2,/,1X,                            
     *            'SEC  =',F6.3)                                                
C                                                                               
C***********************************************************************        
C                                                                               
      J=ND+1461*(NY+4799)/4-3*((NY+4899)/100)/4-2465022                         
C                                                                               
C     COMPUTE ACTUAL TIME IN MINUTES FROM JANUARY 1.0, 1950                     
C                                                                               
      TIMEX=J*1440.D0+NH*60.D0+NM+S/60.D0                                       
      RETURN                                                                    
      END                                                                       
C***********************************************************************        
C**                                                                             
C**   INTEGRAL SYSTEMS, INC.                                                    
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT                     
C**   SYSTEM    : EARTH LOCATION USERS GUIDE                                    
C**   ROUTINE   : SETCONS                                                       
C**   SOURCE    : F.SETCONS                                                     
C**   LOAD NAME : ANY                                                           
C**   PROGRAMMER: IGOR LEVINE                                                   
C**                                                                             
C**   VER.    DATA    BY   COMMENT                                              
C**   ----  --------  ---  ---------------------------------------------        
C**   A     02/16/89  IL   INITIAL CREATION                                     
C**
C**   B     05/19/94  NP   ADDED CALCULATION OF INSTRUMENT ELEVATION AND
C**                        SCAN ANGLE BIASES BASED ON USER INPUT
C***********************************************************************        
C**                                                                             
C**   THIS SUBROUTINE GENERATES CONSTANTS IN COMMON  INSTCOMM                   
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   CALLED BY       : ANY                                                     
C**   COMMONS MODIFIED: INSTCO                                                  
C**   INPUTS          : NONE                                                    
C**   OUTPUTS         : NONE                                                    
C**   ROUTINES CALLED : NONE                                                    
C**                                                                             
C***********************************************************************        
C***********************************************************************        
      SUBROUTINE SETCON(INSTR,NS_NAD_CY,NS_NAD_INC,EW_NAD_CY,EW_NAD_INC)
C                                                                               
C     CALLING PARAMETERS                                                        
C                                                                               
C                                                                               
C     LOCAL VARIABLES

      INTEGER NS_NAD_CY,NS_NAD_INC,EW_NAD_CY,EW_NAD_INC
C                                                                               
C                                                                               
C     INCLUDE FILES                                                             
C                                                                               
      INCLUDE 'elcons.inc'                                                     
        INCLUDE 'instco.inc'                                                          
C***********************************************************************        
      INCMAX(1)=6136                                                            
      INCMAX(2)=2805                                                            
      ELVINC(1)=8.0D-6
      ELVINC(2)=17.5D-6                                                         
      SCNINC(1)=16.D-6                                                          
      SCNINC(2)=35.D-6                                                          
      ELVLN(1)=28.D-6                                                           
      ELVLN(2)=280.D-6                                                          
      SCNPX(1)=16.D-6                                                           
      SCNPX(2)=280.D-6
C     ************************************************************
C     COMMENTED OUT ELEVATION AND SCAN BIAS CONSTANTS SINCE INSTRUMENT
C     EARTH NADIR POSITION IS AVAILABLE IN GVAR DATA AND PERIODICALLY
C     UPDATED
C
C       
C        ELVMAX(1)=0.220896D0                         
C        ELVMAX(2)=0.22089375D0
C        SCNMAX(1)=0.24544D0
C        SCNMAX(2)=0.2454375D0

C     RECOMPUTE ELEVATION AND SCAN BIASES BASED ON USER INPUTS OF                 
C     CYCLES & INCREMENTS OBTAINED FROM GVAR

c     ELVMAX(INSTR) = (NS_NAD_CY*INCMAX(INSTR)+NS_NAD_INC)*ELVINC(INSTR)

      IF(INSTR.EQ.1)THEN                                                                 
        ELVMAX(INSTR)=(NS_NAD_CY*INCMAX(INSTR)+NS_NAD_INC)*ELVINC(INSTR)
      ELSE
        ELVMAX(INSTR)=((9-NS_NAD_CY)*INCMAX(INSTR)-NS_NAD_INC)
     +                *ELVINC(INSTR)
      ENDIF
      
      SCNMAX(INSTR) = (EW_NAD_CY*INCMAX(INSTR)+EW_NAD_INC)*SCNINC(INSTR)
C     ************************************************************
      RETURN                                                                    
      END                                                                       
C***********************************************************************        
C**                                                                             
C**   INTEGRAL SYSTEMS, INC.                                                    
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT                     
C**   SYSTEM    : EARTH LOCATION USERS GUIDE                                    
C**   ROUTINE   : LMODEL
C**   SOURCE    : F.LMODEL                                                      
C**   LOAD NAME : ANY                                                           
C**   PROGRAMMER: IGOR LEVINE                                                   
C**                                                                             
C**   VER.    DATA    BY   COMMENT                                              
C**   ----  --------  ---  ---------------------------------------------        
C**   1     01/09/89  IL   INITIAL CREATION                                     
C**   2     06/02/89  IL   COORDINATE AXES CHANGED ACCORDING TO                 
C**                        FORD'S DEFINITION IN SDAIP, DRL 504-01               
C**   3     08/21/89  IL   CORRECTED ORBIT ANGLE COMPUTATIONS                   
C**   4     03/08/94  SC   S/C COMPENSATION APPLIED UNCONDITIONALLY;
C**                        REFERENCE RADIAL DISTANCE, LATITUDE AND
C**                        ORBIT YAW SET TO ZERO IF IMC DISABLED.
C**   5     03/08/94  SC   ADDED TRAP FOR SLAT=SYAW=0; CORRECTED
C**                        EXPRESSION FOR LAM.
C***********************************************************************        
C**                                                                             
C**   THIS SUBROUTINE COMPUTES THE POSITION OF THE SATELLITE AND THE            
C**   ATTITUDE OF THE IMAGER OR SOUNDER.  THE CALCULATIONS ARE BASED            
C**   ON THE OATS ORBIT AND ATTITUDE MODEL REPRESENTED BY THE O&A               
C**   PARAMETER SET IN GVAR BLOCK 0.                                            
C**        INPUTS:                                                              
C**          TIME, EPOCH TIME, O&A PARAMETER SET, IMC STATUS.                   
C**                                                                             
C**        OUTPUTS:                                                             
C**          THE SPACECRAFT POSITION VECTOR IN EARTH FIXED COORDINATES;         
C**          THE GEOMETRIC ROLL, PITCH, YAW ANGLES AND THE ROLL,                
C**          PITCH MISALIGNMENTS FOR EITHER THE IMAGER OR THE SOUNDER;          
C**          THE EARTH FIXED TO INSTRUMENT FRAME TRANSFORMATION MATRIX;         
C**          GEOGRAPHIC LATITUDE AND LONGITUDE AT SUBSATELLITE POINT.           
C**                                                                             
C**   DESCRIPTION                                                               
C**   LMODEL ACCEPTS AN INPUT DOUBLE PRECISION TIME IN MINUTES FROM             
C**   1950, JAN.1.0 AND AN INPUT SET OF O&A PARAMETERS AND COMPUTES             
C**   POSITION OF THE SATELLITE, THE ATTITUDE ANGLES AND ATTITUDE               
C**   MISALIGNMENTS AND THE INSTRUMENT TO EARTH FIXED COORDINATES               
C**   TRANSFORMATION MATRIX.                                                    
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   CALLED BY       : ANY                                                     
C**   COMMONS MODIFIED: /ELCOMM/ XS,Q3,PITCH,ROLL,YAW,PMA,RMA,BT                
C**   INPUTS          : NONE                                                    
C**   OUTPUTS         : NONE                                                    
C**   ROUTINES CALLED : INST2ER,GATT                                            
C**                                                                             
C***********************************************************************        
C***********************************************************************        
      SUBROUTINE LMODEL(T,TU,REC,IMC,RLAT,RLON)                                 
C                                                                               
C     CALLING ARGUMENTS                                                         
C                                                                               
      REAL*8 T                                                                  
C                                   INPUT TIME FROM 1950, JAN 1.0 (MIN)         
      REAL*8 TU                                                                 
C                                   EPOCH TIME FROM 1950, JAN 1.0 (MIN)         
      REAL*8 REC(336)                                                           
C                                   INPUT O&A PARAMETER SET                     
      INTEGER IMC                                                               
C                                   INPUT IMC STATUS: 0 - ON, 1 - OFF           
      REAL*8 RLAT                                                               
C                                   SUBSATELLITE GEODETIC LATITUDE (RAD)        
      REAL*8 RLON                                                               
C                                   SUBSATELLITE LONGITUDE IN RADIANS           
C                                                                               
C     LOCAL VARIABLES
C                                                                   
      REAL*8  R                                                                 
C                    NORMALIZED SATELLITE DISTANCE (IN UNITS OF KMER9)          
      REAL*8 TS                                                                 
C                    TIME FROM EPOCH IN MINUTES                                 
      REAL*8 B(3,3)                                                             
C                    SPACCRAFT TO EARTH FIXED COORDINATES TRANSFORMATION        
C                    MATRIX                                                     
      REAL*8 TE                                                                 
C                    EXPONENENTIAL TIME DELAY FROM EPOCH (IN MINUTES)           
      REAL*8 PHI                                                                
C                    SUBSATELLITE GEOCENTRIC LATITUDE IN RADIANS                
      REAL*8 DR                                                                 
C                    RADIAL DISTANCE FROM THE NOMINAL (KM)                      
      REAL*8 PSI                                                                
C                    ORBITAL YAW (IN RADIANS)                                   
      REAL*8  LAM                                                               
C                    IMC LONGITUDE (IN RADIANS)                                 
      REAL*8  U                                                                 
C                    ARGUMENT OF LATITUDE (IN RADIANS)                          
      REAL*8 SU,CU                                                              
C                    SIN(U), COS(U)                                             
      REAL*8 SI,CI                                                              
C                    SINE AND COSINE OF THE ORBIT INCLINATION                   
      REAL*8 SLAT                                                               
C                    SINE OF GEOCENTRIC LATITUDE                                
      REAL*8 ASC                                                                
C                    LONGITUDE OF THE ASCENDING NODE IN RADIANS                 
      REAL*8 SA,CA                                                              
C                    SINE AND COSINE OF ASC                                     
      REAL*8 SYAW                                                               
C                    SINE OF THE ORBIT YAW                                      
      REAL*8 WA                                                                 
C                    SOLAR ORBIT ANGLE IN RADIANS                               
      REAL*8 W                                                                  
C                    ORBIT ANGLE IN RADIANS                                     
      REAL*8 SW,CW                                                              
C                    SIN(W),  COS(W)                                            
      REAL*8 S2W,C2W                                                            
C                    SIN(2*W),  COS(2*W)                                        
      REAL*8 SW1,CW1                                                            
C                    SIN(0.927*W),  COS(0.927*W)                                
      REAL*8 SW3,CW3                                                            
C                    SINE AND COSINE OF 1.9268*W                                
      REAL*8 DLAT                                                               
C                    CHANGE IN SINE OF GEOCENTRIC LATITUDE                      
      REAL*8 DYAW                                                               
C                    CHANGE IN SINE OF ORBIT YAW                                
      REAL*8 GATT                                                               
C                    SUBROUTINE FUNCTION                                        
c      REAL*8 A1,A2                                                              
C                    WORK AREAS                                                 
C                                                                               
C     INCLUDE FILES                                                             
C                                                                               
        INCLUDE 'elcons.inc'                                                          
        INCLUDE 'elcomm.inc'                                                          
C                                                                               
C***********************************************************************        
C                                                                               
C     ASSIGN REFERENCE VALUES TO THE SUBSATELLITE LONGITUDE AND                 
C     LATITUDE, THE RADIAL DISTANCE AND THE ORBIT YAW.                          
C                                                                               
      LAM=REC(5)                                                                
      DR=REC(6)                                                                 
      PHI=REC(7)                                                                
      PSI=REC(8)                                                                
C                                                                               
C     ASSIGN REFERENCE VALUES TO THE ATTITUDES AND MISALIGNMENTS                
C                                                                               
      ROLL=REC(9)                                                               
      PITCH=REC(10)                                                             
      YAW=REC(11)                                                               
      RMA=0.0D0
      PMA=0.0D0
C                                                                               
C     IF IMC IS OFF, COMPUTE CHANGES IN THE SATELLITE ORBIT                     
C                                                                               
      IF (IMC.NE.0) THEN                                                        
C
C     SET REFERENCE RADIAL DISTANCE, LATITUDE AND ORBIT YAW TO ZERO
C
      DR=0.0D0
      PHI=0.0D0
      PSI=0.0D0
C
C     COMPUTE TIME SINCE EPOCH (IN MINUTES)                                     
C                                                                               
      TS=T-TU                                                                   
C                                                                               
C     COMPUTES ORBIT ANGLE AND THE RELATED TRIGONOMETRIC FUNCTIONS.             
C     EARTH ROTATIONAL RATE=.729115E-4 (RAD/S)                                  
C                                                                               
      W=0.729115D-4*60.0D0*TS                                                   
      SW=DSIN(W)
      CW=DCOS(W)
      SW1=DSIN(0.927D0*W)
      CW1=DCOS(0.927D0*W)
      S2W=DSIN(2.0D0*W)
      C2W=DCOS(2.0D0*W)
      SW3=DSIN(1.9268D0*W)
      CW3=DCOS(1.9268D0*W)
C                                                                               
C     COMPUTES CHANGE IN THE IMC LONGITUDE FROM THE REFERENCE                   
C                                                                               
      LAM=LAM+REC(18)+(REC(19)+REC(20)*W)*W                                     
     1        +(REC(27)*SW1+REC(28)*CW1+REC(21)*SW+REC(22)*CW
     2        +REC(23)*S2W+REC(24)*C2W + REC(25)*SW3+REC(26)*CW3                
     3        +W*(REC(29)*SW+REC(30)*CW))*2.0D0
C                                                                               
C     COMPUTES CHANGE IN RADIAL DISTANCE FROM THE REFERENCE (KM)                
C                                                                               
      DR=DR + REC(31) + REC(32)*CW+REC(33)*SW                                   
     1        +REC(34)*C2W+REC(35)*S2W + REC(36)*CW3+REC(37)*SW3                
     2        +REC(38)*CW1+REC(39)*SW1 + W*(REC(40)*CW+REC(41)*SW)              
C                                                                               
C     COMPUTES THE SINE OF THE CHANGE IN THE GEOCENTRIC LATITUDE                
C                                                                               
      DLAT=REC(42) + REC(43)*CW+REC(44)*SW                                      
     1        +REC(45)*C2W+REC(46)*S2W                                          
     2        +W*(REC(47)*CW+REC(48)*SW)                                        
     3        +REC(49)*CW1+REC(50)*SW1                                          
C                                                                               
C     COMPUTES GEOCENTRIC LATITUDE BY USING AN EXPANSION FOR ARCSINE            
C                                                                               
      PHI=PHI+DLAT*(1.0D0+DLAT*DLAT/6.0D0)
C                                                                               
C     COMPUTES SINE OF THE CHANGE IN THE ORBIT YAW                              
C                                                                               
      DYAW=REC(51) + REC(52)*SW+REC(53)*CW                                      
     1        +REC(54)*S2W+REC(55)*C2W                                          
     2        +W*(REC(56)*SW+REC(57)*CW)                                        
     3        +REC(58)*SW1+REC(59)*CW1                                          
C                                                                               
C     COMPUTES THE ORBIT YAW BY USING AN EXPANSION FOR ARCSINE.                 
C                                                                               
      PSI=PSI+DYAW*(1.0D0+DYAW*DYAW/6.0D0)
C                                                                               
C     CALCULATION OF CHANGES IN THE SATELLITE ORBIT ENDS HERE                   
C                                                                               
      END IF                                                                    
C                                                                               
C     CONVERSION OF THE IMC LONGITUDE AND ORBIT YAW TO THE SUBSATELLITE         
C     LONGITUDE AND THE ORBIT INCLINATION (REF: GOES-PCC-TM-2473, INPUTS        
C     REQUIRED FOR EARTH LOCATION AND GRIDDING BY SPS,  JUNE 6, 1988)           
C                                                                               
      SLAT=DSIN(PHI)
      SYAW=DSIN(PSI)
      SI=SLAT**2+SYAW**2                                                        
      CI=DSQRT(1.0D0-SI)
      SI=DSQRT(SI)
      IF (SLAT.EQ.0.0D0.AND.SYAW .EQ. 0.0D0) THEN
      U=0.0D0
      ELSE
      U=DATAN2(SLAT,SYAW)
      ENDIF
      SU=DSIN(U)
      CU=DCOS(U)
C                                                                               
C     COMPUTES LONGITUDE OF THE ASCENDING NODE                                  
C                                                                               
      ASC=LAM - U
      SA=DSIN(ASC)
      CA=DCOS(ASC)
C                                                                               
C     COMPUTES THE SUBSATELLITE GEOGRAPHIC LATITUDE                             
C                                                                               
      RLAT=DATAN(AEBE2*DTAN(PHI))
C                                                                               
C     COMPUTES THE SUBSATELLITE LONGITUDE                                       
C                                                                               
      RLON=ASC+DATAN2(CI*SU,CU)

C                                                                               
C     COMPUTES THE SPACECRAFT TO EARTH FIXED COORDINATES TRANSFORMATION         
C     MATRIX:                                                                   
C         (VECTOR IN ECEF COORDINATES) = B * (VECTOR IN S/C COORDINATES)        
C                                                                               
      B(1,2)=-SA*SI                                                             
      B(2,2)= CA*SI                                                             
      B(3,2)=-CI                                                                
      B(1,3)=-CA*CU+SA*SU*CI                                                    
      B(2,3)=-SA*CU-CA*SU*CI                                                    
      B(3,3)=-SLAT                                                              
      B(1,1)=-CA*SU-SA*CU*CI                                                    
      B(2,1)=-SA*SU+CA*CU*CI                                                    
      B(3,1)= CU*SI                                                             
C                                                                               
C     COMPUTES THE NORMALIZED SPACECRAFT POSITION VECTOR IN EARTH FIXED         
C     COORDINATES - XS.                                                         
C                                                                               
      R=(NOMORB+DR)/AE                                                          
      XS(1)=-B(1,3)*R                                                           
      XS(2)=-B(2,3)*R                                                           
      XS(3)=-B(3,3)*R                                                           
C                                                                               
C     PRECOMPUTES Q3 (USED IN LPOINT)                                           
C                                                                               
      Q3=XS(1)**2+XS(2)**2+AEBE2*XS(3)**2-1.0D0
C                                                                               
C     COMPUTES THE ATTITUDES AND MISALIGNMENTS IF IMC IS OFF                    
C                                                                               
      IF (IMC.NE.0) THEN                                                        
C                                                                               
C     COMPUTES THE SOLAR ORBIT ANGLE                                            
C                                                                               
         WA=REC(60)*TS                                                          
C                                                                               
C     COMPUTES THE DIFFERENCE BETWEEN CURRENT TIME, TS, AND THE                 
C     EXPONENTIAL TIME, REC(61). NOTE THAT BOTH TIMES ARE SINCE EPOCH.          
C                                                                               
         TE=TS-REC(61)                                                          
C                                                                               
C     COMPUTES ROLL + ROLL MISALIGNMENT                                         
C                                                                               
         ROLL=ROLL+GATT(62,REC,WA,TE)                                           
C                                                                               
C     COMPUTES PITCH + PITCH MISALIGNMENT                                       
C                                                                               
         PITCH=PITCH+GATT(117,REC,WA,TE)                                        
C                                                                               
C     COMPUTES YAW                                                              
C                                                                               
         YAW=YAW+GATT(172,REC,WA,TE)                                            
C                                                                               
C     COMPUTES ROLL MISALIGNMENT                                                
C                                                                               
         RMA=GATT(227,REC,WA,TE)                                                
C                                                                               
C     COMPUTES PITCH MISALIGNMENT                                               
C                                                                               
         PMA=GATT(282,REC,WA,TE)                                                
C                                                                               
C     APPLY THE SPCECRAFT COMPENSATION
C                                                                               

         ROLL=ROLL+REC(15)
         PITCH=PITCH+REC(16)
         YAW=YAW+REC(17)

      END IF                                                                    
C                                                                               
C     COMPUTES THE INSTRUMENT TO EARTH FIXED COORDINATES TRANSFORMATION         
C     MATRIX - BT                                                               
C                                                                               
      CALL INST2ER(ROLL,PITCH,YAW,B,BT)
      RETURN                                                                    
      END                                                                       
C***********************************************************************        
C***********************************************************************        
C**                                                                             
C**   INTEGRAL SYSTEMS, INC.                                                    
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT                     
C**   SYSTEM    : EARTH LOCATION USERS GUIDE                                    
C**   ROUTINE   : GPOINT                                                        
C**   SOURCE    : F.GPOINT                                                      
C**   LOAD NAME : ANY                                                           
C**   PROGRAMMER: IGOR LEVINE                                                   
C**                                                                             
C**   VER.    DATA    BY   COMMENT                                              
C**   ----  --------  ---  ---------------------------------------------        
C**   A     12/10/87  IL   INITIAL CREATION                                     
C**   A     06/10/88  IL   REPLACED ASIN WITH ATAN TO SAVE TIME                 
C**   A     06/02/89  IL   COORDINATE AXES CHANGED ACCORDING TO                 
C**                        FORD'S DEFINITION IN SDAIP, DRL 504-01               
C**   4     03/08/94  SC   IMPLEMENTED NEW FORMULAE FOR SCAN ANGLE
C**                        CORRECTION DUE TO MISALIGNMENTS
C***********************************************************************        
C**                                                                             
C**   THIS SUBROUTINE CONVERTS GEOGRAPHIC LATITUDE AND LONGITUDE                
C**   TO THE RELATED ELEVATION AND SCAN ANGLES.                                 
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   CALLED BY       : ANY                                                     
C**   COMMONS MODIFIED: NONE                                                    
C**   INPUTS          : NONE                                                    
C**   OUTPUTS         : NONE                                                    
C**   ROUTINES CALLED : NONE                                                    
C**                                                                             
C***********************************************************************        
C***********************************************************************        
      SUBROUTINE GPOINT(RLAT,RLON,ALF,GAM,IERR)                                 
C                                                                               
C     CALLING PARAMETERS                                                        
C                                                                               
      REAL*8   RLAT                                                             
C                             GEOGRAPHIC LATITUDE IN RADIANS (INPUT)            
      REAL*8   RLON                                                             
C                             GEOGRAPHIC LONGITUDE IN RADIANS (INPUT)           
      REAL*8   ALF                                                              
C                             ELEVATION ANGLE IN RADIANS (OUTPUT)               
      REAL*8   GAM                                                              
C                             SCAN ANGLE IN RADIANS (OUTPUT)                    
      INTEGER IERR                                                              
C                             OUTPUT STATUS; 0 - SUCCESSFUL COMPLETION,         
C                             1 - POINT WITH GIVEN LAT/LON IS INVISIBLE         
C                                                                               
C     LOCAL VARIABLES                                                           
C                                                                               
      REAL*8 F(3)                                                               
C                         POINTING VECTOR IN EARTH CENTERED COORDINATES         
      REAL*8 FT(3)                                                              
C                         POINTING VECTOR IN INSTRUMENT COORDINATES             
      REAL*8 U(3)                                                               
C                         COORDINATES OF THE EARTH POINT (KM)                   
      REAL*8 SING,SLAT,W1,W2                                                    
C                                    WORK SPACE                                 
C                                                                               
C     INCLUDE FILES                                                             
C                                                                               
        INCLUDE 'elcons.inc'                                                          
        INCLUDE 'elcomm.inc'                                                          
C***********************************************************************        
C                                                                               
C     COMPUTES SINUS OF GEOGRAPHIC (GEODETIC) LATITUDE                          
C                                                                               
      SING=DSIN(RLAT)
      W1=AEBE4*SING*SING                                                        
C                                                                               
C     SINUS OF THE GEOCENTRIC LATITUDE                                          
C                                                                               
      SLAT=((0.375D0*W1-0.5D0)*W1+1.0D0)*SING/AEBE2
C                                                                               
C     COMPUTES LOCAL EARTH RADIUS AT SPECIFIED POINT                            
C                                                                               
      W2=SLAT*SLAT                                                              
      W1=AEBE3*W2                                                               
      W1=(0.375D0*W1-0.5D0)*W1+1.D0
C                                                                               
C     COMPUTES CARTESIAN COORDINATES OF THE POINT                               
C                                                                               
      U(3)=SLAT*W1                                                              
      W2=W1*DSQRT(1.0D0-W2)
      U(1)=W2*DCOS(RLON)
      U(2)=W2*DSIN(RLON)
C                                                                               
C     POINTING VECTOR FROM SATELLITE TO THE EARTH POINT                         
C                                                                               
      F(1)=U(1)-XS(1)                                                           
      F(2)=U(2)-XS(2)                                                           
      F(3)=U(3)-XS(3)                                                           
      W2=U(1)*SNGL(F(1))+U(2)*SNGL(F(2))+                                       
     1   U(3)*SNGL(F(3))*AEBE2                                                  
C                                                                               
C     VERIFIES VISIBILITY OF THE POINT                                          
C                                                                               
      IF (W2.GT.0.0D0) THEN
C                               INVISIBLE POINT ON THE EARTH                    
                   IERR=1                                                       
                   ALF=99999.0D0
                   GAM=99999.0D0
                   RETURN                                                       
       END IF                                                                   
C                                                                               
C     CONVERTS POINTING VECTOR TO INSTRUMENT COORDINATES                        
C                                                                               
      FT(1)=BT(1,1)*F(1)+BT(2,1)*F(2)+BT(3,1)*F(3)                              
      FT(2)=BT(1,2)*F(1)+BT(2,2)*F(2)+BT(3,2)*F(3)                              
      FT(3)=BT(1,3)*F(1)+BT(2,3)*F(2)+BT(3,3)*F(3)                              
C                                                                               
C     CONVERTS POINTING VECTOR TO SCAN AND ELEVATION ANGLES AND                 
C     CORRECTS FOR THE ROLL AND PITCH MISALIGNMENTS                             
C                                                                               
      GAM=DATAN(FT(1)/SQRT(FT(2)**2+FT(3)**2))
      ALF=-DATAN(FT(2)/FT(3))
      W1=DSIN(ALF)
      W2=DCOS(GAM)
      ALF=ALF+RMA*(1.0D0-DCOS(ALF)/W2)+PMA*W1*(1.0D0/W2+DTAN(GAM))
      GAM=GAM-RMA*W1                                                            
      IERR=0                                                                    
      RETURN                                                                    
      END                                                                       
C***********************************************************************        
C***********************************************************************        
C**                                                                             
C**   INTEGRAL SYSTEMS, INC.                                                    
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT                     
C**   SYSTEM    : EARTH LOCATION USERS GUIDE                                    
C**   ROUTINE   : INST2ER                                                       
C**   SOURCE    : F.INST2ER                                                     
C**   LOAD NAME : ANY                                                           
C**   PROGRAMMER: IGOR LEVINE                                                   
C**                                                                             
C**   VER.    DATA    BY   COMMENT                                              
C**   ----  --------  ---  ---------------------------------------------        
C**   1     08/16/88  IL   INITIAL CREATION                                     
C**   2     11/11/88  IL   TRIGONOMETRIC FUNCTIONS REPLACED WITH                
C**                        SMALL ANGLE APPROXIMATIONS                           
C**   3    06/02/89   IL   COORDINATE AXES CHANGED ACCORDING TO                 
C**                        FORD'S DEFINITION IN SDAIP, DRL 504-01               
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   INST2ER ACCEPTS THE SINGLE PRECISION ROLL, PITCH AND YAW ANGLES           
C**   OF AN INSTRUMENT AND RETURNS THE DOUBLE PRECISION INSTRUMENT TO           
C**   EARTH COORDINATES TRANSFORMATION MATRIX.                                  
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   CALLED BY       : ANY                                                     
C**   COMMONS MODIFIED: NONE                                                    
C**   INPUTS          : NONE                                                    
C**   OUTPUTS         : NONE                                                    
C**   ROUTINES CALLED : NONE                                                    
C**                                                                             
C***********************************************************************        
C***********************************************************************        
      SUBROUTINE INST2ER(R,P,Y,A,AT)
C                                                                               
C     CALLING PARAMETERS                                                        
C                                                                               
      REAL*8 R                                                                  
C                                   ROLL ANGLE IN RADIANS                       
      REAL*8 P                                                                  
C                                   PITCH ANGLE IN RADIANS                      
      REAL*8 Y                                                                  
C                                   YAW ANGLE IN RADIANS                        
      REAL*8 A(3,3)                                                             
C                                   SPACECRAFT TO ECEF COORDINATES              
C                                   TRANSFORMATION MATRIX                       
      REAL*8 AT(3,3)                                                            
C                                   INSTRUMENT TO ECEF COORDINATES              
C                                   TRANSFORMATION MATRIX                       
C                                                                               
C     LOCAL VARIABLES                                                           
C                                                                               
      REAL*8    RPY(3,3)                                                        
C                                   INSTRUMENT TO BODY COORDINATES              
C                                   TRANSFORMATION MATRIX                       
      INTEGER*4 I,J                                                             
C                                   INDICES                                     
C                                                                               
C     INCLUDE FILES                                                             
C                                                                               
C***********************************************************************        
C                                                                               
C     WE COMPUTE INSTRUMENT TO BODY COORDINATES TRANSFORMATION                  
C     MATRIX BY USING A SMALL ANGLE APPROXIMATION OF TRIGONOMETRIC              
C     FUNCTIONS OF THE ROLL, PITCH AND YAW.                                     
C                                                                               
      RPY(1,1)=1.0D0-0.50D0*(P*P+Y*Y)
      RPY(1,2)=-Y                                                               
      RPY(1,3)=P                                                                
      RPY(2,1)=Y+P*R                                                            
      RPY(2,2)=1.0D0-0.50D0*(Y*Y+R*R)
      RPY(2,3)=-R                                                               
      RPY(3,1)=-P+R*Y                                                           
      RPY(3,2)=R+P*Y                                                            
      RPY(3,3)=1.0D0-0.50D0*(P*P+R*R)
C                                                                               
C     MULTIPLICATION OF MATRICES A AND RPY                                      
C                                                                               
      DO 20 I=1,3                                                               
      DO 10 J=1,3                                                               
      AT(I,J)=A(I,1)*RPY(1,J)+A(I,2)*RPY(2,J)+A(I,3)*RPY(3,J)                   
   10 CONTINUE                                                                  
   20 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
C***********************************************************************        
C***********************************************************************        
C**                                                                             
C**   INTEGRAL SYSTEMS, INC.                                                    
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT                     
C**   SYSTEM    : EARTH LOCATION USERS GUIDE                                    
C**   ROUTINE   : LPOINT                                                        
C**   SOURCE    : F.LPOINT                                                      
C**   LOAD NAME : ANY                                                           
C**   PROGRAMMER: IGOR LEVINE                                                   
C**                                                                             
C**   VER.    DATA    BY   COMMENT                                              
C**   ----  --------  ---  ---------------------------------------------        
C**   A     01/09/89  IL   INITIAL CREATION                                     
C**   A     06/02/89  IL   COORDINATE AXES CHANGED ACCORDING TO                 
C**                        FORD'S DEFINITION IN SDAIP, DRL504-01                
C**   3     03/08/94  SC   IMPLEMENTED NEW FORMULAE FOR SCAN ANGLE
C**                        CORRECTIONS DUE TO MISALIGNMENTS
C***********************************************************************        
C**                                                                             
C**   THIS SUBROUTINE CONVERTS THE INSTRUMENT ELEVATION AND SCAN                
C**   ANGLES TO THE RELATED GEOGRAPHIC LATITUDE AND LONGITUDE.                  
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   CALLED BY       : ANY                                                     
C**   COMMONS MODIFIED: NONE                                                    
C**   INPUTS          : NONE                                                    
C**   OUTPUTS         : NONE                                                    
C**   ROUTINES CALLED : NONE                                                    
C**                                                                             
C***********************************************************************        
C***********************************************************************        
      SUBROUTINE LPOINT(ALPHA,ZETA,RLAT,RLON,IERR)                              
C                                                                               
C     CALLING PARAMETERS                                                        
C                                                                               
      REAL*8   ALPHA                                                            
C                             ELEVATION ANGLE (RAD)                             
      REAL*8   ZETA                                                             
C                             SCAN ANGLE (RAD)                                  
      REAL*8   RLAT                                                             
C                             LATITUDE IN RADIANS (OUTPUT)                      
      REAL*8   RLON                                                             
C                             LONGITUDE IN RADIANS (OUTPUT)                     
      INTEGER IERR                                                              
C                             OUTPUT STATUS; 0 - POINT ON THE EARTH             
C                             FOUND, 1 - INSTRUMENT POINTS OFF EARTH            
C                                                                               
C     LOCAL VARIABLES                                                           
C                                                                               
      REAL*8 G1(3)                                                              
C                          POINTING VECTOR IN EARTH CENTERED COORDINATES        
      REAL*8 H                                                                  
C                          SLANT DISTANCE TO THE EARTH POINT (KM)               
      REAL*8 Q1,Q2,D                                                            
C                          WORK SPACE                                           
      REAL*8 G(3)                                                               
C                          POINTING VECTOR IN INSTRUMENT COORDINATES            
      REAL*8 U(3)                                                               
C                          COORDINATES OF THE EARTH POINT (KM)                  
      REAL*8 SA,CA,DA,DZ,D1,CZ                                                  
C                                     WORK SPACE                                
C                                                                               
C     INCLUDE FILES                                                             
C                                                                               
        INCLUDE 'elcons.inc'                                                          
        INCLUDE 'elcomm.inc'                                                          
C***********************************************************************        
      IERR=1                                                                    
C                                                                               
C     COMPUTES TRIGONOMETRIC FUNCTIONS OF THE SCAN AND ELEVATION                
C     ANGLES CORRECTED FOR THE ROLL AND PITCH MISALIGNMENTS                     
C                                                                               
      CA=DCOS(ALPHA)
      SA=DSIN(ALPHA)
      CZ=DCOS(ZETA)
      DA=ALPHA-PMA*SA*(1.0D0/CZ+DTAN(ZETA))-RMA*(1.0D0-CA/CZ)
      DZ=ZETA+RMA*SA                                                            
C                              CORRECTED SCAN ANGLE                             
      CZ=DCOS(DZ)
C                                                                               
C     COMPUTES POINTING VECTOR IN INSTRUMENT COORDINATES                        
C                                                                               
      G(1)=DSIN(DZ)
      G(2)=-CZ*DSIN(DA)
      G(3)=CZ*DCOS(DA)
C                                                                               
C     TRANSFORMS THE POINTING VECTOR TO EARTH FIXED COORDINATES                 
C                                                                               
      G1(1)=BT(1,1)*G(1)+BT(1,2)*G(2)+BT(1,3)*G(3)                              
      G1(2)=BT(2,1)*G(1)+BT(2,2)*G(2)+BT(2,3)*G(3)                              
      G1(3)=BT(3,1)*G(1)+BT(3,2)*G(2)+BT(3,3)*G(3)                              
C                                                                               
C     COMPUTES COEFFICIENTS AND SOLVES A QUADRATIC EQUATION TO                  
C     FIND THE INTERSECT OF THE POINTING VECTOR WITH THE EARTH                  
C     SURFACE                                                                   
C                                                                               
      Q1=G1(1)**2+G1(2)**2+AEBE2*G1(3)**2                                       
      Q2=XS(1)*G1(1)+XS(2)*G1(2)+AEBE2*XS(3)*G1(3)                              
      D=Q2*Q2-Q1*Q3                                                             
      IF (DABS(D).LT.1.D-9) D=0.0D0
C                                                                               
C     IF THE DISCIMINANTE OF THE EQUATION, D, IS NEGATIVE, THE                  
C     INSTRUMENT POINTS OFF THE EARTH                                           
C                                                                               
      IF (D.LT.0.0D0) THEN
         RLAT=999999.0D0
         RLON=999999.0D0
         RETURN                                                                 
      END IF                                                                    
      D=DSQRT(D)                                                                
C                                                                               
C     SLANT DISTANCE FROM THE SATELLITE TO THE EARTH POINT                      
C                                                                               
      H=-(Q2+D)/Q1                                                              
C                                                                               
C     CARTESIAN COORDINATES OF THE EARTH POINT                                  
C                                                                               
      U(1)=XS(1)+H*G1(1)                                                        
      U(2)=XS(2)+H*G1(2)                                                        
      U(3)=XS(3)+H*G1(3)                                                        
C                                                                               
C     SINUS OF GEOCENTRIC LATITUDE                                              
C                                                                               
      D1=U(3)/DSQRT(U(1)**2+U(2)**2+U(3)**2)
C                                                                               
C     GEOGRAPHIC (GEODETIC) COORDINATES OF THE POINT                            
C                                                                               
      RLAT=DATAN(AEBE2*D1/DSQRT(1.0D0-D1*D1))
      RLON=DATAN2(U(2),U(1))
      IERR=0                                                                    
      RETURN                                                                    
      END                                                                       
C***********************************************************************        
C***********************************************************************        
C**                                                                             
C**   INTEGRAL SYSTEMS, INC.                                                    
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT                     
C**   SYSTEM    : EARTH LOCATION USERS GUIDE                                    
C**   ROUTINE   : SNDELOC                                                       
C**   SOURCE    : F.SNDELOC                                                     
C**   LOAD NAME : ANY                                                           
C**   PROGRAMMER: IGOR LEVINE                                                   
C**                                                                             
C**   VER.    DATA    BY   COMMENT                                              
C**   ----  --------  ---  ---------------------------------------------        
C**   A     02/16/89  IL   INITIAL CREATION                                     
C***********************************************************************        
C**                                                                             
C**   SNDELOC ACCEPTS THE MIRROR POSITION IN CYCLES AND INCREMENTS,             
C**   SERVO ERROR VALUES, AND THE POSITIONAL OFFSETS FOR FOUR DETECTORS         
C**   OF A SELECTED SOUNDER CHANNEL AND COMPUTES THE DETECTOR EARTH             
C**   LOCATIONS IN LATITUDE/LONGITUDE COORDINATES.                              
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   CALLED BY       : ANY                                                     
C**   COMMONS MODIFIED: NONE                                                    
C**   INPUTS          : NONE                                                    
C**   OUTPUTS         : NONE                                                    
C**   ROUTINES CALLED : LPOINT                                                  
C**                                                                             
C***********************************************************************        
C***********************************************************************        
      SUBROUTINE SNDELO(CYEW,INCEW,CYNS,INCNS,SVEW,SVNS,DOFF,GEO)               
C                                                                               
C     CALLING PARAMETERS                                                        
C                                                                               
      INTEGER CYEW                                                              
C                         E-W CYCLES                                            
      INTEGER INCEW                                                             
C                         E-W INCREMENTS                                        
      INTEGER CYNS                                                              
C                         N-S CYCLES                                            
      INTEGER INCNS                                                             
C                         N-S INCREMENTS                                        
      REAL*8 SVEW                                                               
C                         E-W SERVO ERROR IN RADIANS                            
      REAL*8 SVNS                                                               
C                         N-S SERVO ERROR IN RADIANS                            
      REAL*8 DOFF(4,2)                                                          
C                         OFFSETS FOR 4 DETECTORS (RADIANS)                     
C                            DOFF(*,1) = E-W OFFSET                             
C                            DOFF(*,2) = N-S OFFSET                             
      REAL*8 GEO(4,2)                                                           
C                         GEOGRAPHIC COORDINATES RELATED TO 4 DETECTORS         
C                            GEO(*,1) = LATITUDE IN RADIANS                     
C                            GEO(*,2) = LONGITUDE IN RADIANS                    
C                                                                               
C     LOCAL VARIABLES                                                           
C                                                                               
c     REAL*8 E,S,H,EV,SC,ALPHA,BETA,SINE,COSE,DE,DS
      REAL*8 E,S,H,EV,SC,SINE,COSE,DE,DS
      INTEGER I,IER                                                             
C                                                                               
C     INCLUDE FILES                                                             
C                                                                               
        INCLUDE 'instco.inc'                                                          
C***********************************************************************        
C                                                                               
C     CONVERT THE MIRROR POSITION, GIVEN IN CYCLES AND INCREMENTS, TO           
C     ELEVATION AND SCAN ANGLES                                                 
C                                                                               
C     E=(CYNS*INCMAX(2)+INCNS)*ELVINC(2)-ELVMAX(2)
      E=((CYNS-9)*INCMAX(2)+INCNS)*ELVINC(2)-ELVMAX(2)                              
      S=(CYEW*INCMAX(2)+INCEW)*SCNINC(2)-SCNMAX(2)                              
C                                                                               
C     CORRECT ELEVATION AND SCAN ANGLES FOR SERVO ERRORS OBTAINING THE          
C     TRUE MIRROR POINTING                                                      
C                                                                               
      E=E+SVNS                                                                  
      S=S+SVEW
      SINE=DSIN(E)
      COSE=DCOS(E)
      H=-2.0D0*SCNPX(2)
C
C     COMPUTE DETECTOR ROTATION OFFSETS FOR EACH DETECTOR
C
C      ALPHA = 0.643501D0 + E
C      BETA = 0.244979D0 - E
C
C      DOFF(1,1) = -0.064976D0
C      DOFF(1,2) = 0.00042D0
C      DOFF(2,1) = 0.00056D0
C      DOFF(2,2) = 0.00014D0 
C      DOFF(3,1) = -0.064976D0
C      DOFF(3,2) = -0.065396D0
C      DOFF(4,1) = 0.00056D0
C      DOFF(4,2) = -0.065116D0
C      DOFF(1,1) = - 700.0D0*DCOS(ALPHA)*1.0D-6
C      DOFF(1,2) =   700.0D0*DSIN(ALPHA)*1.0D-6
C      DOFF(2,1) =   577.23479D0*DCOS(BETA)*1.D-6
C      DOFF(2,2) =   577.23479D0*DSIN(BETA)*1.0D-6
C      DOFF(3,1) = - 577.23479D0*DCOS(BETA)*1.0D-6
C      DOFF(3,2) = - 577.23479D0*DSIN(BETA)*1.0D-6
C      DOFF(4,1) =   700.0D0*DCOS(ALPHA)*1.0D-6
C      DOFF(4,2) = - 700.0D0*DSIN(ALPHA)*1.0D-6
C                                                                               
C     COMPUTE EARTH LOCATIONS FOR FOUR DETECTORS
C                                                                               
      DO 10 I=1,4

C     COMPUTE POSITIONAL OFFSETS OF I-TH DETECTOR

        DE=(2.5-I)*ELVLN(2)+DOFF(I,2)
        DS=H+DOFF(I,1)
C                                                                               
C     COMPUTE ELEVATION AND SCAN ANGLES RELATED TO I-TH DETECTOR                
C     AND CORRECT THEM FOR THE DETECTOR POSITIONAL OFFSETS                      
C                                                                               
C           EV=E+DOFF(I,2)
C           SC=S+DOFF(I,1)
C
C     CONVERT POSITIONAL OFFSETS TO ANGULAR OFFSETS AND
C     CORRECT ELEVATION AND SCAN ANGLES

            EV = E + DE*COSE - DS*SINE
            SC = S + DE*SINE + DS*COSE

C     TRANSFORM DETECTOR'S POINTING ANGLES TO GEOGRAPHIC COORDINATES            
C     OF THE CORRESPONDING POINT ON THE EARTH SURFACE.                          
C     NOTE:  IF A DETECTOR LOOKS OFF THE EARTH, THE RELATED LATITUDE            
C            AND LONGITUDE ARE SET TO 999999.                                   
C                                                                               
           CALL LPOINT(EV,SC,GEO(I,1),GEO(I,2),IER)                             
           H=-H                                                                 
   10 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
C***********************************************************************        
C***********************************************************************        
C***********************************************************************        
C**                                                                             
C**   INTEGRAL SYSTEMS, INC.                                                    
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT                     
C**   SYSTEM    : EARTH LOCATION USERS GUIDE                                    
C**   ROUTINE   : EVSC2LPF                                                      
C**   SOURCE    : F.EVSC2LPF                                                    
C**   LOAD NAME : ANY                                                           
C**   PROGRAMMER: IGOR LEVINE                                                   
C**                                                                             
C**   VER.    DATA    BY   COMMENT                                              
C**   ----  --------  ---  ---------------------------------------------        
C**   A     10/27/88  IL   INITIAL CREATION                                     
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   THIS SUBROUTINE CONVERTS ELEVATION AND SCAN ANGLES                        
C**   TO THE FRACTIONAL LINE AND PIXEL NUMBERS.                                 
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   CALLED BY       : ANY                                                     
C**   COMMONS MODIFIED: NONE                                                    
C**   INPUTS          : NONE                                                    
C**   OUTPUTS         : NONE                                                    
C**   ROUTINES CALLED : NONE                                                    
C**                                                                             
C***********************************************************************        
C***********************************************************************        
      SUBROUTINE EVSC2L(INSTR,ELEV,SCAN,RL,RP)                                  
C                                                                               
C     CALLING PARAMETERS                                                        
C                                                                               
      INTEGER INSTR                                                             
C                       INSTRUMENT CODE (1-IMAGER, 2-SOUNDER)                   
      REAL*8 ELEV                                                                 
C                       ELEVATION ANGLE IN RADIANS                              
      REAL*8 SCAN                                                                 
C                       SCAN ANGLE IN RADIANS                                   
      REAL*8 RL
C                       LINE NUMBER                                             
      REAL*8 RP
C                       PIXEL NUMBER                                            
C                                                                               
C     LOCAL VARIABLES - NONE                                                    
C                                                                               
C                                                                               
C     INCLUDE FILES                                                             
C                                                                               
        INCLUDE 'instco.inc'                                                          
C**************************************************************                 
C                                                                               
C     COMPUTE FRACTIONAL LINE NUMBER                                            
C                                                                               
      RL=(ELVMAX(INSTR)-ELEV)/ELVLN(INSTR)                                      
      IF (INSTR.EQ.1) THEN
           RL=RL+4.5D0
      ELSE
           RL=RL+2.5D0
      END IF
C                                                                               
C     COMPUTE FRACTIONAL PIXEL NUMBER                                           
C                                                                               
      RP=(SCNMAX(INSTR)+SCAN)/SCNPX(INSTR)+1.0D0
      RETURN                                                                    
      END                                                                       
C***********************************************************************        
C**                                                                             
C**   INTEGRAL SYSTEMS, INC.                                                    
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT                     
C**   SYSTEM    : EARTH LOCATION USERS GUIDE                                    
C**   ROUTINE   : EVLN                                                          
C**   SOURCE    : F.EVLN                                                        
C**   LOAD NAME : ANY                                                           
C**   PROGRAMMER: IGOR LEVINE                                                   
C**                                                                             
C**   VER.    DATA    BY   COMMENT                                              
C**   ----  --------  ---  ---------------------------------------------        
C**   A     10/27/88  IL   INITIAL CREATION                                     
C***********************************************************************        
C**                                                                             
C**   THIS FUNCTION CONVERTS FRACTIONAL LINE NUMBER TO ELEVATION ANGLE          
C**   IN RADIANS.                                                               
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   CALLED BY       : ANY                                                     
C**   COMMONS MODIFIED: NONE                                                    
C**   INPUTS          : NONE                                                    
C**   OUTPUTS         : NONE                                                    
C**   ROUTINES CALLED : NONE                                                    
C**                                                                             
C***********************************************************************        
C***********************************************************************        
      FUNCTION EVLN(INSTR,RLINE)                                                
C                                                                               
C     CALLING PARAMETERS                                                        
C                                                                               
      INTEGER INSTR                                                             
C                       INSTRUMENT CODE (1-IMAGER, 2-SOUNDER)                   
      REAL*8  EVLN,RLINE
C                       FRACTIONAL LINE  NUMBER                                 
C                                                                               
C     LOCAL VARIABLES - NONE                                                    
C                                                                               
C                                                                               
C     INCLUDE FILES                                                             
C                                                                               
        INCLUDE 'instco.inc'                                                          
C***********************************************************************        
      IF (INSTR.EQ.1) THEN                                                      
       EVLN=ELVMAX(INSTR)*1.0D0 - (RLINE-4.5 )*(ELVLN(INSTR)*1.0D0)
      ELSE                                                                      
           EVLN=ELVMAX(INSTR)*1.0D0 - (RLINE
     +          -2.5D0)*(ELVLN(INSTR)*1.0D0)
      END IF                                                                    
      RETURN                                                                    
      END                                                                       
C***********************************************************************        
C**                                                                             
C**   INTEGRAL SYSTEMS, INC.                                                    
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT                     
C**   SYSTEM    : EARTH LOCATION USERS GUIDE                                    
C**   ROUTINE   : SCPX                                                          
C**   SOURCE    : F.SCPX                                                        
C**   LOAD NAME : ANY                                                           
C**   PROGRAMMER: IGOR LEVINE                                                   
C**                                                                             
C**   VER.    DATA    BY   COMMENT                                              
C**   ----  --------  ---  ---------------------------------------------        
C**   A     09/22/87  IL   INITIAL CREATION                                     
C***********************************************************************        
C**                                                                             
C**   THIS FUNCTION CONVERTS FRACTIONAL PIXEL NUMBER TO SCAN ANGLE              
C**   IN RADIANS.                                                               
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   CALLED BY       : ANY                                                     
C**   COMMONS MODIFIED: NONE                                                    
C**   INPUTS          : NONE                                                    
C**   OUTPUTS         : NONE                                                    
C**   ROUTINES CALLED : NONE                                                    
C**                                                                             
C***********************************************************************        
C***********************************************************************        
      FUNCTION SCPX(INSTR,PIX)                                                  
C                                                                               
C     CALLING PARAMETERS                                                        
C                                                                               
      INTEGER INSTR                                                             
C                       INSTRUMENT CODE (1-IMAGER, 2-SOUNDER)                   
      REAL*8 SCPX,PIX
C                       FRACTIONAL PIXEL NUMBER                                 
C                                                                               
C     LOCAL VARIABLES                                                           
C                                                                               
C                                                                               
C     INCLUDE FILES                                                             
C                                                                               
        INCLUDE 'instco.inc'                                                    
C***********************************************************************        
      SCPX=((PIX*1.0D0)-1.0D0)*(SCNPX(INSTR)*1.0D0)
     +     -(SCNMAX(INSTR)*1.0D0)
      RETURN                                                                    
      END                                                                       
C***********************************************************************        
C**                                                                             
C**   INTEGRAL SYSTEMS, INC.                                                    
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT                     
C**   SYSTEM    : EARTH LOCATION USERS GUIDE                                    
C**   ROUTINE   : GATT                                                          
C**   SOURCE    : F.GATT                                                        
C**   LOAD NAME : ANY                                                           
C**   PROGRAMMER: IGOR LEVINE                                                   
C**                                                                             
C**   VER.    DATA    BY   COMMENT                                              
C**   ----  --------  ---  ---------------------------------------------        
C**   A     12/01/88  IL   INITIAL CREATION                                     
C**                                                                             
C***********************************************************************        
C**                                                                             
C**    THIS FUNCTION COMPUTES AN ATTITUDE/MISALIGNMENT ANGLE FROM A             
C**    GIVEN SUBSET OF THE O&A PARAMETERS IN GVAR BLOK 0.                       
C**    ARGUMENT K0 INDICATES THE FIRST WORD OF THE SUBSET.                      
C**                                                                             
C***********************************************************************        
C**                                                                             
C**   CALLED BY       : LMODEL                                                  
C**   COMMONS MODIFIED: NONE
C**   INPUTS          : NONE                                                    
C**   OUTPUTS         : NONE                                                    
C**   ROUTINES CALLED : NONE                                                    
C**                                                                             
C***********************************************************************        
C***********************************************************************        
      FUNCTION GATT(K0,REC,WA,TE)
C                                                                               
C     CALLING PARAMETERS                                                        
C                                                                               
      INTEGER K0                                                                
C                    STARTING POSITION OF A PARAMETER SUBSET IN THE             
C                    O&A SET                                                    
      REAL*8 REC(336)                                                           
C                    INPUT O&A PARAMETER SET                                    
      REAL*8 WA
C                    INPUT SOLAR ORBIT ANGLE IN RADIANS                         
      REAL*8 TE
C                    INPUT EXPONENTIAL TIME DELAY FROM EPOCH (MINUTES)          
C                                                                               
C     LOCAL VARIABLES                                                           
C                                                                               
      INTEGER*4 I,J,M,L,LL,K                                                    
      REAL*8 GATT,IR,JR,MR,ATT
      EQUIVALENCE (J,JR),(M,MR)
C                                                                               
C     INCLUDE FILES                                                             
C                                                                               
C***********************************************************************        
C                                                                               
C     CONSTANT COMPONENT                                                        
C                                                                               
      K=K0                                                                      
      ATT=REC(K+2)                                                              
C                                                                               
C     COMPUTES THE EXPONENTIAL TERM                                             
C                                                                               
      IF ((TE.GE.0.0D0).AND.(REC(K+1).GT. 0))THEN
         ATT=ATT+REC(K)*DEXP(-TE/REC(K+1))
      ENDIF
C                                                                               
C     EXTRACTS THE NUMBER OF SINUSOIDS                                          
C                                                                               
      IR=REC(K+3)
      I = INT(IR)
C                                                                               
C     CALCULATION OF SINUSOIDS                                                  
C                                                                               
      DO 10 L=1,I                                                               
           ATT=ATT+REC(K+2*L+2)*DCOS(WA*L+REC(K+2*L+3))
   10 CONTINUE                                                                  
C                                                                               
C     POINTER TO THE NUMBER OF MONOMIAL SINUSOIDS                               
C                                                                               
      K=K+34                                                                    
C                                                                               
C     EXTACTS NUMBER OF MONOMIAL SINUSOIDS                                      
C                                                                               
      IR=REC(K)
      I=INT(IR)
C     KKK=REC(K)
C                                                                               
C     COMPUTES MONOMIAL SINUSOIDS
C                                                                               
      DO 20 L=1,I
           LL=K+5*L
C                                                                               
C          ORDER OF SINUSOID                                                    
C                                                                               
           JR=REC(LL-4)                                                         
C                                                                               
C          ORDER OF MONOMIAL SINUSOID                                           
C                                                                               
           MR=REC(LL-3)                                                         
C                                                                               
           ATT=ATT+REC(LL-2)*((WA-REC(LL))**M)*DCOS(J*WA+REC(LL-1))
   20 CONTINUE                                                                  
      GATT=ATT                                                                  
      RETURN                                                                    
      END                                                                       
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *        
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *        
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *        
C DLAT,DLON...... 50.00000  -150.00000       (15X,F9.5,2X,F10.5)
C LINE...(1),(2). 3584   1230                (15X,I5,2X,I5)
C IPIXEL.(1),(2).10253   1155                (15X,I5,2X,I5)
C KEY............0                           (15X,I1)
C IMC............1                           (15X,I1)
C INSTR..........2                           (15X,I1)
C EPOCH TIME.....1989032062934.567           (15X,I4,I3,2I2,F6.3)
C START TIME.....1989032064934.567           (15X,I4,I3,2I2,F6.3)
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
* *
C      DLAT ---------- LATITUDE  (NEGATIVE IS SOUTH)(F9.5)
C      DLON ---------- LONGITUDE (NEGATIVE IS WEST)(F10.5)

C      LINE(2) ------- LINE NO.  (SEE INSTR FOR LINE(INSTR)(I5)
C      IPIXEL(2) ----- PIXEL NO. (SEE INSTR FOR IPIXEL(INSTR)(I5)

C      KEY ----------- 0 = LINE/PIXEL TO LAT/LON CONVERSION.
C                      1 = LAT/LON TO LINE/PIXEL CONVERSION.

C      IMC ----------- IMC STATUS (0=IMC ON, 1=IMC OFF)

C      INSTR --------- INSTRUMENT (1=IMAGER, 2=SOUNDER)

C      EPOCH TIME----- YEAR, J-DAY, HOUR, MINUTE, SECONDS

C      START TIME----- YEAR, J-DAY, HOUR, MINUTE, SECONDS
C  TEST FOR USER:INPUT DATA BLOCK(UNNUM DATA BLOCK)..........
C  FIFTY(50) LINES AVAILABLE FOR USE IN INPUT BLOCK..........
C***********************************************************************
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

