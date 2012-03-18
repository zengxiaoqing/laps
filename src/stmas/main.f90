! This analysis program is based on the idea of Dr. Xie in Xie et al. (2005) paper         !
! 'A SEQUENTIAL VARIATIONAL ANALYSIS APPROACH FOR MESOSCALE DATA ASSIMILATION', in         !
! which an innovation of data assimilation has been proposed that data assimilation        !
! should be decomposed to two step, with the first step to catch information provided      !
! by observations and the second step to combine background and observations to make       !
! the final optimal analysis. This code is to implement the first step which can provide   !
! reasonable analysis field. For high efficiency, we choose the multi-grid method to       !
! implement this idea. The previous version of this code is the STMAS surface analysis     !
! system coded by Dr. Xie himself. Some applicaiton of this surface STMAS in ocean can     !
! be refered to Li et al. (2007) and He et al. (2007). In this four-dimensional version,   !
! we use some penalty term following the suggestion in Xie et al. (2002) paper 'IMPACT     !
! OF FORMULATION OFCOST FUNCTION AND CONSTRAIONTS ON THREE-DIMENSIONAL VARIATIONAL DATA    !
! ASSIMIALTION'. So far, this four-dimensional version can deal with analysis from 1D to   !
! 4D, multi-variable anaysis using balance and general vertical coordinate.                !


! NAME OF LOCAL VARIABLE               HAS  1 OR 2 CHARACTERS
! PHYSICAL SCALE VARIABLE              HAS     3   CHARACTERS
! NAME OF I/O DEVICE UNIT              HAS     5   CHARACTERS
! NAME OF INTEGER GLOBAL VARIALBE      HAS     7   CHARACTERS
! NAME OF PARAMETER                    HAS     7   CHARACTERS
! NAME OF GLOBAL ARRAY                 HAS     7   CHARACTERS
! NAME OF OBSERVATION OR GRID VARIABLE HAS     8   CHARACTERS
! NAME OF SUBROUTINE                   HAS     9   CHARACTERS
! NAME OF MODULE                       HAS    12   CHARACTERS

PROGRAM MAIN

  USE INPUT_BG_OBS, ONLY : BKGRNDOBS
  USE PRMTRS_STMAS, ONLY : IF_TEST
  USE STMAS4D_CORE, ONLY : MGANALYSS
  USE OUTPUT_ANALS, ONLY : TMPOUTPUT,OUTPTLAPS,OUTPUTANA

  CALL BKGRNDOBS
  CALL MGANALYSS
  IF (IF_TEST .EQ. 0) THEN
    CALL OUTPTLAPS
  ELSE
    CALL OUTPUTANA
  ENDIF
  PRINT*,'END OF ANALYSIS'

END PROGRAM MAIN
