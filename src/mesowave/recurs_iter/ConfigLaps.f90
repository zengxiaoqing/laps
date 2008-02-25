MODULE ConfigLaps

!==========================================================
!  This module sets up LAPS configuration and access to its
!  LSO observation data.
!
!  HISTORY: MAY. 2004 by YUANFU XIE 
!           adapted from Dr. McGinley's codes.
!==========================================================

  USE Definition

  IMPLICIT NONE

  ! LAPS variables:
  INTEGER, PARAMETER :: nvlaps = 6 	!number of variables
  INTEGER, PARAMETER :: ncycles = 6     !number of cycles to carry
  INTEGER   :: nx,ny,nfic 			! Grid dimensions
  INTEGER   :: laps_cycle_time,len,i1,i2,i3,i4,istatus,ilaps,ivar
  INTEGER   :: maxstations 		! maximum number of ob stations
  INTEGER   :: istarttime
  CHARACTER :: ext_s*30,dir_s*200      
  CHARACTER :: units*60,comment*60,name*100
  CHARACTER :: maproj*6,nest7grid*9,var_s*3 ! map projection character string
  REAL      :: stanlat,stanlat2,stanlon,badflag,grid_spacingx,grid_spacingy
  ! gridded lat, lon , observation arrary o (variables, stations*time)
  REAL, ALLOCATABLE, DIMENSION (:,:) :: lat,lon,ldf,olaps
  ! observation lat, lon, ob time and weight
  REAL, ALLOCATABLE, DIMENSION (:) :: olat,olon,otime,wght

  ! Background fields:
  REAL,   ALLOCATABLE, DIMENSION (:,:,:,:) :: bkgd

CONTAINS

  INCLUDE 'GridBarnes.f90'
  INCLUDE 'LapsInfo.f90'
  INCLUDE 'LSO_Data_QC.f90'
  INCLUDE 'WriteAnalysis.f90'

END MODULE
