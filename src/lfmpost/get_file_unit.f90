SUBROUTINE get_file_unit(unit)

  ! Simple subroutine to get the next available file unit number

  IMPLICIT NONE
  INTEGER, PARAMETER   :: min_unit = 7
  INTEGER, PARAMETER   :: max_unit = 1023
  INTEGER, INTENT(OUT) :: unit
  LOGICAL              :: unit_used
  LOGICAL              :: unit_found

  unit_found = .false.
  unit_loop: DO unit = min_unit, max_unit
    INQUIRE(UNIT=unit,OPENED=unit_used)
    IF (.NOT. unit_used) THEN
      unit_found = .TRUE.
      EXIT unit_loop
    ENDIF
  ENDDO unit_loop
  IF (.NOT. unit_found) THEN
    PRINT '(A)', 'No available unit numbers!'
    STOP
  ENDIF
  RETURN
END SUBROUTINE get_file_unit
