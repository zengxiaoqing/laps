

!       LAPS Parameters
!       Note: Do not include this in the same module as 'lapsparms.cmn'
!             The same variable names might still used in a few cases?


!       LAPS Grid Dimensions
!       Notes: The character declaration should equal the length of the string.
!              There should be no trailing blanks in the domain name.
!              The domain name should be in lower case.
        character*9 LAPS_DOMAIN_FILE
        parameter  (LAPS_DOMAIN_FILE = 'nest7grid')

!       integer*4 DOMAIN_NAME_LEN
!       DOMAIN_NAME_LEN should be exactly the length of the domain name
!       parameter (DOMAIN_NAME_LEN = 9)
