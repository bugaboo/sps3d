      program bin2txt
      INTEGER :: n
      CHARACTER(LEN=30) :: NAME, STR
      complex(16), allocatable :: data (:)
      IF (IARGC() .LT. 2) stop
      CALL GETARG(1, NAME)
      CALL GETARG(2, STR)
      read (str, *) n
      allocate (data(n))
      OPEN(UNIT = 1, FILE=NAME, FORM = 'unformatted', 
     &		ACCESS = 'direct', RECL = n * 16)
      read (1, rec = 1) data
      do i = 1, n
	write (*, 77) real(data(i)), aimag(data(i))
      enddo
 77   FORMAT(4(E19.12,1X))     
      end program