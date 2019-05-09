MODULE starting_points
!
USE constants
USE printing_routines
!
IMPLICIT NONE
!
CONTAINS
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
	SUBROUTINE admissible_starting_point ( seed, theta, theta_inf, theta_sup, llcheck )
	!
	! Computes an admissible starting parameters vector for the PPR iterations
	!
	IMPLICIT NONE
	!
	! Declaring dummy variables
	!
	INTEGER, INTENT(INOUT) :: seed(2)				! Seed for r.n. generation
	REAL(8), INTENT(OUT) :: theta(num_theta)		! Admissible value of the parameters vector
	REAL(8), INTENT(IN) :: theta_inf(num_theta)		! Lower bounds for generating completely 
													! random starting points
	REAL(8) , INTENT(IN):: theta_sup(num_theta)		! Upper bounds for generating completely
													! random starting points
	REAL(8), INTENT(OUT) :: llcheck
    !
    ! Declaring external routine
    !
    REAL(8), EXTERNAL :: r8_uniform_01
    !
	! Declaring local variables
	!
	INTEGER :: i									! Index
    CHARACTER(len = 2) :: ichar
	!
	! Beginning execution
	!
    ! Modalità di test
    !
    IF (to0 .EQ. 1) THEN
        !
        CALL open_read_file(unit_theta,file_theta)
        READ(unit_theta,*) (theta(i), i = 1, num_theta)
        CLOSE(UNIT = unit_theta)
        !
    END IF
    !
	! Step 1: Starting points chosen randomly
	!
    IF (to1 .EQ. 1) THEN
		!
        CALL initialize()
        CALL set_initial_seed(seed(1),seed(2))
		DO i = 1, num_theta
            theta(i) = theta_inf(i)+(theta_sup(i)-theta_inf(i))*r8_uniform_01()
		END DO 
		!
    END IF
	!
	! Step 2: Starting point read from the res_to1.txt file
	!
	IF (to2 .EQ. 1) THEN
        !
        CALL GETARG(1,ichar)
        i = INUM(ichar)
        !
        CALL open_read_file(unit_res_to1,file_res_to1)
        IF (i .GT. 1) THEN
            READ(unit_res_to1, 2534) 
2534        FORMAT(<i-1>(/))
            BACKSPACE(UNIT=unit_res_to1)
        END IF
        READ(unit_res_to1,*) llcheck, theta
        CLOSE(UNIT = unit_res_to1)
        !
	END IF 
	!
	! Step 3: Starting point read from the res_to2.txt file
	!
	IF (to3 .EQ. 1) THEN
        !
        CALL GETARG(1,ichar)
        i = INUM(ichar)
        !
        CALL open_read_file(unit_res_to2,file_res_to2)
        IF (i .GT. 1) THEN
            READ(unit_res_to2, 2535) 
2535        FORMAT(<i-1>(/))
            BACKSPACE(UNIT=unit_res_to2)
        END IF
        READ(unit_res_to2,11) llcheck, theta
11      FORMAT(6X,ES25.18,3X,<num_theta>(ES25.18,1X))
        CLOSE(UNIT = unit_res_to2)
        CALL open_write_file(unit_theta,file_theta)
        WRITE(unit_theta,117) theta
117     FORMAT(<num_theta>( ES25.18,1X))
        CLOSE(unit_theta)
        !
    END IF 
    !
END SUBROUTINE admissible_starting_point
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
END MODULE starting_points