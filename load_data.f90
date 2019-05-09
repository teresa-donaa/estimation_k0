MODULE load_data
!
USE constants
USE observations
!
IMPLICIT NONE
!
CONTAINS
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE input_model ( sel_X )
	!
    ! Reads the model's specification from text file
    !
	IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(OUT) :: sel_X(num_tot_X,num_states)
	! 
	! Declaring local variables
	!
	INTEGER :: open_err								! Open error code
	INTEGER :: n, i_x, i_tot_x, i					! Loop index
	REAL(8) :: tot_x(num_N,num_tot_X)				! Matrix of available explanatory variables
	REAL(8) :: tot_xn(num_tot_X-1)
	!
	! Beginning execution
	!
	OPEN(UNIT = unit_model,FILE = file_model,STATUS = 'old',ACTION = 'read',IOSTAT=open_err)
	IF (open_err .NE. 0) THEN
		WRITE(*,1) file_model
1       FORMAT('PROBLEM: Unable to open file ', A20)
		STOP
	END IF
	!
	! First, skip the first row (containing names of states variables)
	!
	READ(unit_model,1001)
	1001 FORMAT()
	!
	! Second, read the following num_tot_K rows, num_states columns
	!
	DO n = 1, num_tot_X
        READ(unit_model,*) sel_X(n,:)
	    PRINT*, 'n = ', n
    END DO
	!
	CLOSE(UNIT = unit_model)
    !
    ! Compute number of covariates in state variables
    !
    num_X_m_s = SUM(sel_X(:,1))
    num_X_m_q = SUM(sel_X(:,2))
    num_X_sigma_s = SUM(sel_X(:,3))
    !
    ! Compute the number of parameters
    !
    num_theta_beta = num_X_m_s+num_X_m_q
    num_theta_sigma_s = num_X_sigma_s
    num_theta = num_theta_beta+num_theta_sigma_s+num_theta_sigma_z+num_theta_rho+num_theta_delta_z
    !
    ! Parameters for POLITOPE
    ! 
    max_iters_politope = (100*num_theta)*to1+(600*num_theta)*to2
    max_repeat_politope = num_theta*to1+8*num_theta*to2
    !
    ! Ending execution and returning control
    !
	END SUBROUTINE input_model
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE input_data ( theta_inf, theta_sup )
	!
	IMPLICIT NONE
	! 
	! Declaring dummy variables
	!
	REAL(8), INTENT(OUT) :: theta_inf(num_theta)	! Lower bounds for generating completely 
													! random starting points
	REAL(8), INTENT(OUT) :: theta_sup(num_theta)	! Upper bounds for generating completely
													! random starting points
	! 
	! Declaring local variables
	!
	INTEGER :: open_err								! Open error code
	INTEGER :: n, i_x, i_tot_x, i					! Loop index
	REAL(8) :: tot_x(num_N,num_tot_X)				! Matrix of available explanatory variables
	REAL(8) :: tot_xn(num_tot_X-1)
    INTEGER :: horizon
	!
	! Beginning execution
	!
	! Opening data file
	!
	OPEN(UNIT = unit_data,FILE = file_data,STATUS = 'old',ACTION = 'read',IOSTAT=open_err)
	IF (open_err .NE. 0) THEN
		WRITE(*,1) file_data
1       FORMAT('PROBLEM: Unable to open file ', A20)
		STOP
	END IF
	!
	! Reading data
	!
	! First, skip the first row (containing names)
	!
	READ (unit_data,1001)
	1001 FORMAT ()
	!
	! Second, read the following num_N rows
	!
    tot_x(:,1) = 1.d0
	DO n = 1, num_N
        !
        READ (unit_data,*) d(n), ab(n), as(n), riskav(n), horizon, tot_xn, peso(n)
        IF ((MOD(n,1000) .EQ. 0) .OR. (n .EQ. num_N)) PRINT*, 'n = ', n, ' ; peso = ', peso(n)
        tot_x(n,2:) = tot_xn
        !
    END DO
	!
	CLOSE(UNIT = unit_data)
    normpeso = peso/DBLE(SUM(peso))
	!
    ! Select explanatory variables in m_s
    !
    i_x = 0
    DO i_tot_x = 1, num_tot_X
        !
        IF (sel_X(i_tot_x,1) .EQ. 1) THEN
            !
            i_x = i_x+1
            x_m_s(:,i_x) = tot_x(:,i_tot_x)
            !
        END IF
        !
    END DO
	!
    ! Select explanatory variables in m_q
    !
    i_x = 0
    DO i_tot_x = 1, num_tot_X
        !
        IF (sel_X(i_tot_x,2) .EQ. 1) THEN
            !
            i_x = i_x+1
            x_m_q(:,i_x) = tot_x(:,i_tot_x)
            !
        END IF
        !
    END DO
	!
    ! Select explanatory variables in sigma_s
    !
    i_x = 0
    DO i_tot_x = 1, num_tot_X
        !
        IF (sel_X(i_tot_x,3) .EQ. 1) THEN
            !
            i_x = i_x+1
            x_sigma_s(:,i_x) = tot_x(:,i_tot_x)
            !
        END IF
        !
    END DO
    !
	! Computing bounds for random starting values of the parameters
	!
    theta_inf = -1.d0
    theta_sup = 1.d0
	!
    ! Defining identity matrix (num_psi x num_psi) to be used to compute derivatives
    !
    DO i = 1, num_psi
        !
        idmat(i,:) = 0.d0
        idmat(i,i) = 1.d0
        !
    END DO    
    !
    ! Ending execution and returning control
    !
	END SUBROUTINE input_data
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
	SUBROUTINE input_names ( )
	!
	IMPLICIT NONE
	! 
	! Declaring local variables
	!
	INTEGER :: open_err								! Open error code
	INTEGER :: i_x, i_tot_x							! Loop index
	CHARACTER(len = 15) :: tot_names(num_tot_X)		! Available explanatory variables' names
	!
	! Beginning execution
	!
	! Opening names file
	!
	OPEN(UNIT = unit_names,FILE = file_names,STATUS = 'old',ACTION = 'read',IOSTAT = open_err)
	IF (open_err .NE. 0) THEN
		WRITE(*,1) file_names
1       FORMAT('PROBLEM: Unable to open file ', A20)
		STOP
	END IF
	!
	! Reading names
	!
	DO i_tot_x = 1, num_tot_X
        !
		READ (unit_names,*) tot_names(i_tot_x)
        !
	END DO
	!
	CLOSE(UNIT = unit_names)
	!
    ! Select names of explanatory variables in m_s
    !
    i_x = 0
    DO i_tot_x = 1, num_tot_X
        !
        IF (sel_X(i_tot_x,1) .EQ. 1) THEN
            !
            i_x = i_x+1
            names_m_s(i_x) = tot_names(i_tot_x)
            !
        END IF
        !
    END DO
	!
    ! Select names of explanatory variables in m_q
    !
    i_x = 0
    DO i_tot_x = 1, num_tot_X
        !
        IF (sel_X(i_tot_x,2) .EQ. 1) THEN
            !
            i_x = i_x+1
            names_m_q(i_x) = tot_names(i_tot_x)
            !
        END IF
        !
    END DO
	!
    ! Select names of explanatory variables in sigma_s
    !
    i_x = 0
    DO i_tot_x = 1, num_tot_X
        !
        IF (sel_X(i_tot_x,3) .EQ. 1) THEN
            !
            i_x = i_x+1
            names_sigma_s(i_x) = tot_names(i_tot_x)
            !
        END IF
        !
    END DO
	!
	! Ending execution and returning control
	!
	END SUBROUTINE input_names
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
END MODULE load_data