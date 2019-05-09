MODULE printing_routines
!
USE constants
USE observations
USE utilities
!
IMPLICIT NONE
!
CONTAINS
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE open_read_file ( unit_number, file_name )
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: unit_number
    CHARACTER(len=20), INTENT(IN) :: file_name
    !
    ! Declaring local variables
    !
    INTEGER :: open_err                 ! Open file error code
    !
    ! Beginning execution
    !
    OPEN ( UNIT=unit_number, FILE=file_name, STATUS='old', IOSTAT=open_err )
    IF (open_err .NE. 0) THEN
        WRITE(*,1234) file_name
        STOP
    END IF
    1234 FORMAT ('Unable to open input file ', A20)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE open_read_file 
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE open_write_file ( unit_number, file_name )
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: unit_number
    CHARACTER(len=20), INTENT(IN) :: file_name
    !
    ! Declaring local variables
    !
    INTEGER :: open_err                 ! Open file error code
    !
    ! Beginning execution
    !
    OPEN ( UNIT=unit_number, FILE=file_name, STATUS='replace', IOSTAT=open_err )
    IF (open_err .NE. 0) THEN
        WRITE(*,1234) file_name
        STOP
    END IF
    1234 FORMAT ('Unable to open input file ', A20)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE open_write_file 
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE print_res ( i_stime, objf, theta, task, grad )
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: i_stime      ! Estimation trial number
    REAL(8), INTENT(IN) :: objf         ! Latent criterion function at the optimum
    REAL(8), INTENT(IN) :: theta(num_theta)
    CHARACTER(len=60), INTENT(IN) :: task
    REAL(8), INTENT(IN) :: grad(num_theta)
    !
    ! Declaring local variables
    !
    REAL(8) :: psi(num_psi)
    REAL(8) :: m_s, m_z, sigma_s, sigma_z, rho_sz
    REAL(8) :: gamma, delta_z(num_delta_z)
    !
    ! Beginning execution
    !
    CALL compute_psi_parameters(x_m_s(1,:),x_m_q(1,:),x_sigma_s(1,:), &
        theta,psi)
    CALL compute_b_parameters(psi,m_s,m_z, &
            sigma_s,sigma_z,rho_sz, &
            gamma,delta_z)
    !
    ! Print res file
    !
    WRITE (unit_res,35) i_stime, objf, theta, &
        sigma_z, rho_sz, &
        delta_z, &
        error_flag, task, grad
35  FORMAT ( I3, " # ", ES25.18, " #", <num_theta>(1X, ES25.18), &
        " #", <2>(1X, ES15.8), &
        " #", <num_delta_z>(1X, ES15.8), &
        " # ", L5, " # ", A60, " # ", <num_theta>(ES15.8,1X))
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE print_res 
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE print_screen_politope ( i_stime, iter, objf )
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: i_stime      ! Estimation trial number
    INTEGER, INTENT(IN) :: iter         ! Iteration number
    REAL(8), INTENT(IN) :: objf         ! Latent criterion function at the optimum
    !
    ! Beginning execution
    !
    WRITE(*,19) i_stime, iter, objf
19  FORMAT('Stima ', I3, ', Iter: ', I5, ' ; -LL/N: ', ES15.8)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE print_screen_politope 
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE print_final_results ( theta, stderr, objf, grad, vmat )
    !
    ! Computes final results and writes them on file
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: theta(num_theta)             ! Estimates of parameters
    REAL(8), INTENT(IN) :: stderr(num_theta)            ! Estimated asymptotic standard errors
    REAL(8), INTENT(IN) :: objf                         ! Optimized latent criterion 
    REAL(8), INTENT(IN) :: grad(num_theta)              ! Gradient at theta_star
    REAL(8), INTENT(IN) :: vmat(num_theta,num_theta)    ! Estimated variance matrix of the parameters
    !
    ! Declaring local variables
    !
    INTEGER :: i_par, i_L, ind                          ! Loop indexes
    CHARACTER(len = 15) :: extra_name
    REAL(8) :: econst, stderr_econst
    !
    ! Beginning execution
    !
    ! Computing statistics for the parameter estimates
    !
    ! m_s
    !
    WRITE (unit_fin_res,3525) 
    3525 FORMAT ( /, &
        '-----------------------------------------------------------------', /, &
        '    Parameter     Estimate        AsySE   AsyT-Ratio     Gradient', /, &
        '-----------------------------------------------------------------', /, &
        '                             M_S                                 ', /, &
        '-----------------------------------------------------------------' )
    
    DO i_par = 1, num_X_m_s
        WRITE (unit_fin_res, 9999) names_m_s(i_par), theta(i_par), stderr(i_par), &
            theta(i_par)/stderr(i_par), grad(i_par)
        9999 FORMAT ( 1X, A12, 4(1X, ES12.5) )
    END DO
    !
    ! m_q
    !
    WRITE (unit_fin_res,35251) 
    35251 FORMAT ( &
        '-----------------------------------------------------------------', /, &
        '                             M_Z                                 ', /, &
        '-----------------------------------------------------------------' )
    !
    ind = num_X_m_s
    DO i_par = ind+1, ind+num_X_m_q
        WRITE (unit_fin_res, 9999) names_m_q(i_par-ind), theta(i_par), stderr(i_par), &
            theta(i_par)/stderr(i_par), grad(i_par)
    END DO
    !
    ! sigma_s
    ! 
    WRITE (unit_fin_res,35261) 
    35261 FORMAT ( &
        '-----------------------------------------------------------------', /, &
        '                          SIGMA_S                                ', /, &
        '-----------------------------------------------------------------' )
    !
    ind = num_X_m_s+num_X_m_q
    DO i_par = ind+1, ind+num_X_sigma_s
        WRITE (unit_fin_res, 9999) names_sigma_s(i_par-ind), theta(i_par), stderr(i_par), &
            theta(i_par)/stderr(i_par), grad(i_par)
    END DO
    !
    ! Add exp(const), in case no covariate is used in sigma_s
    !
    i_par = num_X_m_s+num_X_m_q+1
    extra_name = 'e^const'
    econst = EXP(theta(i_par))
    stderr_econst = econst*stderr(i_par)
    WRITE(unit_fin_res,13131)
    13131 FORMAT('-----------------------------------------------------------------')
    WRITE (unit_fin_res, 9999) extra_name, econst, stderr_econst, &
        econst/stderr_econst, grad(i_par)
    ind = num_X_m_s+num_X_m_q+num_X_sigma_s
    !
    ! sigma_z
    !
    IF (switch_sigma_z .EQ. 1) THEN
        !
        WRITE (unit_fin_res,35262) 
        35262 FORMAT ( &
            '-----------------------------------------------------------------', /, &
            '                           SIGMA_Z                               ', /, &
            '-----------------------------------------------------------------' )
        !
        ind = ind+1
        WRITE (unit_fin_res, 9999) 'const', theta(ind), stderr(ind), &
            theta(ind)/stderr(ind), grad(ind)
        !
    END IF
    !
    ! rho_sz
    !
    IF (switch_rho_sz .EQ. 1) THEN
        !
        WRITE (unit_fin_res,35264) 
        35264 FORMAT ( &
            '-----------------------------------------------------------------', /, &
            '                            RHO_SZ                               ', /, &
            '-----------------------------------------------------------------' )
        !
        ind = ind+1
        WRITE (unit_fin_res, 9999) 'const', theta(ind), stderr(ind), &
            theta(ind)/stderr(ind), grad(ind)
        !
    END IF
    !
    ! delta_z
    ! 
    WRITE (unit_fin_res,35267) 
    35267 FORMAT ( &
        '-----------------------------------------------------------------', /, &
        '                          DELTA_Z                                ', /, &
        '-----------------------------------------------------------------' )
    !
    DO i_par = ind+1, ind+num_theta_delta_z
        WRITE (unit_fin_res, 9991) 'delta_z(', i_par-ind+1, ')  ', theta(i_par), stderr(i_par), &
            theta(i_par)/stderr(i_par), grad(i_par)
        9991 FORMAT ( 1X, A8, I1, A3, 4(1X, ES12.5) )
    END DO
    ind = ind+num_theta_delta_z
    !
    WRITE (unit_fin_res,3524) 
    !
    WRITE (unit_fin_res, 9997) -num_N*objf
    9997 FORMAT ( 1X, 'Total objf  ', 1X, ES12.5 )
    !
    WRITE (unit_fin_res, 9992) -objf
    9992 FORMAT ( 1X, 'Avg. objf   ', 1X, ES12.5 )
    !
    WRITE (unit_fin_res,3524) 
    3524 FORMAT ( &
        '---------------------------------------------------------------------------------------------------------------------' )
    !
    WRITE (unit_fin_res,14121) (vmat(ind,:), ind = 1, num_theta)
    14121 FORMAT (1X, &
        1X, 'Estimated asymptotic variance matrix of theta:', /, &
        <num_theta>(1X,ES12.5))
    !
    WRITE (unit_fin_res,3524) 
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE print_final_results
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
END MODULE printing_routines
