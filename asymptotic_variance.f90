MODULE asymptotic_variance
!
USE constants
USE maxlik
USE printing_routines
USE utilities
!
IMPLICIT NONE
!
CONTAINS
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE compute_asymptotic_variance ( theta, llcheck, thetaF, objf, &
        grad, varQMLmat, stderrQML )
    !
    ! Computes the asymptotic variance matrix of the maximum likelihood estimator 
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: theta(num_theta)         ! Tentative estimates
    REAL(8), INTENT(IN) :: llcheck
    REAL(8), INTENT(OUT) :: thetaF(num_theta)        ! Final parameters
    REAL(8), INTENT(OUT) :: objf                    ! Tentative optimized loglikelihood 
    REAL(8), INTENT(OUT) :: grad(num_theta)         ! Gradient at theta
    REAL(8), INTENT(OUT) :: varQMLmat(num_theta,num_theta)
                                                    ! QML estimate of the variance matrix of the parameters
    REAL(8), INTENT(OUT) :: stderrQML(num_theta)    ! QML estimate of the asymptotic standard errors
    !
    ! Declaring local variables
    !
    REAL(8) :: psi(num_psi), pdpsi(num_psi), p2dpsi(num_psi,num_psi)
    REAL(8) :: p, spdl(num_theta,num_theta)
    REAL(8) :: imat(num_theta,num_theta)            ! I
    INTEGER :: n, i, j                              ! Loop indexes
    REAL(8) :: l(num_N)                             ! Loglikelihood contributions
    REAL(8) :: dl(num_N,num_theta)                  ! Loglikelihood contributions gradient 
    REAL(8) :: d2l(num_theta,num_theta)             ! Loglikelihood contributions hessian
    REAL(8) :: hess(num_theta,num_theta)            ! Loglikelihood hessian
    INTEGER :: ipiv(num_theta)                      ! Used in f07mdf
    INTEGER :: lwork_f07mdf                         ! Used in f07mdf
    INTEGER :: ifail1, ifail2                       ! Used in f07mdf
    REAL(8) :: work_f07mdf(num_theta*100)           ! Used in f07mdf
    REAL(8) :: auxmat(num_theta,num_theta)          ! Used in f07mdf
    REAL(8) :: invjmat(num_theta,num_theta)         ! Used in f07mdf
    REAL(8) :: dtheta_dpsiprime(num_theta,num_psi)
    REAL(8) :: dthetaF_dthetaprime(num_theta,num_theta)
    INTEGER :: iflag
    !
    ! Beginning execution
    !
    WRITE (*,1357)
    1357 FORMAT (//, 'Starting evaluation of asymptotic variance matrix', //)
    !
    ! Check the correctness of the parameter values on input
    !
    CALL loglik_fct(num_theta,theta,objf,grad)
    WRITE (*,1212) theta, objf, llcheck
1212 FORMAT('Parameters on input:', /, <num_theta>(ES9.2,1X), /, &
    'With these parameters, the average loglikelihood equals ', ES20.13, /, &
    'Should have been                                        ', ES20.13)
    IF (ABS(objf-llcheck) .GT. 1.d-10) THEN 
        PRINT*, 'Problem:'
        PRINT*, 'The loglikelihood value of the starting parameters'
        PRINT*, 'is not equal to the previously computed value'
        PRINT*, 'Type any key to stop'
        STOP
    END IF
    !
    ! First, computing matrix I 
    !
    WRITE (*,*) 'Computing gradients and hessian ...'
    !
    ! Beginning loop over observations
    !
    imat = 0.d0
    grad = 0.d0
    hess = 0.d0
    dl = 0.d0
    d2l = 0.d0
    DO n = 1, num_N
        !
        IF ((MOD(n,1000) .EQ. 0) .OR. (n .EQ. num_N)) PRINT*, 'n = ', n
        !
        CALL compute_psi_parameters(x_m_s(n,:),x_m_q(n,:),x_sigma_s(n,:), &
            theta,psi)
        dtheta_dpsiprime = matrix_dtheta_dpsiprime(x_m_s(n,:),x_m_q(n,:),x_sigma_s(n,:), &
            theta)    
        CALL INDIVIDUAL_CONTRIBUTION_DV_DV(psi,idmat,idmat,d(n),as(n),riskav(n), &
            p,pdpsi,p2dpsi,num_psi,num_psi)
        !
        dl(n,:) = MATMUL(dtheta_dpsiprime,pdpsi)
        spdl = SPREAD(dl(n,:),DIM = 2,NCOPIES = num_theta)
        dl(n,:) = normpeso(n)*dl(n,:)/p
        grad = grad+dl(n,:)
        !
        d2l = MATMUL(MATMUL(dtheta_dpsiprime,p2dpsi),TRANSPOSE(dtheta_dpsiprime))
        d2l = normpeso(n)*(d2l/p-spdl*TRANSPOSE(spdl)/(p**2))
        hess = hess+d2l
        !
    END DO
    !
    ! Middle matrix in formula (17.73) on p. 594 of Wooldridge (2002)
    !
    imat = MATMUL(TRANSPOSE(dl),dl)     
    CALL open_write_file(unit_dl,file_dl)
    WRITE(unit_dl,1010) (l(n), dl(n,:), n = 1, num_N)
    CLOSE(UNIT=unit_dl)
    CALL open_write_file(unit_imat,file_imat)
    WRITE(unit_imat,1011) (imat(n,:), n = 1, num_theta)
    CLOSE(UNIT=unit_imat)
    CALL open_write_file(unit_jmat,file_jmat)
    WRITE(unit_jmat,1011) (-hess(n,:), n = 1, num_theta)
    CLOSE(UNIT=unit_jmat)
    !
    1010 FORMAT (<num_theta+1>(ES20.13, 1X))
    1011 FORMAT (<num_theta>(ES20.13, 1X))
    !
    ! Inverting the weighted hessian matrix 
    !
    auxmat = -hess
    lwork_f07mdf = num_theta*100
    CALL dsytrf('U',num_theta,auxmat,num_theta,ipiv,work_f07mdf,lwork_f07mdf,ifail1)
    IF (ifail1 .NE. 0) THEN
        PRINT*, "dsytrf can't be applied to -hess"
        STOP
    END IF
    CALL dsytri('U',num_theta,auxmat,num_theta,ipiv,work_f07mdf,ifail2)
    IF (ifail2 .NE. 0) THEN
        PRINT*, "f07mjf can't be applied to -hess"
        STOP
    END IF
    DO i = 2, num_theta
        DO j = 1, i-1
            auxmat(i,j) = auxmat(j,i)
        END DO
    END DO
    invjmat = auxmat
    !
    CALL open_write_file(unit_invjmat,file_invjmat)
    WRITE(unit_invjmat,1011) (invjmat(n,:), n = 1, num_theta)
    CLOSE(UNIT = unit_invjmat)
    !
    ! Computing the QML variance matrix and standard errors
    !
    CALL compute_thetaF_parameters(theta,thetaF)
    dthetaF_dthetaprime = matrix_dthetaF_dthetaprime(theta)
    varQMLmat = MATMUL(MATMUL(invjmat,imat),invjmat)
    varQMLmat = MATMUL(MATMUL(dthetaF_dthetaprime,varQMLmat),TRANSPOSE(dthetaF_dthetaprime))
    DO i = 1, num_theta
        !
        WRITE(unit_varthetaQML,12345) varQMLmat(i,:)
12345   FORMAT(<num_theta>(ES20.13,1X))
        !
    END DO
    DO i = 1, num_theta
        IF (varQMLmat(i,i) .GE. 0.d0) stderrQML(i) = SQRT(varQMLmat(i,i)/REAL(num_N))
    END DO
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE compute_asymptotic_variance
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
END MODULE asymptotic_variance
