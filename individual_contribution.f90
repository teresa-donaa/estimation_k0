SUBROUTINE individual_contribution ( psi, dn, asn, riskavn, p )
!
USE constants
USE observations
USE utilities
USE sp
!
IMPLICIT NONE
! 
! Declaring dummy variables
!   
REAL(8), INTENT(IN) :: psi(num_psi)         ! Vector of metaparameters
INTEGER, INTENT(IN) :: dn                   ! Observed regime dummy
REAL(8), INTENT(IN) :: asn
INTEGER, INTENT(IN) :: riskavn              ! Observed self-reported risk attitude
REAL(8), INTENT(OUT) :: p                   ! Computed regime probability
!
! Declaring external routines
!
REAL(8), EXTERNAL :: TVTL                   ! Trivariate normal cdf
REAL(8), EXTERNAL :: BVND                   ! Bivariate normal survival sdf
REAL(8), EXTERNAL :: PHID                   ! Univariate normal cdf
! 
! Declaring local variables
! 
REAL(8) :: m_s, m_q, m_z 
REAL(8) :: sigma_s, sigma_z
REAL(8) :: rho_sz
REAL(8) :: gamma
REAL(8) :: delta_z(num_delta_z)
REAL(8) :: m_zs, m_ys, psi_z, psi_y, rho_psi_zy
REAL(8) :: p_as, p_zy
REAL(8) :: limit2(2,2), corr, pcomp(2,2)
! 
! Beginning execution
! 
! Computing the b parameters
! 
CALL compute_b_parameters(psi,m_s,m_z, &
        sigma_s,sigma_z,rho_sz, &
        gamma,delta_z)
pcomp = 0.d0
limit2 = 0.d0
! 
! Computing individual contribution to the loglikelihood
! 
SELECT CASE (dn)
    !
    CASE (1)   ! 0 < as < 1          
        ! 
        ! density of a_s
        ! 
        p_as = pdf_N01((gamma*asn-m_s)/sigma_s)/sigma_s
        !
        ! probability of self-reported risk aversion and planning horizon
        !
        m_zs = m_z+rho_sz*sigma_z/sigma_s*(m_s-gamma*asn)
        psi_z = sigma_z*SQRT(1.d0-rho_sz**2)
        !
        IF (riskavn .EQ. 1) THEN
            !
            limit2(2,1) = (delta_z(1)-m_zs)/psi_z
            pcomp(1,1) = PHID(limit2(2,1))
            !
        ELSE IF (riskavn .EQ. 2) THEN
            !
            limit2(2,1) = (delta_z(2)-m_zs)/psi_z
            limit2(2,2) = (delta_z(1)-m_zs)/psi_z
            pcomp(1,1) = PHID(limit2(2,1)) 
            pcomp(2,1) = PHID(limit2(2,2))
            !
        ELSE IF (riskavn .EQ. 3) THEN
            !
            limit2(2,2) = (delta_z(2)-m_zs)/psi_z
            pcomp(1,1) = 1.d0
            pcomp(2,1) = PHID(limit2(2,2))
            !
        END IF
        !
        p_zy = pcomp(1,1)-pcomp(2,1)
        ! 
    CASE (2)   ! as = 0
        ! 
        ! joint probability
        ! 
        p_zy = 1.d0
        corr = rho_sz
        limit2(1,:) = -m_s/sigma_s
        !
        IF (riskavn .EQ. 1) THEN
            !
            limit2(2,1) = (delta_z(1)-m_z)/sigma_z
            pcomp(1,1) = BVND(-limit2(1,1),-limit2(2,1),corr)
            !
        ELSE IF (riskavn .EQ. 2) THEN
            !
            limit2(2,1) = (delta_z(2)-m_z)/sigma_z
            limit2(2,2) = (delta_z(1)-m_z)/sigma_z
            pcomp(1,1) = BVND(-limit2(1,1),-limit2(2,1),corr)
            pcomp(2,1) = BVND(-limit2(1,1),-limit2(2,2),corr)
            !
        ELSE IF (riskavn .EQ. 3) THEN
            !
            limit2(2,2) = (delta_z(2)-m_z)/sigma_z
            pcomp(1,1) = PHID(limit2(1,1)) 
            pcomp(2,1) = BVND(-limit2(1,1),-limit2(2,2),corr)
            !
        END IF
        !
        p_as = pcomp(1,1)-pcomp(2,1)
        ! 
    CASE (3)   ! as = 1
        ! 
        ! joint probability
        ! 
        p_zy = 1.d0
        corr = rho_sz
        limit2(1,:) = (gamma-m_s)/sigma_s
        !
        IF (riskavn .EQ. 1) THEN
            !
            limit2(2,1) = (delta_z(1)-m_z)/sigma_z
            pcomp(1,1) = PHID(limit2(2,1))
            pcomp(1,2) = BVND(-limit2(1,1),-limit2(2,1),corr)
            !
        ELSE IF (riskavn .EQ. 2) THEN
            !
            limit2(2,1) = (delta_z(2)-m_z)/sigma_z
            limit2(2,2) = (delta_z(1)-m_z)/sigma_z
            pcomp(1,1) = PHID(limit2(2,1))
            pcomp(2,1) = PHID(limit2(2,2))
            pcomp(1,2) = BVND(-limit2(1,1),-limit2(2,1),corr)
            pcomp(2,2) = BVND(-limit2(1,1),-limit2(2,2),corr)
            !
        ELSE IF (riskavn .EQ. 3) THEN
            !
            limit2(2,2) = (delta_z(2)-m_z)/sigma_z
            pcomp(1,1) = 1.d0
            pcomp(2,1) = PHID(limit2(2,2))
            pcomp(1,2) = PHID(limit2(1,1))
            pcomp(2,2) = BVND(-limit2(1,1),-limit2(2,2),corr)
            !
        END IF
        !
        p_as = pcomp(1,1)-pcomp(2,1)-(pcomp(1,2)-pcomp(2,2))
        ! 
END SELECT
!
! Evaluating the integrand
!
p = p_as*p_zy+minimum_p
! 
! Ending subroutine and returning execution
! 
END SUBROUTINE individual_contribution
