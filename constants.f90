MODULE constants
!
IMPLICIT NONE
! 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Selecting optimization method and other switches
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
INTEGER, PARAMETER :: to0 = 0           ! Modalità di test
INTEGER, PARAMETER :: to1 = 0           ! Prime stime
INTEGER, PARAMETER :: to2 = 0           ! Seconde stime
INTEGER, PARAMETER :: to3 = 1           ! Varianza
!
INTEGER, PARAMETER :: num_stime = to0*1+to1*100+to2*1+to3*0
                                        ! Total number of completed estimation trials
INTEGER, PARAMETER :: compute_var_as = to0*0+to1*0+to2*0+to3*1      
                                        ! Switch to compute the asymptotic variance matrix
!
! Selection of optimization algorithm
!
INTEGER, PARAMETER :: switch_lbfgs = 1          ! = 1: L-BFGS optimization ON; = 0 L-BFGS optimization OFF
INTEGER, PARAMETER :: switch_politope = 1       ! = 1: POLITOPE optimization ON; = 0 POLITOPE optimization OFF
! 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declaring sample characteristics and selecting explanatory variables
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! SCF 1995 - 2016 data
!
INTEGER, PARAMETER :: num_N = 24811     ! mv4 data 95-16
INTEGER, PARAMETER :: num_tot_X = 36    ! Total number of available explanatory variables
INTEGER, PARAMETER :: num_L = 3         ! Degrees of self-reported risk aversion
INTEGER, PARAMETER :: num_states = 3    ! Number of state variables: (m_s,m_z,sigma_s)
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Defining the number of parameters
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! variance parameters
!
INTEGER, PARAMETER :: switch_sigma_z = 0        ! = 1: sigma_z unconstrained; = 0: sigma_z = 1
INTEGER, PARAMETER :: num_theta_sigma_z = switch_sigma_z
!
! correlation parameters
!
INTEGER, PARAMETER :: switch_rho_sz = 1         ! = 1: rho_sz unconstrained; = 0: rho_sz = 0
INTEGER, PARAMETER :: num_theta_rho = switch_rho_sz
!
! threshold parameters
!
INTEGER, PARAMETER :: num_delta_z = num_L-1
INTEGER, PARAMETER :: num_theta_delta_z = num_L-2
!
! psi
!
INTEGER, PARAMETER :: num_psi_m = 3       
INTEGER, PARAMETER :: num_psi_sigma = 1+num_theta_sigma_z
INTEGER, PARAMETER :: num_psi_rho = num_theta_rho
INTEGER, PARAMETER :: num_psi = num_psi_m+num_psi_sigma+num_psi_rho+num_theta_delta_z
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declaring optimization switches
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 
REAL(8), PARAMETER :: factr = (to0+to1)*1.d+7+(to2+to3)*1.d+1
REAL(8), PARAMETER :: pgtol = (to0+to1)*1.0d-4+(to2+to3)*1.d-5
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declaring parameters for POLITOPE
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 
REAL(8), PARAMETER :: tol = 1.d-7*to1+1.d-11*to2
                                                ! Tolerance in second round optimization
REAL(8), PARAMETER :: tol_conv = 1.d-6*to1+1.d-10*to2
                                                ! Convergence tolerance in second round optimization
REAL(8), PARAMETER :: tol_politope_p = tol
REAL(8), PARAMETER :: crit_politope_conv_p = tol_conv
REAL(8), PARAMETER :: tol_politope_y = tol
REAL(8), PARAMETER :: crit_politope_conv_y = tol_conv
INTEGER, PARAMETER :: rtol_formula = 2          ! Chooses the formula used to compute rtol
                                                ! See lines 150-190 in simplex_M.f90
INTEGER, PARAMETER :: crit_conv_formula = 1     ! Politope minimizations are restarted looking 
                                                ! at improvements in: 
                                                ! 1 = y
                                                ! 2 = p
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declaring parameters about the output files
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 
CHARACTER(len=30), PARAMETER :: file_data = 'data_2.txt'
CHARACTER(len=30), PARAMETER :: file_names = 'names_data_2.txt'
CHARACTER(len=30), PARAMETER :: file_model = 'model_2.txt'
CHARACTER(len=30), PARAMETER :: file_theta = 'theta.txt'
!
CHARACTER(len=20), PARAMETER :: file_res = 'res.txt'
CHARACTER(len=20), PARAMETER :: file_res_to1 = 'res_to1.txt'
CHARACTER(len=20), PARAMETER :: file_res_to2 = 'res_to2.txt'
CHARACTER(len=20), PARAMETER :: file_fin_res = 'fin_res.txt'
CHARACTER(len=20), PARAMETER :: file_varthetaQML = 'varthetaQML.txt'
CHARACTER(len=20), PARAMETER :: file_imat = 'imat.txt'
CHARACTER(len=20), PARAMETER :: file_jmat = 'jmat.txt'
CHARACTER(len=20), PARAMETER :: file_invjmat = 'invjmat.txt'
CHARACTER(len=20), PARAMETER :: file_dl = 'dl.txt'
CHARACTER(len=20), PARAMETER :: file_loglik = 'loglik.txt'
CHARACTER(len=30), PARAMETER :: file_politope = 'politope.txt'
!
INTEGER, PARAMETER :: unit_data = 1             
INTEGER, PARAMETER :: unit_names = 2            
INTEGER, PARAMETER :: unit_model = 3
INTEGER, PARAMETER :: unit_theta = 4
INTEGER, PARAMETER :: unit_politope = 5
!
INTEGER, PARAMETER :: unit_res = 21             
INTEGER, PARAMETER :: unit_res_to1 = 22         
INTEGER, PARAMETER :: unit_res_to2 = 23         
INTEGER, PARAMETER :: unit_fin_res = 24         
INTEGER, PARAMETER :: unit_varthetaQML = 25
INTEGER, PARAMETER :: unit_imat = 26            
INTEGER, PARAMETER :: unit_jmat = 27            
INTEGER, PARAMETER :: unit_invjmat = 28     
INTEGER, PARAMETER :: unit_dl = 29      
INTEGER, PARAMETER :: unit_loglik = 30
! 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declaring constants
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
REAL(8), PARAMETER :: minimum_p = 1.d-10                                ! Minimum probability value
REAL(8), PARAMETER :: minimum_sigma_s = 1.d-6                           ! Minimum sigma_s value
REAL(8), PARAMETER :: greek_pi = 3.14159265358979323846264338328d0      ! Pi
REAL(8), PARAMETER :: twooverpi = 0.636619772367581343075535053490d0    ! 2/Pi
REAL(8), PARAMETER :: invsqrtpi = 0.564189583547756286948079451561d0    ! 1/Sqrt[Pi]
REAL(8), PARAMETER :: sqrt2 = 1.41421356237309504880168872421d0         ! Sqrt[2]
REAL(8), PARAMETER :: invsqrt2pi = 0.398942280401432677939946059934d0   ! 1/Sqrt[2*Pi]
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
END MODULE constants
