MODULE sp
!
IMPLICIT NONE
!
CONTAINS
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    FUNCTION cdf_N01 ( x, upper )
    !
    ! cdf_N01 computes the cumulative density of the standard normal distribution.
    ! It can be differentiated automatically by TAPENADE
    !
    ! cdf_N01 is based on the ALNORM function
    ! Original FORTRAN77 version by David Hill.
    ! FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, is one endpoint of the semi-infinite interval
    !    over which the integration takes place.
    !
    !    Input, logical UPPER, determines whether the upper or lower
    !    interval is to be integrated:
    !    .TRUE.  => integrate from X to + Infinity;
    !    .FALSE. => integrate from - Infinity to X.
    !
    !    Output, real ( kind = 8 ) ALNORM, the integral of the standard normal
    !    distribution over the desired interval.
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: x
    LOGICAL, INTENT(IN) :: upper
    !
    ! Declaring local variables
    !
    LOGICAL :: up
    REAL(8) :: y, z
    !
    ! Declaring parameters
    !
    REAL(8), PARAMETER :: a1 = 5.75885480458D+00
    REAL(8), PARAMETER :: a2 = 2.62433121679D+00
    REAL(8), PARAMETER :: a3 = 5.92885724438D+00
    REAL(8), PARAMETER :: b1 = -29.8213557807D+00
    REAL(8), PARAMETER :: b2 = 48.6959930692D+00
    REAL(8), PARAMETER :: c1 = -0.000000038052D+00
    REAL(8), PARAMETER :: c2 = 0.000398064794D+00
    REAL(8), PARAMETER :: c3 = -0.151679116635D+00
    REAL(8), PARAMETER :: c4 = 4.8385912808D+00
    REAL(8), PARAMETER :: c5 = 0.742380924027D+00
    REAL(8), PARAMETER :: c6 = 3.99019417011D+00
    REAL(8), PARAMETER :: con = 1.28D+00
    REAL(8), PARAMETER :: d1 = 1.00000615302D+00
    REAL(8), PARAMETER :: d2 = 1.98615381364D+00
    REAL(8), PARAMETER :: d3 = 5.29330324926D+00
    REAL(8), PARAMETER :: d4 = -15.1508972451D+00
    REAL(8), PARAMETER :: d5 = 30.789933034D+00
    REAL(8), PARAMETER :: ltone = 7.0D+00
    REAL(8), PARAMETER :: p = 0.398942280444D+00
    REAL(8), PARAMETER :: q = 0.39990348504D+00
    REAL(8), PARAMETER :: r = 0.398942280385D+00
    REAL(8), PARAMETER :: utzero = 18.66D+00
    !
    ! Declaring function's type
    !
    REAL(8) :: cdf_N01
    !
    ! Beginning execution
    !
    up = upper
    z = x
    !
    IF (z < 0.0D+00) THEN
        up = .not. up
        z = - z
    END IF
    !
    IF ((ltone < z) .AND. ((.NOT. up) .OR. (utzero < z))) THEN
        !
        IF (up) THEN
            cdf_N01 = 0.0D+00
        ELSE
            cdf_N01 = 1.0D+00
        END IF
        RETURN
        !
    END IF
    !
    y = 0.5D+00*z*z
    !
    IF (z <= con) THEN
        cdf_N01 = 0.5D+00-z*(p-q*y/(y+a1+b1/(y+a2+b2/(y+a3))))
    ELSE
        cdf_N01 = r*EXP(-y)/(z+c1+d1/(z+c2+d2/(z+c3+d3/(z+c4+d4/(z+c5+d5/(z+c6))))))
    END IF
    IF (.NOT. up) cdf_N01 = 1.0D+00-cdf_N01
    RETURN
    !
    ! Ending execution
    !
    END FUNCTION cdf_N01
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    FUNCTION pdf_N01 ( x )
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: x
    !
    ! Declaring function's type
    !
    REAL(8) :: pdf_N01
    !
    ! Declaring constants
    !
    REAL(8), PARAMETER :: pi = 3.14159265358979323846264338328d0        
    !
    ! Beginning execution
    !
    pdf_N01 = 1.d0/SQRT(2.d0*pi)*EXP(-x**2/2.d0)
    !
    ! Ending execution and returning control
    !
    END FUNCTION pdf_N01
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    FUNCTION cdf_biv_N01 ( u1, u2, r )
    !
    ! cdf_biv_N01 computes the (upper tail) probability for two normal variates X and Y
    ! whose correlation is R, that X >= u1 and Y >= u2.
    ! It can be differentiated automatically by TAPENADE    
    !
    ! Original FORTRAN77 version by Thomas Donnelly.
    ! FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !    Input, real ( kind = 8 ) u1, u2, the lower limits of integration.
    !    Input, real ( kind = 8 ) R, the correlation between X and Y.
    !    Output, real ( kind = 8 ) cdf_biv_N01, the bivariate normal CDF.
    !
    !  Local Parameters:
    !    Local, integer ( kind = 4 ) IDIG, the number of significant digits
    !    to the right of the decimal point desired in the answer.
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: u1
    REAL(8), INTENT(IN) :: u2
    REAL(8), INTENT(IN) :: r
    !
    ! Declaring local variables
    !
    INTEGER(4) :: i, is
    REAL(8) :: a2, ap, b, cn, con, conex, ex, g2, gh, gk, gw, h2, h4
    REAL(8) :: rr, s1, s2, sgn, sn, sp, sqr, t, w2, wh, wk
    ! 
    ! Declaring function's type
    !
    REAL(8) :: cdf_biv_N01
    !
    ! Declaring parameters
    !
    INTEGER(4), PARAMETER :: idig = 15
    REAL(8), PARAMETER :: twopi = 6.283185307179587D+00
    !
    ! Beginning execution
    !
    b = 0.0D+00
    gh = cdf_N01(-u1,.FALSE.)/2.0D+00
    gk = cdf_N01(-u2,.FALSE.)/2.0D+00
    IF (r == 0.0D+00) THEN
        b = 4.0D+000*gh*gk
        b = MAX(b,0.0D+00)
        b = MIN(b,1.0D+00)
        cdf_biv_N01 = b
        RETURN
    END IF
    !
    rr = (1.0D+00+r)*(1.0D+00-r)
    !
    IF (rr < 0.0D+00) THEN
        WRITE(*,'(a)') ' '
        WRITE(*,'(a)') 'cdf_biv_N01 - Fatal error!'
        WRITE(*,'(a)') '  1 < |R|.'
        STOP
    END IF
    !
    IF (rr == 0.0D+00) THEN
        IF (r < 0.0D+00) THEN
            IF (u1+u2 < 0.0D+00) THEN
                b = 2.0D+00*(gh+gk)-1.0D+00
            END IF
        ELSE
            IF (u1-u2 < 0.0D+00) THEN
                b = 2.0D+00*gk
            ELSE
                b = 2.0D+00*gh
            END IF
        END IF
        !
        b = MAX(b,0.0D+00)
        b = MIN(b,1.0D+00)
        cdf_biv_N01 = b
        RETURN
    END IF
    !
    sqr = SQRT(rr)
    !
    IF (idig == 15) THEN
        con = twopi*1.0D-15/2.0D+00
    ELSE
        con = twopi/2.0D+00
        DO i = 1, idig
            con = con/10.0D+00
        END DO
    END IF
    !
    !  (0,0)
    !
    IF ((u1 == 0.0D+00) .AND. (u2 == 0.0D+00)) THEN
        b = 0.25D+00+ASIN(r)/twopi
        b = MAX(b,0.0D+00)
        b = MIN(b,1.0D+00)
        cdf_biv_N01 = b
        RETURN
    END IF
    !
    !  (0,nonzero)
    !
    IF ((u1 == 0.0D+00) .AND. (u2 /= 0.0D+00)) THEN
        b = gk
        wh = -u2
        wk = (u1/u2-r)/sqr
        gw = 2.0D+00*gk
        is = 1
    !
    !  (nonzero,0)
    !
    ELSE IF ((u1 /= 0.0D+00) .AND. (u2 == 0.0D+00)) THEN
        b = gh
        wh = -u1
        wk = (u2/u1-r)/sqr
        gw = 2.0D+00*gh
        is = -1
    !
    !  (nonzero,nonzero)
    !
    ELSE IF ((u1 /= 0.0) .AND. (u2 /= 0.0)) THEN
        b = gh+gk
        IF (u1*u2 < 0.0D+00) THEN
            b = b-0.5D+00
        END IF
        wh = -u1
        wk = (u2/u1-r)/sqr
        gw = 2.0D+00*gh
        is = -1
    END IF
    !
    DO
        sgn = -1.0D+00
        t = 0.0D+00
        IF (wk /= 0.0D+00) THEN
            IF (ABS(wk) == 1.0D+00) THEN
                t = wk*gw*(1.0D+00-gw)/2.0D+00
                b = b+sgn*t
            ELSE
                IF (1.0D+00 < ABS(wk)) THEN
                    sgn = -sgn
                    wh = wh*wk
                    g2 = cdf_N01(wh,.FALSE.)
                    wk = 1.0D+00/wk
                    IF (wk < 0.0D+00) THEN
                        b = b+0.5D+00
                    END IF
                    b = b-(gw+g2)/2.0D+00+gw*g2
                END IF
                h2 = wh*wh
                a2 = wk*wk
                h4 = h2/2.0D+00
                ex = EXP(-h4)
                w2 = h4*ex
                ap = 1.0D+00
                s2 = ap-ex
                sp = ap
                s1 = 0.0D+00
                sn = s1
                conex = ABS(con/wk)
                !
                DO
                    cn = ap*s2/(sn+sp)
                    s1 = s1+cn
                    IF (ABS(cn) <= conex) THEN
                        EXIT
                    END IF
                    !
                    sn = sp
                    sp = sp+1.0D+00
                    s2 = s2-w2
                    w2 = w2*h4/sp
                    ap = -ap*a2
                END DO
                t = (ATAN(wk)-wk*s1)/twopi
                b = b+sgn*t
            END IF
        END IF
        !
        IF (0 <= is) THEN
            EXIT
        END IF
        !
        IF (u2 == 0.0D+00) THEN
            EXIT
        END IF
        wh = -u2
        wk = (u1/u2-r)/sqr
        gw = 2.0D+00*gk
        is = 1
    END DO
    !
    b = MAX(b,0.0D+00)
    b = MIN(b,1.0D+00)
    cdf_biv_N01 = b
    RETURN
    !
    END FUNCTION cdf_biv_N01
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    FUNCTION pdf_biv_N01 ( u1, u2, r )
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: u1
    REAL(8), INTENT(IN) :: u2
    REAL(8), INTENT(IN) :: r
    !
    ! Declaring local variables
    !
    REAL(8) :: q
    !
    ! Declaring function's type
    !
    REAL(8) :: pdf_biv_N01
    !
    ! Declaring constants
    !
    REAL(8), PARAMETER :: pi = 3.14159265358979323846264338328d0        
    !
    ! Beginning execution
    !
    q = 1.d0-r**2
    pdf_biv_N01 = EXP((-(u2*(u2/q-(u1*r)/q))-u1*(u1/q-(u2*r)/q))/2.d0)/(2.*pi*SQRT(q))
    !
    ! Ending execution and returning control
    !
    END FUNCTION pdf_biv_N01
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    FUNCTION InverseMillsRatio ( x, u )
!    !
!    ! Calcola 
!    !
!    ! phi(x)/[1-Phi(x)]       se u = .TRUE.
!    ! phi(x)/Phi(x)           se u = .FALSE.
!    !
!    ! NB: C'è un errore nel manuale online dell'Intel Fortran 
!    ! alla voce ERFC_SCALED:
!    ! la frazione 2/sqrt(x) dovrebbe essere 2/sqrt(pi)
!    !
!    REAL(8), INTENT(IN) :: x
!    LOGICAL, INTENT(IN) :: u
!    !
!    REAL(8) :: z
!    !
!    REAL(8) :: InverseMillsRatio
!    !
!    REAL(8), PARAMETER :: invsqrt2 = 0.707106781186547524400844362105d0
!    REAL(8), PARAMETER :: sqrtpiover2 = 1.25331413731550025120788264241d0
!    !
!    IF (u) z = invsqrt2*x
!    IF (.NOT.(u)) z = -invsqrt2*x
!    InverseMillsRatio = 1.d0/(sqrtpiover2*ERFC_SCALED(z))
!    !
!    END FUNCTION InverseMillsRatio
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    SUBROUTINE comp_eigenvals_rsym ( task, need_eigenvectors, n, mat, &
!        mat_eigvect, eigval, flag_error)
!    !
!    ! Computes all eigenvalues and, optionally, eigenvectors of a
!    ! real (n x n) symmetric matrix mat
!    ! Eigenvalues are returned in ascending order
!    !
!    IMPLICIT NONE
!    !
!    ! Declaring dummy variables
!    !
!    CHARACTER(len=30), INTENT(IN) :: task
!    LOGICAL, INTENT(IN) :: need_eigenvectors
!    INTEGER, INTENT(IN) :: n
!    REAL(8), INTENT(IN) :: mat(n,n)
!    REAL(8), INTENT(OUT) :: mat_eigvect(n,n)
!    REAL(8), INTENT(OUT) :: eigval(n)
!    LOGICAL, INTENT(OUT) :: flag_error
!    !
!    ! Declaring local variables
!    !
!    CHARACTER(LEN=1) :: jobz
!    INTEGER :: lwork
!    REAL(8) :: work(10*n)
!    INTEGER :: info
!    !
!    ! Beginning execution
!    !
!    IF (need_eigenvectors) jobz = 'V'
!    IF (.NOT.(need_eigenvectors)) jobz = 'N'
!    mat_eigvect = mat
!    lwork = 10*n
!    CALL dsyev(jobz,'U',n,mat_eigvect,n,eigval,work,lwork,info)
!    IF (info .NE. 0) THEN
!        WRITE(*,1) task
!1       FORMAT('Problem with dsyev in task ', A30)
!        PRINT*, 'info = ', info
!        PAUSE 
!        flag_error = .TRUE.
!        RETURN
!    END IF
!    !
!    ! Ending execution and returning control
!    !
!    END SUBROUTINE comp_eigenvals_rsym
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    SUBROUTINE solve_lin_sys ( task, n, mat1, m, mat2, matsol, flag_error )
!    !
!    ! Computes the solution to a real system of linear equations
!    ! mat1 * X = mat2,
!    ! where mat1 is an (n x n) matrix and X and mat2 are (n x m)matrices.
!    !
!    IMPLICIT NONE
!    !
!    ! Declaring dummy variables
!    !
!    CHARACTER(len=30), INTENT(IN) :: task
!    INTEGER, INTENT(IN) :: n
!    REAL(8), INTENT(IN) :: mat1(n,n)
!    INTEGER, INTENT(IN) :: m
!    REAL(8), INTENT(IN) :: mat2(n,m)
!    REAL(8), INTENT(OUT) :: matsol(n,m)
!    LOGICAL, INTENT(OUT) :: flag_error
!    !
!    ! Declaring local variables
!    !
!    REAL(8) :: tmp(n,n)
!    INTEGER :: ipiv(n), info
!    !
!    ! Beginning execution
!    !
!    tmp = mat1
!    matsol = mat2
!    CALL dgesv(n,m,tmp,n,ipiv,matsol,n,info)
!    IF (info .NE. 0) THEN
!        WRITE(*,1) task
!1       FORMAT('Problem with dgesv in task ', A30)
!        PRINT*, 'info = ', info
!        PAUSE 
!        flag_error = .TRUE.
!        RETURN
!    END IF
!    !
!    ! Ending execution and returning control
!    !
!    END SUBROUTINE solve_lin_sys
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    SUBROUTINE dbrent ( fct, ax, bx, cx, tol, itmax, xmin, fmin )
!    !
!    ! Based on NR77, 10.3, p. 400-1
!    ! Given a function f and its derivative function df, and given 
!    ! a bracketing triplet of abscissas ax, bx, cx [such that bx is 
!    ! between ax and cx, and f(bx) is less than both f(ax) and f(cx)], 
!    ! this routine isolates the minimum to a fractional precision of 
!    ! about tol using a modification of Brent's method that uses 
!    ! derivatives. The abscissa of the minimum is returned as xmin, 
!    ! and the minimum function value is returned as fmin.
!    !
!    ! The objective function and its derivative must be evaluated 
!    ! in a subroutine upon request:
!    !
!    ! SUBROUTINE fct ( c, need_f, f, need_df, df )
!    !
!    ! with LOGICAL need_f and need_df
!    !
!    IMPLICIT NONE
!    !
!    ! Declaring dummy variables
!    !
!    REAL(8), INTENT(IN) :: ax
!    REAL(8), INTENT(IN) :: bx
!    REAL(8), INTENT(IN) :: cx
!    REAL(8), INTENT(IN) :: tol
!    INTEGER, INTENT(IN) :: itmax
!    REAL(8), INTENT(OUT) :: xmin
!    REAL(8), INTENT(OUT) :: fmin
!    !
!    ! Declaring local variables
!    !
!    INTEGER :: iter
!    REAL(8) :: a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw,fx,olde,fa,da
!    REAL(8) :: tol1,tol2,u,u1,u2,v,w,x,xm
!    LOGICAL :: ok1,ok2
!    !
!    ! Declaring parameters
!    !
!    REAL(8), PARAMETER :: zeps = 1.0d-10
!    !
!    ! Beginning execution
!    !
!    a = MIN(ax,cx)
!    b = MAX(ax,cx)
!    v = bx
!    w = v
!    x = v
!    e = 0.d0
!    CALL fct(x,.TRUE.,fx,.TRUE.,dx)
!    fv = fx
!    fw = fx
!    dv = dx
!    dw = dx
!    !
!    DO iter = 1, itmax
!        !
!        xm = 0.5d0*(a+b)
!        tol1 = tol*ABS(x)+zeps
!        tol2 = 2.d0*tol1
!        !
!        IF (ABS(x-xm) .LE. (tol2-0.5d0*(b-a))) GOTO 3
!        IF (ABS(e) .GT. tol1) THEN
!            d1 = 2.d0*(b-a) 
!            d2 = d1
!            IF (dw .NE. dx) d1 = (w-x)*dx/(dx-dw)
!            IF (dv .NE. dx) d2 = (v-x)*dx/(dx-dv) 
!            u1 = x+d1
!            u2 = x+d2
!            ok1 = ((a-u1)*(u1-b) .GT. 0.d0) .AND. (dx*d1 .LE. 0.d0)
!            ok2 = ((a-u2)*(u2-b) .GT. 0.d0) .AND. (dx*d2 .LE. 0.d0)
!            olde = e
!            e = d
!            IF (.NOT.(ok1 .OR.ok2)) THEN 
!                GOTO 1
!            ELSE IF (ok1 .AND. ok2) THEN
!                IF (ABS(d1) .LT. ABS(d2)) THEN
!                    d=d1
!                ELSE
!                    d=d2
!                END IF
!            ELSE IF (ok1) THEN
!                d=d1
!            ELSE
!                d=d2
!            END IF
!            IF (ABS(d) .GT. ABS(0.5d0*olde)) GOTO 1
!            u=x+d
!            IF (((u-a) .LT. tol2) .OR. ((b-u) .LT. tol2)) d = SIGN(tol1,xm-x)
!            GOTO 2
!        END IF
!        !
!1       IF (dx .GE. 0.) THEN 
!            e = a-x
!        ELSE
!            e = b-x
!        END IF
!        d = 0.5d0*e
!        !
!2       IF (ABS(d) .GE. tol1) THEN
!            u = x+d
!            CALL fct(u,.TRUE.,fu,.FALSE.,da)
!        ELSE
!            u = x+SIGN(tol1,d)
!            CALL fct(u,.TRUE.,fu,.FALSE.,da)
!            IF (fu .GT. fx) GOTO 3 
!        END IF 
!        CALL fct(u,.FALSE.,fa,.FALSE.,du)
!        IF (fu .LE. fx) THEN
!            IF (u .GE. x) THEN
!                a = x
!            ELSE
!                b = x
!            END IF
!            v = w
!            fv = fw
!            dv = dw
!            w = x
!            fw = fx
!            dw = dx
!            x = u
!            fx = fu
!            dx = du
!        ELSE
!            IF (u .LT. x) THEN
!                a = u
!            ELSE
!                b = u
!            END IF
!            IF ((fu .LE. fw) .OR. (w .EQ. x)) THEN
!                v = w
!                fv = fw
!                dv = dw
!                w = u
!                fw = fu
!                dw = du
!            ELSE IF ((fu .LE. fv) .OR. (v .EQ. x) .OR. (v .EQ. w)) THEN
!                v = u
!                fv = fu
!                dv = du
!            END IF
!        END IF
!    END DO
!    !
!    PAUSE 'dbrent exceeded maximum iterations'
!    !
!3   xmin = x
!    fmin = fx
!    RETURN
!    !
!    ! Ending subroutine and returning control
!    !
!    END SUBROUTINE dbrent
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    SUBROUTINE mnbrak ( ax, bx, cx, fa, fb, fc, fct )
!    !
!    ! Based on NR77, 10.1, p. 393-4
!    ! Given a function func, and given distinct initial points ax and bx, 
!    ! this routine searches in the downhill direction (defined by the function 
!    ! as evaluated at the initial points) and returns new points ax, bx, cx 
!    ! that bracket a minimum of the function. Also returned are 
!    ! the function values at the three points, fa, fb, and fc.
!    ! Parameters: GOLD is the default ratio by which successive intervals 
!    ! are magnified; GLIMIT is the maximum magnification allowed 
!    ! for a parabolic-fit step.
!    !
!    ! The objective function and its derivative must be evaluated 
!    ! in a subroutine upon request:
!    !
!    ! SUBROUTINE fct ( c, need_f, f, need_df, df )
!    !
!    ! with LOGICAL need_f and need_df
!    !
!    IMPLICIT NONE
!    !
!    ! Declaring dummy arguments
!    !
!    REAL(8), INTENT(INOUT) :: ax
!    REAL(8), INTENT(INOUT) :: bx
!    REAL(8), INTENT(OUT) :: cx
!    REAL(8), INTENT(OUT) :: fa
!    REAL(8), INTENT(OUT) :: fb
!    REAL(8), INTENT(OUT) :: fc
!    !
!    ! Declaring parameters
!    !
!    REAL(8), PARAMETER :: GOLD = 1.618034d0
!    REAL(8), PARAMETER :: GLIMIT=100.d0
!    REAL(8), PARAMETER :: TINY=1.d-20
!    !
!    ! Declaring local variables
!    !
!    REAL(8) :: dum,fu,q,r,u,ulim,da
!    !
!    ! Beginning execution
!    !
!    CALL fct(ax,.TRUE.,fa,.FALSE.,da)
!    CALL fct(bx,.TRUE.,fb,.FALSE.,da)
!    IF (fb .GT. fa) THEN
!        dum = ax 
!        ax = bx
!        bx = dum
!        dum = fb
!        fb = fa
!        fa = dum
!    END IF
!    !
!    cx = bx+GOLD*(bx-ax) 
!    CALL fct(cx,.TRUE.,fc,.FALSE.,da)
!    !
!1   IF (fb .GE. fc) THEN
!        r = (bx-ax)*(fb-fc) 
!        q = (bx-cx)*(fb-fa)
!        u = bx-((bx-cx)*q-(bx-ax)*r)/(2.*SIGN(MAX(ABS(q-r),TINY),q-r))
!        ulim = bx+GLIMIT*(cx-bx)
!        IF ((bx-u)*(u-cx) .GT. 0.d0) THEN
!            CALL fct(u,.TRUE.,fu,.FALSE.,da)
!            IF (fu .LT. fc) THEN
!                ax = bx
!                fa = fb
!                bx = u
!                fb = fu
!                RETURN
!            ELSE IF (fu .GT. fb) THEN 
!                cx = u
!                fc = fu
!                RETURN
!            END IF
!            !
!            u = cx+GOLD*(cx-bx) 
!            CALL fct(u,.TRUE.,fu,.FALSE.,da)
!        ELSE IF ((cx-u)*(u-ulim) .GT. 0.d0) THEN 
!            CALL fct(u,.TRUE.,fu,.FALSE.,da)
!            IF (fu .LT. fc) THEN
!                bx = cx
!                cx = u
!                u = cx+GOLD*(cx-bx)
!                fb = fc
!                fc = fu
!                CALL fct(u,.TRUE.,fu,.FALSE.,da)
!            END IF
!        ELSE IF ((u-ulim)*(ulim-cx) .GE. 0.d0) THEN
!            u = ulim
!            CALL fct(u,.TRUE.,fu,.FALSE.,da)
!        ELSE
!            u = cx+GOLD*(cx-bx)
!            CALL fct(u,.TRUE.,fu,.FALSE.,da)
!        END IF
!        !
!        ax = bx
!        bx = cx
!        cx = u
!        fa = fb
!        fb = fc
!        fc = fu
!        GOTO 1
!    END IF
!    !
!    RETURN
!    !
!    ! Ending execution and returning control
!    !
!    END SUBROUTINE mnbrak
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
END MODULE sp
