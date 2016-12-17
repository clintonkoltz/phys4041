!###################################################################################
!#  Module: mod_spline - contains types and procedures for generating a cubic spline
!#      through an arbitrary number of data points (x_i,y_i), using natural spline
!#      for boundary condition closure.
!#
!###################################################################################
MODULE mod_spline

IMPLICIT NONE
!### spline_gen - spline generator object has storage for x,y data pairs and 
!    coefficients of a_i + b_i*t + c_i*t^2 + d_i*t^3 in each interval
TYPE :: spline_gen

    REAL*8, DIMENSION(:), ALLOCATABLE :: x_data, y_data
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: coeffs
    INTEGER :: n_pts !### n_pts - # of data points

END TYPE spline_gen

CONTAINS
!###################################################################################
!#  build - allocates space for spline generator object
!#
!###################################################################################
SUBROUTINE build(dat)

    TYPE(spline_gen) :: dat

    ALLOCATE( dat%x_data(dat%n_pts), dat%y_data(dat%n_pts), dat%coeffs(4, dat%n_pts) )

END SUBROUTINE build
!###################################################################################
!#  compute_spline_coeffs - solves tridiagonal matrix equation A*D = b by
!#  LU decomposition, uses D_i = dy/dx at x_i to compute cubic polynomial 
!#  coefficients in  each interval between data points.
!#
!###################################################################################
SUBROUTINE compute_spline_coeffs(dat)

    IMPLICIT NONE

    TYPE(spline_gen) :: dat
    REAL*8, DIMENSION(dat%n_pts) :: alpha, beta, b, D, y
    INTEGER :: i

    !### b vector as defined for natural splines
    b(1) = 3.d0*( dat%y_data(2) - dat%y_data(1) )
    DO i = 2, dat%n_pts - 1
       b(i) = 3.d0*( dat%y_data(i+1) - dat%y_data(i-1) )
    END DO
    b(dat%n_pts) = 3.d0*( dat%y_data(dat%n_pts) - dat%y_data(dat%n_pts-1) )


    DO i = 1,6
        PRINT *,i,b(i)
    ENDDO

    !### based on Crout's algorithm for LU decomp, see pdf for explicit decomposition
    alpha(1) = 2.d0
    DO i = 2, dat%n_pts - 1
       beta(i) = 1.d0 / alpha(i-1)
       alpha(i) = 4.d0 - beta(i)
    END DO
    beta(dat%n_pts) = 1.d0 / alpha(dat%n_pts-1)
    alpha(dat%n_pts) = 2.d0 - beta(dat%n_pts)

    !### solve L*y = b for y by forward substitution
    y(1) = b(1)
    DO i = 2, dat%n_pts
       y(i) = b(i) - beta(i)*y(i-1)
    END DO

    !### solve U*x = y for x by backward substituion
    D(dat%n_pts) = y(dat%n_pts) / alpha(dat%n_pts)
    print *,d(dat%n_pts)
    DO i = dat%n_pts-1, 1, -1
       D(i) = (y(i) - D(i+1)) / alpha(i)
       print *,i,d(i)
    END DO



    !### generate 4 x (n-1) matrix of cubic polynomial coefficients. for n points,
    !    there are n-1 intervals. coeffs(:,i) are a,b,c,d coefficients in cubic 
    !    polynomial equation, which is portion of the spline in between points i & i+1
    DO i = 1, dat%n_pts-1
       dat%coeffs(1,i) = dat%y_data(i)
       dat%coeffs(2,i) = D(i)
       dat%coeffs(3,i) = 3.d0*( dat%y_data(i+1) - dat%y_data(i) ) - 2.d0*D(i) - D(i+1)
       dat%coeffs(4,i) = 2.d0*( dat%y_data(i) - dat%y_data(i+1) ) + D(i) + D(i+1)
    END DO

END SUBROUTINE compute_spline_coeffs
!###################################################################################
!#  spline_interp - returns y-value of spline given x and spline_gen object
!#
!###################################################################################
FUNCTION spline_interp(dat, x_int) result(y_int)

    TYPE(spline_gen) :: dat
    REAL*8 :: x_int
    REAL*8 :: x_norm
    REAL*8 :: y_int
    INTEGER :: bin

    !### find which interval x_int is in
    DO bin = 1, dat%n_pts - 1
       IF (x_int .gt. dat%x_data(bin) .and. x_int .le. dat%x_data(bin+1)) EXIT
       IF (x_int .eq. dat%x_data(bin) ) EXIT
    END DO

    IF (bin .ge. dat%n_pts) THEN
       PRINT *,'x value outside x range',x_int
       STOP
    END IF

    !### polynomial is defined y = a_i + b_i*t + c_i*t^2 + d_i*t^3, where t
    !    runs from 0 to 1 over the interval, normalize x_int accordingly
    x_norm = (x_int - dat%x_data(bin)) / (dat%x_data(bin+1) - dat%x_data(bin))

    y_int = dat%coeffs(1,bin) + x_norm*dat%coeffs(2,bin) + &
            x_norm**2*dat%coeffs(3,bin) + x_norm**3*dat%coeffs(4,bin)

    RETURN

END FUNCTION spline_interp
!###################################################################################
!#  bisection_meth - using bisection root finder to find x-value (root) where spline 
!#      function returns yval
!###################################################################################
SUBROUTINE bisection_meth(root,dat,yval,bounds)

    REAL*8, PARAMETER :: tol = 1.d-8
    REAL*8, PARAMETER :: imax = 1000
    TYPE(spline_gen),INTENT(IN) :: dat
    REAL*8, INTENT(OUT) :: root
    REAL*8, INTENT(IN) :: yval
    REAL*8, DIMENSION(2) :: sgns = (/1.d0,1.d0/)
    REAL*8, DIMENSION(2),INTENT(INOUT) :: bounds
    REAL*8, DIMENSION(2) :: y_bounds
    REAL*8 :: testpt, yroot
    INTEGER :: i

    y_bounds = (/ spline_interp(dat,bounds(1)) - yval, &
                 spline_interp(dat,bounds(2)) - yval /)
    sgns = sign(sgns,y_bounds)

    yroot = 1.d0
    i = 0

    DO WHILE (abs(yroot) .ge. tol)
       testpt = 0.5d0*(bounds(1) + bounds(2))
       yroot = spline_interp(dat,testpt) - yval
       IF (sign(sgns(1),yroot) .eq. sgns(1)) THEN
          bounds(1) = testpt
       ELSE IF (sign(sgns(2),yroot) .eq. sgns(2) ) THEN 
          bounds(2) = testpt
       END IF
       i = i + 1
    END DO

    root = testpt

END SUBROUTINE bisection_meth

END MODULE mod_spline
!###################################################################################
!#
!#  Main Program: problem2_main tests routines in module mod_spline, 
!#      provides solutions for problem 2
!#
!###################################################################################
PROGRAM problem2_main

USE mod_spline

IMPLICIT NONE

    TYPE(spline_gen) :: myprob
    REAL*8 :: x, y_int
    REAL*8 :: x_root, dx
    REAL*8, DIMENSION(2) :: init_bnd
    INTEGER :: interp_pts
    INTEGER :: i

    OPEN(unit=10, file='exper.dat')
    OPEN(unit=20, file='interp.dat')

!### this format writes out two floating points in exponential format with
!    12 digits following the decimal point
100 FORMAT(2(1pe24.12))

    myprob%n_pts = 6

    !### allocate myprob data structure and input values for P(v)
    CALL build(myprob)
    myprob%x_data = (/ 1.d0, 1.1d0, 1.2d0, 1.3d0, 1.4d0, 1.5d0 /)
    myprob%y_data = (/ 4.7d0, 3.5d0, 3.d0, 2.7d0, 3.d0, 2.4d0 /)

    !### write out data values for plotting
    DO i = 1,myprob%n_pts
       WRITE(10,100) myprob%x_data(i), myprob%y_data(i)
    END DO

    !### calculate polynomials in each interval
    CALL compute_spline_coeffs(myprob)

    !### interpolate spline throughout interval, write out to file for plotting
    interp_pts = 200
    dx = (myprob%x_data(myprob%n_pts)-myprob%x_data(1))/ dble(interp_pts)
    DO i = 0, interp_pts
       x = 1.d0 + dx * dble(i)
       y_int = spline_interp(myprob,x)
       WRITE(20,100) x, y_int
    END DO

    !### from plot, P = 3.25 appears to occur around v = 1.13, use bisection method
    init_bnd = (/ 1.11d0, 1.15d0 /)
    CALL bisection_meth(x_root,myprob,3.25d0,init_bnd)

    PRINT '(2f13.8)', x_root, spline_interp(myprob,x_root)


END PROGRAM
