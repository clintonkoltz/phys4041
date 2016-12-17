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
    REAL*8, DIMENSION(:), ALLOCATABLE :: x_data, y_data, D
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: coeffs
    INTEGER :: n_pts !### n_pts - # of data points
    CONTAINS
    PROCEDURE :: build_spline
    PROCEDURE :: destroy_spline
END TYPE spline_gen

CONTAINS
!###################################################################################
!#  build - allocates space for spline generator object
!#  destroy - deallocates memory
!###################################################################################
SUBROUTINE build_spline(dat,n_pts)
    CLASS(spline_gen) :: dat
    INTEGER, INTENT(IN) :: n_pts

    dat%n_pts = n_pts
    ALLOCATE( dat%x_data(dat%n_pts), dat%y_data(dat%n_pts), dat%D(dat%n_pts),dat%coeffs(4, dat%n_pts) )
END SUBROUTINE build_spline
SUBROUTINE destroy_spline(dat)
    CLASS(spline_gen) :: dat

    DEALLOCATE( dat%x_data,dat%y_data,dat%coeffs )
END SUBROUTINE destroy_spline
!###################################################################################
!#  compute_spline_coeffs - solves tridiagonal matrix equation A*D = b by
!#      LU decomposition, uses D_i = dy/dx at x_i to compute cubic polynomial 
!#      coefficients in  each interval between data points.
!#
!###################################################################################
SUBROUTINE compute_spline_coeffs(dat)
    TYPE(spline_gen) :: dat
    REAL*8, DIMENSION(dat%n_pts) :: alpha, beta, b, y
    INTEGER :: i

    !### b vector as defined for natural splines
    b(1) = 3.d0*( dat%y_data(2) - dat%y_data(1) )
    DO i = 2, dat%n_pts - 1
       b(i) = 3.d0*( dat%y_data(i+1) - dat%y_data(i-1) )
    END DO
    b(dat%n_pts) = 3.d0*( dat%y_data(dat%n_pts) - dat%y_data(dat%n_pts-1) )

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
    dat%D(dat%n_pts) = y(dat%n_pts) / alpha(dat%n_pts)
    DO i = dat%n_pts-1, 1, -1
       dat%D(i) = (y(i) - dat%D(i+1)) / alpha(i)
    END DO
    !PRINT *, dat%D
    !### generate 4 x (n-1) matrix of cubic polynomial coefficients. for n points,
    !    there are n-1 intervals. coeffs(:,i) are a,b,c,d coefficients in cubic 
    !    polynomial equation, which is portion of the spline in between points i & i+1
    DO i = 1, dat%n_pts-1
       dat%coeffs(1,i) = dat%y_data(i)
       dat%coeffs(2,i) = dat%D(i)
       dat%coeffs(3,i) = 3.d0*( dat%y_data(i+1) - dat%y_data(i) ) - 2.d0*dat%D(i) - dat%D(i+1)
       dat%coeffs(4,i) = 2.d0*( dat%y_data(i) - dat%y_data(i+1) ) + dat%D(i) + dat%D(i+1)
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

    bin = locate_bin(dat,x_int)
    !### polynomial is defined y = a_i + b_i*t + c_i*t^2 + d_i*t^3, where t
    !    runs from 0 to 1 over the interval, normalize x_int accordingly
    x_norm = (x_int - dat%x_data(bin)) / (dat%x_data(bin+1) - dat%x_data(bin))

    y_int = dat%coeffs(1,bin) + x_norm*dat%coeffs(2,bin) + &
            x_norm**2*dat%coeffs(3,bin) + x_norm**3*dat%coeffs(4,bin)

END FUNCTION spline_interp
!###################################################################################
!#  locate_bin - returns integer bin # for given input x value
!#
!###################################################################################
FUNCTION locate_bin(dat,x_int) result(bin_x)
    TYPE(spline_gen),INTENT(IN) :: dat
    REAL*8,INTENT(IN) :: x_int
    INTEGER :: bin_x
    INTEGER :: bin

    DO bin = 1, dat%n_pts - 1
       IF (x_int .gt. dat%x_data(bin) .and. x_int .le. dat%x_data(bin+1)) EXIT
       IF (x_int .eq. dat%x_data(bin) ) EXIT
    END DO
    IF (bin .ge. dat%n_pts) THEN
       PRINT *,'x value outside x range',x_int
       STOP
    END IF

    bin_x = bin

END FUNCTION locate_bin
END MODULE mod_spline
