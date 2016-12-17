!#################################################################################
!#  MODULE: mod_finite_diff - functions for approximating f_prime(x) of 
!#      user-defined function f(x) by finite differences. 1st-order forward
!#      difference from x to x + dx. 2nd-order central difference approximating
!#      f_prime(x) using values of f(x) at x-dx & x+dx
!#################################################################################
MODULE mod_finite_diff
IMPLICIT NONE

CONTAINS
REAL*8 FUNCTION forward_diff(x_val,dx)
    REAL*8, INTENT(IN) :: x_val
    REAL*8, INTENT(IN) :: dx

    forward_diff = (f(x_val + dx) - f(x_val) ) / dx
END FUNCTION forward_diff

REAL*8 FUNCTION central_diff(x_val,dx)
    REAL*8, INTENT(IN) :: x_val
    REAL*8, INTENT(IN) :: dx

    central_diff = (f(x_val + dx) - f(x_val - dx) ) / (2.d0*dx)
END FUNCTION central_diff

REAL*8 FUNCTION f(x)
    REAL*8, INTENT(IN) :: x
    f = tan(x)
END FUNCTION

REAL*8 FUNCTION f_prime(x)
    REAL*8, INTENT(IN) :: x
    f_prime = 1.d0 / cos(x)**2
END FUNCTION

END MODULE mod_finite_diff

!#################################################################################
!#  Main Program: approximates derivative f(x) = tan(x) over interval [-2,2] at 
!#      n points for [100,5000]. excludes areas around infinities at -pi/2 and pi/2
!#      using buffer zone. demonstrates different rates of convergence for each method
!#################################################################################
PROGRAM prob1
USE mod_finite_diff

IMPLICIT NONE
    REAL*8, PARAMETER :: pihalf = 2.d0*atan(1.d0)
    REAL*8 :: x,dx
    REAL*8 :: for_err, cent_err, buffer
    REAL*8 :: first_order, second_order, analytic
    REAL*8,DIMENSION(2) :: bounds = (/ -2.d0, 2.d0 /)
    INTEGER :: n,i

!### open data file to write out n & maximum error for each method
    OPEN(UNIT=10, FILE='diff.dat')
100 FORMAT(I6,2(1pe20.10))
!### loop over n points from 100 to 5000 by 100
    DO n = 100, 5000, 100
!### define dx for n points and rezero errors
        dx = (bounds(2) - bounds(1)) / (n-1)
        for_err = 0.d0 
        cent_err = 0.d0
!### loop over range, finding maximum error for each method
        DO i = 1, n
           x = bounds(1) + (i-1)*dx
           analytic = f_prime(x)
           buffer = 0.1d0
           IF ( .NOT.((x.gt.-pihalf-buffer).AND.(x.lt.-pihalf+buffer)) .AND. &
               .NOT.((x.gt. pihalf-buffer).AND.(x.lt. pihalf+buffer)) ) THEN
              first_order = forward_diff(x,dx)
              for_err = max(for_err, &
                 abs((first_order - analytic)/ analytic ))
              second_order = central_diff(x,dx)
              cent_err = max(cent_err, &
                 abs((second_order - analytic)/ analytic ))
           END IF
        END DO
!### write out data
        WRITE(10,100) n, for_err, cent_err
    END DO
END PROGRAM
