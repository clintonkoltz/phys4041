!#################################################################################
!#  MODULE: mod_trap_rule - performs trapezoidal rule to approximate
!#      integral of user-defined functions
!#
!#################################################################################
MODULE mod_trap_rule
IMPLICIT NONE

ABSTRACT INTERFACE
    FUNCTION func(x)
       REAL*8,INTENT(IN) :: x
       REAL*8 :: func
    END FUNCTION
END INTERFACE

CONTAINS
!#################################################################################
!# trap_rule - integrates input function f(x) from limits(1) to limits(2)
!#      using N-1 trapezoids
!#
!#################################################################################
REAL*8 FUNCTION trap_rule(limits,N,f)
    REAL*8, DIMENSION(2), INTENT(IN) :: limits
    INTEGER, INTENT(IN) :: N        !### N - # of points used in sum
    PROCEDURE(func) :: f            !### f is a procedure with in/out similar to func
    REAL*8 :: dx,x
    INTEGER :: i

!### define width of trapezoids
    dx = (limits(2)-limits(1))/dble(N-1)

!###  sum [ dx/2*( f_i + f_{i+1}) ]
!### exterior points have weighting 1/2, interior points have weighting 1
    trap_rule = 0.5d0*f(limits(1))
    DO i = 2,N-1
       x = limits(1) + dx * dble(i-1)
       trap_rule = trap_rule + f(x)
    END DO
    trap_rule = trap_rule + 0.5d0*f(limits(2))

    trap_rule = trap_rule * dx

END FUNCTION trap_rule
END MODULE mod_trap_rule
!#################################################################################
!#  MODULE: mod_prob2 - user defined functions to be integrated, see hw4.pdf
!#
!#################################################################################
MODULE mod_prob2
IMPLICIT NONE

    REAL*8 s        !### s - ind. parameter varies from 0.5 - 3.0

CONTAINS
    FUNCTION f1(y)
       REAL*8 :: f1
       REAL*8, INTENT(IN) :: y
       f1 = exp(y) / sqrt( exp(3.d0*y) + s*exp(y) )
    END FUNCTION f1
    FUNCTION f2(x)
       REAL*8 :: f2
       REAL*8, INTENT(IN) :: x
       f2 = exp(-s*x) / sqrt(x**2 + 1.d0)
    END FUNCTION f2
END MODULE mod_prob2
!#################################################################################
!#  Main program: Uses trapezoidal rule quadrature to approximate integrals in
!#      mod_prob2 for s over range [0.5,3.0] using approximations and
!#      transformations in hw4.pdf.
!#################################################################################
PROGRAM prob2
USE mod_trap_rule
USE mod_prob2

IMPLICIT NONE

    REAL*8 :: x0,x1
    REAL*8 :: partA,partB
    INTEGER :: i

100 FORMAT(3(a20))
200 FORMAT(3(f20.11))
    PRINT 100, 's','part A', 'part B'
!### step through range in s by 0.1
    DO i = 0,25
       s = 0.5d0 + dble(i)/10.d0
       
!### integral in part a.) is divided into three sections:
!      0 to x0 && x0 to x1 && x1 to Infty
       x0 = sqrt(s) * 1.d-5
       x1 = sqrt(s) * 1.d5
       partA = 2.d0 * sqrt(x0/s)
       partA = partA + trap_rule( (/log(x0),log(x1)/), 10000, f1)
       partA = partA + 2.d0/sqrt(x1)       

!### integral in part b.) truncated where exponential goes to 0
       x0 = 10.d0*log(10.d0)/s
       partB = trap_rule( (/0.d0,x0/), 100000, f2)

!### mathematica say:
!   part a.) = (8 *Gamma[5/4]^2)/(Sqrt[Pi] s^(1/4))
!   part b.) =  Pi/2 *(-BesselY[0, s] + StruveH[0, s])
!   these agree to ~8 digits
       PRINT 200, s, partA, partB


    END DO

END PROGRAM
