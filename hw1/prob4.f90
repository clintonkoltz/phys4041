!#################################################################################
!#
!#  Subroutine: single_epsilon
!#
!#  Purpose: find machine epsilon for single precision.
!#           evaluates to 2.98023224E-08 on my machine,
!#           this means there are ~7 digits of precision.
!#################################################################################
SUBROUTINE single_epsilon

    IMPLICIT NONE

    REAL*4 x,y,z,eps

    eps = 1.
    x = 1. + eps
    y = 1. - eps
    z = x - y

    DO WHILE ( z/2. .GT. 0. )
        eps = eps/2.

        x = 1. + eps
        y = 1. - eps
        z = x - y

    END DO

    PRINT *, 'machine precision for 4-byte floating point:'
    PRINT *, eps

END SUBROUTINE single_epsilon
!#################################################################################
!#
!#  Subroutine: double_epsilon
!#
!#  Purpose: find machine epsilon for double precision.
!#           evaluates to 5.55111512312578270E-017 on my
!#           machine, this means there are ~16 digits
!#           of precision.
!#################################################################################
SUBROUTINE double_epsilon

    IMPLICIT NONE

    REAL*8 x,y,z,eps

    eps = 1.d0
    x = 1.d0 + eps
    y = 1.d0 - eps
    z = x - y

    DO WHILE ( z/2.d0 .GT. 0.d0 )
        eps = eps/2.d0

        x = 1.d0 + eps
        y = 1.d0 - eps
        z = x - y

    END DO

    PRINT *, 'machine precision for 8-byte floating point:'
    PRINT *, eps

END SUBROUTINE double_epsilon
!#################################################################################
!#
!#  Main Program: test routines for determining digits of precision
!#
!#################################################################################
PROGRAM problem4_main

    CALL single_epsilon()
    CALL double_epsilon()

END PROGRAM problem4_main




!#################################################################################
!#################################################################################
