!#################################################################################
!#
!#  Subroutines: formula1, formula2, formula3
!#          perform floating point summations for each nominally equivalent formula,
!#          demonstrates loss of precision from underflow for third formula
!#          true answer =  Ln( e/2 ) = 0.30685281944005...
!#
!#################################################################################
SUBROUTINE formula1(N,sum_N)

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: N
    REAL*4,INTENT(OUT) :: sum_N
    INTEGER :: ii

    sum_N = 0.

    DO ii = 1, 2*N

        sum_N = sum_N + (-1)**ii * ii /(ii+1.)

    END DO

END SUBROUTINE formula1
!#################################################################################
SUBROUTINE formula2(N,sum_N)

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: N
    REAL*4,INTENT(OUT) :: sum_N
    INTEGER :: ii

    sum_N = 0.

    DO ii = 1, N

        sum_N = sum_N - (2.*ii - 1.)/(2.*ii) + 2.*ii/(2.*ii + 1 )

    END DO

END SUBROUTINE formula2
!#################################################################################
SUBROUTINE formula3(N,sum_N)

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: N
    REAL*4,INTENT(OUT) :: sum_N
    INTEGER :: ii

    sum_N = 0.

    DO ii = 1, N

        sum_N = sum_N + 1./ (2.*dble(ii) * (2.*dble(ii)+1.))

    END DO 

END SUBROUTINE formula3
!#################################################################################
!### formula3d - demonstrates better convergence for 3rd formula with double precision
!#################################################################################
SUBROUTINE formula3d(N,sum_N)

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: N
    REAL*8,INTENT(OUT) :: sum_N
    INTEGER :: ii

    sum_N = 0.d0

    DO ii = 1, N

        sum_N = sum_N + 1.d0/ (2.d0*dble(ii) * (2.d0*dble(ii)+1.d0))

    END DO

END SUBROUTINE formula3d
!#################################################################################
!### formula3_alt - work around to account for underflow loss of precision
!#################################################################################
SUBROUTINE formula3_alt(N,sum_N)

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: N
    REAL*4,INTENT(OUT) :: sum_N
    REAL*4 :: small_N,smaller_N,tiny_N
    INTEGER :: ii

    sum_N = 0.
    small_N = 0.
    smaller_N = 0.
    tiny_N = 0.

    DO ii = 1, 1000

        sum_N = sum_N + 1./ (2.*dble(ii) * (2.*dble(ii)+1.))

    END DO
    PRINT *,''
    PRINT *, 'big numbers:', sum_N

    DO  ii = 1001,10000

        small_N = small_N + 1. / (2.*dble(ii) * (2.*dble(ii)+1.))

    END DO

    PRINT *, 'small numbers:', small_N

    DO ii = 10001, 100000

        smaller_N = smaller_N + 1. / (2.*dble(ii) * (2.*dble(ii)+1.))

    END DO

    PRINT *, 'smaller numbers:', smaller_N

    DO ii = 100001, N

        tiny_N = tiny_N + 1. / (2.*dble(ii) * (2.*dble(ii)+1.))

    END DO

    PRINT *, 'tiny numbers:', tiny_N
    PRINT *,''

    smaller_N = smaller_N + tiny_N

    small_N = small_N + smaller_N

    sum_N = small_N + sum_N

END SUBROUTINE formula3_alt
!#################################################################################
!#
!#  Main Program: test sumation formula routines
!#
!#################################################################################
PROGRAM problem5_main

    IMPLICIT NONE

    REAL*4 :: var1,var2,var3
    INTEGER :: ii,j
    INTEGER,DIMENSION(3) :: list = (/100,1000,10000/)


    PRINT *,'    N        formula1        formula2         formula3'

    !### test each formula at the requested values of N
    DO j = 1,3
        ii = list(j)
        CALL formula1(ii,var1)
        CALL formula2(ii,var2)
        CALL formula3(ii,var3)
        PRINT *, ii,var1,var2,var3
    END DO

    !### demonstrates how formula 3 fails to converge to the right answer in
    !    in single precision.
    ii = 2500000
    CALL formula3(ii,var1)
    CALL formula3_alt(ii,var2)
    CALL formula3d(ii,var3)

    print *, '  N        4-byte formula3  4-byte formula3_alt  8-byte_formula3'
    print *, ii, var1, var2, var3

END PROGRAM problem5_main

!#################################################################################
