!##############################################################################
!#  uses open source code to generate spline, confirm P(v) = 3.25
!##############################################################################
PROGRAM burkardt_spline
IMPLICIT NONE

    INTEGER, PARAMETER :: N = 6
    REAL*8 :: t(N) = (/ 1.d0, 1.1d0, 1.2d0, 1.3d0, 1.4d0, 1.5d0 /)
    REAL*8 :: y(N) = (/ 4.7d0, 3.5d0, 3.d0, 2.7d0, 3.d0, 2.4d0 /)
    INTEGER :: ibcbeg = 2, ibcend = 2   !### boundary condition ids reqd by spline.f90
    REAL*8 :: ybcbeg = 0.d0, ybcend = 0.d0 !### natural splines
    REAL*8 :: ypp(N) !### 2nd derivs output by spline_cubic_set
    REAL*8 :: tval = 1.13923d0
    REAL*8 :: yval ,ypval ,yppval


    CALL spline_cubic_set( N, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp )

    CALL spline_cubic_val( N, t, y, ypp, tval, yval, ypval, yppval )

    PRINT *, yval

END PROGRAM
