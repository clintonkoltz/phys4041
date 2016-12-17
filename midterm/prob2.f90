!#################################################################################
!#  MODULE: mod_pc_int - integrate simple ODE dy/dt = f(y) using a predictor-
!#      corrector method, called by pre_cor
!#
!#################################################################################
MODULE mod_pc_int
IMPLICIT NONE

    REAL*8 :: t = 0.d0          ! t - ind. var, initialized to 0
    REAL*8 :: dt = 0.1d0        ! dt - stepsize, initialized to 0.1
    REAL*8 :: y_i = 1.d0        ! y_i - current value of y
    REAL*8 :: tol = 1.d-6       ! tol - tolerance for corrector

CONTAINS
!#################################################################################
!#  pre_cor - calculates y(t_i+1) given y(t_i) and user-defined function dydt
!#
!#################################################################################
FUNCTION pre_cor(y_now) RESULT(y_new)
    REAL*8, INTENT(IN) :: y_now
    REAL*8 :: y_new

!### y_i+1 = corrected value output from euler method on y_i
    y_new = corrector(euler(y_now))
!### advance timestep
    t = t + dt

END FUNCTION
!#################################################################################
!#  euler - calculates prediction value 'y_pre' from simple first-order
!#      forward euler method
!#
!#################################################################################
FUNCTION euler(y_now) RESULT(y_pre)
    REAL*8, INTENT(IN) :: y_now
    REAL*8 :: y_pre

    y_pre = y_now + dt * dydt(y_now,t)

END FUNCTION euler
!#################################################################################
!#  corrector - iteratively calculate y_i+1 assuming the finite difference 
!#      [delta y / delta t] equals the average of the analytic function
!#      evaluated at each point [= 1/2 (dydt(y_i,t_i) + dydt(y_i+1,t_i+1)].
!#      repeat until absolute difference between new guess and previous guess
!#      is less than tolerance
!#################################################################################
FUNCTION corrector(y_pre) RESULT(y_new)
    INTEGER, PARAMETER :: imax = 10000000
    REAL*8 :: y_pre     !   y_pre - previous guess
    REAL*8 :: y_new     !   y_new - new guess
    INTEGER :: i 

    i = 0
    !### for first iteration, new guess is value from euler method,
    !   previous guess is fixed point, y_i
    y_new = y_pre;  y_pre = y_i
    !### trapezoidal rule estimation of dydt, iterate to solve implicit equation
    DO WHILE(abs(y_pre-y_new).gt.tol)
       i = i+1
       y_pre = y_new
       y_new = y_i + 0.5d0*dt*(dydt(y_i,t) + dydt(y_pre,t+dt))
       !### check against max iterations, end if convergence not found
       IF (i.eq. imax) THEN
          PRINT '("NO CONVERGENCE",4f14.8)',dt ,y_i, y_pre, y_new; STOP
       END IF
    END DO
END FUNCTION corrector
!#################################################################################
!#  dydt - analytic function to be integrated. y -> Pi for y(0) = 1 as t -> infty
!#    
!#################################################################################
FUNCTION dydt(y,t)
    REAL*8, INTENT(IN) :: y,t
    REAL*8 :: dydt

    dydt = sin(y) / sqrt(y)

END FUNCTION dydt
END MODULE mod_pc_int
!#################################################################################
!#  Main Program: Integrate system out to t = 1, then further to demonstrate
!#      convergence of the system to Pi 
!#################################################################################
PROGRAM prob2
USE mod_pc_int
IMPLICIT NONE

    REAL*8 :: tfinal = 1.d0
    REAL*8 :: tmore = 30.d0

    OPEN( UNIT=10, FILE= 'pre_cor.dat')
    OPEN( UNIT=20, FILE= 'long.dat')
    OPEN( UNIT=30, FILE= 'stab.dat')
    OPEN( UNIT=40, FILE= 'instab.dat')
200 FORMAT(2(1pe28.15))
!### write out data to t = 1
    DO WHILE (t.lt.tfinal)
        WRITE(10,200) t,y_i
        WRITE(20,200) t,y_i
        y_i = pre_cor(y_i)
    END DO
!### demonstrate how solution goes to Pi for large t
    DO WHILE (t.lt.tmore)
        WRITE(20,200) t,y_i
        y_i = pre_cor(y_i)
    END DO
!### demonstrate how this step size still converges...
    t = 0.d0;   dt = 3.337d0
    y_i = 1.d0
    DO WHILE (t.lt.tmore)
       WRITE(30,200) t,y_i
       y_i = pre_cor(y_i)
    END DO
!### ...but this one doesn't 
    t = 0.d0;  dt = 3.338d0
    y_i = 1.d0
    DO WHILE (t.lt.tmore)
       WRITE(40,200) t,y_i
       y_i = pre_cor(y_i)
    END DO
END PROGRAM
