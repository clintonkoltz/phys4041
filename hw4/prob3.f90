!#################################################################################
!#  MODULE: mod_rk4 - data types and procedures to perform 4th-order Runge-Kutta
!#      integration of ODE. requires user-defined functional form for dfdt
!#
!#################################################################################
MODULE mod_rk4
IMPLICIT NONE

TYPE :: rk4
    REAL*8, DIMENSION(:), ALLOCATABLE :: f  ! f - array of dependent variables
    REAL*8 :: t,dt                          ! t - independent variable, dt - stepsize
    INTEGER :: nvar                         ! nvar - # of dependent variables

    CONTAINS
    PROCEDURE :: init
    PROCEDURE :: rk4_step
    PROCEDURE :: dfdt       !### dummy function for derivatives of f w.r.t. t
END TYPE rk4

CONTAINS
!#################################################################################
!#  init - initializes ODE problem for rk4 solver. allocates dependent variable
!#      array. stores initial data in structure
!#
!#################################################################################
SUBROUTINE init(self,nvar,t0,dt,f0)

    CLASS(rk4) :: self
    REAL*8, INTENT(IN) :: t0,dt
    INTEGER, INTENT(IN) :: nvar
    REAL*8, INTENT(IN), DIMENSION(nvar) :: f0

    ALLOCATE(self%f(nvar))
    self%nvar = nvar;   self%t = t0;    self%dt = dt;    self%f = f0
END SUBROUTINE init
!#################################################################################
!#  rk4_step - advances f(t) to f(t+dt) using 4th order Runge Kutta method
!#
!#################################################################################
SUBROUTINE rk4_step(self)

    CLASS(rk4) :: self
    REAL*8, DIMENSION(self%nvar) :: f_p, k1, k2, k3, k4 ! trial variables

!### numerical recipes eqn. 17.1.3
    k1 = self%dt * self%dfdt(self%f, self%t)
    f_p = self%f + 0.5d0*k1

    k2 = self%dt * self%dfdt(f_p, self%t + 0.5d0*self%dt)
    f_p = self%f + 0.5d0*k2

    k3 = self%dt * self%dfdt(f_p, self%t + 0.50*self%dt)
    f_p = self%f + k3

    k4 = self%dt * self%dfdt(f_p, self%t + self%dt)

    self%f = self%f + (1.d0/6.d0) * (k1 + 2.d0*k2 + 2.d0*k3 + k4) 

!### advance independent variable by stepsize
    self%t = self%t + self%dt
END SUBROUTINE rk4_step
!#################################################################################
!#   dfdt - dummy function, provides generic interface.
!#       overridden by dfdt of specific problem
!################################################################################
FUNCTION dfdt(self,dep,ind)
    CLASS(rk4) :: self;    REAL*8, DIMENSION(self%nvar), INTENT(IN) :: dep
    REAL*8, INTENT(IN) :: ind;    REAL*8, DIMENSION(self%nvar) :: dfdt
END FUNCTION dfdt
END MODULE mod_rk4
!#################################################################################
!# MODULE: mod_prob3 - provides growth rates dxdt & dydt for predator-prey model
!#      
!#################################################################################
MODULE mod_prob3;   USE mod_rk4;    IMPLICIT NONE
TYPE, EXTENDS(rk4) :: rk4_prob3
    CONTAINS
    PROCEDURE :: dfdt=>dfdt_prob3 !### override calls to dummy dfdt with dfdt for prob3
END TYPE rk4_prob3
CONTAINS
FUNCTION dfdt_prob3(self,dep,ind)
    CLASS(rk4_prob3) :: self
    REAL*8, DIMENSION(self%nvar), INTENT(IN) :: dep
    REAL*8, INTENT(IN) :: ind
    REAL*8, DIMENSION(self%nvar) :: dfdt_prob3

    REAL*8 :: eps_a,eps_b,gam_a,gam_b

    eps_a = 1.d0; eps_b = 1.d0; gam_a = 1.d0; gam_b = 1.d0

    dfdt_prob3(1) = eps_a*dep(1) - gam_a*dep(1)*dep(2)
    dfdt_prob3(2) = -eps_b*dep(2) + gam_b*dep(1)*dep(2)

END FUNCTION dfdt_prob3;    END MODULE mod_prob3
!#################################################################################
!#  Main Program: solve ODE problem predator-prey model for two populations of fish
!#      type A and type B
!#
!#################################################################################
PROGRAM prob3
USE mod_prob3

IMPLICIT NONE
    TYPE(rk4_prob3) :: pred_prey
    REAL*8, DIMENSION(2) :: f_0 = (/ 0.1d0, 0.1d0/)     !   initial populations of A & B
    REAL*8 :: t = 0.d0                                  !   initial time
    REAL*8 :: tfinal = 9.999d0                          !   final time
    REAL*8 :: dt = 1.d-2                                !   timestep size

    OPEN(UNIT=10,FILE='pred_prey.dat')
100 FORMAT(3(1pe20.11))
!### initialize problem
    CALL pred_prey%init(2,t,dt,f_0)
    WRITE(10,100) pred_prey%t, pred_prey%f(:)
!### advance system by dt until t = 10, write out population at each step
    DO WHILE(pred_prey%t .lt. tfinal)
       CALL pred_prey%rk4_step()
       WRITE(10,100) pred_prey%t, pred_prey%f(:)
    END DO
END PROGRAM
