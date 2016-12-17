!#################################################################################
!#  MODULE: mod_domain2D - 'field2D' stores a scalar quantity on an xy grid with
!#      uniform grid spacing 'del'
!#
!#################################################################################
MODULE mod_domain2D
IMPLICIT NONE
    REAL*8 :: t = 0.d0          ! t - time, initialized to 0
    REAL*8 :: dt                ! dt - timestep
    TYPE :: field2D
       REAL*8, DIMENSION(:,:), ALLOCATABLE :: qty
       REAL*8, DIMENSION(2) :: xbounds,ybounds
       REAL*8 :: del
       INTEGER :: nx,ny
    END TYPE field2D
CONTAINS
!#################################################################################
!#  build_domain - calculate 'nx' & 'ny' (# of grid zones along x and y) and 
!#      allocate memory space
!#################################################################################
SUBROUTINE build_domain(fld,del,xbounds,ybounds)
    CLASS(field2D) :: fld
    REAL*8, DIMENSION(2),INTENT(IN) :: xbounds,ybounds
    REAL*8, INTENT(IN) :: del
    INTEGER  :: nx,ny

    fld%del = del
    fld%xbounds = xbounds  ; fld%ybounds = ybounds
    fld%nx = 1 + nint( (xbounds(2)-xbounds(1)) / del )
    fld%ny = 1 + nint( (ybounds(2)-ybounds(1)) / del )
    ALLOCATE(fld%qty(fld%nx,fld%ny))
END SUBROUTINE build_domain
END MODULE mod_domain2D
!#################################################################################
!#  MODULE: mod_diffus_ADI - procedures to perform alternating direction
!#      implicit method to solve 2D diffusion equation
!#################################################################################
MODULE mod_diffus_ADI
USE mod_domain2D
IMPLICIT NONE
    REAL*8 :: kappa = 1.d0          ! kappa - diffusion coefficient
    REAL*8 :: alpha                 ! alpha - tridag. eqn. parameter, see hw5.pdf

CONTAINS
!#################################################################################
!#  adi_step - performs ADI algorithm by calling x implicit first half of time
!#      step, then y implicit second half, see hw5.pdf
!#################################################################################
SUBROUTINE adi_step(u)
    TYPE(field2D), INTENT(INOUT) :: u

    CALL x_imp(u)
    CALL y_imp(u)
END SUBROUTINE adi_step
!#################################################################################
!#  x_imp - advance half time step with x-Laplacian implicitly taken at t_{n+1/2}
!#      and y-Laplacian taken explicitly at t_n. tridag solver routine returns 
!#      distribution T^{n+1/2}
!#################################################################################
SUBROUTINE x_imp(u)
    TYPE(field2D), INTENT(INOUT) :: u   
!### xmaj_vec - vector transformation of matrix such that consecutive x-values are
!   sequentially arranged ; a,b,c - tridiagonal matrix elements, see hw5.pdf ;
!   rhs - right-hand side of equation ( d in hw5.pdf )
    REAL*8, DIMENSION(u%nx*u%ny) :: xmaj_vec, a, b, c, rhs
    INTEGER :: i,j,l

    DO j = 1, u%ny
       DO i = 1, u%nx
          l = (j-1)*u%nx + i
          a(l) = -0.5d0* alpha
          b(l) = 1.d0 + alpha
          c(l) = -0.5d0* alpha
          rhs(l) = 0.5d0* alpha* u%qty(i,j+1) + (1.d0 - alpha)* u%qty(i,j) &
                 + 0.5d0* alpha* u%qty(i,j-1)
!### Neumann boundary conditions
          IF( j.eq.1 ) THEN
             rhs(i) = (1.d0 - alpha)* u%qty(i,j) + alpha* u%qty(i,j+1)
          END IF
!### Dirichlet boundary conditions
          IF( (i.eq.1) .or. (i.eq.u%nx) .or. (j.eq.u%ny) ) THEN
             a(l) = 0.d0;    b(l) = 1.d0;    c(l) = 0.d0;    rhs(l) = u%qty(i,j)
          END IF
       END DO
    END DO

    CALL tridag(a,b,c,rhs,xmaj_vec,u%nx*u%ny)
!### map solution xmaj_vec back to matrix elements    
    DO j = 1, u%ny
       DO i = 1, u%nx
          l = (j-1)*u%nx + i
          u%qty(i,j) = xmaj_vec(l)
       END DO
    END DO
END SUBROUTINE x_imp
!#################################################################################
!#  y_imp - advance half time step with y-Laplacian implicitly taken at t_{n+1}
!#      and x-Laplacian taken explicitly at t_{n+1/2}. tridag solver routine
!#      returns distribution T^{n+1}
!#################################################################################
SUBROUTINE y_imp(u)
    TYPE(field2D), INTENT(INOUT) :: u
!### ymaj_vec - vector transformation of matrix such that consecutive y-values are
!   sequentially arranged ; a,b,c - tridiagonal matrix elements, see hw5.pdf ;
!   rhs - right-hand side of equation ( d in hw5.pdf )
    REAL*8, DIMENSION(u%nx*u%ny) :: ymaj_vec, rhs, a, b, c
    INTEGER :: i,j,m

    DO i = 1, u%nx
       DO j = 1, u%ny
          m = (i-1)*u%ny + j
          a(m) = -0.5d0* alpha
          b(m) = 1.d0 + alpha
          c(m) = -0.5d0* alpha
          rhs(m) = 0.5d0* alpha* u%qty(i+1,j) + (1.d0 - alpha)* u%qty(i,j) &
                 + 0.5d0* alpha* u%qty(i-1,j)
!### Neumann boundary conditions
          IF( j.eq.1 ) THEN
             a(m) = 0.d0;   c(m) = - alpha
          END IF
!### Dirichlet boundary conditions
          IF( (i.eq.1) .or. (i.eq.u%nx) .or. (j.eq.u%ny) ) THEN
             a(m) = 0.d0;   b(m) = 1.d0;    c(m) = 0.d0;    rhs(m) = u%qty(i,j)
          END IF
       END DO
    END DO

    CALL tridag(a,b,c,rhs,ymaj_vec,u%nx*u%ny)
!### map solution ymaj_vec back to matrix elements    
    DO i = 1, u%nx
       DO j = 1, u%ny
          m = (i-1)*u%ny + j
          u%qty(i,j) = ymaj_vec(m)
       END DO
    END DO
END SUBROUTINE y_imp
END MODULE mod_diffus_ADI
!#################################################################################
!#  MODULE: mod_laxwendroff - procedure to solve advection equation of a
!#      quantity in a 2D field
!#################################################################################
MODULE mod_laxwendroff
USE mod_domain2D
IMPLICIT NONE
    REAL*8 :: a_lw = 0.9d0         ! a_lw - CFL parameter
    REAL*8 :: v = 1.d0             ! v - advection velocity

CONTAINS
!#################################################################################
!#  lw_step_x - advect quantity in the x direction using full timestep value
!#      assuming a constant velocity 'v'. in prob2, this represents the flow of
!#      ground water through porous rock
!#################################################################################
SUBROUTINE lw_step_x(u)
    TYPE(field2D), INTENT(INOUT) :: u
!### F = v * du/dx. F_ij is flux from zone i,j into i+1,j
    REAL*8, DIMENSION(u%nx-1,u%ny) :: F
    INTEGER :: i

!### find fluxes at x-interfaces between zones
    DO i = 1, u%nx-1
       F(i,:) = (v/2.d0) * (u%qty(i+1,:) + u%qty(i,:)) &
              - v**2*dt/(u%del*2.d0) * (u%qty(i+1,:) - u%qty(i,:))
    END DO
!### boundary conditions i = 1 : u = const; i = nx : u = const; j = ny : u = const
!   therefore fluxes for j = ny zones are 0; for i = 1 & i = nx, u doesn't get updated
    F(:,u%ny) = 0.d0
!### update change in quantity
    DO i = 2, u%nx-1
       u%qty(i,:) = u%qty(i,:) - (dt/u%del) * (F(i,:) - F(i-1,:))
    END DO

END SUBROUTINE lw_step_x
END MODULE mod_laxwendroff
!#################################################################################
!#  Main Program: solve combined advection-diffusion equation to find steady
!#      state temperature distribution over 2D domain
!#
!#################################################################################
PROGRAM prob2
USE mod_laxwendroff
USE mod_diffus_ADI
    TYPE(field2D) :: T_ij
!### T_old - needed to check change in T per timestep ; T_1, T_N - boundary values
!   tfinal - end point in case no steady state is reached ; delta - grid spacing
!   maxresid, tol - to check for steady state convergence
    REAL*8,DIMENSION(:,:),ALLOCATABLE :: T_old
    REAL*8,DIMENSION(2) :: xbounds, ybounds
    REAL*8 :: tfinal = 1000.d0, delta = 0.25d0 
    REAL*8 :: T_1 = 65.d0, T_N = 25.d0
    REAL*8 :: maxresid = 1.d0, tol = 1.d-8
    REAL*8 :: x,y
    INTEGER :: i,j

    OPEN( UNIT=10, FILE='T_xy.dat')
    OPEN( UNIT=20, FILE='resid_t.dat')

!### set domain boundaries and build data structures
    xbounds = (/ 0.d0, 30.d0*kappa/v /)
    ybounds = (/ 0.d0, 10.d0*kappa/v /)
    CALL build_domain(T_ij, delta, xbounds, ybounds)
    ALLOCATE( T_old(T_ij%nx,T_ij%ny) )

!### set time step and tridag eqn parameter
    dt = a_lw*T_ij%del/abs(v)
    alpha = kappa * dt/ T_ij%del**2

!### initialize temperature distribution with constant gradient from x = 0 to x = 30kap/v
    DO i = 1, T_ij%nx
       x = T_ij%xbounds(1) + dble(i-1)*T_ij%del
       T_ij%qty(i,:) = T_1 + (T_N - T_1)/(T_ij%xbounds(2)-T_ij%xbounds(1)) * x
    END DO
!### ...but if y = 10kap/v, T = 25 = const
    T_ij%qty(:,T_ij%ny) = T_N

!### integrate forward in time until steady state is reached
    DO WHILE(maxresid .gt. tol)
       maxresid = 0.d0
       T_old = T_ij%qty          ! store T(x,y) from previous step
       CALL lw_step_x(T_ij)
       CALL adi_step(T_ij)
       DO j = 1, T_ij%ny
          DO i = 1, T_ij%nx
             ! check residuals at i,j: resid_ij = |T^n+1 - T^n|/T^n
             maxresid = max( abs((T_ij%qty(i,j)-T_old(i,j))/T_old(i,j)), maxresid )
          END DO
       END DO
       t = t + dt
!### quit if max time reached
       IF ( t.gt.tfinal ) EXIT
!### write out file to demonstrate convergence to steady state
       WRITE(20,'(1p2e22.10)') t, maxresid 
    END DO
!### write out final steady state temperature distribution
    DO j = 1, T_ij%ny
       y = T_ij%ybounds(1) + dble(j-1)*T_ij%del
       DO i = 1, T_ij%nx
          x = T_ij%xbounds(1) + dble(i-1)*T_ij%del
          WRITE(10,'(1p3e22.10)') x, y, T_ij%qty(i,j)
       END DO
    END DO
END PROGRAM
!#################################################################################
