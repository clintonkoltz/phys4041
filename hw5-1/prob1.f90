!#################################################################################
!# MODULE: mod_domain1D - type 'field1D' stores a scalar quantity at discrete
!#      positions, x_i, seperated by uniform intervals 'dx' 
!#
!#################################################################################
MODULE mod_domain1D
IMPLICIT NONE
    REAL*8 :: t = 0.d0          ! t - time, initialized to 0
    REAL*8 :: dt                ! dt - timestep
    TYPE :: field1D
       REAL*8, DIMENSION(:), ALLOCATABLE :: qty
       REAL*8, DIMENSION(2) :: xbounds
       REAL*8 :: dx
       INTEGER :: nx
    CONTAINS
       PROCEDURE :: build_domain
       PROCEDURE :: data_out
    END TYPE field1D
CONTAINS
!#################################################################################
!#  build_domain - set interval dx, allocate array space for domain plus ghost zones
!#      for periodic boundary conditions
!#################################################################################
SUBROUTINE build_domain(self,N,xbounds)
    CLASS(field1D) :: self
    REAL*8, DIMENSION(2) :: xbounds
    INTEGER, INTENT(IN) :: N

    self%nx = N
    self%xbounds = xbounds
    self%dx = (self%xbounds(2)-self%xbounds(1)) / dble(N-1)
    ALLOCATE(self%qty(0:N+1))
END SUBROUTINE build_domain
!#################################################################################
!#  data_out - data is dumped to new file every 'dumpt', file is named based on 
!#      current value of 't' and input value 'a'
!#################################################################################
SUBROUTINE data_out(self,a)
    CLASS(field1D) :: self
    CHARACTER(LEN=8) :: c1,c2       ! c1,c2 - dummy strings
    CHARACTER(LEN=24) :: filename
    REAL*8, INTENT(IN) :: a
    REAL*8 :: x
    INTEGER :: i,id

    id = nint(t)
    WRITE (c1,'(f8.4)') t   ! this is how you convert reals into strings in fortran,
    WRITE (c2,'(f8.2)') a   !    pretend the character string is a file and write to it

!### trim(adjustl('string')) to eliminate blank characters from string.
!    '//' concatenates stings together
    filename = 't-'//trim(adjustl(c1))//'_a-'//trim(adjustl(c2))//'.dat'
    OPEN( UNIT=10+id, FILE=filename )

!### write out quantity vs. x to file
    DO i = 1, self%nx
       x = self%xbounds(1) + self%dx*dble(i-1)
       WRITE(10+id,'(1p2e20.10)') x, self%qty(i)
    END DO

    CLOSE(10+id)

END SUBROUTINE data_out
END MODULE mod_domain1D
!#################################################################################
!#  MODULE: mod_laxwendroff - procedures to solve the advection equation
!#      ( du/dt = -v*dx/dt ) using the Lax-Wendroff method
!#
!#################################################################################
MODULE mod_laxwendroff
USE mod_domain1D
IMPLICIT NONE
    REAL*8 :: a_lw = 0.99d0    ! a_lw - CFL parameter, initialized to stable value
    REAL*8 :: v = 1.d0         ! v - advection velocity

CONTAINS
!#################################################################################
!#  init_prob - initializes sinusoidal distribution of u(x). use input 'CFL'
!#      to calculate size of timestep
!#################################################################################
SUBROUTINE init_prob(u,CFL)
    TYPE(field1D), INTENT(INOUT) :: u
    REAL*8,PARAMETER :: twopi = 8.d0*atan(1.d0)
    REAL*8, INTENT(IN) :: CFL                   ! CFL - stability parameter
    REAL*8 :: x
    INTEGER :: i

!### set time to 0 & set timestep
    t = 0.d0
    dt = CFL*u%dx/abs(v)
!### set values of quantity through
    DO i = 0, u%nx+1
       x = u%xbounds(1) + dble(i-1)*u%dx
       u%qty(i) = sin(twopi*x)
    END DO

END SUBROUTINE init_prob
!#################################################################################
!#  lw_step - advance advected quantity from u^n to u^{n+1} using L-W method
!#
!#################################################################################
SUBROUTINE lw_step(u)
    TYPE(field1D), INTENT(INOUT) :: u
    REAL*8, DIMENSION(0:u%nx) :: F  ! F = v*du/dx , F_i is flux from x_i to x_{i+1}
    INTEGER :: i

!### advance timestep
    t = t + dt
!### apply periodic boundary conditions
    u%qty(0) = u%qty(u%nx)
    u%qty(u%nx+1) = u%qty(1)
!### find fluxes at boundaries
    DO i = 0, u%nx
       F(i) = (v/2.d0) * (u%qty(i+1) + u%qty(i)) &
            - v**2*dt/(u%dx*2.d0) * (u%qty(i+1) - u%qty(i))
!        print *, F(i)
    END DO
!### update quantity values
    DO i = 1, u%nx
       u%qty(i) = u%qty(i) - (dt/u%dx) * (F(i) - F(i-1))
    END DO

END SUBROUTINE lw_step
END MODULE mod_laxwendroff
!#################################################################################
!#  Main Program: advect 1D quantity 'u' with periodic boundary conditions through
!#      a few periods. compare staple timestep size vs. unstable timestep size
!#
!#################################################################################
PROGRAM prob1
USE mod_laxwendroff
    TYPE(field1D) :: u
!### tfinal - end time ; dumpt - interval between datadumps ; tdump - time to
!   write out data
    REAL*8 :: tfinal = 2.d0, dumpt = 1.d0, tdump = 0.1d0
    REAL*8 :: x
    INTEGER :: i

!### build array N = 256. initialize problem for a = 0.99
    CALL u%build_domain(256, (/0.d0,1.d0/) )
    CALL init_prob(u,a_lw)
!### write out initial distribution
    CALL u%data_out(a_lw)
    i = 1
!### step forward in time, write data out at t = 0.1 (for comparison with a_lw > 1
!   system), then write out every t = 1 until final time
    DO WHILE( t .lt. tfinal)
       CALL lw_step(u)
       IF ( t .ge. tdump ) THEN
          tdump = dble(i)*dumpt
          i = i + 1
          CALL u%data_out(a_lw)
       END IF
    END DO
!### demonstrate instability of system with timestep calculated from CFL
!   coefficient a_lw > 1
    a_lw = 1.1d0 ;  tfinal = 0.2d0 ; dumpt = 0.02d0 ; tdump = 0.02d0
    CALL init_prob(u,a_lw)
    CALL u%data_out(a_lw)
    DO WHILE( t .lt. tfinal)
       CALL lw_step(u)
       IF ( t .ge. tdump ) THEN
          tdump = tdump + dumpt
          CALL u%data_out(a_lw)
       END IF
    END DO

END PROGRAM
!#################################################################################
