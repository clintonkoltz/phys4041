!#################################################################################
!#  MODULE: mod_datapoints - data storage and procedures for storing discrete 
!#      data points, dep = func(ind), and their finite difference derivatives
!#
!#################################################################################
MODULE mod_datapoints
IMPLICIT NONE
!### dataset - type contains arrays for discrete data points (ind,dep), and first 
!   and second finite difference derivative of dep w.r.t. ind evaluated at same 
!   data points. n_pts is # of data points, i.e. length of ind & dep
TYPE :: dataset
    REAL*8, DIMENSION(:), ALLOCATABLE :: ind,dep
    REAL*8, DIMENSION(:) ,ALLOCATABLE :: f_prime,f_doubprime
    INTEGER :: n_pts
    CONTAINS
    PROCEDURE :: build_data
    PROCEDURE :: destroy_data
    PROCEDURE :: frst_deriv
    PROCEDURE :: scnd_deriv
END TYPE dataset

CONTAINS
!###################################################################################
!#  build_data - allocates data structures based on input n_pts
!#  destroy_data - frees memory
!###################################################################################
SUBROUTINE build_data(self, n_pts)
    CLASS(dataset) :: self
    INTEGER :: n_pts

    self%n_pts = n_pts
    ALLOCATE(self%ind(n_pts),self%dep(n_pts),self%f_prime(n_pts),self%f_doubprime(n_pts))
END SUBROUTINE build_data
SUBROUTINE destroy_data(self)
    CLASS(dataset) :: self

    DEALLOCATE(self%ind,self%dep,self%f_prime,self%f_doubprime)
END SUBROUTINE destroy_data
!###################################################################################
!#  frst_deriv - approximate first derivative of stored data with forward difference 
!#      at leftmost point, backward difference at rightmost point, and central
!#      differences at all points in between
!###################################################################################
SUBROUTINE frst_deriv(self)
    CLASS(dataset) :: self
    INTEGER :: i
    REAL*8 :: delta

    delta = self%ind(2) - self%ind(1) !###  constant delta for evenly spaced data

    !### forward difference
    self%f_prime(1) = (self%dep(2) - self%dep(1)) / delta
    DO i = 2, self%n_pts-1  !### central differences
        self%f_prime(i) = (self%dep(i+1) - self%dep(i-1)) / (2.d0*delta)
    END DO
    !### backward diff
    self%f_prime(self%n_pts) = (self%dep(self%n_pts) - self%dep(self%n_pts-1)) / delta

END SUBROUTINE frst_deriv
!###################################################################################
!#  scnd_deriv - approximate discrete 2nd derivatives. use central difference 2nd
!#      derivative for interior points. for endpoints, forward difference
!#      2nd derivative at x_i = central difference 2nd derivative at x_i+1 and
!#      backward difference 2nd deriv at x_i = central difference 2nd derivative at x_i-1
!###################################################################################
SUBROUTINE scnd_deriv(self)
    CLASS(dataset) :: self
    INTEGER :: i
    REAL*8 :: delta

    delta = self%ind(2) - self%ind(1) !###  constant delta for evenly spaced data

    DO i = 2, self%n_pts-1  !### central differences
        self%f_doubprime(i) = (self%dep(i+1) - 2.d0*self%dep(i) + self%dep(i-1)) / delta**2
    END DO
!### forward diff at i is central diff at i+1, backward diff at i is central diff at i-1
    self%f_doubprime(1) = self%f_doubprime(2)
    self%f_doubprime(self%n_pts) = self%f_doubprime(self%n_pts-1)

END SUBROUTINE scnd_deriv
END MODULE mod_datapoints
!###################################################################################
!#  MODULE: mod_prob1 - function calc_stress calculates shear stress at stored
!#      data points from finite differences of discrete data points in data set
!#
!###################################################################################
MODULE mod_prob1
USE mod_spline
USE mod_datapoints
IMPLICIT NONE
TYPE, EXTENDS(dataset) :: turbprob
    CONTAINS
    PROCEDURE :: calc_stress
END TYPE turbprob
CONTAINS
FUNCTION calc_stress(self)
    CLASS(turbprob) :: self
    REAL*8, DIMENSION(self%n_pts) :: calc_stress

    CALL self%frst_deriv()
    CALL self%scnd_deriv()
    calc_stress = self%f_doubprime + self%f_prime / self%ind

END FUNCTION calc_stress
END MODULE mod_prob1
!#################################################################################
!#  Main Program: calculate spline for shear stress T(r) from experimental 
!#      measurements of turbulent flow in tube over range r=[0.1,0.9]
!#
!#################################################################################
PROGRAM prob1
USE mod_prob1
    INTEGER, PARAMETER :: n_pts = 9, N = 200 ! array size parameters
    TYPE(turbprob) :: vz_exp                   ! vz_exp - dataset for vz(r)
    TYPE(spline_gen) :: vz_spline            ! vz_spline - spline generator for vz(r)
    TYPE(turbprob) :: vz_int                 ! vz_int - dataset of interpolated points
    REAL*8, DIMENSION(n_pts) :: T_i          ! T_i - shear stress from FDs of experimental data
    REAL*8, DIMENSION(N) :: T_ii             ! T_ii - shear stress from FDs of vz_int
    REAL*8 :: r,dr
    INTEGER :: i

    OPEN(UNIT = 10, FILE='stress.dat')
    OPEN(UNIT = 20, FILE='stress_spline.dat')
200 FORMAT(2(1pe24.11))
!### build dataset memory, store experimental data
    CALL vz_exp%build_data(n_pts)
    vz_exp%ind = (/ 0.1d0,0.2d0,0.3d0,0.4d0,0.5d0,0.6d0,0.7d0,0.8d0,0.9d0 /)
    vz_exp%dep = (/ 1.13d0, 1.10d0, 1.05d0, 0.90d0, 0.82d0, &
            0.65d0, 0.43d0, 0.25d0, 0.12d0 /)

!### calculate shear stress using finite differences of vz(r)
    T_i =  vz_exp%calc_stress()
!### write out data 
    DO i = 1, vz_exp%n_pts
       WRITE(10,200) vz_exp%ind(i), T_i(i)
    END DO
!### build dataset and spline generator, and store data for cubic spline of vz(r)
    CALL vz_int%build_data(N)
    CALL vz_spline%build_spline(n_pts)
    vz_spline%x_data = vz_exp%ind
    vz_spline%y_data = vz_exp%dep
!### compute coefficients of cubic polynomials in intervals between discrete data points
    CALL compute_spline_coeffs(vz_spline)

!### evaluate spline of vz throughout interval, store in dataset
    dr = (vz_exp%ind(vz_exp%n_pts) - vz_exp%ind(1)) / dble(N)
    DO i = 1, N
       r = vz_exp%ind(1) + dr*dble(i)
       vz_int%ind(i) = r
       vz_int%dep(i) = spline_interp(vz_spline,r)
    END DO
!### calculate shear stress from FDs of spline of vz
    T_ii = vz_int%calc_stress()

!### write out values of shear stress
    DO i = 1, N
       r = vz_exp%ind(1) + dr*dble(i)
       WRITE(20,200) r, T_ii(i)
    END DO
!### free memory
    CALL vz_exp%destroy_data()
    CALL vz_spline%destroy_spline()
    CALL vz_int%destroy_data()
END PROGRAM
