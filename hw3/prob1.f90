!#################################################################################
!#  Module: mod_poisson - sturctures and routines for solving 2D poisson's equation
!#      nabla^2 V(x,y) = - rho(x,y)
!#
!#################################################################################
MODULE mod_poisson
IMPLICIT NONE

TYPE :: poi_prob !### poisson problem object
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: rho !### rho - 2D array storing source term values at grid points
    REAL*8, DIMENSION(:), ALLOCATABLE :: V_vec !### V_vec - vector form of 2D potential array 
                                               !    see Numerical Recipes sect. 20.0, (j,l) -> i= (j-1)*L + l
    REAL*8 :: del_x         !### - grid spacing
    REAL*8 :: tol = 1.d-8   !### - tolerance for checking solution
    REAL*8 :: max_resid     !### - maximum residual to check against tolerance 
    INTEGER :: N,M  !### N - grid size along x, M - grid size along y

    CONTAINS !### type bound procedures

    PROCEDURE :: build
    PROCEDURE :: trans2vec
    PROCEDURE :: trans2mat
    PROCEDURE :: iter_step
    PROCEDURE :: poi_solve

END TYPE poi_prob

CONTAINS
!#################################################################################
!#  build - allocates space for arrays
!#
!#################################################################################
SUBROUTINE build(self)

    CLASS(poi_prob) :: self !### passed-object dummy argument (see line 111), must be declared CLASS

    ALLOCATE(self%rho(self%N,self%M),self%V_vec(self%N*self%M))
    self%max_resid = 1.d2 * self%tol

END SUBROUTINE build
!#################################################################################
!#  trans2vec - takes initial guess for potential, transforms to 1D vec
!#  trans2mat - takes solution V_vec, transforms back to 2D matrix form
!#
!################################################################################
SUBROUTINE trans2vec(self,input_array)
    CLASS(poi_prob) :: self
    REAL*8, DIMENSION(:,:), INTENT(IN) :: input_array
    INTEGER :: i,j,k

    DO i = 1, self%N
       DO j = 1, self%M
          k = (i-1)*self%N + j
          self%V_vec(k) = input_array(i,j)
       END DO
    END DO
END SUBROUTINE trans2vec

SUBROUTINE trans2mat(self,output_array)
    CLASS(poi_prob) :: self
    REAL*8, DIMENSION(:,:) ,INTENT(OUT) :: output_array
    INTEGER :: i,j,k

    DO i = 1, self%N
       DO j = 1, self%M
          k = (i-1)*self%N + j
          output_array(i,j) = self%V_vec(k)
       END DO
    END DO
END SUBROUTINE trans2mat
!################################################################################
!#  iterstep - perform one iteration of Jacobi Method for solving Poisson eqn
!#      by relaxation. (Jacobi Method is SUPERSLOW. for large problems always, 
!#      use Gauss-Seidel along with multigrid and/or successive overrelaxation
!#      (Numerical Recipes sect. 20.5-20.6)
!################################################################################
SUBROUTINE iter_step(self)

    CLASS(poi_prob) :: self
    REAL*8, DIMENSION(self%N*self%M) :: V_old
    INTEGER :: i,j,k

    V_old = self%V_vec

    self%max_resid = 0.d0

    DO i = 2, self%N-1
       DO j = 2, self%M-1
          k = (i-1)*self%N + j
          !### Numerical Recipes eqns 20.0.8 & 20.5.5
          self%V_vec(k) = 0.25d0* & 
            ( V_old(k-self%N) + V_old(k-1) + V_old(k+1) + V_old(k+self%N) &
            - self%del_x**2*self%rho(i,j) )
          !### residuals are absolute difference betw. old guess and new guess
          self%max_resid = max( abs(self%V_vec(k) - V_old(k)), self%max_resid )
       END DO
    END DO

END SUBROUTINE iter_step
!#################################################################################
!#  poi_solve - transforms potential 2D array into vector form to build matrix equation
!#      as in Fig. 20.0.3 in Numerical Recipes. calls iterstep until maximum residual 
!#      value in our solution for the potential is less than specified tolerance.
!#      transform solution back into 2D array
!#################################################################################
SUBROUTINE poi_solve(self,potential,density)

    CLASS(poi_prob) :: self
    REAL*8,DIMENSION(:,:),INTENT(IN) :: density
    REAL*8,DIMENSION(:,:),INTENT(INOUT) :: potential
    INTEGER :: i = 0, max_iter = 100000

    self%N = size(potential,1)
    self%M = size(potential,2)
    !###  calling type-bound procedures is similar to accessing attributes of a
    !   derived type [variable_name%procedure_name()]. variable_name is implicitly
    !   passed to procedure_name by default
    CALL self%build() !### allocate data structures

    CALL self%trans2vec(potential)  !### transform 2D array to 1D array
    self%rho = density              !### store source function in object

    DO WHILE (self%max_resid .gt. self%tol .and. i .lt. max_iter)
       CALL self%iter_step()        !### iterate until solution is found
       i = i + 1
    END DO

    CALL self%trans2mat(potential)  !### transform solution back to 2D array

END SUBROUTINE poi_solve

END MODULE mod_poisson

!#################################################################################
!#  Main Program: For fun let's solve for and plot the potential well of a 
!#      point source
!#
!#################################################################################
PROGRAM prob1
USE mod_poisson

IMPLICIT NONE

    INTEGER,PARAMETER :: N = 40 
    TYPE(poi_prob) :: my_poi
    REAL*8, DIMENSION(N,N) :: my_potential = 0.d0   !### initial guess to relax from
    REAL*8, DIMENSION(N,N) :: my_source = 0.d0      !### source function is zero everywhere except...
    REAL*8 :: x,y
    INTEGER :: i,j


    my_source(N/2,N/2) = 10.d0                      !  ...at one point in center
    my_poi%del_x = 0.1d0                            !### grid spacing

    CALL my_poi%poi_solve(my_potential,my_source)   !### call poisson solver

    !### write out potential for plotting
    OPEN(UNIT=10, FILE='potential.dat')
100 FORMAT(3(1pe20.10))
    DO i=1,N
       x =  i*my_poi%del_x
       DO j=1,N
          y = j*my_poi%del_x
          WRITE(10,100) x, y, my_potential(i,j)
       END DO
    END DO

END PROGRAM
