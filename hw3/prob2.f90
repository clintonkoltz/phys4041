!#################################################################################
!#  Module: mod_schrod_eigen - defines type 'eigen_prob' which builds storage
!#      for LAPACK subroutine 'dgeev', general eigenvalue problem solver.
!#      see http://www.netlib.org/lapack/explore-html/d9/d28/dgeev_8f.html
!#      compile with: gfortran prob2.f90 -o prob2 -llapack
!#################################################################################
MODULE mod_eigen
IMPLICIT NONE

TYPE :: eigen_prob

    REAL*8, DIMENSION(:,:), ALLOCATABLE :: A,VR,VL
    REAL*8, DIMENSION(:), ALLOCATABLE :: WR,WI,WORK
    INTEGER :: N,LWORK,INFO

    CONTAINS

    PROCEDURE :: build

END TYPE eigen_prob

CONTAINS
!#################################################################################
!#  build - allocate space for eigen_prob object
!#
!#################################################################################
SUBROUTINE build(self,N)

    CLASS(eigen_prob) :: self
    INTEGER,INTENT(IN) :: N

    self%N = N
    self%LWORK = 4.d0 * self%N
    ALLOCATE(self%A(self%N,self%N), self%VR(self%N,self%N), self%VL(self%N,self%N))
    ALLOCATE(self%WR(self%N), self%WI(self%N), self%WORK(self%LWORK))

END SUBROUTINE build

END MODULE mod_eigen

MODULE sorting_algorithms
IMPLICIT NONE

CONTAINS
!#################################################################################
!#  insertion sort - an outer loop runs over all the elements in the array except the
!#  first. an inner loop swaps the element backwards until in proper sequence.
!#
!#################################################################################
SUBROUTINE insertion_sort(array)
    REAL*8, INTENT(INOUT), DIMENSION(:) :: array
    REAL*8 :: temp_var
    INTEGER :: array_length
    INTEGER :: ii,jj

    array_length = size(array)

    DO ii = 2, array_length
        jj = ii
        DO WHILE ( (jj .gt. 1) .AND.  (array(jj-1) .gt. array(jj)) )
            temp_var = array(jj-1)
            array(jj-1) = array(jj)
            array(jj) = temp_var
            jj = jj - 1
        END DO
    END DO

END SUBROUTINE insertion_sort

END MODULE sorting_algorithms
!#################################################################################
!#  Main Program: test eigen_prob module for truncated harmonic oscillator potential
!#
!#################################################################################
PROGRAM prob2

USE mod_eigen
USE sorting_algorithms

IMPLICIT NONE

    REAL*8, PARAMETER :: pi = 4.d0*atan(1.d0) !### arctan(1) = pi/4
    TYPE(eigen_prob) :: prob
    REAL*8 :: dj,di
    INTEGER :: i,j


    CALL prob%build(10)

    DO i = 1, prob%N
       di = dble(i)
       DO j = 1, prob%N
         dj = dble(j)
         IF (i.ne.j) THEN   !### see  hw3.pdf
            prob%A(i,j) = ( 2.d0*di*dj*(-1)**(i+j) ) / ( (di-dj)**2*(di+dj)**2*pi**2 )
         ELSE
            prob%A(i,j) = ( 2.d0 - 3.d0/(i*pi)**2 + 12.d0*(i*pi)**2 )/ 24.d0
         END IF
       END DO
    END DO

    !### matrix A that goes into eigensolver is 2*H (see hw3.pdf)
    prob%A = 2.d0* prob%A

    CALL dgeev('N','V',prob%N,prob%A,prob%N,prob%WR,prob%WI,prob%VL,prob%N, &
         prob%VR,prob%N,prob%WORK,prob%LWORK,prob%INFO)

    !### dgeev returns eigenvalues unsorted, sort for plotting
    CALL insertion_sort(prob%WR)

    OPEN(UNIT=10, FILE='schrod.dat')
    OPEN(UNIT=20, FILE='schrodlog.dat')

    DO i = 1, prob%N

       !### plot E_n vs. n
       WRITE(10,'(i3,1pe20.10)') i, prob%WR(i)
       !### write log(E_n) vs. log(n) to determine power-law index
       WRITE(20,'(2(1pe20.10))') log(dble(i)), log(prob%WR(i))
       !### linear regression evaluates to m~1.995
       !    E_n ~ n^1.995, similar to scaling for infinite square well (E_n ~ n^2)
    END DO

END PROGRAM prob2
