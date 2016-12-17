!#################################################################################
!#
!#  Module: sorting_algorithms -
!#            routines that sort single precision reals in numerical order
!#
!#################################################################################
MODULE sorting_algorithms
IMPLICIT NONE

CONTAINS

!#################################################################################
!#  insertion sort - an outer loop runs over all the elements in the array except the
!#  first. an inner loop swaps the element backwards until in proper sequence.
!#################################################################################
SUBROUTINE insertion_sort(array)

    REAL*4, INTENT(INOUT), DIMENSION(:) :: array
    REAL*4 :: temp_var
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
!#################################################################################
!#  insertion sort abs, using absolute value instead for swap criteria, output
!#  array should be original elements sorted by abs value, not an array of the 
!#  abs values of the elements
!#################################################################################
SUBROUTINE insertion_sort_abs(array)

    REAL*4, INTENT(INOUT), DIMENSION(:) :: array
    REAL*4 :: temp_var
    INTEGER :: array_length
    INTEGER :: ii,jj

    array_length = size(array)


    DO ii = 2, array_length
        jj = ii
        DO WHILE ( ( jj .gt. 1 ) .AND.  ( abs(array(jj-1)) .gt. abs(array(jj)) ) )
            temp_var = array(jj-1)
            array(jj-1) = array(jj)
            array(jj) = temp_var
            jj = jj - 1
        END DO
    END DO

END SUBROUTINE insertion_sort_abs

END MODULE sorting_algorithms
!#################################################################################
!#
!#   Module: mod_random - 
!#          FORTRAN's intrinsic random number routines are tedious. this provides
!#          rand_seed() subroutine which insures random seeds used by
!#          pseudorandom number generator are different with each execution
!#
!#################################################################################
MODULE mod_random

IMPLICIT NONE

CONTAINS

    SUBROUTINE rand_seed()
       implicit none
       integer,dimension(8) :: time_ints
       integer :: time_sum

       !### returns 8 integers: year, month, day, etc. stored in time_ints
       CALL date_and_time(VALUES = time_ints)

       !### intrinsic random_seed(PUT) requires 12 integer array, use sum of 
       !    time_ints to fill out input
       time_sum = sum(time_ints)
       CALL random_seed( PUT = (/ time_ints, time_sum, time_sum**2, &
          int(sqrt(real(time_sum))), int(real(time_sum)**1.5) /) )

    END SUBROUTINE

END MODULE mod_random
!#################################################################################
!#
!#  Main Program:  test sorting algorithms using real array of pseudorandom numbers
!#
!#################################################################################
program problem3_main

USE mod_random          !### load modules
USE sorting_algorithms
IMPLICIT NONE

    REAL*4, DIMENSION(6) :: original_array, array

    !### this generates pseudorandom numbers between -1 & 1 for each element in array
    CALL rand_seed()
    CALL random_number(original_array)
    original_array = 2. * original_array - 1.

    array = original_array

    print *, 'This is the original array:'
    print *, array

    CALL insertion_sort(array)

    print *, 'This is the array sorted by value:'
    print *, array

    array = original_array

    CALL insertion_sort_abs(array)

    print *, 'This is the array sorted by absolute value:'
    print *, array

end program problem3_main
