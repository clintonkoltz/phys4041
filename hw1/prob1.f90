!#################################################################################
!#  subroutine: test_array_sum
!#
!#  Takes an array of real numbers as argument, returns +1 if the sum of
!#  elements is positive, 0 if the sum is zero, -1 if the sum is negative.
!#  Function uses an automatic array as input, to allow a rank 1 array 
!#  of arbitrary length.
!#
!#################################################################################
INTEGER FUNCTION test_array_sum(input_array)

    IMPLICIT NONE

    REAL*4, INTENT(IN), DIMENSION(:) :: input_array
    REAL*4 :: array_sum

!### intrinsic subroutine 'sum' is used, which sums all the elements of the input
!    array. compiled code from intrinsic processes like this are at least
!    as efficient as writing an explicit loop structure, can be more so
    array_sum = sum(input_array)

    IF(array_sum .GT. 0.) THEN
        test_array_sum = 1
    ELSE IF(array_sum .LT. 0.) THEN
        test_array_sum = -1
    ELSE IF(array_sum .EQ. 0.) THEN
        test_array_sum = 0
    END IF

END FUNCTION test_array_sum
!#################################################################################
!#
!#  Main Program - test subroutine 'testarray' for various arrays
!#
!#################################################################################
PROGRAM problem1_main

    IMPLICIT NONE
    !### tells compiler all variables must be explicitly defined. always include this,
    !    otherwise implicit typing means typos become hard-to-find bugs

    REAL*4 :: pos_array(4),zero_array(5),neg_array(6)
    INTERFACE
        INTEGER FUNCTION test_array_sum(input_array)
            REAL*4, INTENT(IN), DIMENSION(:) :: input_array
        END FUNCTION
    END INTERFACE

    pos_array = (/ 1., 2., 3., 4. /)

    print *, pos_array
    print *, 'test_array_sum(pos_array) = ', test_array_sum(pos_array)
    print *, ''

    zero_array = (/ -2., -1., 0., 1., 2. /)

    print *, zero_array
    print *, 'test_array_sum(zero_array) = ', test_array_sum(zero_array)
    print *, ''

    neg_array = (/ -3., -2., -1., 0., 1., 2. /)

    print *, neg_array
    print *, 'test_array_sum(neg_array) = ', test_array_sum(neg_array)
    print *, ''

END PROGRAM problem1_main
