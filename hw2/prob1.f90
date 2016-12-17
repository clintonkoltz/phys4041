!#################################################################################
!#
!#  Module: mod_bessel - subroutines for generating bessel functions by recursion
!#
!#  Subroutines: build - allocates space for bessel function array
!#               up_recursion - generates bessel function values by upward recursion
!#                  up to supplied value of l over supplied range in x
!#               down_recursion - calculates relative values of lower order bessel
!#                  functions through downward recursion. uses actual value of J0(x)
!#                  to normalize higher order recursively calculated J_l(x)
!#################################################################################
MODULE mod_bessel

IMPLICIT NONE

PUBLIC :: bessel_funcs

TYPE bessel_funcs
    !### J_l(l,ix) is value of l-th order bessel function at x = ix*dx + bounds(1)
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: J_l 
    REAL*8 :: dx = 1.d-2  !### dx - step size over x range
    !### bounds - range of x to plot over, set bounds(1)=bounds(2) if J_l only 
    !        needed at one x value, default is x=0 to x=12
    REAL*8, DIMENSION(2) :: bounds = (/ 0.d0, 12.d0 /)
    INTEGER :: nx !### nx - # of x values J_l calculated at, set in build subroutine
    INTEGER :: l = 8 !### caluclate up to l-th order, default is l=8

END TYPE bessel_funcs

CONTAINS

SUBROUTINE build(bessel)

    TYPE(bessel_funcs) :: bessel

    !### determine # of x values over range in order to allocate enough space
    !    in function value array J_l
    bessel%nx = int ( (bessel%bounds(2)-bessel%bounds(1) ) / bessel%dx )
    ALLOCATE( bessel%J_l( 0:bessel%l, 0:bessel%nx) )

END SUBROUTINE build

SUBROUTINE up_recursion(bessel)

    TYPE(bessel_funcs) :: bessel
    REAL*8 :: x_i
    INTEGER :: ix, il

    !### loop over x values of array
    DO ix = 0, bessel%nx

        x_i = bessel%bounds(1) + bessel%dx * dble(ix)

        !### store 1st & 2nd from explicit bessel functions
        bessel%J_l(0, ix) = sin(x_i) / x_i
        bessel%J_l(1, ix) = ( sin(x_i) - x_i*cos(x_i) ) / x_i**2

        !### recursively generate higher orders of J_l(x) from previous orders
        DO il = 1, bessel%l - 1
           bessel%J_l(il+1, ix) = (2.d0*dble(il) + 1.d0) / x_i &
            * bessel%J_l(il, ix) - bessel%J_l(il-1, ix)
        END DO

    END DO

END SUBROUTINE up_recursion

SUBROUTINE down_recursion(bessel)

    TYPE(bessel_funcs) :: bessel
    REAL*8 :: x_i
    REAL*8 :: J0
    INTEGER :: ix, il

    !### loop over x values of array
    DO ix = 0, bessel%nx

       !### store actual value of J_0(x)
       x_i = bessel%bounds(1) + bessel%dx * dble(ix)
       J0 = sin(x_i) / x_i

       !### random values for J_75(x) & J_74(x)
       bessel%J_l(bessel%l  ,ix) = 2.d0
       bessel%j_l(bessel%l-1,ix) = 1.d0

       !### recursively generate relative values for J_l-1 from random values 
       !    for arbitrarily high order bessel functions.
       DO il = bessel%l, 1, -1
          bessel%J_l(il-1, ix) = (2.d0*dble(il) + 1.d0) / x_i &
           *bessel%J_l(il, ix) - bessel%J_l(il+1, ix)
       END DO

       !### the ratio between actual J_0(x) and recursively calculated J_0(x)
       !    renormalizes J_l(x) values
       DO il = 1, bessel%l
          bessel%J_l(il, ix) = bessel%J_l(il, ix) * ( J0/ bessel%J_l(0,ix) )
       END DO

    END DO


END SUBROUTINE down_recursion

END MODULE mod_bessel


!#################################################################################
!#
!#  Main Program: tests routines in mod_bessel, provides solutions for problem 1
!#
!#################################################################################

PROGRAM problem1_main

USE mod_bessel

IMPLICIT NONE

    TYPE(bessel_funcs) :: up_plt, up_xp1, up_x1, up_x10
    TYPE(bessel_funcs) :: down_xp1, down_x1, down_x10
    INTEGER :: ii
    
    OPEN(unit = 10, file = 'bessel.dat' )

100 FORMAT(6(1pe24.10)) 
200 FORMAT(1f8.2,3(1pe20.10))

    !### calculate for plotting up to l = 4 over default range of x
    up_plt%l = 4
    CALL build(up_plt)
    CALL up_recursion(up_plt)

    !### write out x and all J_l(x) for plotting
    DO ii = 1, up_plt%nx
       WRITE(10,100) dble(ii)*up_plt%dx, up_plt%J_l(:,ii)
    END DO


    !### test upwards routine against values in table
    up_xp1%bounds = 0.1d0 ; up_x1%bounds = 1.d0 ;  up_x10%bounds = 10.d0

    CALL build(up_xp1) ; CALL build(up_x1) ; CALL build(up_x10) 
    CALL up_recursion(up_xp1) ; CALL up_recursion(up_x1) ; CALL up_recursion(up_x10)

    PRINT *, 'table calulated using upward recursion'
    PRINT *, '     x           j_3              j_5               j_8'
    PRINT 200, 0.1 , up_xp1%J_l(3,0), up_xp1%J_l(5,0), up_xp1%J_l(8,0)
    PRINT 200, 1.0 , up_x1%J_l(3,0), up_x1%J_l(5,0), up_x1%J_l(8,0)
    PRINT 200, 10.0, up_x10%J_l(3,0), up_x10%J_l(5,0), up_x10%J_l(8,0)
    PRINT *,''

    !### test downwards routine against values in table
    down_xp1%bounds = 0.1d0 ; down_x1%bounds = 1.d0 ; down_x10%bounds = 10.d0
    down_xp1%l = 75 ;  down_x1%l = 75 ; down_x10%l = 75
    
    CALL build(down_xp1) ; CALL build(down_x1) ; CALL build(down_x10)
    CALL down_recursion(down_xp1) ; CALL down_recursion(down_x1) ; CALL down_recursion(down_x10)

    PRINT *, 'table calculated using downward recurion'
    PRINT *, '     x           j_3              j_5               j_8'
    PRINT 200, 0.1 , down_xp1%J_l(3,0), down_xp1%J_l(5,0), down_xp1%J_l(8,0)
    PRINT 200, 1.0 , down_x1%J_l(3,0), down_x1%J_l(5,0), down_x1%J_l(8,0)
    PRINT 200, 10.0, down_x10%J_l(3,0), down_x10%J_l(5,0), down_x10%J_l(8,0)



END PROGRAM problem1_main
