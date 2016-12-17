! tridag.f
! Numerical Recipes (http://www.nr.com), Section 2.4

      SUBROUTINE tridag(a,b,c,r,u,n)

      INTEGER n
      REAL*8 a(n),b(n),c(n),r(n),u(n)
! INTEGER NMAX
! PARAMETER (NMAX=500)
! Solves for a vector u(1:n) of length n the tridiagonal linear set given by equation (2.4.1).
! a(1:n), b(1:n), c(1:n), and r(1:n) are input vectors and are not modified.
! Parameter: NMAX is the maximum expected value of n.
      INTEGER j
! One vector of workspa!e, gam is needed.
      REAL*8 bet
! REAL gam(NMAX)
      REAL*8 gam(n)

      if (b(1).eq.0.) then 
        print *,'tridag.f: rewrite equations' ; stop
      end if
! If this happens then you should rewrite your equations as a set of order N-1, with u_2
! trivially eliminated.
      bet=b(1)
      u(1)=r(1)/bet
! Decomposition and forward substitution.
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
! Algorithm fails; see below.
        if (bet.eq.0.) then
           print *, 'tridag.f: tridag failed' ; stop
        end if
        u(j)=(r(j)-a(j)*u(j-1))/bet
11 enddo
! Backsubstitution.
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12 enddo

      return
END


