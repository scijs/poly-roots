SUBROUTINE noshft(l1)
! COMPUTES  THE DERIVATIVE  POLYNOMIAL AS THE INITIAL H
! POLYNOMIAL AND COMPUTES L1 NO-SHIFT H POLYNOMIALS.

INTEGER, INTENT(IN)  :: l1

REAL (dp) :: xni, t1, t2
INTEGER   :: i, j, jj, n, nm1

n = nn - 1
nm1 = n - 1
DO  i = 1, n
  xni = nn - i
  hr(i) = xni * pr(i) / n
  hi(i) = xni * pi(i) / n
END DO
DO  jj = 1, l1
  IF (cmod(hr(n), hi(n)) > eta*10.0_dp*cmod(pr(n), pi(n))) THEN
    CALL cdivid(-pr(nn), -pi(nn), hr(n), hi(n), tr, ti)
    DO  i = 1, nm1
      j = nn - i
      t1 = hr(j-1)
      t2 = hi(j-1)
      hr(j) = tr * t1 - ti * t2 + pr(j)
      hi(j) = tr * t2 + ti * t1 + pi(j)
    END DO
    hr(1) = pr(1)
    hi(1) = pi(1)
  ELSE

! IF THE CONSTANT TERM IS ESSENTIALLY ZERO, SHIFT H COEFFICIENTS.
    DO  i = 1, nm1
      j = nn - i
      hr(j) = hr(j-1)
      hi(j) = hi(j-1)
    END DO
    hr(1) = 0.0_dp
    hi(1) = 0.0_dp
  END IF
END DO
RETURN
END SUBROUTINE noshft
