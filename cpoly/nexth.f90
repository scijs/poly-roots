SUBROUTINE nexth(bool)
! CALCULATES THE NEXT SHIFTED H POLYNOMIAL.
! BOOL   -  LOGICAL, IF .TRUE. H(S) IS ESSENTIALLY ZERO

LOGICAL, INTENT(IN)  :: bool

REAL (dp) :: t1, t2
INTEGER   :: j, n

n = nn - 1
IF (.NOT.bool) THEN
  DO  j = 2, n
    t1 = qhr(j-1)
    t2 = qhi(j-1)
    hr(j) = tr * t1 - ti * t2 + qpr(j)
    hi(j) = tr * t2 + ti * t1 + qpi(j)
  END DO
  hr(1) = qpr(1)
  hi(1) = qpi(1)
  RETURN
END IF

! IF H(S) IS ZERO REPLACE H WITH QH.
DO  j = 2, n
  hr(j) = qhr(j-1)
  hi(j) = qhi(j-1)
END DO
hr(1) = 0.0_dp
hi(1) = 0.0_dp
RETURN
END SUBROUTINE nexth
