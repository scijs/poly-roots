SUBROUTINE calct(bool)
! COMPUTES  T = -P(S)/H(S).
! BOOL   - LOGICAL, SET TRUE IF H(S) IS ESSENTIALLY ZERO.

LOGICAL, INTENT(OUT)  :: bool

REAL (dp) :: hvr, hvi
INTEGER   :: n

n = nn - 1

! EVALUATE H(S).
CALL polyev(n, sr, si, hr, hi, qhr, qhi, hvr, hvi)
bool = cmod(hvr,hvi) <= are * 10.0_dp * cmod(hr(n), hi(n))
IF (.NOT.bool) THEN
  CALL cdivid(-pvr, -pvi, hvr, hvi, tr, ti)
  RETURN
END IF
tr = 0.0_dp
ti = 0.0_dp
RETURN
END SUBROUTINE calct
