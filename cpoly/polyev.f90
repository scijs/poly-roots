SUBROUTINE polyev(nn, sr, si, pr, pi, qr, qi, pvr, pvi)
! EVALUATES A POLYNOMIAL  P  AT  S  BY THE HORNER RECURRENCE
! PLACING THE PARTIAL SUMS IN Q AND THE COMPUTED VALUE IN PV.

INTEGER, INTENT(IN)     :: nn
REAL (dp), INTENT(IN)   :: sr
REAL (dp), INTENT(IN)   :: si
REAL (dp), INTENT(IN)   :: pr(:)
REAL (dp), INTENT(IN)   :: pi(:)
REAL (dp), INTENT(OUT)  :: qr(:)
REAL (dp), INTENT(OUT)  :: qi(:)
REAL (dp), INTENT(OUT)  :: pvr
REAL (dp), INTENT(OUT)  :: pvi

REAL (dp) :: t
INTEGER   :: i

qr(1) = pr(1)
qi(1) = pi(1)
pvr = qr(1)
pvi = qi(1)
DO  i = 2, nn
  t = pvr * sr - pvi * si + pr(i)
  pvi = pvr * si + pvi * sr + pi(i)
  pvr = t
  qr(i) = pvr
  qi(i) = pvi
END DO
RETURN
END SUBROUTINE polyev
