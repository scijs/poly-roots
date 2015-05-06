SUBROUTINE cauchy(nn, pt, q, fn_val)
! CAUCHY COMPUTES A LOWER BOUND ON THE MODULI OF THE ZEROS OF A
! POLYNOMIAL - PT IS THE MODULUS OF THE COEFFICIENTS.

INTEGER, INTENT(IN)     :: nn
REAL (dp), INTENT(OUT)  :: pt(:), q(:), fn_val

REAL (dp) :: x, xm, f, dx, df
INTEGER   :: i, n

pt(nn) = -pt(nn)

! COMPUTE UPPER ESTIMATE OF BOUND.
n = nn - 1
x = EXP((LOG(-pt(nn)) - LOG(pt(1))) / n)
IF (pt(n) /= 0.0_dp) THEN

! IF NEWTON STEP AT THE ORIGIN IS BETTER, USE IT.
  xm = -pt(nn) / pt(n)
  IF (xm < x) x = xm
END IF

! CHOP THE INTERVAL (0,X) UNITL F LE 0.
10 xm = x * .1_dp
f = pt(1)
DO  i = 2, nn
  f = f * xm + pt(i)
END DO
IF (f > 0.0_dp) THEN
  x = xm
  GO TO 10
END IF
dx = x

! DO NEWTON ITERATION UNTIL X CONVERGES TO TWO DECIMAL PLACES.
30 IF (ABS(dx/x) > .005_dp) THEN
  q(1) = pt(1)
  DO  i = 2, nn
    q(i) = q(i-1) * x + pt(i)
  END DO
  f = q(nn)
  df = q(1)
  DO  i = 2, n
    df = df * x + q(i)
  END DO
  dx = f / df
  x = x - dx
  GO TO 30
END IF
fn_val = x

RETURN
END SUBROUTINE cauchy

