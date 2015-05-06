FUNCTION scale(nn, pt, eta, infin, smalno, base) RESULT(fn_val)
! RETURNS A SCALE FACTOR TO MULTIPLY THE COEFFICIENTS OF THE POLYNOMIAL.
! THE SCALING IS DONE TO AVOID OVERFLOW AND TO AVOID UNDETECTED UNDERFLOW
! INTERFERING WITH THE CONVERGENCE CRITERION.  THE FACTOR IS A POWER OF THE
! BASE.
! PT - MODULUS OF COEFFICIENTS OF P
! ETA, INFIN, SMALNO, BASE - CONSTANTS DESCRIBING THE FLOATING POINT ARITHMETIC.

INTEGER, INTENT(IN)    :: nn
REAL (dp), INTENT(IN)  :: pt(:)
REAL (dp), INTENT(IN)  :: eta
REAL (dp), INTENT(IN)  :: infin
REAL (dp), INTENT(IN)  :: smalno
REAL (dp), INTENT(IN)  :: base
REAL (dp)              :: fn_val

REAL (dp) :: hi, lo, MAX, MIN, x, sc
INTEGER   :: i, l

! FIND LARGEST AND SMALLEST MODULI OF COEFFICIENTS.
hi = SQRT(infin)
lo = smalno / eta
MAX = 0.0_dp
MIN = infin
DO  i = 1, nn
  x = pt(i)
  IF (x > MAX) MAX = x
  IF (x /= 0.0_dp .AND. x < MIN) MIN = x
END DO

! SCALE ONLY IF THERE ARE VERY LARGE OR VERY SMALL COMPONENTS.
fn_val = 1.0_dp
IF (MIN >= lo .AND. MAX <= hi) RETURN
x = lo / MIN
IF (x <= 1.0_dp) THEN
  sc = 1.0_dp / (SQRT(MAX)*SQRT(MIN))
ELSE
  sc = x
  IF (infin/sc > MAX) sc = 1.0_dp
END IF
l = LOG(sc) / LOG(base) + .500
fn_val = base ** l
RETURN
END FUNCTION scale
