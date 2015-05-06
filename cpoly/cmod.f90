FUNCTION cmod(r, i) RESULT(fn_val)
! MODULUS OF A COMPLEX NUMBER AVOIDING OVERFLOW.

REAL (dp), INTENT(IN)  :: r
REAL (dp), INTENT(IN)  :: i
REAL (dp)              :: fn_val

REAL (dp) :: ar, ai

ar = ABS(r)
ai = ABS(i)
IF (ar < ai) THEN
  fn_val = ai * SQRT(1.0_dp + (ar/ai)**2)
  RETURN
END IF
IF (ar > ai) THEN
  fn_val = ar * SQRT(1.0_dp + (ai/ar)**2)
  RETURN
END IF
fn_val = ar * SQRT(2.0_dp)
RETURN
END FUNCTION cmod
