FUNCTION errev(nn, qr, qi, ms, mp, are, mre) RESULT(fn_val)
! BOUNDS THE ERROR IN EVALUATING THE POLYNOMIAL BY THE HORNER RECURRENCE.
! QR,QI - THE PARTIAL SUMS
! MS    -MODULUS OF THE POINT
! MP    -MODULUS OF POLYNOMIAL VALUE
! ARE, MRE -ERROR BOUNDS ON COMPLEX ADDITION AND MULTIPLICATION

INTEGER, INTENT(IN)    :: nn
REAL (dp), INTENT(IN)  :: qr(:)
REAL (dp), INTENT(IN)  :: qi(:)
REAL (dp), INTENT(IN)  :: ms
REAL (dp), INTENT(IN)  :: mp
REAL (dp), INTENT(IN)  :: are
REAL (dp), INTENT(IN)  :: mre
REAL (dp)              :: fn_val

REAL (dp) :: e
INTEGER   :: i

e = cmod(qr(1), qi(1)) * mre / (are+mre)
DO  i = 1, nn
  e = e * ms + cmod(qr(i), qi(i))
END DO
fn_val = e * (are+mre) - mp * mre
RETURN
END FUNCTION errev
