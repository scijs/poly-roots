SUBROUTINE fxshft(l2, zr, zi, conv)
! COMPUTES L2 FIXED-SHIFT H POLYNOMIALS AND TESTS FOR CONVERGENCE.
! INITIATES A VARIABLE-SHIFT ITERATION AND RETURNS WITH THE
! APPROXIMATE ZERO IF SUCCESSFUL.
! L2 - LIMIT OF FIXED SHIFT STEPS
! ZR,ZI - APPROXIMATE ZERO IF CONV IS .TRUE.
! CONV  - LOGICAL INDICATING CONVERGENCE OF STAGE 3 ITERATION

INTEGER, INTENT(IN)     :: l2
REAL (dp), INTENT(OUT)  :: zr
REAL (dp), INTENT(OUT)  :: zi
LOGICAL, INTENT(OUT)    :: conv

REAL (dp) :: otr, oti, svsr, svsi
LOGICAL   :: test, pasd, bool
INTEGER   :: i, j, n

n = nn - 1

! EVALUATE P AT S.
CALL polyev(nn, sr, si, pr, pi, qpr, qpi, pvr, pvi)
test = .true.
pasd = .false.

! CALCULATE FIRST T = -P(S)/H(S).
CALL calct(bool)

! MAIN LOOP FOR ONE SECOND STAGE STEP.
DO  j = 1, l2
  otr = tr
  oti = ti

! COMPUTE NEXT H POLYNOMIAL AND NEW T.
  CALL nexth(bool)
  CALL calct(bool)
  zr = sr + tr
  zi = si + ti

! TEST FOR CONVERGENCE UNLESS STAGE 3 HAS FAILED ONCE OR THIS
! IS THE LAST H POLYNOMIAL.
  IF (.NOT.(bool.OR..NOT.test.OR.j == l2)) THEN
    IF (cmod(tr-otr,ti-oti) < .5_dp*cmod(zr,zi)) THEN
      IF (pasd) THEN

! THE WEAK CONVERGENCE TEST HAS BEEN PASSED TWICE, START THE THIRD STAGE
! ITERATION, AFTER SAVING THE CURRENT H POLYNOMIAL AND SHIFT.
        DO  i = 1, n
          shr(i) = hr(i)
          shi(i) = hi(i)
        END DO
        svsr = sr
        svsi = si
        CALL vrshft(10,zr,zi,conv)
        IF (conv) RETURN

! THE ITERATION FAILED TO CONVERGE. TURN OFF TESTING AND RESTORE H,S,PV AND T.
        test = .false.
        DO  i = 1, n
          hr(i) = shr(i)
          hi(i) = shi(i)
        END DO
        sr = svsr
        si = svsi
        CALL polyev(nn, sr, si, pr, pi, qpr, qpi, pvr, pvi)
        CALL calct(bool)
        CYCLE
      END IF
      pasd = .true.
    ELSE
      pasd = .false.
    END IF
  END IF
END DO

! ATTEMPT AN ITERATION WITH FINAL H POLYNOMIAL FROM SECOND STAGE.
CALL vrshft(10, zr, zi, conv)
RETURN
END SUBROUTINE fxshft
