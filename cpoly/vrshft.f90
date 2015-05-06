SUBROUTINE vrshft(l3, zr, zi, conv)
! CARRIES OUT THE THIRD STAGE ITERATION.
! L3 - LIMIT OF STEPS IN STAGE 3.
! ZR,ZI   - ON ENTRY CONTAINS THE INITIAL ITERATE, IF THE
!           ITERATION CONVERGES IT CONTAINS THE FINAL ITERATE ON EXIT.
! CONV    -  .TRUE. IF ITERATION CONVERGES

INTEGER, INTENT(IN)        :: l3
REAL (dp), INTENT(IN OUT)  :: zr
REAL (dp), INTENT(IN OUT)  :: zi
LOGICAL, INTENT(OUT)       :: conv

REAL (dp) :: mp, ms, omp, relstp, r1, r2, tp
LOGICAL   :: b, bool
INTEGER   :: i, j

conv = .false.
b = .false.
sr = zr
si = zi

! MAIN LOOP FOR STAGE THREE
DO  i = 1, l3

! EVALUATE P AT S AND TEST FOR CONVERGENCE.
  CALL polyev(nn, sr, si, pr, pi, qpr, qpi, pvr, pvi)
  mp = cmod(pvr,pvi)
  ms = cmod(sr,si)
  IF (mp <= 20.0_dp*errev(nn, qpr, qpi, ms, mp, are, mre)) THEN

! POLYNOMIAL VALUE IS SMALLER IN VALUE THAN A BOUND ON THE ERROR
! IN EVALUATING P, TERMINATE THE ITERATION.
    conv = .true.
    zr = sr
    zi = si
    RETURN
  END IF
  IF (i /= 1) THEN
    IF (.NOT.(b .OR. mp < omp .OR. relstp >= .05_dp)) THEN

! ITERATION HAS STALLED. PROBABLY A CLUSTER OF ZEROS.  DO 5 FIXED
! SHIFT STEPS INTO THE CLUSTER TO FORCE ONE ZERO TO DOMINATE.
      tp = relstp
      b = .true.
      IF (relstp < eta) tp = eta
      r1 = SQRT(tp)
      r2 = sr * (1.0_dp+r1) - si * r1
      si = sr * r1 + si * (1.0_dp+r1)
      sr = r2
      CALL polyev(nn, sr, si, pr, pi, qpr, qpi, pvr, pvi)
      DO  j = 1, 5
        CALL calct(bool)
        CALL nexth(bool)
      END DO
      omp = infin
      GO TO 20
    END IF

! EXIT IF POLYNOMIAL VALUE INCREASES SIGNIFICANTLY.
    IF (mp*.1_dp > omp) RETURN
  END IF
  omp = mp

! CALCULATE NEXT ITERATE.
  20 CALL calct(bool)
  CALL nexth(bool)
  CALL calct(bool)
  IF (.NOT.bool) THEN
    relstp = cmod(tr,ti) / cmod(sr,si)
    sr = sr + tr
    si = si + ti
  END IF
END DO
RETURN
END SUBROUTINE vrshft
