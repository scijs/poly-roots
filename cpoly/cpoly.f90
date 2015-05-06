MODULE Solve_Complex_Poly
!     ALGORITHM 419 COLLECTED ALGORITHMS FROM ACM.

! Code converted using TO_F90 by Alan Miller
! Date: 2000-01-08  Time: 16:02:44

!     ALGORITHM APPEARED IN COMM. ACM, VOL. 15, NO. 02, P. 097.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(14, 60)

! COMMON AREA
! COMMON /global/ pr, pi, hr, hi, qpr, qpi, qhr, qhi, shr, shi, sr,  &
!    si, tr, ti, pvr, pvi, are, mre, eta, infin, nn
! REAL (dp) :: sr, si, tr, ti, pvr, pvi, are, mre, eta, infin,  &
!    pr(50), pi(50), hr(50), hi(50), qpr(50), qpi(50), qhr(50),  &
!    qhi(50), shr(50), shi(50)

REAL (dp), ALLOCATABLE, SAVE  :: pr(:), pi(:), hr(:), hi(:), qpr(:), qpi(:), &
                                 qhr(:), qhi(:), shr(:), shi(:)
REAL (dp), SAVE               :: sr, si, tr, ti, pvr, pvi, are, mre, eta,  &
                                 infin
INTEGER, SAVE                 :: nn

PRIVATE
PUBLIC  :: dp, cpoly


CONTAINS


SUBROUTINE cpoly(opr, opi, degree, zeror, zeroi, fail)
! FINDS THE ZEROS OF A COMPLEX POLYNOMIAL.
! OPR, OPI  -  REAL (dp) VECTORS OF REAL AND IMAGINARY PARTS OF THE
!              COEFFICIENTS IN ORDER OF DECREASING POWERS.
! DEGREE    -  INTEGER DEGREE OF POLYNOMIAL.
! ZEROR, ZEROI  -  OUTPUT REAL (dp) VECTORS OF REAL AND IMAGINARY
!                  PARTS OF THE ZEROS.
! FAIL      -  OUTPUT LOGICAL PARAMETER,  TRUE ONLY IF LEADING COEFFICIENT
!              IS ZERO OR IF CPOLY HAS FOUND FEWER THAN DEGREE ZEROS.
! THE PROGRAM HAS BEEN WRITTEN TO REDUCE THE CHANCE OF OVERFLOW OCCURRING.
! IF IT DOES OCCUR, THERE IS STILL A POSSIBILITY THAT THE ZEROFINDER WILL WORK
! PROVIDED THE OVERFLOWED QUANTITY IS REPLACED BY A LARGE NUMBER.

REAL (dp), INTENT(IN)   :: opr(:)
REAL (dp), INTENT(IN)   :: opi(:)
INTEGER, INTENT(IN)     :: degree
REAL (dp), INTENT(OUT)  :: zeror(:)
REAL (dp), INTENT(OUT)  :: zeroi(:)
LOGICAL, INTENT(OUT)    :: fail

REAL (dp) :: xx, yy, cosr, sinr, smalno, base, xxx, zr, zi, bnd
LOGICAL   :: conv
INTEGER   :: cnt1, cnt2, i, idnn2

! INITIALIZATION OF CONSTANTS
CALL mcon(eta, infin, smalno, base)
are = eta
mre = 2.0_dp * SQRT(2.0_dp) * eta
xx = .70710678
yy = -xx
cosr = -.060756474
sinr = .99756405
fail = .false.
nn = degree + 1

! ALGORITHM FAILS IF THE LEADING COEFFICIENT IS ZERO.
IF (opr(1) == 0.0_dp .AND. opi(1) == 0.0_dp) THEN
  fail = .true.
  RETURN
END IF

! Allocate arrays

IF (ALLOCATED(pr)) DEALLOCATE(pr, pi, hr, hi, qpr, qpi, qhr, qhi, shr, shi)

ALLOCATE(pr(nn), pi(nn), hr(nn), hi(nn), qpr(nn), qpi(nn), qhr(nn), qhi(nn), shr(nn), shi(nn))

! REMOVE THE ZEROS AT THE ORIGIN IF ANY.
10 IF (opr(nn) == 0.0_dp .AND. opi(nn) == 0.0_dp) THEN
  idnn2 = degree - nn + 2
  zeror(idnn2) = 0.0_dp
  zeroi(idnn2) = 0.0_dp
  nn = nn - 1
  GO TO 10
END IF

! MAKE A COPY OF THE COEFFICIENTS.
DO  i = 1, nn
  pr(i) = opr(i)
  pi(i) = opi(i)
  shr(i) = cmod(pr(i),pi(i))
END DO

! SCALE THE POLYNOMIAL.
bnd = scale(nn, shr, eta, infin, smalno, base)
IF (bnd /= 1.0_dp) THEN
  DO  i = 1, nn
    pr(i) = bnd * pr(i)
    pi(i) = bnd * pi(i)
  END DO
END IF


! START THE ALGORITHM FOR ONE ZERO .
40 IF (nn <= 2) THEN

! CALCULATE THE FINAL ZERO AND RETURN.
  CALL cdivid(-pr(2), -pi(2), pr(1), pi(1), zeror(degree), zeroi(degree))
  RETURN
END IF

! CALCULATE BND, A LOWER BOUND ON THE MODULUS OF THE ZEROS.
DO  i = 1, nn
  shr(i) = cmod(pr(i), pi(i))
END DO
CALL cauchy(nn, shr, shi, bnd)

! OUTER LOOP TO CONTROL 2 MAJOR PASSES WITH DIFFERENT SEQUENCES OF SHIFTS.
DO  cnt1 = 1, 2

  ! FIRST STAGE CALCULATION, NO SHIFT.
  CALL noshft(5)

  ! INNER LOOP TO SELECT A SHIFT.
  DO  cnt2 = 1, 9

    ! SHIFT IS CHOSEN WITH MODULUS BND AND AMPLITUDE ROTATED BY
    ! 94 DEGREES FROM THE PREVIOUS SHIFT
    xxx = cosr * xx - sinr * yy
    yy = sinr * xx + cosr * yy
    xx = xxx
    sr = bnd * xx
    si = bnd * yy

    ! SECOND STAGE CALCULATION, FIXED SHIFT.
    CALL fxshft(10*cnt2,zr,zi,conv)
    IF (conv) THEN

      ! THE SECOND STAGE JUMPS DIRECTLY TO THE THIRD STAGE ITERATION.
      ! IF SUCCESSFUL THE ZERO IS STORED AND THE POLYNOMIAL DEFLATED.
      idnn2 = degree - nn + 2
      zeror(idnn2) = zr
      zeroi(idnn2) = zi
      nn = nn - 1
      DO  i = 1, nn
        pr(i) = qpr(i)
        pi(i) = qpi(i)
      END DO
      GO TO 40
    END IF

  ! IF THE ITERATION IS UNSUCCESSFUL ANOTHER SHIFT IS CHOSEN.
  END DO

! IF 9 SHIFTS FAIL, THE OUTER LOOP IS REPEATED WITH ANOTHER SEQUENCE OF SHIFTS.
END DO

! THE ZEROFINDER HAS FAILED ON TWO MAJOR PASSES.
! RETURN EMPTY HANDED.
fail = .true.
RETURN
END SUBROUTINE cpoly
