SUBROUTINE cdivid(ar, ai, br, bi, cr, ci)
! COMPLEX DIVISION C = A/B, AVOIDING OVERFLOW.

REAL (dp), INTENT(IN)   :: ar
REAL (dp), INTENT(IN)   :: ai
REAL (dp), INTENT(IN)   :: br
REAL (dp), INTENT(IN)   :: bi
REAL (dp), INTENT(OUT)  :: cr
REAL (dp), INTENT(OUT)  :: ci

REAL (dp) :: r, d, t, infin

IF (br == 0.0_dp .AND. bi == 0.0_dp) THEN

! DIVISION BY ZERO, C = INFINITY.
  CALL mcon(t, infin, t, t)
  cr = infin
  ci = infin
  RETURN
END IF
IF (ABS(br) < ABS(bi)) THEN
  r = br / bi
  d = bi + r * br
  cr = (ar*r+ai) / d
  ci = (ai*r-ar) / d
  RETURN
END IF
r = bi / br
d = br + r * bi
cr = (ar+ai*r) / d
ci = (ai-ar*r) / d
RETURN
END SUBROUTINE cdivid
