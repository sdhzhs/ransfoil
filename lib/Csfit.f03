SUBROUTINE CSFIT ( N, X, Y, IENDL, DERIVL, IENDR, DERIVR, B, C, D, IER )

!IENDL     = 0: The 3rd derivative of the spline at the
!               left endpoint is to match the 3rd deriv-
!               ative of the cubic passing through the
!               first 4 data points; if N=3 and IENDL=
!               IENDR=0, both parts of the spline become
!               the same quadratic defined by the pts.;
!          = 1: The 1st derivative of the spline at the
!               left endpoint is to be the given DERIVL;
!          = 2: The 2nd derivative of the spline at the
!               left endpoint is to be the given DERIVL;
!          = 3: The 3rd derivative of the spline at the
!               left endpoint is to be the given DERIVL;
!          = 4: The spline and its first 3 derivatives
!               at X(N) are to match the values at X(1).
!               This is the cyclic, or periodic, case.
!DERIVL    Value of derivative used if IENDL=1, 2, or 3;
!          ignored if IENDL=0 or 4.
!IENDR,    As for IENDL, DERIVL, but pertaining to the
!DERIVR    right endpoint; ignored if IENDL=4.
!B,C,D     Spline coefficients (see PURPOSE and NOTES).
!IER       =0: No errors were detected.
!          =1: Too few data points; N < 2.
!          =2: Degenerate cases of N=2 for all end con-
!              ditions and of N=3 with second derivative
!              end condition require that the same order
!              of derivative is applied at both ends.
!          =3: Degenerate case where N=2 with third der-
!              ivative applied must have the same value
!              for the derivative at each endpoint.
!          =4: Cyclic mode -- Y(N) not equal to Y(1).
!          =5: Non-cyclic mode -- IENDL or IENDR is
!              beyond the range [0,4].

IMPLICIT NONE

INTEGER  N, IENDL, IENDR, IER
REAL(8)     X(N), Y(N), B(N), C(N), D(N)

REAL(8)     ZERO, ONE, TWO, THREE, FOUR, FIVE, SIX
PARAMETER ( ZERO=0.E+0, ONE=1.E+0, TWO=2.E+0, THREE=3.E+0, FOUR=4.E+0, FIVE=5.E+0, SIX=6.E+0 )


INTEGER  I, ISYS, NSYS
REAL(8)  DERIVL, DERIVR, DYDX1, DYDX2, DYDX3, H1, H2, H3

EXTERNAL  TRICPS

IER = 0
IF ( N<2 ) THEN
   IER = 1
   RETURN
END IF

IF ( N==2 ) THEN
   IF ( IENDL/=IENDR ) THEN
      IER = 2
      RETURN
   END IF

   H1    = X(2) - X(1)
   DYDX1 = (Y(2) - Y(1)) / H1

   IF ( IENDL==0 ) THEN
      B(1) = DYDX1
      C(1) = ZERO
      D(1) = ZERO
   ELSE IF ( IENDL==1 ) THEN
      B(1) = DERIVL
      C(1) = (DYDX1 * THREE - DERIVL * TWO - DERIVR) / H1
      D(1) = (DERIVL + DERIVR - DYDX1 * TWO) / H1**2
   ELSE IF ( IENDL==2 ) THEN
      B(1) = DYDX1 - (DERIVL * FOUR - DERIVR) * (H1 / SIX)
      C(1) = DERIVL / TWO
      D(1) = (DERIVR - DERIVL) / (H1 * SIX)
   ELSE IF ( IENDL==3 ) THEN
      IF ( DERIVL==DERIVR ) THEN
         B(1) = (Y(2) - Y(1)) / H1 - DERIVL / SIX * H1**2
         C(1) = ZERO
         D(1) = DERIVL / SIX
      ELSE
         IER = 3
         RETURN
      END IF
   ELSE IF ( IENDL==4 ) THEN
      IF ( Y(N)/=Y(1) ) THEN
         IER = 4
         RETURN
      END IF
      B(1) = ZERO
      C(1) = ZERO
      D(1) = ZERO
   END IF
   IF ( IENDL/=4 ) THEN
      H1   = X(N) - X(N-1)
      B(N) = B(N-1) + C(N-1) * H1 * TWO + D(N-1) * H1**2 * THREE
      C(N) = C(N-1) + D(N-1) * H1 * THREE
      D(N) = D(N-1)
   ELSE
      B(N) = B(1)
      C(N) = C(1)
      D(N) = D(1)
   END IF
   RETURN

ELSE IF ( N==3 ) THEN
   H1  = X(2) - X(1)
   H2  = X(3) - X(2)
   DYDX1 = (Y(2) - Y(1)) / H1
   DYDX2 = (Y(3) - Y(2)) / H2

   IF ( IENDL==2 .OR. IENDR==2 ) THEN
      IF ( IENDL/=IENDR ) THEN
         IER = 2
         RETURN
      END IF

      C(1) = DERIVL / SIX
      C(2) = (DYDX2 - DYDX1 - (DERIVL * H1 + DERIVR * H2) / SIX) / ((H1 + H2) * TWO)
      C(3) = DERIVR / SIX
      DO I = 1, N-1
         H1   = X(I+1) - X(I)
         B(I) = (Y(I+1) - Y(I)) / H1 - (C(I+1) + C(I)*TWO) * H1
         D(I) = (C(I+1) - C(I)) / H1
         C(I) = C(I) * THREE
      END DO

      IF ( IENDL/=4 ) THEN
         H1   = X(N) - X(N-1)
         B(N) = B(N-1) + C(N-1) * H1 * TWO + D(N-1) * H1**2 * THREE
         C(N) = C(N-1) + D(N-1) * H1 * THREE
         D(N) = D(N-1)
      ELSE
         B(N) = B(1)
         C(N) = C(1)
         D(N) = D(1)
      END IF
   
      RETURN  
   ELSE IF ( IENDL==4 ) THEN
      IF ( Y(N)/=Y(1) ) THEN
         IER = 4
         RETURN
      END IF
      C(1) = (DYDX1 - DYDX2) / (H1 + H2)
      C(2) =-C(1)
      C(3) = C(1)
      DO I = 1, N-1
         H1   = X(I+1) - X(I)
         B(I) = (Y(I+1) - Y(I)) / H1 - (C(I+1) + C(I)*TWO) * H1
         D(I) = (C(I+1) - C(I)) / H1
         C(I) = C(I) * THREE
      END DO

      IF ( IENDL/=4 ) THEN
         H1   = X(N) - X(N-1)
         B(N) = B(N-1) + C(N-1) * H1 * TWO + D(N-1) * H1**2 * THREE
         C(N) = C(N-1) + D(N-1) * H1 * THREE
         D(N) = D(N-1)
      ELSE
         B(N) = B(1)
         C(N) = C(1)
         D(N) = D(1)
      END IF
   
      RETURN
   END IF
ELSE IF ( N==4 ) THEN
   IF ( IENDL==2 .AND. IENDR==2 ) THEN
      H1  = X(2) - X(1)
      H2  = X(3) - X(2)
      H3  = X(4) - X(3)
      DYDX1 = (Y(2) - Y(1)) / H1
      DYDX2 = (Y(3) - Y(2)) / H2
      DYDX3 = (Y(4) - Y(3)) / H3

      C(1)  = DERIVL / SIX
      C(2)  = (DYDX3 - DYDX2 - DERIVR * H3 / SIX - (H3 / H2 + ONE) * (DYDX2 - DYDX1 - DERIVL * H1 / SIX) * TWO) / &
             (H2 - (H2 + H3) * (H1 / H2 + ONE) * FOUR)
      C(3)  = (DYDX2 - DYDX1 - DERIVL * H1 / SIX - (H1 / H2 + ONE) * (DYDX3 - DYDX2 - DERIVR * H3 / SIX) * TWO) / &
             (H2 - (H2 + H3) * (H1 / H2 + ONE) * FOUR)
      C(4)  = DERIVR / SIX
      DO I = 1, N-1
         H1   = X(I+1) - X(I)
         B(I) = (Y(I+1) - Y(I)) / H1 - (C(I+1) + C(I)*TWO) * H1
         D(I) = (C(I+1) - C(I)) / H1
         C(I) = C(I) * THREE
      END DO

      IF ( IENDL/=4 ) THEN
         H1   = X(N) - X(N-1)
         B(N) = B(N-1) + C(N-1) * H1 * TWO + D(N-1) * H1**2 * THREE
         C(N) = C(N-1) + D(N-1) * H1 * THREE
         D(N) = D(N-1)
      ELSE
         B(N) = B(1)
         C(N) = C(1)
         D(N) = D(1)
      END IF
   
      RETURN
   END IF
END IF

D(1) = X(2) - X(1)

DO I = 2, N-1
   D(I) = X(I+1) - X(I)
   B(I) = ( X(I+1) - X(I-1) ) * TWO
   C(I) = ( Y(I+1) - Y(I) ) / D(I) - ( Y(I) - Y(I-1) ) / ( X(I) - X(I-1) )
END DO

ISYS = 1
IF ( IENDL/=4 ) THEN

   NSYS = N
   D(N) = ZERO

   IF ( IENDL==0 ) THEN
      B(1) = -D(1)
      C(1) = ZERO
      IF ( N>3 ) C(1) = (C(3)/(X(4)-X(2)) - C(2)/(X(3)-X(1)))*D(1)**2/(X(4)-X(1))
   ELSE IF ( IENDL==1 ) THEN
      B(1) = D(1) * TWO
      C(1) = - ( DERIVL - ( Y(2) - Y(1) )/D(1) )
   ELSE IF ( IENDL==2 ) THEN
      NSYS = NSYS - 1
      ISYS = 2
      C(1) = DERIVL/SIX
      C(2) = C(2) - D(1) * C(1)
   ELSE IF ( IENDL==3 ) THEN
      B(1) = -D(1)
      C(1) = ( DERIVL/SIX ) * D(1)**2
   ELSE
      IER  = 5
      RETURN
   END IF

   IF ( IENDR==0 ) THEN
      B(N) = -D(N-1)
      C(N) = ZERO
      IF ( N>3 ) C(N) = (C(N-2)/(X(N-1)-X(N-3)) - C(N-1)/(X(N)-X(N-2)))*D(N-1)**2/(X(N)-X(N-3))
   ELSE IF ( IENDR==1 ) THEN
      B(N) = D(N-1) * TWO
      C(N) = DERIVR - (Y(N)-Y(N-1))/D(N-1)
   ELSE IF ( IENDR==2 ) THEN
      NSYS   = NSYS - 1
      C(N)   = DERIVR/SIX
      C(N-1) = C(N-1) - D(N-1) * C(N)
      D(N-1) = ZERO
   ELSE IF ( IENDR==3 ) THEN
      B(N) = -D(N-1)
      C(N) = -(DERIVR/SIX) * D(N-1)**2
   ELSE
      IER  = 5
      RETURN
   END IF
ELSE

   IF ( Y(N)/=Y(1) ) THEN
      IER = 4
      RETURN
   END IF

   NSYS = N-1
   B(1) = (D(1) + D(NSYS)) * TWO
   C(1) = (Y(2) - Y(1)) / D(1) - (Y(N) - Y(NSYS)) / D(NSYS)

END IF

CALL TRICPS ( NSYS, B(ISYS), D(ISYS), C(ISYS), C(ISYS) )

IF ( IENDL==4 ) C(N) = C(1)

DO I = 1, N-1
   H1   = X(I+1) - X(I)
   B(I) = (Y(I+1) - Y(I)) / H1 - (C(I+1) + C(I)*TWO) * H1
   D(I) = (C(I+1) - C(I)) / H1
   C(I) = C(I) * THREE
END DO

IF ( IENDL/=4 ) THEN
   H1   = X(N) - X(N-1)
   B(N) = B(N-1) + C(N-1) * H1 * TWO + D(N-1) * H1**2 * THREE
   C(N) = C(N-1) + D(N-1) * H1 * THREE
   D(N) = D(N-1)
ELSE
   B(N) = B(1)
   C(N) = C(1)
   D(N) = D(1)
END IF

RETURN
END SUBROUTINE CSFIT

SUBROUTINE CSEVAL(NDATA, X, Y, NU, U, B, C, D, S)

IMPLICIT NONE

REAL(8) ONE
PARAMETER (ONE = 1.0E+0)

INTEGER NDATA, NU
REAL(8) B(*), C(*), D(*), S(NU), U(NU), X(*), Y(*)

INTEGER IU, LEFT, NX
REAL(8) ARROW, DX, PERIOD, XEVAL, XLEFT, XRIGHT

ARROW  = SIGN (ONE, X (2) - X (1))
NX     = ABS (NDATA)
XRIGHT = X (NX)
LEFT   = 1

IF (NDATA > 0) THEN

   DO IU = 1, NU

      XEVAL = U (IU)

      CALL INTERVAL (NX, X, XEVAL, ARROW, LEFT)

      DX = XEVAL - X (LEFT)

      S (IU) = Y (LEFT) + (DX * (B (LEFT) + (DX * (C (LEFT) + (DX * (D (LEFT)))))))

      IF (XEVAL == XRIGHT) S (IU) = Y (NX)

   END DO

ELSE

   IF (ARROW == ONE) THEN
      XLEFT = X (1)
   ELSE
      XLEFT = XRIGHT
      XRIGHT = X (1)
   END IF

   PERIOD = XRIGHT - XLEFT

   DO IU = 1, NU

      XEVAL = U (IU)

      IF (XEVAL < XLEFT) THEN

         XEVAL = XRIGHT - MOD (XLEFT - XEVAL, PERIOD)

      ELSE IF (XEVAL > XRIGHT) THEN

         XEVAL = XLEFT + MOD (XEVAL - XRIGHT, PERIOD)

      END IF

      CALL INTERVAL (NX, X, XEVAL, ARROW, LEFT)

      DX = XEVAL - X (LEFT)

      S (IU) = Y (LEFT) + (DX * (B (LEFT) + (DX * (C (LEFT) + (DX * (D (LEFT)))))))

      IF (XEVAL == X (NX)) S (IU) = Y (NX)

   END DO

END IF

RETURN
END SUBROUTINE CSEVAL

SUBROUTINE TRICPS ( N, D, L, R, S )

IMPLICIT NONE

INTEGER  N
REAL(8)  D(N), L(N), R(N), S(N)

INTEGER  I
REAL(8)  DI, DINV, LIM1

DINV = 1.E+0 / D(1)
D(1) = L(N) * DINV
S(1) = R(1) * DINV
D(N) = D(N) - L(N) * D(1)
S(N) = R(N) - L(N) * S(1)

DO I = 2, N-2
   LIM1   = L(I-1)
   L(I-1) = LIM1 * DINV
   DI     = D(I) - LIM1 * L(I-1)
   DINV   = 1.E+0 / DI
   D(I)   = - LIM1 * D(I-1) * DINV
   S(I)   = ( R(I) - LIM1 * S(I-1) ) * DINV
   D(N)   = D(N) - ( DI * D(I) ) * D(I)
   S(N)   = S(N) - ( DI * D(I) ) * S(I)
END DO

I      = N-1
LIM1   = L(I-1)
L(I-1) = LIM1 * DINV
DI     = D(I) - LIM1 * L(I-1)
D(I)   = ( L(I) - LIM1 * D(I-1) ) / DI
S(I)   = ( R(I) - LIM1 * S(I-1) ) / DI
D(N)   = D(N) - ( DI * D(I) ) * D(I)
S(N)   = ( S(N) - ( DI * D(I) ) * S(I) ) / D(N)
L(I)   = 0.E+0

DO I = N-1, 1, -1
   S(I) = S(I) - L(I) * S(I+1) - D(I)*S(N)
END DO

RETURN
END SUBROUTINE TRICPS

SUBROUTINE INTERVAL (NX, X, XFIND, ARROW, LEFT)

IMPLICIT NONE

REAL(8) ONE
PARAMETER (ONE = 1.0E+0)

INTEGER LEFT, NX
REAL(8) ARROW, X (NX), XFIND

INTEGER LENGTH, NXLESS1, RIGHT, TRIAL
REAL(8) XBYARROW

XBYARROW = XFIND * ARROW

NXLESS1 = NX - 1

IF (XBYARROW >= X (NXLESS1) * ARROW) THEN
   LEFT = NXLESS1
   RETURN
ELSE IF (XBYARROW < X (2) * ARROW) THEN
   LEFT = 1
   RETURN
END IF

LEFT = MIN (MAX (2, LEFT), NX - 2)

IF (XBYARROW >= X (LEFT) * ARROW) THEN
   IF (XBYARROW < X (LEFT + 1) * ARROW) THEN
      RETURN
   ELSE
      RIGHT = NXLESS1
      LEFT  = LEFT + 1
   END IF
ELSE
   RIGHT = LEFT
   LEFT  = 2
END IF

DO WHILE(.TRUE.)
   LENGTH = RIGHT - LEFT
   IF (LENGTH > 1) THEN

      TRIAL = MIN (RIGHT - 1, LEFT + MAX (0, INT (REAL (LENGTH) * (XFIND - X (LEFT)) / (X (RIGHT) - X (LEFT)))))

      IF (XBYARROW >= X (TRIAL + 1) * ARROW) THEN
         LEFT  = TRIAL + 1
      ELSE IF (XBYARROW < X (TRIAL) * ARROW) THEN
         RIGHT = TRIAL
      ELSE
         LEFT  = TRIAL
         EXIT
      END IF
   ELSE
      EXIT
   END IF
END DO 

RETURN
END SUBROUTINE INTERVAL