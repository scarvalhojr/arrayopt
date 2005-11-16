      PROGRAM QAPGLB
C main program for Gilmore Lawler bound

      INTEGER nn
      PARAMETER( nn  = 256)
      DOUBLE PRECISION XS, lowbd, cost
      DOUBLE PRECISION A( nn, nn), B( nn, nn)
      DOUBLE PRECISION C( nn*nn), eps, suma( nn), sumb( nn)
      double precision C1( nn*nn)

      INTEGER INDA(NN,NN),INDB(NN,NN),IWORK1(NN),IWORK2(NN),IPERM(NN)
      DOUBLE PRECISION  WORK1( NN),WORK2( NN), WORK3( NN), WORK4( NN)
      LOGICAL LABEL( NN)

C
      eps = 1.0e-7
C
        call lese( n, a, b, c, nn)
C Zeilensummen von A und B
        do 10 i=1,n
                suma( i) = 0.0
10              sumb( i) = 0.0
        do 20 i=1,n
        do 20 j=1,n
                suma( i) = suma( i) + a( i, j)
                sumb( i) = sumb( i) + b( i, j)
20              continue

c sort columns of A and B in increasing order
        CALL MSORT(N, NN, A, INDA, WORK1, IWORK1)
        CALL MSORT(N, NN, B, INDB, WORK1, IWORK1)

C linear term (minimal scalar products of columns
C from A and B and add to C
        DO 60 I=1,N
        DO 60 J=1,N
c multiply column i from A with column j from B
                XS = 0
                DO 50 K=1,N-1
50              XS = XS + A( INDA(K, I), I) * B( INDB(N-K, J), J)
C add min. scalar product to C and also main diagonal contribution
60      C1( (I-1)*N + J) = C( (I-1)*N + J) + XS + A( I, I) * B( J, J)

C solve linear assignment problem with C1
        CALL LSAPR( N, C1, COST, IWORK1, IPERM, WORK1, WORK2,
     1     WORK3,WORK4, IWORK2, LABEL, EPS)

C lower bound is given by COST
        lowbd = COST

C  check also value of QAP corresponding to opt. permutation of C1
        XS = 0
        DO 65 I=1,N
                XS = XS + C( (I-1)*N + IPERM( I))
                DO 65 J=1,N
65                      XS = XS + A( I, J) * B( IPERM( I), IPERM( J))
        print 991, lowbd, xs
991     format(5x, 'gilmore-lawler:', f20.2,/, 5x,
     *             'feas. solution:',f20.2)

      END


        SUBROUTINE MSORT( N, NMAX, A, INDA, WORK, IWORK)
        DOUBLE PRECISION A( NMAX, NMAX), WORK( NMAX)
        INTEGER INDA( NMAX, NMAX), IWORK( NMAX)
c sort columns of matrix a (without element on main diagonal) in
c increasing order. inda(. ,j) contains order of elements from col. j.
c i.e. a( inda(i,j), j) i=1,...,n-1  produces increasing sequence

        DO 30 J = 1,N
c for all columns
                K = 0
                DO 10 I = 1,N
                        IF (I .EQ. J) GOTO 10
                        K = K + 1
                        WORK( K) = A( I, J)
                        IWORK( K) = I
10                      CONTINUE

                CALL SSORT( WORK, IWORK, N-1)
c store permutation
                DO 20 I = 1,N-1
20                      INDA( I, J) = IWORK( I)
30      CONTINUE
        RETURN
        END

      subroutine ssort (a, b, l)
      implicit integer (a-z)
      dimension b(1)
      double precision a(1), ah
      f=1
      if(l.le.f) return
      n2 =(l-f+1)/2
      s =1023
      do 100 t=1,10
      if(s.gt.n2) goto 90
      ls = l-s
      do 20 i=f,ls
      is = i+s
      ah = a(is)
      bh = b(is)
      j = i
      js = is
    5 if(ah.ge.a(j)) goto 10
      a(js) = a(j)
      b(js) = b(j)
      js =j
      j =j-s
      if(j.ge.f) goto 5
   10 a(js) = ah
      b(js) = bh
   20 continue
   90 s=s/2
  100 continue
      return
      end

      SUBROUTINE LSAPR (N,C,Z,ZEILE,SPALTE,DMINUS,DPLUS,YS,YT,VOR
     *                  ,LABEL,EPS)
C *** ****************************************************************
C     *            LINEAR SUM ASSIGNMENT PROBLEM                     *
C     *                    WITH double precision costs               *
C *** ****************************************************************
C     *      INPUT:                                                  *
C     *         N     DIMENSION OF THE COST MATRIX C                 *
C     *         C(I)  REAL COST MATRIX , ROWWISE STORED  I=1,...,N*N *
C     *         EPS   MACHINE ACCURACY (REAL)                        *
C     *      OUTPUT:                                                 *
C     *         SPALTE(I)   OPTIMAL ASSIGNMENT (I=1,...,N)           *
C     *         Z           OPTIMAL VALUE  (REAL)                    *
C     *         YS(I)    OPTIMAL DUAL (ROW) VARIABLES (REAL)         *
C     *                   I=1,...,N                                  *
C     *         YT(I)    OPTIMAL DUAL (COLUMN) VARIABLES (REAL)      *
C     *                   I=1,...,N                                  *
C     *      REAL ARRAYS OF DIMENSION N                              *
C     *         DMINUS(I), DPLUS(I)                                  *
C     *      INTEGER ARRAYS OF DIMENSION N                           *
C     *         ZEILE(I), VOR(I)                                     *
C     *      LOGICAL ARRAY OF LENGTH N:                              *
C     *        LABEL(I)                                              *
C *** ****************************************************************
C<<
      double precision C(1),YS(1),YT(1),DMINUS(1),DPLUS(1)
      double precision z, eps, sup, cc, ui, vj, d, vgl
      INTEGER ZEILE(1),SPALTE(1),VOR(1)
      INTEGER U,US,USI,W,WS,WSI
      LOGICAL LABEL(1)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** STARTPROCEDURE
C     CONSTRUCTION OF AN INITIAL (PARTIAL) ASSIGNMENT
C
      sup = 1.0d15
      DO 50 I=1,N
      ZEILE(I)=0
      SPALTE(I)=0
      VOR(I)=0
      YS(I)=0.
   50 YT(I)=0.
      IK=0
      DO 2 I=1,N
      DO 3 J=1,N
      IK=IK+1
      CC=C(IK)
      IF(J.EQ.1) GOTO 4
      IF(CC-UI.GE.EPS) GOTO 3
    4 UI=CC
      JO=J
    3 CONTINUE
      YS(I)=UI
      IF(ZEILE(JO).NE.0) GOTO 2
      ZEILE(JO)=I
      SPALTE(I)=JO
    2 CONTINUE
      DO 5 J=1,N
      IF(ZEILE(J).EQ.0) YT(J)=SUP
    5 CONTINUE
      IK=0
      DO 6 I=1,N
      UI=YS(I)
      DO 7 J=1,N
      IK=IK+1
      VJ=YT(J)
      IF(VJ .LE. EPS) GOTO 7
      CC=C(IK)-UI
      IF(CC+EPS .GE. VJ) GOTO 7
      YT(J)=CC
      VOR(J)=I
    7 CONTINUE
    6 CONTINUE
      DO 8 J=1,N
      I=VOR(J)
      IF(I.EQ.0) GOTO 8
      IF(SPALTE(I).NE.0) GOTO 8
      SPALTE(I)=J
      ZEILE(J)=I
    8 CONTINUE
C<<
      DO 9 I=1,N
      IF(SPALTE(I).NE.0) GOTO 9
      UI=YS(I)
      IK=(I-1)*N
      DO 10 J=1,N
      IK=IK+1
      IF(ZEILE(J).NE.0) GOTO 10
      CC=C(IK)
      IF((CC-UI-YT(J)+EPS) .GT. 0.) GOTO 10
      SPALTE(I)=J
      ZEILE(J)=I
      GOTO 9
   10 CONTINUE
    9 CONTINUE
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** CONSTRUCTION OF THE OPTIMAL ASSIGNMENT
C
      DO 1000 U=1,N
      IF(SPALTE(U).GT.0) GOTO 1000
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** SHORTEST PATH COMPUTATION
C
      US=(U-1)*N
      DO 100 I=1,N
      VOR(I)=U
      LABEL(I)=.FALSE.
      DPLUS(I)=SUP
      USI=US+I
  100 DMINUS(I)=C(USI)-YS(U)-YT(I)
      DPLUS(U)=0.
  105 D=SUP
      DO 110 I=1,N
      IF(LABEL(I)) GOTO 110
      IF(DMINUS(I)+EPS .GE. D) GOTO 110
      D=DMINUS(I)
      INDEX=I
  110 CONTINUE
      IF(ZEILE(INDEX).LE.0) GOTO 400
      LABEL(INDEX)=.TRUE.
      W=ZEILE(INDEX)
      WS=(W-1)*N
      DPLUS(W)=D
      DO 130 I=1,N
      IF(LABEL(I)) GOTO 130
      WSI=WS+I
      VGL=D+C(WSI)-YS(W)-YT(I)
      IF(DMINUS(I) .LE. VGL+EPS) GOTO 130
      DMINUS(I)=VGL
      VOR(I)=W
  130 CONTINUE
      GOTO 105
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** AUGMENTATION
C
  400 W=VOR(INDEX)
      ZEILE(INDEX)=W
      IND=SPALTE(W)
      SPALTE(W)=INDEX
      IF(W.EQ.U) GOTO 500
C<<
      INDEX=IND
      GOTO 400
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** TRANSFORMATION
C
  500 DO 510 I=1,N
      IF(DPLUS(I).EQ.SUP) GOTO 505
      YS(I)=YS(I)+D-DPLUS(I)
  505 IF(DMINUS(I)+EPS .GE. D) GOTO 510
      YT(I)=YT(I)+DMINUS(I)-D
  510 CONTINUE
 1000 CONTINUE
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** COMPUTATION OF THE OPTIMAL VALUE
C
      Z=0.
      DO 2000 I=1,N
      IS=(I-1)*N
      J=SPALTE(I)
      ISJ=IS+J
      Z=Z+C(ISJ)
 2000 CONTINUE
      RETURN
      END

        subroutine lese(n, a, b, c, ndim)
        double precision a( ndim, ndim), b( ndim, ndim)
        double precision c(ndim*ndim)

        read (5,*) n
        if (n .gt. ndim) then
                print *,'size too large: n, nmax=',n, ndim
                stop
                        endif
        do 10 i=1, n
10              read (5,*) (a(i, j), j=1, n)
        do 20 i=1, n
20              read (5,*) (b(i, j), j=1, n)

        read (5,*,end=40) (c(i), i = 1, n*n)
        return

40      continue
        print *,'no linear term.'
        do 50 i = 1,n*n
50      c( i) = 0.0
        return
        end
