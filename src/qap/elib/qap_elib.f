        program qapeli
        parameter (ndim = 256)
        double precision a(ndim, ndim), b( ndim, ndim), c(ndim*ndim)
        double precision lbd
        call lese( n, a, b, c, ndim)
        call elibd( n, ndim, a, b, c, lbd)
        print 1,lbd
1       format(2x, 'elimination bound:',2x,f12.3)
        end

        subroutine elibd( n, ndim, a, b, c, lbd)
c compute elimination bound for QAP given by A, B (symmetric), C
        integer n, ndim
        double precision a( ndim, ndim), b( ndim, ndim), c(ndim*ndim)
     *    , lbd

c local variables
        parameter ( nn = 256)
        integer iw1( nn), iw2( nn), iw3( nn)
        logical label( nn)
        double precision v( nn, nn-1), a0( nn-1, nn-1), sum, sum1,
     *   b0( nn-1, nn-1), work( nn-1, nn), work1( nn-1, nn), work2( nn)
     *  , r1(nn), r2( nn), alfa( nn), beta( nn), c0( nn*nn), eps

c set up matrix v and transform a and b into a0 and b0
c alfa and beta contain row sms of a and b
        do 10 i=1,n
                sum = 0.0
                sum1 = 0.0
                do 8 j=1,n
                        sum = sum + a( i, j)
                        sum1 = sum1 + b( i, j)
 8                      v( i, j) = 0.0
                alfa( i) = sum
10              beta( i) = sum1

        sum = 0.0
        sum1 = 0.0
        do 20 i=1,n
                sum = sum + alfa( i)
                sum1 = sum1 + beta( i)
                do 20 j=i,n-1
20                      v( i, j) = 1.0 / sqrt( float(j*j + j) )

c constant contribution to lower bound
        lbd = - sum * sum1 / float( n*n)
        do 30 i=1,n-1
30              v( i+1, i) = -float(i) / sqrt( float(i*i + i) )

c first form work=a*v, work1=b*v
        do 50 i=1,n
        do 50 j=1,n-1
                sum = 0.0
                sum1 = 0.0
                do 40 k=1,n
                        sum = sum + a( i, k) * v( k, j)
40                      sum1 = sum1 + b( i, k) * v( k, j)
                work( i, j) = sum
50              work1( i, j) = sum1

c now form a0=v'*work, b0=v'*work1
        do 70 i=1,n-1
        do 70 j=1,n-1
                sum = 0.0
                sum1 = 0.0
                do 60 k=1,n
                        sum = sum + v( k, i)* work( k, j)
60                      sum1 = sum1 + v( k, i) * work1( k, j)
                a0( i, j) = sum
70              b0( i, j) = sum1

c compute eigenvalues of a0 and b0
        n1 = n-1
        nn1 = nn - 1
        call eigen( n1, nn1, a0, b0, r1, r2, work, work1, work2)

c minimal scalar eigenvalue product bound quadratic term
        do 80 i=1,n1
80              lbd = lbd + r1( i) * r2( n1 + 1 - i)

c setup matrix for linear assignment problem to bound linear term
c c0 = (2/n) * alfa(i)*beta(j) + c(i,j)
        sum = 2.0 / float( n)
        do 90 i=1,n
        do 90 j=1,n
90              c0( (i-1)*n+j) = alfa( i)*beta( j)*sum + c((i-1)*n+j)

c solve lsap with costs c0
        isup =2**30
        eps = 1.0e-7
c r1, r2, alfa, beta are used as workspace, sum contains opt. value
        call lsapr( n, isup, c0, sum, iw1, iw2, r1, r2, alfa, beta, iw3,
     *     label, eps)

c compute final value of bound
        lbd = lbd + sum
        return
        end

        subroutine eigen(n, ndim, a, b, r1, r2, u, v, wkspce)
c computes eigenvalues and vectors of symmetric A and B
c r1, u: decomposition of A; r2,v: of B
        double precision a(ndim, 1), b(ndim, 1), u(ndim, 1), v(ndim, 1)
        double precision r1( n), r2( n), wkspce( n)

        ifail=0
        call f02abf(a, ndim, n, r1, u, ndim, wkspce, ifail)
        do 15 i=1,n-1
        do 15 j=i+1,n
 15           a(j,i)=a(i,j)
        call f02abf(b, ndim, n, r2, v, ndim, wkspce, ifail)
        do 25 i=1,n-1
        do 25 j=i+1,n
 25           b(j,i)=b(i,j)
        return
        end

      subroutine f01ajf(n, atol, a, ia, d, e, z, iz)
c
c     mark 5c revised
c     this subroutine reduces the given lower triangle of a
c     symmetric matrix, a, stored in the array a(n,n), to
c     tridiagonal form using householders reduction. the diagonal
c     of the result is stored in the array d(n) and the
c     sub-diagonal in the last n - 1 stores of the array e(n)
c     (with the additional element e(1) = 0). the transformation
c     matrices are accumulated in the array z(n,n). the array
c     a is left unaltered unless the actual parameters
c     corresponding to a and z are identical.
c     1st august 1971
c
      integer i, ia, ii, iz, j1, j, k, l, n
      double precision atol, f, g, h, hh, a(ia,n), d(n), e(n), z(iz,n)
      do 40 i=1,n
         do 20 j=1,i
            z(i,j) = a(i,j)
   20    continue
   40 continue
      if (n.eq.1) go to 280
      do 260 ii=2,n
         i = n - ii + 2
         l = i - 2
         f = z(i,i-1)
         g = 0.0d0
         if (l.eq.0) go to 80
         do 60 k=1,l
            g = g + z(i,k)*z(i,k)
   60    continue
   80    h = g + f*f
c     if g is too small for orthogonality to be
c     guaranteed the transformation is skipped
         if (g.gt.atol) go to 100
         e(i) = f
         h = 0.0d0
         go to 240
  100    l = l + 1
         g = dsqrt(h)
         if (f.ge.0.0d0) g = -g
         e(i) = g
         h = h - f*g
         z(i,i-1) = f - g
         f = 0.0d0
         do 180 j=1,l
            z(j,i) = z(i,j)/h
            g = 0.0d0
c     form element of a*u
            do 120 k=1,j
               g = g + z(j,k)*z(i,k)
  120       continue
            j1 = j + 1
            if (j1.gt.l) go to 160
            do 140 k=j1,l
               g = g + z(k,j)*z(i,k)
  140       continue
c     form element of p
  160       e(j) = g/h
            f = f + g*z(j,i)
  180    continue
c     form k
         hh = f/(h+h)
c     form reduced a
         do 220 j=1,l
            f = z(i,j)
            g = e(j) - hh*f
            e(j) = g
            do 200 k=1,j
               z(j,k) = z(j,k) - f*e(k) - g*z(i,k)
  200       continue
  220    continue
  240    d(i) = h
  260 continue
  280 e(1) = 0.0d0
      d(1) = 0.0d0
c     accumulation of transformation matrices
      do 400 i=1,n
         l = i - 1
         if (d(i).eq.0.0d0) go to 360
         do 340 j=1,l
            g = 0.0d0
            do 300 k=1,l
               g = g + z(i,k)*z(k,j)
  300       continue
            do 320 k=1,l
               z(k,j) = z(k,j) - g*z(k,i)
  320       continue
  340    continue
  360    d(i) = z(i,i)
         z(i,i) = 1.0d0
         if (l.eq.0) go to 400
         do 380 j=1,l
            z(i,j) = 0.0d0
            z(j,i) = 0.0d0
  380    continue
  400 continue
      return
      end
      subroutine f02abf(a, ia, n, r, v, iv, e, ifail)
c     mark 4.5 revised
c     eigenvalues and eigenvectors of a real symmetrix matrix
c     1st august 1971
c
      integer isave, ifail, n, ia, iv
      double precision tol, a(ia,n), r(n), v(iv,n), e(n)
      isave = ifail
      ifail = 1
      tol = 10.0 **(-99.0) * 2.0 ** (52.0)
      call f01ajf(n, tol, a, ia, r, e, v, iv)
      tol = 2.0 ** (-52.0)
      call f02amf(n, tol, r, e, v, iv, ifail)
      if (ifail .eq. 0) return
c error termination
      print *,'ifail=',ifail, ' in f02abf'
      stop
      end
      subroutine f02amf(n, acheps, d, e, z, iz, ifail)
c
c     mark 9 revised. ier-326 (sep 1981).
c     this subroutine finds the eigenvalues and eigenvectors of a
c     tridiagonal matrix, t, given with its diagonal elements in
c     the array d(n) and its sub-diagonal elements in the last n
c     - 1 stores of the array e(n), using ql transformations. the
c     eigenvalues are overwritten on the diagonal elements in the
c     array d in ascending order. the eigenvectors are formed in
c     the array z(n,n), overwriting the accumulated
c     transformations as supplied by the subroutine f01ajf. the
c     subroutine will fail if all eigenvalues take more than 30*n
c     iterations.
c     1st april 1972
c
      integer isave, ifail, n, i, l, j, m, i1, m1, ii, k, iz
      double precision b, f, h, acheps, g, p, r, c, s, d(n), e(n), z(iz,
     *n)
      isave = ifail
      if (n.eq.1) go to 40
      do 20 i=2,n
         e(i-1) = e(i)
   20 continue
   40 e(n) = 0.0d0
      b = 0.0d0
      f = 0.0d0
      j = 30*n
      do 300 l=1,n
         h = acheps*(dabs(d(l))+dabs(e(l)))
         if (b.lt.h) b = h
c     look for small sub-diag element
         do 60 m=l,n
            if (dabs(e(m)).le.b) go to 80
   60    continue
   80    if (m.eq.l) go to 280
  100    if (j.le.0) go to 400
         j = j - 1
c     form shift
         g = d(l)
         h = d(l+1) - g
         if (dabs(h).ge.dabs(e(l))) go to 120
         p = h*0.5d0/e(l)
         r = dsqrt(p*p+1.0d0)
         h = p + r
         if (p.lt.0.0d0) h = p - r
         d(l) = e(l)/h
         go to 140
  120    p = 2.0d0*e(l)/h
         r = dsqrt(p*p+1.0d0)
         d(l) = e(l)*p/(1.0d0+r)
  140    h = g - d(l)
         i1 = l + 1
         if (i1.gt.n) go to 180
         do 160 i=i1,n
            d(i) = d(i) - h
  160    continue
  180    f = f + h
c     ql transformation
         p = d(m)
         c = 1.0d0
         s = 0.0d0
         m1 = m - 1
         do 260 ii=l,m1
            i = m1 - ii + l
            g = c*e(i)
            h = c*p
            if (dabs(p).lt.dabs(e(i))) go to 200
            c = e(i)/p
            r = dsqrt(c*c+1.0d0)
            e(i+1) = s*p*r
            s = c/r
            c = 1.0d0/r
            go to 220
  200       c = p/e(i)
            r = dsqrt(c*c+1.0d0)
            e(i+1) = s*e(i)*r
            s = 1.0d0/r
            c = c/r
  220       p = c*d(i) - s*g
            d(i+1) = h + s*(c*g+s*d(i))
c     form vector
            do 240 k=1,n
               h = z(k,i+1)
               z(k,i+1) = s*z(k,i) + c*h
               z(k,i) = c*z(k,i) - s*h
  240       continue
  260    continue
         e(l) = s*p
         d(l) = c*p
         if (dabs(e(l)).gt.b) go to 100
  280    d(l) = d(l) + f
  300 continue
c     order eigenvalues and eigenvectors
      do 380 i=1,n
         k = i
         p = d(i)
         i1 = i + 1
         if (i1.gt.n) go to 340
         do 320 j=i1,n
            if (d(j).ge.p) go to 320
            k = j
            p = d(j)
  320    continue
  340    if (k.eq.i) go to 380
         d(k) = d(i)
         d(i) = p
         do 360 j=1,n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  360    continue
  380 continue
      ifail = 0
      return
  400 ifail = 1
      return
      end
      SUBROUTINE LSAPR (N,ISUP,C,Z,ZEILE,SPALTE,DMINUS,DPLUS,YS,YT,VOR
     *                  ,LABEL,EPS)
C *** ****************************************************************
C     *                                                              *
C     *            LINEAR SUM ASSIGNMENT PROBLEM                     *
C     *                    WITH REAL COSTS                           *
C     *                                                              *
C *** ****************************************************************
C     *                                                              *
C     *  1.  CALL:                                                   *
C     *      CALL LSAPR (N,ISUP,C,Z,ZEILE,SPALTE,DMINUS,DPLUS,YS,YT, *
C     *                  VOR,LABEL,EPS)                              *
C     *                                                              *
C     *  2.  COMPUTER CODE:                                          *
C     *      FORTRAN IV                                              *
C     *                                                              *
C     *  3.  METHOD:                                                 *
C     *      SHORTEST AUGMENTING PATH METHOD                         *
C     *                                                              *
C     *  4.  PARAMETERS:                                             *
C     *      INPUT:                                                  *
C     *         N     DIMENSION OF THE COST MATRIX C                 *
C     *         ISUP  LARGE MACHINE NUMBER (INTEGER)                 *
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
C     *                                                              *
C     *  5.  EXTERNAL SUBROUTINES:                                   *
C     *      NONE                                                    *
C     *                                                              *
C     *  6.  AUTHOR:                                                 *
C     *      U.DERIGS                                                *
C     *                                                              *
C *** ****************************************************************
C
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
      SUP = FLOAT(ISUP)
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
