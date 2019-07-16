*> \brief \b ZUNCSD
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZUNCSD + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zuncsd.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zuncsd.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zuncsd.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       RECURSIVE SUBROUTINE ZUNCSD( JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS,
*                                    SIGNS, M, P, Q, X11, LDX11, X12,
*                                    LDX12, X21, LDX21, X22, LDX22, THETA,
*                                    U1, LDU1, U2, LDU2, V1T, LDV1T, V2T,
*                                    LDV2T, WORK, LWORK, RWORK, LRWORK,
*                                    IWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBU1, JOBU2, JOBV1T, JOBV2T, SIGNS, TRANS
*       INTEGER            INFO, LDU1, LDU2, LDV1T, LDV2T, LDX11, LDX12,
*      $                   LDX21, LDX22, LRWORK, LWORK, M, P, Q
*       ..
*       .. Array Arguments ..
*       INTEGER            IWORK( * )
*       DOUBLE PRECISION   THETA( * )
*       DOUBLE PRECISION   RWORK( * )
*       COMPLEX*16         U1( LDU1, * ), U2( LDU2, * ), V1T( LDV1T, * ),
*      $                   V2T( LDV2T, * ), WORK( * ), X11( LDX11, * ),
*      $                   X12( LDX12, * ), X21( LDX21, * ), X22( LDX22,
*      $                   * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZUNCSD computes the CS decomposition of an M-by-M partitioned
*> unitary matrix X:
*>
*>                                 [  I  0  0 |  0  0  0 ]
*>                                 [  0  C  0 |  0 -S  0 ]
*>     [ X11 | X12 ]   [ U1 |    ] [  0  0  0 |  0  0 -I ] [ V1 |    ]**H
*> X = [-----------] = [---------] [---------------------] [---------]   .
*>     [ X21 | X22 ]   [    | U2 ] [  0  0  0 |  I  0  0 ] [    | V2 ]
*>                                 [  0  S  0 |  0  C  0 ]
*>                                 [  0  0  I |  0  0  0 ]
*>
*> X11 is P-by-Q. The unitary matrices U1, U2, V1, and V2 are P-by-P,
*> (M-P)-by-(M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. C and S are
*> R-by-R nonnegative diagonal matrices satisfying C^2 + S^2 = I, in
*> which R = MIN(P,M-P,Q,M-Q).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOBU1
*> \verbatim
*>          JOBU1 is CHARACTER
*>          = 'Y':      U1 is computed;
*>          otherwise:  U1 is not computed.
*> \endverbatim
*>
*> \param[in] JOBU2
*> \verbatim
*>          JOBU2 is CHARACTER
*>          = 'Y':      U2 is computed;
*>          otherwise:  U2 is not computed.
*> \endverbatim
*>
*> \param[in] JOBV1T
*> \verbatim
*>          JOBV1T is CHARACTER
*>          = 'Y':      V1T is computed;
*>          otherwise:  V1T is not computed.
*> \endverbatim
*>
*> \param[in] JOBV2T
*> \verbatim
*>          JOBV2T is CHARACTER
*>          = 'Y':      V2T is computed;
*>          otherwise:  V2T is not computed.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER
*>          = 'T':      X, U1, U2, V1T, and V2T are stored in row-major
*>                      order;
*>          otherwise:  X, U1, U2, V1T, and V2T are stored in column-
*>                      major order.
*> \endverbatim
*>
*> \param[in] SIGNS
*> \verbatim
*>          SIGNS is CHARACTER
*>          = 'O':      The lower-left block is made nonpositive (the
*>                      "other" convention);
*>          otherwise:  The upper-right block is made nonpositive (the
*>                      "default" convention).
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows and columns in X.
*> \endverbatim
*>
*> \param[in] P
*> \verbatim
*>          P is INTEGER
*>          The number of rows in X11 and X12. 0 <= P <= M.
*> \endverbatim
*>
*> \param[in] Q
*> \verbatim
*>          Q is INTEGER
*>          The number of columns in X11 and X21. 0 <= Q <= M.
*> \endverbatim
*>
*> \param[in,out] X11
*> \verbatim
*>          X11 is COMPLEX*16 array, dimension (LDX11,Q)
*>          On entry, part of the unitary matrix whose CSD is desired.
*> \endverbatim
*>
*> \param[in] LDX11
*> \verbatim
*>          LDX11 is INTEGER
*>          The leading dimension of X11. LDX11 >= MAX(1,P).
*> \endverbatim
*>
*> \param[in,out] X12
*> \verbatim
*>          X12 is COMPLEX*16 array, dimension (LDX12,M-Q)
*>          On entry, part of the unitary matrix whose CSD is desired.
*> \endverbatim
*>
*> \param[in] LDX12
*> \verbatim
*>          LDX12 is INTEGER
*>          The leading dimension of X12. LDX12 >= MAX(1,P).
*> \endverbatim
*>
*> \param[in,out] X21
*> \verbatim
*>          X21 is COMPLEX*16 array, dimension (LDX21,Q)
*>          On entry, part of the unitary matrix whose CSD is desired.
*> \endverbatim
*>
*> \param[in] LDX21
*> \verbatim
*>          LDX21 is INTEGER
*>          The leading dimension of X11. LDX21 >= MAX(1,M-P).
*> \endverbatim
*>
*> \param[in,out] X22
*> \verbatim
*>          X22 is COMPLEX*16 array, dimension (LDX22,M-Q)
*>          On entry, part of the unitary matrix whose CSD is desired.
*> \endverbatim
*>
*> \param[in] LDX22
*> \verbatim
*>          LDX22 is INTEGER
*>          The leading dimension of X11. LDX22 >= MAX(1,M-P).
*> \endverbatim
*>
*> \param[out] THETA
*> \verbatim
*>          THETA is DOUBLE PRECISION array, dimension (R), in which R =
*>          MIN(P,M-P,Q,M-Q).
*>          C = DIAG( COS(THETA(1)), ... , COS(THETA(R)) ) and
*>          S = DIAG( SIN(THETA(1)), ... , SIN(THETA(R)) ).
*> \endverbatim
*>
*> \param[out] U1
*> \verbatim
*>          U1 is COMPLEX*16 array, dimension (LDU1,P)
*>          If JOBU1 = 'Y', U1 contains the P-by-P unitary matrix U1.
*> \endverbatim
*>
*> \param[in] LDU1
*> \verbatim
*>          LDU1 is INTEGER
*>          The leading dimension of U1. If JOBU1 = 'Y', LDU1 >=
*>          MAX(1,P).
*> \endverbatim
*>
*> \param[out] U2
*> \verbatim
*>          U2 is COMPLEX*16 array, dimension (LDU2,M-P)
*>          If JOBU2 = 'Y', U2 contains the (M-P)-by-(M-P) unitary
*>          matrix U2.
*> \endverbatim
*>
*> \param[in] LDU2
*> \verbatim
*>          LDU2 is INTEGER
*>          The leading dimension of U2. If JOBU2 = 'Y', LDU2 >=
*>          MAX(1,M-P).
*> \endverbatim
*>
*> \param[out] V1T
*> \verbatim
*>          V1T is COMPLEX*16 array, dimension (LDV1T,Q)
*>          If JOBV1T = 'Y', V1T contains the Q-by-Q matrix unitary
*>          matrix V1**H.
*> \endverbatim
*>
*> \param[in] LDV1T
*> \verbatim
*>          LDV1T is INTEGER
*>          The leading dimension of V1T. If JOBV1T = 'Y', LDV1T >=
*>          MAX(1,Q).
*> \endverbatim
*>
*> \param[out] V2T
*> \verbatim
*>          V2T is COMPLEX*16 array, dimension (LDV2T,M-Q)
*>          If JOBV2T = 'Y', V2T contains the (M-Q)-by-(M-Q) unitary
*>          matrix V2**H.
*> \endverbatim
*>
*> \param[in] LDV2T
*> \verbatim
*>          LDV2T is INTEGER
*>          The leading dimension of V2T. If JOBV2T = 'Y', LDV2T >=
*>          MAX(1,M-Q).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the work array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is DOUBLE PRECISION array, dimension MAX(1,LRWORK)
*>          On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK.
*>          If INFO > 0 on exit, RWORK(2:R) contains the values PHI(1),
*>          ..., PHI(R-1) that, together with THETA(1), ..., THETA(R),
*>          define the matrix in intermediate bidiagonal-block form
*>          remaining after nonconvergence. INFO specifies the number
*>          of nonzero PHI's.
*> \endverbatim
*>
*> \param[in] LRWORK
*> \verbatim
*>          LRWORK is INTEGER
*>          The dimension of the array RWORK.
*>
*>          If LRWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the RWORK array, returns
*>          this value as the first entry of the work array, and no error
*>          message related to LRWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (M-MIN(P,M-P,Q,M-Q))
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          > 0:  ZBBCSD did not converge. See the description of RWORK
*>                above for details.
*> \endverbatim
*
*> \par References:
*  ================
*>
*>  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer.
*>      Algorithms, 50(1):33-65, 2009.
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date June 2017
*
*> \ingroup complex16OTHERcomputational
*
*  =====================================================================
      RECURSIVE SUBROUTINE zuncsd( JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS,
     $                             SIGNS, M, P, Q, X11, LDX11, X12,
     $                             LDX12, X21, LDX21, X22, LDX22, THETA,
     $                             U1, LDU1, U2, LDU2, V1T, LDV1T, V2T,
     $                             LDV2T, WORK, LWORK, RWORK, LRWORK,
     $                             IWORK, INFO )
*
*  -- LAPACK computational routine (version 3.7.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2017
*
*     .. Scalar Arguments ..
      CHARACTER          JOBU1, JOBU2, JOBV1T, JOBV2T, SIGNS, TRANS
      INTEGER            INFO, LDU1, LDU2, LDV1T, LDV2T, LDX11, LDX12,
     $                   ldx21, ldx22, lrwork, lwork, m, p, q
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   THETA( * )
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         U1( ldu1, * ), U2( ldu2, * ), V1T( ldv1t, * ),
     $                   v2t( ldv2t, * ), work( * ), x11( ldx11, * ),
     $                   x12( ldx12, * ), x21( ldx21, * ), x22( ldx22,
     $                   * )
*     ..
*
*  ===================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      parameter( one = (1.0d0,0.0d0),
     $                     zero = (0.0d0,0.0d0) )
*     ..
*     .. Local Scalars ..
      CHARACTER          TRANST, SIGNST
      INTEGER            CHILDINFO, I, IB11D, IB11E, IB12D, IB12E,
     $                   ib21d, ib21e, ib22d, ib22e, ibbcsd, iorbdb,
     $                   iorglq, iorgqr, iphi, itaup1, itaup2, itauq1,
     $                   itauq2, j, lbbcsdwork, lbbcsdworkmin,
     $                   lbbcsdworkopt, lorbdbwork, lorbdbworkmin,
     $                   lorbdbworkopt, lorglqwork, lorglqworkmin,
     $                   lorglqworkopt, lorgqrwork, lorgqrworkmin,
     $                   lorgqrworkopt, lworkmin, lworkopt, p1, q1
      LOGICAL            COLMAJOR, DEFAULTSIGNS, LQUERY, WANTU1, WANTU2,
     $                   wantv1t, wantv2t
      INTEGER            LRWORKMIN, LRWORKOPT
      LOGICAL            LRQUERY
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zbbcsd, zlacpy, zlapmr, zlapmt,
     $                   zunbdb, zunglq, zungqr
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. Intrinsic Functions
      INTRINSIC          int, max, min
*     ..
*     .. Executable Statements ..
*
*     Test input arguments
*
      info = 0
      wantu1 = lsame( jobu1, 'Y' )
      wantu2 = lsame( jobu2, 'Y' )
      wantv1t = lsame( jobv1t, 'Y' )
      wantv2t = lsame( jobv2t, 'Y' )
      colmajor = .NOT. lsame( trans, 'T' )
      defaultsigns = .NOT. lsame( signs, 'O' )
      lquery = lwork .EQ. -1
      lrquery = lrwork .EQ. -1
      IF( m .LT. 0 ) THEN
         info = -7
      ELSE IF( p .LT. 0 .OR. p .GT. m ) THEN
         info = -8
      ELSE IF( q .LT. 0 .OR. q .GT. m ) THEN
         info = -9
      ELSE IF ( colmajor .AND.  ldx11 .LT. max( 1, p ) ) THEN
        info = -11
      ELSE IF (.NOT. colmajor .AND. ldx11 .LT. max( 1, q ) ) THEN
        info = -11
      ELSE IF (colmajor .AND. ldx12 .LT. max( 1, p ) ) THEN
        info = -13
      ELSE IF (.NOT. colmajor .AND. ldx12 .LT. max( 1, m-q ) ) THEN
        info = -13
      ELSE IF (colmajor .AND. ldx21 .LT. max( 1, m-p ) ) THEN
        info = -15
      ELSE IF (.NOT. colmajor .AND. ldx21 .LT. max( 1, q ) ) THEN
        info = -15
      ELSE IF (colmajor .AND. ldx22 .LT. max( 1, m-p ) ) THEN
        info = -17
      ELSE IF (.NOT. colmajor .AND. ldx22 .LT. max( 1, m-q ) ) THEN
        info = -17
      ELSE IF( wantu1 .AND. ldu1 .LT. p ) THEN
         info = -20
      ELSE IF( wantu2 .AND. ldu2 .LT. m-p ) THEN
         info = -22
      ELSE IF( wantv1t .AND. ldv1t .LT. q ) THEN
         info = -24
      ELSE IF( wantv2t .AND. ldv2t .LT. m-q ) THEN
         info = -26
      END IF
*
*     Work with transpose if convenient
*
      IF( info .EQ. 0 .AND. min( p, m-p ) .LT. min( q, m-q ) ) THEN
         IF( colmajor ) THEN
            transt = 'T'
         ELSE
            transt = 'N'
         END IF
         IF( defaultsigns ) THEN
            signst = 'O'
         ELSE
            signst = 'D'
         END IF
         CALL zuncsd( jobv1t, jobv2t, jobu1, jobu2, transt, signst, m,
     $                q, p, x11, ldx11, x21, ldx21, x12, ldx12, x22,
     $                ldx22, theta, v1t, ldv1t, v2t, ldv2t, u1, ldu1,
     $                u2, ldu2, work, lwork, rwork, lrwork, iwork,
     $                info )
         RETURN
      END IF
*
*     Work with permutation [ 0 I; I 0 ] * X * [ 0 I; I 0 ] if
*     convenient
*
      IF( info .EQ. 0 .AND. m-q .LT. q ) THEN
         IF( defaultsigns ) THEN
            signst = 'O'
         ELSE
            signst = 'D'
         END IF
         CALL zuncsd( jobu2, jobu1, jobv2t, jobv1t, trans, signst, m,
     $                m-p, m-q, x22, ldx22, x21, ldx21, x12, ldx12, x11,
     $                ldx11, theta, u2, ldu2, u1, ldu1, v2t, ldv2t, v1t,
     $                ldv1t, work, lwork, rwork, lrwork, iwork, info )
         RETURN
      END IF
*
*     Compute workspace
*
      IF( info .EQ. 0 ) THEN
*
*        Real workspace
*
         iphi = 2
         ib11d = iphi + max( 1, q - 1 )
         ib11e = ib11d + max( 1, q )
         ib12d = ib11e + max( 1, q - 1 )
         ib12e = ib12d + max( 1, q )
         ib21d = ib12e + max( 1, q - 1 )
         ib21e = ib21d + max( 1, q )
         ib22d = ib21e + max( 1, q - 1 )
         ib22e = ib22d + max( 1, q )
         ibbcsd = ib22e + max( 1, q - 1 )
         CALL zbbcsd( jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q,
     $                theta, theta, u1, ldu1, u2, ldu2, v1t, ldv1t,
     $                v2t, ldv2t, theta, theta, theta, theta, theta,
     $                theta, theta, theta, rwork, -1, childinfo )
         lbbcsdworkopt = int( rwork(1) )
         lbbcsdworkmin = lbbcsdworkopt
         lrworkopt = ibbcsd + lbbcsdworkopt - 1
         lrworkmin = ibbcsd + lbbcsdworkmin - 1
         rwork(1) = lrworkopt
*
*        Complex workspace
*
         itaup1 = 2
         itaup2 = itaup1 + max( 1, p )
         itauq1 = itaup2 + max( 1, m - p )
         itauq2 = itauq1 + max( 1, q )
         iorgqr = itauq2 + max( 1, m - q )
         CALL zungqr( m-q, m-q, m-q, u1, max(1,m-q), u1, work, -1,
     $                childinfo )
         lorgqrworkopt = int( work(1) )
         lorgqrworkmin = max( 1, m - q )
         iorglq = itauq2 + max( 1, m - q )
         CALL zunglq( m-q, m-q, m-q, u1, max(1,m-q), u1, work, -1,
     $                childinfo )
         lorglqworkopt = int( work(1) )
         lorglqworkmin = max( 1, m - q )
         iorbdb = itauq2 + max( 1, m - q )
         CALL zunbdb( trans, signs, m, p, q, x11, ldx11, x12, ldx12,
     $                x21, ldx21, x22, ldx22, theta, theta, u1, u2,
     $                v1t, v2t, work, -1, childinfo )
         lorbdbworkopt = int( work(1) )
         lorbdbworkmin = lorbdbworkopt
         lworkopt = max( iorgqr + lorgqrworkopt, iorglq + lorglqworkopt,
     $              iorbdb + lorbdbworkopt ) - 1
         lworkmin = max( iorgqr + lorgqrworkmin, iorglq + lorglqworkmin,
     $              iorbdb + lorbdbworkmin ) - 1
         work(1) = max(lworkopt,lworkmin)
*
         IF( lwork .LT. lworkmin
     $       .AND. .NOT. ( lquery .OR. lrquery ) ) THEN
            info = -22
         ELSE IF( lrwork .LT. lrworkmin
     $            .AND. .NOT. ( lquery .OR. lrquery ) ) THEN
            info = -24
         ELSE
            lorgqrwork = lwork - iorgqr + 1
            lorglqwork = lwork - iorglq + 1
            lorbdbwork = lwork - iorbdb + 1
            lbbcsdwork = lrwork - ibbcsd + 1
         END IF
      END IF
*
*     Abort if any illegal arguments
*
      IF( info .NE. 0 ) THEN
         CALL xerbla( 'ZUNCSD', -info )
         RETURN
      ELSE IF( lquery .OR. lrquery ) THEN
         RETURN
      END IF
*
*     Transform to bidiagonal block form
*
      CALL zunbdb( trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21,
     $             ldx21, x22, ldx22, theta, rwork(iphi), work(itaup1),
     $             work(itaup2), work(itauq1), work(itauq2),
     $             work(iorbdb), lorbdbwork, childinfo )
*
*     Accumulate Householder reflectors
*
      IF( colmajor ) THEN
         IF( wantu1 .AND. p .GT. 0 ) THEN
            CALL zlacpy( 'L', p, q, x11, ldx11, u1, ldu1 )
            CALL zungqr( p, p, q, u1, ldu1, work(itaup1), work(iorgqr),
     $                   lorgqrwork, info)
         END IF
         IF( wantu2 .AND. m-p .GT. 0 ) THEN
            CALL zlacpy( 'L', m-p, q, x21, ldx21, u2, ldu2 )
            CALL zungqr( m-p, m-p, q, u2, ldu2, work(itaup2),
     $                   work(iorgqr), lorgqrwork, info )
         END IF
         IF( wantv1t .AND. q .GT. 0 ) THEN
            CALL zlacpy( 'U', q-1, q-1, x11(1,2), ldx11, v1t(2,2),
     $                   ldv1t )
            v1t(1, 1) = one
            DO j = 2, q
               v1t(1,j) = zero
               v1t(j,1) = zero
            END DO
            CALL zunglq( q-1, q-1, q-1, v1t(2,2), ldv1t, work(itauq1),
     $                   work(iorglq), lorglqwork, info )
         END IF
         IF( wantv2t .AND. m-q .GT. 0 ) THEN
            CALL zlacpy( 'U', p, m-q, x12, ldx12, v2t, ldv2t )
            IF( m-p .GT. q) THEN
               CALL zlacpy( 'U', m-p-q, m-p-q, x22(q+1,p+1), ldx22,
     $                      v2t(p+1,p+1), ldv2t )
            END IF
            IF( m .GT. q ) THEN
               CALL zunglq( m-q, m-q, m-q, v2t, ldv2t, work(itauq2),
     $                      work(iorglq), lorglqwork, info )
            END IF
         END IF
      ELSE
         IF( wantu1 .AND. p .GT. 0 ) THEN
            CALL zlacpy( 'U', q, p, x11, ldx11, u1, ldu1 )
            CALL zunglq( p, p, q, u1, ldu1, work(itaup1), work(iorglq),
     $                   lorglqwork, info)
         END IF
         IF( wantu2 .AND. m-p .GT. 0 ) THEN
            CALL zlacpy( 'U', q, m-p, x21, ldx21, u2, ldu2 )
            CALL zunglq( m-p, m-p, q, u2, ldu2, work(itaup2),
     $                   work(iorglq), lorglqwork, info )
         END IF
         IF( wantv1t .AND. q .GT. 0 ) THEN
            CALL zlacpy( 'L', q-1, q-1, x11(2,1), ldx11, v1t(2,2),
     $                   ldv1t )
            v1t(1, 1) = one
            DO j = 2, q
               v1t(1,j) = zero
               v1t(j,1) = zero
            END DO
            CALL zungqr( q-1, q-1, q-1, v1t(2,2), ldv1t, work(itauq1),
     $                   work(iorgqr), lorgqrwork, info )
         END IF
         IF( wantv2t .AND. m-q .GT. 0 ) THEN
            p1 = min( p+1, m )
            q1 = min( q+1, m )
            CALL zlacpy( 'L', m-q, p, x12, ldx12, v2t, ldv2t )
            IF( m .GT. p+q ) THEN
               CALL zlacpy( 'L', m-p-q, m-p-q, x22(p1,q1), ldx22,
     $                      v2t(p+1,p+1), ldv2t )
            END IF
            CALL zungqr( m-q, m-q, m-q, v2t, ldv2t, work(itauq2),
     $                   work(iorgqr), lorgqrwork, info )
         END IF
      END IF
*
*     Compute the CSD of the matrix in bidiagonal-block form
*
      CALL zbbcsd( jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, theta,
     $             rwork(iphi), u1, ldu1, u2, ldu2, v1t, ldv1t, v2t,
     $             ldv2t, rwork(ib11d), rwork(ib11e), rwork(ib12d),
     $             rwork(ib12e), rwork(ib21d), rwork(ib21e),
     $             rwork(ib22d), rwork(ib22e), rwork(ibbcsd),
     $             lbbcsdwork, info )
*
*     Permute rows and columns to place identity submatrices in top-
*     left corner of (1,1)-block and/or bottom-right corner of (1,2)-
*     block and/or bottom-right corner of (2,1)-block and/or top-left
*     corner of (2,2)-block
*
      IF( q .GT. 0 .AND. wantu2 ) THEN
         DO i = 1, q
            iwork(i) = m - p - q + i
         END DO
         DO i = q + 1, m - p
            iwork(i) = i - q
         END DO
         IF( colmajor ) THEN
            CALL zlapmt( .false., m-p, m-p, u2, ldu2, iwork )
         ELSE
            CALL zlapmr( .false., m-p, m-p, u2, ldu2, iwork )
         END IF
      END IF
      IF( m .GT. 0 .AND. wantv2t ) THEN
         DO i = 1, p
            iwork(i) = m - p - q + i
         END DO
         DO i = p + 1, m - q
            iwork(i) = i - p
         END DO
         IF( .NOT. colmajor ) THEN
            CALL zlapmt( .false., m-q, m-q, v2t, ldv2t, iwork )
         ELSE
            CALL zlapmr( .false., m-q, m-q, v2t, ldv2t, iwork )
         END IF
      END IF
*
      RETURN
*
*     End ZUNCSD
*
      END
