      SUBROUTINE SLAHD2( IOUNIT, PATH )
*
*  -- LAPACK auxiliary test routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER*3        PATH
      INTEGER            IOUNIT
*     ..
*
*  Purpose
*  =======
*
*  SLAHD2 prints header information for the different test paths.
*
*  Arguments
*  =========
*
*  IOUNIT  (input) INTEGER.
*          On entry, IOUNIT specifies the unit number to which the
*          header information should be printed.
*
*  PATH    (input) CHARACTER*3.
*          On entry, PATH contains the name of the path for which the
*          header information is to be printed.  Current paths are
*
*             SHS, CHS:  Non-symmetric eigenproblem.
*             SST, CST:  Symmetric eigenproblem.
*             SBD, CBD:  Singular Value Decomposition (SVD)
*
*          These paths also are supplied in double precision (replace
*          leading S by D and leading C by Z in path names).
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            CORZ, SORD
      CHARACTER*2        C2
      INTEGER            J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, LSAMEN
      EXTERNAL           LSAME, LSAMEN
*     ..
*     .. Executable Statements ..
*
      IF( IOUNIT.LE.0 )
     $   RETURN
      SORD = LSAME( PATH, 'S' ) .OR. LSAME( PATH, 'D' )
      CORZ = LSAME( PATH, 'C' ) .OR. LSAME( PATH, 'Z' )
      IF( .NOT.SORD .AND. .NOT.CORZ ) THEN
         WRITE( IOUNIT, FMT = 9999 )PATH
 9999    FORMAT( 1X, A3, ':  no header available' )
      END IF
      C2 = PATH( 2: 3 )
*
      IF( LSAMEN( 2, C2, 'HS' ) ) THEN
         IF( SORD ) THEN
*
*           Real Non-symmetric Eigenvalue Problem:
*
            WRITE( IOUNIT, FMT = 9998 )PATH
*
*           Matrix types
*
            WRITE( IOUNIT, FMT = 9992 )
            WRITE( IOUNIT, FMT = 9991 )
            WRITE( IOUNIT, FMT = 9990 )'pairs ', 'pairs ', 'prs.',
     $         'prs.'
            WRITE( IOUNIT, FMT = 9989 )
*
*           Tests performed
*
            WRITE( IOUNIT, FMT = 9988 )'orthogonal', '''=transpose',
     $         ( '''', J = 1, 6 )
*
         ELSE
*
*           Complex Non-symmetric Eigenvalue Problem:
*
            WRITE( IOUNIT, FMT = 9997 )PATH
*
*           Matrix types
*
            WRITE( IOUNIT, FMT = 9992 )
            WRITE( IOUNIT, FMT = 9991 )
            WRITE( IOUNIT, FMT = 9990 )'e.vals', 'e.vals', 'e.vs',
     $         'e.vs'
            WRITE( IOUNIT, FMT = 9989 )
*
*           Tests performed
*
            WRITE( IOUNIT, FMT = 9988 )'unitary', '*=conj.transp.',
     $         ( '*', J = 1, 6 )
         END IF
*
      ELSE IF( LSAMEN( 2, C2, 'ST' ) ) THEN
*
         IF( SORD ) THEN
*
*           Real Symmetric Eigenvalue Problem:
*
            WRITE( IOUNIT, FMT = 9996 )PATH
*
*           Matrix types
*
            WRITE( IOUNIT, FMT = 9987 )
            WRITE( IOUNIT, FMT = 9986 )
            WRITE( IOUNIT, FMT = 9985 )'Symmetric'
*
*           Tests performed
*
            WRITE( IOUNIT, FMT = 9984 )'orthogonal', '''=transpose',
     $         ( '''', J = 1, 6 )
*
         ELSE
*
*           Complex Hermitian Eigenvalue Problem:
*
            WRITE( IOUNIT, FMT = 9995 )PATH
*
*           Matrix types
*
            WRITE( IOUNIT, FMT = 9987 )
            WRITE( IOUNIT, FMT = 9986 )
            WRITE( IOUNIT, FMT = 9985 )'Hermitian'
*
*           Tests performed
*
            WRITE( IOUNIT, FMT = 9984 )'unitary', '*=conj.transp.',
     $         ( '*', J = 1, 6 )
         END IF
*
      ELSE IF( LSAMEN( 2, C2, 'BD' ) ) THEN
*
         IF( SORD ) THEN
*
*           Real Singular Value Decomposition:
*
            WRITE( IOUNIT, FMT = 9994 )PATH
*
*           Matrix types
*
            WRITE( IOUNIT, FMT = 9983 )
*
*           Tests performed
*
            WRITE( IOUNIT, FMT = 9982 )'orthogonal'
         ELSE
*
*           Complex Singular Value Decomposition:
*
            WRITE( IOUNIT, FMT = 9993 )PATH
*
*           Matrix types
*
            WRITE( IOUNIT, FMT = 9983 )
*
*           Tests performed
*
            WRITE( IOUNIT, FMT = 9982 )'unitary   '
         END IF
*
      ELSE
*
         WRITE( IOUNIT, FMT = 9999 )PATH
         RETURN
      END IF
*
      RETURN
*
 9998 FORMAT( / 1X, A3, ' -- Real Non-symmetric eigenvalue problem' )
 9997 FORMAT( / 1X, A3, ' -- Complex Non-symmetric eigenvalue problem' )
 9996 FORMAT( / 1X, A3, ' -- Real Symmetric eigenvalue problem' )
 9995 FORMAT( / 1X, A3, ' -- Complex Hermitian eigenvalue problem' )
 9994 FORMAT( / 1X, A3, ' -- Real Singular Value Decomposition' )
 9993 FORMAT( / 1X, A3, ' -- Complex Singular Value Decomposition' )
*
 9992 FORMAT( ' Matrix types (see xCHKHS for details): ' )
*
 9991 FORMAT( / ' Special Matrices:', / '  1=Zero matrix.             ',
     $      '           ', '  5=Diagonal: geometr. spaced entries.',
     $      / '  2=Identity matrix.                    ', '  6=Diagona',
     $      'l: clustered entries.', / '  3=Transposed Jordan block.  ',
     $      '          ', '  7=Diagonal: large, evenly spaced.', / '  ',
     $      '4=Diagonal: evenly spaced entries.    ', '  8=Diagonal: s',
     $      'mall, evenly spaced.' )
 9990 FORMAT( ' Dense, Non-Symmetric Matrices:', / '  9=Well-cond., ev',
     $      'enly spaced eigenvals.', ' 14=Ill-cond., geomet. spaced e',
     $      'igenals.', / ' 10=Well-cond., geom. spaced eigenvals. ',
     $      ' 15=Ill-conditioned, clustered e.vals.', / ' 11=Well-cond',
     $      'itioned, clustered e.vals. ', ' 16=Ill-cond., random comp',
     $      'lex ', A6, / ' 12=Well-cond., random complex ', A6, '   ',
     $      ' 17=Ill-cond., large rand. complx ', A4, / ' 13=Ill-condi',
     $      'tioned, evenly spaced.     ', ' 18=Ill-cond., small rand.',
     $      ' complx ', A4 )
 9989 FORMAT( ' 19=Matrix with random O(1) entries.    ', ' 21=Matrix ',
     $      'with small random entries.', / ' 20=Matrix with large ran',
     $      'dom entries.   ' )
 9988 FORMAT( / ' Tests performed:   ', '(H is Hessenberg, T is Schur,',
     $      ' U and Z are ', A, ',', / 20X, A, ', W is a diagonal matr',
     $      'ix of eigenvalues,', / 20X, 'L and R are the left and rig',
     $      'ht eigenvector matrices)', / '  1 = | A - U H U', A1, ' |',
     $      ' / ( |A| n ulp )         ', '  2 = | I - U U', A1, ' | / ',
     $      '( n ulp )', / '  3 = | H - Z T Z', A1, ' | / ( |H| n ulp ',
     $      ')         ', '  4 = | I - Z Z', A1, ' | / ( n ulp )',
     $      / '  5 = | A - UZ T (UZ)', A1, ' | / ( |A| n ulp )     ',
     $      '  6 = | I - UZ (UZ)', A1, ' | / ( n ulp )', / '  7 = | T(',
     $      'e.vects.) - T(no e.vects.) | / ( |T| ulp )', / '  8 = | W',
     $      '(e.vects.) - W(no e.vects.) | / ( |W| ulp )', / '  9 = | ',
     $      'TR - RW | / ( |T| |R| ulp )     ', ' 10 = | LT - WL | / (',
     $      ' |T| |L| ulp )', / ' 11= |HX - XW| / (|H| |X| ulp)  (inv.',
     $      'it)', ' 12= |YH - WY| / (|H| |Y| ulp)  (inv.it)' )
*
*     Symmetric/Hermitian eigenproblem
*
 9987 FORMAT( ' Matrix types (see xCHKST for details): ' )
*
 9986 FORMAT( / ' Special Matrices:', / '  1=Zero matrix.             ',
     $      '           ', '  5=Diagonal: clustered entries.', / '  2=',
     $      'Identity matrix.                    ', '  6=Diagonal: lar',
     $      'ge, evenly spaced.', / '  3=Diagonal: evenly spaced entri',
     $      'es.    ', '  7=Diagonal: small, evenly spaced.', / '  4=D',
     $      'iagonal: geometr. spaced entries.' )
 9985 FORMAT( ' Dense ', A, ' Matrices:', / '  8=Evenly spaced eigen',
     $      'vals.            ', ' 12=Small, evenly spaced eigenvals.',
     $      / '  9=Geometrically spaced eigenvals.     ', ' 13=Matrix ',
     $      'with random O(1) entries.', / ' 10=Clustered eigenvalues.',
     $      '              ', ' 14=Matrix with large random entries.',
     $      / ' 11=Large, evenly spaced eigenvals.     ', ' 15=Matrix ',
     $      'with small random entries.' )
*
 9984 FORMAT( / ' Tests performed:   ', '(S is Tridiag, D is diagonal,',
     $      ' U and Z are ', A, ',', / 20X, A, ', W is a diagonal matr',
     $      'ix of eigenvalues)', / '  1= | A - U S U', A1, ' | / ( |A',
     $      '| n ulp )     ', '  2= | I - U U', A1, ' | / ( n ulp )',
     $      / '  3= | S - Z D Z', A1, ' | / ( |S| n ulp )     ', '  4=',
     $      ' | I - Z Z', A1, ' | / ( n ulp )',
     $      / '  5= | A - UZ D (UZ)', A1, ' | / ( |A| n ulp ) ', '  6=',
     $      ' | I - UZ (UZ)', A1, ' | / ( n ulp )', / '  7= |D(with Z)',
     $      ' - D(w/o Z)| / (|D| ulp) ', '  8= | D(PWK) - D(QR) | / (|',
     $      'D| ulp)', / '  9=   Sturm sequence test on W         ',
     $      ' 10= | Z(inv it.) - Z(QR) | / (|Z| ulp)' )
*
*     Singular Value Decomposition
*
 9983 FORMAT( ' Matrix types (see xCHKBD for details):',
     $      / ' Diagonal matrices:', / '   1: Zero', 28X,
     $      ' 5: Clustered entries', / '   2: Identity', 24X,
     $      ' 6: Large, evenly spaced entries',
     $      / '   3: Evenly spaced entries', 11X,
     $      ' 7: Small, evenly spaced entries',
     $      / '   4: Geometrically spaced entries',
     $      / ' General matrices:', / '   8: Evenly spaced sing. vals.',
     $      7X, '12: Small, evenly spaced sing vals',
     $      / '   9: Geometrically spaced sing vals  ',
     $      '13: Random, O(1) entries', / '  10: Clustered sing. vals.',
     $      11X, '14: Random, scaled near overflow',
     $      / '  11: Large, evenly spaced sing vals  ',
     $      '15: Random, scaled near underflow' )
*
 9982 FORMAT( / ' Test ratios:  ',
     $      '(B: bidiagonal, S: diagonal, Q, P, U, and V: ', A10, / 16X,
     $      'X: m x nrhs, Y = Q'' X, and Z = U'' Y)',
     $      / '   1: norm( A - Q B P'' ) / ( norm(A) max(m,n) ulp )',
     $      / '   2: norm( I - Q'' Q )   / ( m ulp )',
     $      / '   3: norm( I - P'' P )   / ( n ulp )',
     $      / '   4: norm( B - U S V'' ) / ( norm(B) min(m,n) ulp )', /
     $      '   5: norm( Y - U Z )    / ( norm(Z) max(min(m,n),k) ulp )'
     $      , / '   6: norm( I - U'' U )   / ( min(m,n) ulp )',
     $      / '   7: norm( I - V'' V )   / ( min(m,n) ulp )',
     $      / '   8: Test ordering of S  (0 if nondecreasing, 1/ulp ',
     $      'otherwise)', / '   9: Sturm sequence test ',
     $      '(0 if sing. vals of B within THRESH of S)',
     $      / '  10: norm( S - S2 )     / ( norm(S) ulp ),',
     $      ' where S2 is computed'/ 44X, 'without computing U and V''',
     $      / '  11: norm( A - (QU) S (V'' P'') ) / ',
     $      '( norm(A) max(m,n) ulp )', /
     $      '  12: norm( X - (QU) Z )         / ( |X| max(M,k) ulp )',
     $      / '  13: norm( I - (QU)''(QU) )      / ( M ulp )',
     $      / '  14: norm( I - (V'' P'') (P V) )  / ( N ulp )' )
*
*     End of SLAHD2
*
      END
