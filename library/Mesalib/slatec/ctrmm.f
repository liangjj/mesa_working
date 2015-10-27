*deck ctrmm
      subroutine ctrmm (side, uplo, transa, diag, m, n, alpha, a, lda,
     $   b, ldb)
c***begin prologue  ctrmm
c***purpose  multiply a complex general matrix by a complex triangular
c            matrix.
c***library   slatec (blas)
c***category  d1b6
c***type      complex (strmm-s, dtrmm-d, ctrmm-c)
c***keywords  level 3 blas, linear algebra
c***author  dongarra, j., (anl)
c           duff, i., (aere)
c           du croz, j., (nag)
c           hammarling, s. (nag)
c***description
c
c  ctrmm  performs one of the matrix-matrix operations
c
c     b := alpha*op( a )*b,   or   b := alpha*b*op( a )
c
c  where  alpha  is a scalar,  b  is an m by n matrix,  a  is a unit, or
c  non-unit,  upper or lower triangular matrix  and  op( a )  is one  of
c
c     op( a ) = a   or   op( a ) = a'   or   op( a ) = conjg( a' ).
c
c  parameters
c  ==========
c
c  side   - character*1.
c           on entry,  side specifies whether  op( a ) multiplies b from
c           the left or right as follows:
c
c              side = 'l' or 'l'   b := alpha*op( a )*b.
c
c              side = 'r' or 'r'   b := alpha*b*op( a ).
c
c           unchanged on exit.
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the matrix a is an upper or
c           lower triangular matrix as follows:
c
c              uplo = 'u' or 'u'   a is an upper triangular matrix.
c
c              uplo = 'l' or 'l'   a is a lower triangular matrix.
c
c           unchanged on exit.
c
c  transa - character*1.
c           on entry, transa specifies the form of op( a ) to be used in
c           the matrix multiplication as follows:
c
c              transa = 'n' or 'n'   op( a ) = a.
c
c              transa = 't' or 't'   op( a ) = a'.
c
c              transa = 'c' or 'c'   op( a ) = conjg( a' ).
c
c           unchanged on exit.
c
c  diag   - character*1.
c           on entry, diag specifies whether or not a is unit triangular
c           as follows:
c
c              diag = 'u' or 'u'   a is assumed to be unit triangular.
c
c              diag = 'n' or 'n'   a is not assumed to be unit
c                                  triangular.
c
c           unchanged on exit.
c
c  m      - integer.
c           on entry, m specifies the number of rows of b. m must be at
c           least zero.
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the number of columns of b.  n must be
c           at least zero.
c           unchanged on exit.
c
c  alpha  - complex         .
c           on entry,  alpha specifies the scalar  alpha. when  alpha is
c           zero then  a is not referenced and  b need not be set before
c           entry.
c           unchanged on exit.
c
c  a      - complex          array of dimension ( lda, k ), where k is m
c           when  side = 'l' or 'l'  and is  n  when  side = 'r' or 'r'.
c           before entry  with  uplo = 'u' or 'u',  the  leading  k by k
c           upper triangular part of the array  a must contain the upper
c           triangular matrix  and the strictly lower triangular part of
c           a is not referenced.
c           before entry  with  uplo = 'l' or 'l',  the  leading  k by k
c           lower triangular part of the array  a must contain the lower
c           triangular matrix  and the strictly upper triangular part of
c           a is not referenced.
c           note that when  diag = 'u' or 'u',  the diagonal elements of
c           a  are not referenced either,  but are assumed to be  unity.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program.  when  side = 'l' or 'l'  then
c           lda  must be at least  max( 1, m ),  when  side = 'r' or 'r'
c           then lda must be at least max( 1, n ).
c           unchanged on exit.
c
c  b      - complex          array of dimension ( ldb, n ).
c           before entry,  the leading  m by n part of the array  b must
c           contain the matrix  b,  and  on exit  is overwritten  by the
c           transformed matrix.
c
c  ldb    - integer.
c           on entry, ldb specifies the first dimension of b as declared
c           in  the  calling  (sub)  program.   ldb  must  be  at  least
c           max( 1, m ).
c           unchanged on exit.
c
c***references  dongarra, j., du croz, j., duff, i., and hammarling, s.
c                 a set of level 3 basic linear algebra subprograms.
c                 acm toms, vol. 16, no. 1, pp. 1-17, march 1990.
c***routines called  lsame, xerbla
c***revision history  (yymmdd)
c   890208  date written
c   910605  modified to meet slatec prologue standards.  only comment
c           lines were modified.  (bks)
c***end prologue  ctrmm
c     .. scalar arguments ..
      character*1        side, uplo, transa, diag
      integer            m, n, lda, ldb
      complex            alpha
c     .. array arguments ..
      complex            a( lda, * ), b( ldb, * )
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg, max
c     .. local scalars ..
      logical            lside, noconj, nounit, upper
      integer            i, info, j, k, nrowa
      complex            temp
c     .. parameters ..
      complex            one
      parameter        ( one  = ( 1.0e+0, 0.0e+0 ) )
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c***first executable statement  ctrmm
c
c     test the input parameters.
c
      lside  = lsame( side  , 'l' )
      if( lside )then
         nrowa = m
      else
         nrowa = n
      end if
      noconj = lsame( transa, 't' )
      nounit = lsame( diag  , 'n' )
      upper  = lsame( uplo  , 'u' )
c
      info   = 0
      if(      ( .not.lside                ).and.
     $         ( .not.lsame( side  , 'r' ) )      )then
         info = 1
      else if( ( .not.upper                ).and.
     $         ( .not.lsame( uplo  , 'l' ) )      )then
         info = 2
      else if( ( .not.lsame( transa, 'n' ) ).and.
     $         ( .not.lsame( transa, 't' ) ).and.
     $         ( .not.lsame( transa, 'c' ) )      )then
         info = 3
      else if( ( .not.lsame( diag  , 'u' ) ).and.
     $         ( .not.lsame( diag  , 'n' ) )      )then
         info = 4
      else if( m  .lt.0               )then
         info = 5
      else if( n  .lt.0               )then
         info = 6
      else if( lda.lt.max( 1, nrowa ) )then
         info = 9
      else if( ldb.lt.max( 1, m     ) )then
         info = 11
      end if
      if( info.ne.0 )then
         call xerbla( 'ctrmm ', info )
         return
      end if
c
c     quick return if possible.
c
      if( n.eq.0 )
     $   return
c
c     and when  alpha.eq.zero.
c
      if( alpha.eq.zero )then
         do 20, j = 1, n
            do 10, i = 1, m
               b( i, j ) = zero
   10       continue
   20    continue
         return
      end if
c
c     start the operations.
c
      if( lside )then
         if( lsame( transa, 'n' ) )then
c
c           form  b := alpha*a*b.
c
            if( upper )then
               do 50, j = 1, n
                  do 40, k = 1, m
                     if( b( k, j ).ne.zero )then
                        temp = alpha*b( k, j )
                        do 30, i = 1, k - 1
                           b( i, j ) = b( i, j ) + temp*a( i, k )
   30                   continue
                        if( nounit )
     $                     temp = temp*a( k, k )
                        b( k, j ) = temp
                     end if
   40             continue
   50          continue
            else
               do 80, j = 1, n
                  do 70 k = m, 1, -1
                     if( b( k, j ).ne.zero )then
                        temp      = alpha*b( k, j )
                        b( k, j ) = temp
                        if( nounit )
     $                     b( k, j ) = b( k, j )*a( k, k )
                        do 60, i = k + 1, m
                           b( i, j ) = b( i, j ) + temp*a( i, k )
   60                   continue
                     end if
   70             continue
   80          continue
            end if
         else
c
c           form  b := alpha*b*a'   or   b := alpha*b*conjg( a' ).
c
            if( upper )then
               do 120, j = 1, n
                  do 110, i = m, 1, -1
                     temp = b( i, j )
                     if( noconj )then
                        if( nounit )
     $                     temp = temp*a( i, i )
                        do 90, k = 1, i - 1
                           temp = temp + a( k, i )*b( k, j )
   90                   continue
                     else
                        if( nounit )
     $                     temp = temp*conjg( a( i, i ) )
                        do 100, k = 1, i - 1
                           temp = temp + conjg( a( k, i ) )*b( k, j )
  100                   continue
                     end if
                     b( i, j ) = alpha*temp
  110             continue
  120          continue
            else
               do 160, j = 1, n
                  do 150, i = 1, m
                     temp = b( i, j )
                     if( noconj )then
                        if( nounit )
     $                     temp = temp*a( i, i )
                        do 130, k = i + 1, m
                           temp = temp + a( k, i )*b( k, j )
  130                   continue
                     else
                        if( nounit )
     $                     temp = temp*conjg( a( i, i ) )
                        do 140, k = i + 1, m
                           temp = temp + conjg( a( k, i ) )*b( k, j )
  140                   continue
                     end if
                     b( i, j ) = alpha*temp
  150             continue
  160          continue
            end if
         end if
      else
         if( lsame( transa, 'n' ) )then
c
c           form  b := alpha*b*a.
c
            if( upper )then
               do 200, j = n, 1, -1
                  temp = alpha
                  if( nounit )
     $               temp = temp*a( j, j )
                  do 170, i = 1, m
                     b( i, j ) = temp*b( i, j )
  170             continue
                  do 190, k = 1, j - 1
                     if( a( k, j ).ne.zero )then
                        temp = alpha*a( k, j )
                        do 180, i = 1, m
                           b( i, j ) = b( i, j ) + temp*b( i, k )
  180                   continue
                     end if
  190             continue
  200          continue
            else
               do 240, j = 1, n
                  temp = alpha
                  if( nounit )
     $               temp = temp*a( j, j )
                  do 210, i = 1, m
                     b( i, j ) = temp*b( i, j )
  210             continue
                  do 230, k = j + 1, n
                     if( a( k, j ).ne.zero )then
                        temp = alpha*a( k, j )
                        do 220, i = 1, m
                           b( i, j ) = b( i, j ) + temp*b( i, k )
  220                   continue
                     end if
  230             continue
  240          continue
            end if
         else
c
c           form  b := alpha*b*a'   or   b := alpha*b*conjg( a' ).
c
            if( upper )then
               do 280, k = 1, n
                  do 260, j = 1, k - 1
                     if( a( j, k ).ne.zero )then
                        if( noconj )then
                           temp = alpha*a( j, k )
                        else
                           temp = alpha*conjg( a( j, k ) )
                        end if
                        do 250, i = 1, m
                           b( i, j ) = b( i, j ) + temp*b( i, k )
  250                   continue
                     end if
  260             continue
                  temp = alpha
                  if( nounit )then
                     if( noconj )then
                        temp = temp*a( k, k )
                     else
                        temp = temp*conjg( a( k, k ) )
                     end if
                  end if
                  if( temp.ne.one )then
                     do 270, i = 1, m
                        b( i, k ) = temp*b( i, k )
  270                continue
                  end if
  280          continue
            else
               do 320, k = n, 1, -1
                  do 300, j = k + 1, n
                     if( a( j, k ).ne.zero )then
                        if( noconj )then
                           temp = alpha*a( j, k )
                        else
                           temp = alpha*conjg( a( j, k ) )
                        end if
                        do 290, i = 1, m
                           b( i, j ) = b( i, j ) + temp*b( i, k )
  290                   continue
                     end if
  300             continue
                  temp = alpha
                  if( nounit )then
                     if( noconj )then
                        temp = temp*a( k, k )
                     else
                        temp = temp*conjg( a( k, k ) )
                     end if
                  end if
                  if( temp.ne.one )then
                     do 310, i = 1, m
                        b( i, k ) = temp*b( i, k )
  310                continue
                  end if
  320          continue
            end if
         end if
      end if
c
      return
c
c     end of ctrmm .
c
      end
