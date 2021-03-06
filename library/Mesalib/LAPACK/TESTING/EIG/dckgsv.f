      SUBROUTINE DCKGSV( NM, MVAL, PVAL, NVAL, NMATS, ISEED, THRESH,
     $                   NMAX, A, AF, B, BF, U, V, Q, ALPHA, BETA, R,
     $                   IWORK, WORK, RWORK, NIN, NOUT, INFO )
*
*  -- LAPACK test routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      INTEGER            INFO, NIN, NM, NMATS, NMAX, NOUT
      DOUBLE PRECISION   THRESH
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 ), IWORK( * ), MVAL( * ), NVAL( * ),
     $                   PVAL( * )
      DOUBLE PRECISION   A( * ), AF( * ), ALPHA( * ), B( * ), BETA( * ),
     $                   BF( * ), Q( * ), R( * ), RWORK( * ), U( * ),
     $                   V( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DCKGSV tests DGGSVD:
*         the GSVD for M-by-N matrix A and P-by-N matrix B.
*
*  Arguments
*  =========
*
*  NM      (input) INTEGER
*          The number of values of M contained in the vector MVAL.
*
*  MVAL    (input) INTEGER array, dimension (NM)
*          The values of the matrix row dimension M.
*
*  PVAL    (input) INTEGER array, dimension (NP)
*          The values of the matrix row dimension P.
*
*  NVAL    (input) INTEGER array, dimension (NN)
*          The values of the matrix column dimension N.
*
*  NMATS   (input) INTEGER
*          The number of matrix types to be tested for each combination
*          of matrix dimensions.  If NMATS >= NTYPES (the maximum
*          number of matrix types), then all the different types are
*          generated for testing.  If NMATS < NTYPES, another input line
*          is read to get the numbers of the matrix types to be used.
*
*  ISEED   (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator.  The array
*          elements should be between 0 and 4095, otherwise they will be
*          reduced mod 4096, and ISEED(4) must be odd.
*          On exit, the next seed in the random number sequence after
*          all the test matrices have been generated.
*
*  THRESH  (input) DOUBLE PRECISION
*          The threshold value for the test ratios.  A result is
*          included in the output file if RESULT >= THRESH.  To have
*          every test ratio printed, use THRESH = 0.
*
*  NMAX    (input) INTEGER
*          The maximum value permitted for M or N, used in dimensioning
*          the work arrays.
*
*  A       (workspace) DOUBLE PRECISION array, dimension (NMAX*NMAX)
*
*  AF      (workspace) DOUBLE PRECISION array, dimension (NMAX*NMAX)
*
*  B       (workspace) DOUBLE PRECISION array, dimension (NMAX*NMAX)
*
*  BF      (workspace) DOUBLE PRECISION array, dimension (NMAX*NMAX)
*
*  U       (workspace) DOUBLE PRECISION array, dimension (NMAX*NMAX)
*
*  V       (workspace) DOUBLE PRECISION array, dimension (NMAX*NMAX)
*
*  Q       (workspace) DOUBLE PRECISION array, dimension (NMAX*NMAX)
*
*  ALPHA   (workspace) DOUBLE PRECISION array, di €€  [  b@ ж   U €€   ј  l ~`  U €€  ]  l@ ж(  U €€  F  x Џ0  U €€  _  x@ ж0  U €€  8  Ъ ~x  U €€  a  Ъ@ ж@  U €€  J  Ђ №ƒ  U €€  c  Ђ@ жH  U €€  ^  Є@ ЁЎ  U €€  P  њ    U €€  f  њ@ жP  U €€  D  …   U €€  h  …@ жX  U €€  i  џ@ №   U €€  H  в и8  U €€  k  в@ ж`  U €€  b  т ~Р  U €€  m  т@ жh  U €€  @   ~®  U €€  o  @ жp  U €€  Z  ' ~ј  U €€  q  '@ жx  U €€  p  3 wј  U €€  s  3@ жА  U €€  e  | ј   U €€  u  |@ жШ  U €€  j  С   U €€  w  С@ ж†  U €€  `  І    U €€  y  І@ ж®  U €€  n  є 8  U €€  {  є@ ж∞  U €€  z  – P  U €€  }  –@ жЄ  U €€  ~  з ј†  U €€    з@ жј  U €€  t  т ≈H  U €€  Б  т@ ж»  U €€  X  ъ ƒ  U €€  Г  ъ@ ж–  U €€  B   h  U €€  Е  @ жЎ  U €€  r   А  U €€  З  @ жа  U €€  |  < Ш  U €€  Й  <@ жи  U €€  V  O ∞  U €€  Л  O@ жр  U €€  Ж  b »  U €€  Н  b@ жш  U €€ getuid stat memset close setsockopt sqrt main _$_t7dVector1Zi checkpoint__Ct6dArray1ZiR9dofstream loadbalance__t7dVector1Zi lb_info__t6dArray1ZiRiT1Rd isisomorphic__t6dArray1ZiP7domeobj _$_t6dArray1Zi loadbalance__t6dArray1Zi _$_t7dVector1Zd checkpoint__Ct6dArray1ZdR9dofstream loadbalance__t7dVector1Zd lb_info__t6dArray1ZdRiT1Rd isisomorphic__t6dArray1ZdP7domeobj _$_t6dArray1Zd loadbalance__t6dArray1Zd _$_t6dArray1Zc checkpoint__Ct6dArray1ZcR9dofstream loadbalance__t6dArray1Zc lb_info__t6dArray1ZcRiT1Rd isisomorphic__t6dArray1ZcP7domeobj _$_t5stack1Zi _$_9difstream checkpoint__C13dfstream_baseR9dofstream loadbalance__13dfstream_base lb_info__13dfstream_baseRiT1Rd isisomorphic__13dfstream_baseP7domeobj _$_9dofstream _$_13dfstream_base _$_9difstream loadbalance__13dfstream_base lb_info__13dfstream_baseRiT1Rd isisomorphic__13dfstream_baseP7domeobj _$_9dofstream _$_9difstream loadbalance__13dfstream_base lb_info__13dfstream_baseRiT1Rd isisomorphic__13dfstream_baseP7domeobj _$_9dofstream _$_9difstream loadbalance__13dfstream_base lb_info__13dfstream_baseRiT1Rd isisomorphic__13dfstream_baseP7domeobj _$_9dofstream _$_9difstream loadbalance__13dfstream_base lb_info__13dfstream_baseRiT1Rd isisomorphic__13dfstream_baseP7domeobj _$_9dofstream _$_9difstream loadbalance__13dfstream_base lb_info__13dfstream_baseRiT1Rd isisomorphic__13dfstream_baseP7domeobj _$_9dofstream _$_7domeobj checkpoint__C7domeobjR9dofstream __pure_virtual _$_22_IO_ostream_withassign _$_22_IO_istream_withassign _$_8iostream _$_7istream _$_7ostream underflow__7filebuf overflow__7filebufi doallocate__7filebuf seekoff__7filebuflQ23ios8seek_diri seekpos__9streambufli _$_7filebuf sync__7filebuf pbackfail__9streambufi setbuf__7filebufPci xsputn__7filebufPCci xsgetn__7filebufPci get_column__9streambuf set_column__9streambufi sys_read__7filebufPci sys_seek__7filebuflQ23ios8seek_dir sys_write__7filebufPCci sys_stat__7filebufPv sys_close__7filebuf underflow__9streambuf overflow__9streambufi doallocate__9streambuf seekoff__9streambuflQ23ios8seek_diri _$_9streambuf sync__9streambuf setbuf__9streambufPci xsputn__9streambufPCci xsgetn__9streambufPci sys_read__9streambufPci sys_seek__9streambuflQ23ios8seek_dir sys_write__9streambufPCci sys_stat__9streambufPv sys_close__9streambuf _$_3ios underflow__10builtinbuf overflow__10builtinbufi doallocate__10builtinbuf seekoff__10builtinbuflQ23ios8seek_diri seekpos__10builtinbufli _$_10builtinbuf sync__10builtinbuf pbackfail__10builtinbufi setbuf__10builtinbufPci xsputn__10builtinbufPCci xsgetn__10builtinbufPci sys_read__10builtinbufPci sys_seek__10builtinbuflQ23ios8seek_dir sys_write__10builtinbufPCci sys_stat__10builtinbufPv sys_close__10builtinbuf _$_12ostdiostream _$_12istdiostream overflow__8stdiobufi _$_8stdiobuf sync__8stdiobuf xsputn__8stdiobufPCci sys_read__8stdiobufPci sys_seek__8stdiobuflQ23ios8seek_dir sys_write__8stdiobufPCci sys_close__8stdiobuf _GLOBAL_$I$dcout _GLOBAL_$I$__11_ios_fields _GLOBAL_$D$dcout _GLOBAL_$D$__11_ios_fields __main def_match enc_raw_init dec_raw_init enc_raw_any dec_raw_any enc_xdr_init dec_xdr_init enc_xdr_byte dec_xdr_byte enc_xdr_short dec_xdr_short enc_xdr_int dec_xdr_int enc_xdr_long dec_xdr_long dec_xdr_ushort dec_xdr_uint dec_xdr_ulong enc_xdr_float dec_xdr_float enc_xdr_double dec_xdr_double enc_xdr_cplx dec_xdr_cplx enc_xdr