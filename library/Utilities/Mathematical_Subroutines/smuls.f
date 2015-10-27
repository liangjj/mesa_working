*deck @(#)smuls.f	1.1  11/30/90
      subroutine smuls(v,w,sa,sb,n)
c***begin prologue     smul
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, scalar, multiply.
c***author             saxe, paul (lanl)
c***source             @(#)smul.f	1.1   11/30/90
c***purpose            vectorized scalar multiply v=w*sa + sb, where v,w are 
c                      vectors and sa, sb are scalars.
c***description
c                      call smul(v,w,s,n)
c                        v        output vector dimensioned (n).
c                        w        input vector dimensioned (n).
c                        sa       scalar.
c                        sb       scalar
c
c***references
c***routines called    (none)
c***end prologue       smul
      real*8 v(n),w(n),sa, sb
      do 1 i=1,n
         v(i)=w(i)*sa + sb
    1 continue
      return
      end
