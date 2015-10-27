*deck dpcoef
      subroutine dpcoef (l, c, tc, a)
c***begin prologue  dpcoef
c***purpose  convert the dpolft coefficients to taylor series form.
c***library   slatec
c***category  k1a1a2
c***type      double precision (pcoef-s, dpcoef-d)
c***keywords  curve fitting, data fitting, least squares, polynomial fit
c***author  shampine, l. f., (snla)
c           davenport, s. m., (snla)
c***description
c
c     abstract
c
c     dpolft  computes the least squares polynomial fit of degree  l  as
c     a sum of orthogonal polynomials.  dpcoef  changes this fit to its
c     taylor expansion about any point  c , i.e. writes the polynomial
c     as a sum of powers of (x-c).  taking  c=0.  gives the polynomial
c     in powers of x, but a suitable non-zero  c  often leads to
c     polynomials which are better scaled and more accurately evaluated.
c
c     the parameters for  dpcoef  are
c
c     input -- all type real variables are double precision
c         l -      indicates the degree of polynomial to be changed to
c                  its taylor expansion.  to obtain the taylor
c                  coefficients in reverse order, input  l  as the
c                  negative of the degree desired.  the absolute value
c                  of l  must be less than or equal to ndeg, the highest
c                  degree polynomial fitted by  dpolft .
c         c -      the point about which the taylor expansion is to be
c                  made.
c         a -      work and output array containing values from last
c                  call to  dpolft .
c
c     output -- all type real variables are double precision
c         tc -     vector containing the first ll+1 taylor coefficients
c                  where ll=abs(l).  if  l.gt.0 , the coefficients are
c                  in the usual taylor series order, i.e.
c                    p(x) = tc(1) + tc(2)*(x-c) + ... + tc(n+1)*(x-c)**n
c                  if l .lt. 0, the coefficients are in reverse order,
c                  i.e.
c                    p(x) = tc(1)*(x-c)**n + ... + tc(n)*(x-c) + tc(n+1)
c
c***references  l. f. shampine, s. m. davenport and r. e. huddleston,
c                 curve fitting by polynomials in one variable, report
c                 sla-74-0270, sandia laboratories, june 1974.
c***routines called  dp1vlu
c***revision history  (yymmdd)
c   740601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891006  cosmetic changes to prologue.  (wrb)
c   891006  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dpcoef
c
      integer i,l,ll,llp1,llp2,new,nr
      double precision a(*),c,fac,save,tc(*)
c***first executable statement  dpcoef
      ll = abs(l)
      llp1 = ll + 1
      call dp1vlu (ll,ll,c,tc(1),tc(2),a)
      if (ll .lt. 2) go to 2
      fac = 1.0d0
      do 1 i = 3,llp1
        fac = fac*(i-1)
 1      tc(i) = tc(i)/fac
 2    if (l .ge. 0) go to 4
      nr = llp1/2
      llp2 = ll + 2
      do 3 i = 1,nr
        save = tc(i)
        new = llp2 - i
        tc(i) = tc(new)
 3      tc(new) = save
 4    return
      end
