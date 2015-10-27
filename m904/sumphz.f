*deck @(#)sumphz.f	5.1  11/6/94
      subroutine sumphz(ms,ifiga,ifigb,ifig,isuma,isumb,nphase)
c
c  calculate the sums of the orbital numbers for the alpha and beta
c  lists, and the phase of the determinants ((-1)**n stored as mod(n,2))
c
      common /c1/ nbf,nel,ntotfg,mxopn,maxdet,ndbgu,mxocc,maxne0,mel
      common /c2/ msbas,jsym,nopen,kdbl,ndeti,nsefi,ncpi
c
      dimension ms(*), ifiga(mxocc,*), ifigb(mxocc,*), ifig(nbf,2),
     1  isuma(*), isumb(*), nphase(*)
c
c  construct the determinants for configuration i.
c
      do 80 ii=1,ndeti
c
c  calculate nalpha and nbeta for this determinant
c
      nbeta = mxocc - ms(ii)
      nalpha  = nel - nbeta
c
      isign = 0
c
c  clear the array ifig
c
      do 10 kk=1,2
      do 10 ia=1,nbf
   10 ifig(ia,kk) = 0
c
c  put a 1 in ifig(i,1) if the ith mo is occupied with alpha spin
c
      lsum = 0
      do 30 ka=1,nalpha
      ia = ifiga(ka,ii)
      lsum = lsum + ia
      ifig(ia,1) = 1
      if(ia.eq.nbf) go to 30
      ia = ia + 1
c
c  count the number of orbitals between the present orbital and
c  the end of the list to determine the sign of the determinant
c
      do 20 kk=ia,nbf
   20 isign = isign + ifig(kk,1)
c
   30 continue
c
      isuma(ii) = lsum
c
c  process the list of beta spin orbitals in a similar manner
c
      lsum = 0
      do 50 kb=1,nbeta
      ia = ifigb(kb,ii)
      lsum = lsum + ia
      ifig(ia,2) = 1
      if(ia.eq.nbf) go to 50
      ia = ia + 1
c
      do 40 kk=ia,nbf
   40 isign = isign + ifig(kk,2)
c
   50 continue
c
      isumb(ii) = lsum
c
      nphase(ii) = mod(isign,2)
c
c  rearrange the list of alpha and beta orbitals according to
c  increasing orbital number
c
      ka = 0
      kb = 0
      do 70 ia=1,nbf
      if( ifig(ia,1).ne.0 ) then
        ka = ka + 1
        ifiga(ka,ii) = ia
      endif
      if( ifig(ia,2).ne.0 ) then
        kb = kb + 1
        ifigb(kb,ii) = ia
      endif
   70 continue
c
   80 continue
c
      return
      end
