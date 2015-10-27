*deck rgamma.f
c***begin prologue     rgamma
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            R-matrix amplitudes.
c***                   
c***references         
c
c***routines called    
c***end prologue       rgamma
      subroutine rgamma(pham,sym,n,dim,prn)
      implicit integer (a-z)
      integer*8 pham, pgam
      character*80 title
      character*(*) sym
      logical prn
      dimension n(4,2)
      dimension pham(4,2)
      common/io/inp, iout
      if(dim.eq.1) then
         call gam1(pham(1,1),n(1,1),prn)
      elseif(dim.eq.2) then
c
c        the one particle functions and the two particle functions are
c        expanded in the DVR basis.  Therefore
c
c          Chi(1) = Sum d(k,p) phi(1) 
c             p      k            k
c  
c          Chi(2) = Sum d(k,p) psi(2) ; note that for the symmetric case 
c             p      k            k     psi = phi
c  
c        and
c
c          Psi(1,2) = Sum c(ji,q) [ phi(1) * psi(2) ]  unsymmetrized case
c             q        ij              i        j    
c
c          Psi(1,2) = Sum c(ij,q) [ phi(1) * phi(2) + phi(2) * phi(1) ] 
c             q        i>j             i        j        i        j    /sqrt(2)
c                                                      symmetrized case 
c                       
c
c        unsymmetrized case:
c   
c                   project onto Chi(1) to get
c
c         A(2) =  <Chi(1) | Psi(1,2) >  = Sum c(ji,q) * phi(2) * d(i,p)
c          pq         p        q      1    j,i             j
c
c                      project onto phi(2) to get
c
c         B(1) =  <Chi(2) | Psi(1,2) >  = Sum c(ji,q) * psi(1) * d(j,p)
c          pq         p        q      2    j,i             i
c
c
c                since the coefficients c(ij,q) are not
c                symmetric in (ij) the two sums are not
c                not necessarily the same, even if the same basis
c                is used to expand the two coordinate functions.  
c
c        in a DVR basis the sums simplify considerably, since we only need 
c        the value at the surface and at the surface only ONE, the last 
c        DVR function is non-zero.
c        thus:
c
         if(sym.eq.'unsymmetric') then
c
            call gkqunsy(pham,pgam,ngot,n,dim,prn)
c
c         A(a) =  <Chi(1) | Psi(1,a) >  = phi(a) * Sum c(ni,q) * d(i,p)
c          pq         p        q      1      n      i                 
c                           
c                           and similarly for B(a)
c
c         B(a) =  <Chi(2) | Psi(a,2) >  = psi(a) * Sum c(in,q) * d(i,p)
c          pq         p        q      1      n      i                 
c
c        symmetrized case:
c
         elseif(sym.eq.'symmetric') then
c
            call gkqsy(pham,pgam,ngot,n,dim,prn)
c
c
c        <Chi(1) |Psi(1,2)> = Sum c(ji,q) * 
c            p       q        i>j                                    
c                              k
c                                 [ delta(i,k) * d(k,p) * phi(2) 
c                                                         j
c                                                + 
c                                   delta(j,k) * d(k,p) * phi(2) ] / sqrt(2) 
c                                                            i
c
c        now set coordinate 2 = a which requires i = j = n
c        however, since i>j, only the second term survives.  
c
c                                        n-1
c      <Chi(1) | Psi(1,a)> = phi(a) * [ Sum c(jn,q) * d(j,p) * phi(a) 
c          p        q           n        j                           / sqrt(2) 
c                                                   +
c                                           c(nn,q) * d(n,p) ]
c
         else
            call lnkerr('error in symmetry')
         endif
      endif
      return      
      end       






