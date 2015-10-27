*deck cpzero
      subroutine cpzero (in, a, r, t, iflg, s)
c***begin prologue  cpzero
c***purpose  find the zeros of a polynomial with complex coefficients.
c***library   slatec
c***category  f1a1b
c***type      complex (rpzero-s, cpzero-c)
c***keywords  polynomial roots, polynomial zeros, real roots
c***author  kahaner, d. k., (nbs)
c***description
c
c      find the zeros of the complex polynomial
c         p(z)= a(1)*z**n + a(2)*z**(n-1) +...+ a(n+1)
c
c    input...
c       in = degree of p(z)
c       a = complex vector containing coefficients of p(z),
c            a(i) = coefficient of z**(n+1-i)
c       r = n word complex vector containing initial estimates for zeros
c            if these are known.
c       t = 4(n+1) word array used for temporary storage
c       iflg = flag to indicate if initial estimates of
c              zeros are input.
c            if iflg .eq. 0, no estimates are input.
c            if iflg .ne. 0, the vector r contains estimates of
c               the zeros
c       ** warning ****** if estimates are input, they must
c                         be separated, that is, distinct or
c                         not repeated.
c       s = an n word array
c
c    output...
c       r(i) = ith zero,
c       s(i) = bound for r(i) .
c       iflg = error diagnostic
c    error diagnostics...
c       if iflg .eq. 0 on return, all is well
c       if iflg .eq. 1 on return, a(1)=0.0 or n=0 on input
c       if iflg .eq. 2 on return, the program failed to converge
c                after 25*n iterations.  best current estimates of the
c                zeros are in r(i).  error bounds are not calculated.
c
c***references  (none)
c***routines called  cpevl
c***revision history  (yymmdd)
c   810223  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  cpzero
c
      real  s(*)
      complex r(*),t(*),a(*),pn,temp
c***first executable statement  cpzero
      if( in .le. 0 .or. abs(a(1)) .eq. 0.0 ) go to 30
c
c       check for easily obtained zeros
c
      n=in
      n1=n+1
      if(iflg .ne. 0) go to 14
    1 n1=n+1
      if(n .gt. 1) go to 2
         r(1)=-a(2)/a(1)
         s(1)=0.0
         return
    2 if( abs(a(n1)) .ne. 0.0 ) go to 3
         r(n)=0.0
         s(n)=0.0
         n=n-1
         go to 1
c
c          if initial estimates for zeros not given, find some
c
    3 temp=-a(2)/(a(1)*n)
      call cpevl(n,n,a,temp,t,t,.false.)
      imax=n+2
      t(n1)=abs(t(n1))
      do 6 i=2,n1
         t(n+i)=-abs(t(n+2-i))
         if(real(t(n+i)) .lt. real(t(imax))) imax=n+i
    6 continue
      x=(-real(t(imax))/real(t(n1)))**(1./(imax-n1))
    7 x=2.*x
         call cpevl(n,0,t(n1),cmplx(x,0.0),pn,pn,.false.)
      if (real(pn).lt.0.) go to 7
      u=.5*x
      v=x
   10 x=.5*(u+v)
         call cpevl(n,0,t(n1),cmplx(x,0.0),pn,pn,.false.)
         if (real(pn).gt.0.) v=x
         if (real(pn).le.0.) u=x
         if((v-u) .gt. .001*(1.+v)) go to 10
      do 13 i=1,n
         u=(3.14159265/n)*(2*i-1.5)
   13    r(i)=max(x,.001*abs(temp))*cmplx(cos(u),sin(u))+temp
c
c          main iteration loop starts here
c
   14 nr=0
      nmax=25*n
      do 19 nit=1,nmax
         do 18 i=1,n
            if(nit .ne. 1 .and. abs(t(i)) .eq. 0.) go to 18
               call cpevl(n,0,a,r(i),pn,temp,.true.)
               if(abs(real(pn))+abs(aimag(pn)) .gt. real(temp)+
     1              aimag(temp)) go to 16
                  t(i)=0.0
                  nr=nr+1
                  go to 18
   16          temp=a(1)
               do 17 j=1,n
   17             if(j .ne. i) temp=temp*(r(i)-r(j))
               t(i)=pn/temp
   18    continue
         do 15 i=1,n
   15       r(i)=r(i)-t(i)
         if(nr .eq. n) go to 21
   19 continue
      go to 26
c
c          calculate error bounds for zeros
c
   21 do 25 nr=1,n
         call cpevl(n,n,a,r(nr),t,t(n+2),.true.)
         x=abs(cmplx(abs(real(t(1))),abs(aimag(t(1))))+t(n+2))
         s(nr)=0.0
         do 23 i=1,n
            x=x*real(n1-i)/i
            temp=cmplx(max(abs(real(t(i+1)))-real(t(n1+i)),0.0),
     1           max(abs(aimag(t(i+1)))-aimag(t(n1+i)),0.0))
   23       s(nr)=max(s(nr),(abs(temp)/x)**(1./i))
   25    s(nr)=1./s(nr)
      return
c        error exits
   26 iflg=2
      return
   30 iflg=1
      return
      end
