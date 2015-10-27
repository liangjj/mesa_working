*deck evnodd
      subroutine evnodd(p,dp,ddp,n,m,dir)
c***begin prologue     evnodd
c***date written       940504   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            utility routine for producing even and odd
c***                   lagrange polynomials from unsymmetric polynomials.
c
c***description        note that the point distribution must be taken to
c***                   symmetrically about zero for proper results.
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       evnodd
c
      implicit integer (a-z)
      real*8 p, dp, ddp
      real*8 sym, asym
      character*(*) dir
      dimension p(n,n), dp(n,n), ddp(n,n)
      common /io/ inp, iout
      upper=n/2
      do 10 i=1,n
         bkwrd=n
         do 20 j=1,upper
            sym = p(i,j) + p(i,bkwrd)
            asym = p(i,j) - p(i,bkwrd)
            p(i,j) = sym
            p(i,bkwrd) = asym
            sym = dp(i,j) + dp(i,bkwrd)
            asym = dp(i,j) - dp(i,bkwrd)
            dp(i,j) = sym
            dp(i,bkwrd) = asym 
            sym = ddp(i,j) + ddp(i,bkwrd)
            asym = ddp(i,j) - ddp(i,bkwrd)
            ddp(i,j) = sym
            ddp(i,bkwrd) = asym
            bkwrd = bkwrd -1
 20      continue
 10   continue
c
c     the even polynomials are in the first vector slots and the odd
c     polynomials follow.
c
      check=mod(n,2)
      if(check.eq.0) then
         even = n/2
      else
         even = n/2 +1
      endif
      odd = n - even
      write(iout,1) even, odd
      if(dir.eq.'even') then
         m = even
      else
         m = odd
         count=0
         do 30 i=even+1,n
            count=count+1
            call copy(p(1,i),p(1,count),n)
            call copy(dp(1,i),dp(1,count),n)
            call copy(ddp(1,i),ddp(1,count),n)
 30      continue
      endif
 1    format(/,5x,'number of even polynomials = ',i3,/,5x,
     1            'number of odd polynomials  = ',i3)                   
      return
      end















