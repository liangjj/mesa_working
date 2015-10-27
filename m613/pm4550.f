*deck @(#)pm4550.f	5.1 11/6/94
      subroutine pm4550(z,a)
c***begin prologue     pm4550.f
c***date written       920324   (yymmdd)
c***revision date      @(#)pm4550.f	5.1
c
c
c***keywords
c***author             
c***source             @(#)pm4550.f	5.1 11/6/94 
c***purpose
c                       
c***description
c***references
c
c***routines called
c
c***end prologue       pm4550.f
c
      implicit integer (a-z)
c
      real*8 z, alpha(100), r(100), fn(100)
      character*1600 card
      character*8 cpass
      character*12 fptoc
      character*80 title
      dimension z(*), a(*)
c     
      common/io/inp,iout
c
      call posinp('$input',cpass)
      call cardin(card)
      lmax=intkey(card,'maximum-l-value',1,' ')
      nalpha=intkey(card,'number-of-alpha-values',1,' ')
      call fparr(card,'alpha-values',alpha,nalpha,' ')
      nr=intkey(card,'number-of-r-values',1,' ')
      call fparr(card,'r-values',r,nr,' ')
      do 10 i=1,nalpha
         title='f(l,r) for alpha = '//fptoc(alpha(i))
         do 20 j=1,nr
            call gfunct(fn,alpha(i),r(j),lmax,.true.)
   20    continue
   10 continue             
c
      call chainx(0)
c
c
      stop
      end
