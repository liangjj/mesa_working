*deck chk.f
c***begin prologue     chk
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           orthogonal polynomials
c***author             schneider, barry (nsf)
c***source             math
c***purpose            check orthonormality of polynomials.
c***                   
c***references         
c
c***routines called    
c***end prologue       chk
      subroutine chk(pn,wts,scr,work,n)
      implicit integer (a-z)
      real*8 pn, wts, scr
      character*80 title
      dimension pn(n,n), wts(n), scr(n,n), work(n,n)
      common/io/inp, iout 
      call vsqrt(work,wts,n)
      call copy(pn,scr,n*n)
      call vmmul(work,scr,scr,n,n)
      call ebtc(work,scr,scr,n,n,n)
      title='overlap matrix'
      call prntrm(title,work,n,n,n,n,iout)
      call ebct(work,scr,scr,n,n,n)
      title='completeness matrix'
      call prntrm(title,work,n,n,n,n,iout)
      return
      end       
