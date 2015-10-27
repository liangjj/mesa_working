*deck ovmat.f
c***begin prologue     ovmat
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            lobatto overlap ovrix elements.
c***                   
c***references         
c
c***routines called    
c***end prologue       ovmat
      subroutine ovmat(ov,fa,fb,wt,n,prn,sec)
      implicit integer (a-z)
      real*8 ov, fa, fb, wt
      character*80 title
      character*3 itoc
      logical prn
      dimension ov(n,n), fa(n,n), fb(n,n), wt(n)
      common/io/inp, iout
      call rzero(ov,n*n)
      do 10 i=1,n
         ov(i,i) = ov(i,i) + fa(i,i)*wt(i)*fb(i,i)
 10   continue   
      if(prn) then
         title='unnormalized overlap matrix for sector = '//itoc(sec)
         call prntrm(title,ov,n,n,n,n,iout)
      endif
      return
      end       
