*deck dslv.f
c***begin prologue     dslv
c***date written       970503   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           precondition
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***routines called    
c***end prologue       dslv
      subroutine dslv(dinv,ham,vecold,vecnew,ndvr,nvc,prnt)
      implicit integer (a-z)
      character*80 title
      logical prnt
      complex*16 dinv, ham, vecold, vecnew
      dimension dinv(ndvr), ham(ndvr,ndvr), vecold(ndvr,nvc)
      dimension vecnew(ndvr,nvc) 
      common/io/inp, iout
      call cebc(vecnew,ham,vecold,ndvr,ndvr,nvc)
      do 10 i=1,ndvr
         do 20 j=1,nvc
            vecnew(i,j) = dinv(i)*vecnew(i,j)
   20    continue
   10 continue
      if(prnt) then
         title='new trial vectors'
         call prntcm(title,vecnew,ndvr,nvc,ndvr,nvc,iout)
      endif
      return
      end       
