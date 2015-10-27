*deck preslv.f
c***begin prologue     preslv
c***date written       970503   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           precondition
c***author             schneider, barry (nsf)
c***source             
c***purpose            create a preconditioned trial vector based
c***                   on a soluble zeroth order hamiltonian.
c***references 
c                      dinv(n) = inverse of ( energy - eigenvalue )        
c                      v(n,n) = diagonal matrix of potential in the dvr
c                               representation.
c                      vecold(n,nvc) = vector from previous iteration in
c                                      dvr representation
c                      vecnew(n,nvc) = new vector in dvr rerpresentation
c                      ul(n,n) = left transformation matrix from dvr to h0
c                                representation.
c                      ur(n,n) = right transformation matrix from dvr to h0
c                                representation.
c***routines called    
c***end prologue       preslv
      subroutine preslv(dinv,v,vecold,vecnew,ul,ur,t1,t2,
     1                  ndvr,nh0,nvc,prnt)
      implicit integer (a-z)
      character*80 title
      logical prnt
      complex*16 dinv, v, vecold, vecnew, ul, ur
      complex*16 t1, t2
      dimension dinv(ndvr), v(ndvr,ndvr), vecold(ndvr,nvc)
      dimension vecnew(ndvr,nvc), ul(ndvr,nh0), ur(ndvr,nh0)
      dimension t1(ndvr,nvc), t2(nh0,nvc) 
      common/io/inp, iout
      do 10 i=1,ndvr
         do 20 j=1,nvc
            t1(i,j) = v(i,i)*vecold(i,j)
   20    continue
   10 continue
      call cebtc(t2,ul,t1,nh0,ndvr,nvc)
      do 30 i=1,nh0
         do 40 j=1,nvc
            t2(i,j) = dinv(i)*t2(i,j)
 40      continue
 30   continue   
      call cebc(vecnew,ur,t2,ndvr,nh0,nvc)
      if(prnt) then
         title='new trial vectors'
         call prntrm(title,vecnew,ndvr,nvc,ndvr,nvc,iout)
      endif
      return
      end       
