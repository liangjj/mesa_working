      subroutine genvec (gr1,gr2,varr,veco,vecn,t1,wts,msize)
      implicit integer(a-z)
      real *8 gr1, gr2, varr, veco, vecn, t1, wts, sumf, sumb
      dimension varr(msize), veco(msize), vecn(msize), wts(msize)
      dimension gr1(msize), gr2(msize), t1(msize)
*
*          subroutine to generate a new iteration in the vector
*          sequence to solve the integral equations.
*          the method used is a direct technique in which the
*          matrix is not actually formed. its effect on the old vector
*          is computed directly from the regular and
*          irregular solutions , the potential matrix and the
*          previous vector in the chain.
*
*
*
      call rzero(vecn,msize)
      call rzero(t1,msize)
*
      do 20 i=1,msize
         t1(i)=-varr(i)*veco(i)
   20 continue
*
*          now begin main loops.
*
*         do forward recursion
*
      sumf=0.d+00
      do 30 i=1,msize
      sumf=sumf+gr1(i)*wts(i)*t1(i)
      vecn(i)=vecn(i)+gr2(i)*sumf
   30 continue
*
*         do backward recursion
*
      ptmnus=msize-1
      sumb=0.d+00
      do 60 i=ptmnus,1,-1
      i1=i+1
      sumb=sumb+gr2(i1)*wts(i1)*t1(i1)
      vecn(i)=vecn(i)+gr1(i)*sumb
   60 continue
      return
      end
