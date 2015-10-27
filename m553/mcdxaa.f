*deck @(#)mcdxaa.f	5.1  11/6/94
      subroutine mcdxaa(den,xjk,nij,mrs,nstep,itran,thrsh,hess,
     $     nstart,ntot)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcdxaa.f	5.1   11/6/94
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      implicit real*8(a-h,o-z)
cc
      dimension den(2),xjk(2),hess(2)
c
      common /io/ inp,iout
c
      common / number / zero,pt5,one,two,four,eight
c---------------------------------------------------------
c   this program multiplies a coulomb or exchange integral
c   block times a row or column of density matrix to build
c   the orbital hessian
c----------------------------------------------------------
c
c  den    density matrix elements
c  xjk    block of integrals
c  hess   the orbital hessian in the ao basis
c  nij    length of a column of the density matrix (dm)
c  mrs    the number integrals in xjk
c  thrsh  threshhold for neglecting an integral or dm element
c  itran  transformation flag which dictates the order of the loops
c---------
c         the next three parameters are used to read across the
c         dm rather than down a column
c---------
c  nstep  the size of the step used in moving thru the dm
c  ntot   the number of elements in this dm
c  nstart the starting row index in the dm
c
c     write(iout,9901) mrs
c9901 format('  entry  mcdxaa  mrs ',i8)
c
      mstep=iabs(nstep)
c
      if(nstep.ne.1)go to 100
c
      if(itran.ne.0)go to 40
      nx=0
      do 20 ij=1,nij
         denm=den(ij)
         if(abs(denm).lt.thrsh)go to 15
         do 10 lrs=1,mrs
            hess(nx+lrs)=hess(nx+lrs)+denm*xjk(lrs)
 10      continue
 15      continue
         nx=nx+mrs
 20   continue
c
      return
c
 40   continue
c
      do 60 lrs=1,mrs
         xx=xjk(lrs)
         if(abs(xx).lt.thrsh)go to 55
         nx=0
         do 50 ij=1,nij
            hess(nx+lrs)=hess(nx+lrs)+xx*den(ij)
            nx=nx+mrs
 50      continue
 55      continue
 60   continue
c
      return
c
 100  continue
c
      if(itran.ne.0)go to 140
      nx=0
      do 120 ij=nstart,ntot,mstep
         denm=den(ij)
c     write(iout,9906) denm
c9906 format('  cross mcdxaa  denm  ',f18.10)
         if(abs(denm).lt.thrsh)go to 115
         do 110 lrs=1,mrs
            hess(nx+lrs)=hess(nx+lrs)+denm*xjk(lrs)
 110     continue
 115     continue
         nx=nx+mrs
 120  continue
c
      return
c
 140  continue
c
      do 160 lrs=1,mrs
         xx=xjk(lrs)
         if(abs(xx).lt.thrsh)go to 155
         nx=0
         do 150 ij=nstart,ntot,mstep
            hess(nx+lrs)=hess(nx+lrs)+xx*den(ij)
            nx=nx+mrs
 150     continue
 155     continue
 160  continue
      return
      end
