*deck @(#)hbound.f
c***begin prologue     hbound
c***date written       920417   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           basis, link 6080
c***author             schneider, barry (lanl)
c***source             m6090
c***purpose            kohn matrix elements for numerical functions.
c***description        matrix elements of  h and s are computed for
c***                   bound components by radial quadrature. then a
c***                   transformation to an orthonormal set is performed
c***                   and the hamiltonian transformed to that basis.
c***                   a copy of the new hamiltonian is placed in the
c***                   input matrix.
c***references       
c
c***routines called
c***end prologue       hbound
      subroutine hbound(fns,ddfns,v,wt,hbbp,sbbp,hbbt,eig,work,temp,
     1                  scrb,nbfn,npts,count,prntbb)
      implicit integer (a-z)
      real *8 fns, ddfns, v, wt, hbbp, sbbp, hbbt, eig, work
      real *8 scrb, temp, tol, sqeig
      character *80 title
      logical prntbb
      dimension fns(npts,nbfn), ddfns(npts,nbfn)
      dimension v(npts), wt(npts), hbbp(nbfn,nbfn)
      dimension sbbp(nbfn,nbfn), hbbt(*), eig(nbfn), work(nbfn)
      dimension scrb(npts,nbfn), temp(nbf,nbf)
      common /io/ inp, iout
c
c               bound-bound matrix elements of  h 
c
      do 10 i=1,nbfn
         do 20 j=1,npts
            scrb(j,i)=(-.5d0*ddfns(j,i)+v(j)*fns(j,i))*wt(j)
   20    continue
   10 continue
      call ebtc(hbbp,fns,scrb,nbfn,npts,nbfn)
c
      if(prntbb) then
         title='bound-bound primitive hamiltonian'
         call prntrm(title,hbbp,nbfn,nbfn,nbfn,nbfn,iout)
      endif
c
c               bound-bound overlaps
c
      do 30 i=1,nbfn
         do 40 j=1,npts
            scrb(j,i)=fns(j,i)*wt(j)
   40    continue
   30 continue
      call ebtc(sbbp,fns,scrb,nbfn,npts,nbfn)
      if(prntbb) then
         title='bound-bound primitive overlap'
         call prntrm(title,sbbp,nbfn,nbfn,nbfn,nbfn,iout)
      endif
c               transform to orthogonal bound basis
c
      call tred2(nbfn,nbfn,sbbp,eig,work,sbbp)
      call tql2(nbfn,nbfn,eig,work,sbbp,ierr)
      count=0
      do 70 i=1,nbfn
         if (abs(eig(i)).gt.tol) then
             sqeig=1.d0/sqrt(eig(i))
             count=count+1
             do 80 j=1,nbfn
                sbbp(j,count)=sqeig*sbbp(j,i)
   80        continue
         endif
   70 continue     
      write(iout,90) count
c
c               bound-bound
c
      call ebc(temp,hbbp,sbbp,nbfn,nbfn,count)
      call ebtc(hbbt,sbbp,temp,count,nbfn,count)
      if(prntbb) then
         title='bound-bound final hamiltonian'
         call prntrm(title,hbbt,count,count,count,count,iout)
      endif
      call copy(hbbt,hbbp,count*count) 
      return
   90 format(//,5x,'number eigenvectors kept = ',1x,i4)  
      end






