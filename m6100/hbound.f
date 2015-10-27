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
     1                  scrb,nbf,npts,count,prntbb)
      implicit integer (a-z)
      real *8 fns, ddfns, v, wt, hbbp, sbbp, hbbt, eig, work
      real *8 scrb, temp, tol, sqeig
      character *80 title
      logical prntbb
      dimension fns(npts,nbf), ddfns(npts,nbf)
      dimension v(npts), wt(npts), hbbp(nbf,nbf)
      dimension sbbp(nbf,nbf), hbbt(*), eig(nbf), work(nbf)
      dimension scrb(npts,nbf), temp(nbf,nbf)
      common /io/ inp, iout
      data tol/1.e-06/
c
c               bound-bound matrix elements of  h 
c
      do 10 i=1,nbf
         do 20 j=1,npts
            scrb(j,i)=(-.5d0*ddfns(j,i)+v(j)*fns(j,i))*wt(j)
   20    continue
   10 continue
      call ebtc(hbbp,fns,scrb,nbf,npts,nbf)
c
      if(prntbb) then
         title='bound-bound primitive hamiltonian'
         call prntrm(title,hbbp,nbf,nbf,nbf,nbf,iout)
      endif
c
c               bound-bound overlaps
c
      do 30 i=1,nbf
         do 40 j=1,npts
            scrb(j,i)=fns(j,i)*wt(j)
   40    continue
   30 continue
      call ebtc(sbbp,fns,scrb,nbf,npts,nbf)
      if(prntbb) then
         title='bound-bound primitive overlap'
         call prntrm(title,sbbp,nbf,nbf,nbf,nbf,iout)
         call copy(sbbp,hbbt,nbf*nbf)
      endif
c               transform to orthogonal bound basis
c
      call tred2(nbf,nbf,sbbp,eig,work,sbbp)
      call tql2(nbf,nbf,eig,work,sbbp,ierr)
      write(iout,100) (eig(i),i=1,nbf)
      count=0
      do 70 i=1,nbf
         if (abs(eig(i)).gt.tol) then
             sqeig=1.d0/sqrt(eig(i))
             count=count+1
             do 80 j=1,nbf
                sbbp(j,count)=sqeig*sbbp(j,i)
   80        continue
         endif
   70 continue     
      write(iout,90) count
      if (prntbb) then
          call ebc(temp,hbbt,sbbp,nbf,nbf,count)
          call ebtc(hbbt,sbbp,temp,count,nbf,count)
          title='transformation matrix'
          call prntrm(title,sbbp,nbf,count,nbf,nbf,iout)
          title='overlap test'
          call prntrm(title,hbbt,count,count,count,count,iout)
      endif
c
c               bound-bound
c
      call ebc(temp,hbbp,sbbp,nbf,nbf,count)
      call ebtc(hbbt,sbbp,temp,count,nbf,count)
      if(prntbb) then
         title='bound-bound final hamiltonian'
         call prntrm(title,hbbt,count,count,count,count,iout)
      endif
      call rzero(hbbp,nbf*nbf)
      call copy(hbbt,hbbp,count*count) 
      return
   90 format(//,5x,'number eigenvectors kept = ',1x,i4,/)  
  100 format(/,5x,'eigenvalues',/,(/,5(1x,e15.8))) 
      end






