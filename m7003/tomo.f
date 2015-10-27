*deck @(#)onebb.f	1.1 9/8/91
c***begin prologue     m7001
c***date written       920525   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m7001, link 7001, spline
c***author             schneider, b. (nsf)
c***source             m7001
c***purpose            bound-bound one electron integrals
c***
c***description        compute all standard one electron integrals
c***                   between real orbitals by numerical quadrature.
c***                   the routine is vectorized for efficiency.
c***                 
c***                 
c***references       
c
c***routines called    ebtc(mylib),vmul(math),iosys(io),sscal(clams)
c***end prologue       m7001
      subroutine onebb(fns,delfns,scr,onemat,mask,pt,z,n,nbf,prnt)
      implicit integer (a-z)
      real *8 fns, delfns, scr, pt, z, onemat, mask
      character *80 title
      logical prnt
      dimension fns(n,nbf), delfns(n,nbf), pt(n), onemat(nbf,nbf)
      dimension mask(nbf,nbf), scr(n,nbf)
      common /io/ inp, iout
c**********************************************************************c
c                       overlap                                        c
c**********************************************************************c
      call ebtc(onemat,fns,fns,nbf,n,nbf)
      call vmul(onemat,onemat,mask,nbf*nbf)
      call iosys ('write real "real bb atomic overlap integrals" '//
     1            'to atomci',nbf*nbf,onemat,0,' ')
      if (prnt) then
          title='real-real overlap integrals'
          call prntrm(title,onemat,nbf,nbf,nbf,nbf,iout)
      endif
c**********************************************************************c
c                   orthonormalize the bound orbitals                  c
c**********************************************************************c
      call copy(onemat,smat,nbf*nbf)
      call tred2(nbf,nbf,smat,eig,work,onemat)
      call tql2(nbf,nbf,eig,work,smat,ierr) 
      if (prnt) then
          write(iout,30) (eig(i),i=1,nbf)
      endif
      nbfo=0
      do 100 i=1,nbf      
         if (abs(eig(i)).gt.tol) then
             nbfo=nbfo+1
             eig(i)=1.d0/sqrt(eig(i))
             do 110 j=1,nbf
                smat(j,nbfo)=eig(i)*smat(j,i)
  110        continue
         endif
  100 continue    
      if (prnt) then
          write(iout,40) tol, nbfo
      endif
      call iosys('write integer "number of orthogonal bound '//
     1           'functions" to atomci',1,nbfo,0,' ')
      call iosys('write real "bound transformation matrix" to atomci',
     1            nbf*nbfo,smat,0,' ')
c**********************************************************************c
c                       kinetic energy                                 c
c**********************************************************************c
      call ebtc(onemat,fns,delfns,nbf,n,nbf)
      call vmul(onemat,onemat,mask,nbf*nbf)
      call iosys ('write real "real bb atomic kinetic energy '//
     1            'integrals" to atomci',nbf*nbf,onemat,0,' ')
      if (prnt) then
          title='real-real kinetic energy integrals'
          call prntrm(title,onemat,nbf,nbf,nbf,nbf,iout)
      endif
c**********************************************************************c
c                       potential energy                               c
c**********************************************************************c
      do 10 i=1,nbf
         do 20 j=1,n
            scr(j,i)=fns(j,i)/pt(j)
   20    continue
   10 continue     
      call ebtc(onemat,fns,scr,nbf,n,nbf)
      call vmul(onemat,onemat,mask,nbf*nbf)
      call sscal(nbf*nbf,-z,onemat,1)
      call iosys ('write real "real bb atomic potential energy '//
     1            'integrals" to atomci',nbf*nbf,onemat,0,' ')
      if (prnt) then
          title='real-real potential energy integrals'
          call prntrm(title,onemat,nbf,nbf,nbf,nbf,iout)
      endif
      return
   30 format(//,5x,'eigenvalues of bound overlap matrix',
     1              (/,5x,5(e15.8,1x)))
   40 format(//,5x,'number of eigenvectors above tolerance of',1x,
     1              e15.8,1x,'is',1x,i3)  
      end







