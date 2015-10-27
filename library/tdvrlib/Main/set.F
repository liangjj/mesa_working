*deck set.f 
c***begin prologue     set
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           time, dvr, orthogonal polynomial
c***                   bose-einstein
c***author             schneider, b. i.(nsf)
c***source             tprop
c***purpose            set pointers to first word of the dvr arrays 
c***                   needed in the tdvr code.
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       set
      subroutine set(pham,pv,h0,eigv0,eig0,q0,q,qwt,pq,dpq,n,
     1               prnt,notim,spac)
c
      implicit integer (a-z)
#ifdef DECPOINTER
      integer*8 pham, pv
#endif DECPOINTER
#ifdef SGIPOINTER
      integer*4 pham, pv
#endif SGIPOINTER
      real*8 ham, v
      logical prnt, notim, spac
      character*80 title
      common /io/ inp, iout
      pointer (pham,ham(1))
      pointer (pv,v(1))
c
c     one body hamiltonian, potential, transformation matrix, points, weights,
c     and basis function location, follow for spatial variables.
c
      if(notim) then
         h0=1
         v0=h0+n*n
         hx=v0+n
         srf=hx+n*n
         q0=srf+2
         q=q0+1
         qwt=q+n
         pq=qwt+n
         dpq=pq+n*n
         ddpq=dpq+n*n
         eigv0=ddpq+n*n
         rgama=eigv0+n*n
         eig0=rgama+2*n
         eigvx=eig0+n
         rgamax=eigvx+n*n
         eigx=rgamax+2*n
         if(prnt) then
            title='hamiltonian'
            call prntrm(title,ham(h0),n,n,n,n,iout)
            title='potential'
            call prntrm(title,ham(v0),n,1,n,1,iout)
            title='eigenvectors'
            call prntrm(title,ham(eigv0),n,n,n,n,iout)
            title='eigenvalues'
            call prntrm(title,ham(eig0),n,1,n,1,iout)
            write(iout,1) ham(q0)
            title='coordinates'
            call prntrm(title,ham(q),n,1,n,1,iout)
            title='weights'
            call prntrm(title,ham(qwt),n,1,n,1,iout)
            title='polynomials'
            call prntrm(title,ham(pq),n,n,n,n,iout)
            title='first derivative of polynomials'
            call prntrm(title,ham(dpq),n,n,n,n,iout)
         endif
         if(spac) then
            call rzero(ham(h0),n*n)
         endif
      else
         h0=1
         q0=h0+n*n
         q=q0+1
         qwt=q+n
         pq=qwt+n
         dpq=pq+n*n
         if(prnt) then
            title='hamiltonian'
            call prntrm(title,ham(h0),n,n,n,n,iout)
            title='potential'
            call prntrm(title,v,n,1,n,1,iout)
            write(iout,1) ham(q0)
            title='coordinates'
            call prntrm(title,ham(q),n,1,n,1,iout)
            title='weights'
            call prntrm(title,ham(qwt),n,1,n,1,iout)
            title='polynomials'
            call prntrm(title,ham(pq),n,n,n,n,iout)
            title='first derivative of polynomials'
            call prntrm(title,ham(dpq),n,n,n,n,iout)
         endif
      endif
      return
 1    format(/,5x,'first coordinate value = ',e15.8)
      end

