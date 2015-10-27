*deck toorth
c***begin prologue     toorth
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           schmidt, orthogonal
c***author             schneider, barry (nsf)
c***source             
c***purpose            orthonormal open channel orbitals
c***description        orthogonalize open to closed channel functions.
c***                   form open channel orthonormal set by diagonalizing
c***                   open channel overlap matrix.  finally transform 
c***                   functions to new open channel basis.
c***references       
c
c***routines called
c***end prologue       toorth
      subroutine toorth(fi,dfi,ddfi,fj,dfj,ddfj,flsti,flstj,dflsti,
     1                  dflstj,s,eig,dum,scr,ni,nj,npts,n,cnt2,prnt)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 fi, dfi, ddfi, fj, dfj, ddfj
      real*8 flsti, dflsti, flstj, dflstj, s, eig, dum, tmp, scr
      real*8 tol
      logical prnt
      character*80 title
      dimension fi(npts,ni), dfi(npts,ni), ddfi(npts,ni)
      dimension fj(npts,nj), dfj(npts,nj), ddfj(npts,nj)
      dimension flsti(ni), flstj(nj), dflsti(ni), dflstj(nj)
      dimension eig(nj), dum(nj), scr(*)
      dimension s(n,n)
      parameter (tol=1.d-07)
      nbeg=ni+1
      nfinal=ni+nj
c     form new open channel function orthogonal to closed channel space
c     using schmidt procedure.
c     then do same for derivatives and second derivatives.      
      call ambcxx(fj,fi,s(1,nbeg),npts,ni,nj,npts,npts,n)
      call ambcxx(dfj,dfi,s(1,nbeg),npts,ni,nj,npts,npts,n)
      call ambcxx(ddfj,ddfi,s(1,nbeg),npts,ni,nj,npts,npts,n)
c     do the same for the functions and derivatives at last point.
      call ambtcxx(flstj,s(1,nbeg),flsti,nj,ni,1,nj,n,ni)
      call ambtcxx(dflstj,s(1,nbeg),dflsti,nj,ni,1,nj,n,ni)
c     calculate new overlap matrix
      call ebtcxx(s,fi,fi,ni,npts,ni,n,npts,npts)
      call ebtcxx(s(1,nbeg),fi,fj,ni,npts,nj,n,npts,npts)
      call ebtcxx(s(nbeg,1),fj,fi,nj,npts,ni,n,npts,npts)
      call ebtcxx(s(nbeg,nbeg),fj,fj,nj,npts,nj,n,npts,npts)
      if (prnt) then
          title='final overlap matrix before open channel '//
     1          'diagonalization'
          call prntrm(title,s,n,n,n,n,iout)
      endif          
c     diagonalize open channel overlap matrix
      call tred2(n,nj,s(nbeg,nbeg),eig,dum,s(nbeg,nbeg))
      call tql2(n,nj,eig,dum,s(nbeg,nbeg),ierr)
      if (prnt) then
          write(iout,1) (eig(i),i=1,nj)
      endif          
      cnt1=0
      cnt2=ni 
      do 10 i=nbeg,nfinal
         cnt1=cnt1+1
         if (abs(eig(cnt1)).ge.tol) then
             cnt2=cnt2+1
             tmp=1.d0/sqrt(eig(cnt1))
             do 20 j=nbeg,nfinal
                s(j,cnt2)=tmp*s(j,i)
   20        continue
         endif
   10 continue
      cnt2=cnt2-ni
c     form new open orthonormal basis set and same for
c     derivatives and second derivatives
      call ebcxx(scr,fj,s(nbeg,nbeg),npts,nj,cnt2,npts,npts,n)
      call copy(scr,fj,npts*cnt2)
      call ebtcxx(scr,s(nbeg,nbeg),flstj,cnt2,nj,1,cnt2,n,nj)
      call copy(scr,flstj,cnt2) 
      call ebcxx(scr,dfj,s(nbeg,nbeg),npts,nj,cnt2,npts,npts,n)
      call copy(scr,dfj,npts*cnt2)
      call ebtcxx(scr,s(nbeg,nbeg),dflstj,cnt2,nj,1,cnt2,n,nj)
      call copy(scr,dflstj,cnt2)       
      call ebcxx(scr,ddfj,s(nbeg,nbeg),npts,nj,cnt2,npts,npts,n)
      call copy(scr,ddfj,npts*cnt2)
c     we are now in a completely orthonormal basis
    1 format(/,5x,'overlap eigenvalues of open channel subspace',
     1         (/,5x,5e15.8))                                    
      return
      end



