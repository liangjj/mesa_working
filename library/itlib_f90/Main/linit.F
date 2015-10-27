*deck linit.f
c***begin prologue     linit
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source                           
c***purpose            initialize davidson
c***                    
c***description        
c                                                                  
c***references         
c***routines called    
c***end prologue       linit
      subroutine linit(ht,hx,hy,hz,v,driver,vec,hvec,h,htmp,b,btmp,
     1                 n,m,nt,nx,ny,nz,nbeg,nend,nout,spdim,maxvec)
      implicit integer (a-z)
      real*8 ht, hx, hy, hz, v, driver, vec, hvec, h, htmp
      real*8 b, btmp
      dimension ht(nt,nt), hx(nx,nx), hy(ny,ny), hz(nz,nz)
      dimension v(n,nt), driver(m), vec(m,maxvec), hvec(m,maxvec)
      dimension h(maxvec,maxvec), htmp(maxvec,maxvec), b(maxvec)
      dimension btmp(maxvec)
      common/io/inp, iout
c
c
c    initialize the effect of the hamiltonian on these vectors.
c
      call hvdvr(ht,hx,hy,hz,v,vec(1,nbeg),hvec(1,nbeg),n,nt,nx,ny,nz,
     1           nout,spdim)
c
c     initialize the small hamiltonian matrix and right hand side.
c
      call hinit(h,htmp,b,btmp,vec,hvec,driver,m,1,nend,maxvec)
      return
      end       



