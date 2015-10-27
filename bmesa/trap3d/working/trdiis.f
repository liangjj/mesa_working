*deck trdiis.f
c***begin prologue     trdiis
c***date written       961209   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           pulay, extrapolation
c***author             schneider, barry (nsf)
c***source             math
c***purpose            truncate the diis space.
c***                   
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       trdiis
      subroutine trdiis(vnl,psi,hpsi,b,hh,pp,ph,hp,old,new,n,
     1                  maxit,prnt)
      implicit integer (a-z)
      real*8 vnl, psi, hpsi, b, hh, pp, ph, hp
      logical prnt
      dimension vnl(n,*), psi(n,*), hpsi(n,*)
      dimension b(maxit+1,maxit+1), hh(maxit,maxit) 
      dimension pp(maxit,maxit), ph(maxit,maxit), hp(maxit,maxit)
      common/io/inp, iout
      drop=old-new
      drop=1
      write(iout,1) drop
      ci = 0
      do 10 i=drop+1,old
         ci = ci + 1
         cj = 0
         call copy(vnl(1,i),vnl(1,ci),n)
         call copy(psi(1,i),psi(1,ci),n)
         call copy(hpsi(1,i),hpsi(i,ci),n)
         do 20 j=drop+1,old
            cj = cj+1
            b(ci,cj)  = b(i,j)
            hh(ci,cj) = hh(i,j)
            hp(ci,cj) = hp(i,j)
            ph(ci,cj) = ph(i,j)
            pp(ci,cj) = pp(i,j)
 20      continue
 10   continue
      return
 1    format(/,5x,'truncating diis subspace by ',i3)
      end       
