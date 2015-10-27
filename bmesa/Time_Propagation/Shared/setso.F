c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Ch3d}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck setso.f 
***begin prologue     setso
***date written       960718   (yymmdd)
***revision date               (yymmdd)
***keywords           time, dvr, fd
***author             schneider, b. i.(nsf)
***source             
***purpose            set pointers to first word of the arrays 
***                   needed in the time propagation code.
***routines called    iosys, util and mdutil
***end prologue       setso
      subroutine setso(pham,h0,q0,q,qwt,pq,dpq,n,ndiff,prnt,typke,spac)
      implicit integer (a-z)
#ifdef DECPOINTER
      integer*8 pham
#endif DECPOINTER
#ifdef SGIPOINTER
      integer*4 pham
#endif SGIPOINTER
      real*8 ham
      logical prnt, spac
      character*(*) typke
      character*80 title
      common /io/ inp, iout
      pointer (pham,ham(1))

c     One body hamiltonian, potential, transformation matrix, points, weights,
c     and basis function location, follow for spatial variables.

      if(typke.eq.'dvr'.or.typke.eq.'packed') then
         h=1
         vphy=h+n*n
         h0=vphy+n
         srf=h0+n*n
         q0=srf+2
         q1=q0+1
         qwt=q1+n
         pq=qwt+n
         dpq=pq+n*n
         ddpq=dpq+n*n            
         h0=1
         v0=h0+n*n
         hx=v0+n
         srf=hx+n*n
         q0=srf+2
         q=q0+1
         qwt=q+n
         pq=qwt+n
         dpq=pq+n*n
      elseif(typke.eq.'fd') then
         row=ndiff/2+1
         h0=1
         v0=h0+n*row
         q0=v0+n
         q=q0+1
         qwt=q+n
      else
         call lnkerr('bad ke call')
      endif
      if(typke.eq.'dvr'.or.typke.eq.'packed') then
         if(spac) then
            call rzero(ham(h0),n*n)
         endif
         if(prnt) then
            title='hamiltonian'
            call prntrm(title,ham(h0),n,n,n,n,iout)
            title='coordinates'
            call prntrm(title,ham(q0),n+1,1,n+1,1,iout)
            title='weights'
            call prntrm(title,ham(qwt),n,1,n,1,iout)
            title='polynomials'
            call prntrm(title,ham(pq),n,n,n,n,iout)
            title='first derivative of polynomials'
            call prntrm(title,ham(dpq),n,n,n,n,iout)
         endif
      else
         if(spac) then
            call rzero(ham(h0),row*n)
         endif
         if(prnt) then
            title='hamiltonian columns'
            call prntrm(title,ham(h0),row,n,row,n,iout)
            title='coordinates'
            call prntrm(title,ham(q0),n+1,1,n+1,1,iout)
            title='weights'
            call prntrm(title,ham(qwt),n,1,n,1,iout)
         endif
      endif
      return
      end

