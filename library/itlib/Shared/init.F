c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{INIT: Driver for Initialization of Davidson Vectors/Matrix}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck init.f
c***begin prologue     init
c***date written       010829   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           initialize, davidson
c***purpose            driver for initialization of davidson vectors
c***                   and matrices.   
c***description
c***author             schneider, barry (nsf)
c***source             
c***                   
c***references         
c
c***routines called    
c***end prologue       init
      subroutine init(hx,hy,hz,diag,vec,hvec,b,bwrk,scr,thresh,
     1                drctv,dim,n,nin,begin,size,nx,ny,nz,rtdone,
     2                code,maxvec,orth,prnt)
      implicit integer (a-z)
      real*8 hx, hy, hz, diag, vec, hvec, b, bwrk, scr
      real*8 thresh
      character*80 title
      character*4 itoc
      character*(*) code, drctv
      logical orth, prnt
      dimension hx(nx,nx), hy(ny,ny), hz(nz,nz), diag(n)
      dimension vec(n,*), hvec(n,*), b(maxvec,*), bwrk(maxvec,*)
      dimension scr(*)
      dimension prnt(3)
      common/io/inp, iout
c
c     get the additional set of vectors by orthogonalizing the trials
c     to the set of roots which have been converged.  this will ensure
c     that the subspace is orthogonal to previously converged vectors. 
c
      nout = nin
      if(orth) then
         if(rtdone.ne.0) then
            call invec(scr,code,n,rtdone,prnt(1))
            call abschm(scr,vec(1,begin),thresh,n,rtdone,nin,
     1                  nout,.true.,.false.)
         endif
      endif
      end = begin + nout - 1
      if(nout.ne.0) then
         call gschmt(vec,thresh,n,begin,end,nout,.true.,prnt(2))
      endif
      if(nout.eq.0) then
         call lnkerr('quit davidson. no more trial vectors '//
     1               'possible')
      endif
      size = begin + nout - 1                  
c
c     initialize the effect of the hamiltonian on these vectors.
c
      title='h on initial vectors'
      call hvdvr(hx,hy,hz,diag,vec(1,begin),hvec(1,begin),
     1           n,nx,ny,nz,nout,dim,prnt(3))
c
c     initialize the small hamiltonian matrix.
c            
      call hsmall(b,bwrk,vec,hvec,n,begin,size,maxvec,drctv,.false.)
      return
      end

