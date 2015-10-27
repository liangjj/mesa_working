c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Generalized Kinetic Energy Matrix Elements}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle    

*deck kemat.f
c***begin prologue     kemat
c***date written       022202   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           coordinate eigenfunctions
c***author             schneider, barry (nsf)
c***source
c***purpose            generate kinetic energy matrix elements with s
c***                   singularities removed and bloch operators added.
c***           
c***description
c***references         
c
c***routines called    
c***end prologue       kemat

c\begin{eqnarray}
c\end{eqnarray} 

      subroutine kemat(ke,bl,fa,dfb,ddfb,pt,wt,mass,typarg,
     #                 n,prn,sec)
      implicit integer (a-z)
      real*8 ke, bl, fa, dfb, ddfb, pt, wt, mass, scale
      character*(*) typarg
      character*80 title
      character*3 itoc
      logical prn
      dimension ke(n,n), bl(n,n), fa(n,n), dfb(n,n), ddfb(n,n)
      dimension pt(n), wt(n)
      common/io/inp, iout 
      call rzero(ke,n*n)
      if(typarg.eq.'linear') then
         do 10 i=1,n
            do 20 j=1,n
               ke(i,j) = ke(i,j) + fa(i,i)*ddfb(i,j)*wt(i)
 20         continue
 10      continue
      elseif(typarg.eq.'quadratic') then

c        In the cylindrical case, the functions are defined wrt the
c        argument $\rho^{2}$ and we need to use the chain rule on the 
c        derivatives to get them wrt $\rho$.  This is important in order to
c        get the limiting property at $\rho = 0$ correctly.

         do 30 i=1,n
            do 40 j=1,n
               ke(i,j) = ke(i,j)  + fa(i,i)*wt(i)*4.d0*
     1                            ( pt(i)*pt(i)*ddfb(i,j) + dfb(i,j) )
 40         continue
 30      continue 
      else
         call lnkerr('quit')
      endif
      scale= -.5d0/mass
      call vscale(ke,ke,scale,n*n)
      call vadd(ke,ke,bl,n*n)
      if(prn) then
         title='unnormalized kinetic energy matrix for '//
     #         'sector = '//itoc(sec)
         call prntrm(title,ke,n,n,n,n,iout)
      endif
      return
      end       
