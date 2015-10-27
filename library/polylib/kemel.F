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

*deck kemel.f
c***begin prologue     kemel
c***date written       022202   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           coordinate eigenfunctions
c***author             schneider, barry (nsf)
c***source
c***purpose            generate kinetic energy matrix elements with singularities
c***                   removed and bloch operators added.
c***           
c***description
c***references         
c
c***routines called    
c***end prologue       kemel
      subroutine kemel(p,dp,ddp,kmat,pt,wt,wtfn,typarg,edge,n)
      implicit integer (a-z)
      real*8 p, dp, ddp, pt, wt, wtfn, kmat, edge
      character*(*) typarg
      character*80 title
      dimension p(n,n), dp(n,n), ddp(n,n), kmat(n,n)
      dimension pt(n), wt(n), wtfn(n), edge(2)
      common/io/inp, iout 
      call rzero(kmat,n*n)
      if(typarg.eq.'linear') then
         do 10 i=1,n
            do 20 j=1,n
               kmat(i,j) = kmat(i,j) + 
     1                     p(i,i)*ddp(i,j)*wt(i)
 20         continue
 10      continue
         if(pt(1).eq.edge(1)) then
            do 30 j=1,n
               kmat(1,j) = kmat(1,j) + p(1,1)*dp(1,j)
 30         continue
         endif   
         if(pt(n).eq.edge(2)) then
            do 40 j=1,n
               kmat(n,j) = kmat(n,j) - p(n,n)*dp(n,j)
 40         continue   
         endif
      elseif(typarg.eq.'quadratic') then

c        In the cylindrical case, the functions are defined wrt the
c        argument $\rho^{2}$ and we need to use the chain rule on the 
c        derivatives to get them wrt $\rho$.  This is important in order to
c        get the limiting property at $\rho = 0$ correctly.

         do 50 i=1,n
            do 60 j=1,n
               kmat(i,j) = kmat(i,j)  + 
     1                     p(i,i)*wt(i)*4.d0*
     2                   ( pt(i)*pt(i)*ddp(i,j) 
     3                                + dp(i,j) )
 60         continue
 50      continue 
         title='before bloch'
         call prntrm(title,kmat,n,n,n,n,iout)        
         if(pt(1).eq.edge(1)) then
            do 70 j=1,n
               kmat(1,j) = kmat(1,j) + 
     1                           2.d0*pt(1)*pt(1)*p(1,1)*dp(1,j)
 70         continue
         endif   
         if(pt(n).eq.edge(2)) then
            do 80 j=1,n
               kmat(n,j) = kmat(n,j) - 
     1                           2.d0*pt(n)*pt(n)*p(n,n)*dp(n,j)
 80         continue   
         endif
         title='after bloch'
         call prntrm(title,kmat,n,n,n,n,iout)        
      else
         call lnkerr('quit')
      endif
      return
      end       
