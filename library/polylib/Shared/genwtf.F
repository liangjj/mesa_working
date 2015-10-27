c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Generalized Weight Functions for Orthogonal Functions}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle    

*deck genwtf.f
c***begin prologue     genwtf
c***date written       022202   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           coordinate eigenwtfctions
c***author             schneider, barry (nsf)
c***source
c***purpose            generate weight functions and derivatives on grid.
c***           
c***description
c***references         
c
c***routines called    
c***end prologue       genwtf

      subroutine genwtf(pt,arg,wt,dwt,ddwt,alpha,beta,edge,npt,
     1                  wtfn,prn)
      implicit integer (a-z)
      real*8 pt, arg, wt, dwt, ddwt, alpha, beta, edge
      character*(*) wtfn
      character*80 title
      logical prn
      dimension pt(npt), arg(npt), wt(npt), dwt(npt), ddwt(npt)
      dimension edge(2)
      common/io/inp, iout 

c     Get the weight functions and their derivatives

      write(iout,1) wtfn
      call genrwt(wt,dwt,ddwt,arg,wtfn,alpha,beta,.true.,
     1            edge,npt)
      if(prn) then
         title='weight function'
         call prntfm(title,wt,npt,1,npt,1,iout)   
         title='first derivative of weight function'
         call prntfm(title,dwt,npt,1,npt,1,iout)   
         title='second derivative of weight function'
         call prntfm(title,ddwt,npt,1,npt,1,iout)   
      endif
      return
 1    format(/,1x,'type weight function = ',a32)
      end       
