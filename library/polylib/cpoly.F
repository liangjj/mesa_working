c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{ Coordinate Functions and Derivatives}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle    

*deck cpoly.f
c***begin prologue     cpoly
c***date written       022202   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           coordinate eigenfunctions
c***author             schneider, barry (nsf)
c***source
c***purpose            generate coordinate functions. 
c***           
c***description
c***references         
c
c***routines called    
c***end prologue       cpoly
      subroutine cpoly(cp,dcp,ddcp,pt,arg,n,npt,typ,prn)
      implicit integer (a-z)
      real*8 cp, dcp, ddcp, pt
      real*8 arg, fac, fac2, tmp
      character*(*) typ
      character*80 title
      logical prn
      dimension cp(npt,0:n), dcp(npt,0:n), ddcp(npt,0:n)
      dimension pt(npt), arg(npt)
      common/io/inp, iout 
      if(typ.eq.'quadratic') then
        call vmul(arg,pt,pt,npt)
      else
        call copy(pt,arg,npt)
      endif
      call lgngr(cp,dcp,ddcp,arg,arg,npt,npt,.false.,'all')
      if(prn) then
         title='coordinate function'
         call prntfm(title,cp(1,0),npt,n+1,npt,n+1,iout)   
         title='first derivative of coordinate function'
         call prntfm(title,dcp(1,0),npt,n+1,npt,n+1,iout)   
         title='second derivative of coordinate function'
         call prntfm(title,ddcp(1,0),npt,n+1,npt,n+1,iout)   
      endif 
      return
      end       
