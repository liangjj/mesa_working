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
      subroutine kemel(kmat,p,dp,ddp,pt,wt,edge,n,
     1                 cordsys,coord,parity)
      implicit integer (a-z)
      real*8 p, dp, ddp, pt, wt, kmat, edge
      character*(*) cordsys, coord, parity
      character*80 title
      dimension p(n,n), dp(n,n), ddp(n,n), kmat(n,n)
      dimension pt(n), wt(n), edge(2)
      common/io/inp, iout 
      call rzero(kmat,n*n)
      if(cordsys.eq.'cartesian') then
         if(coord.eq.'x'.or.
     1      coord.eq.'y'.or.
     2      coord.eq.'z') then

c The kinetic energy operator is assumed to be,
c \begin{equation}
c         T = - \frac{1}{2M} \frac{d^2}{dx^2}
c \end{equation}
c $x$ can be the transformed radial equation as well as cartesian coordinates.

            call kinxyz(kmat,p,ddp,wt,n)
            call blxyz(kmat,p,dp,pt(1),edge(1),pt(n),edge(2),n)
         elseif(coord.eq.'r') then
            call kinxyz(kmat,p,ddp,wt,n)
            call blxyz(kmat,p,dp,pt(1),edge(1),pt(n),edge(2),n)
         else
            call lnkerr('error in axis type')
         endif
      elseif(cordsys.eq.'spherical') then

c The kinetic energy operator is assumed to be,
c \begin{equation}
c    T =- \frac{1}{2M} \frac{1}{r^2} \frac{d}{dr} \big (  r^2 \frac{d} {dr} \big )
c \end{equation}

         call kinrad(kmat,p,dp,ddp,pt,wt,n,parity)
         call blrad(kmat,p,dp,pt(1),edge(1),pt(n),edge(2),n,parity)
      elseif(cordsys.eq.'cylindrical') then

c The kinetic energy operator is assumed to be,
c \begin{equation}
c    T =- \frac{1}{2M} \frac{1}{\rho} \frac{d}{d \rho} \big (  \rho \frac{d} {d \rho} \big )
c \end{equation}

         call kincyl(kmat,p,dp,ddp,pt,wt,n,parity)
         call blcyl(kmat,p,dp,pt(1),edge(1),pt(n),edge(2),n,parity)
      else
         call lnkerr('quit')
      endif
      return
      end       
