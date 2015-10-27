c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Pert}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck pert.f
***begin prologue     pert
***date written       960723   (yymmdd)
***revision date      yymmdd   (yymmdd)
***keywords           time development
***author             schneider, barry (nsf)
***source             
***purpose            time-dependent potential
***description        calculate the time and space dependent potential
***                   matrix elements in the dvr representation.
***references         
***routines called    
***end prologue       pert 
      subroutine pert(v,vt,t,q1,q2,q3,typ,scale,omega,
     1                width,shift,dim,n,nd,tim,prn)
      implicit integer (a-z)
      real*8 v, vt, t, q1, q2, q3, tfac
      real*8 omega, scale, width, shift
      character*80 typ, title
      character*8 itoc
      logical prn
      dimension nd(3)      
      dimension v(n)
      dimension q1(nd(1)), q2(nd(2)), q3(nd(3))
      common/io/inp, iout
      write(iout,1)

c     Deal with pure time potentials first.  They just need to be added it as
c     they were computed in a previous call.

      call sadd(v,v,vt,n)

c     The current routine has a number of options for potentials.  They include
c     sine, cosine and dipole-like terms.  Others could be added by the user.

      if(typ.eq.'cosine') then
         tfac = scale*cos(omega*t)
         call sadd(v,v,tfac,n)         
      elseif(typ.eq.'sine') then
         tfac = scale*sin(omega*t)
         call sadd(v,v,tfac,n)         
      elseif(typ.eq.'t') then
         call sadd(v,v,t,n)         
      elseif(typ.eq.'cosine-dipole-field') then
         tfac = scale*cos(omega*t)
         call scdip(v,tfac,q1,q2,q3,dim,n,nd)
      elseif(typ.eq.'sine-dipole-field') then
         tfac = scale*sin(omega*t)
         call scdip(v,tfac,q1,q2,q3,dim,n,nd)
      elseif(typ.eq.'cosine-gaussian-pulse') then
         tfac= scale*cos(omega*t)*exp(-width*(t-shift)*(t-shift))
         call scdip(v,tfac,q1,q2,q3,dim,n,nd)
      elseif(typ.eq.'sine-gaussian-pulse') then
         tfac= scale*cos(omega*t)*exp(-width*(t-shift)*(t-shift))
         call scdip(v,tfac,q1,q2,q3,dim,n,nd)
      endif       
      if(prn) then
         title='time-dependent perturbation'
         call prntrm(title,v,n,1,n,1,iout)
      endif   
      return
 1    format(/,'calculating linear perturbing potential')
      end       


























