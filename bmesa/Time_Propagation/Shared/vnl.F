c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Vnl}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck vnl.f
***begin prologue     vnl
***date written       960723   (yymmdd)
***revision date      yymmdd   (yymmdd)
***keywords           copy
***author             schneider, barry (nsf)
***source             
***purpose            mean time dependent hartree field
***references         
***routines called    
***end prologue       vnl
      subroutine vnl(v,psi,gamma,p1,p2,p3,n,nd,dim,prn)
      implicit integer (a-z)
      real*8 v, gamma, p1, p2, p3, fi, fj, fk
      complex*16 psi
      character*80 title
      logical prn
      dimension nd(dim)
      dimension v(n), psi(n)
      dimension p1(nd(1),nd(1)), p2(nd(2),nd(2)), p3(nd(3),nd(3))
      common/io/inp, iout
      write(iout,1)
      if(dim.eq.1) then
         do 10 i=1,nd(1)
            v(i) = v(i) + gamma * psi(i)*conjg( psi(i) )*
     1                            p1(i,i)*p1(i,i)
 10      continue
      elseif(dim.eq.2) then
         cnt=0
         do 20 i=1,nd(1)
            fi = p1(i,i)*p1(i,i)
            do 30 j=1,nd(2)
               fj = p2(j,j)*p2(j,j)
               cnt=cnt+1 
               v(cnt) = v(cnt) + gamma * psi(cnt)*conjg( psi(cnt) ) *
     1                                          fi *fj
 30         continue   
 20      continue   
      elseif(dim.eq.3) then
         cnt=0 
         do 40 i=1,nd(1)
            fi = p1(i,i)*p1(i,i)
            do 50 j=1,nd(2)
               fj = p2(j,j)*p2(j,j)
               do 60 k=1,nd(3)
                  fk = p3(k,k)*p3(k,k)
                  cnt = cnt + 1
                  v(cnt) = v(cnt) + gamma * psi(cnt)*conjg( psi(cnt) ) *
     1                                          fi *fj *fk
 60            continue
 50         continue
 40      continue
      else
         call lnkerr('dimension error')  
      endif
      if(prn) then
         title='non-linear potential'
	 call prntrm(title,v,n,1,n,1,iout)
      endif	 
      return
 1    format(/,'calculating non-linear perturbing potential')
      end       


