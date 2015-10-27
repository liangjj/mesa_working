c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Lanczos}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle    

*deck lancz.f
c***begin prologue     lancz
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           eigenvalues, eigenvectors
c***author             schneider, barry (nsf)
c***source
c***purpose            generate orthogonal polynomials using recursion.
c***                   using the generated alpha and beta coefficients find
c***                   the points and weights of the generalized gauss
c***                   quadrature.
c***                   
c***description
c***references
c***routines called    
c***end prologue       lancz

c\begin{eqnarray}
c  \beta_{j} P_{j}(x) &=& \big ( x - \alpha_{j} \big ) P_{j-1}(x) - \beta_{j-1} P_{j-2}(x) 
c                             \\ \nonumber
c  \alpha_{j} &=& \langle P_{j-1} \mid x \mid P_{j-1} \rangle \\ \nonumber
c  \beta_{j} &=& \langle P_{j-1} \mid x \mid P_{j} \rangle \\ \nonumber
c  \beta_{j} P_{j}^{\prime}(x) &=& \big ( x - \alpha_{j} \big ) P_{j-1}^{\prime}(x) 
c                             - \beta_{j-1} P_{j-2}^{\prime}(x) + P_{j-1}(x)
c                                 \\ \nonumber
c  \beta_{j} P_{j}^{\prime \prime}(x) &=& \big ( x - \alpha_{j} \big ) 
c                                      P_{j-1}^{\prime \prime}(x) 
c                            - \beta_{j-1} P_{j-2}^{\prime \prime}(x) 
c                            + 2 P_{j-1}^{\prime}(x)
c\end{eqnarray}

      subroutine lancz(v,arg,a,b,wt,refwt,mu,scr,n,iter)
      implicit integer (a-z)
      real*8 v, arg, a, b, wt, refwt, mu, scr, anorm, scaprd
      dimension v(n,0:iter), arg(n), scr(n), wt(n), refwt(n)
      dimension a(1:iter), b(1:iter)
      common/io/inp, iout

c     fill scr with unit vector and compute the weight integral

      call vfill(scr,1.d0,n )
      mu=scaprd(scr,scr,wt,refwt,n)
      write(iout,1) mu

c     normalize the input function

      anorm=scaprd(v(1,0),v(1,0),wt,refwt,n)
      anorm=sqrt(1.d0/anorm) 
      call sscal(n,anorm,v(1,0),1)

c     we now have the first function

      if (iter.gt.1) then

c         form argument times the first function

          call vmul(scr,arg,v(1,0),n)

c         calculate a(1)

          a(1)=scaprd(v(1,0),scr,wt,refwt,n)

c         form mat times the first function - a(1) times the first function
c         and store it in the next polynomial

          call vwxs(v(1,1),scr,v(1,0),a(1),-1,n)

c         calculate b(1)

          b(1)=sqrt( scaprd(v(1,1),v(1,1),wt,refwt,n) )

c         normalize the second polynomial

          call sscal(n,1.d0/b(1),v(1,1),1)  
      endif
      if (iter.gt.2) then
          do 10 i=2,iter

c            multiply the last calculated polynomial by mat

             call vmul(scr,arg,v(1,i-1),n)

c            orthogonalize it to the two previous polynomials
c            calculating a(i) as we go

             a(i)=scaprd(v(1,i-1),scr,wt,refwt,n)
             call vwxs(v(1,i),scr,v(1,i-1),a(i),-1,n)
             call vwxs(v(1,i),v(1,i),v(1,i-2),b(i-1),-1,n)

c            calculate b(i)

             b(i) = sqrt (scaprd(v(1,i),v(1,i),wt,refwt,n ) )
 
c            normalize the polynomial and we are done

             call sscal(n,1.d0/b(i),v(1,i),1)
 10       continue

c         calculate the last a

      endif
      return
 1    format(/,1x,'weight integral = ',e15.8)
      end       


