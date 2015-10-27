c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Eigen}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck eigen.f
***begin prologue     eigen
***date written       951229   (yymmdd)
***revision date      yymmdd   (yymmdd)
***keywords           
***author             schneider, barry (nsf)
***source             
***purpose            eigenvalues and eigenvectors
***                   
***references         
***routines called    
***end prologue       eigen
      subroutine eigen(pham,peigen,n,ndiff,prn)
      implicit integer (a-z)
      real*8 ham, eig, work
      logical prn
      dimension prn(2)
      character*80 title
      common/io/inp, iout
      pointer (pham,ham(1))
      pointer (peigen,eig(1))
      pointer (pwork,work(1))
      row=ndiff/2 + 1
      e=1
      ev=e+n
      need=wpadti(ev+n*n)
      call getmem(need,peigen,ngot,'eigen',0)

c     diagonalize saving unitary transformation matrix in eigvec

      if(ndiff.eq.3) then
         h=1
         scr=h+row*n
         need=wpadti(scr+5*n)
         call getmem(need,pwork,nwd,'work',0)
         call cpy3(ham,eig(e),work(h),n)
         call dstev('v',n,eig(e),work(h),eig(ev),n,work(scr),info)
      elseif(ndiff.eq.5) then
         h=1
         scr=h+row*n
         need=wpadti(scr+5*n)
         call getmem(need,pwork,nwd,'work',0)
         call copy(ham,work(h),n*row)
         call dsbev('v','l',n,2,work(h),3,eig(e),eig(ev),n,
     1               work(scr),info)
      elseif(ndiff.eq.7) then
         h=1
         scr=h+row*n
         need=wpadti(scr+5*n)
         call getmem(need,pwork,nwd,'work',0)
         call copy(ham,work(h),n*row)
         call dsbev('v','l',n,3,work(h),4,eig(e),eig(ev),n,
     1               work(scr),info)
      else
         need=wptoin(5*n)
         call getmem(need,pwork,nwd,'work',0)
         call copy(ham,eig(ev),n*n)
         title='hamiltonian in eigen'
         call prntrm(title,eig(ev),n,n,n,n,iout)             
         call dsyev('v','l',n,eig(ev),n,eig(e),work,5*n,info)
      endif
      if(info.eq.0) then
         if(prn(1)) then
            title='eigenvalues'
            call prntrm(title,eig(e),n,1,n,1,iout)
         endif 
         if(prn(2)) then
            title='eigenvectors'
            call prntrm(title,eig(ev),n,n,n,n,iout)
         endif
      else
         call lnkerr('error in diagonalization')
      endif 
      call getmem(-nwd,pwork,idum,'work',idum)
      return
      end       




