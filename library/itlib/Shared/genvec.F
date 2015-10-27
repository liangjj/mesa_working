c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{GENVEC: Generate New Trial Vectors based on Diagonal Preconditioning}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck genvec.f
c***begin prologue     genvec
c***date written       010829   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           davidson, trial, vector
c***author             schneider, barry (nsf)
c***source             
c***purpose            generate a new trial vector based on some zeroth
c***                   order model.  this is equivalent to a preconditioning
c***                   of the matrix.
c***references         
c
c***routines called    
c***end prologue       genvec
      subroutine genvec(vec,resid,diag,eig,work,n1,n2,n3,n,dim,
     1                  p,add,iter,precon,maxblk,prnt)
      implicit integer (a-z)
      real*8 vec, resid, diag, eig, work
      real*8 test, zero, nrzero, one, array
#ifdef DECPOINTER
      integer*8 p
#endif DECPOINTER
#ifdef SGIPOINTER
      integer*4 p
#endif SGIPOINTER
      logical prnt
      character*4 itoc
      character*80 title  
      character*(*) precon
      data zero, nrzero, one / 0.d0, 1.0d-06, 1.d0 /
      dimension vec(n,add), resid(n,add), diag(n), eig(add), work(*)
      common/io/inp, iout
      pointer (p,array(1)), (p,iarray(1))
      if(prnt) then
         title='input residuals to genvec'
         call prntrm(title,resid,n,add,n,add,iout)
      endif
      if(precon.eq.'diagonal') then	 
         do 10 i=1,add
            do 20 j=1,n
               test=eig(i) - diag(j)
               if(abs(test).ge.nrzero) then
                  vec(j,i) = resid(j,i)/test
               else
                  vec(j,i) = one
               endif
 20         continue
 10      continue   
c         title='vectors'
c         call prntrm(title,vec,n1,add,n1,add,iout)            
      elseif(precon.eq.'separable') then
         u1=1
         eig1=u1+n1*n1
         u2=eig1+n1
         eig2=u2+n2*n2
         u3=eig2+n2
         eig3=u3+n3*n3
         title='eigenvalues'
         call prntrm(title,array(eig1),n1,1,n1,1,iout)
	 title='eigenvalues'
	 call prntrm(title,eig,add,1,add,1,iout)
         title='transformation matrix'
         call prntrm(title,array(u1),n1,n1,n1,n1,iout)
         title='residuals'
         call prntrm(title,resid,n1,add,n1,add,iout)	 
         if(dim.eq.1) then
            call h12h0(resid,work,array(u1),n1,add)
	    title='h12h0 work'
	    call prntrm(title,work,n1,add,n1,add,iout)
            call h1e(work,array(eig1),eig,n1,add)
	    title='h1e work'
	    call prntrm(title,work,n1,add,n1,add,iout)
            call h12dvr(work,vec,array(u1),n1,add)
            title='vectors'
            call prntrm(title,vec,n1,add,n1,add,iout)            
         elseif(dim.eq.2) then
            call h22h0(resid,work,array(u1),array(u2),n1,n2,add)
            call h2e(work,array(eig1),array(eig2),eig,n1,n2,add)
            call h22dvr(work,vec,array(u1),array(u2),n1,n2,add)
         elseif(dim.eq.3) then
            call h32h0(resid,work,array(u1),array(u2),array(u3),
     1                 n1,n2,n3,add)
            call h3e(work,array(eig1),array(eig2),array(eig3),eig,
     1                n1,n2,n3,,add)
            call h32dvr(work,vec,array(u1),array(u2),array(u3),
     1                   n1,n2,n3,add)
         endif
      elseif(precon.eq.'block') then
         u0=1
         eig0=u0+maxblk*maxblk
         call rblvec(resid,vec,array(eig0),eig,array(u0),n,add)
      endif
      if(prnt) then
         title='new trial vectors iteration = '//itoc(iter)
         call prntrm(title,vec,n,add,n,add,iout)
      endif
      return
      end       





