c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{FRMRES: Form Residuals and Transform to Optimal Vectors}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck frmres.f
c***begin prologue     frmres
c***date written       010829   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           residual calculation
c***author             schneider, barry (nsf)
c***source             
c***                   
c***references         
c
c***routines called    
c***end prologue       frmres
      subroutine frmres(vec,hvec,resid,trmat,eig,eigwrk,b,bwrk,rep,
     1                  cnverg,maxerr,n,m,nroot,con,uncon,
     2                  maxvec,it,prnt,point,code)
      implicit integer (a-z)
      real*8 vec, hvec, resid, trmat, eig, eigwrk, b, bwrk
      real*8 rep, cnverg, sdot, maxerr, err, temp
      character*16 status
      character*80 title
      character*4 itoc
      character*(*) code
      logical prnt
      dimension vec(n,*), hvec(n,*), resid(n,*), trmat(maxvec,*)
      dimension eig(*), eigwrk(*), b(maxvec,*), bwrk(maxvec,*)
      dimension prnt(4), code(2)
      common/io/inp, iout
c
c        first transform the vectors to the new basis. 
c       
      call ebcxx(resid,vec,trmat,n,m,m,n,n,maxvec)
      call copy(resid,vec,n*m)
c
c        second transform the effect of the hamiltonian on the basis
c        to the new basis.  
c
      call ebcxx(resid,hvec,trmat,n,m,m,n,n,maxvec)
      call copy(resid,hvec,n*m)
      if(prnt(1)) then
         title='information for iteration = '//itoc(it)
         write(iout,1) title
      endif
c
c     form residuals for the desired targeted nroots.
c     
      n2chk=min(m,nroot)
      do 10 i=1,m
         do 20 j=1,n
            resid(j,i) = hvec(j,i) - eigwrk(i)*vec(j,i)
 20      continue
 10   continue
      if(prnt(2)) then
         title='transformed vectors iteration = '//itoc(it)
         call prntfm(title,vec,n,n2chk,n,maxvec,iout)
      endif
      if(prnt(3)) then
          title='hamiltonian on transformed vectors iteration = '
     1          //itoc(it)
         call prntfm(title,hvec,n,n2chk,n,maxvec,iout)
      endif
      if(prnt(4)) then
         title='residuals iteration = '//itoc(it)
         call prntfm(title,resid,n,n2chk,n,maxvec,iout)
      endif      
c
c     re-constitute the small matrix
c
      call rzero(b,maxvec*maxvec)
      call rzero(bwrk,maxvec*maxvec)
      do 30 i=1,m
         b(i,i) = eigwrk(i)
         bwrk(i,i) = eigwrk(i)
 30   continue   
c
c     check converged and unconverged roots
c
      uncon=0
      con=0
      maxerr=0.d0
      do 40 i=1,n2chk
         err = sqrt (sdot(n,resid(1,i),1,
     1                      resid(1,i),1) )
         temp=eigwrk(i) + rep
         maxerr=max(err,maxerr)
         if(err.le.cnverg) then
            status='converged'
            con=con+1
            eig(i+point)=eigwrk(i)
c
c           write out converged eigenpairs to rwf
c
            call iosys('write real "'//code(1)//itoc(point+i)
     1                              //'" to rwf',1,temp,0,' ')
            call iosys('write real "'//code(2)//itoc(point+i)
     1                              //'" to rwf',n,vec(1,i),0,' ')
            write(iout,2) i, temp, err, status
         else
c
c           move eigenvalues and residuals of unconverged eigenpairs
c           
            status='unconverged'
            uncon=uncon+1
            call copy(resid(1,i),resid(1,uncon),n)
            eigwrk(uncon) = eigwrk(i)
            write(iout,2) i, temp, err, status
         endif
 40   continue
      return
 1    format(/,5x,a80)         
 2    format(/,5x,'root            = ',i4,/,5x,
     1            'davidson energy = ',f15.8,/,5x,
     2            'rms error       = ',f15.8,/,5x,
     3            'status          = ',a16)
      end

