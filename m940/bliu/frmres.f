*deck frmres.f
c***begin prologue     frmres
c***date written       980420   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           residual calculation
c***author             schneider, barry (nsf)
c***source             
c***                   
c***references         
c
c***routines called    
c***end prologue       frmres
      subroutine frmres(vec,hvec,trmat,eig,etmp,b,btmp,rep,cnverg,resid,
     1                 maxerr,n,m,nroot,ncon,addvec,
     2                 maxvec,it,prnt,point,code,ilen)
      implicit integer (a-z)
      real*8 vec, hvec, trmat, eig, etmp, resid, b, btmp
      real*8 rep, cnverg, sdot, maxerr, err, temp
      character*16 status
      character*80 title
      character*4 itoc
      character*(*) code
      logical prnt
      dimension vec(n,*), hvec(n,*), trmat(maxvec,*)
      dimension eig(*), etmp(*), resid(n,*)
      dimension b(maxvec,*), btmp(maxvec,*), prnt(4)
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
      do 10 i=1,nroot
         do 20 j=1,n
            resid(j,i) = hvec(j,i) - etmp(i)*vec(j,i)
 20      continue
 10   continue
      if(prnt(2)) then
         title='residuals iteration = '//itoc(it)
         call prntfm(title,resid,n,nroot,n,maxvec,iout)
      endif
      if(prnt(3)) then
         title='transformed vectors iteration = '//itoc(it)
         call prntfm(title,vec,n,nroot,n,maxvec,iout)
      endif
      if(prnt(4)) then
          title='hamiltonian on transformed vectors iteration = '
     1          //itoc(it)
         call prntfm(title,hvec,n,nroot,n,maxvec,iout)
      endif
c
c     re-constitute the small matrix
c
      call rzero(b,maxvec*maxvec)
      call rzero(btmp,maxvec*maxvec)
      do 30 i=1,m
         b(i,i) = etmp(i)
         btmp(i,i) = etmp(i)
 30   continue   
      addvec=0
      ncon=0
      maxerr=0.d0
      do 40 i=1,min(nroot,m)
         err = sqrt (sdot(n,resid(1,i),1,
     1                      resid(1,i),1) )
         temp=etmp(i) + rep
         maxerr=max(err,maxerr)
         if(err.le.cnverg) then
            status='converged'
            eig(i+point)=etmp(i)
            ncon=ncon+1
            call iosys('write real "'//code(1:ilen)//' energy '
     1                              //itoc(point+i)
     2                             //'" to rwf',1,temp,0,' ')
            call iosys('write real "'//code(1:ilen)//' root '
     1                              //itoc(point+i)
     2                             //'" to rwf',
     3                                  n,vec(1,i),0,' ')
            write(iout,2) i, temp, err, status
         else
            status='unconverged'
            addvec=addvec+1
            call copy(resid(1,i),resid(1,addvec),n)
            etmp(addvec) = etmp(i)
            write(iout,2) i, temp, err, status
         endif
 40   continue
      write(iout,3) addvec
      return
 1    format(/,5x,a80)         
 2    format(/,5x,'root            = ',i4,/,5x,
     1            'davidson energy = ',f15.8,/,5x,
     2            'rms error       = ',f15.8,/,5x,
     3            'status          = ',a16)
 3    format(/,5x,'number of unconverged vectors = ',i5)
      end

