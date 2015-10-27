*deck crmres.f
c***begin prologue     crmres
c***date written       980420   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           residual calculation
c***author             schneider, barry (nsf)
c***source             
c***                   
c***references         
c
c***routines called    
c***end prologue       crmres
      subroutine crmres(vec,hvec,trmat,eig,etmp,scale,cnverg,resid,
     1                 maxerr,t1,t2,n,m,nroot,ncon,addvec,
     2                 maxvec,it,prnt)
      implicit integer (a-z)
      complex*16 vec, hvec, trmat, eig, etmp, resid, t1, t2
      complex*16 cdotc, temp
      real*8 scale, cnverg, maxerr, err
      character*16 status
      character*80 title
      character*4 itoc
      logical prnt
      dimension vec(n,*), hvec(n,*), trmat(maxvec,*), eig(*), etmp(*)
      dimension resid(n,*), t1(n,*), t2(n,*), prnt(4)
      common/io/inp, iout
c
c        first transform the vectors to the new basis. put them in t1
c       
      call cebcx(t1,n,vec,n,trmat,maxvec,n,m,m)
c
c        second transform the effect of the hamiltonian on the basis
c        to the new basis.  put them in t2
c
      call cebcx(t2,n,hvec,n,trmat,maxvec,n,m,m)
      if(prnt(1)) then
         title='information for iteration = '//itoc(it)
         write(iout,1) title
      endif
      if(prnt(3)) then
         title='transformed vectors iteration = '//itoc(it)
         call prntcm(title,t1,n,nroot,n,maxvec,iout)
      endif
      if(prnt(4)) then
         title='hamiltonian on transformed vectors iteration = '
     1          //itoc(it)
         call prntcm(title,t2,n,nroot,n,maxvec,iout)
      endif
      do 10 i=1,nroot
         do 20 j=1,n
            resid(j,i) = t2(j,i) - eig(i)*t1(j,i)
 20      continue
 10   continue
      if(prnt(2)) then
         title='residuals iteration = '//itoc(it)
         call prntcm(title,resid,n,nroot,n,maxvec,iout)
      endif
      addvec=0
      ncon=0
      maxerr=0.d0
      do 30 i=1,nroot                         
         err = scale*sqrt (cdotc(n,resid(1,i),1,
     1                             resid(1,i),1) )
         temp=scale*eig(i)
         maxerr=max(err,maxerr)
         if(err.le.cnverg) then
            status='converged'
            ncon=ncon+1
            call iosys('write real "energy for root = '//itoc(i)
     1                 //'" to ham',2,temp,0,' ')
            call iosys('write real "vector for root = '//itoc(i)
     1                 //'" to ham',2*n,t1(1,i),0,' ')
            if(prnt(1)) then
               write(iout,2) i, temp, err, status
            endif                        
         else
            status='unconverged'
            addvec=addvec+1
            call cc2opy(resid(1,i),resid(1,addvec),n)
            etmp(addvec) = eig(i)
            if(prnt(1)) then
               write(iout,2) i, temp, err, status
            endif     
         endif
 30   continue
      return
 1    format(/,5x,a80)         
 2    format(/,5x,'root            = ',i4,/,5x,
     1            'davidson energy = ',e15.8,1x,e15.8,/,5x,
     2            'rms error       = ',e15.8,/,5x,
     3            'status          = ',a16)
      end

