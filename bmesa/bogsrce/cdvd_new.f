*deck cdvd.f
c***begin prologue     cdvd
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           iterative eigenvalue
c***author             schneider, barry (nsf)
c***source             
c***purpose            iterative davidson eigenvalue for unsymmetric
c***                   real eigenvalue equation.
c***references         
c
c***routines called    
c***end prologue       cdvd
      subroutine cdvd(trials,etrial,ntrials,cnverg,thresh,n,nroots,
     1                maxit,drctv,nonsym,sblen,ibvec,jbvec,
     2                sbvec,sclen,icvec,jcvec,scvec,diag)
      implicit integer (a-z)
      complex*16 hbuf, diag
      complex*16 trials, etrial, pvec, hpvec, vec, svecl, svecr 
      complex*16 b, btmp, eig, eig01, eig02, eig03
      complex*16 u01, u02, u03
      complex*16 resid, t1, t2, cdotc, tmp
      real*8 error, cnverg, thresh, scale
      real*8 zero, one, nrzero
			integer sblen, sclen, ibvec, jbvec, icvec, jcvec
			real*8 sbvec, scvec
      logical prnt, drctv, nonsym, incore
      character*16 status
      character*3 itoc
      character*1 it
      character*80 title
      dimension diag(n), trials(n,ntrials)
      dimension etrial(ntrials), pvec(n,maxit), hpvec(n,maxit)
      dimension vec(n,maxit), b(maxit,maxit), btmp(maxit,maxit)
      dimension eig(maxit)
C     dimension eig01(n1), eig02(n2), eig03(n3)
C     dimension u01(n1,n1), u02(n2,n2), u03(n3,n3)
      dimension resid(n,maxit), t1(n,maxit), t2(n,maxit)
      dimension svecl(maxit,maxit), svecr(maxit,maxit)
      dimension prnt(11)
			dimension ibvec(sblen),jbvec(sblen),sbvec(sblen)
			dimension icvec(sclen),jcvec(sclen),scvec(sclen)
      common/io/inp, iout
      data zero, one/ 0.d0, 1.d0/
      data nrzero / 1.0d-06 /
      it=itoc(dim)
C     if(.not.incore) then
C        call iosys('read integer "'//it//
C    1              'd number of elements" from ham',1,ntot,0,' ')
C     endif

      scale=1.0
      if(drctv) then
         write(*,1) nroots, ntrials, maxit, cnverg
      else
         write(*,2) nroots, ntrials, maxit, cnverg
      endif
      if(prnt(1)) then
         title='trial eigenvalues' 
         call prntcm(title,etrial,ntrials,1,ntrials,1)
         title='trial vectors'
         call prntcm(title,trials,n,ntrials,n,ntrials)
      endif
      call cc2opy(trials,pvec,n*ntrials)
      write(*,3) scale
      nbeg=1
      nend=ntrials

      do 10 iter=1,maxit
c      
c        orthonormalize the new trials to the old vectors
c
         call cschmt(pvec,thresh,n,nbeg,nend,nout,.true.)

         if(nout.ne.0) then
            nend=nbeg+nout-1
            if(prnt(11)) then
               write(*,4) iter, nend
            endif
            if(prnt(3)) then
               call tstovl(pvec,pvec,n,nend,nonsym)
            endif
            if(prnt(4)) then
               title='vectors iteration = '//itoc(iter)
               call prntcm(title,pvec(1,nbeg),n,nout,n,nout)
            endif
c        operate with hamiltonian on these vectors
c
            print*,n,nbeg,nout,nend
            call honv2(n,sblen,ibvec,jbvec,sbvec,sclen,
     1                icvec,jcvec,scvec,pvec,hpvec,nbeg,nend)
            if(prnt(5)) then
               title='h on vectors iteration = '//itoc(iter)
               call prntcm(title,hpvec(1,nbeg),n,nout,n,nout)
            endif
c
c        update the current small hamiltonian matrix.
c
            do 20 i=1,nend
               do 30 j=nbeg,nend
                  b(i,j) = cdotc(n,pvec(1,i),1,hpvec(1,j),1)
                  b(j,i) = cdotc(n,pvec(1,j),1,hpvec(1,i),1)
 30            continue   
 20         continue
            do 40 i=1,nend   
               do 50 j=1,nend
                  btmp(i,j) = b(i,j)
 50            continue
 40         continue
            if(prnt(6)) then
               title='small matrix iteration = '//itoc(iter)
               call prntcm(title,b,nend,nend,maxit,maxit)
            endif   
c
c        diagonalize small matrix
c
            call rdiag(btmp,btmp,eig,eig,svecl,svecl,svecr,
     1                 maxit,nend,nonsym)
            if(prnt(11)) then
               do 60 i=1,nroots
                  tmp=scale*eig(i)
                  write(*,6) i, tmp
 60            continue
            endif               
c                
c
c        transform vectors and hamiltonian on vectors to new basis
c
            call cebcx(t1,n,pvec,n,svecr,maxit,n,nend,nend)
            if(prnt(7)) then
               title='transformed vectors iteration = '//itoc(iter)
               call prntcm(title,t1,n,nend,n,maxit)
            endif
c
c        calculate new hamiltonian on vectors
c
            call cebcx(t2,n,hpvec,n,svecr,maxit,n,nend,nend)
            if(prnt(8)) then
               title='transformed  h on vectors iteration = '
     1                //itoc(iter)
               call prntcm(title,t2,n,nend,n,maxit)
            endif
c   
c        calculate residuals for desired roots and test for convergence
c
            addvec=0
            ncon=0
            do 500 i=1,min(nroots,nend)
               do 600 j=1,n
                  resid(j,i) = t2(j,i) - eig(i)*t1(j,i) 
 600           continue
               error = scale*sqrt(cdotc(n,resid(1,i),1,resid(1,i),1) )
               if(error.le.cnverg) then
                  status='converged'
                  ncon=ncon+1
               else
                  status='unconverged'
                  addvec=addvec+1
                  call cc2opy(resid(1,i),resid(1,addvec),n)
                  eig(addvec) = eig(i)
               endif
               if(prnt(11)) then
                  write(*,5) i, error, status
               endif
 500        continue
            if(ncon.eq.nroots) then
               write(*,8) nend
               do 2000 i=1,min(nroots,nend)
                  error = scale*sqrt(cdotc(n,resid(1,i),1,
     1                                      resid(1,i),1) )
                  status='unconverged'
                  if(error.le.cnverg) then
                     status='converged'
                  endif
                  write(*,9) i, scale*eig(i), error, status
 2000          continue   
C              call iosys('write integer "size of davidson vector '//
C    1                    'space" to ham',1,nend,0,' ')
C              call cvscal(eig,eig,scale,nend)
C              call iosys('write real "davidson eigenvalues" to ham',
C    1                     2*nend,eig,0,' ')
C              call cvscal(eig,eig,1.d0/scale,nend)
C              call iosys('write real "davidson vectors" to ham',
C    1                     2*nend*n,pvec,0,' ')
               return
            endif
            if(addvec.eq.0) then
               write(*,7)
               return
            else               
               if(prnt(9)) then
                  title='residuals iteration = '//itoc(iter)
                  call prntcm(title,resid,n,addvec,n,maxit)
               endif
c
c        make sure we will not exceed maxit
c
               newnum = nend + addvec
               newnum = min(newnum,maxit)
               addvec = newnum - nend
c
c            prepare new trial vectors from residuals for unconverged
c            roots
c
               nbeg = nend + 1         
               call cnvec(eig,diag,resid,pvec(1,nend+1),eig01,eig02,
     1                    eig03,u01,u02,u03,t1,t2,dim,n1,n2,n3,n,
     2                    addvec,drctv)
               nend=nend+addvec
               if(prnt(10)) then
                  title='new trial vectors iteration = '//itoc(iter)
                  call prntcm(title,pvec(1,nbeg),n,addvec,n,maxit)
               endif
            endif
      else
         write(*,11)
         do 3000 i=1,nroots
            error = scale*sqrt(cdotc(n,resid(1,i),1,
     1                                resid(1,i),1) )
            status='unconverged'
            if(error.le.cnverg) then
               status='converged'
            endif
            write(*,9) i, scale*eig(i), error, status
 3000    continue   
         return
      endif
c           
 10   continue
      return
 1    format(/,1x,'davidson eigenvalue solver using preconditioning',
     1                                                  /,10x,
     2            'number of roots               = ',i4,/,10x,
     3            'number of trials              = ',i4,/,10x,
     4            'maximum number of iterations  = ',i4,/,10x,
     5            'convergence criterion for rms = ',e15.8)      
 2    format(/,1x,'davidson eigenvalue solver',/,10x,
     1            'number of roots               = ',i4,/,10x,
     2            'number of trials              = ',i4,/,10x,
     3            'maximum number of iterations  = ',i4,/,10x,
     4            'convergence criterion for rms = ',e15.8)      
 3    format(/,5x,'hamiltonian scale factor = ',e15.8)
 4    format(/,1x,'cycle = ',i3,5x,'size of vector space = ',i3) 
 5    format(/,5x,'root = ',i3,2x,'rms error = ',e15.8,2x,
     1            'status = ',a16)
 6    format(/,5x,'root = ',i3,2x,'energy = ',e15.8,1x,e15.8)
 7    format(/1x,'cannot add any more vectors.  will quit')
 8    format(/,1x,'all roots are converged.',/,1x,
     1            'number of davidson vectors = ',i5)
 9    format(/,5x,'root = ',i3,2x,'energy = ',e15.8,1x,e15.8,/,5x,
     1            'rms error = ',e15.8,1x,'status = ',a16)
 11   format(/1x,'no more orthonormal vectors can be added',/,1x,
     1           'summary of unconverged results' )
      end       
