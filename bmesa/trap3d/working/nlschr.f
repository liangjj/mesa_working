*deck nlschr.f
c***begin prologue     nlschr
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           non-linear schroedinger equation, gross-pitaevski
c***author             schneider, barry (nsf)
c***source             timeprp
c***purpose            calculate eigenvalues of gross-pitaevski equation
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       nlschr
      subroutine nlschr(eig,psi,hpsi,hbuf,ibuf,diag,f,trials,etrial,
     1                  pvec,hpvec,vec,vnl,vnlcur,p1,p2,p3,
     2                  eigr1,eigr2,eigr3,eig01,eig02,eig03,
     3                  u01,u02,u03,b,btmp,hh,pp,ph,hp,sol,
     4                  dvm,dvmtmp,t1,t2,
     5                  ipvt,natmin,natmax,norm0,
     6                  u0,scale,todiis,diiscn,mindiis,cnverg,thresh,n,
     7                  diiscy,diisit,maxit,nroots,ntrials,nstep,
     8                  ndim,n1,n2,n3,nc1,nc2,nc3,lenbuf,incore,
     9                  prdiis,prnt,drctv,davdir,itdiag,trunc,diisze,
     x                  restrt,nrestr,atom,usetf,plot,becin,becout)
      implicit integer (a-z)
      real*8 eig, psi, hpsi, hbuf, diag, f, trials, etrial, pvec
      real*8 hpvec, vec, vnl, vnlcur, p1, p2, p3
      real*8 eigr1, eigr2, eigr3, eig01, eig02, eig03
      real*8 u01, u02, u03
      real*8 b, btmp, hh, pp, ph, hp, sol
      real*8 dvm, dvmtmp, t1, t2, natmin, natmax, natoms
      real*8 todiis, diiscn, mindiis, cnverg
      real*8 thresh, error, test
      real*8 gamma, norm0, pi, zero, one, delat
      real*8 eig0, crit, eigtr, scale, iscale, u0, dum, nrestr
      character*80 title
      character*120 pltfle
      character*(*) atom, becin, becout
      character*3 itoc, ans
      character*1 it
      logical drctv, prdiis, prnt, incore, itdiag, restrt, trunc
      logical usetf, plot, davdir, flag
      dimension eig(n), psi(n,*), hpsi(n,*)
      dimension hbuf(*), ibuf(*), diag(*), f(*), trials(*), etrial(*)
      dimension pvec(*), hpvec(*), vec(*), vnl(n,*), vnlcur(*)
      dimension p1(n1,0:nc1-1), p2(n2,0:nc2-1)
      dimension p3(n3,0:nc3-1)
      dimension eigr1(nc1), eigr2(nc2), eigr3(nc3)
      dimension eig01(nc1), eig02(nc2), eig03(nc3)    
      dimension u01(nc1,nc1), u02(nc2,nc2), u03(nc3,nc3)
      dimension b(diisit+1,diisit+1), btmp(diisit+1,diisit+1)
      dimension sol(diisit+1)
      dimension hh(diisit,diisit), pp(diisit,diisit), ph(diisit,diisit)
      dimension hp(diisit,diisit)
      dimension dvm(maxit,maxit), dvmtmp(maxit,maxit)
      dimension t1(*), t2(*), ipvt(n)
      common/io/inp, iout
      data pi/3.14159265358979323844d0/
      data zero, one/ 0.d0, 1.d0/
      if(itdiag) then
         it=itoc(ndim)
         call iosys('read integer "'//it//'d number of elements" '//
     1              'from ham',1,nel,0,' ')
      endif
      call iosys('write real scale to bec',1,scale,0,' ')
      nobj=n*n
      nerr=nobj
      iscale=1.d0/scale
      write(iout,1)
      write(iout,2) u0
      write(iout,3) diiscy, diisit, diiscn
      natoms=natmin
      delat=(natmax - natmin)/nstep
      if(restrt) then
c        
c           get last converged result
c      
         call iosys('read integer "number of files" from bec',1,
     1               nfiles,0,' ')
         offset = 2*(n+1)*(nfiles-1)
         call iosys('read real '//becin//' from bec '//
     1              'without rewinding',1,natoms,offset,' ')
         call iosys('read real '//becin//' from bec '//
     1              'without rewinding',1,eig(1),0,' ')
         call iosys('read real '//becin//' from bec '//
     1              'without rewinding',n,vnlcur,0,' ')
         call iosys('read real '//becin//' from bec '//
     1              'without rewinding',n,psi,0,' ')
         delat=(natmax-natoms)/nstep
         write(iout,100) natoms, eig(1)/iscale, delat
      else
         call copy(trials,psi,n)
         eig(1)=etrial(1)
      endif
      if(usetf) then
         write(iout,90)
         call tfsol(atom,dum,dum,dum,dum,dum,dum,n)
      endif         
      if(plot) then
         call iosys('does '//becin//' exist on bec',0,0,0,ans)
         if(ans.eq.'no') then
            call iosys('create real '//becin//' on bec',-1,0,0,' ')
            pltfle=becin
            dum=0.d0
            call iosys('write real '//pltfle//' to bec without '//
     1                  'rewinding',1,dum,0,' ')   
            write(iout,110) becin, pltfle
         else
            call iosys('does '//becout//' exist on bec',0,0,0,ans)
            if(ans.eq.'no') then
               call iosys('create real '//becout//' on bec',-1,0,0,' ')
               dum=0.d0
               call iosys('write real '//becout//' to bec without '//
     1                    'rewinding',1,dum,0,' ') 
            endif
            pltfle=becout
            write(iout,120) becin, becout
         endif               
      endif
      count=0
      do 1000 num=1,nstep
         crit = todiis
         test = crit
         gamma=natoms*u0
         gamma=gamma*norm0*norm0
         write(iout,15)
         if(usetf) then
            call tfsol(atom,psi,p1,eig(1),eigr1,natoms,norm0,n)
         endif   
         write(iout,4) natoms, gamma
         write(iout,5) eig(1)/iscale
         eig0=eig(1)
         do 10 cyc=1,diiscy
            begin=1
            do 20 iter=1,diisit
               call potnl(vnl(1,begin),eigr1,eigr2,eigr3,p1,p2,p3,
     1                    psi(1,begin),gamma,n,n1,n2,n3,
     2                    nc1,nc2,nc3,ndim)
               call fpsi(hbuf,ibuf,diag,vnl(1,begin),psi(1,begin),
     1                   hpsi(1,begin),test,n,n1,n2,n3,
     2                   lenbuf,nel,ndim,incore,itdiag,prdiis)
               error=test/iscale
               if(prnt) then
                  write(iout,80) iter, error
               endif
               call errfac(psi,hpsi,hh,pp,ph,hp,n,begin,diisit)
               if(drctv) then
c
c                     diis procedure when appropriate
c
c                    check current size of diis space and truncate
c                    if that option is turned on.
                     if(begin.gt.diisze) then
                        if(trunc) then
                           call trdiis(vnl,psi,hpsi,b,hh,pp,ph,hp,
     1                                 begin,diisze,n,diisit,prdiis)
                           begin=begin-1                           
                        endif
                     endif               
                 if(test.le.crit) then
                     call diis(vnlcur,vnl,psi,hpsi,
     1                         b,btmp,hh,pp,ph,hp,sol,test,
     2                         ipvt,begin,n,diisit,prdiis)
c                     error=test/iscale
                 endif
               else
                  call copy(vnl(1,begin),vnlcur,n)
               endif
               begin=begin+1
c
c              diagonalize fock matrix for the current iteration
c              and find the new wavefunction
c
               if(itdiag) then
                  call nldvd(hbuf,ibuf,diag,vnlcur,trials,etrial,
     1                       eig01,eig02,eig03,u01,u02,u03,
     2                       pvec,hpvec,vec,dvm,dvmtmp,vec,
     3                       eig,t1,t2,iscale,cnverg,thresh,n1,n2,n3,
     4                       ndim,n,nroots,ntrials,maxit,lenbuf,
     5                       nel,davdir,incore)
                  call copy(pvec,psi(1,begin),n)
                  call copy(eig,etrial,ntrials)
                  call copy(pvec,trials,n*ntrials)
               else 
                 call newfn(f,hbuf,vnlcur,eig,psi(1,begin),n)
               endif   
               if(prnt) then
                  write(iout,6) iter, eig(1)/iscale, error
               endif
               if(error.le.diiscn) then
                  eigtr=eig(1)/iscale
                  count=count+1
                  if(plot) then
                     call iosys('write real '//pltfle//' to bec '//
     1                          'without rewinding',1,natoms,0,' ')
                     call iosys('write real '//pltfle//' to bec '//
     1                          'without rewinding',1,eigtr,0,' ')
                     call iosys('write real '//pltfle//' to bec '//
     1                          'without rewinding',n,vnlcur,0,' ')
                     call iosys('write real '//pltfle//' to bec '//
     1                          'without rewinding',n,
     2                           psi(1,begin),0,' ')
                  endif
c                  table(count,1)=natoms
c                  table(count,2)=eigtr
c                  table(count,3)=error
                  write(iout,7) iter, natoms, eig(1), eigtr
                  if(prdiis) then
                     title='full spectrum'
                     call prntfm(title,eig,n,1,n,1,iout)
                     title='final wavefunction'
                     call prntfm(title,t1,n,1,n,1,iout)
                  endif
                  go to 1001
               endif
               eig0=eig(1)
 20         continue
            if(prnt) then
               write(iout,8) cyc, iter
            endif
            if(cyc.eq.diiscy) then
               if(error.le.mindiis) then
                  eigtr=eig(1)/iscale
                  count=count+1
                  if(plot) then
                     call iosys('write real '//pltfle//' to bec '//
     1                          'without rewinding',1,natoms,0,' ')
                     call iosys('write real '//pltfle//' to bec '//
     1                          'without rewinding',1,eigtr,0,' ')
                     call iosys('write real '//pltfle//' to bec '//
     1                          'without rewinding',n,vnlcur,0,' ')
                     call iosys('write real '//pltfle//' to bec '//
     1                          'without rewinding',n,
     2                           psi(1,begin),0,' ')
                  endif
c                  table(count,1)=natoms
c                  table(count,2)=eigtr
c                  table(count,3)=error
                  write(iout,7) iter, natoms, eig(1), eigtr
                  if(prdiis) then
                     title='full spectrum'
                     call prntfm(title,eig,n,1,n,1,iout)
                     title='final wavefunction'
                     call prntfm(title,t1,n,1,n,1,iout)
                  endif
                  go to 1001
               endif      
               if(prnt) then
                  write(iout,9) cyc
               endif
               go to 30
            else
               call copy(psi(1,begin),psi,n)   
            endif
 10      continue   
 1001    continue
         natoms = natoms + delat
         write(iout,15)
 1000 continue   
 30   continue
c      do 2000 i=1,count
c         write(iout,70) table(i,1), table(i,2), table(i,3)
c 2000 continue   
      if(plot) then
         call iosys('endfile '//pltfle//' on bec',0,0,0,' ')
         call iosys('rewind '//pltfle//' on bec read-and-write',
     1               0,0,0,' ')
         dum=count
         call iosys('write real '//pltfle//' to bec',1,dum,0,' ')
         call iosys('read real '//pltfle//' from bec',1,dum,0,' ')
         count=dum
         write(iout,130) pltfle, count
c         call iosys('write real "plot table" to bec',3*count,
c     1               table,0,' ')
      endif   
      return
    1 format(/,5x,'diagonalize time-independent non-linear schroedinger'
     1            ' equation')
    2 format(/,5x,'u0 = ',e15.8)
    3 format(/,5x,'number of diis cycles           = ',i4,/,5x,
     1            'number of diis iterations/cycle = ',i4,/,5x,
     2            'diis convergence criterion      = ',e15.8)
    4 format(/,5x,'number of atoms = ',e15.8,/,5x,
     1            'gamma           = ',e15.8)
    5 format(/,10x,'starting approximation',2x,'diis energy = ',e15.8)
    6 format(/,10x,'diis iteration  = ',i3,1x,'diis energy = ',e15.8,
     1       /,10x,'diis rms error = ',e15.8)         
    7 format(/,5x,'converged after ',i4,' diis iterations',/,5x,
     1            'number of atoms = ',e15.8,1x,'energy = ',e15.8,
     2       /,5x,'energy in trap units = ',e15.8)
    8 format(/,10x,'cycle = ',i3,/,5x,
     1             'no convergence after ',i4,' diis iterations')    
 9    format(/,5x,'quit calculation.  no convergence after cycle = ',i3)
 15   format(/,'*******************************************************'
     1         '****************')
 60   format(//,25x,'summary table',//,10x,'number of atoms',
     1           7x,'energy',12x,'convergence')
 70   format(7x,e15.8,5x,e15.8,5x,e15.8)
 80   format(/,5x,'iteration  = ',i3,1x,'rms error = ',e15.8)
 90   format(/,5x,'using thomas-fermi wavefunction as initial guess')
100   format(/,5x,'restarting calculation with:',
     1       /,5x,'number of atoms = ',i7,
     2       /,5x,'energy          = ',e15.8,
     3       /,5x,'atom step       = ',i5)
110   format(/,5x,'plot file          = ',a16,' does not exist',/,5x,
     1            'new plot file name = ',a16)
120   format(/,5x,'plot file          = ',a16,' exists',/,5x,
     1            'new plot file name = ',a16)
130   format(/,5x,'plot file = ',a60,/,5x,'number of files = ',i10)                         
      end       
