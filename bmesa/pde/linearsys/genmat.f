*deck genmat.f
c***begin prologue     genmat
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***references         
c
c***routines called    
c***end prologue       genmat
      subroutine genmat(dim,itsolv,mattyp)
      implicit integer (a-z)
      real*8 z, scale, cnverg, thresh, eps, kappa
      character*80 cpass, title, chrkey
      character*128 filham
      character*1600 card
      character*1 itoc
      logical logkey, prlin, drctv, incore, dollar, itsolv
      logical prnton, prital, prbufh, cgrid
      character*8 tycalc
      character*(*) mattyp
      dimension ia(1), prlin(11)
      common z(1)
      equivalence (ia(1),z(1))
      common/io/inp, iout
      common/memory/ioff
      scale=1.d0
      incore=.true.
      if ( dollar('$genmat',card,cpass,inp) ) then
         n=intkey(card,'matrix-size',1,' ')
         tycalc=chrkey(card,'matrix-type','input',' ')
         prnton=logkey(card,'print=on',.false.,' ')
      endif
      if(itsolv) then
         call lindat(card,cpass,n,nrhs,nattim,cnverg,thresh,
     1               niter,nvec,lenbuf,cgrid,prbufh,drctv,
     2               prlin,prital,filham)
         write(iout,1) n, nrhs, thresh, cnverg, niter,
     1                 nvec, mattyp, drctv         
      else
         write(iout,2) n, nrhs, mattyp
      endif
 1    format(/,1x,'inputting a matrix of dimension = ',i5,/,1x,
     1            'number of right hand sides      = ',i3,/,1x,
     2            'convergence                     = ',e15.8,/,1x,
     3            'overlap criterion               = ',e15.8,/,1x,
     4            'number of iterations            = ',i3,/,1x,
     5            'number of vectors               = ',i3,/,1x,
     6            'type matrix                     = ',a24,/,1x,
     7            'preconditioner                  = ',l1)
 2    format(/,1x,'inputting a matrix of dimension = ',i5,/,1x,
     1            'number of right hand sides      = ',i5,/,1x,
     2            'type matrix                     = ',a24)
      fac=1
      if(mattyp.eq.'complex'.or.mattyp.eq.'real-unsymmetric') then
         fac=2
      endif   
      ioff=1
      do 10 i=1,2
         if(itsolv) then
            hbuf=ioff
            ibuf=wpadti(hbuf+2*fac*lenbuf)
            diag=iadtwp(ibuf+2*lenbuf)
            ind=wpadti(diaga+fac*n)
            etrial=iadtwp(ind+n)
            trials=etrial+ntrial*fac
            pvec=trials+n*ntrial*fac
            hpvec=pvec+n*nvec*fac
            vec=hpvec+nvec*n*fac
            bmat=vec+n*nvec*fac
            bmatm=bmat+nvec*nvec*fac
            vecl=bmatm+nvec*nvec*fac
            vecr=vecl+nvec*nvec*fac
            eigtot=vecr+nvec*nvec*fac
            etmp=eigtot+n*fac
            work=etmp+n*fac
            mwork=10*nvec*fac
            work1=work+max(n*fac*nvec,mwork)
            words=work1+max(n*fac*nvec,2*fac*nvec)            
         else
            ham=add
            hamp=ham+n*n*fac
            hamm=hamp+n*n*fac
            hbufa=ham
            hbufb=ham
            diaga=ham
            ibufa=wpadti(ham)
            ibufb=wpadti(ham)
            eig=ham+n*n*fac
            rwork=eig+n*fac
            vecl=work+2*n
            if(mattyp.eq.'complex'.or.mattyp.eq.
     1                   'real-unsymmetric') then
               vecr=vecl+n*n*fac
               work=vecr+n*n*fac
               words=rwork+10*n*fac   
            else
               vecr=vecl
               rwork=vecl
               words=rwork+n*n*fac
            endif
         endif            
         words=wpadti(words)
         if (i.eq.1) then
             call iosys ('read integer maxsiz from rwf',1,
     1                    canget,0,' ')
             if (words.gt.canget) then
                 call lnkerr('not enough memory. will quit')
             endif
             call iosys ('write integer maxsiz to rwf',1,
     1                    words,0,' ')
         else
             call getscm(words,z,ngot,'poly',0)
         endif
 10   continue         
      if(tycalc.eq.'input') then
         call fromin(z(ham),z(ham),z(work),n,mattyp)    
      elseif(tycalc.eq.'feder') then
         call federin(z(ham),z(ham),n,mattyp)
      elseif(tycalc.eq.'rpa') then
         call genrpa(eps,kappa,z(q),n)
         call filrpa(z(ham),z(hamp),z(hamm),z(hbufa),z(diaga),
     1               ia(ibufa),z(hbufb),ia(ibufb),eps,
     2               kappa,z(q),lenbuf,isolv,n,nela,nelb)
         if(.not.itsolv) then
            call ebc(z(ham),z(hamm),z(hamp),n,n,n) 
         endif
      endif
      if(itsolv) then
         if(tycalc.eq.'rpa') then
            call gentrl(z(diaga),z(etrial),z(trials),z(work),
     1                  ia(ind),n,ntrial)
            call lnkerr('')
         endif
c            call rpa(z(hbuf),ia(ibuf),z(diag),z(trials),
c     1               z(etrial),z(vecl),z(vecr),
c     2               z(pvec),z(hpvec),z(vec),
c     3               z(bmat),z(bmatm),z(etmp),z(svec),z(sveca),
c     4               z(vec),z(eigtot),z(work),z(work1),
c     5               lwork,scale,cnverg,thresh,
c     6               n,nroots,ntrial,niter,nvec,
c     7               lenbuf,hamd,incore,nel,prdvd,mattyp)
      else
         call rdiag(z(ham),z(ham),z(eig),z(eig),z(vecl),z(vecl),
     1              z(vecr),z(work),lwork,z(rwork),
     2              n,n,mattyp)
         title='eigenvalues'
         if(mattyp.eq.'complex'.or.mattyp.eq.'real-unsymmetric') then
            call prntcm(title,z(eig),n,1,n,1,iout)
            title='left eigenvectors'
            call prntcm(title,z(vecl),n,n,n,n,iout)
            title='right eigenvectors'
            call prntcm(title,z(vecr),n,n,n,n,iout)
            title='overlap matrix'
            call cehbtc(z(ham),z(vecl),z(vecr),n,n,n)
            call prntcm(title,z(ham),n,n,n,n,iout)
         else
            call prntrm(title,z(eig),n,1,n,1,iout)
            title='eigenvectors'
            call prntcm(title,z(vecl),n,n,n,n,iout)
         endif
      endif
      return
      end       
