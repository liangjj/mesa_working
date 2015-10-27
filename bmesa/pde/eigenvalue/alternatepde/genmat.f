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
      subroutine genmat(dim,itdiag,mattyp)
      implicit integer (a-z)
      real*8 z, scale, cnverg, thresh, eps, kappa
      character*80 cpass, title, chrkey
      character*128 filham
      character*1600 card
      character*1 itoc
      logical logkey, prdvd, drctv, incore, dollar, itdiag
      logical prnton, dvdall, prbufh, cgrid, keeprt, addvd
      character*8 tycalc
      character*(*) mattyp
      dimension prdvd(11)
      pointer(p,z(1)),(p,ia(1))
      common/io/inp, iout
      scale=1.d0
      incore=.true.
      if ( dollar('$genmat',card,cpass,inp) ) then
         n=intkey(card,'matrix-size',1,' ')
         tycalc=chrkey(card,'matrix-type','input',' ')
         prnton=logkey(card,'print=on',.false.,' ')
         keeprt=logkey(card,'keep-lower-roots',.false.,' ') 
      endif
      if(itdiag) then
         call dvddat(card,cpass,n,nroots,ntrial,nattim,cnverg,thresh,
     1               niter,nvec,lenbuf,cgrid,prbufh,drctv,addvd,n0,
     2               prdvd,dvdall,filham)
         write(iout,1) n, ntrial, nroots, cnverg, thresh, niter,
     1                 nvec, mattyp, drctv, addvd         
      else
         nroots=intkey(ops,'number-of-roots',n,' ')
         nroots=min(nroots,n)
         write(iout,2) n, nroots, mattyp
      endif
      str=1
      if(mattyp.eq.'complex'.or.mattyp.eq.'real-unsymmetric') then
         str=2
      endif   
      ioff=1
      q=ioff
      if(tycalc.eq.'rpa') then
         add=q+n
      else
         add=q
      endif                                       
      if(itdiag) then
         hbufa=add
         ibufa=wpadti(hbufa+2*str*lenbuf)
         diaga=iadtwp(ibufa+2*lenbuf)
         eigtot=diaga+str*n
         wds=eigtot+str*n
         if(tycalc.eq.'rpa') then
            hbufb=wds
            ibufb=wpadti(hbufb+2*str*lenbuf)
            er=iadtwp(ibufb+2*lenbuf)
            ei=er+n 
            resid=ei+n
            wds=resid+n*nvec
         endif
         ind=wpadti(wds)
         etrial=iadtwp(ind+n)
         trials=etrial+ntrial*str
         psi=trials+n*ntrial*str
         pvec=psi+n*nroots*str
         hpvec=pvec+n*nvec*str
         vec=hpvec+nvec*n*str
         bmat=vec+n*nvec*str
         bmatm=bmat+nvec*nvec*str
         vecl=bmatm+nvec*nvec*str
         vecr=vecl+nvec*nvec*str
         etmp=vecr+nvec*nvec*str
         work=etmp+n*str
         lwork=16*nvec*str
         work1=work+max(n*str*nvec,lwork)
         words=work1+max(n*str*nvec,2*str*nvec)            
      else
         ham=add
         hamp=ham+n*n*str
         hamm=hamp+n*n*str
         hbufa=ham
         hbufb=ham
         diaga=ham
         ibufa=wpadti(ham)
         ibufb=wpadti(ham)
         if(tycalc.eq.'rpa') then
            er=ham+n*n*str
            ei=er+n
            lwork=16*n 
            work=ei+n
            vecl=work+lwork
            vecr=vecl+n*n
            words=vecr+n*n
         else
            eig=ham+n*n*str
            lwork=32*n
            work=eig+n*str
            vecl=work+lwork
            vecr=vecl+n*n*str
            rwork=vecr+n*n*str
            words=rwork+n
         endif
      endif            
      words=wpadti(words)
      call manmem(0,idum,idum,'genmat',idum)
      call manmem(words,p,ngot,'genmat',0)
      if(tycalc.eq.'input') then
         call fromin(z(ham),z(ham),z(work),n,mattyp)    
      elseif(tycalc.eq.'feder') then
         call federin(z(ham),z(ham),n,mattyp)
      elseif(tycalc.eq.'rpa') then
         call genrpa(eps,kappa,z(q),n)
         call filrpa(z(hamp),z(hamm),z(hbufa),z(diaga),
     1               ia(ibufa),z(hbufb),ia(ibufb),z(etrial),eps,
     2               kappa,z(q),lenbuf,itdiag,n,nela,nelb,.false.)
         if(itdiag) then
            call gentrl(z(etrial),z(trials),z(work),
     1                  ia(ind),n,ntrial)
            call rpa(z(hbufa),ia(ibufa),z(diaga),z(hbufb),ia(ibufb),
     1               z(trials),z(etrial),z(vecl),z(vecr),z(psi),
     2               z(pvec),z(hpvec),z(vec),z(bmat),z(bmatm),z(etmp),
     3               z(resid),z(er),z(ei),z(work),z(work1),
     4               lwork,cnverg,thresh,n,nroots,ntrial,nattim,
     5               niter,nvec,lenbuf,hamd,incore,nela,nelb,
     6               prdvd,mattyp)
         else
            call ebc(z(ham),z(hamm),z(hamp),n,n,n) 
            call rpadiag(z(ham),z(er),z(ei),z(vecl),z(vecr),
     1                   z(work),n,n,lwork)
            title='real part of eigenvalues'
            call prntrm(title,z(er),n,1,n,1,iout)
            title='imaginary part of eigenvalues'
            call prntrm(title,z(ei),n,1,n,1,iout)
            title='left eigenvectors'
            call prntrm(title,z(vecl),n,n,n,n,iout)
            title='right eigenvectors'
            call prntrm(title,z(vecr),n,n,n,n,iout)
            title='overlap matrix'
            call ebtc(z(ham),z(vecl),z(vecr),n,n,n)
            call prntrm(title,z(ham),n,n,n,n,iout)
         endif
      endif
      call manmem(-ngot,p,idum,'genmat',idum)
      return
 1    format(/,1x,'inputting a matrix of dimension = ',i5,/,1x,
     1            'number of trial vectors         = ',i3,/,1x,
     2            'number of roots                 = ',i3,/,1x,
     3            'convergence                     = ',e15.8,/,1x,
     4            'overlap criterion               = ',e15.8,/,1x,
     5            'number of iterations            = ',i3,/,1x,
     6            'number of vectors               = ',i3,/,1x,
     7            'type matrix                     = ',a24,/,1x,
     8            'preconditioner                  = ',l1,/,1x,
     9            'add diagonal of v               = ',l1)
 2    format(/,1x,'inputting a matrix of dimension = ',i5,/,1x,
     1            'number of roots                 = ',i5,/,1x,
     2            'type matrix                     = ',a24)
      end       
