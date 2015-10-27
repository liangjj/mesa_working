*deck @(#)diagh.f	1.4  8/3/91
      subroutine diagh(calc,rep,fzcore,lnbuf,nroots,nwks,ops,headr)
c
c***begin prologue     diagh
c***date written       000811   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           ci, diagonalize
c***author             saxe, paul (lanl)
c***source             @(#)m940.f	1.4   8/3/91
c
c***purpose            to diagonalize either directly or iteratively
c***                   a hamiltonian matrix
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       m940
c
c
      implicit integer (a-z)
c
      character*(*) calc, ops, headr
      character*4 code
      character*8 prtflg
      integer a
      logical logkey, dvd, bliu, givens, prdvd, dvdall 
      real*8 z, davd
      real*8 rep, thresh, cnverg
      real*8 fzcore
      real*8 cutoff
      real*8 fpkey
      dimension headr(3), prdvd(11)
      pointer(p,a(1)),(p,z(1))
      pointer(pdavd,davd(1)), (pdavd,idavd(1))
c
      common /io/ inp,iout
      if(headr(1).eq.'"p space buffers"') then
         code='p-ci'
      elseif(headr(1).eq.'"q space buffers"') then
         code='q-ci'
      else
         code='ci'
      endif
c
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
c
      dvd=logkey(ops,'ci=davidson',.false.,' ')
      bliu=logkey(ops,'ci=bliu',.false.,' ')
      givens=logkey(ops,'ci=givens',.true.,' ')
c
c     ----- allocate core -----
c
c
c
c     ----- try for in-core diagonalization if possible -----
c
      h=1
      eigvec=h
      eigval=h+nwks*nwks
      t1=eigval+nwks
      ibuf=wpadti(t1+5*nwks)       
      rbuf=iadtwp(ibuf+2*lnbuf)
      need=wpadti(rbuf+lnbuf)
      nnpwks=nwks*(nwks+1)/2
c
      call iosys('read integer mxcore from rwf',1,canget,0,' ')
      if (need.gt.canget) then
          givens=.false.  
          write(iout,1) canget, need
      endif
c
      if (givens.and.need.le.canget.and.(.not.dvd).and.
     1                                  (.not.bliu)) then
          call getmem(need,p,ngot,'incore',0)
c
c           ----- n**3 diagonalization -----
c
          call incore(lnbuf,a(ibuf),z(rbuf),nwks,nroots,
     $                z(eigvec),z(eigval),z(t1),z(h),nnpwks,
     $                ops,rep,fzcore,prtflg,calc,dagtyp,headr,code)
          call getmem(-ngot,p,idum,'incore',idum)
      else
c
c           ----- out-of-core davidson diagonalization -----
c           ----- use either older routine or newer davidson -----
c
c          basic core allocation needed to read in matrix elements
c
          ibuf=1
          rbuf=iadtwp(ibuf+2*lnbuf)
          diag=rbuf+lnbuf
          need=wpadti(diag+nwks)
          left = canget - need
          call getmem(need,p,ngot,'matrix',0)
c
c         get data for iterative eigenvalue calculation
c
          call dvddat(ops,nwks,nwksg,nroots,nattim,cnverg,thresh,
     1                mxiter,left,prdvd,dvdall)     
c
c         we now know how many walks we can fit in the available getmem.
c         go get it.
c
          nnpg=nwksg*(nwksg+1)/2
          ptgues=1
          hguess=iadtwp(ptgues+nwks)
          smvec=hguess
          smeig=hguess+nwksg*nwksg
          t1=smeig+nwksg
          needg=wpadti(t1+max(5*nwksg,nwks))
          call getmem(needg,pdavd,ngotg,'guess',0)  
          nguess=intkey(ops,'ci=nguess',nroots,' ')
          nguess=max(nguess,nroots)
          nguess=min(nwksg,nguess) 
          call iosys('open guess as scratch',0,0,0,' ')
c
          call iosys('read integer '//headr(2)//' from hamiltonian',1,
     1                ntotal,0,' ')
          call guess2(idavd(ptgues),nwks,a(ibuf),z(rbuf),lnbuf,
     $                davd(hguess),nwksg,z(diag),davd(t1),
     $                davd(smvec),davd(smeig),ops,ntotal,rep+fzcore,
     $                prtflg,dagtyp,headr)
          call getmem(-ngotg,pdavd,idum,'guess',idum)  
c
          if(dvd) then          
c
c         now use the davidson routine to get the exact eigenpairs
c
             root=1
             dvdvec=root+mxiter
             dvdmat=dvdvec+mxiter**2
             c=dvdmat+mxiter*(mxiter+1)/2
c
c            this routine does a lot of disk storage and does not use much
c            memory.  we need the buffers for the non-zero hamiltonian matrix
c            elements and a few small matrices.  the remaining memory holds
c            the vectors and the effect of the hamiltonian on the vectors.
c            we now calculate how many we can store, get the memory and go.
c
             avail = iadtwp(left) - c
             mxvec=avail/(nwks*2)
             mxvec=min(mxvec,nwks)
             s=c+nwks*mxvec
             needvd=wpadti(s+nwks*mxvec)
c
             call getmem(needvd,pdavd,ngotd,'dvdson',0)
c
c
             call david(davd(c),davd(s),nwks,nguess,mxvec,davd(root),
     $                  nroots,davd(dvdvec),davd(dvdmat),mxiter,
     $                  mxiter*(mxiter+1)/2,rep,fzcore,ops,a(ibuf),
     $                  z(rbuf),lnbuf,headr,prtflg,calc,code)
             call getmem(-ngotd,pdavd,idum,'dvdson',idum)
          elseif(bliu) then
c
c            we absolutely need in memory tonow words
c
             mxvec = mxiter
             vec=1
             hvec=vec+nwks*mxiter
             b = hvec + mxiter*nwks
             btmp = b + mxiter*mxiter
             etmp = btmp + mxiter*mxiter
             work = etmp + mxiter
             tonow = work + 5*mxiter
c
c            whats left after buffers and diag
c
             avail = iadtwp(left) - tonow
c
c            whats needed
c
             remain = mxiter*mxiter + mxiter*nwks + nroots 
             if(remain.le.avail) then
                svec = tonow
                resid = svec + mxiter*mxiter
                eig = resid + nwks*mxiter
                needvb=wpadti(eig+nroots) 
                call getmem(needvb,pdavd,ngotd,'bliu',0)           
             else
                call lnkerr('not enough memory for rsdvd')
             endif
c
c
             call rsdvd(a(ibuf),z(rbuf),z(diag),davd(eig),
     1                  davd(vec),davd(hvec),davd(resid),
     2                  davd(b),davd(btmp),davd(etmp),davd(work),
     3                  davd(svec),rep,cnverg,thresh,nwks,
     4                  nroots,nguess,nattim,mxiter,mxvec,
     5                  lnbuf,headr,prdvd,code) 
             call getmem(-ngotd,pdavd,idum,'bliu',idum)
          else
             call lnkerr('error in iterative matrix call')
          endif
c
          call iosys('destroy guess',0,0,0,' ')
          call getmem(-ngot,p,idum,'matrix',idum)
      endif
      return
 1    format(/,5x,'in core diagonalization not possible.'/,5x,
     1            'canget = ',i8,1x,'need = ',i8,/,5x,
     2            'try davidson')
      end

















