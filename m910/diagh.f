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
          write(iout,1)
      endif
c

      if (givens.and.need.le.canget.and.(.not.dvd).and.
     1                                  (.not.bliu)) then
          call memory(need,p,ngot,'incore',0)
c
c           ----- n**3 diagonalization -----
c
          call incore(lnbuf,a(ibuf),z(rbuf),nwks,nroots,
     $                z(eigvec),z(eigval),z(t1),z(h),nnpwks,
     $                ops,rep,fzcore,prtflg,calc,dagtyp,headr)
          call memory(-ngot,p,idum,'incore',idum)
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
          call memory(need,p,ngot,'matrix',0)
c
c         get data for iterative eigenvalue calculation
c
          call dvddat(ops,nwks,nwksg,nroots,nattim,cnverg,thresh,
     1                mxiter,mxvec,left,prdvd,dvdall)     
c
c         we now know how many walks we can fit in the available memory.
c         go get it.
c
          nnpg=nwksg*(nwksg+1)/2
          ptgues=1
          hguess=iadtwp(ptgues+nwks)
          smvec=hguess
          smeig=hguess+nwksg*nwksg
          t1=smeig+nwksg
          needg=wpadti(t1+max(5*nwksg,nwks))
          call memory(needg,pdavd,ngotg,'guess',0)  
          nguess=intkey(ops,'ci=nguess',nwksg,' ')
          call iosys('open guess as scratch',0,0,0,' ')
c
          call iosys('read integer '//headr(2)//' from hamiltonian',1,
     1                ntotal,0,' ')
          call guess2(idavd(ptgues),nwks,a(ibuf),z(rbuf),lnbuf,
     $                davd(hguess),nwksg,z(diag),davd(t1),
     $                davd(smvec),davd(smeig),ops,ntotal,rep+fzcore,
     $                prtflg,dagtyp,headr)
          call memory(-ngotg,pdavd,idum,'guess',idum)  
c
          if(dvd) then          
c
c         now use the davidson routine to get the exact eigenpairs
c
             root=1
             dvdvec=root+mxiter
             dvdmat=dvdvec+mxiter**2
             c=dvdmat+mxiter*(mxiter+1)/2
             s=c+nwks*mxvec
             needvd=wpadti(s+nwks*mxvec)
c
             call memory(needvd,pdavd,ngotd,'dvdson',0)
c
c
             call david(davd(c),davd(s),nwks,nguess,mxvec,davd(root),
     $                  nroots,davd(dvdvec),davd(dvdmat),mxiter,
     $                  mxiter*(mxiter+1)/2,rep,fzcore,ops,a(ibuf),
     $                  z(rbuf),lnbuf,headr,prtflg,calc)
             call memory(-ngotd,pdavd,idum,'dvdson',idum)
          elseif(bliu) then
             trials = 1
             psi = trials + nguess*nwks
             vec = psi + nroots*nwks
             hvec = vec + mxvec*nwks
             b = hvec + mxvec*nwks
             btmp = b + mxvec*mxvec
             etmp = btmp + mxvec*mxvec
             svec = etmp + mxvec
             resid = svec + mxvec*mxvec
             eig = resid + max(mxvec*nwks,5*mxvec)
             t1 = eig + nroots
             t2 = t1 + mxvec*nwks
             needvb=wpadti(t2+mxvec*nwks) 
             call memory(needvb,pdavd,ngotd,'bliu',0)           
c
c
             call rsdvd(a(ibuf),z(rbuf),z(diag),davd(trials),davd(psi),
     1                  davd(vec),davd(hvec),davd(b),davd(btmp),
     2                  davd(etmp),davd(svec),davd(resid),davd(eig),
     3                  davd(t1),davd(t2),rep,cnverg,thresh,nwks,nroots,
     4                  nguess,nattim,mxiter,mxvec,lnbuf,headr,prdvd) 
             call memory(-ngotd,pdavd,idum,'bliu',idum)
          else
             call lnkerr('error in iterative matrix call')
          endif
c
          call iosys('destroy guess',0,0,0,' ')
          call memory(-ngot,p,idum,'matrix',idum)
        end if
      return
 1    format(/,5x,'in core diagonalization not possible.  try davidson')
      end

















