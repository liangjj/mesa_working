*deck @(#)directj.f	1.3  11/28/95
      subroutine directj(jmat,jprim,ptprim,noprim,nbtype,ex,c,
     $                   nx,ny,nz,lenxyz,nocart,mintyp,maxmom,
     $                   z,left,nat,npf,nnprim,nprim,ops,cutexp,
     $                   rhotest,a,prnt,ndmat,
     $                   pstart,prtoao,dpr,
     $                   nbf,nnp,ntotal,nonzer,nnshl,qint,dijmax,
     $                   qtest,calc,dok,kmat,kprim,doiijk,doiikk)
c
c***begin prologue     directj.f
c***date written       931215   (yymmdd)
c***revision date      4/18/95
c
c  22 June 1995        russo at lanl
c     add kmatrix functionality
c  18 december, 1993   rlm at lanl
c     revising the driver routine in m712 for direct j-matrix.
c***keywords           direct, j-matrix
c***author             martin, saxe, lengsfield, russo
c***source             @(#)directj.f	1.3   11/28/95
c***purpose            calculates two-electron integrals and returns
c                      the coulomb matrix
c***description
c***references         
c
c***routines called    
c
c***end prologue       directj.f
c
      implicit none
c     --- input variables -----
      integer nat,nbtype,nprim,lenxyz,left
      integer npf,nnprim,ndmat,nbf,nnp,nnshl
      real*8 cutexp,qtest
      logical prnt,rhotest,dok
      character*(*) calc
      logical doiijk,doiikk
c     --- input arrays (unmodified) ---
      integer ptprim(nat,nbtype),noprim(nat,nbtype)
      integer nx(lenxyz),ny(lenxyz),nz(lenxyz)
      integer nocart(nbtype),mintyp(nbtype)
      integer maxmom(nbtype)
      integer pstart(nat,nbtype)
      real*8 prtoao(npf*nbf)
      real*8 ex(nprim),c(3,nat)
      real*8 dpr(nnprim,ndmat)
      real*8 qint(nnshl),dijmax(nnshl)
      character*(*) ops
c     --- input arrays (scratch) ---
      real*8 z(left)
      integer a(*)
c     --- output arrays ---
      real*8 jmat(nnp,ndmat)
      real*8 kmat(nnp,ndmat)
c     --- output variables ---
      integer ntotal,nonzer
c     --- scratch arrays ---
      real*8 jprim(nnprim,ndmat)
      real*8 kprim(nnprim,ndmat)
c     --- local variables ---
      integer inp,iout
      integer iadtwp,wpadti
      integer symcen(4),angmom(4)
      integer dmat
      integer iatom,jatom,katom,latom
      integer itype,jtype,ltype,ktype
      integer imax,jmax,kmax,lmax
      integer nfi,nfj,nfk,nfl
      integer nprimi,nprimj,nprimk,npriml
      integer jtmax,ktmax,latmax,ltmax
      integer ishell,jshell,kshell,lshell
      integer nmax,mmax
      integer ipstrt,jpstrt,kpstrt,lpstrt
      integer ijblk,klblk
      integer nroots,nij,nkl,lenblk
      integer dij,dkl,jij,jkl,ijindx,klindx
      integer dik,dil,djk,djl,kik,kil,kjk,kjl
      integer nik,nil,njk,njl
      integer aij,bkl
      integer ar,xyza,xyzam1,xyzam3,br,xyzb,xyzbm1,xyzbm3
      integer tmptop,len,index,ab,wt,urho
      integer f00,denom,b00,b10,bp01,c00,cp00,g,h,i2
      integer top
      integer expon,aa,bb,rho,t0,t1,t2,t3,t4,t5,t6,t7,t8
      integer kl,nv
      integer aplusb
      integer ijp,klp
      integer isamax
c
      real*8 pi252
      real*8 dpijmx,dpklmx,zero,one,thresh,toosmall
      real*8 estmij,estmkl,estm
      real*8 dum1,dum2,dum3,dum4,dum5,dum6
      real*8 tim1dm,timexp,timpre,timprm,timghi,timfmj
      real*8 timptj,timtrn
c
      logical debug,ijsh,klsh,ijeqkl,timeit
      logical iksh,ilsh,jksh,jlsh,iandj,iandk,iandl,jandl,kandl
      logical opensh
      save tim1dm,timexp,timpre,timprm,timghi,timfmj
      save timptj,timtrn
c
      parameter (debug=.false.,timeit=.false.)
      parameter (zero=0.0d+00,one=1.0d+00,toosmall=1.0d-50)
      data tim1dm,timexp,timpre,timprm,timghi,timfmj/6*0.0d0/
      data timptj,timtrn /2*0.0d0/
c      
      integer nodeid,mynodeid,mdtob,nxtask,nnodes,nproc,next,i,j
      integer nxtval
      integer ierr,stderr
      include 'msgtypesf.h'
      logical ispar
      common /tcgmesa/ ispar
c
      common /io/inp,iout
c
 9000 format(5x,'# primitive integrals possible',14x,i15,
     $      /5x,'# primitive integrals computed',14x,i15,'(',f5.1,'%)')
c
c
      mynodeid=nodeid()
      nproc=nnodes()
      ierr=stderr()
      opensh=(calc.eq.'open')
c
c     ----- get 2*pi**(5/2) -----
      if (mynodeid .eq. 0) then
         call iosys('read real pi from rwf',1,pi252,0,' ')
      endif
      if (ispar) call brdcst(600+MSGDBL,pi252,mdtob(1),0)
      pi252=2.0d+00*pi252**2.5d+00
c
c     --- initialize counters and set primitive j-matrix to zero ---
      ntotal=0
      nonzer=0
      call rzero(jprim,nnprim*ndmat)
      call rzero(kprim,nnprim*ndmat)
c
c
      ishell=0
      ijblk=0
      next=nxtval(nproc)+1
      do 8000 iatom=1,nat
         symcen(1)=iatom
         do 7000 itype=1,nbtype
            nprimi=noprim(iatom,itype)
            if (nprimi.le.0) go to 7000
            nfi=nocart(itype)
            ishell=ishell+1
            angmom(1)=itype
            imax=maxmom(itype)
            ipstrt=pstart(symcen(1),angmom(1))
            jshell=0
            do 6000 jatom=1,iatom
               symcen(2)=jatom
               if (jatom.eq.iatom) then
                  jtmax=itype
                  iandj=.true.
               else
                  jtmax=nbtype
                  iandj=.false.
               end if
               do 5000 jtype=1,jtmax
                  nprimj=noprim(jatom,jtype)
                  if (nprimj.le.0) go to 5000
                  nfj=nocart(jtype)
                  jshell=jshell+1
                  angmom(2)=jtype
                  jmax=maxmom(angmom(2))
                  jpstrt=pstart(symcen(2),angmom(2))
c                 --- set shell coincidence flags ---
                  ijsh=.false.
                  if(ishell.eq.jshell) ijsh=.true.
                  ijblk=ijblk+1
c
c Here it is.  PARALLELISM.  Like, wow!
c
                  if (ijblk .ne. next) goto 5000
                  nij=nprimi*nprimj
c
c                 --- allocate some core ---
                  dij=1
                  if(timeit) then
                     call timing(dum1,dum2,dum3)
                  endif
c
c                 --- get the local blocks of the density matrices ---
                  call get1dm(dpr,nnprim,ndmat,z(dij),nprimi,nprimj,
     $                        nfi,nfj,ipstrt,jpstrt)
                  if(timeit) then
                     call timing(dum4,dum5,dum6)
                     tim1dm=tim1dm+dum4-dum1
                  endif
c                 --- find the biggest one ---
                  dpijmx=z(dij+isamax(nij*nfi*nfj*ndmat,z(dij),1)-1)
c
c
                  kshell=0
                  klblk=0
                  do 4000 katom=1,iatom
                     symcen(3)=katom
                     if (katom.eq.iatom) then
                        ktmax=itype
                        iandk=.true.
                     else
                        ktmax=nbtype
                        iandk=.false.
                     end if
                     do 3000 ktype=1,ktmax
                        nprimk=noprim(katom,ktype)
                        if (nprimk.le.0) go to 3000
                        nfk=nocart(ktype)
                        kshell=kshell+1
                        angmom(3)=ktype
                        kmax=maxmom(angmom(3))
                        kpstrt=pstart(symcen(3),angmom(3))
                        if (katom.eq.iatom.and.ktype.eq.itype) then
                           latmax=jatom
                        else
                           latmax=katom
                        end if
                        lshell=0
                        do 2000 latom=1,latmax
                           symcen(4)=latom
                           iandl=(iatom.eq.latom)
                           jandl=(jatom.eq.latom)
                           kandl=(katom.eq.latom)
                           if (katom.eq.iatom.and.ktype.eq.itype.and.
     $                          latom.eq.jatom) then
                              ltmax=jtype
                           else if (latom.eq.katom) then
                              ltmax=ktype
                           else
                              ltmax=nbtype
                           end if
                           do 1000 ltype=1,ltmax 
                              npriml=noprim(latom,ltype)
                              if (npriml.le.0) go to 1000
                              klblk=klblk+1
                              lshell=lshell+1
                              ntotal=ntotal+1
c handle case where this call is to do only (II|KK) or (II|KL)/(IJ/KK) types
                              if (doiikk) then
                                 if(iandj.and.kandl) then
                                    continue
                                 else
                                    goto 1000
                                 endif
                              else if (doiijk) then
                                 if (iandj.or.kandl) then
                                    continue
                                 else
                                    goto 1000
                                 endif
                              endif
                              nfl=nocart(ltype)
                              angmom(4)=ltype
                              lmax=maxmom(angmom(4))
                              lpstrt=pstart(symcen(4),angmom(4))
c                             --- set shell coincidence flags ---
                              klsh=(kshell.eq.lshell)
c
c                             --- coarse integral screening
c                                 form an estimate of the largest
c                                 integral contribution from this
c                                 set of shell blocks. if it is smaller
c                                 than the threshhold, skip out. 
cSHOULD WE WORRY ABOUT FACTORS OF TWO HERE
                              estmij=qint(ijblk)*qint(klblk)
     $                               *dijmax(klblk)
                              estmkl=qint(klblk)*qint(ijblk)
     $                               *dijmax(ijblk)
                              estm=max(estmij,estmkl)
                              if(estm.le.qtest) go to 1000
c                              write(ierr,*) 'doing (',
c     $                             ishell,jshell,'|',kshell,lshell,')'
c
c                             --- we compute the contribution.
                              nonzer=nonzer+1
                              ijeqkl=.false.
                              if(ishell.eq.kshell.and.jshell.eq.lshell)
     $                           ijeqkl=.true.
c
                              nkl=nprimk*npriml
c
c
               if(debug.and.(mynodeid.eq.0)) then
                  write(iout,*) 'ish,jsh,ksh,lsh',ishell,jshell,
     $                           kshell,lshell
                  write(iout,*) 'symcen',
     $                           symcen(1),symcen(2),symcen(3),symcen(4)
                  write(iout,*) 'c',
     $               c(1,symcen(1)),c(2,symcen(1)),c(3,symcen(1)),
     $               c(1,symcen(2)),c(2,symcen(2)),c(3,symcen(2)),
     $               c(1,symcen(3)),c(2,symcen(3)),c(3,symcen(3)),
     $               c(1,symcen(4)),c(2,symcen(4)),c(3,symcen(4))
               endif
               nmax=imax+jmax
               mmax=kmax+lmax
               nroots=(nmax+mmax)/2+1
               lenblk=nfi*nfj*nfk*nfl
c
c              --- allocate some more core ---
               dkl=dij+nij*nfi*nfj*ndmat
               jij=dkl+nkl*nfk*nfl*ndmat
               jkl=jij+nij*nfi*nfj*ndmat
c
               aij=jkl+nkl*nfk*nfl*ndmat
               ar=aij+nij
               xyza=ar+nij
               xyzam1=xyza+3*nij
               xyzam3=xyzam1+3*nij
               ijindx=wpadti(xyzam3+3*nij)
c
               bkl=iadtwp(ijindx+2*nij)
               br=bkl+nkl
               xyzb=br+nkl
               xyzbm1=xyzb+3*nkl
               xyzbm3=xyzbm1+3*nkl
               klindx=wpadti(xyzbm3+3*nkl)
               t0=iadtwp(klindx+2*nkl)
               tmptop=t0
               top=t0+max(nprimi,nprimk)
c
c              --- get local kl density matrix ---
               if(timeit) then
                  call timing(dum1,dum2,dum3)
               endif
               call get1dm(dpr,nnprim,ndmat,z(dkl),nprimk,npriml,
     $                     nfk,nfl,kpstrt,lpstrt)
               if(timeit) then
                  call timing(dum4,dum5,dum6)
                  tim1dm=tim1dm+dum4-dum1
               endif
c              --- and the largest one
               dpklmx=z(dkl+isamax(nkl*nfk*nfl*ndmat,z(dkl),1)-1)
               if(debug.and.(mynodeid.eq.0)) then
                  call prntld(symcen,angmom,nprimi,nprimj,nprimk,
     $                        npriml,nfi,nfj,nfk,nfl,ipstrt,jpstrt,
     $                        kpstrt,lpstrt,ndmat,z(dij),z(dkl))
               end if
c
c              --- get primitive information ---
               if(timeit) then
                  call timing(dum1,dum2,dum3)
               endif
               call ldexp(z(aij),z(ar),ex,nprim,nij,
     $                    ptprim(symcen(1),angmom(1)),
     $                    ptprim(symcen(1),angmom(1))+nprimi-1,
     $                    ptprim(symcen(2),angmom(2)),
     $                    ptprim(symcen(2),angmom(2))+nprimj-1,
     $                    symcen(1),symcen(3),symcen(1),symcen(2),c,
     $                    z(t0),z(xyza),z(xyzam1),z(xyzam3),nat,nprimi,
     $                    a(ijindx))
               call ldexp(z(bkl),z(br),ex,nprim,nkl,
     $                    ptprim(symcen(3),angmom(3)),
     $                    ptprim(symcen(3),angmom(3))+nprimk-1,
     $                    ptprim(symcen(4),angmom(4)),
     $                    ptprim(symcen(4),angmom(4))+npriml-1,
     $                    symcen(1),symcen(3),symcen(3),symcen(4),c,
     $                    z(t0),z(xyzb),z(xyzbm1),z(xyzbm3),nat,nprimk,
     $                    a(klindx))
               if(timeit) then
                  call timing(dum4,dum5,dum6)
                  timexp=timexp+dum4-dum1
               endif
c
c              --- clear the primitive j-matrix blocks ---
               call rzero(z(jij),nij*nfi*nfj*ndmat)
               call rzero(z(jkl),nkl*nfk*nfl*ndmat)
c
c              --- work out how long the vectors can be ---
c              we need a number of matrices which are some length len.
c              before we allocate them, see how much core we
c              have and determine how big len can be to fit them all in.
c
c              prefac needs 6*len integers and 1*len working precision
c              prims needs (4+2*3)*nroots*len + scratch((13+5*nroots)*len)
c              vmkghi needs 3*(imax+1)*(jmax+1)*(mmax+1)*(lmax+1)*nroots*len
c                           +scratch(3*(nmax+1)*(mmax+1)*nroots
c                                   +3*(nmax+1)*(mmax+1)*(jmax+1)*nroots)*len
               len=nij*nkl
c
c              --- allocate even more core ---
               index=wpadti(tmptop)
               expon=iadtwp(index+6*len)
               t0=expon+len
               tmptop=t0
               top=t0+nij
c
c              --- fine integral screening. tests individual primitives.
c                  form the exponential prefactor and toss out
c                  small integrals. since this integral block is going to be
c                  multiplied by density matrix elements, the appropriate
c                  threshhold to test against is cutexp/max(dpijmx,dpklmx) 
               thresh=cutexp
               if(rhotest) then
                  thresh=max(abs(dpijmx),abs(dpklmx))
                  if(thresh.gt.toosmall) then
                     thresh=cutexp/thresh
                  else
                     thresh=one/toosmall
                  endif
               endif 
               kl=0
               if(timeit) then
                  call timing(dum1,dum2,dum3)
               endif
               call prefac(z(ar),z(aij),nij,z(br),z(bkl),nkl,
     $                     z(expon),a(index),a(ijindx),
     $                     a(klindx),kl,len,nv,thresh,
     $                     pi252,z(t0))
               if(timeit) then
                  call timing(dum4,dum5,dum6)
                  timpre=timpre+dum4-dum1
               endif
               if(kl.ne.nkl) then
                  call plnkerr('kl.ne.nkl after prefac',501)
               endif
               if(debug.and.(mynodeid.eq.0)) then
                  write(iout,*) 'nv of total',nv,nij*nkl
               endif
c
c              --- if no integrals pass the prescreening, skip out ---
               if(nv.ne.0) then
c                 --- allocate more core ---
                  f00=tmptop
                  b00=f00+nv*nroots
                  b10=b00+nv*nroots
                  bp01=b10+nv*nroots
                  c00=bp01+nv*nroots
                  cp00=c00+3*nv*nroots
c
                  ab=cp00+3*nv*nroots
                  aplusb=ab+nv
                  urho=aplusb+nv
                  wt=urho+nv*nroots
                  denom=wt+nv*nroots
                  aa=denom+nv*nroots
                  bb=aa+nv
                  rho=bb+nv
                  t1=rho+nv
                  t2=t1+nv
                  t3=t2+nv
                  t4=t3+nv
                  t5=t4+nv*nroots
                  t6=t5+nv*nroots
                  t7=t6+nv
                  t8=t7+nv
c
c                 --- form auxiliary arrays from primitive info ---
                  if(timeit) then
                     call timing(dum1,dum2,dum3)
                  endif
                  call prims(z(aij),z(xyza),z(xyzam1),z(xyzam3),
     $                       nij,z(bkl),z(xyzb),z(xyzbm1),
     $                       z(xyzbm3),nkl,z(f00),z(b00),
     $                       z(b10),z(bp01),z(c00),z(cp00),
     $                       z(ab),z(aplusb),z(urho),z(wt),z(denom),
     $                       z(aa),z(bb),z(rho),z(expon),
     $                       z(t1),z(t2),z(t3),z(t4),z(t5),
     $                       z(t6),z(t7),z(t8),nv,nv*nroots,nmax,
     $                       mmax,nroots,a(index),len)
                  if(timeit) then
                     call timing(dum4,dum5,dum6)
                     timprm=timprm+dum4-dum1
                  endif
c
c                 --- allocate last bit of core ---
c                     note that g and i2 are overlapped
                  g=ab
                  i2=g
                  h=i2+
     $                 3*(imax+1)*(jmax+1)*(mmax+1)*(lmax+1)*nv*nroots
                  top=h+3*(nmax+1)*(jmax+1)*(mmax+1)*nv*nroots
                  if(top.gt.left) then
                     if (mynodeid.eq.0)
     $                    write(iout,*) 'top,left',top,left
                     call plnkerr('need more memory in directj',502)
                  endif
c
c                 --- form two-dimensional integrals ---
                  if(timeit) then
                     call timing(dum1,dum2,dum3)
                  endif
                  call vmkghi(z(g),z(h),z(i2),z(f00),z(b00),z(b10),
     $                        z(bp01),z(c00),z(cp00),nmax,mmax,
     $                        imax,jmax,kmax,lmax,c,nat,symcen,
     $                        nv*nroots)
                  if(timeit) then
                     call timing(dum4,dum5,dum6)
                     timghi=timghi+dum4-dum1
                  endif
c
c                 --- form coulomb matrices ---
                  do 100 dmat=1,ndmat
                     ijp=(dmat-1)*nij*nfi*nfj
                     klp=(dmat-1)*nkl*nfk*nfl
                     if(timeit) then
                        call timing(dum1,dum2,dum3)
                     endif
                     call vfmj(z(i2),z(dij+ijp),z(dkl+klp),
     $                         z(jij+ijp),z(jkl+klp),
     $                         a(index),nij,nkl,
     $                         nfi,nfj,nfk,nfl,
     $                         ijsh,klsh,
     $                         nv*nroots,angmom,imax,jmax,mmax,lmax,
     $                         nroots,nx,ny,nz,lenxyz,mintyp,nocart,
     $                         nbtype,nv,len)
                     if(timeit) then
                        call timing(dum4,dum5,dum6)
                        timfmj=timfmj+dum4-dum1
                     endif
  100             continue
c
c                 --- scatter back into full primitive array
                  if(timeit) then
                     call timing(dum1,dum2,dum3)
                  endif
                  call putj(jprim,nnprim,ndmat,z(jij),nprimi,nprimj,
     $                      nfi,nfj,ipstrt,jpstrt,ijsh)
                  if (mynodeid .eq. 0) then
                     if(debug) then
                        write(iout,*) 'after ijputj,jprim:'
                        call print(jprim,nnprim,npf,iout)
                     endif
                  endif
c                 don't accumulate if klshell=ijshell (already did it above)
                  if(.not.ijeqkl) then
                     call putj(jprim,nnprim,ndmat,z(jkl),nprimk,npriml,
     $                         nfk,nfl,kpstrt,lpstrt,klsh)
                     if (mynodeid .eq. 0) then
                        if(debug) then
                           write(iout,*) 'after klputj,jprim:'
                           call print(jprim,nnprim,npf,iout)
                        endif
                     endif
                  endif
                  if(timeit) then
                     call timing(dum4,dum5,dum6)
                     timptj=timptj+dum4-dum1
                  endif
c now do exchange matrices if asked to
                  if (dok) then
                     iksh=(ishell.eq.kshell)
                     ilsh=(ishell.eq.lshell)
                     jksh=(jshell.eq.kshell)
                     jlsh=(jshell.eq.lshell)
                     nik=nprimi*nprimk
                     nil=nprimi*npriml
                     njk=nprimj*nprimk
                     njl=nprimj*npriml
c
c                 --- form exchange matrices
c                     funky logic to handle alpha/beta K instead of open/closed
                     if (.not.opensh) then
c at this point we shouldn't need h's storage anymore, so:
                        top=h
                        dik=top
                        dil=dik+nik*nfi*nfk
                        djk=dil+nil*nfi*nfl
                        djl=djk+njk*nfj*nfk
                        kik=djl+njl*nfj*nfl
                        kil=kik+nik*nfi*nfk
                        kjk=kil+nil*nfi*nfl
                        kjl=kjk+njk*nfj*nfk
                        top=kjl+njl*nfj*nfl
                        if(top.gt.left) then
                           if (mynodeid.eq.0)
     $                          write(iout,*) 'top,left',top,left
                           call plnkerr('need more memory in directj',
     $                          503)
                        endif
                        call rzero(z(kik),nprimi*nprimk*nfi*nfk*ndmat)
                        call rzero(z(kil),nprimi*npriml*nfi*nfl*ndmat)
                        call rzero(z(kjk),nprimj*nprimk*nfj*nfk*ndmat)
                        call rzero(z(kjl),nprimj*npriml*nfj*nfl*ndmat)
                        call get1dm(dpr,nnprim,ndmat,z(dik),nprimi,
     $                       nprimk,nfi,nfk,ipstrt,kpstrt)
                        call get1dm(dpr,nnprim,ndmat,z(dil),nprimi,
     $                       nprimk,nfi,nfl,ipstrt,lpstrt)
                        call get1dm(dpr,nnprim,ndmat,z(djk),nprimj,
     $                       nprimk,nfj,nfk,jpstrt,kpstrt)
                        call get1dm(dpr,nnprim,ndmat,z(djl),nprimj,
     $                       npriml,nfj,nfl,jpstrt,lpstrt)
                        call vfmk(z(i2),z(dik),z(dil),z(djk),z(djl),
     $                       z(kik),z(kil),z(kjk),z(kjl),a(index),
     $                       nprimi,nprimj,nprimk,npriml,nfi,nfj,nfk,
     $                       nfl,ijsh,klsh,iksh,ilsh,jksh,jlsh,
     $                       nv*nroots,angmom,imax,jmax,mmax,lmax,
     $                       nroots,nx,ny,nz,
     $                       lenxyz,mintyp,nocart,nbtype,nv,len)
c
c and put them into the proper places in K
c
                        call putj(kprim,nnprim,ndmat,z(kik),nprimi,
     $                       nprimk,nfi,nfk,ipstrt,kpstrt,iksh)
                        call putj(kprim,nnprim,ndmat,z(kil),nprimi,
     $                       npriml,nfi,nfl,ipstrt,lpstrt,ilsh)

c j<k is one of the cases we have to throw away --- that's the
c wrong triangle.  The other cases of wrong-triangle are already taken care of.
                        if (jshell .ge. kshell) then
                           call putj(kprim,nnprim,ndmat,z(kjk),nprimj,
     $                          nprimk,nfj,nfk,jpstrt,kpstrt,jksh)
                        endif
                        if (jshell .ge. lshell) then
                           call putj(kprim,nnprim,ndmat,z(kjl),nprimj,
     $                          npriml,nfj,nfl,jpstrt,lpstrt,jlsh)
                        endif
                     else
c form open and closed versions, add up to make alpha and beta
c first form "open" guy in f(1), then form closed guy in f(2), add 2 to 1
c and viola!

                     endif
                  endif
               endif
c
c
 1000                      continue
 2000                   continue
 3000                continue
 4000             continue
                  next=nxtval(nproc)+1

 5000          continue
 6000       continue
 7000    continue
 8000 continue
c
c     --- transform the j-matrix in the primitive matrix to
c         the contracted matrix.
      t1=1
      t2=t1+npf*npf
      top=t2+npf*nbf
c
      if (mynodeid .eq. 0) then
         if(debug) then
            write(iout,*) 'before final transformation'
            do 8105 dmat=1,ndmat
               write(iout,*) 'primitive j-matrix:',dmat
               call print(jprim(1,dmat),nnprim,npf,iout)
 8105       continue
         endif
      endif
c
      if(timeit) then
         call timing(dum1,dum2,dum3)
      endif
      do 8100 dmat=1,ndmat
         call trtosq(z(t1),jprim(1,dmat),npf,nnprim)
         call ebc(z(t2),z(t1),prtoao,npf,npf,nbf)
         call ebtc(z(t1),prtoao,z(t2),nbf,npf,nbf)
         call sqtotr(jmat(1,dmat),z(t1),nbf,nnp)
 8100 continue
      if(timeit) then
         call timing(dum4,dum5,dum6)
         timtrn=timtrn+dum4-dum1
      endif
      if (mynodeid .eq. 0) then
         if(debug) then
            do 8110 dmat=1,ndmat
               write(iout,*) 'contracted j-matrix:',dmat
               call print(jmat(1,dmat),nnp,nbf,iout)
 8110       continue
         endif
      endif
      if (dok) then
         if (mynodeid .eq. 0) then
            if (debug) then
               write(iout,*) 'before final transformation'
               do 8115 dmat=1,ndmat
                  write(iout,*) 'primitive k-matrix:',dmat
                  call print(kprim(1,dmat),nnprim,npf,iout)
 8115          continue
            endif
         endif

         do 8200 dmat=1,ndmat
            call trtosq(z(t1),kprim(1,dmat),npf,nnprim)
            call ebc(z(t2),z(t1),prtoao,npf,npf,nbf)
            call ebtc(z(t1),prtoao,z(t2),nbf,npf,nbf)
            call sqtotr(kmat(1,dmat),z(t1),nbf,nnp)
 8200    continue
         if (mynodeid .eq. 0) then
            if(debug) then
               do 8210 dmat=1,ndmat
                  write(iout,*) 'contracted k-matrix:',dmat
                  call print(jmat(1,dmat),nnp,nbf,iout)
 8210          continue
            endif
         endif
      endif
c
c
      if (mynodeid .eq. 0) then
         if(prnt) write (iout,9000) ntotal,nonzer,
     $        float(nonzer)*100.0d0/float(ntotal)
         if(timeit) then
            write(iout,*) '1dm,exp,pre,prm,ghi,fmj,ptj,trn',tim1dm,
     $           timexp,
     $           timpre,timprm,timghi,timfmj,timptj,timtrn
         endif
      endif
c
c BARRIER:
c
      next=nxtask(-nproc)
      return
      end
