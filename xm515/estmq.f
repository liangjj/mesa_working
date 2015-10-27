*deck @(#)estmq.f	1.2  11/28/95
      subroutine estmq(ptprim,noprim,nbtype,nocont,ptcont,cont,
     $                 ncont,ex,c,nx,ny,nz,lenxyz,nocart,mintyp,
     $                 maxmom,z,left,nat,nprim,ops,
     $                 a,pstart,nbf,nnp,nnshl,qint)
c
c***begin prologue     estmq.f
c***date written       931215   (yymmdd)
c***revision date      4/18/95
c
c  23 march, 1995      rlm at lanl
c     revising the directj routine to generate estimates for integrals.
c     see haser and ahlrichs, j.comp.chem. 10, 104(1989).
c
c  18 december, 1993   rlm at lanl
c     revising the driver routine in m712 for direct j-matrix.
c***keywords           direct, j-matrix
c***author             martin, saxe, lengsfield
c***source             @(#)estmq.f	1.2   11/28/95
c***purpose            calculates two-electron integrals and returns
c                      the coulomb matrix
c***description
c***references         
c
c***routines called    
c
c***end prologue       estmq.f
c
      implicit none
c     --- input variables -----
      integer nat,nbtype,nprim,lenxyz,left,ncont
      integer nbf,nnp,nnshl
c     --- input arrays (unmodified) ---
      integer ptprim(nat,nbtype),noprim(nat,nbtype)
      integer ptcont(nat,nbtype),nocont(nat,nbtype)
      integer nx(lenxyz),ny(lenxyz),nz(lenxyz)
      integer nocart(nbtype),mintyp(nbtype)
      integer maxmom(nbtype)
      integer pstart(nat,nbtype)
      real*8 cont(ncont)
      real*8 ex(nprim),c(3,nat)
      character*(*) ops
c     --- input arrays (scratch) ---
      real*8 z(left)
      integer a(*)
c     --- output arrays ---
      real*8 qint(nnshl)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer iadtwp,wpadti
      integer symcen(4),angmom(4)
      integer iatom,jatom,katom,latom
      integer itype,jtype,ltype,ktype
      integer imax,jmax,kmax,lmax
      integer nfi,nfj,nfk,nfl
      integer nprimi,nprimj,nprimk,npriml
      integer nconti,ncontj,ncontk,ncontl
      integer jtmax,ktmax,latmax,ltmax
      integer ishell,jshell,kshell,lshell
      integer nmax,mmax
      integer ipstrt,jpstrt,kpstrt,lpstrt
      integer nroots,nij,nkl,lenblk
      integer ijindx,klindx
      integer aij,bkl
      integer ar,xyza,xyzam1,xyzam3,br,xyzb,xyzbm1,xyzbm3
      integer tmptop,len,index,ab,wt,urho
      integer f00,denom,b00,b10,bp01,c00,cp00,g,h,i2
      integer top
      integer expon,aa,bb,rho,t0,t1,t2,t3,t4,t5,t6,t7,t8
      integer kl,nv
      integer aplusb
      integer isamax,i
      integer half,test,temp,temp1,minkl,maxkl,conint,ncint,ijshl
      integer ntotal,nonzer
c
      real*8 pi252
      real*8 zero,one,toosmall
c
      logical debug,ijsh,klsh,ijeqkl
c
      parameter (debug=.false.)
      parameter (zero=0.0d+00,one=1.0d+00,toosmall=1.0d-50)
c
      integer nodeid,mdtob,mynodeid
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
c
c     ----- get 2*pi**(5/2) -----
      if (mynodeid.eq.0) then
         call iosys('read real pi from rwf',1,pi252,0,' ')
      endif
      if (ispar) call brdcst(700+MSGDBL,pi252,mdtob(1),0)
      pi252=2.0d+00*pi252**2.5d+00
c
c     --- initialize counters.
      ntotal=0
      nonzer=0
c
c     --- we compute only a subblock of two-electron integrals. those of the
c         form (mu,nu|mu,nu).
      ijshl=0
      ishell=0
      do 8000 iatom=1,nat
         symcen(1)=iatom
         do 7000 itype=1,nbtype
            nprimi=noprim(iatom,itype)
            if (nprimi.le.0) go to 7000
            nfi=nocart(itype)
            nconti=nocont(iatom,itype)
            ishell=ishell+1
            angmom(1)=itype
            imax=maxmom(itype)
            ipstrt=pstart(symcen(1),angmom(1))
            jshell=0
            do 6000 jatom=1,iatom
               symcen(2)=jatom
               if (jatom.eq.iatom) then
                  jtmax=itype
               else
                  jtmax=nbtype
               end if
               do 5000 jtype=1,jtmax
                  nprimj=noprim(jatom,jtype)
                  if (nprimj.le.0) go to 5000
                  nfj=nocart(jtype)
                  ncontj=nocont(jatom,jtype)
                  jshell=jshell+1
                  angmom(2)=jtype
                  jmax=maxmom(angmom(2))
                  jpstrt=pstart(symcen(2),angmom(2))
c                 --- set shell coincidence flags ---
                  ijsh=.false.
                  if(ishell.eq.jshell) ijsh=.true.
                  nij=nprimi*nprimj
                  ijshl=ijshl+1
c
                  kshell=0
                  katom=iatom
                     symcen(3)=katom
                     if (katom.eq.iatom) then
                        ktmax=itype
                     else
                        ktmax=nbtype
                     end if
                     ktype=itype
                        nprimk=noprim(katom,ktype)
                        if (nprimk.le.0) go to 3000
                        nfk=nocart(ktype)
                        ncontk=nocont(katom,ktype)
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
                        latom=jatom
                           symcen(4)=latom
                           if (katom.eq.iatom.and.ktype.eq.itype.and.
     $                          latom.eq.jatom) then
                              ltmax=jtype
                           else if (latom.eq.katom) then
                              ltmax=ktype
                           else
                              ltmax=nbtype
                           end if
                           ltype=jtype
                              npriml=noprim(latom,ltype)
                              if (npriml.le.0) go to 1000
                              nfl=nocart(ltype)
                              ncontl=nocont(latom,ltype)
                              lshell=lshell+1
                              angmom(4)=ltype
                              lmax=maxmom(angmom(4))
                              lpstrt=pstart(symcen(4),angmom(4))
c                             --- set shell coincidence flags ---
                              klsh=.false.
                              if(kshell.eq.lshell) klsh=.true.
                              ijeqkl=.false.
                              if(ishell.eq.kshell.and.jshell.eq.lshell)
     $                           ijeqkl=.true.
c
                              nkl=nprimk*npriml
c
c
                              
               if(debug .and. (mynodeid .eq. 0)) then
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
               ntotal=ntotal+nij*nkl*lenblk
               ncint=nconti*ncontj*ncontk*ncontl
c
c              --- allocate some core ---
c              aij=1
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
               half=iadtwp(klindx+2*nkl)
               test=wpadti(half+nconti*ncontj*nkl*lenblk)
               t0=iadtwp(test+nij*nkl)
               tmptop=t0
               top=t0+max(nprimi,nprimk)
c
c              --- get primitive information ---
               call rzero(z(half),nconti*ncontj*nkl*lenblk)
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
c              --- form the exponential prefactor and toss out
c                  small integrals. the threshhold value is set negative
c                  so that no integrals should be tossed out.
               kl=0
               call prefac(z(ar),z(aij),nij,z(br),z(bkl),nkl,
     $                     z(expon),a(index),a(ijindx),
     $                     a(klindx),kl,len,nv,-1.0d0,
     $                     pi252,z(t0))
               if(kl.ne.nkl) then
                  call plnkerr('kl.ne.nkl after prefac',700)
               endif
               if(debug .and. (mynodeid.eq.0)) then
                  write(iout,*) 'nv of total',nv,nij*nkl
               endif
c
c              --- if no integrals pass the prescreening, skip out ---
c                  except for this estimation purpose, be sure we don't
c                  throw any out.
               if(nv.ne.len) then
                  call plnkerr('estmq: integrals thrown away',701)
               endif
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
                  nonzer=nonzer+nv*lenblk
                  call prims(z(aij),z(xyza),z(xyzam1),z(xyzam3),
     $                       nij,z(bkl),z(xyzb),z(xyzbm1),
     $                       z(xyzbm3),nkl,z(f00),z(b00),
     $                       z(b10),z(bp01),z(c00),z(cp00),
     $                       z(ab),z(aplusb),z(urho),z(wt),z(denom),
     $                       z(aa),z(bb),z(rho),z(expon),
     $                       z(t1),z(t2),z(t3),z(t4),z(t5),
     $                       z(t6),z(t7),z(t8),nv,nv*nroots,nmax,
     $                       mmax,nroots,a(index),len)
c
c                 --- allocate last bit of core ---
c                     note that g and i2 are overlapped
                  g=ab
                  i2=g
                  h=i2+
     $                 3*(imax+1)*(jmax+1)*(mmax+1)*(lmax+1)*nv*nroots
                  temp=h
                  temp1=h+len
                  top=h+3*(nmax+1)*(jmax+1)*(mmax+1)*nv*nroots
                  top=max(top,temp1+len)
                  if(top.gt.left) then
                     if (mynodeid .eq.0)
     $                    write(iout,*) 'top,left',top,left
                     call plnkerr('need more memory in directj',702)
                  endif
c
c                 --- form two-dimensional integrals ---
                  call vmkghi(z(g),z(h),z(i2),z(f00),z(b00),z(b10),
     $                        z(bp01),z(c00),z(cp00),nmax,mmax,
     $                        imax,jmax,kmax,lmax,c,nat,symcen,
     $                        nv*nroots)
c
c                 --- form primitive integrals ---
c                     must also trick the test array since we're using
c                     the new version of prims.
                  minkl=1
                  maxkl=nkl
                  do 110 i=1,nij*nkl
                     a(test+i-1)=i
  110             continue
                  call vfmint(z(i2),z(half),z(temp),lenblk,nv*nroots,
     $                        angmom,imax,jmax,mmax,lmax,0,nroots,
     $                        nx,ny,nz,lenxyz,mintyp,nocart,nbtype,
     $                        minkl,len,a(test),nv,
     $                        cont(ptcont(symcen(1),angmom(1))),
     $                        cont(ptcont(symcen(2),angmom(2))),
     $                        nprimi,nconti,nprimj,ncontj,nkl,nkl,
     $                        z(temp1),nij,len)
c
c                 --- contract the primitive integrals ---
                  conint=half+nconti*ncontj*nkl*lenblk
                  temp1=conint+ncint*lenblk
                  top=temp1+ncontk*npriml
                  call trans(z(half),z(conint),nprimk,npriml,
     $                       nconti*ncontj*lenblk,ncontk,ncontl,
     $                       z(temp1),left-temp1+1,
     $                       cont(ptcont(symcen(3),angmom(3))),
     $                       cont(ptcont(symcen(4),angmom(4))))
               endif
c
c              --- find the largest integral in the block, and store
c                  the square root in the output array.
               i=isamax(ncint*lenblk,z(conint),1)
               qint(ijshl)=sqrt(z(conint+i-1))
c
c
 1000                      continue
 3000                continue
 5000          continue
 6000       continue
 7000    continue
 8000 continue
c
c
      return
      end
