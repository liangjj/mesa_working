*deck @(#)salc.f	5.1  11/6/94
      subroutine salc(atprmt,symat,bfchar,bftrac,lambda,char,gampt,
     #     gamma,bftran,ptbftr,p,coeffs,maxfnc,
     #     natoms,nbf,nop,maxmom,nirrep,lengam,nfunc,
     #     sc,lnsc,t,labels,maxlam,dump,ncart,nsymat,ns,relatm,mcu,
     $     momatm,ptsc,aords,maxsao,naords,numso,nocont,bfstrt,temp)
c
c***begin prologue     salc
c***date written       871004   (yymmdd)
c***revision date      900109   (yymmdd)
c
c     9 january 1990   rlm at lanl
c          fixing bug with the bfchar array.
c     4 october 1987   pws at lanl
c          adding 'numso' to keep track of the number of symmetry
c          orbitals of each irrep.
c
c***keywords           symmetry adapted linear combinations
c***author             saxe, paul (lanl)
c***source             @(#)salc.f	5.1   11/6/94
c
c***purpose            to form symmetry adapted linear combinations
c                      of atomic orbitals.
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       salc
c
      implicit integer (a-z)
c
      integer atprmt(natoms,nop),symat(natoms),lambda(nirrep)
      integer numso(nirrep)
      integer nocont(natoms,0:maxmom),bfstrt(natoms,0:maxmom)
      integer gampt(nirrep),nfunc(0:maxmom),labels(nbf)
      integer ptbftr(0:maxmom)
      integer nsymat(ns)
      integer relatm(mcu,ns)
      integer momatm(natoms)
      integer ptsc(0:maxmom,ns)
      integer aords(maxsao,naords)
      logical dump
      real*8 sc(lnsc)
      real*8 t(maxfnc,maxlam)
      real*8 trace
      real*8 bfchar(0:maxmom,nop),bftrac(0:maxmom,nop)
      real*8 char(nirrep,nop),gamma(lengam,nop)
      real*8 bftran(*),temp(nbf)
      real*8 p(maxfnc*natoms,maxfnc*natoms,maxlam)
      real*8 coeffs(nbf,nbf)
      real*8 rep
      integer ncart(0:maxmom)
c
      common/io/inp,iout
c
c     zero the coefficient array which will end up with the basis function
c     to salc transformation matrix.
      call rzero(coeffs,nbf*nbf)
c
c     ----- get the traces of the s,p,d... transformation matrices
      call rzero(bftrac,(maxmom+1)*nop)
      call izero(numso,nirrep)
c
      do 50 op=1,nop
         bftrac(0,op)=1.0d+00
         do 49 angmom=1,maxmom
            pt=ptbftr(angmom)+ncart(angmom)**2*(op-1)
            bftrac(angmom,op)=trace(bftran(pt),ncart(angmom))
 49      continue
 50   continue
c
c     ----- loop over atoms, procrastinating till we get to the
c     the last of a symmetry related set
c
      rds=0
      norbs=0
      do 1000 iatom=1,natoms
         do 1 op=1,nop
            if (atprmt(iatom,op).gt.iatom) go to 1000
    1    continue
c
c     ----- find the characters of the reducible representations
c
         call rzero(bfchar,(maxmom+1)*nop)
c
         do 10 jatom=1,natoms
            if (symat(jatom).eq.symat(iatom)) then
c
c     ----- this equivalent center only contributes to
c     the character of operations which map it into
c     itself.
c
               do 5 op=1,nop
                  if (atprmt(jatom,op).eq.jatom) then
                     do 4 angmom=0,momatm(iatom)
                        bfchar(angmom,op)=bfchar(angmom,op)+
     #                       bftrac(angmom,op)
    4                continue
                  end if
    5          continue
            end if
 10      continue
c
c     ----- decompose to irreducible representations -----
c
         pt=ptsc(0,symat(iatom))
         do 800 angmom=0,momatm(iatom)
            rds=rds+1
            nso=0
            nsf=ncart(angmom)*nsymat(symat(iatom))
            ptsave=pt
            if (ptsave.ne.ptsc(angmom,symat(iatom))) then
               call lnkerr('bad symmetry contraction pointer')
            end if
            do 700 irrep=1,nirrep
               l=lambda(irrep)
               rep=0.0d+00
               do 20 op=1,nop
                  rep=rep+bfchar(angmom,op)*char(irrep,op)
 20            continue
c
               nrep=(rep+0.5d+00)/nop
c
c
c     dump debug information?
               if(dump) then
                  write(iout,234) iatom,irrep,nrep,angmom
 234              format(' atom',i2,' irrep',i2,
     $                 ' number of salcs',i3,' momentum',i2)
               endif
               if (nrep.le.0) go to 700
c
c     ----- form this part of the ao-so matrix using
c     projection operators.
c
               call projct(gampt(irrep),atprmt,gamma,
     $              bftran(ptbftr(angmom)),
     #              p,sc(pt+1),nirrep,natoms,nop,
     #              lengam,nfunc(angmom),lambda(irrep),
     #              iatom,angmom,nrep,symat,nsymat(symat(iatom)),
     #              irrep,relatm(1,symat(iatom)),dump)
c
c     ----- put this block of coefficients into the main transformation
c           matrix.
c
               if (pt.gt.lnsc) then
                  call lnkerr('symmetry contraction matrices longer'//
     $                 ' than expected !?')
               end if
c
               call bftoso(nfunc(angmom),nsymat(symat(iatom)),
     $                     lambda(irrep),
     $                     nrep,sc(pt+1),coeffs(1,norbs+1),
     $                     natoms,relatm(1,symat(iatom)),
     $                     nbf,bfstrt(1,angmom),nocont(iatom,angmom),
     $                     dump)
c
               do 650 i=1,nrep*lambda(irrep)
                  aords(nso+i,rds)=irrep
  650          continue
               do 660 j=1,nrep*lambda(irrep)*nocont(iatom,angmom)
                  labels(norbs+j)=irrep
  660          continue
               norbs=norbs+nrep*lambda(irrep)*nocont(iatom,angmom)
               nso=nso+nrep*lambda(irrep)
               pt=pt+nsf*nrep*lambda(irrep)
               numso(irrep)=numso(irrep)+nrep*nocont(iatom,angmom)
c
  700       continue
c
c
            if(dump) then
               write (iout,710) iatom,angmom
  710          format(//,t5,'general symmetry contraction scheme: atom '
     $               ,i3,'.  angular momentum ',i3)
               write (iout,720) (aords(i,rds),i=1,nsf)
  720          format(/,8i10,/)
               call matout(sc(ptsave+1),nsf,nsf,nsf,nsf,iout)
            end if
c
  800    continue
 1000 continue
c
      if (rds.ne.naords) then
         call lnkerr('fault in number of ao reduction sets')
      end if
      if (norbs.ne.nbf) then
         call lnkerr('fault in number of symmetry orbitals.')
      end if
c
c     reorder the basis function to symmetry orbital transformation
c     matrix so that all functions of symmetry 1 come first, etc.
      k=0
      do 1100 irrep=1,nirrep
         do 1090 j=k+1,norbs
            if(labels(j).eq.irrep) then
               k=k+1
               tmplbl=labels(j)
               labels(j)=labels(k)
               labels(k)=tmplbl
               call vmove(temp,coeffs(1,j),nbf)
               call vmove(coeffs(1,j),coeffs(1,k),nbf)
               call vmove(coeffs(1,k),temp,nbf)
            end if
 1090    continue
 1100 continue
c
c
      return
      end
