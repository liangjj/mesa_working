*deck @(#)soint.f	5.1   11/6/94
      subroutine soint(c,ex,z,iz,cont,s,ptprim,noprim,nocont,ptcont,
     #                  nat,nprim,maxcor,ntypes,nbtype,nnp,ncont,
     #                  start,nbasis,zan,nocart,nobf,maxmom,mintyp,
     #                  nx,ny,nz,minmom,dolp,bflabl,ops,zeff)
c
c***module to form the one-electron spin-orbit integrals 
c   integrals over generally contracted gaussian basis
c   sets.
c
c matt braunstein              dec. 1991   lanl 
c
      implicit integer (a-z)
c
      character*(*) ops, bflabl(*)
      real*8 c(3,nat),ex(nprim),z(maxcor),s(nbasis*nbasis),cont(ncont)
      real*8 zan(nat),zeff(nat,ntypes)
c
      real*8 alpha
c
      integer iz(*)
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(ntypes)
      integer nobf(ntypes),maxmom(ntypes),mintyp(ntypes),minmom(ntypes)
      integer nx(*),ny(*),nz(*)
      logical dolp,logkey
      character refops*8
c
      common/io/inp,iout
c
      data alpha/7.29720e-3/
c
      write(iout,*)
      write(iout,*)'Evaluating one-electron spin-orbit integrals'
      write(iout,*)
c
c     ----- form the one-electron spin-orbit integrals -----
c
c The spin-orbit operator has the form 
c
c     V^eff_s-o = alpha^2/2  sum_{i_k,k} (l_k,i_k * s_i_k) * Z^eff_k /r^3, 
c
c   where alpha is the fine structure constant, k indexes the nuclei and
c   i_k indexes the electrons on center k. Z^eff_k is adjusted to reproduce the
c   atomic spin-orbit splitting. Therefore only one-center one-electron terms
c   are included in this operator
c
c   See, for example, W. R. Wadt, Chem. Phys. Lett. 89, 245 (1982).
c
c
      call rzero(zeff,nat*ntypes)
      call fparr(ops,'spin-orbit',zeff,nat*ntypes,' ')
c
c print zeff
c
      if(logkey(ops,'print=int=so',.false.,refops)) then
        write(iout,*)'   zeff(nat,ntypes)   '
        do 5 i=1,nat
          write(iout,*)(zeff(i,j),j=1,ntypes)
5       continue
      end if
c
c loop over lx, ly, and lz operators
c
      do 10 lxyz=1,3
c
c loop over centers
c
        do 9 iatom=1,nat
          do 8 jatom=1,iatom
c
c restrict to one-center integrals
c
            if(jatom.ne.iatom)go to 8

c
c loop over angular momentum type
c
            do 7 itype=1,nbtype
               if (noprim(iatom,itype).le.0) go to 7
               if (iatom.ne.jatom) then
                  jtypmx=nbtype
               else
                  jtypmx=itype
               end if
               do 6 jtype=1,jtypmx
                  if (noprim(jatom,jtype).le.0) go to 6
c
c restrict to integrals of the same angular momentum type
c
                  if(jtype.ne.itype) go to 6
c
                  imax=maxmom(itype)
                  jmax=maxmom(jtype)
                  nprimi=noprim(iatom,itype)
                  nprimj=noprim(jatom,jtype)
                  nconti=nocont(iatom,itype)
                  ncontj=nocont(jatom,jtype)
                  npint=nprimi*nprimj
                  lenblk=nocart(itype)*nocart(jtype)
c
c     ----- allocate core for temporary vectors, etc. -----
c
                  prmint=1
                  xyz=lenblk*npint
c
                  conint=prmint+npint*lenblk
                  tmp1=conint+nconti*ncontj*lenblk
                  len1=nconti*nprimj
                  top2=tmp1+len1
c
                  if (top2.gt.maxcor) then
                     call lnkerr('m304: soint..not enough core for so')
                  end if
c
c     ----- form the primitive spin-orbit integrals -----
c
c  lx block
c
                  if(lxyz.eq.1) then
                    call sox(nprim,nprimi,nprimj,ptprim(iatom,itype),
     #                     ptprim(jatom,jtype),ex,xyz,itype,jtype,
     #                     nocart(itype),nocart(jtype),z(prmint))
                  end if
c
c  ly block
c
                  if(lxyz.eq.2) then
                    call soy(nprim,nprimi,nprimj,ptprim(iatom,itype),
     #                     ptprim(jatom,jtype),ex,xyz,itype,jtype,
     #                     nocart(itype),nocart(jtype),z(prmint))
                  end if
c
c  lz block
c
                  if(lxyz.eq.3) then
                    call soz(nprim,nprimi,nprimj,ptprim(iatom,itype),
     #                     ptprim(jatom,jtype),ex,xyz,itype,jtype,
     #                     nocart(itype),nocart(jtype),z(prmint))
                  end if
c
c At this point multiply by constants
c
                  do 20 i=1,xyz
                    z(prmint+i-1)=(alpha*alpha/2.d0)*zeff(iatom,itype)
     $                            *z(prmint+i-1)
20                continue

c
c     ----- transform to contracted functions -----
c
                  call trans1(z(prmint),z(conint),nprimi,nprimj,nconti,
     #                        ncontj,cont(ptcont(iatom,itype)),
     #                        cont(ptcont(jatom,jtype)),z(tmp1),len1,
     #                        lenblk,minmom(itype),maxmom(itype),
     #                        minmom(jtype),maxmom(jtype),nocart)
c
c
                 
c
c     ----- transfer integrals to total array -----
c             this array is square (s)
c
                call put1so(s,z(conint),start,iatom,jatom,itype,jtype,
     #                      nconti,ncontj,nbasis,lenblk,nat,nbtype,nobf)
c
    6           continue
    7         continue
    8      continue
    9    continue
c
c     ----- print the integrals -----
c
         if(logkey(ops,'print=int=so',.false.,refops)) then
           write(iout,*)'atomic spin-orbit integrals'
           call wmat(s,nbasis,nbasis,bflabl,bflabl)
         end if
c
c     ----- put atomic spin-orbit integrals on the read-write file -----
c
c lx block
c
         if (lxyz.eq.1) then
           call iosys('write real "asox integrals" on rwf',
     #                 nbasis*nbasis,s,0,' ')
         end if
c
c ly block
c
         if (lxyz.eq.2) then
           call iosys('write real "asoy integrals" on rwf',
     #                 nbasis*nbasis,s,0,' ')
         end if
c
c lz block
c
         if (lxyz.eq.3) then
           call iosys('write real "asoz integrals" on rwf',
     #                 nbasis*nbasis,s,0,' ')
         end if
10    continue
c
      return
      end
