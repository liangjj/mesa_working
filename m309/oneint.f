*deck @(#)oneint.f	5.3 11/28/95
      subroutine oneint(c,ex,z,iz,cont,s,ptprim,noprim,nocont,ptcont,
     #                  nat,nprim,maxcor,ntypes,nbtype,nnp,ncont,
     #                  start,nbasis,zan,nocart,nobf,maxmom,mintyp,
     #                  nx,ny,nz,minmom,dis,dolp,bflabl,atnam,ops)
c
c***module to form the overlap, kinetic-energy and potential-energy
c   one-electron integrals over generally contracted gaussian basis
c   sets.
c
c paul saxe               23 july 1984                  lanl
c
      implicit integer (a-z)
c
      character*(*) ops,bflabl(*),atnam(nat)
      real*8 c(3,nat),ex(nprim),z(maxcor),s(nnp),cont(ncont)
      real*8 zan(nat)
      real*8 h(21),wt(21)
      real*8 dis(nat,nat)
      real*8 eh,ehe,ep,t,onsite,rmax,fpkey,zero,one,delta
      real*8 tp,td,tpd,scale
      integer iz(*)
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(ntypes)
      integer nobf(ntypes),maxmom(ntypes),mintyp(ntypes),minmom(ntypes)
      integer nx(*),ny(*),nz(*)
      logical dolp,logkey
      logical drop
      character refops*8
      character card*80
      logical positn
c
      common/io/inp,iout
      parameter (zero=0.0d+00,one=1.0d+00)

c
      data h   /0.0d+00
     2,         -.707106781186548d+00,  0.707106781186548d+00
     3,         -1.22474487139159d+00,0.0d+00,1.22474487139159d+00
     4,         -1.65068012388578d+00, -0.524647623275290d+00
     4,          0.524647623275290d+00, 1.65068012388578d+00
     5,   -2.02018287045609d+00,-0.958572464613819d+00,0.0d+00
     5,          0.958572464613819d+00, 2.02018287045609d+00
     6,         -2.350604973674d+00  , -1.335849074014d+00
     6,         -0.436077411928d+00  ,  0.436077411928d+00
     6,          1.335849074014d+00  ,  2.350604973674d+00/
      data wt  /1.77245385090552d+00
     2,         0.8862269254528d+00  ,  0.8862269254528d+00
     3,         0.2954089751509d+00  ,  1.181635900604d+00
     3,         0.2954089751509d+00
     4,         8.131283544725d-02   ,  8.049140900055d-01
     4,         8.049140900055d-01   ,  8.131283544725d-02
     5,         1.995324205905d-02   ,  3.936193231522d-01
     5,         9.453087204829d-01   ,  3.936193231522d-01
     5,         1.995324205905d-02
     6,         4.530009905509d-03   ,  1.570673203229d-01
     6,         7.246295952244d-01   ,  7.246295952244d-01
     6,         1.570673203229d-01   ,  4.530009905509d-03/
      data mxpts /21/
      save h,wt,mxpts
c
c
   10 format(5x,'overlap integrals:')
   20 format(5x,'kinetic energy integrals:')
   30 format(5x,'potential energy integrals:')
   40 format(5x,'effective core potential integrals:')
c
c     are we to drop function components?
      drop=logkey(ops,'drop',.false.,' ')
c
c     get the on-site and kinetic energy for ppp models.
      if(logkey(ops,'ppp',.false.,' ')) then
         t=fpkey(ops,'ppp=t',0.0d+00,' ')
         eh=fpkey(ops,'ppp=eh',0.0d+00,' ')
         if(eh.eq.0.0d+00) eh=fpkey(ops,'ppp=e',0.0d+00,' ')
         ehe=fpkey(ops,'ppp=ehe',0.0d+00,' ')
         ep=fpkey(ops,'ppp=ep',0.0d+00,' ')
         rmax=fpkey(ops,'ppp=rmax',0.0d+00,' ')
      endif
c     perhaps we are doing the two-site model.
      if(logkey(ops,'twosit',.false.,' ')) then
         tp=fpkey(ops,'twosit=tp',0.0d+00,' ')
         td=fpkey(ops,'twosit=td',0.0d+00,' ')
         tpd=fpkey(ops,'twosit=tpd',0.0d+00,' ')
    
         delta=fpkey(ops,'twosit=delta',0.0d+00,' ')
         rmax=fpkey(ops,'twosit=rmax',0.0d+00,' ')
      endif
c
c     compute the distance matrix
      call dismat(nat,c,dis,one)
c
c     ----- form the overlap integrals -----
      call rzero(s,nnp)
      do 9 iatom=1,nat
         do 8 jatom=1,iatom
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
                  a=1
                  b=a+npint
                  alpha=b+npint
                  ainv=alpha+2*npint
                  xyza=ainv+npint
                  expon=xyza+3*npint
                  prmint=1
                  xyz=max(expon+npint,prmint+lenblk*npint)
                  top1=xyz+3*npint*(imax+1)*(jmax+1)
c
                  conint=prmint+npint*lenblk
                  tmp1=conint+nconti*ncontj*lenblk
                  len1=nconti*nprimj
                  top2=tmp1+len1
c
                  if (top1.gt.maxcor.or.top2.gt.maxcor) then
                     call lnkerr('m302: oneint..not enough core for s')
                  end if
c
c     ----- form the two-dimensional integrals -----
c
                  call sints(iatom,jatom,imax,jmax,nprimi,nprimj,
     #                       ptprim(iatom,itype),ptprim(jatom,jtype),
     #                       ex,z(alpha),c,z(ainv),z(xyza),
     #                       z(expon),z(xyz),npint,nprim,nat,
     #                       nbtype,h,wt,mxpts,z(a),z(b))
c
c     ----- form the primitive integrals -----
c
                  call fmonel(z(prmint),z(xyz),npint,lenblk,imax,jmax,
     #                   mintyp(itype),mintyp(itype)+nocart(itype)-1,
     #                   mintyp(jtype),mintyp(jtype)+nocart(jtype)-1,
     #                   nx,ny,nz)
c
c     ----- transform to contracted functions -----
c
                  call trans1(z(prmint),z(conint),nprimi,nprimj,nconti,
     #                        ncontj,cont(ptcont(iatom,itype)),
     #                        cont(ptcont(jatom,jtype)),z(tmp1),len1,
     #                        lenblk,minmom(itype),maxmom(itype),
     #                        minmom(jtype),maxmom(jtype),nocart)
c
c     ----- transfer integrals to total array -----
c
                  call put1el(s,z(conint),start,iatom,jatom,itype,jtype,
     #                        nconti,ncontj,nnp,lenblk,nat,nbtype,nobf)
c
    6          continue
    7       continue
    8    continue
    9 continue
c
c     ----- print the integrals -----
c
      if(logkey(ops,'print=int=s',.false.,refops)) then
         write(iout,10)
         call trtosq(z(prmint),s,nbasis,nnp)
         call wlmat(z(prmint),nbasis,nbasis,bflabl,bflabl)
      end if
c
c     ----- put overlap integrals on the read-write file -----
c
      call iosys('write real "overlap integrals" on rwf',nnp,s,0,' ')
c
c     ----- form the kinetic-energy integrals -----
c
      call rzero(s,nnp)
      do 19 iatom=1,nat
         do 18 jatom=1,iatom
            do 17 itype=1,nbtype
               if (noprim(iatom,itype).le.0) go to 17
               if (iatom.ne.jatom) then
                  jtypmx=nbtype
               else
                  jtypmx=itype
               end if
               do 16 jtype=1,jtypmx
                  if (noprim(jatom,jtype).le.0) go to 16
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
                  a=1
                  b=a+npint
                  bp2=b+npint
                  bm2=bp2+npint
                  alpha=bm2+npint
                  ainv=alpha+2*npint
                  xyza=ainv+npint
                  expon=xyza+3*npint
                  prmint=1
                  xyz=max(expon+npint,prmint+lenblk*npint)
                  t1=xyz+3*npint*(imax+1)*(jmax+1)*2
                  top1=t1+npint
c
                  conint=prmint+npint*lenblk
                  tmp1=conint+nconti*ncontj*lenblk
                  len1=nconti*nprimj
                  top2=tmp1+len1
c
                  if (top1.gt.maxcor.or.top2.gt.maxcor) then
                     call lnkerr('m302: oneint..not enough core for t')
                  end if
                  call rzero(z(conint),nconti*ncontj*lenblk)
                  if(iatom.ne.jatom) then
                     if(dis(iatom,jatom).le.rmax) then
                        if(logkey(ops,'ppp',.false.,' ')) then
                           ij=0
                           do 15 i=1,nocart(itype)
                              do 14 j=1,nocart(jtype)
                                ij=ij+1
                                z(conint+ij-1)=t
   14                         continue
   15                      continue
                        end if
                        if(logkey(ops,'twosit',.false.,' ')) then
                           z(conint)=td
                           z(conint+1)=tpd
                           z(conint+2)=tpd
                           z(conint+3)=tp
                        endif
                     end if
                  end if
c
c     ----- transfer integrals to total array -----
c
                  call put1el(s,z(conint),start,iatom,jatom,itype,jtype,
     #                        nconti,ncontj,nnp,lenblk,nat,nbtype,nobf)
c
   16          continue
   17       continue
   18    continue
   19 continue
c
c     read back in overlap matrix and fix up sign of t.
      call iosys('read real "overlap integrals" from rwf',nnp,
     $            z(a),0,' ')
      do 21 i=1,nnp
         if(s(i).ne.zero) then
            if(z(i).gt.zero) then
            else if(z(i).eq.zero) then
               s(i)=zero
            else
               s(i)=-s(i)
            end if
         end if
   21 continue
c
c     for a temporary fix for the cuo problem, read in a phasing matrix
c     and use it to fix the signs of t.
      if(positn('$t-matrix',card,inp)) then
         write(iout,*) 'reading t-matrix from cards'
c        we are now positioned to read it.
         read(inp,*) (s(i),i=1,nnp)
         write(iout,*) 't-matrix',(s(i),i=1,nnp)
      endif
      
c
c     ----- print the integrals -----
c
      if(logkey(ops,'print=int=t',.false.,refops)) then
         write(iout,20)
         call trtosq(z(prmint),s,nbasis,nnp)
         call wlmat(z(prmint),nbasis,nbasis,bflabl,bflabl)
      end if
c
c     ----- put kinetic-energy integrals on the rwf -----
c
      call iosys('write real "kinetic integrals" on rwf',nnp,s,0,' ')
c
c     prepare and write a unit s-matrix.
      call rzero(s,nnp)
      do 23 i=1,nbasis
         ij=i*(i+1)/2
         s(ij)=one
   23 continue
c
c     ----- print the integrals -----
c
      if(logkey(ops,'print=int=s',.false.,refops)) then
         write(iout,10)
         call trtosq(z(prmint),s,nbasis,nnp)
         call wlmat(z(prmint),nbasis,nbasis,bflabl,bflabl)
      end if
c
c     ----- put overlap and kinetic-energy integrals on the rwf -----
c
      call iosys('write real "overlap integrals" on rwf',nnp,s,0,' ')
c
c     ----- form the potential-energy integrals -----
c
      call rzero(s,nnp)
      do 29 iatom=1,nat
         do 28 jatom=1,iatom
            do 27 itype=1,nbtype
               if (noprim(iatom,itype).le.0) go to 27
               if (iatom.ne.jatom) then
                  jtypmx=nbtype
               else
                  jtypmx=itype
               end if
               do 26 jtype=1,jtypmx
                  if (noprim(jatom,jtype).le.0) go to 26
c
                  imax=maxmom(itype)
                  jmax=maxmom(jtype)
                  nroots=(imax+jmax)/2+1
                  nprimi=noprim(iatom,itype)
                  nprimj=noprim(jatom,jtype)
                  nconti=nocont(iatom,itype)
                  ncontj=nocont(jatom,jtype)
                  npint=nprimi*nprimj
                  lenblk=nocart(itype)*nocart(jtype)
c
c     ----- allocate core for temporary vectors, etc. -----
c
                  a=1
                  b=a+npint
                  aiaj=b+npint
                  xyz0=aiaj+npint
                  xyz1=xyz0+npint*3
                  rysrt=xyz1+npint*3
                  ryswt=rysrt+npint*nroots
                  alpha=ryswt+npint*nroots
                  ainv=alpha+2*npint
                  xyza=ainv+npint
                  expon=xyza+3*npint
                  prmint=expon+npint
                  xyz=max(expon+npint,prmint+lenblk*npint)
                  t1=xyz+3*npint*(imax+1)*(jmax+1)*2
                  top1=t1+npint
c
                  conint=prmint+npint*lenblk
                  tmp1=conint+nconti*ncontj*lenblk
                  len1=nconti*nprimj
                  top2=tmp1+len1
c
                  if (top1.gt.maxcor.or.top2.gt.maxcor) then
                     call lnkerr('m302: oneint..not enough core for v')
                  end if
c     semiempirical fix ... load the diagonal energies.
                  call rzero(z(conint),nconti*ncontj*lenblk)
                  if(logkey(ops,'ppp',.false.,' ')) then
                     if(iatom.eq.jatom) then
                        if(atnam(iatom)(1:2).eq.'he') then
                           onsite=ehe
                        else if(atnam(iatom)(1:1).eq.'h') then
                           onsite=eh
                        else if(atnam(iatom)(1:1).eq.'c') then
                           onsite=ep
                        end if
                        do 25 i=1,nocart(itype)
c                          scale the p-z components to core-like energies.
                           if(i.eq.3) then
                              scale=10.0
                           else 
                              scale=1.0
                           endif
                           ij=(i-1)*nocart(itype)+i
                           z(conint+ij-1)=scale*onsite
  25                    continue
                     end if
                  endif
                  if(logkey(ops,'twosit',.false.,' ')) then
                     if(iatom.eq.jatom) then
                        z(conint+3)=delta
                     endif
                  endif
c
c
c     ----- transfer integrals to total array -----
c
                  call put1el(s,z(conint),start,iatom,jatom,itype,jtype,
     #                        nconti,ncontj,nnp,lenblk,nat,nbtype,nobf)
c
   26          continue
   27       continue
   28    continue
   29 continue
      if(positn('$e-matrix',card,inp)) then
         write(iout,*) 'reading e-matrix from cards'
c        we are now positioned to read it.
         read(inp,*) (s(i),i=1,nnp)
         write(iout,*) 'e-matrix',(s(i),i=1,nnp)
      endif
c
c     ----- print the integrals -----
c
      if(logkey(ops,'print=int=v',.false.,refops)) then
         write(iout,30)
         call trtosq(z(prmint),s,nbasis,nnp)
         call wlmat(z(prmint),nbasis,nbasis,bflabl,bflabl)
      end if
c
c     ----- put nuclear-potential integrals on the read-write file -----
c
      call iosys('write real "potential integrals" on rwf',nnp,s,0,' ')
c
c     are there any ecp centers?
      call rzero(s,nnp)
      if(dolp) then
c
c        ----- form the effective core potential integrals -----
c
         call ldata
         call ztab
         do 39 iatom=1,nat
            do 38 jatom=1,iatom
               do 37 itype=1,nbtype
                  if (noprim(iatom,itype).le.0) go to 37
                  if (iatom.ne.jatom) then
                     jtypmx=nbtype
                  else
                     jtypmx=itype
                  end if
                  do 36 jtype=1,jtypmx
                     if (noprim(jatom,jtype).le.0) go to 36
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
c        ----- allocate core for temporary vectors, etc. -----
c
                     prmint=1
                     ntpse=prmint+npint*lenblk
                     nlp=ntpse+10
                     zlp=nlp+100
                     clp=zlp+100
                     qq=clp+100
                     top1=qq+343*npint
c
                     conint=prmint+npint*lenblk
                     tmp1=conint+nconti*ncontj*lenblk
                     len1=nconti*nprimj
                     top2=tmp1+len1
c
                     if (top1.gt.maxcor.or.top2.gt.maxcor) then
                        call lnkerr('m302: oneint..not enough core'
     $                            //' for lp')
                     end if
c
c        ----- form the primitive integrals -----
c
                     call lpints(iatom,jatom,imax,jmax,nprimi,nprimj,
     #                       nprim,ncont,nat,nbtype,ntypes,npint,lenblk,
     #                       ptprim(iatom,itype),ptprim(jatom,jtype),
     #                       c,ex,noprim,ptprim,ptcont,cont,maxmom,
     #                       mintyp(itype),
     #                       mintyp(itype)+nocart(itype)-1,
     #                       mintyp(jtype),
     #                       mintyp(jtype)+nocart(jtype)-1,
     #                       nx,ny,nz,z(prmint),z(ntpse),z(nlp),z(zlp),
     #                       z(clp),z(qq))
c
c        ----- transform to contracted functions -----
c
                     call trans1(z(prmint),z(conint),nprimi,nprimj,
     #                        nconti,ncontj,cont(ptcont(iatom,itype)),
     #                        cont(ptcont(jatom,jtype)),z(tmp1),len1,
     #                        lenblk,minmom(itype),maxmom(itype),
     #                        minmom(jtype),maxmom(jtype),nocart)
c
c        ----- transfer integrals to total array -----
c
                     call put1el(s,z(conint),start,iatom,jatom,itype,
     #                        jtype,nconti,ncontj,nnp,lenblk,nat,nbtype,
     #                        nobf)
c
   36             continue
   37          continue
   38       continue
   39    continue
c
c        ----- print the integrals -----
c
         if(logkey(ops,'print=int=ecp',.false.,refops)) then
            write(iout,40)
            call trtosq(z(prmint),s,nbasis,nnp)
            call wlmat(z(prmint),nbasis,nbasis,bflabl,bflabl)
         end if
c
c        ----- write nuclear + effective core integrals to read-write file.
c
         call iosys('read real "potential integrals" from rwf',
     $        nnp,z,0,' ')
         call vadd(s,s,z,nnp)
         call iosys('write real "potential integrals" on rwf',
     $        nnp,s,0,' ')
      endif
c
c     ----- re-pack and write out one electron integrals -----
c           if drop occurs
c 
      pkindx=1
      if (drop) then
         call iosys('read integer "packing index vector" from rwf',
     $               nbasis,iz(pkindx),0,' ')
         call iosys('read integer "truncated number of basis'
     $              //' functions" from rwf',1,newnbf,0,' ')
         newnnp=(newnbf+1)*newnbf/2
c
         call iosys('read real "overlap integrals" from rwf',
     $               nnp,s,0,' ')
         write(iout,*)' packing the overlap integrals'
         call fixone(s,nbasis,iz(pkindx))
         call iosys('write real "overlap integrals" to rwf',
     $               newnnp,s,0,' ')
c
         call iosys('read real "kinetic integrals" from rwf',
     $               nnp,s,0,' ')
         write(iout,*)' packing the kinetic energy integrals'
         call fixone(s,nbasis,iz(pkindx))
         call iosys('write real "kinetic integrals" to rwf',
     $               newnnp,s,0,' ')
c
         call iosys('read real "potential integrals" from rwf',
     $               nnp,s,0,' ')
         write(iout,*)' packing the potential energy integrals'
         call fixone(s,nbasis,iz(pkindx))
         call iosys('write real "potential integrals" to rwf',
     $               newnnp,s,0,' ')
      endif
c
c
      return
      end
