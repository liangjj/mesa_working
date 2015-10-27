*deck @(#)oneint.f	5.2  2/5/95
      subroutine oneint(c,ex,z,iz,cont,s,ptprim,noprim,nocont,ptcont,
     $                  nat,nprim,maxcor,ntypes,nbtype,nnp,ncont,
     $                  start,nbasis,zan,nocart,nobf,maxmom,mintyp,
     $                  nx,ny,nz,minmom,dolp,bflabl,ops,ncharge)
c***begin prologue     oneint.f
c***date written       840723  
c***revision date      2/5/95      
c   january 8, 1994    rlm at lanl
c      adding ncharge variable which flags solvent surface charges
c      which should be included in the one-electron integrals
c
c***keywords           
c***author             saxe, paul and martin, richard (lanl) 
c***source             @(#)oneint.f	5.2   2/5/95
c***purpose            
c***description
c      module to form the overlap, kinetic-energy and potential-energy
c      one-electron integrals over generally contracted gaussian basis
c      sets.
c     
c***references
c
c***routines called
c
c***end prologue       oneint.f
      implicit none
c     --- input variables -----
      integer nat,nprim,maxcor,ntypes,nbtype,nnp,ncont
      integer nbasis,ncharge
      character*(*) ops
      logical dolp
c     --- input arrays (unmodified) ---
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(ntypes)
      integer nobf(ntypes),maxmom(ntypes),mintyp(ntypes),minmom(ntypes)
      integer nx(*),ny(*),nz(*)
      real*8 c(3,nat+ncharge),ex(nprim),cont(ncont)
      real*8 zan(nat+ncharge)
      character*(*) bflabl(*)
c     --- input arrays (scratch) ---
      integer iz(*)
      real*8 z(maxcor),s(nnp)
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer mxpts
      integer iatom,jatom,itype,jtype,jtypmx,imax,jmax,nprimi,nprimj
      integer nconti,ncontj,lenblk
      integer a,alpha,b,npint,prmint,xyz,len1,ainv,xyza,t1,xyz0
      integer expon,conint,bp2,bm2,aiaj,rysrt,ryswt,nroots,xyz1
      integer nlp,ntpse,clp,zlp,qq,pkindx,newnbf,newnnp
      integer top1,top2,tmp1
      integer len,imin,jj
      logical drop,logkey
      character refops*8
      real*8 h(21),wt(21)
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
      common/io/inp,iout
c
   10 format(5x,'overlap integrals:')
   20 format(5x,'kinetic energy integrals:')
   30 format(5x,'potential energy integrals:')
   40 format(5x,'effective core potential integrals:')
c
c     --- are we to drop function components?
      drop=logkey(ops,'drop',.false.,' ')
c
c     --- form the overlap integrals ---
      do 9 iatom=1,nat
         do 9999 itype=1,nbtype
         write(iout,*) 'iatom,itype',iatom,itype
         write(iout,*) 'noprim(iatom,itype)',
     $                  (noprim(iatom,itype))
         write(iout,*) 'nocont(iatom,itype)',
     $                  (nocont(iatom,itype))
         write(iout,*) 'ptprim(iatom,itype)',
     $                  (ptprim(iatom,itype))
         write(iout,*) 'ptcont(iatom,itype)',
     $                  (ptcont(iatom,itype))
         write(iout,*) 'start(iatom,itype)',
     $                  (start(iatom,itype))
         len=noprim(iatom,itype)*nocont(iatom,itype)
         imin=minmom(itype)
         imax=maxmom(itype)
c        each of these blocks is dimensioned cont(nprimi,nconti,imin:imax)
c        where nprimi=noprim(iatom,itype),nconti=nocont(iatom,itype), and imin and imax are
c        the miniumum and maximum values of l in the shell -- i.e. 0-0 for s; 1-1 for p, 0-1 for sp, etc.
         write(iout,*) 'cont(ptcont(iatom,itype))',
     $           (cont(ptcont(iatom,itype)+jj-1),jj=1,len*(imin-imax+1))
 9999    continue
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
c                 --- allocate core for temporary vectors, etc. ---
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
c                 --- form the two-dimensional integrals ---
                  call sints(iatom,jatom,imax,jmax,nprimi,nprimj,
     $                       ptprim(iatom,itype),ptprim(jatom,jtype),
     $                       ex,z(alpha),c,z(ainv),z(xyza),
     $                       z(expon),z(xyz),npint,nprim,nat,
     $                       nbtype,h,wt,mxpts,z(a),z(b))
c
c                 --- form the primitive integrals ---
                  call fmonel(z(prmint),z(xyz),npint,lenblk,imax,jmax,
     $                   mintyp(itype),mintyp(itype)+nocart(itype)-1,
     $                   mintyp(jtype),mintyp(jtype)+nocart(jtype)-1,
     $                   nx,ny,nz)
c
c                 --- transform to contracted functions ---
                  call trans1(z(prmint),z(conint),nprimi,nprimj,nconti,
     $                        ncontj,cont(ptcont(iatom,itype)),
     $                        cont(ptcont(jatom,jtype)),z(tmp1),len1,
     $                        lenblk,minmom(itype),maxmom(itype),
     $                        minmom(jtype),maxmom(jtype),nocart)
c
c                 --- transfer integrals to total array ---
                  call put1el(s,z(conint),start,iatom,jatom,itype,jtype,
     $                        nconti,ncontj,nnp,lenblk,nat,nbtype,nobf)
    6          continue
    7       continue
    8    continue
    9 continue
c
c     --- print the integrals ---
      if(logkey(ops,'print=int=s',.false.,refops)) then
         write(iout,10)
         call trtosq(z(prmint),s,nbasis,nnp)
         call wlmat(z(prmint),nbasis,nbasis,bflabl,bflabl)
      end if
c
c     --- put overlap integrals on the read-write file ---
      call iosys('write real "overlap integrals" on rwf',
     $            nnp,s,0,' ')
c
c     --- form the kinetic-energy integrals ---
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
c                 --- allocate core for temporary vectors, etc. ---
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
c
c                 --- form the two-dimensional integrals ---
                  call tints(iatom,jatom,imax,jmax,nprimi,nprimj,
     $                       ptprim(iatom,itype),ptprim(jatom,jtype),
     $                       ex,z(alpha),c,z(ainv),z(xyza),
     $                       z(expon),z(xyz),npint,nprim,nat,
     $                       nbtype,h,wt,mxpts,z(a),z(b),z(bp2),z(bm2))
c
c                 --- form the primitive integrals ---
                  call fmt(z(prmint),z(xyz),npint,lenblk,
     $                     imax,jmax,z(t1),
     $                   mintyp(itype),mintyp(itype)+nocart(itype)-1,
     $                   mintyp(jtype),mintyp(jtype)+nocart(jtype)-1,
     $                   nx,ny,nz)
c
c                 --- transform to contracted functions ---
                  call trans1(z(prmint),z(conint),nprimi,nprimj,nconti,
     $                        ncontj,cont(ptcont(iatom,itype)),
     $                        cont(ptcont(jatom,jtype)),z(tmp1),len1,
     $                        lenblk,minmom(itype),maxmom(itype),
     $                        minmom(jtype),maxmom(jtype),nocart)
c
c                 --- transfer integrals to total array ---
                  call put1el(s,z(conint),start,iatom,jatom,itype,jtype,
     $                        nconti,ncontj,nnp,lenblk,nat,nbtype,nobf)
c
   16          continue
   17       continue
   18    continue
   19 continue
c
c     --- print the integrals ---
      if(logkey(ops,'print=int=t',.false.,refops)) then
         write(iout,20)
         call trtosq(z(prmint),s,nbasis,nnp)
         call wlmat(z(prmint),nbasis,nbasis,bflabl,bflabl)
      end if
c
c     --- put kinetic-energy integrals on the read-write file ---
      call iosys('write real "kinetic integrals" on rwf',
     $            nnp,s,0,' ')
c
c     --- form the potential-energy integrals ---
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
c                 --- allocate core for temporary vectors, etc. ---
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
c
c                 --- form the primitive integrals ---
                  call vints(iatom,jatom,imax,jmax,nprimi,nprimj,
     $                       ptprim(iatom,itype),ptprim(jatom,jtype),
     $                       ex,z(alpha),c,z(ainv),z(xyza),
     $                       z(expon),z(xyz),npint,nprim,nat,
     $                       nbtype,h,wt,mxpts,z(a),z(b),z(aiaj),
     $                       z(xyz0),z(rysrt),z(ryswt),nroots,
     $                       z(prmint),lenblk,zan,z(xyz1),z(t1),
     $                 mintyp(itype),mintyp(itype)+nocart(itype)-1,
     $                 mintyp(jtype),mintyp(jtype)+nocart(jtype)-1,
     $                 nx,ny,nz,ncharge)
c
c                 --- transform to contracted functions ---
                  call trans1(z(prmint),z(conint),nprimi,nprimj,nconti,
     $                        ncontj,cont(ptcont(iatom,itype)),
     $                        cont(ptcont(jatom,jtype)),z(tmp1),len1,
     $                        lenblk,minmom(itype),maxmom(itype),
     $                        minmom(jtype),maxmom(jtype),nocart)
c
c                 --- transfer integrals to total array ---
                  call put1el(s,z(conint),start,iatom,jatom,itype,jtype,
     $                        nconti,ncontj,nnp,lenblk,nat,nbtype,nobf)
c
   26          continue
   27       continue
   28    continue
   29 continue
c
c     --- print the integrals ---
      if(logkey(ops,'print=int=v',.false.,refops)) then
         write(iout,30)
         call trtosq(z(prmint),s,nbasis,nnp)
         call wlmat(z(prmint),nbasis,nbasis,bflabl,bflabl)
      end if
c
c     --- put nuclear-potential integrals on the read-write file ----
      call iosys('write real "potential integrals" on rwf',
     $            nnp,s,0,' ')
c
c     --- are there any ecp centers?
      if(dolp) then
c
c        --- form the effective core potential integrals ---
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
c                    --- allocate core for temporary vectors, etc. ---
                     prmint=1
                     ntpse=prmint+npint*lenblk
                     nlp=ntpse+10
                     zlp=nlp+100
                     clp=zlp+100
                     qq=clp+100
                     top1=qq+11*9*9*npint
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
c                    --- form the primitive integrals ---
                     call lpints(iatom,jatom,imax,jmax,nprimi,nprimj,
     $                       nprim,ncont,nat,nbtype,ntypes,npint,lenblk,
     $                       ptprim(iatom,itype),ptprim(jatom,jtype),
     $                       c,ex,noprim,ptprim,ptcont,cont,maxmom,
     $                       mintyp(itype),
     $                       mintyp(itype)+nocart(itype)-1,
     $                       mintyp(jtype),
     $                       mintyp(jtype)+nocart(jtype)-1,
     $                       nx,ny,nz,z(prmint),z(ntpse),z(nlp),z(zlp),
     $                       z(clp),z(qq))
c
c                    --- transform to contracted functions ---
                     call trans1(z(prmint),z(conint),nprimi,nprimj,
     $                        nconti,ncontj,cont(ptcont(iatom,itype)),
     $                        cont(ptcont(jatom,jtype)),z(tmp1),len1,
     $                        lenblk,minmom(itype),maxmom(itype),
     $                        minmom(jtype),maxmom(jtype),nocart)
c
c                    --- transfer integrals to total array ---
                     call put1el(s,z(conint),start,iatom,jatom,itype,
     $                        jtype,nconti,ncontj,nnp,lenblk,nat,nbtype,
     $                        nobf)
   36             continue
   37          continue
   38       continue
   39    continue
c
c        --- print the integrals ---
      if(logkey(ops,'print=int=ecp',.false.,refops)) then
         write(iout,40)
         call trtosq(z(prmint),s,nbasis,nnp)
         call wlmat(z(prmint),nbasis,nbasis,bflabl,bflabl)
         end if
c
c        --- write nuclear + effective core integrals to read-write file.
         call iosys('read real "potential integrals" from rwf',
     $               nnp,z,0,' ')
         call vadd(s,s,z,nnp)
         call iosys('write real "potential integrals" on rwf',
     $               nnp,s,0,' ')
      endif
c
c
c     --- re-pack and write out one electron integrals ---
c         if drop occurs
      pkindx=1
      if (drop) then
         call iosys('read integer "packing index vector" from rwf',
     $               nbasis,iz(pkindx),0,' ')
         call iosys('read integer "truncated number of basis'
     $              //' functions" from rwf',1,newnbf,0,' ')
         newnnp=(newnbf+1)*newnbf/2
c
c
         call iosys('read real "overlap integrals" from rwf',
     $              nnp,s,0,' ')
         call fixone(s,nbasis,iz(pkindx))
         call iosys('write real "overlap integrals" to rwf',
     $              newnnp,s,0,' ')
c
c
         call iosys('read real "kinetic integrals" from rwf',
     $               nnp,s,0,' ')
         call fixone(s,nbasis,iz(pkindx))
         call iosys('write real "kinetic integrals" to rwf',
     $               newnnp,s,0,' ')
c
c
         call iosys('read real "potential integrals" from rwf',
     $               nnp,s,0,' ')
         call fixone(s,nbasis,iz(pkindx))
         call iosys('write real "potential integrals" to rwf',
     $               newnnp,s,0,' ')
      endif
c
c
      return
      end
