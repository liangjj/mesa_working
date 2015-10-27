*deck @(#)sizer1.f	5.1  11/6/94
      subroutine sizer1(ptprim,noprim,nbtype,nocont,ptcont,ex,cont,
     #                  ncont,c,ints,labels,lenbuf,start,nx,ny,nz,
     #                  lenxyz,nobf,nocart,mintyp,minmom,maxmom,
     #                  z,maxcor,nat,num,nnp,nprim,ops,cutexp,maxc)
c
      implicit integer(a-z)
c
      character*(*) ops
      real*8 ex(nprim),cont(ncont),c(3,nat),ints(lenbuf),z(maxcor)
      real*8 cutexp
      integer ptprim(nat,nbtype),noprim(nat,nbtype),nocont(nat,nbtype)
      integer ptcont(nat,nbtype),labels(lenbuf),nobf(nbtype)
      integer start(nat,nbtype),nocart(nbtype),mintyp(nbtype)
      integer minmom(nbtype),maxmom(nbtype)
      integer nx(lenxyz),ny(lenxyz),nz(lenxyz)
      integer symcen(4),angmom(4)
c
      common /io/     inp,iout
c
c     ----- timing -----
c
c
c     ----- initialize counters -----
c
      maxc=0
      maxc1=0
      maxc2=0
      maxc3=0
      maxc4=0
      numbuf=1
      ntotal=0
      nonzer=0
      nactul=0
c
      do 8000 iatom=1,nat
         symcen(1)=iatom
         do 7000 itype=1,nbtype
            nprimi=noprim(iatom,itype)
            if (nprimi.le.0) go to 7000
            angmom(1)=itype
            imax=maxmom(itype)
            nconti=nocont(iatom,itype)
            do 6000 jatom=1,iatom
               symcen(2)=jatom
               if (jatom.eq.iatom) then
                  jtmax=itype
               else
                  jtmax=nbtype
               end if
               do 5000 jtype=1,jtmax
                  nprimj=noprim(jatom,jtype)
                  if (nprimj.eq.0) go to 5000
                  angmom(2)=jtype
                  jmax=maxmom(jtype)
                  ncontj=nocont(jatom,jtype)
                  do 4000 katom=1,iatom
                     symcen(3)=katom
                     if (katom.eq.iatom) then
                        ktmax=itype
                     else
                        ktmax=nbtype
                     end if
                     do 3000 ktype=1,ktmax
                        nprimk=noprim(katom,ktype)
                        if (nprimk.le.0) go to 3000
                        angmom(3)=ktype
                        kmax=maxmom(ktype)
                        ncontk=nocont(katom,ktype)
                        if (katom.eq.iatom.and.ktype.eq.itype) then
                           latmax=jatom
                        else
                           latmax=katom
                        end if
                        do 2000 latom=1,latmax
                           symcen(4)=latom
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
                              angmom(4)=ltype
                              lmax=maxmom(ltype)
                              ncontl=nocont(latom,ltype)
c
               nmax=imax+jmax
               mmax=kmax+lmax
               nroots=(nmax+mmax)/2+1
               npint=nprimi*nprimj*nprimk*npriml
               ncint=nconti*ncontj*ncontk*ncontl
               nij=nprimi*nprimj
               nkl=nprimk*npriml
               lenblk=nocart(angmom(1))*nocart(angmom(2))*
     #                nocart(angmom(3))*nocart(angmom(4))
               numint=ncint*lenblk
c
c     ----- core allocation. nb there is a great deal of overlapping --
c
               a=1
               ar=a+nij
               xyza=ar+nij
               xyzam1=xyza+nij*3
               xyzam3=xyzam1+nij*3
               b=xyzam3+nij*3
               br=b+nkl
               xyzb=br+nkl
               xyzbm1=xyzb+nkl*3
               xyzbm3=xyzbm1+nkl*3
               half=xyzbm3+nkl*3
               test=half+nconti*ncontj*nkl*lenblk
c
               numkl=1
               npass=(nkl+numkl-1)/numkl
               nijkl=nij*numkl
               lenv=nijkl*nroots
c
               g=test+nijkl
               ab=g
               aplusb=ab+nijkl
               urho=aplusb+nijkl
               wt=urho+lenv
               denom=wt+lenv
               t1=denom+lenv
               t2=t1+nijkl
               t3=t2+nijkl
               t4=t3+nijkl
               t5=t4+lenv
               t6=t5+lenv
               t7=t6+nijkl
               t8=t7+nijkl
               t9=t8+nijkl
c
               if (nmax+mmax.eq.0) then
                  f00=denom
                  b00=1
                  b10=1
                  bp01=1
                  c00=1
                  cp00=1
                  top1=t9+nijkl
               else
                  f00=max(g+3*(nmax+1)*(mmax+1)*lenv,t9+nijkl)
                  b00=f00+lenv
                  b10=b00+lenv
                  bp01=b10+lenv
                  c00=bp01+lenv
                  cp00=c00+lenv*3
                  top1=cp00+lenv*3
               end if
c
               if (nmax+mmax.eq.0) then
                  temp1=ab
                  i2=1
                  h=1
                  temp=1
                  tp1=1
                  top2=f00+lenv
               else
                  i2=g
                  h=i2+3*(imax+1)*(jmax+1)*(mmax+1)*(lmax+1)*lenv
                  lentmp=nijkl
                  if (lentmp.lt.nijkl) call lnkerr('internal error '//
     #                 'with length of temps for vfmint')
                  temp=h
                  temp1=temp+lentmp
                  top2=max(h+3*(nmax+1)*(jmax+1)*(mmax+1)*lenv,
     #                    temp1+lentmp)
               end if
c
               conint=half+nconti*ncontj*nkl*lenblk
               tmp1=conint+ncint*lenblk
               lent=maxcor-tmp1+1
               top3=tmp1+ncontk*npriml
c
               numint=ncint*lenblk
               plbls=half
               flbls=plbls+numint
               if (flbls.gt.conint) call lnkerr('internal error in '//
     #                                  'driver with flbls')
               if (flbls+lenblk*4.le.conint) then
                  top4=conint+numint
               else
                  flbls=conint+numint
                  top4=flbls+lenblk*4
               end if
c
               top=max(top1,top2,top3,top4)
               maxc=max(top,maxc)
               maxc1=max(maxc1,top1)
               maxc2=max(maxc2,top2)
               maxc3=max(maxc3,top3)
               maxc4=max(maxc4,top4)
 1000                      continue
 2000                   continue
 3000                continue
 4000             continue
 5000          continue
 6000       continue
 7000    continue
 8000 continue
c
c     ----- timing -----
c
c
c
      return
      end
