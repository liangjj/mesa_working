*deck %W%  %G%
      subroutine sizer(noprim,nbtype,
     $                  nx,ny,nz,
     $                  lenxyz,nobf,nocart,mintyp,minmom,maxmom,
     $                  nat,npf,nnprim,nprim,ops,cutexp,
     $                  nderiv,prnt,ndmat,
     $                  nshell,mincor)
c***begin prologue     sizer.f
c***date written       851113   (yymmdd)
c***revision date      11/6/94
c  22 march 1987    pws at lanl
c     modifying 'driver' to run through and calculate core needed.
c***keywords
c***author             saxe, paul    (lanl)
c***source             %W%   %G%
c***purpose
c***description
c***references         (none)
c
c***routines called    (none)
c
c***end prologue       sizer.f
      implicit none
c     --- input variables -----
      integer lenxyz,nat,npf,nnprim,nprim,nderiv,ndmat,nshell,mincor
      integer nbtype
      real*8 cutexp
c     --- input arrays (unmodified) ---
      character*(*) ops
      logical prnt
      integer noprim(nat,nbtype)
      integer nobf(nbtype)
      integer nocart(nbtype),mintyp(nbtype)
      integer minmom(nbtype),maxmom(nbtype)
      integer nx(lenxyz),ny(lenxyz),nz(lenxyz)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer iatom,itype,jatom,jtype,katom,ktype,latom,ltype
      integer jtmax,ltmax,latmax,ktmax
      integer imax,jmax,kmax,lmax,nmax,mmax,nroots
      integer nprimi,nprimj,nprimk,npriml,nfi,nfj,nfk,nfl
      integer npint,lenblk,numint,ntotal
      integer ndcen
      integer dij,dkl,ijindx,klindx,aij,ar,xyza,xyzam1,xyzam3
      integer bkl,br,xyzb,xyzbm1,xyzbm3,tmptop
      integer fac1,fac2,len,lenv,index,twopdm
      integer dg,g,ab,aplusb,urho,wt,denom,expon,a,b,rho
      integer t1,t2,t3,t4,t5,t6,t7,t8
      integer rhotsq,expnts,camcb
      integer f00,b00,b10,bp01,c00,cp00,top1,dc00,dcp00
      integer temp1,i2,h,temp,tp1,top2
      integer di,dh,d1exp,maxc,maxc1,maxc2,top
      integer numbuf,nonzer,nactul
      integer symcen(4),angmom(4),dercen(4),dermom(4)
      integer npass,flip,nij,nkl
      integer wptoin,wpadti,iadtwp
c
      common /io/     inp,iout
c
c     --- initialize counters ---
c
      maxc=0
      maxc1=0
      maxc2=0
      numbuf=1
      ntotal=0
      nonzer=0
      nactul=0
c
      do 8000 iatom=1,nat
         symcen(1)=iatom
         do 7000 itype=1,nbtype
            if (noprim(iatom,itype).le.0) go to 7000
            angmom(1)=itype
            do 6000 jatom=1,iatom
               symcen(2)=jatom
               if (jatom.eq.iatom) then
                  jtmax=itype
               else
                  jtmax=nbtype
               end if
               do 5000 jtype=1,jtmax
                  if (noprim(jatom,jtype).le.0) go to 5000
                  angmom(2)=jtype
                  do 4000 katom=1,iatom
                     symcen(3)=katom
                     if (katom.eq.iatom) then
                        ktmax=itype
                     else
                        ktmax=nbtype
                     end if
                     do 3000 ktype=1,ktmax
                        if (noprim(katom,ktype).le.0) go to 3000
                        angmom(3)=ktype
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
                              if (noprim(latom,ltype).le.0) go to 1000
                              angmom(4)=ltype
c
c              --- order centres and number of derivatives ---
               call redund(symcen,angmom,dercen,dermom,npass,flip)
c
               imax=maxmom(dermom(1))
               jmax=maxmom(dermom(2))
               kmax=maxmom(dermom(3))
               lmax=maxmom(dermom(4))
               nmax=imax+jmax
               mmax=kmax+lmax
               nroots=(nmax+mmax+nderiv)/2+1
               nprimi=noprim(dercen(1),dermom(1))
               nprimj=noprim(dercen(2),dermom(2))
               nprimk=noprim(dercen(3),dermom(3))
               npriml=noprim(dercen(4),dermom(4))
               nfi=nocart(dermom(1))
               nfj=nocart(dermom(2))
               nfk=nocart(dermom(3))
               nfl=nocart(dermom(4))
               nij=nprimi*nprimj
               nkl=nprimk*npriml
               npint=nij*nkl
               lenblk=nfi*nfj*nfk*nfl
               numint=npint*lenblk
               ntotal=ntotal+numint
c
c              --- for derivatives, the number of centres to
c                  differentiate
               if (nderiv.eq.0.or.npass.eq.0) then
                  ndcen=0
               else
                  if (npass.eq.1) then
                     ndcen=3
                  else if (npass.eq.4) then
                     ndcen=1
                  else
                     ndcen=2
                  end if
               end if
c
c              --- core allocation. nb there is a great deal of overlapping ---
               dij=1
               dkl=dij+nprimi*nprimj*nfi*nfj*ndmat
               ijindx=wpadti(dkl+nprimk*npriml*nfk*nfl*ndmat)
               klindx=ijindx+nij*2
               aij=iadtwp(klindx+nkl*2)
               ar=aij+nij
               xyza=ar+nij
               xyzam1=xyza+nij*3
               xyzam3=xyzam1+nij*3
               bkl=xyzam3+nij*3
               br=bkl+nkl
               xyzb=br+nkl
               xyzbm1=xyzb+nkl*3
               xyzbm3=xyzbm1+nkl*3
               tmptop=wpadti(xyzbm3+nkl*3)
c
c              --- work out how long the vectors can be ---
               if (nmax+mmax.lt.0) then
                  fac1=10+5*nroots
                  fac2=3+3*nroots
               else
                  if (nderiv.eq.0.or.npass.eq.0) then
                     fac1=6+wptoin(lenblk+10*nroots+
     $                   max(3*(nmax+1)*(mmax+1)*nroots,
     $                         12+5*nroots))
                     fac2=6+wptoin(lenblk+
     $                   nroots*3*(jmax+1)*(mmax+1)*((imax+1)*(lmax+1)+
     $                                               (nmax+1)))
                  else
                     fac1=6+wptoin(lenblk+nroots*(10+2*ndcen)+
     $                   max(3*(nmax+1)*(mmax+1)*nroots*(1+ndcen),
     $                         12+5*nroots))
                     fac2=max(6+wptoin(lenblk+
     $                   nroots*3*(imax+1)*(jmax+1)*(mmax+1)*(lmax+1)*
     $                                             (ndcen+1)+
     $                   nroots*3*(nmax+1)*(jmax+1)*(mmax+1)*(ndcen+1)),
     $                   fac1)+wptoin(nroots*3*ndcen)
                  end if
               end if
c
               len=nij
               lenv=len*nroots
c
               index=tmptop
               twopdm=iadtwp(index+len*6)
               if (nderiv.eq.0.or.npass.eq.0) then
                  dg=twopdm+len*lenblk
                  g=dg
               else
                  dg=twopdm+len*lenblk
                  g=dg+lenv*3*(nmax+1)*(mmax+1)*ndcen
               end if
               ab=dg
               aplusb=ab+len
               urho=aplusb+len
               wt=urho+lenv
               denom=wt+lenv
               expon=denom+lenv
               a=expon+len
               b=a+len
               rho=b+len
               t1=rho+len
               t2=t1+len
               t3=t2+len
               t4=t3+len
               t5=t4+lenv
               t6=t5+lenv
               t7=t6+len
               t8=t7+len
c
               if (nderiv.eq.0.or.npass.eq.0) then
                  rhotsq=1
                  expnts=1
               else
                  rhotsq=t1
                  camcb=rhotsq+len*nroots
                  expnts=camcb+len*3
               end if
c
               if (nmax+mmax.lt.0) then
                  f00=denom
                  b00=1
                  b10=1
                  bp01=1
                  c00=1
                  cp00=1
                  top1=t8+len
               else
                  f00=max(g+3*(nmax+1)*(mmax+1)*lenv,t8+len)
                  b00=f00+lenv
                  b10=b00+lenv
                  bp01=b10+lenv
                  c00=bp01+lenv
                  cp00=c00+lenv*3
                  top1=wpadti(cp00+lenv*3)
                  if (nderiv.ne.0.and.npass.ne.0) then
                     dc00=iadtwp(top1)
                     dcp00=dc00+nroots*len*ndcen
                     top1=wpadti(dcp00+nroots*len*ndcen)
                  end if
               end if
c
               if (abs(top1-(tmptop+fac1*len)).gt.4) then
                  call plnkerr('internal error in sizer with top1',201)
               end if
c
               if (nmax+mmax.lt.0) then
                  temp1=ab
                  i2=1
                  h=1
                  temp=1
                  tp1=1
                  top2=f00+lenv
               else
                  if (nderiv.eq.0.or.npass.eq.0) then
                     i2=g
                     h=i2+3*(imax+1)*(jmax+1)*(mmax+1)*(lmax+1)*lenv
                     top2=wpadti(h+3*(nmax+1)*(jmax+1)*(mmax+1)*lenv)
                  else
                     di=dg
                     i2=di+nroots*len*3*(imax+1)*(jmax+1)*(mmax+1)*
     $                                  (lmax+1)*ndcen
                     dh=i2+nroots*len*3*(imax+1)*(jmax+1)*(mmax+1)*
     *                                  (lmax+1)
                     h=dh+nroots*len*3*(nmax+1)*(jmax+1)*(mmax+1)*ndcen
                     d1exp=max(iadtwp(top1),
     $                     h+nroots*len*3*(nmax+1)*(jmax+1)*(mmax+1))
                     top2=wpadti(d1exp+nroots*len*3*ndcen)
                  end if
               end if
c
               if (abs(top2-(tmptop+fac2*len)).gt.4) then
                  call plnkerr('internal error in sizer with top2',202)
               end if
c
               top=max(top1,top2)
               maxc=max(top,maxc)
               maxc1=max(maxc1,top1)
               maxc2=max(maxc2,top2)
c
 1000                      continue
 2000                   continue
 3000                continue
 4000             continue
 5000          continue
 6000       continue
 7000    continue
 8000 continue
c
      mincor=maxc+200
c
c
      return
      end
