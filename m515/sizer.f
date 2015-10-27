*deck @(#)sizer.f	5.1  11/28/95
      subroutine sizer(noprim,nbtype,nocart,maxmom,
     $                 nat,npf,nbf,ndmat,biggy)
c***begin prologue     sizer.f
c***date written       940115      (yymmdd)
c***revision date      11/6/94      
c
c***keywords           direct, j-matrix
c***author             martin, richard(lanl)
c***source             @(#)sizer.f	5.1   11/28/95
c***purpose            runs through the direct integrals loop
c                      to find maximum core needed.
c***description
c***references         
c
c***routines called    
c
c***end prologue       sizer.f
c
      implicit none
c     --- input variables -----
      integer nbtype,nat,npf,nbf,ndmat
c     --- input arrays (unmodified) ---
      integer noprim(nat,nbtype)
      integer nocart(nbtype)
      integer maxmom(nbtype)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
      integer biggy
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer iadtwp,wpadti
      integer iatom,jatom,katom,latom
      integer itype,jtype,ltype,ktype
      integer imax,jmax,kmax,lmax
      integer nfi,nfj,nfk,nfl
      integer nprimi,nprimj,nprimk,npriml
      integer jtmax,ktmax,latmax,ltmax
      integer nmax,mmax
      integer nroots,nij,nkl
      integer dij,dkl,jij,jkl,ijindx,klindx
      integer aij,bkl
      integer ar,xyza,xyzam1,xyzam3,br,xyzb,xyzbm1,xyzbm3
      integer tmptop,len,index,ab,wt,urho
      integer f00,denom,b00,b10,bp01,c00,cp00,g,h,i2
      integer top
      integer expon,aa,bb,rho,t0,t1,t2,t3,t4,t5,t6,t7,t8
      integer kl,nv
      integer aplusb
c
c
      logical debug
c
      parameter (debug=.false.)
c
      common /io/inp,iout
c
c     ----- make sure we have enough core -----
      biggy=0
      do 8000 iatom=1,nat
         do 7000 itype=1,nbtype
            nprimi=noprim(iatom,itype)
            if (nprimi.le.0) go to 7000
            nfi=nocart(itype)
            imax=maxmom(itype)
            do 6000 jatom=1,iatom
               if (jatom.eq.iatom) then
                  jtmax=itype
               else
                  jtmax=nbtype
               end if
               do 5000 jtype=1,jtmax
                  nprimj=noprim(jatom,jtype)
                  if (nprimj.le.0) go to 5000
                  nfj=nocart(jtype)
                  jmax=maxmom(jtype)
c
                  nij=nprimi*nprimj
                  dij=1
                  do 4000 katom=1,iatom
                     if (katom.eq.iatom) then
                        ktmax=itype
                     else
                        ktmax=nbtype
                     end if
                     do 3000 ktype=1,ktmax
                        nprimk=noprim(katom,ktype)
                        if (nprimk.le.0) go to 3000
                        nfk=nocart(ktype)
                        kmax=maxmom(ktype)
                        if (katom.eq.iatom.and.ktype.eq.itype) then
                           latmax=jatom
                        else
                           latmax=katom
                        end if
                        do 2000 latom=1,latmax
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
                              nfl=nocart(ltype)
                              lmax=maxmom(ltype)
c
                              nkl=nprimk*npriml
c
c
               nmax=imax+jmax
               mmax=kmax+lmax
               nroots=(nmax+mmax)/2+1
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
               kl=nkl
               nv=nij*nkl
c
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
               top=t8+nv
               biggy=max(biggy,top)
c
               g=ab
               i2=g
               h=i2+
     $              3*(imax+1)*(jmax+1)*(mmax+1)*(lmax+1)*nv*nroots
               top=h+3*(nmax+1)*(jmax+1)*(mmax+1)*nv*nroots
               biggy=max(biggy,top)
c
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
      t1=1
      t2=t1+npf*npf
      top=t2+npf*nbf
      biggy=max(biggy,top)
c
c
      return
      end
