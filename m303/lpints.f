*deck @(#)lpints.f	5.1  11/6/94
      subroutine lpints(iatom,jatom,imax,jmax,nprimi,nprimj,
     $                   nprim,ncont,nat,nbtype,ntypes,maxp2,maxi2,
     $                   i1,j1,c,ex,noprim,ptprim,ptcont,cont,
     $                   maxmom,mini,maxi,minj,maxj,nx,ny,nz,sp,
     $                   ntpse,nlp,zlp,clp,qq,
     $                   scr,dmaxi2,nderiv,dmini,dmaxi,dminj,dmaxj,
     $                   conint,nconti,ncontj,
     $                   itype,jtype,acf,bcf,tmp1,len1,imin,jmin,nocart,
     $                   ds,nnp,start,nobf)
c
      implicit real*8(a-h,o-z)
      integer ptprim,ptcont,dmaxi2,start,dmini,dmaxi,dminj,dmaxj
      integer atom,coord,junk
      logical skip
c
c     controls computation of pseudopotential integrals.
c     details of the method are presented in
c         l.e. mcmurchie and e.r. davidson, j. comp. phys.,44,289(1981).
c         r.l. martin, to be published.
c
      common/center/xa,ya,za,xb,yb,zb,xc,yc,zc
      common/dist/cax,cay,caz,ca,ca2,cbx,cby,cbz,cb,cb2
      common/pifac/pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
      common/const/zero,one,two,three,four,five,six,ten
      common/prims/igbegn,igend,jgbegn,jgend
      common/limit/istart,jstart,iend,jend
      common/angmax/lamax,lbmax
      dimension ex(nprim),c(3,nat),noprim(nat,ntypes)
      dimension ptprim(nat,ntypes),ptcont(nat,ntypes),cont(ncont)
      dimension maxmom(ntypes)
      dimension nx(*),ny(*),nz(*),nlp(*),zlp(*),clp(*)
      dimension ntpse(*),ds(nnp,3,nat)
c
      dimension qq(*),sp(maxp2,maxi2,7)
      dimension conint(nconti,ncontj,maxi2),tmp1(len1),nocart(nbtype)
      real*8 acf(nprimi,nprimj,imin:imax),bcf(nprimj,ncontj,jmin:jmax)
      dimension scr(maxp2,dmaxi2),start(nat,nbtype)
      dimension nobf(nbtype)
      integer centre(3)
c
      parameter (half=5.0d-01)
c
      common/io/inp,iout
c
      centre(1)=iatom
      centre(2)=jatom
      call rzero(sp,maxp2*maxi2*7)
      lamax=imax+1+nderiv
      lbmax=jmax+1+nderiv
      istart=mini
      iend=maxi
      jstart=minj
      jend=maxj
      if(nderiv.ne.0) then
         istart=dmini
         iend=dmaxi
         jstart=dminj
         jend=dmaxj
      endif
      igbegn=i1
      igend=i1+nprimi-1
      jgbegn=j1
      jgend=j1+nprimj-1
      xa=c(1,iatom)
      ya=c(2,iatom)
      za=c(3,iatom)
      xb=c(1,jatom)
      yb=c(2,jatom)
      zb=c(3,jatom)
c
      do 100 katom=1,nat
         centre(3)=katom
         call rzero(scr,maxp2*dmaxi2)
         if (nderiv .eq.2) then
            call lnkerr('m303: nderiv should not be 2 here')
c            call rzero(sp(1,1,2),maxp2*maxi2*27)
         else
            call rzero(sp(1,1,2),maxp2*maxi2*6)
         endif
         xc=c(1,katom)
         yc=c(2,katom)
         zc=c(3,katom)
         cax=xc-xa
         cay=yc-ya
         caz=zc-za
         ca2=cax*cax+cay*cay+caz*caz
         cbx=xc-xb
         cby=yc-yb
         cbz=zc-zb
         cb2=cbx*cbx+cby*cby+cbz*cbz
         ca=sqrt(ca2)
         cb=sqrt(cb2)
c
         nt=0
         kk=0
         lmax=0
         skip=.true.
         do 10 ktype=nbtype+1,ntypes
            if(noprim(katom,ktype).ne.0) then
               skip=.false.
               lmax=max(lmax,maxmom(ktype))
               nt=nt+1
               ntpse(nt)=noprim(katom,ktype)
               do 5 jj=1,noprim(katom,ktype)
                  kk=kk+1
                  nlp(kk)=int(cont(ptcont(katom,ktype)
     $                             +noprim(katom,ktype)+jj-1)+half)
                  clp(kk)=cont(ptcont(katom,ktype)+jj-1)
                  zlp(kk)=ex(ptprim(katom,ktype)+jj-1)
    5          continue
            endif
   10    continue
         if(skip) goto 100
c
         call pseud1 (ntpse,nlp,zlp,clp,
     1                scr,ex,dmaxi2,maxp2,nx,ny,nz)
         call pseud2 (lmax,ntpse(2),nlp(ntpse(1)+1),zlp(ntpse(1)+1),
     1                clp(ntpse(1)+1),scr,ex,qq,dmaxi2,maxp2,nx,ny,nz)
         call smul(scr,scr,fpi,maxp2*dmaxi2)
c
c        move from the expanded basis function array into the
c        original one.
         call lpd1(nprim,ex,igbegn,igend,jgbegn,jgend,maxp2,
     $        sp,mini,maxi,minj,maxj,
     $        scr,dmini,dmaxi,dminj,dmaxj,nx,ny,nz,nderiv,
     $        iatom,jatom,katom)
         if (nderiv .gt.0) then
c
c        ----- transform the primitive derivative integrals, and
c              then sum into the correct place in the final array
            do 65 junk=1,6
               if (junk.le.3) then
                  coord=junk
                  atom=iatom
               else
                  coord=junk-3
                  atom=jatom
               end if
c     
               if (atom.eq.katom) go to 65
c     
               call trans1(sp(1,1,junk+1),conint,nprimi,nprimj,
     #              nconti,ncontj,acf,bcf,tmp1,len1,lenblk,imin,
     #              imax,jmin,jmax,nocart)
c     
c     ----- put this angular momentum block in the right places
c     
               numi=nobf(itype)
               numj=nobf(jtype)
               istrt=start(iatom,itype)
               jstrt=start(jatom,jtype)
               intgrl=0
c     
               do 64 if=1,numi
                  do 63 jf=1,numj
                     intgrl=intgrl+1
                     do 62 jc=1,ncontj
                        j=jstrt+(jc-1)*numj+jf
                        do 61 ic=1,nconti
                           i=istrt+(ic-1)*numi+if
                           if (i.ge.j) then
                              ij=i*(i-1)/2+j
                           else
                              if (iatom.eq.jatom.and.itype.eq.jtype)
     #                                   go to 61
                              ij=j*(j-1)/2+i
                           end if
c     
                           ds(ij,coord,atom)=ds(ij,coord,atom)+
     #                          conint(ic,jc,intgrl)
                           ds(ij,coord,katom)=ds(ij,coord,katom)-
     #                          conint(ic,jc,intgrl)
 61                     continue
 62                  continue
 63               continue
 64            continue
 65         continue
         endif
 100  continue
c     
c
      return
      end
