*deck @(#)lpints.f	5.1  11/6/94
      subroutine lpints(iatom,jatom,imax,jmax,nprimi,nprimj,
     $                   nprim,ncont,nat,nbtype,ntypes,maxp2,maxi2,
     $                   i1,j1,c,ex,noprim,ptprim,ptcont,cont,
     $                   maxmom,mini,maxi,minj,maxj,nx,ny,nz,sp,
     $                   ntpse,nlp,zlp,clp,qq)
c***begin prologue     lpints.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           effective core potentials
c***author             martin, richard(lanl)
c***source             @(#)lpints.f	5.1   11/6/94
c***purpose            driver rotuine for ecp integrals 
c***description
c     controls computation of pseudopotential integrals.
c     
c    
c
c***references
c         l.e. mcmurchie and e.r. davidson, j. comp. phys.,44,289(1981).
c         r.l. martin, to be published.
c***routines called
c
c***end prologue       lpints.f
      implicit none
c     --- input variables -----
      integer iatom,jatom,imax,jmax,nprimi,nprimj
      integer nprim,ncont,nat,nbtype,ntypes
      integer maxp2,maxi2,i1,j1
      integer mini,maxi,minj,maxj
c     --- input arrays (unmodified) ---
      integer noprim(nat,ntypes)
      integer ptprim(nat,ntypes),ptcont(nat,ntypes)
      integer maxmom(ntypes),nx(*),ny(*),nz(*)
      real*8 ex(nprim),c(3,nat),cont(ncont)
c     --- input arrays (scratch) ---
      integer nlp(*)
      integer ntpse(*)
      real*8 zlp(*),clp(*)
      real*8 qq(*)
c     --- output arrays ---
      real*8 sp(*)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer igbegn,igend,jgbegn,jgend
      integer istart,jstart,iend,jend
      integer lamax,lbmax
      integer katom,kk,nt,ktype,jj,lmax
      logical skip
      real*8 xa,ya,za,xb,yb,zb,xc,yc,zc
      real*8 cax,cay,caz,ca,ca2,cbx,cby,cbz,cb,cb2
      real*8 pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
      real*8 zero,one,two,three,four,five,six,ten
      real*8 half
c
      parameter (half=5.0d-01)
c
      common/prims/igbegn,igend,jgbegn,jgend
      common/limit/istart,jstart,iend,jend
      common/angmax/lamax,lbmax
      common/center/xa,ya,za,xb,yb,zb,xc,yc,zc
      common/dist/cax,cay,caz,ca,ca2,cbx,cby,cbz,cb,cb2
      common/pifac/pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
      common/const/zero,one,two,three,four,five,six,ten
c
c     --- initialize the output integral array
      call rzero(sp,maxp2*maxi2)
      lamax=imax+1
      lbmax=jmax+1
      istart=mini
      iend=maxi
      jstart=minj
      jend=maxj
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
      write(6,*) 'new basis function shell pair'
      write(6,*) 'lamax,lbmax,istart,iend,jstart,jend,igbegn,'
      write(6,*) 'igbegn,igend,jgbgn,jgend,nprimi,nprimj'
      write(6,*) lamax,lbmax,istart,iend,jstart,jend,igbegn,igend,
     $           jgbegn,jgend,nprimi,nprimj
    
c
      write(6,*) 'now looping over centers trying to find an ecp'
      do 100 katom=1,nat
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
         write(6,*) 'katom,xc,yc,zc',katom,xc,yc,zc
         nt=0
         kk=0
         lmax=0
         skip=.true.
         write(6,*) 'looping over basis function types',nbtype+1,'to',ntypes
         write(6,*) 'if the center has an ecp, there will be info here'
         do 10 ktype=nbtype+1,ntypes
            write(6,*) 'ktype',ktype
            if(noprim(katom,ktype).ne.0) then
               skip=.false.
               lmax=max(lmax,maxmom(ktype))
               nt=nt+1
               ntpse(nt)=noprim(katom,ktype)
               write(6,*) 'lmax,nt,ntpse(nt)',lmax,nt,ntpse(nt)
               do 5 jj=1,noprim(katom,ktype)
                  kk=kk+1
                  nlp(kk)=int(cont(ptcont(katom,ktype)
     $                             +noprim(katom,ktype)+jj-1)+half)
                  clp(kk)=cont(ptcont(katom,ktype)+jj-1)
                  zlp(kk)=ex(ptprim(katom,ktype)+jj-1)
                  write(6,*) 'jj,nlp,clp,zlp',jj,nlp(kk),clp(kk),zlp(kk)
    5          continue
            endif
   10    continue
         if(skip) goto 100
c
         call pseud1 (ntpse,nlp,zlp,clp,
     1                sp,ex,maxi2,maxp2,nx,ny,nz)
         call pseud2 (lmax,ntpse(2),nlp(ntpse(1)+1),zlp(ntpse(1)+1),
     1                clp(ntpse(1)+1),sp,ex,qq,maxi2,maxp2,nx,ny,nz)
c
  100 continue
      call smul(sp,sp,fpi,maxi2*maxp2)
c
c
      return
      end
