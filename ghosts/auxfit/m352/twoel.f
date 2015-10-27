*deck %W% %G%
      subroutine twoel(c,ex,z,iz,cont,s,ptprim,noprim,nocont,ptcont,
     $     nat,nprim,maxcor,ntypes,nbtype,nnp,ncont,
     $     start,nbasis,nocart,nobf,maxmom,mintyp,
     $     nx,ny,nz,minmom,ops,amaxm,fna0,tmpsiz,lval,
     $     mycart,yn0,mindx)
c
c***begin prologue     %M%
c***date written       840723  
c***revision date      %G%
c
c***keywords           
c***author             saxe, paul and martin, richard (lanl) 
c***source             %W% %G%
c***purpose            
c***description
c      module to form the auxiliary two-electron integrals over generally 
c      contracted gaussian basis sets.
c
c***references
c
c***routines called
c
c***end prologue       twoel.f
      implicit integer (a-z)
c
c     ----- input arrays(unmodified) -----
      character*(*) ops
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(ntypes)
      integer nobf(ntypes),maxmom(ntypes),mintyp(ntypes),minmom(ntypes)
      integer nx(*),ny(*),nz(*)
      real*8 c(3,nat),ex(nprim),cont(ncont)
c
c     ----- scratch working space -----
      real*8 z(maxcor),s(nnp)
      integer iz(*)
c
c     ----- input flags -----
c
c     ----- local variables -----
      real*8 pi,pi252,one,two,four
      logical debug
c
      integer fna0(0:amaxm),tmpsiz(0:amaxm)
      integer lval(0:amaxm,mindx,3),mycart(0:amaxm),yn0(0:amaxm)
c
      parameter (one=1.0d+00,two=2.0d+00,four=4.0d+00,debug=.false.)
      common/io/inp,iout
c
c fna0(m):  total number of integrals (a,0) to be
c generated, this is just sum from l=0..mmax of (l+1)(l+2)/2, the number 
c of unique integrals of each angular momentum l.
c
c tmpsiz(m):  how much temporary space are we going to need for the (a,0)
c     calculation?  This is just max m=1,mmax ( (mmax-m+1)(m+1)(m+2)/2)
c
c      data fna0/1,4,10,20,35,56,84,120,165/
c      data tmpsiz/1,3,6,12,20,30,45,63,84/
c
c
c     ----- form the auxiliary two-electron integrals -----
      pi=four*atan(one)
      pi252=two*pi*pi*(sqrt(pi))
      do 9 ia=1,nat
         do 8 ja=1,ia
            do 7 it=1,nbtype
               if (noprim(ia,it).le.0) go to 7
               if (ia.ne.ja) then
                  jtypmx=nbtype
               else
                  jtypmx=it
               end if
               do 6 jt=1,jtypmx
                  if (noprim(ja,jt).le.0) go to 6
c
                  imax=maxmom(it)
                  jmax=maxmom(jt)
                  if (imax .lt. jmax) then
                     iatom=ja
                     jatom=ia
                     itype=jt
                     jtype=it
                     imax=maxmom(itype)
                     jmax=maxmom(jtype)
                  else
                     iatom=ia
                     jatom=ja
                     itype=it
                     jtype=jt
                  endif
                  mmax=imax+jmax
                  nprimi=noprim(iatom,itype)
                  nprimj=noprim(jatom,jtype)
                  nconti=nocont(iatom,itype)
                  ncontj=nocont(jatom,jtype)
                  npint=nprimi*nprimj
                  lenblk=nocart(itype)*nocart(jtype)
c
c     ----- allocate core for temporary vectors, etc. -----
c     core needed in gena0int
c
                  t=1
                  f=t+npint
                  kappa=f+npint*(mmax+1)
                  alpha=kappa+npint
                  ainv=alpha+2*npint
                  prctr=ainv+5*npint
                  expon=prctr+3*npint
c                  prmint=1
c                  xyz=max(expon+npint,prmint+lenblk*npint)
                  a0out=ainv+5*npint
                  wma=a0out+npint*fna0(mmax)
                  tmparry=wma+npint*3
                  top1=tmparry+npint*tmpsiz(mmax)*3
c
c     core needed in genacint
c
c     someplace to put the integrals and the temp array.  Note that
c     the "t" array when used clobbers the a0out array, which is fine as 
c     long as genacint copies the a0out stuff into tmparry2, which it does. 
c     
c calculate how big tmparry2 needs to be.  This should eventually go into
c an array like tmpsiz does, but for now generate on-the-fly
c
                  tmpsiz2=0
                  lmn=mod(mmax,2)
                  lmx=mmax
                  do 666 lc=0,jmax
                     jnk=0
                     do 6666 la=lmn,lmx,2
                        jnk=jnk+mycart(la)*mycart(lc)
 6666                continue 
                     tmpsiz2=max(jnk,tmpsiz2)
 666              continue 
                  acout=a0out+npint*fna0(mmax)
                  tmparry2=acout+npint*lenblk
c
c     core needed by trans1
c
                  conint=acout+npint*lenblk
                  tmp1=conint+nconti*ncontj*lenblk
                  len1=nconti*nprimj                  
                  top2=max(tmp1+len1,tmparry2+npint*tmpsiz2*3)
c
                  if (top1.gt.maxcor.or.top2.gt.maxcor) then
                     call lnkerr('m352: overlp..not enough core for s')
                  end if
c
c     ----- form the (a,0) integrals
c
                  call gena0int(iatom,jatom,imax,jmax,nprimi,nprimj,
     $                 ptprim(iatom,itype),ptprim(jatom,jtype),
     $                 ex,z(alpha),c,z(ainv),z(prctr),z(kappa),z(f),
     $                 z(expon),z(wma),z(a0out),z(tmparry),
     $                 tmpsiz(mmax),z(t),npint,nprim,nat,
     $                 mmax,pi252,lval,mycart,yn0,amaxm,mindx)
                  if (debug) then
                     write(iout,*) "returned from gena0int"
                  endif
c     
c     ----- form the (a,c) integrals from the (a,0)'s
c
                  call genacint(iatom,jatom,imax,jmax,nprimi,nprimj,
     $                 z(alpha),z(ainv),z(a0out),z(acout),z(tmparry2),
     $                 tmpsiz2,z(t),npint,mmax,pi252,lval,mycart,yn0,
     $                 amaxm,mindx,mintyp(itype),mintyp(jtype),
     $                 nx,ny,nz)
c
c     ----- transform to contracted functions -----
c
                  call trans1(z(acout),z(conint),nprimi,nprimj,nconti,
     $                        ncontj,cont(ptcont(iatom,itype)),
     $                        cont(ptcont(jatom,jtype)),z(tmp1),len1,
     $                        lenblk,minmom(itype),maxmom(itype),
     $                        minmom(jtype),maxmom(jtype),nocart)
c
c     ----- transfer integrals to total array -----
c
                  call put1el(s,z(conint),start,iatom,jatom,itype,jtype,
     $                        nconti,ncontj,nnp,lenblk,nat,nbtype,nobf)
c
    6          continue
    7       continue
    8    continue
    9 continue
c
c
      return
      end
