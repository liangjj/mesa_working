*deck @(#)vfmj.f	5.1  11/6/94
      subroutine vfmj(i2,dij,dkl,fij,fkl,index,npij,npkl,
     $               nfi,nfj,nfk,nfl,ijsh,klsh,
     $               lenv,angmom,imax,jmax,
     $               mmax,lmax,nroots,px,py,pz,lenxyz,mintyp,
     $               nocart,nbtype,nv,len)
c***begin prologue     vfmj.f
c***date written       931210  (yymmdd)  
c***revision date      11/6/94      
c
c   this code originated from routines vfmint in m312 and m712.
c
c   june 3, 1993       rlm at lanl
c     adding loops to handle higher angular momentum. this choked on
c     second derivatives of f functions.
c   11 november 1985   pws  at lanl
c     modifying original version from m312 to sum up energy by directly
c     multiplying primitive integrals by two-particle density matrix.
c***keywords           direct, j-matrix
c***author             martin, saxe, and lengsfield
c***source             @(#)vfmj.f	5.1   11/6/94
c***purpose            forms the coulomb matrix from the 1-dimensional
c                      components of the two-electron integrals and the
c                      density matrix.  
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       vfmj.f
      implicit none
c     --- input variables -----
      integer lenv,angmom(4),imax,jmax,mmax,lmax
      integer nroots,lenxyz,nbtype,nv,len
      integer npij,npkl,nfi,nfj,nfk,nfl
      integer px(lenxyz),py(lenxyz),pz(lenxyz)
      integer mintyp(nbtype),nocart(nbtype)
      real*8 i2(*)
      logical ijsh,klsh
c     --- input arrays (unmodified) ---
      integer index(len,6)
      real*8 dij(npij,nfi*nfj),dkl(npkl,nfk*nfl)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 fij(npij,nfi*nfj),fkl(npkl,nfk*nfl)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer nindex
      parameter (nindex=56)
      integer ipx(nindex),ipy(nindex),ipz(nindex)
      integer jpx(nindex),jpy(nindex),jpz(nindex)
      integer kpx(nindex),kpy(nindex),kpz(nindex)
      integer lpx(nindex),lpy(nindex),lpz(nindex)
      integer mini,maxi,minj,maxj,mink,maxk,minl,maxl
      integer junk,i,j,k,l,m
      integer ix,iy,iz,jx,jy,jz,kx,ky,kz,lx,ly,lz
      integer prim,nij,nkl
      integer inp,iout
      real*8 facij,fackl
      real*8 zero,one,two
      parameter (zero=0.0d+00,one=1.0d+00,two=2.0d+00)
c
c     real*8 i2(lenv,3,0:imax,0:jmax,0:mmax,0:lmax)
c
      common/io/inp,iout
c
c     --- set up minimum and maximum functions
c         in each momentum group ---
      mini=mintyp(angmom(1))
      maxi=mini+nocart(angmom(1))-1
      minj=mintyp(angmom(2))
      maxj=minj+nocart(angmom(2))-1
      mink=mintyp(angmom(3))
      maxk=mink+nocart(angmom(3))-1
      minl=mintyp(angmom(4))
      maxl=minl+nocart(angmom(4))-1
c
c     --- set up offsets ---
      junk=3*lenv
      do 1 i=mini,maxi
         ipx(i)=       junk*px(i)
         ipy(i)=  lenv+junk*py(i)
         ipz(i)=2*lenv+junk*pz(i)
    1 continue
      junk=junk*(imax+1)
      do 2 j=minj,maxj
         jpx(j)=junk*px(j)
         jpy(j)=junk*py(j)
         jpz(j)=junk*pz(j)
    2 continue
      junk=junk*(jmax+1)
      do 3 k=mink,maxk
         kpx(k)=junk*px(k)
         kpy(k)=junk*py(k)
         kpz(k)=junk*pz(k)
    3 continue
      junk=junk*(mmax+1)
      do 4 l=minl,maxl
         lpx(l)=junk*px(l)
         lpy(l)=junk*py(l)
         lpz(l)=junk*pz(l)
    4 continue
c
c     --- loop over integrals in this block; multiply by the
c         appropriate density matrix elements and accumulate
c         the j matrix as we go.
c
c         note that the routine prefac has stripped out integrals which
c         are too small to bother with. the dimension nv therefore refers
c         only to the effective number of non-zero npi*npj*npk*npl. the
c         array index maps the contracted list back into the full matrix.
c            index(nv,1)=i;  index(nv,2)=j; index(nv,3)=k; index(nv,4)=l
c            index(nv,5)=ij; index(nv,6)=kl
c   
c         the outer loop has made sure that the only shell blocks which
c         make it here have ishell.ge.jshell, kshell.ge.lshell.
c         thus,the integral (ij,kl) contributes to two jmatrix entries:
c            j(ij)=(ij,kl)*2d(kl); j(kl)=(ij,kl)*2d(ij)
c         where the factor of two accounts for the kl=lk equivalence
c
c         note,however,that if ishell=jshell(ijsh.eq.true),
c         or kshell=lshell(klsh.eq.true.), we have the full square
c         block, and so we should limit the sums to i>=j, or k>=l to
c         avoid double counting. instead of this, we do the full unrestricted 
c         but multiply by one to account for the fact that we summed the 
c         entire block. 
c
c         kshell=lshell
c            j(ij)=(ij,kl)*d(kl); 
c         ishell=jshell
c            j(kl)=(ij,kl)*d(ij)
c
      facij=two
      fackl=two
      if(ijsh) then
         facij=one
      endif
      if(klsh) then
         fackl=one
      endif
c           
      do 400 i=mini,maxi
         ix=ipx(i)
         iy=ipy(i)
         iz=ipz(i)
         do 300 j=minj,maxj
            jx=ix+jpx(j)
            jy=iy+jpy(j)
            jz=iz+jpz(j)
            nij=(j-minj)*nfi+(i-mini+1)
c
            do 200 k=mink,maxk
               kx=jx+kpx(k)
               ky=jy+kpy(k)
               kz=jz+kpz(k)
               do 100 l=minl,maxl
                  lx=kx+lpx(l)
                  ly=ky+lpy(l)
                  lz=kz+lpz(l)
                  nkl=(l-minl)*nfk+(k-mink+1)
c
                  go to (5,6,7,8,9,10,11,12,13),nroots
                     write(iout,*) 'vfmj:512,nroots',nroots
                     call lnkerr('bad number of rys roots')
c
    5             continue
c
c  a better way might be to pack the dij and dkl arrays.
c  we could then fill fij,fkl and then scatter into the full matrix.
c  this might vectorize better.
                  do 50 prim=1,nv 
                     fij(index(prim,5),nij)=fij(index(prim,5),nij)
     $                            +i2(prim+lx)*i2(prim+ly)*
     $                             i2(prim+lz)*
     $                             fackl*dkl(index(prim,6),nkl)
                     fkl(index(prim,6),nkl)=fkl(index(prim,6),nkl)
     $                            +i2(prim+lx)*i2(prim+ly)*
     $                             i2(prim+lz)*
     $                             facij*dij(index(prim,5),nij)
   50             continue
                  go to 15
c
    6             continue
                  m=1
                  do 60 prim=1,nv
                     fij(index(prim,5),nij)=fij(index(prim,5),nij)
     $                              +(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                              i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz))*
     $                              fackl*dkl(index(prim,6),nkl)
                     fkl(index(prim,6),nkl)=fkl(index(prim,6),nkl)
     $                              +(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                              i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz))*
     $                              facij*dij(index(prim,5),nij)
                     m=m+2
   60             continue
                  go to 15
c
    7             continue
                  m=1
                  do 70 prim=1,nv
                     fij(index(prim,5),nij)=fij(index(prim,5),nij)
     $                             +(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                              i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     $                              i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz))*
     $                              fackl*dkl(index(prim,6),nkl)
                     fkl(index(prim,6),nkl)=fkl(index(prim,6),nkl)
     $                             +(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                              i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     $                              i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz))*
     $                              facij*dij(index(prim,5),nij)
                     m=m+3
   70             continue
                  go to 15
c
    8             continue
                  m=1
                  do 80 prim=1,nv
                     fij(index(prim,5),nij)=fij(index(prim,5),nij)
     $                             +(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                              i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     $                              i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz)+
     $                              i2(m+3+lx)*i2(m+3+ly)*i2(m+3+lz))*
     $                              fackl*dkl(index(prim,6),nkl)
                     fkl(index(prim,6),nkl)=fkl(index(prim,6),nkl)
     $                             +(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                              i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     $                              i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz)+
     $                              i2(m+3+lx)*i2(m+3+ly)*i2(m+3+lz))*
     $                              facij*dij(index(prim,5),nij)
                     m=m+4
   80             continue
                  go to 15
c
    9             continue
                  m=1
                  do 90 prim=1,nv
                     fij(index(prim,5),nij)=fij(index(prim,5),nij)
     $                             +(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                              i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     $                              i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz)+
     $                              i2(m+3+lx)*i2(m+3+ly)*i2(m+3+lz)+
     $                              i2(m+4+lx)*i2(m+4+ly)*i2(m+4+lz))*
     $                              fackl*dkl(index(prim,6),nkl)
                     fkl(index(prim,6),nkl)=fkl(index(prim,6),nkl)
     $                             +(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                              i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     $                              i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz)+
     $                              i2(m+3+lx)*i2(m+3+ly)*i2(m+3+lz)+
     $                              i2(m+4+lx)*i2(m+4+ly)*i2(m+4+lz))*
     $                              facij*dij(index(prim,5),nij)
                     m=m+5
   90             continue
                  go to 15
c
   10             continue
                  m=1
                  do 110 prim=1,nv
                     fij(index(prim,5),nij)=fij(index(prim,5),nij)
     $                             +(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                              i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     $                              i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz)+
     $                              i2(m+3+lx)*i2(m+3+ly)*i2(m+3+lz)+
     $                              i2(m+4+lx)*i2(m+4+ly)*i2(m+4+lz)+
     $                              i2(m+5+lx)*i2(m+5+ly)*i2(m+5+lz))*
     $                              fackl*dkl(index(prim,6),nkl)
                     fkl(index(prim,6),nkl)=fkl(index(prim,6),nkl)
     $                             +(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                              i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     $                              i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz)+
     $                              i2(m+3+lx)*i2(m+3+ly)*i2(m+3+lz)+
     $                              i2(m+4+lx)*i2(m+4+ly)*i2(m+4+lz)+
     $                              i2(m+5+lx)*i2(m+5+ly)*i2(m+5+lz))*
     $                              facij*dij(index(prim,5),nij)
                     m=m+6
  110             continue
                  go to 15
c
   11             continue
                  m=1
                  do 120 prim=1,nv
                     fij(index(prim,5),nij)=fij(index(prim,5),nij)
     $                             +(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                              i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     $                              i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz)+
     $                              i2(m+3+lx)*i2(m+3+ly)*i2(m+3+lz)+
     $                              i2(m+4+lx)*i2(m+4+ly)*i2(m+4+lz)+
     $                              i2(m+5+lx)*i2(m+5+ly)*i2(m+5+lz)+
     $                              i2(m+6+lx)*i2(m+6+ly)*i2(m+6+lz))*
     $                              fackl*dkl(index(prim,6),nkl)
                     fkl(index(prim,6),nkl)=fkl(index(prim,6),nkl)
     $                             +(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                              i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     $                              i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz)+
     $                              i2(m+3+lx)*i2(m+3+ly)*i2(m+3+lz)+
     $                              i2(m+4+lx)*i2(m+4+ly)*i2(m+4+lz)+
     $                              i2(m+5+lx)*i2(m+5+ly)*i2(m+5+lz)+
     $                              i2(m+6+lx)*i2(m+6+ly)*i2(m+6+lz))*
     $                              facij*dij(index(prim,5),nij)
                     m=m+7
  120             continue
                  go to 15
c
   12             continue
                  m=1
                  do 130 prim=1,nv
                     fij(index(prim,5),nij)=fij(index(prim,5),nij)
     $                             +(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                              i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     $                              i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz)+
     $                              i2(m+3+lx)*i2(m+3+ly)*i2(m+3+lz)+
     $                              i2(m+4+lx)*i2(m+4+ly)*i2(m+4+lz)+
     $                              i2(m+5+lx)*i2(m+5+ly)*i2(m+5+lz)+
     $                              i2(m+6+lx)*i2(m+6+ly)*i2(m+6+lz)+
     $                              i2(m+7+lx)*i2(m+7+ly)*i2(m+7+lz))*
     $                              fackl*dkl(index(prim,6),nkl)
                     fkl(index(prim,6),nkl)=fkl(index(prim,6),nkl)
     $                             +(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                              i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     $                              i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz)+
     $                              i2(m+3+lx)*i2(m+3+ly)*i2(m+3+lz)+
     $                              i2(m+4+lx)*i2(m+4+ly)*i2(m+4+lz)+
     $                              i2(m+5+lx)*i2(m+5+ly)*i2(m+5+lz)+
     $                              i2(m+6+lx)*i2(m+6+ly)*i2(m+6+lz)+
     $                              i2(m+7+lx)*i2(m+7+ly)*i2(m+7+lz))*
     $                              facij*dij(index(prim,5),nij)
                     m=m+8
  130             continue
                  go to 15
c
   13             continue
                  m=1
                  do 140 prim=1,nv
                     fij(index(prim,5),nij)=fij(index(prim,5),nij)
     $                             +(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                              i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     $                              i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz)+
     $                              i2(m+3+lx)*i2(m+3+ly)*i2(m+3+lz)+
     $                              i2(m+4+lx)*i2(m+4+ly)*i2(m+4+lz)+
     $                              i2(m+5+lx)*i2(m+5+ly)*i2(m+5+lz)+
     $                              i2(m+6+lx)*i2(m+6+ly)*i2(m+6+lz)+
     $                              i2(m+7+lx)*i2(m+7+ly)*i2(m+7+lz)+
     $                              i2(m+8+lx)*i2(m+8+ly)*i2(m+8+lz))*
     $                              fackl*dkl(index(prim,6),nkl)
                     fkl(index(prim,6),nkl)=fkl(index(prim,6),nkl)
     $                             +(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                              i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     $                              i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz)+
     $                              i2(m+3+lx)*i2(m+3+ly)*i2(m+3+lz)+
     $                              i2(m+4+lx)*i2(m+4+ly)*i2(m+4+lz)+
     $                              i2(m+5+lx)*i2(m+5+ly)*i2(m+5+lz)+
     $                              i2(m+6+lx)*i2(m+6+ly)*i2(m+6+lz)+
     $                              i2(m+7+lx)*i2(m+7+ly)*i2(m+7+lz)+
     $                              i2(m+8+lx)*i2(m+8+ly)*i2(m+8+lz))*
     $                              facij*dij(index(prim,5),nij)
                     m=m+9
  140             continue
c
   15             continue
  100          continue
  200       continue
  300    continue
  400 continue
c
c
      return
      end
