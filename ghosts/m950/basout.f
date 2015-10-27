*deck basout
      subroutine basout(coords,ex,cont,ptprim,noprim,nocont,ptcont,
     #           nat,nprim,maxcor,ntypes,nbtype,nnp,ncont,
     #           start,nbasis,zan,nocart,nobf,maxmom,mintyp,
     #   nx,ny,nz,minmom,igeobas,idebug,iprt,index,idrop,oldnbf)
c
c
      implicit integer (a-z)
c
      character*1 itoc
      character*4 funcnm
      logical logkey
      real*8 ex(nprim),cont(ncont)
      real*8 zan(nat)
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(ntypes)
      integer nobf(ntypes),maxmom(ntypes),mintyp(ntypes),minmom(ntypes)
      integer nx(*),ny(*),nz(*)
      integer index(oldnbf)
c
c  dimensions for the arrays being built for quad codes
c
      real*8 coords(3*nat),eta(1000,5)
      integer nstart(1000), nstop(1000)
      integer lnew(1000),mnew(1000),nnew(1000)
c
c
      common/io/inp,iout
c
      nstart(1) = 1
      ikount = 0
      iloc = 0
c
c.. bhl 9/13/89         using index array to drop functions
c
      if(idrop.ne.0) then
      write(iout,*)' reading the drop index vector from rwf'
      call iosys('read integer "packing index vector" from kohn',
     $ oldnbf,index,0,' ')
      else
      do 999 i=1,oldnbf
      index(i)=i
 999  continue
      end if
c 
      kbf=0
c
c.. bhl 9/13/89             kbf counter
c
      do 9 iatom=1,nat
c        write(iout,*)' '
c        write(iout,*)' atom number ',iatom
c        write(iout,*)' charge      ',zan(iatom)
             do 7 itype=1,nbtype
               if (noprim(iatom,itype).le.0) go to 7
c
                   mini=mintyp(itype)
                   maxi=mintyp(itype)+nocart(itype)-1
                   imax=maxmom(itype)
c        write(iout,*)'  primitive type s=1 p=2 d=3 etc. ',itype
c        write(iout,*)'  maximum momemtum ',imax
c
                   nprimi=noprim(iatom,itype)
                   nconti=nocont(iatom,itype)
c        write(iout,*)'    number of primitives functions ',nprimi
c        write(iout,*)'    number of contracted functions ',nconti
c
c  build arrays for quad codes: nstart, nstop, lnew, mnew, nnew, eta
c
c
c  the way the  indexing for quad codes works is that there are
c  ncontra contracted functions.  for the ith contracted function nstart(i)
c  and nstop(i) give the first and last locations of the relevant quantities
c  in the arrays lnew(i), mnew(i), nnew(i), and eta(i,j=1,5)
c   lnew, mnew, and new are the powers of x, y, and z of the primitive
c   and eta(i,j=1,5) contains the center, exponent and contraction coefficient
c
      do 6 icont=1,nconti
       do 3 m=mini,maxi
c..bhl 9/13/89                  kbf counter 
      kbf=kbf+1
      if(index(kbf).ne.0) then
c..bhl 9/13/89                  if statement
      ikount = ikount + 1
      do 4 iii=1,nprimi
      iiirel = iii-1+ (icont-1)*nprimi
      if (cont(ptcont(iatom,itype)+iiirel).gt.1.e-20) then
       iloc = iloc+1
       eta(iloc,1) = coords(3*(iatom-1)+1)
       eta(iloc,2) = coords(3*(iatom-1)+2)
       eta(iloc,3) = coords(3*(iatom-1)+3)
       eta(iloc,4) = ex(ptprim(iatom,itype)+iii - 1)
       eta(iloc,5) = cont(ptcont(iatom,itype)+iiirel)
       lnew(iloc) = nx(m)
       mnew(iloc) = ny(m)
       nnew(iloc) = nz(m)
      endif
   4  continue
      nstop(ikount) = iloc
      nstart(ikount+1) = nstop(ikount)+ 1
      end if
c..bhl 9/13/89                  end if
  3    continue
   6  continue
c
c        write(iout,*)' '
c        write(iout,*)'    primitive exponents '
      if(iprt.ne.0) then
        call matout(ex(ptprim(iatom,itype)),nprimi,1,
     $           nprimi,1,iout)
      end if
c        write(iout,*)' '
c        write(iout,*)'    contraction coefficients '
      if(iprt.ne.0) then
        call matout(cont(ptcont(iatom,itype)),nprimi,nconti,
     $           nprimi,nconti,iout)
      end if
c        write(iout,*)' '
c        write(iout,*)'    xyx convention '
        do 5 m=mini,maxi
        if(iprt.ne.0)write(iout,*)nx(m),ny(m),nz(m)
    5   continue
    7       continue
    9 continue
      neta = iloc
      ncontra = ikount
c
c write the newly indexed basis set information to binary file igeobas
c
c
      write(igeobas) ncontra
      write(igeobas) (nstart(i),nstop(i),i=1,ncontra)
      write(iout,*)' '
      write(iout,*)' geometry and basis written to file ',igeobas
      if(idebug.ne.0) then
      write(iout,*)' nbf ',ncontra
      do 12 i=1,ncontra
      nlow = nstart(i)
      nhi = nstop(i)
      do 11 j=nlow,nhi
      write(igeobas) lnew(j),mnew(j),nnew(j),(eta(j,jj),jj=1,5)
      if(iprt.ne.0) then
      write(iout,*) lnew(j),mnew(j),nnew(j),(eta(j,jj),jj=1,5)
      end if
  11   continue
  12   continue
      endif
c
c
c
      return
      end
