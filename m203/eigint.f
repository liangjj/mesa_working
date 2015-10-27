*deck @(#)eigint.f	1.1  11/30/90
      subroutine eigint(eigvec,fre,natoms3,nvar,vname,b,evcint,
     #                  zpe,icfx)
c
      implicit integer (a-z)
      character*(*) vname(nvar)
      real*8 eigvec(natoms3,natoms3),evcint(nvar,nvar)
      real*8 b(natoms3,nvar),fre(natoms3)
      real*8 fac,zpe
      common /io/ inp,iout
c
      write(iout,1099)
 1099 format(/'    cartesian normal modes '/)
      call frqprt(eigvec,natoms3,natoms3,natoms3,natoms3,0,0,
     #            vname,collab,0,fre,.true.)
c
      if(icfx.eq.1) return
c
c     transform the eigenvectors to internal coordinates
c     first crunch the columns of the eigenvector matrix
c
      nrt=6
      count=0
      do 10 icol=1,natoms3
         if(abs(fre(icol)).gt.5.0) then
c
            count=count+1
            fre(count)=fre(icol)
            do 5 j=1,natoms3
               eigvec(j,count)=eigvec(j,icol)
    5       continue
         endif
   10 continue
      if(count.ne.natoms3-nrt) then
c         call lnkerr(' wrong number of zero frequencies ')
      write(iout,2250)
 2250 format(/
     #,' ** warning:  different # of zero frequencies than expected'/)
      endif
      call ebtc(evcint,b,eigvec,nvar,natoms3,nvar)
c
c     normalize the normal modes expressed in z-matrix coodinates
      do 200 j=1,nvar
         fac=0.0
         do 100 i=1,nvar
            fac=fac+evcint(i,j)*evcint(i,j)
  100    continue
         do 150 i=1,nvar
            evcint(i,j)=evcint(i,j)/sqrt(fac)
  150    continue
  200 continue
c
      write(iout,1098)
 1098 format(/'  internal coordinate normal modes '/)
      call frqprt(evcint,nvar,nvar,nvar,nvar,1,0,vname,collab,
     $            0,fre,.true.)
      write(iout,350) zpe
  350 format(/' ** zero point vibrational energy: ',f5.2,' kcal/mole')
      return
      end
