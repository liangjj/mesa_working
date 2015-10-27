*deck @(#)projct.f	5.1  11/6/94
      subroutine projct(gampt,atprmt,gamma,tr,
     #     p,coeffs,nirrep,natoms,nop,lengam,
     #     nfunc,nlamda,iatom,angmom,nrep,symat,nsymat,
     #     irrep,relatm,dump)
c
      implicit integer (a-z)
c
      integer gampt,atprmt(natoms,nop),symat(natoms),relatm(nsymat)
      logical dump
      real*8 gamma(lengam,nop)
      real*8 tr(nfunc,nfunc,nop)
      real*8 p(nfunc*nsymat,nfunc*nsymat,nlamda)
      real*8 coeffs(nfunc*nsymat,nlamda,nrep)
      real*8 ovrlap,norm,sdot,small
c
      parameter (small=1.0d-04)
c
      common /io/     inp,iout
 1000 format(1x,'  raw p-matrices:')
 1010 format(1x,'  orthonormalized p-matrices:')
 1020 format(/,1x,i5,' vectors of ',i2,'th irrep, degeneracy ',i1,
     $       /,' atoms:',(t10,20i3))
c
c     ----- zero space to accumulate coefficients -----
c
c     nfunc is the number of functions in the shell type under consideration.
c     nsymat is the number of atoms which go into one another under one
c     of the operations of the point group. nsf is the dimension of
c     the resulting salc.

      nsf=nfunc*nsymat
      call rzero(p,nsf**2*nlamda)
      call rzero(coeffs,nsf*nlamda*nrep)
c
      pt=gampt
      do 200 l=1,nlamda
         do 100 op=1,nop
            do 90 atom=1,nsymat
               iatom=relatm(atom)
               jatom=atprmt(iatom,op)
c              find the relative index of jatom in the related atom set.
               do 80 ptatom=1,nsymat
                  if(relatm(ptatom).eq.jatom) go to 81
   80          continue
               call lnkerr('error in symmetry related atom list')
   81          continue
c              apply the projection operator. note that we use the
c              transpose of the basis function transformation matrix.
               do 2 j=1,nfunc
                  do 1 i=1,nfunc
                     indi=(atom-1)*nfunc
                     indj=(ptatom-1)*nfunc
                     p(indi+i,indj+j,l)=p(indi+i,indj+j,l)
     #                                  +tr(j,i,op)*gamma(pt,op)
 1                continue
 2             continue
  90        continue
 100     continue
c
c     ----- step to the next diagonal in the representation matrix
c
         pt=pt+nlamda+1
 200  continue
c
c
      if(dump) then
         do 206 l=1,nlamda
            write(iout,1000)
            call matout(p(1,1,l),nsf,nsf,nsf,nsf,iout)
 206     continue
      endif
c
c     orthonormalize the p-matrix.
      do 212 l=1,nlamda
         do 210 vec=1,nsf
            do 208 old=1,vec-1
               ovrlap=sdot(nsf,p(1,old,l),1,p(1,vec,l),1)
               call saxpy(nsf,-ovrlap,p(1,old,l),1,p(1,vec,l),1)
 208        continue
            norm=sqrt(sdot(nsf,p(1,vec,l),1,p(1,vec,l),1))
            if(norm.gt.small) then
               call sscal(nsf,1.0d+00/norm,p(1,vec,l),1)
            else
               call rzero(p(1,vec,l),nsf)
            endif
 210     continue
 212  continue
c
c
      if(dump) then
         write(iout,1010)
         do 214 l=1,nlamda
            call matout(p(1,1,l),nsf,nsf,nsf,nsf,iout)
 214     continue
      endif
c
c     ----- find the non-negligable combinations and transfer to
c     coefficient array
c
      do 300 n=1,nrep
         do 250 l=1,nlamda
            do 220 i=1,nsf
               norm=sdot(nsf,p(1,i,l),1,p(1,i,l),1)
               if (norm.ne.0.0d+00) then
                  norm=sqrt(abs(norm))
               end if
               if (norm.gt.small) go to 230
 220        continue
c
            call lnkerr('not enough salcs found')
c
 230        continue
            call vmove(coeffs(1,l,n),p(1,i,l),nsf)
            call rzero(p(1,i,l),nsf)
 250     continue
 300  continue
c
      junk=nsf*nlamda
c
c
      if(dump) then
         write(iout,1020) nrep,irrep,nlamda,(relatm(i),i=1,nsymat)
         call matout(coeffs,nsf,nlamda*nrep,nsf,nlamda*nrep,iout)
      endif
c
c     ----- check that we did indeed use all salcs formed -----
c
      do 400 l=1,nlamda
         do 390 j=1,nsf
            do 380 i=1,nsf
               if (abs(p(i,j,l)).gt.small) then
                  call lnkerr('remnants left in p-matrix')
               end if
 380        continue
 390     continue
 400  continue
c
c
      return
      end
