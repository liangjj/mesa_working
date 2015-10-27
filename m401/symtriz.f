*deck @(#)symtriz.f	5.1 11/6/94
      subroutine symtriz(z,a,maxcor,ops,prnt,nbf,nnp,s,smhalf,c,
     $                   bftoso,numso,salc,vecrep,nirrep,
     $                   lambda,lirrep)
c***begin prologue     symtriz.f
c***date written       940119  
c***revision date      11/6/94      
c
c***keywords           symmetry, orbitals, projection
c***author             martin, richard(lanl)
c***source             @(#)symtriz.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       symtriz.f
      implicit none
c     --- input variables -----
      integer nbf,nnp,nirrep,maxcor
      logical prnt
c     --- input arrays (unmodified) ---
      character*4096 ops
      character*(*) lirrep(nirrep)
      integer lambda(nirrep)
      real*8 s(nnp),smhalf(nnp)
c     --- input arrays (scratch) ---
      integer a(maxcor)
      real*8 z(*)
      real*8 salc(nbf,nbf)
c     --- output arrays ---
      integer vecrep(nirrep),numso(nirrep)
      real*8 bftoso(nbf,nbf),c(nbf,nbf)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer u,usym,eigval,t1,t2,t3,t4,t5
      integer norm,triang
      integer top,wpadti,iprint
      integer nmo,i
c
      logical debug
      parameter (debug=.false.)
c
      common /io/ inp,iout
c
 1000 format(1x,'raw salcs')
 1010 format(1x,'salc overlap')
 1020 format(1x,'orthogonalized salcs')
 1025 format(1x,'projected vectors')
 1030 format(5x,'project symmetry vectors')
 1040 format(5x,'orbital  index:',20i4,(/20x,20i4))
 1050 format(5x,'symmetry index:',20i4,(/20x,20i4))
 1060 format(5x,'symmetry label:',20(1x,a3),(/20x,20(1x,a3)))
 1070 format(5x,'calculation performed in reduced basis')
c
c
      write (iout,1030)
c
c     --- allocate core ---
      u=1
      usym=u+nbf*nbf
      eigval=usym+nbf*nbf
      t1=eigval+nbf
      t2=t1+nbf*nbf
      t3=t2+nbf*nbf
      t4=t3+nbf*nbf
      t5=t4+nbf
      norm=t5+nbf*nbf
      triang=norm+nbf*nirrep
      top=wpadti(triang+nnp)
c
      if(top.gt.maxcor) then
         call lnkerr('need more core for symtriz in m401')
      endif
c
      if(debug) then
         write(iout,1000)
         call matout(salc,nbf,nbf,nbf,nbf,iout)
      end if
c
c     --- symmetrically orthogonalize the basis function to salc matrix ---
      call trtosq(z(t1),s,nbf,nnp)
      call ebc(z(t2),z(t1),salc,nbf,nbf,nbf)
      call ebtc(z(t1),salc,z(t2),nbf,nbf,nbf)
      if(debug) then
         write(iout,1010)
         call matout(z(t1),nbf,nbf,nbf,nbf,iout)
      end if
      call sqtotr(z(usym),z(t1),nbf,nnp)
c
c     --- usym now holds the salc x salc overlap matrix ---
      iprint=0
      call sinv(z(usym),smhalf,z(u),z(eigval),z(t1),z(t2),
     $          nbf,nnp,z(triang),iprint)
      call trtosq(z(u),smhalf,nbf,nnp)
      call ebc(bftoso,salc,z(u),nbf,nbf,nbf)
c
      if(debug) then
         write(iout,1020)
         call matout(bftoso,nbf,nbf,nbf,nbf,iout)
      end if
c
c     --- project the vectors, guess goes in c,
c         good stuff comes back in usym ---
      call trtosq(z(t5),s,nbf,nnp)
      nmo=nbf
      call symprj(nirrep,nbf,nbf,numso,lambda,lirrep,c,nmo,bftoso,
     $            z(usym),vecrep,z(t5),z(t1),z(t2),z(t3),z(t4),
     $            z(norm),ops)
      call vmove(c,z(usym),nbf*nbf)
c
c
      if(prnt) then
         write (iout,1025)
         call matout(c,nbf,nbf,nbf,nbf,iout)
      end if
      if (prnt) then
         write(iout,1040) (i,i=1,nbf)
         write(iout,1050) (vecrep(i),i=1,nbf)
         write(iout,1060) (lirrep(vecrep(i)),i=1,nbf)
      end if
c
c
      return
      end
