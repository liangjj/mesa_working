*deck @(#)symtriz.f	5.1 11/6/94
      subroutine symtriz(ops,prnt,lirrep,symlabl,bflabl,nbf,nnp,
     $                   maxrep,maxnbf,s,smhalf,c)
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
      integer nbf,nnp,maxrep,maxnbf,oldnbf,nirrep,numsn,numso
      integer salc,salct,bftoso,lambda,vecrep,bfns,indx
      logical prnt
c     --- input arrays (unmodified) ---
      character*4096 ops
      character*4 subgrp
      real*8 s(nnp),smhalf(nnp)
c     --- input arrays (scratch) ---
      real*8 c(nbf,nbf)
      character*(*) lirrep(maxrep), symlabl(maxnbf), bflabl(maxnbf)
      character*32 xform
      character*8 chrkey, scatyp
c     --- dynamically allocated or scratch ---
      integer a, ngot, idum
      real*8 z
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer u,usym,eigval,t1,t2,t3,t4,t5
      integer norm,triang
      integer top,wpadti,iadtwp,iprint
      integer nmo,i
c
      logical logkey,debug,scat,drop
      parameter (debug=.false.)
c
      common /io/ inp,iout
      pointer (p,z(1)), (p,a(1))
c
 1000 format(1x,'raw salcs')
 1010 format(1x,'salc overlap')
 1020 format(1x,'orthogonalized salcs')
 1025 format(1x,'projected vectors')
 1030 format(5x,'project symmetry vectors')
 1040 format(5x,'orbital  index:',20i4,(/20x,20i4))
 1050 format(5x,'symmetry index:',20i4,(/20x,20i4))
 1060 format(5x,'symmetry label:',20(1x,a3),(/20x,20(1x,a3)))
 1070 format(5x,'calculation performed in reduced basis: '/,15x,
     $          'old basis size = ',i5,1x,'new basis size = ',i5)
c
c
      xform='"guess vector"'
      scatyp=chrkey(ops,'scattering','none',' ')
      if(scatyp(1:4).eq.'kohn'.or.scatyp(1:8).eq.'r-matrix') then
         scat=.true.
      end if
      drop=logkey(ops,'drop',.false.,' ')
c----------------------------------------------------------------------c
c         drop is a flag which made m330 drop function components      c
c         oldnbf is the original basis set with no functions deleted   c
c                nbf is the size of the deleted set                    c
c----------------------------------------------------------------------c
      if (drop) then
         call iosys('read integer "old nbf" from rwf',1,oldnbf,0,' ')
         write(iout,1070) oldnbf, nbf
      else
         oldnbf=nbf
      endif   
      write (iout,1030)
c
c        --- get the dimensions, symmetry information, etc ---
c
      call iosys('read character "group symbol" from rwf',0,0,0,                  
     $            subgrp)                                                         
      call grpsym(subgrp,lirrep,maxrep) 
      call iosys('read integer "number of irreducible representations"'
     $          //' from rwf',1,nirrep,0,' ')
c      call iosys('read character "labels of irreducible '//
c     $           'representations" from rwf',nirrep*len(lirrep(1)),
c     $            0,0,lirrep)
c
c     --- allocate core ---
c
      numso=1
      numsn=numso+maxrep
      u=iadtwp(numsn+maxrep)
      usym=u+nbf*nbf
      eigval=usym+nbf*nbf     
      salc=eigval+nbf
      salct=salc
      if(drop) then
         salct=salc+nbf*nbf
      endif
      bftoso=salct+oldnbf*oldnbf
      lambda=wpadti(bftoso+nbf*nbf)
      vecrep=lambda+maxrep
      t1=iadtwp(vecrep+nbf)
      t2=t1+nbf*nbf
      t3=t2+nbf*nbf
      t4=t3+nbf*nbf
      t5=t4+nbf
      norm=t5+nbf*nbf
      triang=norm+nbf*nirrep
      bfns=wpadti(triang+nnp)
      indx=bfns
      if (drop) then
          indx=bfns+nbf
      endif
      top=indx+oldnbf  
      call getmem(top,p,ngot,'symtriz',0)
c
      call iosys('read integer "number of symmetry orbitals" from rwf',
     $           nirrep,a(numso),0,' ')
      call icopy(a(numso),a(numsn),nirrep)
      call iosys('read integer "degeneracies of irreducible'//
     $           ' representations" from rwf',nirrep,a(lambda),0,' ')
c
      call iosys('read real "salc transformation matrix" from rwf',
     $            -1,z(salct),0,' ')
      if (drop) then
          call nsalc(z(salct),z(salc),a(indx),nirrep,a(numso),a(numsn),
     1               oldnbf,nbf,prnt)
      endif
      if (scat) then
         call aosym(z(salc),a(numsn),a(bfns),nirrep,nbf)
      endif
      if(debug) then
         write(iout,1000)
         call matout(z(salc),nbf,nbf,nbf,nbf,iout)
      end if
c
c     --- symmetrically orthogonalize the basis function to salc matrix ---
      call trtosq(z(t1),s,nbf,nnp)
      call ebc(z(t2),z(t1),z(salc),nbf,nbf,nbf)
      call ebtc(z(t1),z(salc),z(t2),nbf,nbf,nbf)
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
      call ebc(z(bftoso),z(salc),z(u),nbf,nbf,nbf)
c
      if(debug) then
         write(iout,1020)
         call matout(z(bftoso),nbf,nbf,nbf,nbf,iout)
      end if
      call iosys('write real "orthogonal salc transformation matrix"'
     $          //' to rwf',nbf*nbf,z(bftoso),0,' ')
c
c     --- project the vectors, guess is in c,
c         good stuff comes back in usym ---
      call trtosq(z(t5),s,nbf,nnp)
      nmo=nbf
      call symprj(nirrep,nbf,nbf,a(numsn),a(lambda),lirrep,c,nmo,
     $            z(bftoso),z(usym),a(vecrep),z(t5),z(t1),z(t2),
     $            z(t3),z(t4),z(norm),ops)
      if (drop) then
c----------------------------------------------------------------------c
c         overwrite original salc matrix with new salc matrix          c
c----------------------------------------------------------------------c
         call iosys('write integer "number of symmetry orbitals" to '
     1            //'rwf',nirrep,a(numsn),0,' ')
         call iosys('write real "salc transformation matrix" to rwf',
     1               nbf*nbf,z(salc),0,' ')
      end if
c
c     save the projected vectors on the rwf, and store their
c     associated symmetry indices.
      call iosys('write real '//xform//' to rwf',nbf*nbf,z(usym),
     $            0,' ')
      call iosys('write integer "guess vector symmetries" '
     $           //'to rwf',nbf,a(vecrep),0,' ')
      call vmove(c,z(usym),nbf*nbf)
c
c
      if(prnt) then
         write (iout,1025)
         call matout(c,nbf,nbf,nbf,nbf,iout)
      end if
      if (prnt) then
         write(iout,1040) (i,i=1,nbf)
         write(iout,1050) (a(vecrep+i-1),i=1,nbf)
         write(iout,1060) (lirrep(a(vecrep+i-1)),i=1,nbf)
      end if
      if (logkey(ops,'print=guess=salc',.false.,' ')) then
          write(iout,*)'Symmetrized guess vector'
          do 10 i=1,nbf
             symlabl(i)=lirrep(a(vecrep+(i-1)))
 10       continue 
          call matprt(c,nbf,nbf,nbf,nbf,1,2,bflabl,
     $                symlabl,0,z(eigval),.true.)
      endif
c
c     used as a test for unit matrix
c      
c      call trtosq(z(t1),s,nbf,nnp)
c      call ebc(z(t2),z(t1),c,nbf,nbf,nbf)
c      call ebtc(z(t1),c,z(t2),nbf,nbf,nbf)
c      call matout(z(t1),nbf,nbf,nbf,nbf,iout)
      call getmem(-ngot,p,idum,'symtriz',idum)
c
c
      return
      end
