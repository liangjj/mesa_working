*deck @(#)pm1031.f	5.1  11/6/94
      subroutine pm1031(z,a)
c***begin prologue     m1031
c***date written       861211  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           ci lagrangian, ci gradients, ci derivatives
c***author             saxe, paul (lanl)
c***source             @(#)pm1031.f	5.1   11/6/94
c***purpose            form the ci lagrangian.
c***description
c
c***references         "generalization of analytic configuration
c   interaction (ci) gradient techniques for potential energy
c   hypersurfaces, including a solution to the coupled perturbed
c   hartree-fock equations for multiconfiguration scf molecular
c   wave functions", yoshihiro osamura, yukio yamaguchi and henry f.
c   schaefer iii, j. chem. phys. 77(1), p.383-390, 1 july 1982.
c
c***routines called
c
c***end prologue       m1031
c
      implicit integer (a-z)
c
      real*8 z(*)
      character*4096 ops
      character*128 tints,moden
      integer a(*)
c
c     ----- lenval is the maximum size of the buffer to read mo integrals
c
      parameter (lenval=50000)
c
      common /io/ inp,iout
c
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- recover dimensions for core allocation -----
c
      call iosys('read integer "number of basis functions" from rwf',
     $     1,num,0,' ')
      nnp=num*(num+1)/2
c
c
c     ----- allocate core -----
c
      call getscm(0,a,maxcor,'?',0)
c
      lag=1
      t1=lag+num**2
      t2=t1+num**2
      ints=t2+num**2
      left=iadtwp(maxcor)-ints
      ntriang=min(left/2,nnp**2,lenval)/nnp
      if (ntriang.lt.1) call lnkerr('cant hold one triangle of mo ints')
      tpdm=ints+nnp*ntriang
      need=wpadti(tpdm+nnp*ntriang)
c
      call getscm(need,a,junk,' ',0)
c
c     ----- open the integrals unit -----
c
      call iosys('read character "transformed integral filename"'
     $         //' from rwf',0,0,0,tints)
      call iosys('open tints as old',0,0,0,tints)
c
c     open the mo density matrix unit
      call iosys('read character "mo density filename" from rwf',
     $           0,0,0,moden)
      call iosys('open moden as old',0,0,0,moden)
c
c     ----- and form the lagrangian -----
c
      call cilag(num,nnp,z(lag),z(ints),z(tpdm),ntriang,z(t1),z(t2),
     $     ops)
c
c
c     ----- and exit this link -----
c
      call iosys('close tints',0,0,0,' ')
      call iosys('close moden',0,0,0,' ')
      call chainx(0)
c
c
      stop
      end
