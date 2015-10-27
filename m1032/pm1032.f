*deck @(#)pm1032.f	5.1  11/6/94
      subroutine pm1032(z,a)
c
c***begin prologue     m1032
c***date written       861211  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           ci lagrangian, ci gradients, ci derivatives
c***author             saxe, paul (lanl)
c***source             @(#)pm1032.f	5.1   11/6/94
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
c***end prologue       m1032
c
      implicit integer (a-z)
c
      real*8 z(*)
      character*4096 ops
      integer a(*)
c
c..bhl..unicos      common //   a
      common /io/ inp,iout
c
c..bhl..unicos      equivalence (a,z)
c
c     ----- open the read-write file -----
c
c..bhl..unicos      call drum
c
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- recover dimensions for core allocation -----
c
      call iosys('read integer "number of atoms" from rwf',1,nat,0,' ')
      call iosys('read integer "number of basis functions" from rwf',
     $     1,num,0,' ')
      call iosys('read integer "number of shells" from rwf',
     $     1,nshell,0,' ')
      nder=3*nat
c
      maxshl=1
      lag=iadtwp(maxshl+nshell)
      u=lag+num**2
      intd1e=u+num*num*nder
      d1e=intd1e+nder
      need=wpadti(d1e+nder)
c
      call getscm(need,a,junk,' ',0)
c
c     ----- and form the total ci derivative -----
c
      call cider(num,nder,nat,z(lag),z(u),z(intd1e),z(d1e),ops,
     $     a(maxshl),nshell)
c
c     ----- and exit this link -----
c
      call chainx(0)
c
c
      stop
      end
