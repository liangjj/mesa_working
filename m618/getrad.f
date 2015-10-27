*deck @(#)getrad.f	5.1 11/6/94
      subroutine getrad(c,ian,nat,xcnt,ycnt,zcnt,rad,
     $                 ncnt,prnt,toang)
c***begin prologue     getrad.f
c***date written       yymmdd
c***revision date      2/7/94
c   january 25, 1994   rlm at lanl
c      converting to atomic units throughout
c***keywords
c***author             tawa, greg(lanl)
c***source             @(#)getrad.f	5.1 11/6/94
c***purpose
c***description
c
c
c
c***references
c                      a.rashin and b. honig, j.phys.chem. 89,5588(1985).
c
c***routines called
c
c***end prologue       getrad.f
      implicit none
c     --- input variables -----
      integer nat,ncnt
      logical prnt
      real*8 toang
c     --- input arrays (unmodified) ---
      integer ian(nat)
      real*8 c(3,nat)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 xcnt(ncnt), ycnt(ncnt), zcnt(ncnt), rad(ncnt)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i
      integer mxatm
      integer inp,iout
      parameter (mxatm=103)
      real*8 radius(mxatm)
      real*8 zero
c
      parameter (zero=0.0d0)
c
      common/io/inp,iout
c
c     --- the atomic radii (Angstroms) ---
c         many of these come from rashin and honig
c            li+,na+,k+,rb+,cs+
c            mg2+,ca2+,sr2+,ba2+
c            al3+,sc3+,y3+,la3+,ce3/4+
c            zn2+,cd2+,hg2+
c            ga3+,in3+
c            f-,cl-,br-,i-
c     --- the others are from ...
c         Van der Waals radii, 
c            Hydrogen 1.375(explicit hydrogen attached to Carbon)
      data radius/1.375d0,0.0d0,
     $            1.316d0,2*0.0d0,1.80d0,1.75d0,1.65d0,1.423d0,0.0d0,
     $            1.680d0,1.455d0,1.338d0,3*0.0d0,1.937d0,0.0d0,
     $            2.172d0,1.862d0,
     $                    1.541d0,8*0.0d0,1.338d0,
     $                            1.338d0,3*0.0d0,2.087d0,0.0d0,
     $            2.311d0,2.054d0,
     $                    1.733d0,8*0.0d0,1.509d0,
     $                            1.605d0,3*0.0d0,2.343d0,0.0d0,
     $            2.514d0,2.119d0,
     $                    1.808d0,
     $                            1.761d0,13*0.0d0,
     $                            8*0.0d0,1.541d0,
     $            23*0.d0/
c
c         Van der Waals radii, Hydrogen 1.000(amino hydrogen), 
c                              Nitrogen 1.85 (sp3, 4 substituents)
c     data radius/1.000d0,4*0.0d0,1.80d0,1.85d0,1.65d0,2*0.0d0,6*0.0d0,
c    $            1.937d0,3*0.0d0,83*0.d0/
c         Experiment
c     data radius/1.200d0,4*0.0d0,1.80d0,2.05d0,1.65d0,2*0.0d0,6*0.0d0,
c    $            1.937d0,3*0.0d0,83*0.d0/
     
c
c     --- define centers and radii for the bubble ---
c         the atomic coordinates enter in c(3,nat),
c         atomic numbers in ian(nat), solvent in rad(nat).
c         if the solvent radii have been previously defined, they will
c         enter with a non-zero radius.
      do 10 i=1,nat
         xcnt(i)=c(1,i)
         ycnt(i)=c(2,i)
         zcnt(i)=c(3,i)
         if(rad(i).eq.zero) then
            rad(i)=radius(ian(i))/toang
         endif
   10    continue
c
c
      return
      end
