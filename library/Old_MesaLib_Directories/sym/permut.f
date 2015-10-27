*deck @(#)permut.f	5.1  11/6/94
      subroutine permut(c,natoms,t,nop,temp,atprmt,atmchg)
c***begin prologue     permut.f
c***date written       850101  (yymmdd) 
c***revision date      900507  (yymmdd)    
c
c     7 may     bhl at llnl
c       allowing scattering and ghost centers to coincide with atomic
c       centers
c***keywords           symmetry, permutation 
c***author             saxe, paul 
c***source             @(#)permut.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c       start here
c
c***end prologue      permut.f
      implicit integer (a-z)
c
      real*8 fuzz,zap
      real*8 c(3,natoms),t(3,3,nop),temp(3,natoms)
      integer atprmt(natoms,nop)
      real*8 atmchg(natoms)
      common/io/inp,iout
c
      parameter (fuzz=1.0d-06)
 1000 format(1x,i3,3f10.5)
c
c     ----- apply each operation to the coordinates, then see which
c     atom maps to which
c
      do 100 op=1,nop
         call ebc(temp,t(1,1,op),c,3,3,natoms)
         do 20 iatom=1,natoms
            do 10 jatom=1,natoms
               zap=abs(c(1,iatom)-c(1,jatom))
     $            +abs(c(2,iatom)-c(2,jatom))
     $            +abs(c(3,iatom)-c(3,jatom))
               if ((iatom.ne.jatom).and.(zap.le.fuzz)) then
c                 these two distinct centers are coincident.
c                 must be a scattering site, disregard test.
               else
                  if (abs(temp(1,jatom)-c(1,iatom)).lt.fuzz.and.
     #                abs(temp(2,jatom)-c(2,iatom)).lt.fuzz.and.
     #                abs(temp(3,jatom)-c(3,iatom)).lt.fuzz.and.
     $                atmchg(iatom).eq.atmchg(jatom)) then
                     atprmt(iatom,op)=jatom
                     go to 15
                  end if
               end if
 10         continue
c
            write(iout,*) ' problems in permut'
            write(iout,*) 'original coordinates'
            do 12 i=1,natoms
               write(iout,1000) i,(c(j,i),j=1,3)
   12       continue
            write(iout,*) 'transformed coordinates,operation:',op
            do 13 i=1,natoms
               write(iout,1000) i,(temp(j,i),j=1,3)
   13       continue
            call lnkerr('not symmetric atom found in permut')
c
 15         continue
 20      continue
 100  continue
c
c
      return
      end
