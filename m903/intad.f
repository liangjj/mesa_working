*deck @(#)intad.f	5.1  11/6/94
      subroutine intad(orbsym,norbs,nsym,numsym,offsym,ijpt,ngsym,nnp)
c
c***begin prologue     intad
c***date written       860920  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           integral indexing
c***author             saxe, paul (lanl)
c***source
c***purpose            to calculate indices for a symmetry packed integral list
c***description
c
c***references
c***routines called    (none)
c***end prologue       intad
c
      implicit integer (a-z)
c
      integer orbsym(norbs),ijpt(nnp),numsym(0:nsym-1),offsym(0:nsym-1)
c
c     ----- form the pointer array for ij pairs -----
c
      call izero(numsym,nsym)
c
      do 1 i=1,norbs
         orbsym(i)=orbsym(i)-1
    1 continue
c
      do 10 i=1,norbs
         ia=i*(i-1)/2
         is=orbsym(i)
         do 5 j=1,i
            ijs=xor(is,orbsym(j))
            ij=ia+j
            numsym(ijs)=numsym(ijs)+1
            ijpt(ij)=numsym(ijs)
    5    continue
   10 continue
c
c     ----- and the global symmetry offsets -----
c
      offsym(0)=0
      do 15 sym=1,nsym-1
         offsym(sym)=offsym(sym-1)+numsym(sym-1)**2
   15 continue
c
      ngsym=offsym(nsym-1)+numsym(nsym-1)**2
c
      return
      end
