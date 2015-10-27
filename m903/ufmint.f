*deck @(#)ufmint.f	5.1  11/6/94
      subroutine ufmint(orbsym,norbs,nsym,numsym,offsym,ijpt,nnp0,
     #                  nnp,ngsym,hsym,h,gsym,g)
c
c***begin prologue     ufmint
c***date written       871022  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords
c***author             saxe, paul (lanl)
c***source
c***purpose            to un-symmetry block the integrals.
c***description
c
c***references
c***routines called
c***end prologue       ufmint
c
      implicit integer (a-z)
c
      real*8 small
      real*8 h(nnp),hsym(nnp0),g(nnp,nnp),gsym(ngsym)
      integer orbsym(norbs),ijpt(nnp),numsym(0:nsym-1),offsym(0:nsym-1)
c
      parameter (small=1.0d-05)
c
c     ----- compact the one-electron integrals -----
c
      do 25 i=1,norbs
         ia=i*(i-1)/2
         is=orbsym(i)
         do 20 j=1,i
            ij=ia+j
            ijs=xor(is,orbsym(j))
            if (ijs.eq.0) then
               h(ij)=hsym(ijpt(ij))
            end if
   20    continue
   25 continue
c
c     ----- and the two-electron integrals -----
c
      do 45 i=1,norbs
         ia=i*(i-1)/2
         is=orbsym(i)
         do 40 j=1,i
            ij=ia+j
            ijs=xor(is,orbsym(j))
            ija=offsym(ijs)+(ijpt(ij)-1)*numsym(ijs)
            do 35 k=1,norbs
               ka=k*(k-1)/2
               ks=orbsym(k)
               do 30 l=1,k
                  kl=ka+l
                  kls=xor(ks,orbsym(l))
                  if (kls.eq.ijs) then
                     g(ij,kl)=gsym(ija+ijpt(kl))
                  end if
   30          continue
   35       continue
   40    continue
   45 continue
c
      return
      end
