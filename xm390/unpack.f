*deck @(#)unpack.f	2.1  10/10/91
      subroutine unpack (pword,i,j,k,l,m)
c     compile ftn with trace off and opt=1 or opt=2
c     pword is not altered
      kfn(i)=i+4
      i=and(compl(mask(kfn(50))),shift(pword,kfn(10)))
      j=and(compl(mask(kfn(50))),shift(pword,kfn(20)))
      k=and(compl(mask(kfn(50))),shift(pword,kfn(30)))
      l=and(compl(mask(kfn(50))),shift(pword,kfn(40)))
      m=and(compl(mask(kfn(55))),shift(pword,kfn(45)))
      return
      end
