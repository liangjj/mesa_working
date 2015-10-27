*deck @(#)noshel.f	5.1  11/6/94
      subroutine noshel(noprim,shelno,nat,nbtype,nshell)
      implicit integer (a-z)
c
      integer noprim(nat,nbtype),shelno(nat,nbtype)
c
c     ----- give each shell present a unique number from 1 to nshell
c
      nshell=0
      do 2 atom=1,nat
         do 1 type=1,nbtype
            if (noprim(atom,type).gt.0) then
               nshell=nshell+1
               shelno(atom,type)=nshell
            else
               shelno(atom,type)=0
            end if
    1    continue
    2 continue
c
c
      return
      end
