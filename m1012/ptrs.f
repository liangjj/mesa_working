*deck @(#)ptrs.f	5.1  11/6/94
      subroutine ptrs(minshl,maxshl,numshl,nshell,pt,nnp,nindep,
     #                orbshl,num)
c
c***begin prologue     ptrs
c***date written       861210  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           independent scf orbital rotations
c***author             saxe, paul (lanl)
c***source             @(#)ptrs.f	5.1   11/6/94
c***purpose            to enumerate the independent rotations twixt orbitals
c***description
c
c***references
c***routines called
c***end prologue       ptrs
c
      implicit integer (a-z)
c
      integer numshl(nshell),minshl(nshell),maxshl(nshell)
      integer pt(nnp),orbshl(num)
c
c     ----- the pointers to independent rotations -----
c
      call izero(pt,nnp)
      nindep=0
      do 5 ishell=2,nshell
         do 4 jshell=1,ishell-1
            do 3 i=minshl(ishell),maxshl(ishell)
               ia=i*(i-1)/2
               do 2 j=minshl(jshell),maxshl(jshell)
                  nindep=nindep+1
                  pt(ia+j)=nindep
    2          continue
    3       continue
    4    continue
    5 continue
c
      call iosys('write integer "number of independent rotations" '//
     $           'to rwf',1,nindep,0,' ')
      call iosys('write integer "independent rotations" to rwf',
     #            nnp,pt,0,' ')
c
c     ----- and the shell of each orbital -----
c
      do 7 shell=1,nshell
         do 6 i=minshl(shell),maxshl(shell)
            orbshl(i)=shell
    6    continue
    7 continue
c
c
      return
      end
