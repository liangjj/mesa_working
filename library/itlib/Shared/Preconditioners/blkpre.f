*deck blkpre.f
c***begin prologue     blkpre
c***date written       970531   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           davidson, preconditioner
c***author             schneider, barry (nsf)
c***source             
c***purpose            prepare what is needed for a preconditioner
c***                   based on a blocked hamiltonian.
c***
c***description        
c                      
c***references         
c
c***routines called    
c***end prologue       blkpre
      subroutine blkpre(p,h1,h2,h3,v,n1,n2,n3,n,nblck,dim)
      implicit integer (a-z)
      real*8 h1, h2, h3, h, v, scr
      character*80 title
      dimension h1(n1,n1), h2(n2,n2), h3(n3,n3), v(n)
      common/io/inp, iout
      pointer (p,h(1))
      pointer (pscr,scr(1)), (pscr,iscr(1))
      write(iout,1)
c
c     get the memory
c 
      u=1
      eig=u+nblck*nblck
c
c     hold all eigenvalues
c
      need=wpadti(eig+n)
      call getmem(need,p,ngot,'blkpre',0)
      work=1
      ind=wpadti(work+5*nblck)
      scrwd=ind+n*dim
      call getmem(scrwd,pscr,nscr,'scr',0)
      call rbldag(h(u),h(eig),h1,h2,h3,v,scr(work),iscr(ind),n1,n2,n3,
     1            n,nblck,dim,.false.)
      call getmem(-nscr,pscr,dum,'scr',dum)
      return
 1    format(/,5x,'setting up for preconditioning based on '
     1            'block diagonal hamiltonian')
      end       





