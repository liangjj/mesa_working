*deck seppre.f
c***begin prologue     seppre
c***date written       970531   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           davidson, preconditioner
c***author             schneider, barry (nsf)
c***source             
c***purpose            prepare what is needed for a preconditioner
c***                   based on a hamiltonian which is a sum of
c***                   separate, commuting operators.
c***
c***description        
c                      
c***references         
c
c***routines called    
c***end prologue       seppre
      subroutine seppre(p,h1,h2,h3,n1,n2,n3,n,dim)
      implicit integer (a-z)
      real*8 h1, h2, h3, h, work
      character*80 title
      dimension h1(n1,n1), h2(n2,n2), h3(n3,n3)
      common/io/inp, iout
      pointer (p,h(1))
      pointer (pw,work(1))
      write(iout,1)      
c
c     get the separate eigenvalues and transformation matrices
c 
      maxn=n1
      u1=1
      eig1=u1+n1*n1
      need=eig1+n1
      if(dim.gt.1) then
         maxn=max(maxn,n2)
         u2=need
         eig2=need+n2*n2
         need=eig2+n2
      endif
      if(dim.gt.2) then
         maxn=max(maxn,n3)
         u3=need
         eig3=u3+n3*n3
         need=eig3+n3
      endif
      need=wpadti(need)
      call getmem(need,p,ngot,'seppre',0)
      scrwd=wptoin(5*maxn)
      call getmem(scrwd,pw,nscr,'scr',0)
      call copy(h1,h(u1),n1*n1)
      call diddle(h(u1),n1)
      call dsyev('v','l',n1,h(u1),n1,h(eig1),work,5*n1,info)
      title='eigenvalues for dimension 1'
      call prntrm(title,h(eig1),n1,1,n1,1,iout)
      if(dim.gt.1) then
         call copy(h2,h(u2),n2*n2)
         call dsyev('v','l',n2,h(u2),n2,h(eig2),work,5*n2,info)
         title='eigenvalues for dimension 2'
         call prntrm(title,h(eig2),n2,1,n2,1,iout)
      endif
      if(dim.gt.2) then
         call copy(h3,h(u3),n3*n3)
         call dsyev('v','l',n3,h(u3),n3,h(eig3),work,5*n3,info)
         title='eigenvalues for dimension 3'
         call prntrm(title,h(eig3),n3,1,n3,1,iout)
      endif
      call getmem(-nscr,pw,dum,'scr',dum)
      return
 1    format(/,5x,'setting up for preconditioning based on '
     1            'separable hamiltonian')
      end       





