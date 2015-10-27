*deck @(#)mcfock.f	5.1  11/6/94
      subroutine mcfock(dab,f,nco,nao,mrs,hess,lhaa,lhcc,deg)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcfock.f	5.1   11/6/94
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      implicit real*8(a-h,o-z)
cc
cmp   extended dummy dab,f,hess
cc
      common / number / zero,pt5,one,two,four,eight
      dimension dab(2),f(2),hess(2)
c
c---------------------------------------------------
c    fock operator contributions to the diagonal
c    core-core blocks of the hessian is now in mccorjk
c---------------------------------------------------
c
c---------------------------------------------------c
c    add fock operator contributions to the diagonal
c    active-active blocks of the hessian
c---------------------------------------------------c
c
      ia=0
      lh=lhaa
      if(nao.eq.0) go to 31
c
      do 30 i=1,nao
         do 20 j=1,i
            factab=dab(ia+j)
            if(i.eq.j) factab=factab*pt5
cc
            call mcstv(hess(lh),f,factab,mrs)
cc
            lh=lh+mrs
 20      continue
         ia=ia+nao
 30   continue
 31   continue
c
      return
      end
