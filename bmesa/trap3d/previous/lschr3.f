*deck lschr3.f
c***begin prologue     lschr3
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           3-d schroedinger equation
c***author             schneider, barry (nsf)
c***source             trap3d
c***purpose            calculate eigenvalues of three dimensional
c***                   schreodinger equation.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       lschr3
      subroutine lschr3(h01,h02,h03,v,ham,eig,ind,n,nd,numdim,prnt,ops)
      implicit integer (a-z)
      real*8 h01, h02, h03, v, ham, eig
      character*80 title
      dimension nd(3)
      dimension h01(nd(1),nd(1)), h02(nd(2),nd(2)), h03(nd(3),nd(3))
      dimension v(n), ham(n,n), eig(n), ind(n,*)
      common/io/inp, iout
c
      write(iout,1)
      call rzero(ham,n*n)      
      do 10 i=1,n
         ham(i,i) = v(i)
 10   continue
      cnti=0
      do 20 i=1,n
         cnti=cnti+1
         i1=ind(cnti,1)
         j1=ind(cnti,2)
         k1=ind(cnti,3)
         cntj=0
         do 30 j=1,i
            cntj=cntj+1
            i2=ind(cntj,1)
            j2=ind(cntj,2)
            k2=ind(cntj,3)
            di1i2=0
            dj1j2=0
            dk1k2=0
            if(i1.eq.i2) then
               di1i2=1
            endif
            if(j1.eq.j2) then
               dj1j2=1
            endif
            if(k1.eq.k2) then
               dk1k2=1
            endif
            ham(cnti,cntj) = ham(cnti,cntj) + h01(i1,i2)*dj1j2*dk1k2
     1                                      + h02(j1,j2)*di1i2*dk1k2
     2                                      + h03(k1,k2)*di1i2*dj1j2
 30      continue
 20   continue   
      do 40 i=1,n
         do 50 j=1,i
            ham(j,i) = ham(i,j)
 50      continue
 40   continue   
      call rdiag(ham,eig,v,n)
    1 format(/,5x,'diagonalize time-independent three dimensional '
     1            'schroedinger equation')
      end       
