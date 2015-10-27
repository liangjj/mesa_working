*deck trnham.f
      subroutine trnham(ham,eig,u,scr,n,m)
      implicit integer (a-z)
      real*8 ham, eig, u, scr
      character*80 title
      dimension ham(n,*), eig(m), u(n,*), scr(n,*)
      common/io/inp, iout 
      title='untransformed hamiltonian'
c      call prntrm(title,ham,n,n,n,n,iout)
      size=n-m
      call ebcxx(scr,ham(m+1,1),u,size,m,m,n,n,n)
      do 10 i=1,m
         do 20 j=1,m
            ham(i,j)=0.d0
 20      continue
         ham(i,i)=eig(i)
 10   continue
      count=0
      do 30 i=m+1,n
         count=count+1
         do 40 j=1,m
            ham(i,j)=scr(count,j)
            ham(j,i)=ham(i,j)
 40      continue
 30   continue   
      title='transformed hamiltonian'
c      call prntrm(title,ham,n,n,n,n,iout)
      return
      end       

