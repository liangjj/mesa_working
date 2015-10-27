*deck oldin.f
      subroutine oldin(vec,eig,n,m,nroots,prnt)
      implicit integer (a-z)
      real*8 vec, eig
      character*80 title
      logical prnt
      dimension vec(n,*), eig(*)
      common/io/inp, iout
       call iosys('read real "davidson iterates" from lamdat',
     1             n*m,vec,0,' ')            
      if(prnt) then
         title='guess eigenvectors'
         call prntrm(title,vec,n,m,n,n,iout)
      endif 
      return
      end       

