*deck orthbs
      subroutine orthbs(smat,eig,vec,scr,rbox,power,nprim,ncon)
      implicit integer (a-z)
      real *8 smat, vec, scr, eig, rbox, tol, sqeig
      dimension smat(nprim,nprim), power(nprim), vec(nprim,*)
      dimension scr(*), eig(nprim)
      parameter ( tol=1.d-10)
      common/io/ inp, iout
c     overlap matrix elements
c
      do 10 i=1,nprim
         do 20 j=1,i
            lsum=power(i)+power(j)+1
            smat(i,j)=rbox**lsum/lsum
 20      continue
         eig(i)=1.d0/sqrt(smat(i,i))   
 10   continue
      do 15 i=1,nprim
         do 16 j=1,i
            vec(i,j)=smat(i,j)*eig(i)*eig(j)
            vec(j,i)=vec(i,j)
   16    continue
         smat(i,i)=eig(i)
   15 continue             
c     diagonalize
c   
      call tred2(nprim,nprim,vec,eig,scr,vec)
      call tql2(nprim,nprim,eig,scr,vec,ierr)
      write(iout,1) (eig(i),i=1,nprim)
      ncon=0
      do 30 i=1,nprim
         if (eig(i).gt.0.d0) then
             if (abs(eig(i)).gt.tol) then
                 sqeig=1.d0/sqrt(eig(i))
                 ncon=ncon+1
                 eig(ncon)=eig(i)
                 do 40 j=1,nprim
                    vec(j,ncon)=sqeig*vec(j,i)
 40              continue
             endif    
         endif
 30   continue
      write(iout,2) ncon
 1    format(/,1x,'eigenvalues of overlap matrix',(/,5e15.8))
 2    format(/,1x,'number contracted functions = ',i4)
      return
      end
