*deck kinmat
      subroutine kinmat(ply,dply,ddply,pt,wt,ekin,rbox,l,nprim,nq)
      implicit integer (a-z)
      real *8 ply, dply, ddply, pt, wt, ekin, rbox
      dimension ply(nq,0:nprim), dply(nq,0:nprim), ddply(nq,0:nprim)
      dimension ekin(0:nprim,0:nprim), wt(nq), pt(nq)
      common/io/ inp, iout
c     kinetic energy matrix elements with bloch operator added.
c
      lfac=l*(l+1)
      do 10 i=0,nprim
         do 20 j=0,i
            ekin(i,j)=0.d0
            do 30 k=1,nq
                ekin(i,j)=ekin(i,j)+ply(k,i)*wt(k)*( ddply(k,j) -
     1                         lfac*ply(k,j)/(pt(k)*pt(k)))
   30       continue
            ekin(i,j)=ekin(i,j)-ply(nq,i)*dply(nq,j)
            ekin(i,j)=-.5d0*ekin(i,j)
            ekin(j,i)=ekin(i,j)                              
   20    continue   
   10 continue
      return
      end
