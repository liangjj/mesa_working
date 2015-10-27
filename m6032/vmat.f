*deck vmat
      subroutine vmat(ply,pt,wt,v,v0,vl,nprim,nq)
      implicit integer (a-z)
      real*8 ply, pt, wt, v, v0, vl
      dimension ply(nq,0:nprim), v(0:nprim,0:nprim), pt(nq), wt(nq)
      common/io/ inp, iout
      do 5 i=1,nq
         if (pt(i).ge.vl) then
             icut=i
             go to 6
         endif
    5 continue
      call lnkerr('screwed up cutoff')
    6 write(iout,*) ' cutoff point = ', icut   
c     potential energy matrix elements
c
      do 10 i=0,nprim
         do 20 j=0,i
            v(i,j)=0.d0
            do 30 k=1,icut
               v(i,j) = v(i,j) + v0*ply(k,i)*ply(k,j)*wt(k)
   30       continue     
            v(j,i)=v(i,j)
   20    continue   
   10 continue
      return
      end


