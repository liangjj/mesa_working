*deck kinmat
      subroutine kinmat(ekp,ekin,vec,s,scr,l,rbox,power,nprim,ncon)
      implicit integer (a-z)
      real *8 ekp, ekin, vec, s, scr, rbox
      dimension ekp(nprim,nprim), ekin(ncon,ncon), power(nprim)
      dimension vec(nprim,*), scr(*), s(nprim,nprim)
      common/io/ inp, iout
c     kinetic energy matrix elements with bloch operator added.
c
      lfac=l*(l+1)
      do 10 i=1,nprim
         do 20 j=1,i
            ekp(i,j)=0.d0
            pre=power(j)*(power(j)-1)
            lsum=power(i)+power(j)-1
            if (pre.ne.0) then
                ekp(i,j)=ekp(i,j)-.5d0*pre*(rbox**lsum/lsum)
            endif
            if (power(j).ne.0) then
                ekp(i,j)=ekp(i,j)+.5d0*power(j)*(rbox**lsum)
            endif
            if (l.ne.0) then
                ekp(i,j)=ekp(i,j)+.5d0*lfac*(rbox**lsum/lsum)
            endif
            ekp(i,j)=ekp(i,j)*s(i,i)*s(j,j)
            ekp(j,i)=ekp(i,j)
 20      continue   
 10   continue
c     transform
c   
      call ebtc(scr,vec,ekp,ncon,nprim,nprim)
      call ebc(ekin,scr,vec,ncon,nprim,ncon)
      return
      end
