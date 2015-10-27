*deck bsabox
      subroutine bsabox(bfbox,dbfbox,vec,s,scr,rbox,power,nprim,ncon)
      implicit integer (a-z)
      real *8 bfbox, dbfbox, vec, s, scr, rbox
      dimension bfbox(ncon), power(nprim), vec(nprim,*), scr(nprim)
      dimension dbfbox(ncon), s(nprim,nprim)
      common/io/ inp, iout
      do 10 i=1,nprim
         scr(i)=s(i,i)*rbox**power(i)
 10   continue
      call ebtc(bfbox,vec,scr,ncon,nprim,1)
      do 20 i=1,nprim
         scr(i)=s(i,i)*power(i)*rbox**(power(i)-1)
 20   continue
      call ebtc(dbfbox,vec,scr,ncon,nprim,1)
      return
      end
