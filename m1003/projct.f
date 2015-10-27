*deck @(#)projct.f	5.1  11/6/94
      subroutine projct(v,x,nm,nvm,nspinm)
      implicit real*8(a-h,o-z)
c
      dimension  v(2),x(2),nm(2),nvm(2)
      real*8 sdot
c
      kx=1
      do 30 k=1,nspinm
         nv=nvm(k)
         n=nm(k)
         jx=kx
         do 20 j=1,nv
            ix=kx
            do 10 i=1,nv
c
               xx=sdot(n,v(jx),1,x(ix),1) 
               rr=-xx
               call saxpy(n,rr,x(ix),1,v(jx),1)
c
               ix=ix+n
  10        continue
            jx=jx+n
 20      continue
         kx=kx+nv*n
30    continue
c
      return
      end
