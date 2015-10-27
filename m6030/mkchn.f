      subroutine mkchn(pii,t,ncon,ngi,lgi)
      implicit integer (a-z)
      real *8 pii, t
      dimension pii(ncon,ncon), t(ncon,ncon)
      dimension lgi(ngi)
      do 10 i=1,ngi
	 ipt=lgi(i)
         do 20 j=1,ngi
	    jpt=lgi(j)
            pii(ipt,jpt)=pii(ipt,jpt) + t(ipt,jpt)
   20    continue
   10 continue
      return
      end
