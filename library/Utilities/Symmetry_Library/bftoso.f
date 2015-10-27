*deck @(#)bftoso.f	5.1  11/6/94
      subroutine bftoso(nfunc,nsymat,nlamda,nrep,c,coeffs,
     #                  natoms,relatm,nbf,bfstrt,
     #                  nocont,dump)
c
      implicit integer (a-z)
c
      real*8 c(nfunc,nsymat,nlamda,nrep)
      real*8 coeffs(nbf,*)
      integer relatm(nsymat),bfstrt(natoms)
      logical dump
c
      common/io/inp,iout
c
   90 format(1x,'bftoso')
c
c     expand this block of salc coefficients into the full salc
c     transformation matrix.
c
c
      so=0
      do 300 rep=1,nrep
         do 200 l=1,nlamda
            do 100 j=1,nocont
               so=so+1
               do 50 atom=1,nsymat
                  iatom=relatm(atom)
                  bf=bfstrt(iatom)+(j-1)*nfunc
                  do 40 func=1,nfunc
                     bf=bf+1
                     coeffs(bf,so)=c(func,atom,l,rep)
   40             continue
   50          continue
  100       continue
  200    continue
  300 continue
c
c
c     if(dump) then
c        write(iout,90) 
c        call matout(coeffs,nbf,so,nbf,so,iout)
c     endif
c
c
      return
      end
