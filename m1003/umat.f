*deck @(#)umat.f	5.1  11/6/94
      subroutine umat(um,del,mixh,nco,nao,nob)
      implicit integer(a-z)
c
      real*8 um(nob,2),del(2)
      integer mixh(nob,2)
c
      nco1=nco+1
      if(nco.eq.0)go to 25
      do 10 i=1,nco
         do 20 j=nco1,nob
            idel=mixh(j,i)
            if(idel.eq.0)go to 20
            um(j,i)=-del(idel)
            um(i,j)=del(idel)
  20     continue
  10  continue
c
  25  continue
c
      if(nao.eq.0) return
c
      is=nco1
      ie=nco+nao
      do 30 i=is,ie
         do 40 j=i,nob
            idel=mixh(j,i)
            if(idel.eq.0)go to 40
            um(j,i)=-del(idel)
            um(i,j)=del(idel)
   40    continue
   30 continue
c
      return
      end
