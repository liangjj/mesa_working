      subroutine bbpot(pij,grid,basis,pot,ncon,npt,ntri,ngi,
     1                                    lgi,ngj,lgj,itri)
      implicit integer (a-z)
      real *8 grid, basis
      real *8 pij, pot
      dimension pij(ncon,ncon), grid(4,npt), basis(npt,ncon)
      dimension pot(npt,ntri)
      dimension lgi(ngi), lgj(ngj)
      do 10 i=1,ngi
	 ipt=lgi(i)
         do 20 j=1,ngj
	    jpt=lgj(j)
            do 30 k=1,npt
               pij(ipt,jpt)=pij(ipt,jpt)+basis(k,ipt)*pot(k,itri)*
     1                                  grid(4,k)*basis(k,jpt)
   30       continue
   20    continue
   10 continue
      return
      end
