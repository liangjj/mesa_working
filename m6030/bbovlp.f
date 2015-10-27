      subroutine bbovlp(s,grid,basis,ncon,npt)
      implicit integer (a-z)
      real *8 s, grid, basis
      dimension s(ncon,ncon), grid(4,npt), basis(npt,ncon)
      do 10 i=1,ncon
         do 20 j=1,i
            do 30 k=1,npt
               s(i,j)=s(i,j)+basis(k,i)*grid(4,k)*basis(k,j)
   30       continue
   20    continue
   10 continue
      return
      end
