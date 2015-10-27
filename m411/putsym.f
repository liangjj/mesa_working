*deck @(#)putsym.f	5.1  11/6/94
      subroutine putsym(sym,n,val)
      implicit integer (a-z)
      dimension sym(n)
      do 10 i=1,n
         sym(i)=val
   10 continue
      return
      end
