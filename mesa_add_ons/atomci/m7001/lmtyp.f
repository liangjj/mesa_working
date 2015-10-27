*deck @(#)lmtyp.f
c***begin prologue     lmtyp
c***date written       920417   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           lmtyp, link 7001
c***author             schneider, barry (nsf)
c***source             m7001
c***purpose            return l and m value for a given
c***                   function type.
c***
c***description        fills arrays to correspond to function
c***                   description. dimensioned to i orbitals. 
c***references       
c
c***routines called
***end prologue       lmtyp
      subroutine lmtyp(n,l,m,ylm)
      implicit integer (a-z)
      parameter (dimsym=49)
      character *8 types, ylm
      dimension types(dimsym), ltyp(dimsym), mtyp(dimsym)
      data types/'s','p(-1)','p(0)','p(1)','d(-2)','d(-1)','d(0)',
     1           'd(1)','d(2)','f(-3)','f(-2)','f(-1)','f(0)','f(1)',
     2           'f(2)','f(3)','g(-4)','g(-3)','g(-2)','g(-1)','g(0)',
     3           'g(1)','g(2)','g(3)','g(4)','h(-5)','h(-4)','h(-3)',
     4           'h(-2)','h(-1)','h(0)','h(1)','h(2)','h(3)','h(4)',
     5           'h(5)','i(-6)','i(-5)','i(-4)','i(-3)','i(-2)','i(-1)',
     6           'i(0)','i(1)','i(2)','i(3)','i(4)','i(5)','i(6)'/
      data ltyp/0,1,1,1,2,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,5,5,
     1          5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6/
      data mtyp/0,-1,0,1,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,
     1          0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,
     2         -1,0,1,2,3,4,5,6/
      l=ltyp(n)
      m=mtyp(n)
      ylm=types(n)
      return
      end
