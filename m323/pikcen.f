*deck @(#)pikcen.f	5.1  11/6/94
      function pikcen(symcen,nsymcn,cencod)
c
c*** module to step through a canonical four-centre list, putting
c    the four centres in the array symcen. nsymcen is the number of
c    centres, and cencod reports on equalities between the centres:
c         if 1=2 bit 1 is set
c         if 3=4 bit 2 is set
c         if 1=3 and 2=4 bit 3 is set
c
c    pikcen returns a .false. until the list is exhausted, when it
c    returns a .true.
c
c    initialisation is accomplished by setting cencod to a negative value.
c
c       paul saxe       28 june 1984                     lanl
c
      implicit integer (a-z)
c
      integer symcen(4)
      logical pikcen
c
      pikcen=.false.
c
c     ----- initialise if so requested -----
c
      if (cencod.lt.0) then
         symcen(1)=1
         symcen(2)=1
         symcen(3)=1
         symcen(4)=1
         cencod=7
         return
      end if
c
c     ----- increment last centre and check -----
c
      symcen(4)=symcen(4)+1
      if (symcen(4).le.symcen(3)) go to 10
c
      symcen(4)=1
      symcen(3)=symcen(3)+1
      if (symcen(3).le.symcen(1)) go to 10
c
      symcen(3)=1
      symcen(2)=symcen(2)+1
      if (symcen(2).le.symcen(1)) go to 10
c
      symcen(2)=1
      symcen(1)=symcen(1)+1
      if (symcen(1).le.nsymcn) go to 10
c
c     ----- if reached here, have finished the canonical list -----
c
      pikcen=.true.
      return
c
c
   10 continue
      cencod=0
      if (symcen(1).eq.symcen(2)) cencod=1
      if (symcen(3).eq.symcen(4)) cencod=cencod+2
      if (symcen(1).eq.symcen(3).and.symcen(2).eq.symcen(4))
     #                                     cencod=cencod+4
      return
      end
