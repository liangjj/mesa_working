*deck @(#)balpar.f	5.1  11/6/94
      function balpar(string,par)
c
c***begin prologue     balpar
c***date written       860814  (yymmdd)
c***revision date      870131  (yymmdd)
c
c   31 january 1987  pws at lanl
c      removing the declaration of an unused variable 'index'.
c
c***keywords           balanced parenthesis
c***author             saxe, paul (lanl)
c***source             @(#)balpar.f	5.1   11/6/94
c***purpose            to find a balancing parenthesis in a string.
c***description        integer=balpar(string,par)
c        when passed an input string starting with a parenthesis, this routine
c        will return the location of the balanced parenthesis.
c
c  on input:
c            string    character*(*)
c                      the input string.
c
c  on return:
c            par       integer
c                      the location of the balancing parenthesis in the string.
c
c            balpar=par
c
c      special conditions: if the string does not start with a left parenthesis,
c           par will be returned as 0. if there is no balanced right
c           parenthesis, par is set to -1.
c
c***references
c***routines called
c***end prologue      balpar
c
      implicit integer (a-z)
c
      integer balpar
      character*(*) string
      character*1 left,right
      integer par
c
      parameter (left='(',right=')')
c
c     ----- watch out for zero length string -----
c
      if (len(string).lt.1) then
         balpar=0
         par=0
         return
      end if
c
c     ----- or string not starting with a left parenthesis -----
c
      if (string(1:1).ne.left) then
         balpar=0
         par=0
         return
      end if
c
c     ----- start with a count of 1 for left parenthesis, increment
c           for each left, decrement for right, till get to zero
c
      depth=1
      do 10 i=2,len(string)
         if (string(i:i).eq.left) then
            depth=depth+1
         else if (string(i:i).eq.right) then
            depth=depth-1
            if (depth.eq.0) then
               balpar=i
               par=i
               return
            end if
         end if
   10 continue
c
c     ----- we have reached the end of the string without match -----
c
      balpar=-1
      par=-1
c
c
      return
      end
