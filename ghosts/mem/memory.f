*deck %W% %G%
      function memory(str,i1,base)
c
c***begin prologue     memory
c***date written       910815   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           memory, pointers 
c
c***author             martin, richard (lanl)
c***source             %W% %G% 
c***purpose            memory manager 
c
c***description        memory is an attempt at a flexible, user-friendly
c memory management package which hides the machine dependent aspects.
c 
c the form of a call to memory is
c
c          pointer=memory(directive,nwords,base)
c
c where 'directive' is a character string giving (in english?) the
c action to be taken. the second argument gives the number of words
c requested, and the third argument is the base address to which the pointer
c refers.
c
c the function returns either the offset relative to base which
c can be used to access core, or a zero in cases where core
c is being released, etc.
c
c      in the following description of the directives <handle> is an
c internal identifier used to refer to the allocated space.
c parenthesis () indicate a choice of words, one of which must be used;
c and brackets [] indicate an optional keyword or phrase. also, all
c (hopefully) non-essential words such as 'to', 'on' and 'of' may be
c omitted wherever they appear; however, their inclusion is suggested
c as an aid to readability.
c      the possible directives are:
c        'initialize'
c        'get (integer) (real) pointer <handle>'
c        'release (integer) (real) pointer <handle>'
c        'define block'
c        'release block'
c
c  'initialize'
c      this must called at the beginning of the job and it initializes
c      the common blocks, etc.  in addition, it allocates the original
c      space which the subsequent calls will divvy up.  the original
c      space is allocated relative to the beginning of blank common.
c
c      nwords ... the size of the main memory area in which to work. 
c                 this is in "integer" words.
c
c      base   ... the address to which all else is referred.
c
c      example:
c         integer a
c         common//a(1)
c
c         ioff=memory('initialize',1000000,a(1))
c         call pm1(a(ioff),a(ioff))
c
c         stop
c         end
c
c         subroutine pm1(a,z)
c         integer a(*)
c         real*8 z(*)
c
c      this will allocate 1000000 "integer" words of core. it begins at
c      a(ioff+1).  in the subroutine pointers to either real or integer
c      variables can be generated by subsequent calls to memory.
c
c  'get (integer) (real) pointer <handle>'
c      this returns a pointer to a region of core. 
c
c      nwords ... the number of words needed.
c                 the length of core allocated depends on integer
c                 vs. real.  
c
c      base   ... the address to which the pointer refers.
c  'release (integer) (real) pointer <handle>'
c      this releases the space allocated to <handle>
c
c      nwords ... unused.
c      base   ... unused.
c    
c      example:
c         subroutine pm1(a,z)
c         implicit integer(a-z)
c         integer a(*)
c         real*8 z(*)
c
c         nint=10
c         idum=memory('get integer pointer idum',nint,a(1))
c         nreal=100
c         rdum=memory('get real pointer rdum',nreal,z(1))
c
c         call sub1(a(idum),nint,z(rdum),nreal)
c         write(6,*) (a(idum+i),i=0,nint-1)
c         write(6,*) (z(rdum+i),i=0,nreal-1)
c
c         idum=memory('release integer pointer idum',0,0)
c         rdum=memory('release real pointer rdum',0,0)
c
c         return
c         end
c
c
c         subroutine sub1(ints,nints,reals,nreals)
c         implicit integer(a-z)
c         integer ints(nints)
c         real*8 reals(nreals)
c
c         do 10 i=1,nints
c            ints(i)=i
c   10    continue
c         do 20 i=1,nreals
c            reals(i)=float(i)
c   20    continue
c
c         return
c         end 
c      
c      note that we explicitly released the core that was allocated.
c      sometimes it is more convenient to allocate a series of pointers
c      which can then all be released with one call to memory.
c      this is the function of the 'block' directive:
c  'define block'
c      this associates all the pointers subsequently allocated  with
c      a lock of pointers which can all be released with a single call.
c
c      nwords ... unused.
c      base   ... unused.
c  'release block'
c      this releases all the pointers  associated with the block.
c
c      nwords ... unused.
c      base   ... unused.
c         
c***references
c
c***routines called    gettok (io)
c                      inimem (mem)
c                      getpt  (mem)
c                      relpt  (mem)
c                      lnkerr (mdutil)
c
c   common blocks:     mem1, mem2, and io.
c
c***end prologue       memory
c
      implicit integer (a-z)
      integer memory
      character*(*) str
      integer i1
      integer base(1)
c
      character*16 token
c
      character*240 string
      character*80 tmplin
      character*8  handle
      integer narray,addr,length
      logical called
      logical block
c
c     ----- functions -----
c
      character*240 dcaptl
c
      parameter (maxarr=1000)
c
      common /mem1/ narray,addr(0:maxarr),length(0:maxarr)
      common /mem2/ handle(0:maxarr)
c
c
      data called/.false./
      save called
      save block
c
c     set up some common block information.
      if(.not.called) then
         called=.true.
      end if
c
c     ----- transfer the command string to tmplin in case we need
c           to print the string for an error
c
      tmplin=str
c
c     ----- convert the command string to lower case -----
c
      string=dcaptl(str)
c
c     ----- get the first token from 'string' -- the operation -----
c
      pos=0
      call gettok(token,string,pos)
c
c     ----- decide what to do based on 'oper' -----
c
      if (token.eq.'initialize') then
         memory=inimem(base,nwords)
      else if (token.eq.'get') then
         memory=getpt(string,pos,i1,array,array)
      else if (token.eq.'release') then
         memory=relpt(string,pos,array,i1)
      else if (token.eq.'length') then
         call ioopen(string,pos,i1)
      else if (token.eq.'dump') then
         call memdat()
      else if (token.eq.'delete') then
      else
         call lnkerr('memory: unrecognised operation: '//
     #                tmplin)
      end if
c
      return
      end