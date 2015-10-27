*deck @(#)comand.f	5.1  11/6/94
      subroutine comand(num,key)
c***begin prologue     comand
c***date written       850601  (yymmdd)
c***revision date      900130  (yymmdd)
c
c     30 january 1990   rlm at lanl
c         rewriting for bsd 4.3 unix / stardent titan.
c     7 february 1987   pws at lanl
c         rewriting functionality for a bsd 4.2 unix version on the
c         sun 3/50 and 3/160 workstations.
c
c***keywords           keyword, reassignment, execute line
c***author             saxe, paul  (lanl)
c***source             @(#)comand.f	5.1   11/6/94
c***purpose            parses the execute line.
c***description
c                      call comand(num,key)
c                        num     the number of keywords for which to search.
c                        key     character array containing the keywords.
c
c***references
c***routines called    
c***end prologue       comand
c
      implicit integer(a-z)
c
      character*(*) key(num)
      character*128 arg,temp,chrkey,value,getarg
c
c     ----- find the number of arguments on the execute line -----
c
      nargs=iargc()
c
c     ----- if there are none, then leave the keywords alone -----
c
      if (nargs.eq.0) return
c
c     ----- otherwise, pick up each keyword in turn and see if it 
c           is to be replaced
c
      do 100 i=1,nargs
         value=getarg(i,arg)
c
c        ----- loop through the list of keywords to see if this is it
c
         do 90 j=1,num
            temp=chrkey(arg,key(j),' ',' ')
            if (temp.ne.' ') then
               key(j)=temp
               go to 100
            end if
   90    continue
  100 continue
c
c
      return
      end
