*deck %W%  %G%
      subroutine comand(num,key)
c***begin prologue     comand
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           keyword, reassignment, execute line
c***author             norton, john (lanl)
c***source             %W%   %G%
c***purpose            parses the execute line.
c***description
c                      call comand(num,key)
c                        num     the number of keywords for which to search.
c                        key     hollerith array containing the keywords.
c
c***references
c***routines called    qgetsym(ctss), syserr(ctss), msgtty(ctss), ljust(ctss)
c***end prologue       comand
      implicit integer(a-z)
c
c  key is really character*128, let's try and fool it with an array 16xn
c
      dimension key(16,*),idlist(4)
      common /qrunbfc/ mess(64), irep(101)
c
c       perform keyword reassignment from execution line.
c
c       input array -key(num)- contains the keywords which are to be
c       matched. the keyword list on the execution line should be of
c       form -
c               key1=val1 key2=val2 .....
c       the -key- element containing -key1- is then replaced by the
c       hollerith -val1-, etc. if -key1- is not found, the -key-
c       element containing it is left unchanged
c
c       an example:
c       l1c inp=water,out=owater,rwf=rwf1
c       would execute link1, reading input from the local file water,
c       writing to owater, and naming the read/write file rwf1.
c       the other files default to int,chk,d2e.
c       the line continuation character is an ampersand,"&".
c
c       this routine is essentially lifted from the innards of
c       the cftlib routine filerep.
c     note that /qrunbfc/ is a cftlib common block.
c       cray/ctss version: calls cftlib routines qgetsym,syserr,msgtty,ljust.
c
      data idlist/1h=,1h,,1h ,-1/
c       get the execute line
      nn=1
      isub=1
      nw=63
      nsym=0
      call qgetsym (nn,irep(isub),nw,n,2,idlist)
c       if nothing is on the execute line,just return
      if (n.lt.0) return
c
c        see if the line is to be continued
10    isubp=isub+n-1
      nsym=nsym+n
      if (irep(isubp).ne.1r&) go to 20
c       yes. see if there is any room.
      isub=isubp
      nsym=nsym-1
      nw=min(63,101-nsym)
      if (nw.lt.3) call syserr (53,16htoo many symbols,16)
c       yes. send the prompt and read the next line.
      itemp=or(shift(2r? ,48),shift(04b,40))
      ntemp=3
      call msgtty(itemp,ntemp)
      nn=0
      call qgetsym (nn,irep(isub),nw,n,2,idlist)
      go to 10
c
c       all the symbols have been found. do the keyword reassignments.
20    n=1
30    call ljust (irep(n),ifa)
      if (irep(n+1).ne.1r=) go to 60
      n=n+2
      if (n.gt.nsym) go to 60
      call ljust (irep(n),ifb)
      do 40 i=1,num
40    if (key(1,i).eq.ifa) key(1,i)=ifb
50    n=n+1
      if (n.gt.nsym) return
      if (irep(n).eq.1r,) go to 50
      go to 30
c
c       error in syntax
60    call syserr (54,35herror in keyword replacement syntax,35)
      return
      end
