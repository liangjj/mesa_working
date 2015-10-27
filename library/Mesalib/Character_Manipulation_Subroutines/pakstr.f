*deck @(#)pakstr.f	5.1  11/6/94
      subroutine pakstr(string,length)
c***begin prologue     pakstr
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           character, string, pack
c***author             martin, richard (lanl)
c***source             @(#)pakstr.f	5.1   11/6/94
c***purpose            pack down a string. pakstr removes all leading and
c                      trailing blanks, and all embedded blanks but one.
c***description
c                      call pakstr(string,length)
c                        string    the string to pack.
c                        length    the length of the string after packing.
c
c***references
c***routines called    cskipf(chr), cskipb(chr)
c***end prologue       pakstr
      implicit integer(a-z)
      character*(*) string
      character*1 char,lastch
c
c
c
      length=0
      bcur=cskipf(string,' ')
      ecur=cskipb(string,' ')
      if(ecur.eq.0) return
 
c     now pack.
      k=0
      lastch=' '
      do 10 i=bcur,ecur
         char=string(i:i)
         if(char.eq.' '.and.lastch.eq.' ') goto 10
         k=k+1
         string(k:k)=char
         lastch=char
   10 continue
      length=k
c
c     ----- now fill out with blanks -----
c
      string(k+1:)=' '
c
c
      return
      end
