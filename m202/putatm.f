*deck @(#)putatm.f	5.1  11/6/94
      subroutine putatm(i,iel,j,jel,k,kel,width,string,cursor)
c***begin prologue     putatm.f
c***date written       850601  yymmdd
c***revision date      11/6/94
c***keywords           formatting
c***author             binkley, et al. (g82)
c                      martin, richard (lanl)
c***source             @(#)putatm.f	5.1   11/6/94
c***purpose            formats information about two atoms in a character
c                      string for printing.
c***description
c     call putatm(i,iel,j,jel,k,kel,width,string,cursor)
c***references
c***routines called    (none)
c***end prologue       putatm.f
      implicit none
c     --- input variables -----
      integer i,j,k,width,cursor
c     --- input arrays (unmodified) ---
      character*(*)iel,jel,kel,string
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer icur
      character substr*24, ci*16, cj*16, ck*16, itoc*16
c
c
      substr=iel
      icur=index(substr,' ')
      ci=itoc(i)
      substr(icur:)=ci
      icur=index(substr,' ')
      substr(icur:)='-'
      icur=icur+1
c
      substr(icur:)=jel
      icur=index(substr,' ')
      cj=itoc(j)
      substr(icur:)=cj
      icur=index(substr,' ')
      substr(icur:)='-'
      icur=icur+1
c
      substr(icur:)=kel
      icur=index(substr,' ')
      ck=itoc(k)
      substr(icur:)=ck
      icur=index(substr,' ')
c
c     --- right justify it into string.
      string(cursor+width-icur+1:cursor+width)=substr(1:icur)
      cursor=cursor+width
c
c
      return
      end
