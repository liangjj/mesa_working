*deck @(#)initx2.f	5.1  11/6/94
      subroutine initx2(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
      implicit real*8 (a-h,o-z)
c
      integer and
      integer arr,tr1,tr2,asm,aos,os,wtw,wtx,wty,wab,ss,symorb
      integer bmax,orbfrm
      integer ijxx(numij ),nklxx(nsym,orbfrm),ijww(numij )
      integer nklww(nsym,orbfrm),klxx(1),klww(1)
      dimension kadd(symorb),ladd(symorb),ijadd(numij),wtw(orbfrm,nsym)
      dimension wtx(orbfrm,nsym),wty(orbfrm),wab(orbfrm),ss(norbs)
      real*8 ci(1),cj(1),si(1),sj(1),h1(1),h2(1)
      real*8 int(nmax),c(nwks),s(nwks)
      real*8 val1,val2,val3
c
      common /dims/ nbf,nsym,norbs,nrows,nrows4,nwks,nwks2,nlevs
     *,             nrowoc,levfrm
     *,             nwksmx,nlwkmx,nuwkmx,bmax,nroots,orbfrm
      common /optns/  iguess,irstrt,irooti,irootf,i34x
      common /ints/   nmax,nmax2,ngroup,nblkoc,numij,symorb,intsrt
c     universal identity of the objects in these common
c     common /all/ val1,val2,val3,arr,itr1,itr2,ia,ja,itype,isegt
c    *,lvfrm1,nlwki,nlwkj,imax,imin
c     common /symq/ jsm,jfsym,ifsym,maxsym(8),minsym(8),ismoff(8)
c    #,             numsym(8)
      common /symq/ asm,js,is,mx(8),mn(8),os(8),numsym(8)
      common /io/     itape5,itape6
      common /all/val1,val2,val3,arr,tr1,tr2,ia,ja,m,iseg,n,n1,n2
     *,           imax,imin
      common /count/  icount,ixx4,ixx5,ixx6,ixx8,iww4,iww5,iww6,iww8
     #,               iwx7,ixw9,ixy3,ixy16,ixy18,ixy22,iwy2,iwy15,iwy17
     #,               iwy21,iyx20,iyw19
      common /count2/ ientry(20),time(20)
c
      icount=0
      ixx4=0
      ixx5=0
      ixx6=0
      ixx8=0
      iww4=0
      iww5=0
      iww6=0
      iww8=0
      iwx7=0
      ixw9=0
      return
c
c--------------------------------------------------------external
c--------------------------------------------------------fourx---
c
      entry external(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss
     #,                 ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww
     #,                 klxx,klww)
      entry fourx(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss
     #,                 ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww
     #,                 klxx,klww)
c
c
      if (i34x.ne.0) then
         do 11 iq=1,nsym
            mins=mn(iq)
            if (mins.gt.n) then
               numsym(iq)=0
            else
               maxs=mx(iq)
               if (maxs.le.n) then
                  numsym(iq)=maxs-mins+1
               else
                  numsym(iq)=n-mins+1
               end if
            end if
   11    continue
      end if
c
c
      ientry(m)=ientry(m)+1
      icount=icount+1
c
      aos=os(asm+1)
      go to (313,305,304,312,311,310,307,308,303,306,309,302,301,314
c     go to (301,302,303,304,305,306,307,308,309,310,311,312,313,314
c     go to ( yz, xz, xy, yy, xx, zy, yx, yw, wz, wy, wx, xw, ww,---yy
c
     1,315,316,317,318,319),m
c       xx, ww, zz,---4 internals),m
c
      stop 'funny m'
  301 continue
c     ----- yz -----
      go to 1
  302 continue
c     ----- xz -----
      go to 1
  303 continue
c     ----- xy -----
      if (and(i34x,1).ne.0.and.iseg.eq.3) go to 1
      if (and(i34x,1).eq.0.and.n.lt.orbfrm.and.iseg.eq.18) go to 9
      call xy(int,c,s,ijxx,nklxx,h1,ci,si,wtx,ladd)
      go to 9
  304 continue
c     ----- yy -----
      go to 1
  305 continue
c     ----- xx -----
      call xx(int,c,s,h1,h2,ci,cj,si,sj,kadd,wtx)
      go to 9
  306 continue
c     ----- zy -----
      go to 1
  307 continue
c     ----- yx -----
      call yx(int,c,s,h1,ci,si,wtx,ladd)
      go to 9
  308 continue
c     ----- yw -----
      call yw(int,c,s,h1,ci,si,wtw,wab,ladd)
      go to 9
  309 continue
c     ----- wz -----
      go to 1
  310 continue
c     ----- wy -----
      if (and(i34x,2).ne.0.and.iseg.eq.2) go to 1
      if (and(i34x,2).eq.0.and.n.lt.orbfrm.and.iseg.eq.17) go to 9
      call wy(int,c,s,ijww,nklww,h1,wtw,wab,ci,si,ladd)
      go to 9
  311 continue
c     ----- wx -----
      call wx(int,c,s,kadd,wtw,wab,wtx,h1,h2,ci,cj,si,sj)
      go to 9
  312 continue
c     ----- xw -----
      call xw(int,c,s,kadd,wtw,wab,wtx,h1,h2,ci,cj,si,sj)
      go to 9
  313 continue
c     ----- ww -----
      call ww(int,c,s,kadd,h1,h2,ci,cj,si,sj,wtw,wab)
      go to 9
  314 continue
c     ----- yy 4x -----
cpws  call shape4
      call shape4(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
      go to 9
  315 continue
c     ----- xx 4x -----
c     if (i34x.eq.0) then
      call xx4x(int,c,s,ijxx,nklxx,wtx,ss)
c     else
c        call shape4
c     end if
      go to 9
  316 continue
c     ----- ww 4x -----
c     if (i34x.eq.0) then
      call ww4x(int,c,s,ijww,nklww,wtw,wab,ss)
c     else
c        call shape4
c     end if
      go to 9
  317 continue
c     ----- zz entries ----
      go to 1
  318 continue
c     ----- 4 internal -----
      go to 1
  319 continue
c     ----- 4 internal, (3,2,1) case -----
      go to 1
c
c
c
    1 continue
      call shapes(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
    9 continue
cpsctemp
cps      if (npt.ge.2*998) call lnkerr('end of iteration')
cps      npt=npt+1
cps      if (npt.lt.999) go to 3004
cps      if (qtemp.eq.s(26)) go to 3004
cps      write (itape6,3002) npt,m,arr,val1,val2,val3,ia,ja,asm,is,js,s(26)
cps 3002 format (1x,3i5,3f10.6,5i3,f10.6)
cpscps      write (itape6,3003) s
 3003 format (1x,10f10.6)
 3004 continue
      qtemp=s(26)
cend
      return
      end