*deck @(#)initx2.f	5.1  11/6/94
      subroutine initx2(int,int1,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,
     $                  ss,ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww)
c
c
c
      implicit real*8 (a-h,o-z)
      integer arr,tr1,tr2,asm,aos,os,wtw,wtx,wty,wab,ss,symorb
      integer bmax,orbfrm
      real*8 int(nmax),int1(norbs*norbs),c(nwks),s(nwks)
      real*8 ci(1),cj(1),si(1),sj(1),h1(1),h2(1)
      integer ijxx(numij),nklxx(nsym),ijww(numij),nklww(nsym)
      dimension kadd(symorb),ladd(symorb),ijadd(numij),wtw(orbfrm,nsym)
      dimension wtx(orbfrm,nsym),wty(orbfrm),wab(orbfrm),ss(norbs)
      common /dims/ nbf,nsym,norbs,nrows,nrows4,nwks,nwks2,nlevs
     *,             nrowoc,nrow4o,nwksoc,nlevoc,norboc,levfrm
     *,             nwksmx,nlwkmx,nuwkmx,bmax,nroots,orbfrm
      common /optns/  iguess,irstrt,irooti,irootf,i34x
      common /ints/   nmax,nmax2,ngroup,nblkoc,numij,symorb,intsrt
c  universal identity of the objects in these common
c     common /all/ val1,val2,val3,arr,itr1,itr2,ia,ja,itype,isegt
c    *,lvfrm1,nlwki,nlwkj,imax,imin
c     common /symm/ jsm,jfsym,ifsym,maxsym(8),minsym(8),ismoff(8)
c    #,             numsym(8)
      common /symm/ asm,js,is,mx(8),mn(8),os(8),numsym(8)
      common /tapes/itap20,itape5,itape6,itape8,itap12,itap03,itap04
     *,             itape3,itap05,itap06
      common /all/val1,val2,val3,arr,tr1,tr2,ia,ja,m,iseg,n,n1,n2
     *,           imax,imin
      common /count/  icount,ixx4,ixx5,ixx6,ixx8,iww4,iww5,iww6,iww8
     #,               iwx7,ixw9,ixy3,ixy16,ixy18,ixy22,iwy2,iwy15,iwy17
     #,               iwy21,iyx20,iyw19
      common /count2/ ientry(20),time(20)
      real*8 val1,val2,val3
      sqrt2=sqrt(2.0d+00)
      sqrt3=sqrt(3.0d+00)
      sqt1p5=sqrt(1.5d+00)
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
cbl   itest=0
cbl   return
c
c--------------------------------------------------------external
c--------------------------------------------------------fourx---
c
      entry fourx(int,int1,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $            ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww)
      entry external(int,int1,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $               ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww)
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
c     if (i34x.ne.0.and.iseg.eq.3) go to 1
      if((iseg.eq.3).or.(iseg.eq.16).or.(iseg.eq.18))go to 1
      call xy(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $        ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww)
      go to 9
  304 continue
c     ----- yy -----
      go to 1
  305 continue
c     ----- xx -----
      call xx(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $        ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww)
      go to 9
  306 continue
c     ----- zy -----
      go to 1
  307 continue
c     ----- yx -----
      call yx(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $        ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww)
      go to 9
  308 continue
c     ----- yw -----
      call yw(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $        ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww)
      go to 9
  309 continue
c     ----- wz -----
      go to 1
  310 continue
c     ----- wy -----
c     if (i34x.ne.0.and.iseg.eq.2) go to 1
      if (iseg.eq.2) go to 1
      if((iseg.eq.15).or.(iseg.eq.17))go to 1
      call wy(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $        ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww)
      go to 9
  311 continue
c     ----- wx -----
      call wx(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $        ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww)
      go to 9
  312 continue
c     ----- xw -----
      call xw(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $        ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww)
      go to 9
  313 continue
c     ----- ww -----
      call ww(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $        ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww)
      go to 9
  314 continue
c     ----- yy 4x -----
      call shape4(int,int1,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
      go to 9
  315 continue
c     ----- xx 4x -----
c     if (i34x.eq.0) then
c        call xx4x(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
c    $        ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww)
c     else
      call shape4(int,int1,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
c     end if
      go to 9
  316 continue
c     ----- ww 4x -----
c     if (i34x.eq.0) then
c        call ww4x(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
c    $        ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww)
c     else
      call shape4(int,int1,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
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
c..rlm     call shapes
      call shapes(int,int1,c,s,ijadd,kadd,ladd,wtw,wtx,wty,
     $            wab,ss)
    9 continue
      return
      end
