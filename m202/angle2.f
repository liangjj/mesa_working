*deck @(#)angle2.f	5.1  11/6/94
      subroutine angle2(iout,n,ian,c)
c***begin prologue     angle2.f
c***date written       850601  yymmdd
c***revision date      11/6/94
c***keywords           bond angles, cartesian coordinates
c***author             binkley, et al., (g82)
c***source             @(#)angle2.f	5.1   11/6/94
c***purpose            computes and prints bond angles between atoms.
c                      this is used when no z-matrix is available.
c***description
c     call angle2(iout,n,ian,c)
c       iout    output file for printing.
c       n       the number of atoms.
c       ian     the atomic number vector(n).
c       c       the cartesian coordinates array(3,n).
c***references
c***routines called    fillel(util), vsub(math), sdot(math), putatm(m202),
c                      putfp(chr)
c***end prologue       angle2.f
      implicit none
c     --- input variables -----
      integer iout,n
c     --- input arrays (unmodified) ---
      integer ian(n)
      real*8 c(3,n)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer maxel,iwidth,idecim
      integer i,j,k,i1,j1,k1
      integer ifield,icur,limj,limk,len,icurm1,jj
      character el(106)*2, line*80, blank*1, equal*1
      real*8 cut,one,f10,f45,f100
      real*8 rij,rik,rjk,r1,r2,cosa,angle
      real*8 sdot,sqrt
      real*8 v1(3), v2(3)
c
      data maxel/104/, el(1)/'x '/, iwidth/14/, idecim/4/, equal/'='/
     $     blank/' '/, cut/5.7d0/, one/1.0d0/, f10/10.0d0/,
     $     f45/45.0d0/, f100/100.0d0/
      save maxel,el,iwidth,idecim,equal
      save blank,cut,one,f10,f45,f100
c
 1000 format(1x,'interatomic angles:')
 1010 format(5x,80a1)
c
c     --- loop over triples of atoms within cut of each other.
      if(n.lt.3) return
c
      write(iout,1000)
      call fillel(0,maxel,el(2))
      ifield = iwidth + idecim + 6
      icur = 0
      line=blank
      do 40 i = 3, n
          limj = i - 1
          do 40 j = 2, limj
              limk = j - 1
              do 40 k = 1, limk
                  call vsub(v1,c(1,i),c(1,j),3)
                  rij = sqrt(sdot(3,v1,1,v1,1))
                  call vsub(v1,c(1,i),c(1,k),3)
                  rik = sqrt(sdot(3,v1,1,v1,1))
                  call vsub(v1,c(1,j),c(1,k),3)
                  rjk = sqrt(sdot(3,v1,1,v1,1))
                  if(rij.ge.rjk.and.rij.ge.rik) then
                      i1 = j
                      j1 = k
                      k1 = i
                  else if(rik.ge.rjk) then
                      i1 = k
                      j1 = j
                      k1 = i
                  else
                      i1 = k
                      j1 = i
                      k1 = j
                  endif
                  call vsub(v1,c(1,i1),c(1,j1),3)
                  call vsub(v2,c(1,k1),c(1,j1),3)
                  r1 = sqrt(sdot(3,v1,1,v1,1))
                  r2 = sqrt(sdot(3,v2,1,v2,1))
                  if(r1.gt.cut.or.r2.gt.cut) goto 40
                  cosa = sdot(3,v1,1,v2,1) / (r1*r2)
                  if(abs(cosa).gt.one) cosa = sign(one,cosa)
                  angle = f45 * acos(cosa) / atan(one)
                  call putatm(i1,el(ian(i1)+2),j1,el(ian(j1)+2),k1,
     $                        el(ian(k1)+2),iwidth,line,icur)
                  icur=icur+1
                  line(icur:icur)=equal
                  if(angle.lt.f100) icur=icur+1
                  if(angle.lt.f10)  icur=icur+1
                  call putfp(angle,idecim,line,icur)
                  len = mod(icur,ifield)
                  if(len.ne.0) then
                      len = ifield - len
                      icur=icur+len
                  endif
                  if((icur+ifield).gt.72) then
                      icurm1=icur-1
                      write(iout,1010) (line(jj:jj),jj=1,icurm1)
                      icur = 0
                      line=blank
                  endif
   40 continue
      icurm1=icur-1
      if(icurm1.gt.0) write(iout,1010) (line(i:i),i=1,icurm1)
c
c
      return
      end
