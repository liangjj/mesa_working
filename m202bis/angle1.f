*deck @(#)angle1.f	1.2  8/19/91
      subroutine angle1(iout,nz,iz,ian,c)
c***begin prologue     angle1
c***date written       850601  yymmdd
c***revision date      yymmdd  yymmdd
c***keywords           bond angles
c***author             gauss82
c***source             @(#)angle1.f	1.2   8/19/91
c***purpose            computes and prints angles between adjacent bonds.
c                      used when a z-matrix is available.
c***description
c     call angle1(iout,nz,iz,ian,c)
c       iout    output file for printing.
c       nz      the number of z-matrix variables.
c       iz
c       ian     the atomic number vector(nz).
c       c       the cartesian coordinates array(3,nz).
c***references
c***routines called    fillel(util), asub(math), sdot(math), putatm(m202),
c                      putfp(chr)
c***end prologue       angle1
      implicit real*8(a-h,o-z)
c
      character el(106)*2, line*80, blank*1, equal*1
      integer iz(4,nz), ian(nz)
      real*8 c(3,nz), v1(3), v2(3)
      data maxel/104/, el(1)/'x ' /, iwidth/14/, idecim/4/,
     $     equal/'='/, blank/' '/, one/1.0d0/, f10/10.0d0/,
     $     f45/45.0d0/, f100/100.0d0/
 1000 format(1x,'interatomic angles:')
 1010 format(5x,80a1)
c
c     loop over pairs of bonds.
c
      if(nz.lt.3) return
      write(iout,1000)
      call fillel(0,maxel,el(2))
      ifield = iwidth + idecim + 6
      icur = 0
      line=blank
      do 40 i = 3, nz
          limj = i - 1
          do 40 j = 2, limj
              if(iz(1,j).eq.iz(1,i)) then
                  i1 = j
                  j1 = iz(1,i)
                  k1 = i
              else if(iz(1,j).eq.i) then
                  i1 = min(j,iz(1,i))
                  j1 = i
                  k1 = max(j,iz(1,i))
              else if(j.eq.iz(1,i)) then
                  i1 = iz(1,j)
                  j1 = j
                  k1 = i
              else
                  goto 40
              endif
              call asub(3,c(1,i1),c(1,j1),v1)
              call asub(3,c(1,k1),c(1,j1),v2)
              r1 = sqrt(sdot(3,v1,1,v1,1))
              r2 = sqrt(sdot(3,v2,1,v2,1))
              cosa = sdot(3,v1,1,v2,1) / (r1*r2)
              if(abs(cosa).gt.one) cosa = sign(one,cosa)
              angle = f45 * acos(cosa) / atan(one)
              call putatm(i1,el(ian(i1)+2),j1,el(ian(j1)+2),k1,
     $                    el(ian(k1)+2),iwidth,line,icur)
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
