*deck @(#)trmat.f	5.1  11/6/94
      subroutine trmat(bftran,t,ptbftr,maxmom,nop,ncart,nx,ny,nz)
c
      implicit integer (a-z)
c
      integer ptbftr(0:maxmom)
      integer ncart(0:maxmom)
      integer nx(*),ny(*),nz(*)
      integer dxyz(6,2),fxyz(10,2)
      integer gxyz(15,2)
      real*8 bftran(*)
      real*8 cosine(3,3),origin(3),vector(3),t1(3),t2(3)
      real*8 t(3,3,nop),norm,temp
c
c     the values which define the order of the cartesian components
c     are reproduced here.  note that they are set in genbas(m102)
c     and if that routine is changed, so should the arrays dxyz,etc.
c
c
c     data nx/0, 1,0,0, 2,0,0,1,1,0, 3,0,0,1,2,2,1,0,0,1,
c    $        4,3,2,1,0,3,2,1,0,2,1,0,1,0,0/
c     data ny/0, 0,1,0, 0,2,0,1,0,1, 0,3,0,2,1,0,0,1,2,1,
c    $        0,1,2,3,4,0,1,2,3,0,1,2,0,1,0/
c     data nz/0, 0,0,1, 0,0,2,0,1,1, 0,0,3,0,0,1,2,2,1,1,
c    $        0,0,0,0,0,1,1,1,1,2,2,2,3,3,4/
c
c
c                 x  y  z  x  x  y
c                 x  y  z  y  z  z
c
      data dxyz / 1, 2, 3, 1, 1, 2,
     #            1, 2, 3, 2, 3, 3 /
c
c                 x  y  z  x  x  x  x  y  y  x
c                 x  y  z  y  x  x  z  z  y  y
c                 x  y  z  y  y  z  z  z  z  z
c
      data fxyz / 1, 2, 3, 4, 1, 1, 5, 6, 2, 4,
     #            1, 2, 3, 2, 2, 3, 3, 3, 3, 3 /
c
c                 x  x  x  x  y  x  x  x  y  x  x  y  x  y  z
c                 x  x  x  y  y  x  x  y  y  x  y  y  z  z  z
c                 x  x  y  y  y  x  y  y  y  z  z  z  z  z  z
c                 x  y  y  y  y  z  z  z  z  z  z  z  z  z  z
c
      data gxyz / 1, 1, 5, 4, 2, 1, 5, 4, 2, 6,10, 9, 7, 8, 3,
     $            1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 /
c
c     ----- pointers into the nx,ny,nz array
c
      data poff/1/, doff/4/, foff/10/, goff/20/
      save dxyz,fxyz,gxyz
      save poff,doff,foff,goff
c
c     ----- statement functions to address elements of p, d, f,...
c           transformation matrices
c
      ppt(i,j)=ptbftr(1)+(i-1)+ncart(1)*(j-1)+ncart(1)**2*(op-1)
      dpt(i,j)=ptbftr(2)+(i-1)+ncart(2)*(j-1)+ncart(2)**2*(op-1)
      fpt(i,j)=ptbftr(3)+(i-1)+ncart(3)*(j-1)+ncart(3)**2*(op-1)
      gpt(i,j)=ptbftr(4)+(i-1)+ncart(4)*(j-1)+ncart(4)**2*(op-1)
      hpt(i,j)=ptbftr(5)+(i-1)+ncart(5)*(j-1)+ncart(5)**2*(op-1)
      ipt(i,j)=ptbftr(6)+(i-1)+ncart(6)*(j-1)+ncart(6)**2*(op-1)
c
c     ----- cosine is the rotation matrix to go from the standard
c     orientation to the local frame, origin is the location of
c     the origin of the local frame
c
      call runit(cosine,3)
      call rzero(origin,3)
c
c     ----- fill out pointers to the transformation matrices -----
c
      ptbftr(0)=1
      do 1 angmom=1,maxmom
         ptbftr(angmom)=ptbftr(angmom-1)+ncart(angmom-1)**2*nop
 1    continue
c
c     initialize the transformation matrices
      call rzero(bftran,ptbftr(maxmom)+ncart(maxmom)**2*nop)
c
c
c     ----- calculate the transformation matrices for the p functions
c     (which are the same as for the three axis vectors)
c
      do 2 op=1,nop
         bftran(op)=1.0d+00
 2    continue
c
      if (maxmom.lt.1) return
c
      do 20 xyz=1,3
         call vmove(vector,origin,3)
         vector(xyz)=vector(xyz)+1.0d+00
c
c     ----- rotate and translate to the local frame -----
c
         call tolocl(t1,vector,cosine,origin)
c
c     ----- operate with the operations of the group -----
c
         do 10 op=1,nop
            call ebc(t2,t1,t(1,1,op),1,3,3)
c
c     ----- and transform back to master frame -----
c
            call tomast(bftran(ppt(1,xyz)),t2,cosine,origin)
c
            do 5 i=1,3
               bftran(ppt(i,xyz))=bftran(ppt(i,xyz))-origin(i)
    5       continue
 10      continue
 20   continue
c
c     ----- produce the d function transformations from p's -----
c
      if (maxmom.lt.2) return
c
      do 50 op=1,nop
         do 40 i=1,6
            i1=dxyz(i,1)
            i2=dxyz(i,2)
c           these two functions, i1 and i2, will produce the transformation
c           for function i.
c           loop over the columns of each of the "generating" functions.
            do 30 j=1,3
               do 25 k=1,3
                  temp=bftran(ppt(i1,j))*bftran(ppt(i2,k))
c
c                 figure out which column of the d-matrix this
c                 corresponds to by determining the polynomial powers
c                 and finding it in the list.
                  ntot=nx(poff+j) +nx(poff+k)
                  ltot=ny(poff+j) +ny(poff+k)
                  mtot=nz(poff+j) +nz(poff+k)
                  do 23 mu=1,6
                     if(ntot.eq.nx(doff+mu).and.
     $                  ltot.eq.ny(doff+mu).and.
     $                  mtot.eq.nz(doff+mu)) goto 24
   23             continue
                  call lnkerr('trmat: problem finding polynomial')
   24             continue
c                 mu is the index of the appropriate column in the d-matrix.
                  bftran(dpt(i,mu))=bftran(dpt(i,mu))+temp
   25          continue
   30       continue
   40    continue
   50 continue
c
c     ---- and the f from d and p transformations -----
c
      if (maxmom.lt.3) return
c
ctemp
c     the normalization stuff remains to be worked out
cend
      do 80 op=1,nop
         do 70 i=1,10
            i1=fxyz(i,1)
            i2=fxyz(i,2)
            do 60 j=1,6
               do 59 k=1,3
                  temp=bftran(dpt(i1,j))*bftran(ppt(i2,k))
                  ntot=nx(doff+j) +nx(poff+k)
                  ltot=ny(doff+j) +ny(poff+k)
                  mtot=nz(doff+j) +nz(poff+k)
                  do 57 mu=1,10
                     if(ntot.eq.nx(foff+mu).and.
     $                  ltot.eq.ny(foff+mu).and.
     $                  mtot.eq.nz(foff+mu)) go to 58
   57             continue
                  call lnkerr('trmat: problem finding polynomial')
   58             continue
                  bftran(fpt(i,mu))=bftran(fpt(i,mu))+temp
   59          continue
   60       continue
   70    continue
   80 continue
c
c     ---- and the g from f and p transformations -----
c
      if (maxmom.lt.4) return
c
ctemp
c     the normalization stuff remains to be worked out
cend
      do 110 op=1,nop
         do 100 i=1,15
            i1=gxyz(i,1)
            i2=gxyz(i,2)
            do 90 j=1,10
               do 89 k=1,3
                  temp=bftran(fpt(i1,j))*bftran(ppt(i2,k))
                  ntot=nx(foff+j) +nx(poff+k)
                  ltot=ny(foff+j) +ny(poff+k)
                  mtot=nz(foff+j) +nz(poff+k)
                  do 85 mu=1,15
                     if(ntot.eq.nx(goff+mu).and.
     $                  ltot.eq.ny(goff+mu).and.
     $                  mtot.eq.nz(goff+mu)) go to 86
   85             continue
                  call lnkerr('trmat: problem finding polynomial')
   86             continue
                  bftran(gpt(i,mu))=bftran(gpt(i,mu)) +temp
   89          continue
   90       continue
  100    continue
  110 continue
c
c
      if(maxmom.ge.5) then
         call lnkerr('symmetry cannnot handle h-functions.')
      end if
c
c

      return
      end
