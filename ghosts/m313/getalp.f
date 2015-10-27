*deck @(#)getalp.f	1.1  11/30/90
      subroutine getalp(ai,aj,ak,al,index,alpha,nv,len,camcb,ca,cb,
     #                  nij,nkl)
c
      implicit integer (a-z)
c
      real*8 ai(*),aj(*),ak(*),al(*),alpha(nv,4),camcb(nv,3)
      real*8 ca(nij,3),cb(nkl,3)
      integer index(len,6)
c
      do 1 i=1,nv
         alpha(i,1)=ai(index(i,1))
    1 continue
      do 2 i=1,nv
         alpha(i,2)=aj(index(i,2))
    2 continue
      do 3 i=1,nv
         alpha(i,3)=ak(index(i,3))
    3 continue
      do 4 i=1,nv
         alpha(i,4)=al(index(i,4))
    4 continue
c
      do 6 coord=1,3
         do 5 i=1,nv
            camcb(i,coord)=ca(index(i,5),coord)-cb(index(i,6),coord)
    5    continue
    6 continue
c
      return
      end
