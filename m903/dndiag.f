*deck @(#)dndiag.f	5.1  11/6/94
      subroutine dndiag(dnarc,dnwt,a,b,dnnwks,nrows,nlevs,cfs,
     #                  acoef,bcoef,irowsv,iwtsv,pagesv,segsv,
     #                  diag,nwks,h,g,nnp,traksv)
c
      implicit integer (a-z)
c
c     ----- external arrays modified -----
c
      real*8 diag(nwks)
c
c     ----- external arrays not modified -----
c
      real*8 cfs(420),h(nnp),g(nnp,nnp)
      integer dnarc(4,nrows),dnwt(4,nrows),a(nrows),b(nrows)
      integer dnnwks(nrows)
c
c     ----- external arrays used for scratch -----
c
      real*8 acoef(nlevs),bcoef(nlevs)
      integer irowsv(nlevs),iwtsv(nlevs),traksv(nlevs)
      integer pagesv(nlevs),segsv(nlevs)
c
c     ----- local arrays -----
c
      integer segmax(3),segmin(3),nxtpag(20),iarc(20),trak(20)
c
c     ----- local scalars -----
c
      real*8 loop,crite,root2,rootn2,toor2,toorn2,d,dx
c
c     ----- common blocks -----
c
      common /io/ inp,iout
c
      equivalence (segmin(2),segmax(1))
c
      data segmin(1) /0/
      data segmax / 6, 13, 20 /
      data nxtpag
     #  /  0, 3, 0, 3, 0, 2,
     #     2, 0, 2, 0, 2, 0, 2,
     #     3, 0, 3, 0, 3, 0, 3/
      data   iarc
     #  /  2, 2, 3, 3, 4, 4,
     #     1, 2, 2, 3, 3, 4, 4,
     #     1, 2, 2, 3, 3, 4, 4/
      data   trak
     #  /  3, 1, 3, 1, 2, 1,
     #     0, 0, 0, 0, 0, 0, 0,
     #     0, 0, 0, 0, 0, 0, 0/
c
      save segmin,segmax,trak,nxtpag,iarc
c
      crite=1.0d-05
      root2=sqrt(2.0d+00)
      rootn2=-root2
      toor2=1.0d+00/root2
      toorn2=-toor2
c
c     ----- descend partial walks to the opening level -----
c
      do 1000 top=nlevs,2,-1
         if (top.eq.nlevs) then
c
c           ----- special situation for loops opening at top of graph
c
            harc=4
            hlev=nlevs
            level=top
            irow=1
            iwtsv(top)=0
            acoef(top)=1.0d+00
            page=1
            seg=segmin(page)
            segmx=segmax(page)
            go to 201
         end if
c
c        ----- initialize head section ------
c
         hrow=1
         hlev=nlevs
         iwtsv(hlev)=0
         harc=0
c
c        ----- start tree search of head sections -----
c
  101    continue
            harc=harc+1
            if (harc.gt.4) then
               hlev=hlev+1
               if (hlev.gt.nlevs) go to 1000
               hrow=irowsv(hlev-1)
               harc=segsv(hlev-1)
               go to 101
            end if
c
c           ----- check that this segment exists on the graph -----
c
            hnxt=dnarc(harc,hrow)
            if (hnxt.eq.0) go to 101
c
c           ----- hop out of search if descended to desired level -----
c
            if (hlev-1.eq.top) go to 102
c
c           ----- otherwise, descend another level -----
c
            hlev=hlev-1
            irowsv(hlev)=hrow
            segsv(hlev)=harc
            iwtsv(hlev)=iwtsv(hlev+1)+dnwt(harc,hrow)
            harc=0
            hrow=hnxt
         go to 101
c
c     ----- initialize stacks for tree search -----
c
  102    continue
         level=top
         irow=hnxt
         iwtsv(level)=iwtsv(level+1)+dnwt(harc,hrow)
         acoef(level)=1.0d+00
         page=1
         seg=segmin(page)
         segmx=segmax(page)
c
c        ----- start tree search in downward direction until we have
c              found all loops
c
  201    continue
c
            seg=seg+1
            if (seg.gt.segmx) then
c
c              ----- exhausted this page of the table, so back up a level
c
               level=level+1
               if (level.gt.top) then
                  go to 101
               end if
c
               seg=segsv(level-1)
               page=pagesv(level-1)
               segmx=segmax(page)
               irow=irowsv(level-1)
               go to 201
            end if
c
c           ----- check that this segment is valid on the graph -----
c
            inxt=dnarc(iarc(seg),irow)
            if (inxt.eq.0) go to 201
c
c           ----- and set the track value -----
c
            if (trak(seg).ne.0) then
               traksv(level-1)=trak(seg)
            else
               traksv(level-1)=traksv(level)
            end if
c
c           ----- evaluate this segment -----
c
            go to ( 1, 44,  1, 45, 50, 51,  1, 77,  1, 77,
     #              1, 78,  1, 71, 87, 75, 68, 76, 82, 71), seg
c
    1          acoef(level-1)=acoef(level)
               goto 120
   44          junk1=b(irow)+221
               acoef(level-1) = acoef(level) * toor2
               bcoef(level-1) = acoef(level) * cfs(junk1)
               go to 120
   45          junk1 = b(irow) + 163
               acoef(level-1) = acoef(level) * toor2
               bcoef(level-1) = acoef(level) * cfs(junk1)
               go to 120
   50          acoef(level-1) = acoef(level) + acoef(level)
               d=0.5d+00
               go to 120
   51          acoef(level-1)=acoef(level)*root2
               go to 120
   68          junk1 = b(irow) + 222
               dx=acoef(level)*toorn2
               d=dx+bcoef(level)*cfs(junk1)
               if(abs(d).lt.crite) go to 111
               acoef(level-1) = d
               d=-(dx+dx)/d
               go to 120
   71          acoef(level-1) = acoef(level)
               bcoef(level-1) = bcoef(level)
               go to 120
   75          junk1 = b(irow) + 302
               acoef(level-1) = acoef(level)
               bcoef(level-1) = bcoef(level) * cfs(junk1)
               go to 120
   76          junk1 = b(irow) + 303
               acoef(level-1) = acoef(level)
               bcoef(level-1) = bcoef(level) * cfs(junk1)
               go to 120
   77          acoef(level-1)=acoef(level)*toorn2
               d=-2.0d+00
               go to 120
   78          acoef(level-1)=acoef(level)*rootn2
               d=-2.0d+00
               go to 120
   82          acoef(level-1) = acoef(level) * rootn2
               d=-2.0d+00
               go to 120
   87          junk1 = b(irow) + 162
               dx=acoef(level)*toorn2
               d=dx+bcoef(level)*cfs(junk1)
               if(abs(d).lt.crite) go to 111
               acoef(level-1)=d
               d=-(dx+dx)/d
               go to 120
  111          continue
               traksv(level-1)=4
               acoef(level-1)=-(dx+dx)
  120          continue
c
c           ----- update stacks -----
c
            newpag=nxtpag(seg)
            if (newpag.eq.0) then
c
c              ----- have finished a loop -----
c
               j=level-1
               i=top-1
c
ctemp            if (orbsym(j).ne.is) go to 201
c
               ii=i*(i-1)/2+i
               ij=i*(i-1)/2+j
               jj=j*(j-1)/2+j
c
               go to (211,212,213,214), traksv(level-1)
                  call lnkerr('bad track in dndiag')
c
  211             continue
                     loop=acoef(level-1)*(g(ij,ij)+d*g(ii,jj))
                     go to 219
  212             continue
                     loop=acoef(level-1)*(h(ii)+d*g(ii,ii))
                     go to 219
  213             continue
                     loop=acoef(level-1)*h(ii)
                     go to 219
  214             continue
                     loop=acoef(level-1)*g(ii,jj)
  219             continue
c
               iwt=iwtsv(level)+dnwt(iarc(seg),irow)
c
               do 220 iwk=iwt+1,iwt+dnnwks(inxt)
                  diag(iwk)=diag(iwk)+loop
  220          continue
c
c              ----- go back up and try next segment -----
c
               go to 201
            end if
c
c           ----- update stacks and then go down one level -----
c
            if (level.le.2) go to 201
c
            irowsv(level-1)=irow
            pagesv(level-1)=page
            iwtsv(level-1)=iwtsv(level)+dnwt(iarc(seg),irow)
            segsv(level-1)=seg
            irow=inxt
            page=newpag
c
            level=level-1
            seg=segmin(page)
            segmx=segmax(page)
         go to 201
 1000 continue
c
c     ----- the return is down here for convenience -----
c
 9000 continue
c
      return
      end
