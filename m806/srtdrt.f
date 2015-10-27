*deck @(#)srtdrt.f	5.1  11/6/94
      subroutine srtdrt(levpt,levpt1,levnr,levnr1,a,a1,b,b1,s,s1,
     #                  x,x1,arc,arc1,nlwks,nlwks1,nlevs,nrows,nsym)
c
      implicit integer (a-z)
c
      integer levpt(nlevs),levpt1(nlevs),levnr(nlevs),levnr1(nlevs)
      integer a(nrows),a1(nrows),b(nrows),b1(nrows),s(nrows),s1(nrows)
      integer x(nrows),x1(nrows),arc(4,nrows),arc1(4,nrows)
      integer nlwks(nrows),nlwks1(nrows)
c
      do 1 i=1,nlevs
         levpt(i)=levpt1(i)
    1 continue
      do 2 i=1,nlevs
         levnr(i)=levnr1(i)
    2 continue
      do 20 level=nlevs,1,-1
         pt=levpt(level)
         nr=levnr(level)
         if (level.lt.nlevs) then
            ptp1=levpt(level+1)
            nrp1=levnr(level+1)
         end if
         amax=0
         amin=9999
         bmax=0
         bmin=9999
         do 3 row=pt+1,pt+nr
            amax=max(amax,a1(row))
            amin=min(amin,a1(row))
            bmax=max(bmax,b1(row))
            bmin=min(bmin,b1(row))
    3    continue
         newrow=pt
         do 19 sym=nsym-1,0,-1
            do 18 ajunk=amax,amin,-1
               do 17 bjunk=bmax,bmin,-1
                  do 16 row=pt+1,pt+nr
                     if (a1(row).eq.ajunk.and.b1(row).eq.bjunk.and.
     #                   s1(row).eq.sym) then
                        newrow=newrow+1
                        a(newrow)=a1(row)
                        b(newrow)=b1(row)
                        s(newrow)=s1(row)
                        x(newrow)=x1(row)
                        arc(1,newrow)=arc1(1,row)
                        arc(2,newrow)=arc1(2,row)
                        arc(3,newrow)=arc1(3,row)
                        arc(4,newrow)=arc1(4,row)
                        nlwks(newrow)=nlwks1(row)
                        if (level.lt.nlevs) then
                           do 10 rowp1=ptp1+1,ptp1+nrp1
                              do 9 case=1,4
                                 if (arc(case,rowp1).eq.row-pt)
     #                              arc(case,rowp1)=-newrow+pt
    9                         continue
   10                      continue
                        end if
                     end if
   16             continue
   17          continue
   18       continue
   19    continue
   20 continue
c
      do 22 case=1,4
         do 21 row=1,nrows
            arc(case,row)=-arc(case,row)
   21    continue
   22 continue
c
c
      return
      end
