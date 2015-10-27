*deck @(#)movdrt.f	5.1  11/6/94
      subroutine movdrt(levpt,levpt1,levnr,levnr1,a,a1,b,b1,s,s1,
     #                  x,x1,arc,arc1,nlwks,nlwks1,nlevs,nrows)
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
      do 3 i=1,nrows
         a(i)=a1(i)
    3 continue
      do 4 i=1,nrows
         b(i)=b1(i)
    4 continue
      do 5 i=1,nrows
         s(i)=s1(i)
    5 continue
      do 6 i=1,nrows
        x(i)=x1(i)
    6 continue
      do 7 i=1,nrows
         arc(1,i)=arc1(1,i)
         arc(2,i)=arc1(2,i)
         arc(3,i)=arc1(3,i)
         arc(4,i)=arc1(4,i)
    7 continue
      do 8 i=1,nrows
         nlwks(i)=nlwks1(i)
    8 continue
c
      return
      end
