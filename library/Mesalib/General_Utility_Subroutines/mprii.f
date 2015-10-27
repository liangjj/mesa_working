*deck mprii
      subroutine mprii(a,rowv,colv,n,m,ia,ja,title,rowt,colt,iout)
      implicit integer (a-z)
      character *80 title
      character *8 rowt, colt
      real *8 a
      dimension a(ia,ja), rowv(*), colv(*)
      write(iout,1) title
    1 format(a80)
      ibeg=0
      do 10 i=1,m,5
         iend=min(ibeg+5,m)
         ibeg=ibeg+1
         write (iout,30) (ii,ii=ibeg,iend)
         if (colv(1).ne.-99) then
             write (iout,20) colt,(colv(ii),ii=ibeg,iend)
         else
             write (iout,200) colt
         endif 
         if (rowv(1).ne.-99) then
             write (iout,50) rowt,rowv(1),(a(1,k),k=ibeg,iend)
         else
             write (iout,500) rowt,(a(1,k),k=ibeg,iend)
         endif            
         if (n.gt.1) then
             if (rowv(1).ne.-99) then
                 do 40 j=2,n
                    write (iout,60) rowv(j),(a(j,k),k=ibeg,iend)
   40            continue
             else
                 do 400 j=2,n
                    write (iout,600) (a(j,k),k=ibeg,iend)
  400            continue
             endif
         endif
         ibeg=iend
   10 continue
   20 format(a8,10x,5(3x,i4,3x))
   30 format(/,18x,5(3x,i4,3x))
  200 format(a8)
   50 format(a8,i10,5f10.5)
  500 format(a8,10x,5f10.5)
   60 format(8x,i10,5f10.5)
  600 format(18x,5f10.5)
      return
      end
