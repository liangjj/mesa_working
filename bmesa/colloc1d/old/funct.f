*deck funct.f
c***begin prologue     funct
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           functions, derivatives
c***author             schneider, barry(nsf)
c***source             
c***purpose            calculate functions and their first and second 
c***                   derivatives on grid
c***
c***
c***references
c
c***routines called
c***end prologue
      subroutine funct(r,f,df,ddf,n,npt,prnt)
      implicit integer (a-z)
      dimension r(npt), f(npt,n), df(npt,n), ddf(npt,n)
      real*8 r, f, df, ddf, alfa, fpkey, fac
      character*16 type, chrkey
      character*80 card, title
      logical prnt
      common /io/ inp, iout
      do 10 i=1,n
         read(inp,1) card
         type=chrkey(card,'type','power',' ')
         if (type.eq.'power') then
             np=intkey(card,'power',1,' ')
             do 20 j=1,npt
                f(j,i)=r(j)**np
 20          continue
             call rzero(df(1,i),npt)
             if(np.gt.0) then   
                do 30 j=1,npt
                   df(j,i)=np*r(j)**(np-1)
 30             continue
             endif 
             call rzero(ddf(1,i),npt)
             if (np.gt.1) then
                 do 40 j=1,npt               
                    ddf(j,i)=np*(np-1)*r(j)**(np-2)
 40              continue
             endif 
         elseif (type.eq.'exponential') then
             np=intkey(card,'power',1,' ')
             alfa=fpkey(card,'type=exponential',1.d0,' ')
             do 50 j=1,npt
                f(j,i)=r(j)**np*exp(-alfa*r(j))
 50          continue
             if (np.gt.0) then
                 do 60 j=1,npt
                    fac=r(j)**(np-1)
                    df(j,i)=( np*fac - alfa*r(j)*fac )
     1                       * exp(-alfa*r(j))
 60              continue
             else
                 do 70 j=1,npt
                    df(j,i)=-alfa*exp(-alfa*r(j))
 70              continue
             endif
             if (np.gt.1) then
                 do 80 j=1,npt
                    fac=r(j)**(np-2)
                    ddf(j,i)=( np*(np-1)*fac -2.d0*alfa*np*fac*r(j) +
     1                         alfa*alfa*fac*r(j)*r(j) )
     2                        * exp(-alfa*r(j))
 80              continue
             elseif (np.eq.1) then
                 do 90 j=1,npt
                    ddf(j,i)= ( -2.d0*alfa +alfa*alfa*r(j) )*
     1                               exp(-alfa*r(j))
 90              continue
             elseif( np.eq.0) then
                 do 100 j=1,npt
                    ddf(j,i)=alfa*alfa*exp(-alfa*r(j))
 100             continue
             endif
         else
             call lnkerr('error in function type')
         endif
 10   continue
      if (prnt) then
          title='basis functions at grid points'
          call prntrm(title,f,npt,n,npt,n,iout)
          title='1 derivative of basis functions at grid points'
          call prntrm(title,df,npt,n,npt,n,iout)
          title='2 derivative of basis functions at grid points'
          call prntrm(title,ddf,npt,n,npt,n,iout)
      endif
      return
 1    format(a80)
      end
