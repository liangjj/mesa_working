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
      dimension r(0:npt-1), f(0:npt-1,n), df(0:npt-1,n), ddf(0:npt-1,n)
      real*8 r, f, df, ddf, alfa, fpkey, fac, center, aux
      character*16 type, chrkey
      character*80 card, title
      logical prnt
      common /io/ inp, iout
      do 10 i=1,n
         read(inp,1) card
         type=chrkey(card,'type','power',' ')
         if (type.eq.'power') then
             np=intkey(card,'power',1,' ')
             write(iout,2) i, type, np
             do 20 j=0,npt-1
                f(j,i)=r(j)**np
 20          continue
             call rzero(df(0,i),npt-1)
             if(np.gt.0) then   
                do 30 j=0,npt-1
                   df(j,i)=np*r(j)**(np-1)
 30             continue
             endif 
             call rzero(ddf(0,i),npt-1)
             if (np.gt.1) then
                 do 40 j=0,npt-1               
                    ddf(j,i)=np*(np-1)*r(j)**(np-2)
 40              continue
             endif 
         elseif (type.eq.'exponential') then
             np=intkey(card,'power',1,' ')
             alfa=fpkey(card,'exponent',1.d0,' ')
             write(iout,3) i, type, np, alfa
             do 50 j=0,npt-1
                f(j,i)=r(j)**np*exp(-alfa*r(j))
 50          continue
             if (np.gt.0) then
                 do 60 j=0,npt-1
                    fac=r(j)**(np-1)
                    df(j,i)=( np*fac - alfa*r(j)*fac )
     1                       * exp(-alfa*r(j))
 60              continue
             else
                 do 70 j=0,npt-1
                    df(j,i)=-alfa*exp(-alfa*r(j))
 70              continue
             endif
             if (np.gt.1) then
                 do 80 j=0,npt-1
                    fac=r(j)**(np-2)
                    ddf(j,i)=( np*(np-1)*fac -2.d0*alfa*np*fac*r(j) +
     1                         alfa*alfa*fac*r(j)*r(j) )
     2                        * exp(-alfa*r(j))
 80              continue
             elseif (np.eq.1) then
                 do 90 j=0,npt-1
                    ddf(j,i)= ( -2.d0*alfa +alfa*alfa*r(j) )*
     1                               exp(-alfa*r(j))
 90              continue
             elseif( np.eq.0) then
                 do 100 j=0,npt-1
                    ddf(j,i)=alfa*alfa*exp(-alfa*r(j))
 100             continue
             endif
         elseif (type.eq.'gaussian') then
             np=intkey(card,'power',1,' ')
             alfa=fpkey(card,'exponent',1.d0,' ')
             center=fpkey(card,'center',0.d0,' ')
             write(iout,4) i, type, np, alfa, center
             do 200 j=0,npt-1
                aux=(r(j)-center)
                fac=r(j)**np
                f(j,i)=fac*exp(-alfa*aux*aux)
 200         continue
             if (np.gt.0) then
                 do 210 j=0,npt-1
                    aux=(r(j)-center)
                    fac=r(j)**(np-1)
                    df(j,i)=( np*fac - 2.d0*alfa*aux*fac*r(j) )
     1                       * exp(-alfa*aux*aux)
 210             continue
             else
                 do 220 j=0,npt-1
                    aux=r(j)-center
                    df(j,i)=-2.d0*alfa*aux*exp(-alfa*aux*aux)
 220             continue
             endif
             if (np.gt.1) then
                 do 230 j=0,npt-1
                    aux=r(j)-center 
                    fac=r(j)**(np-2)
                    ddf(j,i)=( np*(np-1)*fac -4.d0*alfa*np*aux*fac*r(j)
     1                         -2.d0*alfa*fac*r(j)*r(j)
     2                         + 4.d0*alfa*alfa*fac*aux*aux*r(j)*r(j) )
     3                         * exp(-alfa*aux*aux)
 230             continue
             elseif (np.eq.1) then
                 do 240 j=0,npt-1
                    aux=r(j)-center 
                    ddf(j,i)= ( -4.d0*alfa*aux - 2.d0*alfa*r(j)
     1                          + 4.d0*alfa*alfa*r(j)*aux*aux )
     2                            *exp(-alfa*aux*aux)
 240             continue
             elseif( np.eq.0) then
                 do 250 j=0,npt-1
                    aux=r(j)-center 
                    ddf(j,i)=( -2.d0*alfa +4.d0*alfa*alfa*aux*aux )
     1                        *exp(-alfa*aux*aux)
 250             continue
             endif
         else
             call lnkerr('error in function type')
         endif
 10   continue
      if (prnt) then
          title='basis functions at grid points'
          call prntrm(title,f(0,1),npt,n,npt,n,iout)
          title='1 derivative of basis functions at grid points'
          call prntrm(title,df(0,1),npt,n,npt,n,iout)
          title='2 derivative of basis functions at grid points'
          call prntrm(title,ddf(0,1),npt,n,npt,n,iout)
      endif
      return
 1    format(a80)
 2    format('basis function =',i3,1x,'type =',a8,1x,'power = ',i3)
 3    format('basis function =',i3,1x,'type =',a16,1x,'power = ',i3,1x,
     1       'exponent =',e15.8)
 4    format('basis function =',i3,1x,'type =',a16,'power = ',i3,1x,
     1       'exponent =',e15.8,1x,'center = ',e15.8)
      end




