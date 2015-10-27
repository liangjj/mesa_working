      double precision function fun(x)
      implicit real*8(a-h,o-z)
      common/intg/alpa,beta,a,b,n
      arg=alpa*x*x+beta*(x-a)**2+beta*(x-b)**2
      fun=exp(-arg)*x**n
      return
      end
