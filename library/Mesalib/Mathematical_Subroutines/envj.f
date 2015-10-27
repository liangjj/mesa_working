*deck envj
        function envj(n,x)
        real*8 x, envj
        envj=0.5d0*log10(6.28d0*n)-n*log10(1.36d0*x/n)
        return
        end
