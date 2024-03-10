       SUBROUTINE midsqu(funk,aa,bb,s,n)
        IMPLICIT NONE
        integer n,NMAX
        PARAMETER (NMAX=50)
        real*8 aa,bb,s,funk
        EXTERNAL funk
        integer it,j
        real*8 ddel,del,ssum,tnm,x,func,a,b
        !This routine is an exact replacement for midpnt , except that it allows for an inverse square-root singularity in the integrand at the upper limit aa .
        func(x)=2.d0*x*funk(bb-x**2)
        b=dsqrt(bb-aa)
        a=0.d0
        if (n==1) then
            s=(b-a)*func(0.5d0*(a+b))
        else
            it=3**(n-2)
            tnm=it
            del=(b-a)/(3.d0*tnm)
            ddel=del+del
            x=a+0.5d0*del
            ssum=0.d0
            do j=1,it
                ssum=ssum+func(x)
                x=x+ddel
                ssum=ssum+func(x)
                x=x+del
            enddo
            s=(s+(b-a)*ssum/tnm)/3.d0
        endif
        return
        END SUBROUTINE midsqu
