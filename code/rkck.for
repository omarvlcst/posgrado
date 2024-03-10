      SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr)  
       INTEGER n,NMAX  
       REAL*8 h,x,dydx(0:n),y(0:n),yerr(0:n),yout(0:n)  
       PARAMETER (NMAX=50)  
C      USES derivs  
       INTEGER i  
       REAL*8 ak2(0:NMAX),ak3(0:NMAX),ak4(0:NMAX),ak5(0:NMAX),  
     X ak6(0:NMAX),ytemp(0:NMAX),
     X A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53, 
     X B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6  
       PARAMETER (A2=.2d0,A3=.3d0,A4=.6d0,A5=1.d0,A6=.875d0,B21=.2d0,  
     X B31=3.d0/40.d0,B32=9.d0/40.d0,B41=.3d0,B42=-.9d0,B43=1.2d0, 
     X B51=-11.d0/54.d0,B52=2.5d0, B53=-70.d0/27.d0,B54=35.d0/27.d0,
     X B61=1631.d0/55296.d0,B62=175.d0/512.d0,  B63=575.d0/13824.d0,  
     X B64=44275.d0/110592.d0,B65=253.d0/4096.d0,C1=37.d0/378.d0,  
     X C3=250.d0/621.d0,C4=125.d0/594.d0,C6=512.d0/1771.d0,
     X DC1=C1-2825.d0/27648.d0,DC3=C3-18575.d0/48384.d0,
     X DC4=C4-13525.d0/55296.d0,DC5=-277.d0/14336.d0,  DC6=C6-.25d0)  
        do i=1,n  
         ytemp(i)=y(i)+B21*h*dydx(i)  
        enddo
        call derivs(x+A2*h,ytemp,ak2)  
        do i=1,n  
         ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))  
        enddo
        call derivs(x+A3*h,ytemp,ak3)  
        do i=1,n  
         ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))  
        enddo
        call derivs(x+A4*h,ytemp,ak4)  
        do i=1,n  
         ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))  
        enddo
        call derivs(x+A5*h,ytemp,ak5)  
        do i=1,n  
         ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+  
     X   B65*ak5(i))  
        enddo  
        call derivs(x+A6*h,ytemp,ak6)  
        do i=1,n  
         yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))  
        enddo  
        do i=1,n  
         yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*  
     X   ak6(i))  
        enddo 
       return  
      END  