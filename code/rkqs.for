      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext) 
        INTEGER n,NMAX  
        REAL*8 eps,hdid,hnext,htry,x,dydx(0:n),y(0:n),yscal(0:n)   
        PARAMETER (NMAX=50)  
C	      USES derivs,rkck  
        INTEGER i  
        REAL*8 errmax,h,htemp,xnew,yerr(0:NMAX),ytemp(0:NMAX),SAFETY,
     X  PGROW,PSHRNK,ERRCON  
        PARAMETER (SAFETY=0.9d0,PGROW=-.2d0,PSHRNK=-.25d0,
     X  ERRCON=1.89d-4)  
	 h=htry  
1        call rkck(y,dydx,n,x,h,ytemp,yerr)  
         errmax=0.d0  
         do i=1,n  
	   errmax=max(errmax,abs(yerr(i)/yscal(i)))  
	 enddo
         errmax=errmax/eps  
         if(errmax.gt.1.d0)then  
          htemp=SAFETY*h*(errmax**PSHRNK)  
          h=sign(max(abs(htemp),0.1d0*abs(h)),h)  
          xnew=x+h  
          if(xnew.eq.x)pause 'stepsize underflow in rkqs'  
	   goto 1  
	  else  
	   if(errmax.gt.ERRCON)then  
	    hnext=SAFETY*h*(errmax**PGROW)  
	   else  
	    hnext=5.d0*h  
	   endif  
	   hdid=h  
	   x=x+h  
	   do i=1,n  
	    y(i)=ytemp(i)  
           enddo  
          return  
         endif  
	END