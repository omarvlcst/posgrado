      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,noend,
     x                  kmax,kount,dxsav,xp,yp)
       IMPLICIT NONE
       INCLUDE 'size.for' 
       INTEGER nbad,nok,nvar,noend
       REAL*8 eps,h1,hmin,x1,x2,ystart(0:NMAX),TINY    
       PARAMETER (TINY=1.d-30)  
       INTEGER i,kmax,kount,nstp 
       REAL*8 dxsav,h,hdid,hnext,x,xsav,dydx(0:NMAX),y(0:NMAX),  
     X yscal(0:NMAX),ysol(0:MAXSTP,0:NMAX), Umax, x0
       REAL*8 xp(0:KMAXX),yp(0:NMAX,0:KMAXX)
ccc       COMMON /path/ kmax,kount,dxsav,xp,yp
ccc!$OMP   THREADPRIVATE (/path/)
       COMMON /parameters/ x0, Umax
!$OMP   THREADPRIVATE (/parameters/)
c HERE DANY
c        print *,'Odeint0: u1,u2,eps,h1,dxsav=',x1,x2,eps,h1,dxsav 
c        print *,'Odeint0: x0, Umax=',x0,Umax
c        print *,'Odeint0: ystart =',ystart
        x=x1
        h=sign(h1,x2-x1)
        nok=0
        nbad=0
        kount=0
        noend=0
c         print *,'odeint: nvar=',nvar
        do i=1,nvar
         y(i)=ystart(i)
c         print *,'odeint: i,y(i),ystart(i) =',i,y(i),ystart(i)
        enddo 
        if (kmax.gt.0) xsav=x-2.d0*dxsav  
        do nstp=1,MAXSTP  
         call derivs(x,y,dydx)
         do i=1,nvar  
           yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY  
	 enddo   
         if(kmax.gt.0)then  
           if(abs(x-xsav).gt.abs(dxsav)) then  
             if(kount.lt.kmax-1)then  
               kount=kount+1  
               xp(kount)=x  
               do i=1,nvar  
                 yp(i,kount)=y(i)  
	       enddo
c HERE DANY
c               print *,'Odeint1: k,x,y,dydx=',
c     x                  kount,xp(kount),yp(1,kount),dydx(1)
               xsav=x 
             endif  
           endif  
         endif  
         if((x+h-x2)*(x+h-x1).gt.0.d0) h=x2-x  
         call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext) 
c         write(6,*) nstp, Umax, x, y(1) 
         if(hdid.eq.h)then  
           nok=nok+1  
         else  
           nbad=nbad+1  
         endif  
         if((x-x2)*(x2-x1).ge.0.d0)then  
           do i=1,nvar  
             ystart(i)=y(i)
	   enddo  
           if(kmax.ne.0)then  
             kount=kount+1  
             xp(kount)=x 
             do i=1,nvar  
               yp(i,kount)=y(i)
	     enddo
           endif  
c HERE DANY
c            print *,'Odeint2: k,x,y=',kount,xp(kount),yp(1,kount)
           return  
         endif  
         if(abs(hnext).lt.hmin) pause  
     X	     'stepsize smaller than minimum in odeint'  
         h=hnext  
        enddo  
c       pause 'too many steps in odeint'  
        noend=1
       return  
       END
