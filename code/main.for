       PROGRAM Main
        IMPLICIT NONE
        real*8 Ms,R,alpha
        write(6,*) 'Masas solares, R [km]= ,alpha ='
        read(5,*) Ms, R, alpha
        CALL COMPUTE(Ms,R,alpha)
       END PROGRAM Main

       SUBROUTINE COMPUTE(Ms,R,alpha)
cf2py  intent(in) Ms,R,alpha       
       IMPLICIT NONE
       integer nvar,nok,nbad,NMAX, KMAXX, kmax, kount, MAXSTP,k
       integer jx, jth, Nx, Nth
       PARAMETER (nvar=1,Nx=10,Nth=10)
       real*8 Ms, R, Masas_solares,TINY, hx, hth, pi,alpha
       INCLUDE 'size.for'
       PARAMETER (TINY=1.d-30,pi=4.d0*datan(1.d0))
       real*8 G,c,Msun,M,Umax,xc,xang(0:Nx),phistart(0:NMAX)    
       real*8 dusav, usol(0:KMAXX), phisol(0:NMAX,0:KMAXX)
       real*8 uaout(0:KMAXX), phiaout(NMAX,0:KMAXX)
       real*8 uax(0:KMAXX,Nx),rax(0:KMAXX,0:Nx)
       real*8 theta0(0:Nth),phiax(0:NMAX,0:KMAXX,0:Nx)
       real*8 xt(0:KMAXX,0:Nx,0:Nth),yt(0:KMAXX,0:Nx,0:Nth),
     X        zt(0:KMAXX,0:Nx,0:Nth) 
       real*8 xtr(0:KMAXX,0:Nx,0:Nth),ytr(0:KMAXX,0:Nx,0:Nth),
     X        ztr(0:KMAXX,0:Nx,0:Nth) 
       integer kounta(0:Nx,0:Nth),conteo(0:Nx,0:Nth)
       real*8 u1,u2,eps,h1,hmin,deltaphi
       COMMON /puntos/ xt,yt,zt,kounta
       COMMON /rays/ xtr,ytr,ztr,conteo
        print*, Ms,R
        G = 6.6743015d-20
        c = 299792.458d0
        Msun = 1.98847d30
        Masas_solares = G*Msun/c**2
        M = Ms*Masas_solares
        Umax = M/R
        write(6,*) 'Umax = ', Umax
        dusav=1.d-10
        kmax=KMAXX
        alpha=alpha*pi/180.d0
        xc=3.d0*Umax*dsqrt(3.d0*(1.d0-2.d0*Umax))
        if (Umax<1.d0/3.d0)  xc=1.d0
        hx = 0.999999d0*xc/dfloat(Nx)
        do jx=1,Nx
          xang(jx)=dfloat(jx)*hx
          call Trayectoria(xang(jx),Umax,kount,usol,phisol)
          do k=1,kount
            uax(k,jx) = usol(k)
            phiax(nvar,k,jx)=phisol(nvar,k)
            rax(k,jx) = M/uax(k,jx)
            xt(k,jx,1) = rax(k,jx)*dcos(phiax(nvar,k,jx))
            yt(k,jx,1) = 0.d0
            zt(k,jx,1) = rax(k,jx)*dsin(phiax(nvar,k,jx))
            if (k==kount) kounta(jx,1) = kount
            xtr(k,jx,1) = xt(k,jx,1)*dsin(alpha)-zt(k,jx,1)*dcos(alpha)
            ytr(k,jx,1) = yt(k,jx,1)
            ztr(k,jx,1) = -xt(k,jx,1)*dcos(alpha)+
     X                     zt(k,jx,1)*dsin(alpha)
            do jth=2,Nth
              theta0(jth) = 2.d0*pi/dfloat(Nth)*dfloat(jth)
              xt(k,jx,jth) = xt(k,jx,1)
              yt(k,jx,jth) = zt(k,jx,1)*dsin(theta0(jth))
              zt(k,jx,jth) = zt(k,jx,1)*dcos(theta0(jth))
              if (k==kount) kounta(jx,jth) = kount
              xtr(k,jx,jth) = xt(k,jx,jth)*dsin(alpha)-
     X                        zt(k,jx,jth)*dcos(alpha)
              ytr(k,jx,jth) = yt(k,jx,jth)
              ztr(k,jx,jth) = -xt(k,jx,jth)*dcos(alpha)+
     X                         zt(k,jx,jth)*dsin(alpha)
            enddo
          enddo
        enddo
c        write(6,*) 'xc = ', xc
c        write(6,*) 'Da un x'
c        read(5,*) x
c        do k=kount-3,kount
c        write(6,*) k, usol(k), phisol(nvar,k)
c        end do
c        print*, 'Umax, phimax = ', Umax, phisol(nvar,kount)
c        go to 22
        RETURN
        END SUBROUTINE COMPUTE

       SUBROUTINE Trayectoria(xin,Umaxin,kountout,uaout,phiaout)
        IMPLICIT NONE
        integer nvar,nok,nbad,NMAX,KMAXX,kmax,kount,kountout,n
        integer MAXSTP,k,noend
        real*8 Ms, R, Masas_solares,TINY,s,deltaphi
        EXTERNAL deltaphi, midsqu
        INCLUDE 'size.for'
        PARAMETER (TINY=1.d-30)
        real*8 G,c,Msun,M,Umaxin,Umax,x0,xin,phistart(0:NMAX)
        real*8 dxsav, ua(0:KMAXX), phia(0:NMAX,0:KMAXX)
        real*8 uaout(0:KMAXX), phiaout(0:NMAX,0:KMAXX)
        real*8 u1,u2,eps,h1,hmin
        COMMON /path/ kmax,kount,dxsav,ua,phia
        common /parameters/ x0, Umax
        x0=xin
        Umax=Umaxin
        nvar=1
        dxsav = 1.d-10
        phistart(nvar) = 0.d0
        u1 = 0.d0
        u2 = 0.999999d0*Umax     !0.99*Umax
        eps = 1.d-16 !minima = 1.d-16 para MAXSTP=1000000000
        h1=1.d-5
        hmin=0.d0
        kmax = KMAXX
        call odeint(phistart,nvar,u1,u2,eps,h1,hmin,nok,nbad,noend)
c        print*, 'nok = ',nok, 'nbad = ',nbad, 'noend = ', noend
c        if (noend==1) then
                call qromo(deltaphi,ua(kount),Umax,s,midsqu)
                kount=kount+1
                ua(kount) = Umax
                phia(nvar,kount) = phia(nvar,kount-1)+s
c        end if
c        if (noend==0) then
c                call qromo(deltaphi,ua(kount),Umax,s,midsqu)
c                print*, 'prueba : ', phia(nvar,kount), 
c     x                               phia(nvar,kount)+s 
c        end if
        do k=1,kount
                uaout(k) = ua(k)
                phiaout(nvar,k) = phia(nvar,k)
c                write(6,*) 'subrutina trayectoria ', k, uaout(k), 
c     x          phiaout(nvar,k)
        enddo
        kountout=kount
c        write(6,*) 'noend = ', noend
c        write(6,*) '+', kountout
c        write(6,*) 'phimax = ', phiaout(nvar,kountout)
        RETURN
        END SUBROUTINE Trayectoria
       
       SUBROUTINE derivs(u,phi,dphidu)
        integer NMAX
        INCLUDE 'size.for'
        real*8 Umax,x0,u,phi(0:NMAX),dphidu(0:NMAX)
        common /parameters/ x0, Umax
      
        dphidu(1) = x0/dsqrt((1.d0-2.d0*Umax)*Umax**2-
     X                      (1.d0-2.d0*u   )*u**2*x0**2)
        RETURN
        END SUBROUTINE derivs

        FUNCTION deltaphi(u)
        real*8 Umax,x0,u,deltaphi
        INCLUDE 'size.for'
        common /parameters/ x0, Umax
      
        deltaphi = x0/dsqrt((1.d0-2.d0*Umax)*Umax**2-
     X                      (1.d0-2.d0*u   )*u**2*x0**2)
        RETURN
        END FUNCTION deltaphi

c       INCLUDE 'odeint.for'
c       INCLUDE 'rkck.for'
c       INCLUDE 'rkqs.for'
c       INCLUDE 'ultimo.for'
c       INCLUDE 'polint.for'
c       INCLUDE 'midsqu.for'
