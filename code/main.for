c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************
       PROGRAM Main
        IMPLICIT NONE
        integer N1,N2
        PARAMETER (N1=800,N2=800)
        INCLUDE 'size.for'
        !include 'errno.h' 
        real*8 R_max,Ms,R,alpha,beta

        real*8 Rmax, Rinf, xd(0:Nx), thetad(0:Nth),
     x         det(0:2*d1-1,0:2*d2-1,0:2)
cc        COMMON /detector/ Rmax, Rinf, xd, thetad, det
        COMMON /Detector/Rmax,xd,thetad,det
        real*8 phie(0:Nx,0:Nth),thetae(0:Nx,0:Nth)
        integer px(0:2)
        real*8 imagen(0:3000-1,0:6000-1,0:2)
        Ms = 1.5
        R = 12.
        alpha = 60.
        beta = 0.
        R_max=14.
        print *,'MAIN: calling COLOR'
        CALL COLOR(Rinf,Rmax,xd,thetad,thetae,phie,imagen,N1,N2,det)
        print *,'MAIN: COLOR done'
       END PROGRAM Main
c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************
      SUBROUTINE INITIALIZE(Msol,Radius,R_max)
       ! Msol is in solar mass and Radius in km
       IMPLICIT NONE
       INCLUDE 'size.for'
       INTEGER i1,i2,jx,jth
       REAL*8 G,c,Msun,pi
       PARAMETER(G=6.6743015d-20,c=299792.458d0,Msun=1.98847d30)
       PARAMETER (pi=4.d0*datan(1.d0))
       REAL*8 Msol,Radius,R_max,M,Umax,Rinf,xc
       REAL*8 Masas_solares
       COMMON /Star/M,Umax,Rinf,xc
       REAL*8 Rmax,xd(0:Nx),thetad(0:Nth),
     x        det(0:2*d1-1,0:2*d2-1,0:2)
       COMMON /Detector/Rmax,xd,thetad,det
       integer ith_max,ichunk
       COMMON /par/ith_max,ichunk
        print *,'INITIALIZE: Msol Radius =',Msol,Radius
        ! Get some basic quantities:
        Masas_solares = G*Msun/c**2
        M = Msol*Masas_solares
        Umax = M/Radius
        print *,'INITIALIZE: Umax=',Umax
        if (Umax<1.d0/3.d0) then
         xc=1.d0
         Rinf = Radius/dsqrt(1.d0-2.d0*Umax)
        else
         xc=3.d0*Umax*dsqrt(3.d0*(1.d0-2.d0*Umax))
         Rinf = dsqrt(3.d0)*Radius
        end if
        ! Initialize detector:
        Rmax=max(R_max,Rinf)
        do i1=0,2*d1-1
         do i2=0,2*d2-1
          det(i1,i2,0) = 0.d0
          det(i1,i2,1) = 0.d0
          det(i1,i2,2) = 0.d0
         end do
        end do
        xd(0)=0.d0
        do jx=1,Nx
         xd(jx)=dfloat(jx)*0.999999d0*xc/dfloat(Nx)
        end do
        do jth=1,Nth
         thetad(jth) = 2.d0*pi/dfloat(Nth)*dfloat(jth-1)
        end do
       RETURN
      END SUBROUTINE INITIALIZE
c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************
      SUBROUTINE COLOR(Rinf,Rmax,xd,thetad,thetae,phie,imagen,N1,N2,det)
cf2py   intent(in) Rinf,Rmax,xd,thetad,thetae,phie,imagen,N1,N2
cf2py   intent(out) det
        IMPLICIT NONE
        INCLUDE 'size.for'
        integer jx, jth, n, N1, N2
        EXTERNAL INCLINE,TRANSFORM
        real*8 beta,phi0,pi,x,y
        integer xx,yy
        PARAMETER (pi=4.d0*datan(1.d0))
        real*8 phie(0:Nx,0:Nth),thetae(0:Nx,0:Nth)
        integer px(0:2)
        real*8 imagen(0:3000-1,0:6000-1,0:2)
        real*8 Rmax, Rinf, xd(0:Nx), thetad(0:Nth),
     x         det(0:2*d1-1,0:2*d2-1,0:2)
        print *,'COLOR, ENTERING'
        beta=beta*pi/180.d0
        do jx=1,Nx
         do jth=1,Nth
          phi0 = phie(jx,jth)+beta
          n = int(phi0/(2.d0*pi))
          if (phi0.gt.2.d0*pi) phi0 = phi0-2.d0*n*pi
          px(1) = int(thetae(jx,jth)/pi*(N1)) 
          px(2) = int(phi0/2./pi*(N2))
          x = -xd(jx)*dcos(thetad(jth))/xd(Nx)
          y = +xd(jx)*dsin(thetad(jth))/xd(Nx)
          xx = int(float(d1)*(1.d0+x * Rinf/Rmax))
          yy = int(float(d2)*(1.d0+y * Rinf/Rmax))
          det(xx,yy,0) = imagen(px(1),px(2),0)
          det(xx,yy,1) = imagen(px(1),px(2),1)
          det(xx,yy,2) = imagen(px(1),px(2),2)
         end do
        end do
        print *,'COLOR: Done'
        RETURN
      END SUBROUTINE COLOR
c***********************************************************************
c***********************************************************************
       SUBROUTINE NMAGDIP(R,D,Tsup,Nph)
cf2py       intent(in) R,D,Tsup
cf2py       intent(out) Nph
       IMPLICIT NONE
        INCLUDE 'size.for'
        integer jx,jth
        REAL*8 xd(0:Nx),thetad(0:Nth),det(0:2*d1-1,0:2*d2-1,0:2)
        REAL*8 const,R,D,Tsup(0:Nx,0:Nth)
        REAL*8, PARAMETER :: kb=1.3807d-16
        REAL*8, PARAMETER :: hplanck=6.6261d-27
        REAL*8, PARAMETER :: zeta3=1.2020569d0
        REAL*8, PARAMETER :: c=2.99792458d10
        REAL*8, PARAMETER :: pi = 4.d0*datan(1.d0)
        REAL*8 dx,fx,dthetad,ft,sumax,sumathetad,Nph,Rmax
        COMMON /Detector/ Rmax,xd,thetad,det
        const=4.d0*kb**3*zeta3*R**2/(hplanck**3*c**2*D**2) !** 4.82d-20/ cm**2 * s * Kelvin**3
        sumax=0.d0
        do jx=0,Nx-1
          dx = xd(jx+1)-xd(jx)
          fx = 0.5d0*(xd(jx+1)+xd(jx))
          sumathetad=0.d0
          do jth=0,Nth-1
            dthetad=thetad(jth+1)-thetad(jth)
            ft = 0.25d0*(Tsup(jx,jth+1)+Tsup(jx,jth)+
     X                   Tsup(jx+1,jth)+Tsup(jx+1,jth+1))
            sumathetad=sumathetad+dthetad*ft**3
          end do
          sumax=sumax+dx*fx*sumathetad
        end do
        Nph=const*sumax
        print*, 'Flujo de fotones por segundo por cm2 = ', Nph
        RETURN
        END SUBROUTINE NMAGDIP
c***********************************************************************
c***********************************************************************
       SUBROUTINE NSPECTENERGY(R,D,Tsup,E,Nspec)
cf2py       intent(in) R,D,Tsup,E
cf2py       intent(out) Nspec
       IMPLICIT NONE
        INCLUDE 'size.for'
        integer jx,jth
        REAL*8 xd(0:Nx),thetad(0:Nth),det(0:2*d1-1,0:2*d2-1,0:2)
        REAL*8 const,R,D,Tsup(0:Nx,0:Nth),E
        REAL*8, PARAMETER :: kb=1.3807d-16
        REAL*8, PARAMETER :: kboltzev=8.6173d-5
        REAL*8, PARAMETER :: hplanck=6.6261d-27
        REAL*8, PARAMETER :: zeta3=1.2020569d0
        REAL*8, PARAMETER :: c=2.99792458d10
        REAL*8, PARAMETER :: pi = 4.d0*datan(1.d0)
        REAL*8 dx,fx,dthetad,ft,sumax,sumathetad,Nspec,Rmax
        COMMON /Detector/ Rmax,xd,thetad,det
        const=4.d0*kb**3*zeta3*R**2/(hplanck**3*c**2*D**2) !** 4.82d-20/ cm**2 * s * Kelvin**3
        sumax=0.d0
        do jx=0,Nx-1
          dx = xd(jx+1)-xd(jx)
          fx = 0.5d0*(xd(jx+1)+xd(jx))
          sumathetad=0.d0
          do jth=0,Nth-1
            dthetad=thetad(jth+1)-thetad(jth)
            ft = 0.25d0*(Tsup(jx,jth+1)+Tsup(jx,jth)+
     X                   Tsup(jx+1,jth)+Tsup(jx+1,jth+1))
            sumathetad=sumathetad+
     X                 dthetad*E**2/(dexp(E/kboltzev*ft)-1.d0)
          end do
          sumax=sumax+dx*fx*sumathetad
        end do
        Nspec=const*sumax
        print*, 'Flujo de fotones por segundo por cm2 por eV = ', Nspec
        RETURN
        END SUBROUTINE NSPECTENERGY
c***********************************************************************
c***********************************************************************
       SUBROUTINE TEMP(thetae,phie,gama,beta,Tpar,Tperp,Tsup)
cf2py  intent(in) thetae,phie,gama,beta,Tpar,Tperp
cf2py  intent(out) Tsup
        IMPLICIT NONE
        INCLUDE 'size.for'
        integer jx,jth
        REAL*8 gama,Tpar,Tperp,anggamma,chi04,pi,ang,tfac,phi0,beta
        REAL*8 Tsup(0:Nx,0:Nth),phie(0:Nx,0:Nth),thetae(0:Nx,0:Nth)
        PARAMETER (pi=4.d0*datan(1.d0))
        anggamma = gama*pi/180.d0
        beta=beta*pi/180.d0
        chi04 = (Tperp/Tpar)**4.d0
        Tsup(0,0) = Tpar
        do jx=1,Nx
          do jth=1,Nth
            phi0 = phie(jx,jth)+beta
            tfac = dsin(thetae(jx,jth))*dcos(phi0)*
     X             dsin(anggamma)+dcos(thetae(jx,jth))*
     X             dcos(anggamma)
            Tsup(jx,jth) = Tpar*(4.d0*tfac**2*(1.d0-chi04)/
     X                     (1.d0+3.d0*tfac**2)+chi04)**(0.25d0)
          end do
        end do
        RETURN
       END SUBROUTINE TEMP
c***********************************************************************
c***********************************************************************
       SUBROUTINE INCLINE(Radius,alpha,xt,yt,zt,thetae,phie)
cf2py  intent(in) Radius,alpha,xt,yt,zt
cf2py  intent(out) thetae,phie    
        IMPLICIT NONE
        INCLUDE 'size.for'
        integer l
        integer jx, jth, ult
        real*8 Radius,alpha,pi,a,phi,theta,sina,cosa
        PARAMETER (pi = 4.d0*datan(1.d0))
        real*8 xta(0:2),yta(0:2),zta(0:2) 
        real*8 xt(0:2,0:Nx,0:Nth),yt(0:2,0:Nx,0:Nth),
     X         zt(0:2,0:Nx,0:Nth) 
        real*8 phie(0:Nx,0:Nth),thetae(0:Nx,0:Nth)
c!!!        COMMON /angulos/ thetaeout,phieout
         print *,'INCLINE, ENTERING: Radius, alpha =',Radius,alpha
         a=alpha*pi/180.d0
         sina=dsin(a)
         cosa=dcos(a)
         do jx=1,Nx
           do jth=1,Nth
               do l=1,2
               xta(l) = xt(l,jx,jth)*sina
     X                        - zt(l,jx,jth)*cosa
               yta(l) = yt(l,jx,jth)
               zta(l) = xt(l,jx,jth)*cosa
     X                        + zt(l,jx,jth)*sina
             end do
             call TRANSFORM(Radius,xta(1),yta(1),zta(1),thetae(jx,jth),
     X                                                  phie(jx,jth))
           end do
         end do
         print *,'INCLINE: DONE !'
       RETURN
       END SUBROUTINE INCLINE
c**************************************************************************
c**************************************************************************
       SUBROUTINE TRANSFORM(Radius,xta,yta,zta,thetae,phie)
cf2py  intent(in) Radius,xta,yta,zta
cf2py  intent(out) thetae,phie 
       IMPLICIT NONE
       real*8 Radius,pi,phi,theta,tol
       real*8 xta,yta,zta
       real*8 phie,thetae
       parameter (pi=4.d0*datan(1.d0), tol=1.d-8)
             ! Definicion de theta:
             if (abs(zta).lt.Radius) then
              theta = dacos(zta/Radius)
             else if (zta.ge.Radius) then
              theta = 0.d0
             else if (zta.le.-Radius) then
              theta = pi
             end if
             ! Definicion de phi:
c!!             if (xta.gt.tol) then
c!!                 phi=datan(yta/xta)
c!!             else if (xta.lt.-tol) then
c!!                 if (yta.gt.tol) then
c!!                     phi=datan(yta/xta)+pi
c!!                 else if (yta.lt.tol) then
c!!                     phi=datan(yta/xta)-pi
c!!                 else if (abs(yta).lt.tol) then
c!!                     if (yta.gt.0.d0) then
c!!                         phi=pi
c!!                     else if (yta.lt.0.d0) then
c!!                         phi=-pi
c!!                     end if
c!!                 end if
c!!             else if (abs(xta).lt.tol) then
c!!                 if (yta.gt.tol) then
c!!                    phi=0.5d0*pi
c!!                 else if (yta.lt.tol) then
c!!                     phi=-0.5d0*pi
c!!                 else if (abs(yta).lt.tol) then
c!!                     phi=pi
c!!                 end if
c!!             end if
             if (theta.eq.0.d0) then
              phi=0.d0
             else if (theta.eq.pi) then
              phi=0.d0
             else
              if (abs(xta).lt.abs(Radius*dsin(theta))) then
               phi=dacos(xta/(Radius*dsin(theta)))
c!                 phi = datan2(xta,yta)
              else
               if (xta.ge.0.d0) then
                phi=0.d0
               else 
                phi=pi
               end if
              end if
              if (yta.le.0.d0) then
                 phi=2.d0*pi-phi
              end if
             end if
             thetae = theta
             phie = phi
         RETURN
       END SUBROUTINE TRANSFORM
c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************
      SUBROUTINE COMPUTE(M,Umaxin,xd,thetad,xt,yt,zt)
cf2py  intent(in) M,Umax,xd,thetad
cf2py  intent(out) xt,yt,zt
       IMPLICIT NONE
       integer nvar, kmax, kount, k
       integer jx, jth, l
       PARAMETER (nvar=1)
       real*8 TINY,pi,dusav
       INCLUDE 'size.for'
       PARAMETER (TINY=1.d-30,pi=4.d0*datan(1.d0))
       real*8 M,Umax
       real*8 usol(0:KMAXX), phisol(0:NMAX,0:KMAXX)
       real*8 uax,rax,phiax
       real*8 xt(0:2,0:Nx,0:Nth),yt(0:2,0:Nx,0:Nth),
     X        zt(0:2,0:Nx,0:Nth) 
       real*8 Rmax, xd(0:Nx), thetad(0:Nth),
     x         det(0:2*d1-1,0:2*d2-1,0:2)
       real*8 x0, Umaxin
       COMMON /parameters/ x0, Umax
       real*8 rate,time1,time2
       integer c1,c2,cr,cm
        ! First initialize the system_clock
        CALL system_clock(count_rate=cr)
        CALL system_clock(count_max=cm)
        CALL SYSTEM_CLOCK(c1)
        CALL CPU_TIME(time1)
        rate = REAL(cr)
c        WRITE(*,*) "system_clock rate ",rate
        Umax=Umaxin
        x0=0.d0   ! Will be overwritten in the loop
        print *,'COMPUTE, ENTERING: M, Umax =', M, Umax
        dusav=1.d-10
        kmax=KMAXX
        do jx=1,Nx
         x0=xd(jx)
         call Trayectoria(kount,usol,phisol)
         do k=kount-1,kount
          l = kount - k + 1
          uax = usol(k)
          phiax = phisol(nvar,k)
          rax = M/uax
          xt(l,jx,1) = rax*dcos(phiax)
          yt(l,jx,1) = 0.d0
          zt(l,jx,1) = rax*dsin(phiax)
          do jth=1,Nth
           xt(l,jx,jth) = xt(l,jx,1)
           yt(l,jx,jth) = zt(l,jx,1)*dsin(thetad(jth))
           zt(l,jx,jth) = zt(l,jx,1)*dcos(thetad(jth))
          enddo
         enddo
        enddo
        CALL SYSTEM_CLOCK(c2)
        CALL CPU_TIME(time2)
c        WRITE(*,*) "system_clock : ",(c2 - c1)/rate
c        WRITE(*,*) "cpu_time     : ",(time2-time1)
        print *,'COMPUTE: DONE !'
       RETURN
      END SUBROUTINE COMPUTE
c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************
      SUBROUTINE Trayectoria(kountout,uaout,phiaout)
cf2py intent(out) kountout,uaout,phiaout
       IMPLICIT NONE
       integer nvar,nok,nbad,kmax,kount,kountout
       integer k,noend
       real*8 TINY,suma,deltaphi
       EXTERNAL deltaphi, derivs, midsqu, qromo, polint
       INCLUDE 'size.for'
       PARAMETER (TINY=1.d-30)
       real*8 Umax,x0,phistart(0:NMAX)
       real*8 dxsav, ua(0:KMAXX), phia(0:NMAX,0:KMAXX)
       real*8 uaout(0:KMAXX), phiaout(0:NMAX,0:KMAXX)
       real*8 u1,u2,eps,h1,hmin
       common /parameters/ x0, Umax
       common /puntos/ kount,ua,phia
        nvar=1
        dxsav = 1.d-5
        phistart(nvar) = 0.d0
        u1 = 0.d0
        u2 = 0.99999d0*Umax     !0.99*Umax
        eps = 1.d-16 !minima = 1.d-16 para MAXSTP=1000000000
        h1=1.d-3
        hmin=0.d0
        suma=0.d0
        kmax = KMAXX
        call odeint(phistart,nvar,u1,u2,eps,h1,hmin,nok,nbad,noend,
     x              kmax,kount,dxsav,ua,phia)
c!!        call qromo(deltaphi,ua(kount),Umax,suma,midsqu)
        call midsqu(deltaphi,ua(kount),Umax,suma,5)
        PRINT*, 'kount, Umax,ua(kount), suma:  ',kount,Umax,
     x                                           ua(kount),suma
        kount=kount+1
        ua(kount) = Umax
        phia(nvar,kount) = phia(nvar,kount-1)+suma
        print*, phia(nvar,kount)
        do k=1,kount
         uaout(k) = ua(k)
         phiaout(nvar,k) = phia(nvar,k)
        enddo
        kountout=kount
       RETURN
      END SUBROUTINE Trayectoria
c***********************************************************************
c***********************************************************************
       SUBROUTINE derivs(u,phi,dphidu)
        IMPLICIT NONE
        INCLUDE 'size.for'
        real*8 Umax,x0,u,phi(0:NMAX),dphidu(0:NMAX)
        common /parameters/ x0, Umax
        dphidu(1) = x0/dsqrt((1.d0-2.d0*Umax)*Umax**2-
     X                       (1.d0-2.d0*u   )*u**2*x0**2)
        RETURN
        END SUBROUTINE derivs
c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************
        FUNCTION deltaphi(u)
        IMPLICIT NONE
        real*8 Umax,x0,u,deltaphi
        INCLUDE 'size.for'
        common /parameters/ x0, Umax
        deltaphi = x0/dsqrt((1.d0-2.d0*Umax)*Umax**2-
     X                      (1.d0-2.d0*u   )*u**2*x0**2)
c!!!	deltaphi=x0*(Umax-u)
        RETURN
        END FUNCTION deltaphi
c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************
      subroutine xterm
       implicit real*8(a-h,k-z)
       character*40 device
c***********************************************************************
c BETTER DO IT FROM THE COMMAND LINE SO THE XTERM CAN BE OPENED REMOTELY
c Open an Xterm for output:
c        call system(' xterm -T "Python NSCool"
c     x  -geometry 145x40 -sl 100000
c     x  -e "tty > Xterm_device ; bash --init-file Xterm_prompt " & ')
c        call sleep(1)
c***********************************************************************
c Read the Xterm tty file:
        open(unit=100,file='Xterm_device',status='old')
        read(100,'(40a1)',end=999)(device(i:i),i=1,40)
 999    close(unit=100,status='keep')
c        call system('rm -f Xterm_device')
c Just print out the tty file name:
c        print *,'Xterm tty file is: ',trim(device)
c***********************************************************************
c Assign stdout, unit 6, to the Xterm:
        open(unit=6,file=trim(device),status='old')
c Assign stderr, unit 0 in gfortran, to the Xterm:
        open(unit=0,file=trim(device),status='old')
c Warning: that stderr is unit 0 is compiler dependent !
c***********************************************************************
        write(6,*)'Trayectoria in Python mode: let''s go !'
c***********************************************************************
      end subroutine xterm
c***********************************************************************
c***********************************************************************

c       INCLUDE 'odeint.for'
c       INCLUDE 'rkck.for'
c       INCLUDE 'rkqs.for'
c       INCLUDE 'ultimo.for'
c       INCLUDE 'polint.for'
c       INCLUDE 'midsqu.for'
