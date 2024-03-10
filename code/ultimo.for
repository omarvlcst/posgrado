       SUBROUTINE qromo(func,a,b,suma,midsqu)
        IMPLICIT NONE
        INTEGER JMAX,JMAXP,K,KM,NMAX
        REAL*8 a,b,func,ss,EPS,suma
        EXTERNAL func,midsqu
        PARAMETER (EPS=1.d-6, JMAX=14, JMAXP=JMAX+1, K=5, KM=K-1)
c        USES polint
        INTEGER j
        REAL*8 dss,h(0:JMAXP),s(0:JMAXP)
        h(1)=1.d0
c!!!!!        print*, 'JMAX=15'
        do j=1,JMAX
            call midsqu(func,a,b,s(j),j)
            if (j.ge.K) then
                call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
c!!!!!                print*, 'ss antes de IF', ss
                if (abs(dss).le.EPS*abs(ss)) then
c!!!!!                    print*, 'ss dentro de IF', j,ss
                    print*, 'qromo IF j,s(j),ss, dss: ', j,s(j),ss,dss
                    suma=s(j)
                    return
                end if
            endif
            s(j+1)=s(j)
            h(j+1)=h(j)/9.d0
            print*, 'qromo OUT j,h,s :', j, h(j), s(j)
        enddo
c!!!!       pause 'too many steps in qromo'
       END SUBROUTINE qromo
