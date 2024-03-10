        PROGRAM Prueba
         IMPLICIT NONE
         INTEGER i,n,x(0:100)
         x(0)=0
         n=10
c         do n=10,100
            do i=1,n,2
                x(i)=x(i-1)+2
                print*, x(i)
            enddo
c            print*, x(n)
c         enddo
         RETURN
        END PROGRAM Prueba