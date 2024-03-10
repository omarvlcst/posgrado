       PROGRAM MAIN
        IMPLICIT NONE
        REAL*8 a,b,funcion,s
        EXTERNAL midsqu,qromo
        a = 0.d0
        b = 1.d0
        s = 0.d0
        call qromo(funcion,a,b,s,midsqu)
        PRINT*, s 
        END PROGRAM MAIN
        
        FUNCTION funcion(x)
        IMPLICIT NONE
        real*8 x, funcion
        funcion=dsqrt(1.d0-x**2)
        END FUNCTION
