
! $Header: /home/cvsroot/LMDZ4/libf/filtrez/jacobi.F,v 1.1.1.1 2004/05/19
! 12:53:09 lmdzadmin Exp $

SUBROUTINE jacobi(a, n, np, d, v, nrot)
  PARAMETER (nmax=400)
  DIMENSION a(np, np), d(np), v(np, np), b(nmax), z(nmax)

  IF (n>nmax) THEN
    PRINT *, 'n, nmax=', n, nmax
    PRINT *, 'Surdimensionnement insuffisant dans jacobi'
    STOP 1
  END IF
  DO ip = 1, n
    DO iq = 1, n
      v(ip, iq) = 0.
    END DO
    v(ip, ip) = 1.
  END DO
  DO ip = 1, n
    b(ip) = a(ip, ip)
    d(ip) = b(ip)
    z(ip) = 0.
  END DO
  nrot = 0
  DO i = 1, 50
    sm = 0.
    DO ip = 1, n - 1
      DO iq = ip + 1, n
        sm = sm + abs(a(ip,iq))
      END DO
    END DO
    IF (sm==0.) RETURN
    IF (i<4) THEN
      tresh = 0.2*sm/n**2
    ELSE
      tresh = 0.
    END IF
    DO ip = 1, n - 1
      DO iq = ip + 1, n
        g = 100.*abs(a(ip,iq))
        IF ((i>4) .AND. (abs(d(ip))+g==abs(d(ip))) .AND. (abs(d(iq))+g==abs( &
            d(iq)))) THEN
          a(ip, iq) = 0.
        ELSE IF (abs(a(ip,iq))>tresh) THEN
          h = d(iq) - d(ip)
          IF (abs(h)+g==abs(h)) THEN
            t = a(ip, iq)/h
          ELSE
            theta = 0.5*h/a(ip, iq)
            t = 1./(abs(theta)+sqrt(1.+theta**2))
            IF (theta<0.) t = -t
          END IF
          c = 1./sqrt(1+t**2)
          s = t*c
          tau = s/(1.+c)
          h = t*a(ip, iq)
          z(ip) = z(ip) - h
          z(iq) = z(iq) + h
          d(ip) = d(ip) - h
          d(iq) = d(iq) + h
          a(ip, iq) = 0.
          DO j = 1, ip - 1
            g = a(j, ip)
            h = a(j, iq)
            a(j, ip) = g - s*(h+g*tau)
            a(j, iq) = h + s*(g-h*tau)
          END DO
          DO j = ip + 1, iq - 1
            g = a(ip, j)
            h = a(j, iq)
            a(ip, j) = g - s*(h+g*tau)
            a(j, iq) = h + s*(g-h*tau)
          END DO
          DO j = iq + 1, n
            g = a(ip, j)
            h = a(iq, j)
            a(ip, j) = g - s*(h+g*tau)
            a(iq, j) = h + s*(g-h*tau)
          END DO
          DO j = 1, n
            g = v(j, ip)
            h = v(j, iq)
            v(j, ip) = g - s*(h+g*tau)
            v(j, iq) = h + s*(g-h*tau)
          END DO
          nrot = nrot + 1
        END IF
      END DO
    END DO
    DO ip = 1, n
      b(ip) = b(ip) + z(ip)
      d(ip) = b(ip)
      z(ip) = 0.
    END DO
  END DO
  STOP '50 iterations should never happen'
  RETURN
END SUBROUTINE jacobi
