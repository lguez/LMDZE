SUBROUTINE xmist(imr, iml, p, dc)
  PARAMETER (r24=1./24.)
  DIMENSION p(-iml:imr+1+iml), dc(-iml:imr+1+iml)

  DO i = 1, imr
    tmp = r24*(8.*(p(i+1)-p(i-1))+p(i-2)-p(i+2))
    pmax = max(p(i-1), p(i), p(i+1)) - p(i)
    pmin = p(i) - min(p(i-1), p(i), p(i+1))
    dc(i) = sign(min(abs(tmp),pmax,pmin), tmp)
  END DO
  RETURN
END SUBROUTINE xmist

