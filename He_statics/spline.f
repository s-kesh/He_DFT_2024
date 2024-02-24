c---------------------------------------------------------------------------------
c                       subroutine findi
c---------------------------------------------------------------------------------
      subroutine findi2(x,xc,n,i)
c
c      finds the interval (x(i),x(i+1)) containing the value xc.
c   input:
c     x(i) (i=1, ...,n) ........ grid points.
c                    (the x values must be in increasing order).
c     xc ....................... point to be located.
c     n ........................ number of grid points.
c   output:
c     i ........................ interval index.
c
      implicit double precision (a-h,o-z)
      dimension x(*)
      if(xc.gt.x(n)) then
         i=n-1
         return
      endif
      if(xc.lt.x(1)) then
         i=1
         return
      endif
      i=1
      i1=n
    1 it=(i+i1)/2
      if(xc.gt.x(it)) i=it
      if(xc.le.x(it)) i1=it
      if(i1-i.gt.1) go to 1
      return
      end
c---------------------------------------------------------------------------------
c                       subroutine spline
c---------------------------------------------------------------------------------
      subroutine spline2(x,y,a,b,c,d,s1,sn,n)
c
c      cubic spline interpolation between tabulated data.
c   input:
c     x(i) (i=1, ...,n) ........ grid points.
c                    (the x values must be in increasing order).
c     y(i) (i=1, ...,n) ........ corresponding function values.
c     s1,sn ..... second derivatives at x(1) and x(n).
c            (the natural spline corresponds to taking s1=sn=0).
c     n ........................ number of grid points.
c      the interpolating polynomial in the i-th interval, from
c   x(i) to x(i+1), is
c            pi(x) = a(i)+x*(b(i)+x*(c(i)+x*d(i)))
c   output:
c     a(i),b(i),c(i),d(i) ...... spline coefficients.
c
c      ref.: m.j. maron, 'numerical analysis: a practical
c            approach', macmillan publ. co., new york 1982.
c
      implicit double precision (a-h,o-z)
      dimension x(*),y(*),a(*),b(*),c(*),d(*)
      if(n.lt.4) then
         write(*,10) n
         stop '0001'
      endif
      n1=n-1
      n2=n-2
c  ****  auxiliary arrays h(=a) and delta(=d).
      do i=1,n1
         if(x(i+1)-x(i).lt.1.0d-10) then
            write(*,11)
            write(*,7777) (i1,x(i1),i1=1,n1)
            stop '0002'
         endif
         a(i)=x(i+1)-x(i)
         d(i)=(y(i+1)-y(i))/a(i)
      end do
c  ****  symmetric coefficient matrix (augmented).
      do i=1,n2
         b(i)=2.0d0*(a(i)+a(i+1))
         k=n1-i+1
         d(k)=6.0d0*(d(k)-d(k-1))
      end do
      d(2)=d(2)-a(1)*s1
      d(n1)=d(n1)-a(n1)*sn
c  ****  gauss solution of the tridiagonal system.
      do i=2,n2
         r=a(i)/b(i-1)
         b(i)=b(i)-r*a(i)
         d(i+1)=d(i+1)-r*d(i)
      end do
c  ****  the sigma coefficients are stored in array d.
      d(n1)=d(n1)/b(n2)
      do i=2,n2
         k=n1-i+1
         d(k)=(d(k)-a(k)*d(k+1))/b(k-1)
      end do
      d(n)=sn
c  ****  spline coefficients.
      si1=s1
      do i=1,n1
         si=si1
         si1=d(i+1)
         h=a(i)
         hi=1.0d0/h
         a(i)=(hi/6.0d0)*(si*x(i+1)**3-si1*x(i)**3)
     1        +hi*(y(i)*x(i+1)-y(i+1)*x(i))
     2        +(h/6.0d0)*(si1*x(i)-si*x(i+1))
         b(i)=(hi/2.0d0)*(si1*x(i)**2-si*x(i+1)**2)
     1        +hi*(y(i+1)-y(i))+(h/6.0d0)*(si-si1)
         c(i)=(hi/2.0d0)*(si*x(i+1)-si1*x(i))
         d(i)=(hi/6.0d0)*(si1-si)
      end do
c
      return
c
   10 format(5x,'Spline interpolation cannot be performed with',
     1i4,' points. stop.')
   11 format(5x,'Spline x values not in increasing order. stop.')
 7777 format((1h ,0p,i4,1p,d28.17))
      end
