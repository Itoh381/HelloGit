!***********************************************************************
!     Constants for simulation setup                
!***********************************************************************
subroutine condition
    implicit none
    real(lbkind) :: pi

    nstep = 50000

    igrid = 251
    jgrid = 101
    ngrid = jgrid

  !
  !  Reynolds number
  !
    re = 500.0

  !
  !  uc : free-stream velocity normalized by particle velocity
  !
    uc = 0.05

  !
  !  calculate relaxation parameter, tau
  !
    tau = 0.5 * (6.0 * uc * real(ngrid-1) / re + 1.0)

    write(6,'(a,f15.7)') '  Reynolds number  = ', re
    write(6,'(a,f15.7)') '  relaxation time  = ', tau

    pi = acos(-1.0d0)
    omg = 2.0*pi/5000.0
    xc = 50.0
    yc = 50.0

end subroutine

!***********************************************************************
!     Implement Bounce-back on upper/lower boundaries                          
!***********************************************************************
subroutine boundary(n)
    implicit none
    integer, parameter :: nz = 24
    integer, intent(in) :: n
    integer      :: i, j, k, l, m, id, jd
    real(lbkind) :: fa, uw, vw, pi, x, y, theta
    real(lbkind), save :: p(2,5*nz)
    real(lbkind) :: xt(2,2), yt(2,2)

    i = igrid
    do j = 1, jgrid
      q(i, j, 1) = 1.0
      q(i, j, 2) = q(i-1, j, 2)
      q(i, j, 3) = q(i-1, j, 3)
    enddo
    call equilibrium(i, i, 1, jgrid, f)

    i = 1
    do j = 1, jgrid
      q(i, j, 1) = q(i+1, j, 1)
      q(i, j, 2) = uc
      q(i, j, 3) = 0.0
    enddo
    call equilibrium(i, i, 1, jgrid, f)

    j = 1
    f(i, j, 1) = f(i, j, 3)
    f(i, j, 2) = f(i, j, 4)
    f(i, j, 5) = f(i, j, 7)
    fa = 0.5d0*( f(i, j, 6) + f(i, j, 8) )
    f(i, j, 6) = fa
    f(i, j, 8) = fa
    !
    j = jgrid
    f(i, j, 1) = f(i, j, 3)
    f(i, j, 4) = f(i, j, 2)
    f(i, j, 8) = f(i, j, 6)
    fa = 0.5d0*( f(i, j, 5) + f(i, j, 7) )
    f(i, j, 5) = fa
    f(i, j, 7) = fa

    j = 1
    do i = 2, igrid-1
      f(i, j, 2) = f(i, j, 4)
      f(i, j, 5) = f(i, j, 7)
      f(i, j, 6) = f(i, j, 8)
      fa = 0.5d0*( f(i, j, 1) + f(i, j, 3) )
      f(i, j, 1) = fa
      f(i, j, 3) = fa
    enddo

    j = jgrid
    do i = 2, igrid-1
      f(i, j, 4) = f(i, j, 2)
      f(i, j, 7) = f(i, j, 5)
      f(i, j, 8) = f(i, j, 6)
      fa = 0.5d0*( f(i, j, 1) + f(i, j, 3) )
      f(i, j, 1) = fa
      f(i, j, 3) = fa
    enddo

    do j = 1, jgrid
		do i = 1, igrid
			ibk(i,j) = 0
			s(i,j,:) = -1.0
			r(i,j,:) = -1.0
		enddo
	enddo

    pi = acos(-1.0d0)
    uw = 0.0
    vw = 0.05 * cos(2.0d0*pi*real(n)/1000.0d0)
    if( n == 1 )then
      !  definition of square cylinder
      p(1,1) = 20.3;  p(2,1) = 49.1!5.1
      p(1,2) = 20.3;  p(2,2) = 51.7!4.7
      p(1,3) = 30.9;  p(2,3) = 51.7!4.7
      p(1,4) = 30.9;  p(2,4) = 49.1!5.1
      p(1,5) = 20.3;  p(2,5) = 49.1!5.1
      do l = 1, 5
        x = p(1,l)-25.0
        y = p(2,l)-50.0
        theta = -pi/6.0
        p(1,l) = 25.0 + x*cos(theta) - y*sin(theta)
        p(2,l) = 50.0 + x*sin(theta) + y*cos(theta)
      enddo
      do l = 6, 5*nz
        x = p(1,l-5)-xc
        y = p(2,l-5)-yc
        theta = -2.0*pi/real(nz)
        p(1,l) = xc + x*cos(theta) - y*sin(theta)
        p(2,l) = yc + x*sin(theta) + y*cos(theta)
      enddo
      !
      open(11,file='fan.xyz',status='unknown')
      write(11,*) nz
      write(11,*) (2, 2,l=1,nz)
      do l=1,nz
        xt(1,1) = p(1,1+5*(l-1));  yt(1,1) = p(2,1+5*(l-1))
        xt(1,2) = p(1,2+5*(l-1));  yt(1,2) = p(2,2+5*(l-1))
        xt(2,2) = p(1,3+5*(l-1));  yt(2,2) = p(2,3+5*(l-1))
        xt(2,1) = p(1,4+5*(l-1));  yt(2,1) = p(2,4+5*(l-1))
        write(11,*) ((xt(i,j)/real(ngrid-1),i=1,2),j=1,2) &
                   ,((yt(i,j)/real(ngrid-1),i=1,2),j=1,2)
      enddo
      close(11)
      !
    else
      !  movement of square cylinder
      do l = 1, 5*nz
        uw =  (p(2,l)-yc)*omg
        vw = -(p(1,l)-xc)*omg
        p(1,l) = p(1,l) + uw
        p(2,l) = p(2,l) + vw
      enddo
    endif

    do m = 1, nz
		do l = 1, 4
			k = 5*(m-1) + l
			call crossx(p(1,k),p(1,k+1))
			call crossy(p(1,k),p(1,k+1))
		enddo
	enddo
    call cross

    do j = 1, jgrid
		do i = 1, igrid
			if( ibk(i,j) == 1 )then
				do l = 1, 8
					if( s(i,j,l) >= 0.0 .and. s(i,j,l) <= 1.0 )then
						call ibb2(i,j,l,s(i,j,l))
					endif
				enddo
			endif
		enddo
	enddo

end subroutine
 
subroutine crossx(a,b)!,u,v)
    implicit none
    integer      :: i, j
    integer      :: jst, jen
    real(lbkind) :: a(2), b(2)
    real(lbkind) :: t, xc, qq
!    real(lbkind) :: u, v
    !
    if( a(2) < b(2) )then
      jst = int(a(2)+1.0d0)
      jen = int(b(2))
      do j = jst, jen
        t = (real(j)-a(2))/(b(2)-a(2))
        xc = a(1) + t*(b(1)-a(1))
        i = int(xc)
        qq = xc - real(i)
        r(i,j,1) = qq
        ibk(i,j) = 1
      enddo
    else
      jst = int(a(2))
      jen = int(b(2)+1.0d0)
      do j = jst, jen, -1
        t = (real(j)-a(2))/(b(2)-a(2))
        xc = a(1) + t*(b(1)-a(1))
        i = int(xc)
        qq = xc - real(i)
        r(i,j,1) = qq
        ibk(i+1,j) = 1
      enddo
    endif
    !
end subroutine

subroutine crossy(a,b)!,u,v)
    implicit none
    integer      :: i, j
    integer      :: ist, ien
    real(lbkind) :: a(2), b(2)
    real(lbkind) :: t, yc, qq
!    real(lbkind) :: u, v
    !
    if( a(1) < b(1) )then
      ist = int(a(1)+1.0d0)
      ien = int(b(1))
      do i = ist, ien
        t = (real(i)-a(1))/(b(1)-a(1))
        yc = a(2) + t*(b(2)-a(2))
        j = int(yc)
        qq = yc - real(j)
        r(i,j,2) = qq
        ibk(i,j+1) = 1
      enddo
    else
      ist = int(a(1))
      ien = int(b(1)+1.0d0)
      do i = ist, ien, -1
        t = (real(i)-a(1))/(b(1)-a(1))
        yc = a(2) + t*(b(2)-a(2))
        j = int(yc)
        qq = yc - real(j)
        r(i,j,2) = qq
        ibk(i,j) = 1
      enddo
    endif
    !
end subroutine

subroutine cross
    implicit none
    integer :: i, j
    real(lbkind) :: a, b

    do j = 1, jgrid-1
    do i = 1, igrid-1
      if( r(i,j,1) >= 0.0 )then
        !  pattern 1
        if( r(i,j,2) >= 0.0 )then
          r(i,j,3) = r(i,j,1)*r(i,j,2)/(r(i,j,1)+r(i,j,2))
        endif
        !  pattern 2
        if( r(i,j+1,1) >= 0.0 )then
          r(i,j,  3) = r(i,j,  1)/(1.0+r(i,j,1)-r(i,j+1,1))
          r(i,j+1,4) = r(i,j+1,1)/(1.0-r(i,j,1)+r(i,j+1,1))
        endif
        !  pattern 3
        if( r(i+1,j,2) >= 0.0 )then
          r(i,j+1,4) = 1.0 - (1.0-r(i,j,1))*r(i,j,2)/(1.0-r(i,j,1)+r(i,j,2))
        endif
      endif
      !
      if( r(i,j,2) >= 0.0 )then
        !  pattern 4
        if( r(i,j+1,1) >= 0.0 )then
          r(i,j+1,4) = (1.0-r(i,j,2))*r(i,j+1,1)/(1.0-r(i,j,2)+r(i,j+1,1))
        endif
        !  pattern 5
        if( r(i+1,j,2) >= 0.0 )then
          r(i,j,  3) =      r(i,j,2) /(1.0+r(i,j,2)-r(i+1,j,2))
          r(i,j+1,4) = (1.0-r(i,j,2))/(1.0-r(i,j,2)+r(i+1,j,2))
        endif
      endif
      !
      if( r(i,j+1,1) >= 0.0 .and. r(i+1,j,2) >= 0.0 )then
        !  pattern 6
        a = 1.0 - r(i,j+1,1)
        b = 1.0 - r(i+1,j,2)
        r(i,j,3) = 1.0 - a*b/(a+b)
      endif
    enddo; enddo

    do j = 2, jgrid-1
    do i = 2, igrid-1
      s(i,j,1) = 1.0 - r(i-1,j,1)
      s(i,j,2) = 1.0 - r(i,j-1,2)
      s(i,j,3) = r(i,j,1)
      s(i,j,4) = r(i,j,2)
      s(i,j,5) = 1.0 - r(i-1,j-1,3)
      s(i,j,6) = r(i,j,4)
      s(i,j,7) = r(i,j,3)
      s(i,j,8) = 1.0 - r(i-1,j+1,4)
    enddo; enddo

end subroutine
 
subroutine ibb2(i,j,l,t)
    implicit none
    integer     , intent(in) :: i, j, l
    real(lbkind) :: u, v, x, y, pi
    integer      :: id, jd
    real(lbkind) :: t, a, fa
    !
    id = c(1,l)
    jd = c(2,l)
    !
    pi = acos(-1.0d0)
    x = i - t*c(1,l)
    y = j - t*c(2,l)
    u = (y-yc)*omg
    v =-(x-xc)*omg
    !
    fa = 6.0*w(l)*( c(1,l)*u + c(2,l)*v )
    if( t < 0.5 )then
      a = 2.0*t
      f(i, j, l) = a*f(i-id, j-jd, lb(l)) + (1.0-a)*f(i, j, lb(l)) + fa
    else
      a = 0.5/t
      f(i, j, l) = a*f(i-id, j-jd, lb(l)) + (1.0-a)*f(i+id, j+jd, l) + a*fa
    endif
    !
end subroutine