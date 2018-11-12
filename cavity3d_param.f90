!***********************************************************************
!     Constants for simulation setup
!***********************************************************************
subroutine condition
    implicit none
    real(lbkind) :: pi

    nstep = 5000

    igrid = 129
    jgrid = 129
	  kgrid = 129
    ngrid = jgrid

  !
  !  Reynolds number
  !
    re = 200.0

  !
  !  uc : free-stream velocity normalized by particle velocity
  !
    uc = 0.1

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

!*********************************************************************
  subroutine boundary_a(n)
!*********************************************************************
    implicit none
    integer :: i,j,k,l,m,id,jd,kd
    real(4) :: resi
    public :: bounce_back

    do l=1,pdm

      id = c(1,l)
      jd = c(2,l)
      kd = c(3,l)
      m = bounce_back(l)

      i = 1
      do k = 1, kgrid
      do j = 1, jgrid
        if( (i-id)<1 .or. (i-id)>igrid .or. &
            (j-jd)<1 .or. (j-jd)>jgrid .or. &
            (k-kd)<1 .or. (k-kd)>kgrid )then
          f(i,j,k,l) = f(i,j,k,m)
        endif
      enddo; enddo
      i = igrid
      do k = 1, kgrid
      do j = 1, jgrid
        if( (i-id)<1 .or. (i-id)>igrid .or. &
            (j-jd)<1 .or. (j-jd)>jgrid .or. &
            (k-kd)<1 .or. (k-kd)>kgrid )then
          f(i,j,k,l) = f(i,j,k,m)
        endif
      enddo; enddo

      j = 1
      do k = 1, kgrid
      do i = 1, igrid
        if( (i-id)<1 .or. (i-id)>igrid .or. &
            (j-jd)<1 .or. (j-jd)>jgrid .or. &
            (k-kd)<1 .or. (k-kd)>kgrid )then
          f(i,j,k,l) = f(i,j,k,m)
        endif
      enddo; enddo

      k = 1
      do j = 1, jgrid
      do i = 1, igrid
        if( (i-id)<1 .or. (i-id)>igrid .or. &
            (j-jd)<1 .or. (j-jd)>jgrid .or. &
            (k-kd)<1 .or. (k-kd)>kgrid )then
          f(i,j,k,l) = f(i,j,k,m)
        endif
      enddo; enddo
      k = kgrid
      do j = 1, jgrid
      do i = 1, igrid
        if( (i-id)<1 .or. (i-id)>igrid .or. &
            (j-jd)<1 .or. (j-jd)>jgrid .or. &
            (k-kd)<1 .or. (k-kd)>kgrid )then
          f(i,j,k,l) = f(i,j,k,m)
        endif
      enddo; enddo

    enddo

    i = 1
    do k = 1, kgrid
    do j = 1, jgrid
       q(i,j,k,1:5) = macro_(f(i,j,k,0:pdm))
       q(i,j,k,2) = 0.0d0
       q(i,j,k,3) = 0.0d0
       q(i,j,k,4) = 0.0d0
    enddo; enddo
    i = igrid
    do k = 1, kgrid
    do j = 1, jgrid
       q(i,j,k,1:5) = macro_(f(i,j,k,0:pdm))
       q(i,j,k,2) = 0.0d0
       q(i,j,k,3) = 0.0d0
       q(i,j,k,4) = 0.0d0
    enddo; enddo

    j = 1
    do k = 1, kgrid
    do i = 1, igrid
       q(i,j,k,1:5) = macro_(f(i,j,k,0:pdm))
       q(i,j,k,2) = 0.0d0
       q(i,j,k,3) = 0.0d0
       q(i,j,k,4) = 0.0d0
    enddo; enddo

    k = 1
    do j = 1, jgrid
    do i = 1, igrid
       q(i,j,k,1:5) = macro_(f(i,j,k,0:pdm))
       q(i,j,k,2) = 0.0d0
       q(i,j,k,3) = 0.0d0
       q(i,j,k,4) = 0.0d0
    enddo; enddo
    k = kgrid
    do j = 1, jgrid
    do i = 1, igrid
       q(i,j,k,1:5) = macro_(f(i,j,k,0:pdm))
       q(i,j,k,2) = 0.0d0
       q(i,j,k,3) = 0.0d0
       q(i,j,k,4) = 0.0d0
    enddo; enddo

    j = jgrid
    do k = 1, kgrid
    do i = 1, igrid
!       q(i,j,k,1) = 1.0d0
       q(i,j,k,1) = q(i,j-1,k,1)
       q(i,j,k,2) = uc
       q(i,j,k,3) = 0.0d0
       q(i,j,k,4) = 0.0d0
       f(i,j,k,0:pdm) = equi_(q(i,j,k,1:5))
    enddo; enddo

  end subroutine boundary_a
!***********************************************************
  module bb

  use precision_module
  implicit none

  integer,parameter :: pdm = 14
  public :: macro_
  public :: equi_
  public :: bounce_back

  contains

  function bounce_back(l) result(m)
    implicit none
    integer :: l,m
    !
    select case(l)
    case( 1); m =  3
    case( 2); m =  4
    case( 3); m =  1
    case( 4); m =  2
    case( 5); m =  6
    case( 6); m =  5
    case( 7); m = 13
    case( 8); m = 14
    case( 9); m = 11
    case(10); m = 12
    case(11); m =  9
    case(12); m = 10
    case(13); m =  7
    case(14); m =  8
    end select
    !
  end function

  !*******************************************************************
  function macro_(f) result(q)
  !*******************************************************************
    integer      :: l
    real(lbkind) :: f(0:pdm)
    real(lbkind) :: q(5)
      q = 0.0d0
      do l = 0, pdm
        q(1) = q(1) + f(l)
        q(2) = q(2) + f(l)*c(1,l)
        q(3) = q(3) + f(l)*c(2,l)
        q(4) = q(4) + f(l)*c(3,l)
      enddo
      q(2) = q(2)/q(1)
      q(3) = q(3)/q(1)
      q(4) = q(4)/q(1)
  end function

  !*******************************************************************
  function equi_(q) result(feq)
  !*******************************************************************
    integer      :: l
    real(lbkind) :: uu, cu
    real(lbkind) :: q(5)
    real(lbkind) :: feq(0:pdm)
      uu  = q(2)*q(2) + q(3)*q(3) + q(4)*q(4)
      do l = 0, pdm
        cu = q(2)*c(1,l) + q(3)*c(2,l) + q(4)*c(3,l)
        feq(l) = wf(l)*q(1)*( 1.0d0 + 3.0d0*cu + 4.5d0*cu*cu - 1.5d0*uu )
      enddo
  end function
