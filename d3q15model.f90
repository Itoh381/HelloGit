!*********************************************************************
  module D3Q15model
!*********************************************************************
  use precision_module
  implicit none

  integer, parameter :: pdm = 14    ! number of particle -1
  real(lbkind)       :: c(3,0:pdm)  ! particle velocity
  real(lbkind)       :: wf(0:pdm)   ! weight coefficient

  public :: init_D3Q15
  public :: macro_
  public :: equi_
  public :: bounce_back

  contains

!*******************************************************************
  subroutine init_D3Q15
!*******************************************************************
    integer :: l

    c = 0.0d0

    !> define the particle velocity c

    !>    D3Q15model
    !//   8,12     2       7,11
    !//            |
    !//   3 -      0,5,6 - 1
    !//            |
    !//   9,13     4       10,14

    !> c(1 = x direction  2 = y direction, number of direction)

    c(1, 0) =  0.0;  c(2, 0) =  0.0;  c(3, 0) =  0.0
    c(1, 1) =  1.0;  c(2, 1) =  0.0;  c(3, 1) =  0.0
    c(1, 2) =  0.0;  c(2, 2) =  1.0;  c(3, 2) =  0.0
    c(1, 3) = -1.0;  c(2, 3) =  0.0;  c(3, 3) =  0.0
    c(1, 4) =  0.0;  c(2, 4) = -1.0;  c(3, 4) =  0.0
    c(1, 5) =  0.0;  c(2, 5) =  0.0;  c(3, 5) =  1.0
    c(1, 6) =  0.0;  c(2, 6) =  0.0;  c(3, 6) = -1.0
    c(1, 7) =  1.0;  c(2, 7) =  1.0;  c(3, 7) =  1.0
    c(1, 8) = -1.0;  c(2, 8) =  1.0;  c(3, 8) =  1.0
    c(1, 9) = -1.0;  c(2, 9) = -1.0;  c(3, 9) =  1.0
    c(1,10) =  1.0;  c(2,10) = -1.0;  c(3,10) =  1.0
    c(1,11) =  1.0;  c(2,11) =  1.0;  c(3,11) = -1.0
    c(1,12) = -1.0;  c(2,12) =  1.0;  c(3,12) = -1.0
    c(1,13) = -1.0;  c(2,13) = -1.0;  c(3,13) = -1.0
    c(1,14) =  1.0;  c(2,14) = -1.0;  c(3,14) = -1.0

     !> define weight coefficient
    do l = 0, pdm
      select case(l)
      case(0)    ; wf(l) = 2.0d0/9.0d0
      case(1:6)  ; wf(l) = 1.0d0/9.0d0
      case(7:14) ; wf(l) = 1.0d0/72.0d0
      end select
    enddo

  end subroutine init_D3Q15

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

end module D3Q15model
!*******************************************************************
