!***********************************************************************
!     LATTICE BOLTZMANN METHOD                               03/12/04
!***********************************************************************
module precision_module
  implicit none
  integer, parameter :: dp = kind(0.0d0)
  private :: dp
  integer, parameter :: lbkind = dp
  public :: lbkind
end module precision_module

!***********************************************************************
!     The main program                              
!***********************************************************************
program lattice_Boltzmann_method
  use precision_module
  implicit none
  integer :: nstep
  integer :: igrid, jgrid, ngrid, ith, itm, its
  integer :: lb(8)
  integer, pointer :: ibk(:,:)
  real(lbkind), pointer :: f(:,:,:)
  real(lbkind), pointer :: q(:,:,:)
  real(lbkind), pointer :: feq(:,:,:)
  real(lbkind), pointer :: x(:,:)
  real(lbkind), pointer :: y(:,:)
  real(lbkind), pointer :: s(:,:,:)
  real(lbkind), pointer :: r(:,:,:)
  real(lbkind) :: c(2,0:8), w(0:8)
  real(lbkind) :: tau, uc, re, omg, xc, yc
  real(4) :: t0, t1, t2
  integer :: n

  call define_c
  call condition

  allocate( f(igrid,jgrid,0:8) )
  allocate( q(igrid,jgrid,4) )
  allocate( feq(igrid,jgrid,0:8) )
  allocate( x(igrid,jgrid) )
  allocate( y(igrid,jgrid) )
  allocate( ibk(igrid,jgrid) )
  allocate( s(igrid,jgrid,8) )
  allocate( r(igrid,jgrid,4) )

  call initialize
  
 ! call cpu_time(t0)
 ! t1 = t0
  
  do n = 1, nstep

    call equilibrium(1, igrid, 1, jgrid, feq)
    call collision(1, igrid, 1, jgrid)
    call translation
    call boundary(n)
    call macro(n)

!    if(mod(n,1) == 0) then
!		write(10,*) n, re
!		call cpu_time(t2)
!		write(6,'(a,i6,a,e14.7,a,f10.5,a)') ' step = ', n,', residual = ', re,', cpu time = ',t2-t1,'[s]'
!		t1 = t2
!	endif
	
	if( n >= 40000 .and. mod(n,100) == 0 )then
		call outdat(n)
	endif
	
  enddo
  call outdat(n)
  
!  call cpu_time(t1)
!  ith = aint((t1-t0)/3600.0)
!  itm = aint(amod((t1-t0),3600.0)/60.0)
!  its = amod(amod((t1-t0),3600.0),60.0)
	
!  write(6,'(a,3(i3,a))') 'total cpu time=  ',ith,'h',itm,'m',its,'s'
  
  contains

!***********************************************************************
!     Lattice constants for the D2Q9 lattice                              
!***********************************************************************
  subroutine define_c
    implicit none
    integer :: l

  ! D2Q9 model (for incompressible fluid)
  !
  !   6   2   5
  !       |
  !   3 - 0 - 1
  !       |
  !   7   4   8
  !

  ! define particle velocity c
    c(1, 0) =  0.0
    c(1, 1) =  1.0
    c(1, 2) =  0.0
    c(1, 3) = -1.0
    c(1, 4) =  0.0
    c(1, 5) =  1.0
    c(1, 6) = -1.0
    c(1, 7) = -1.0
    c(1, 8) =  1.0

    c(2, 0) =  0.0
    c(2, 1) =  0.0
    c(2, 2) =  1.0
    c(2, 3) =  0.0
    c(2, 4) = -1.0
    c(2, 5) =  1.0
    c(2, 6) =  1.0
    c(2, 7) = -1.0
    c(2, 8) = -1.0

    w(0) = 4.0/ 9.0
    do l = 1, 4
      w(l) = 1.0/ 9.0
    enddo
    do l = 5, 8
      w(l) = 1.0/36.0
    enddo

    lb(1) = 3
    lb(2) = 4
    lb(3) = 1
    lb(4) = 2
    lb(5) = 7
    lb(6) = 8
    lb(7) = 5
    lb(8) = 6

  end subroutine


!***********************************************************************
!    Initialize the simulation at an equilibrium distribution                              
!***********************************************************************
  subroutine initialize
    implicit none
    integer :: i, j

    do j = 1, jgrid
    do i = 1, igrid
      x(i, j) = real(i-1)/real(ngrid-1)
      y(i, j) = real(j-1)/real(ngrid-1)
    enddo; enddo

  !++++ cavity flow ++++++++++++++++++++++++++++++++
    do j = 1, jgrid
    do i = 1, igrid
       q(i, j, 1) = 1.0
       q(i, j, 2) = 0.0
       q(i, j, 3) = 0.0
       q(i, j, 4) = 0.0
    enddo; enddo

    do i = 1, igrid
       q(i, jgrid, 1) = 1.0
       q(i, jgrid, 2) = uc
       q(i, jgrid, 3) = 0.0
       q(i, jgrid, 4) = 0.33
    enddo

  ! initialize distribution function
    call equilibrium(1, igrid, 1, jgrid, f)

  end subroutine

!***********************************************************************
!    Compute equilibrium distribution                           
!***********************************************************************
  subroutine equilibrium(ist, ien, jst, jen, f)
    implicit none
    integer :: ist, ien, jst, jen
    real(lbkind), pointer :: f(:,:,:)
    integer      :: i, j, l
    real(lbkind) :: uu, cu

    do i = ist, ien
    do j = jst, jen
      uu = q(i, j, 2) * q(i, j, 2) + q(i, j, 3) * q(i, j, 3)

      do l = 0, 8
        cu = c(1, l) * q(i, j, 2) + c(2, l) * q(i, j, 3)
        f(i, j, l) = w(l) * q(i, j, 1) * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu)
      enddo

    enddo; enddo

  end subroutine

!***********************************************************************
!    Collision step
!***********************************************************************
  subroutine collision(ist, ien, jst, jen)
    implicit none
    integer :: ist, ien, jst, jen
    integer :: i, j, l

    do i = ist, ien
    do j = jst, jen
    do l = 0, 8
      f(i, j, l) = f(i, j, l) - 1.0/tau * (f(i, j, l) - feq(i, j, l))
    enddo; enddo; enddo

  end subroutine

!***********************************************************************
!    Translation step
!***********************************************************************
  subroutine translation
    implicit none
    integer :: i, j, l, ii, jj, id, jd
    real(lbkind), allocatable :: ft(:,:)

    allocate(ft(igrid,jgrid))

    do l = 0, 8

      id = c(1,l)
      jd = c(2,l)

      do j = 1, jgrid
      do i = 1, igrid
        ii = max(min(i-id,igrid),1)
        jj = max(min(j-jd,jgrid),1)
        ft(i,j) = f(ii,jj,l)
      enddo; enddo

      do j = 1, jgrid
      do i = 1, igrid
        f(i,j,l) = ft(i,j)
      enddo; enddo

    enddo

    deallocate(ft)

  end subroutine
  
!***********************************************************************
!    Print out simulation parameters to screen
!***********************************************************************
  subroutine macro(n)
    implicit none
    integer      :: i, j, l, n
    real(lbkind) :: rho, u, v, resi, qref

    resi = 0.0d0

    do j = 1, jgrid
    do i = 1, igrid

      rho = 0.0d0
      u   = 0.0d0
      v   = 0.0d0

      do l = 0, 8
        rho = rho + f(i, j, l)
        u   = u   + f(i, j, l) * c(1, l)
        v   = v   + f(i, j, l) * c(2, l)
      enddo

      resi = resi + abs(q(i, j, 2) - u / rho)**2.0

      q(i, j, 1) = rho
      q(i, j, 2) = u / rho
      q(i, j, 3) = v / rho
      q(i, j, 4) = rho / 3.0

      qref = qref + abs(q(i, j, 2))**2.0

    enddo; enddo

    resi = sqrt( resi / real(igrid * jgrid) )
    qref = sqrt( qref / real(igrid * jgrid) )
    resi = resi / qref

    !if (mod(n, 10).eq.0) write(10,*) n, log10(resi)

    if( mod(n,100) == 0 )then
      write(6,'(a,i6,a,e15.7)') '  step = ', n, '  residual = ', resi
    endif

  end subroutine

!***********************************************************************
!    Write the components of the velocity to a file with indices (x,y)
!***********************************************************************
  subroutine outdat(n)
    implicit none
    integer, intent(in) :: n
    integer :: i, j, l
    real(lbkind) :: fsmach, alpha, time
    character(len=7) :: cn

    open(1, file='res.xyz', form='unformatted', status='unknown')
    write(1) igrid, jgrid
    write(1) ((real(x(i, j)), i = 1, igrid), j = 1, jgrid) &
            ,((real(y(i, j)), i = 1, igrid), j = 1, jgrid)
    close(1)

    do i = 1, igrid
    do j = 1, jgrid
      q(i, j, 2) = q(i, j, 1) * q(i, j, 2)
      q(i, j, 3) = q(i, j, 1) * q(i, j, 3)
    enddo; enddo

    write(cn,'(i7.7)') n
    open(2, file='res'//cn//'.q', form='unformatted', status='unknown')
    write(2) igrid, jgrid
    write(2) fsmach, alpha, re, time
    write(2) (((real(q(i, j, l)), i = 1, igrid), j = 1, jgrid), l = 1, 4)
    close(2)

    do i = 1, igrid
    do j = 1, jgrid
      q(i, j, 2) = q(i, j, 2) / q(i, j, 1)
      q(i, j, 3) = q(i, j, 3) / q(i, j, 1)
    enddo; enddo

  end subroutine

  include'lbm_param.f90'
	
 end
