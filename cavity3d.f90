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
  integer :: igrid, jgrid, kgrid, ngrid, ith, itm, its
  integer :: lb(14)
  ! 14 = number of particle - 1
  integer, pointer :: ibk(:,:)
  real(lbkind), pointer :: f(:,:,:,:)
  real(lbkind), pointer :: q(:,:,:,:)
  real(lbkind), pointer :: feq(:,:,:,:)
  real(lbkind), pointer :: x(:,:)
  real(lbkind), pointer :: y(:,:)
  real(lbkind), pointer :: s(:,:,:,:)
  real(lbkind), pointer :: r(:,:,:,:)
  real(lbkind) :: c(3,0:14), w(0:14) ! particle velocity and weight coefficient
  real(lbkind) :: tau, uc, re, omg, xc, yc
  real(4) :: t0, t1, t2
  integer :: n

  call define_c15
  call condition

  allocate( f(igrid,jgrid,kgrid,0:14) )
  allocate( q(igrid,jgrid,kgrid,5) )
  allocate( feq(igrid,jgrid,kgrid,0:14) )
  allocate( x(igrid,jgrid,kgrid) )
  allocate( y(igrid,jgrid,kgrid) )
  allocate( ibk(igrid,jgrid,kgrid) )
  allocate( s(igrid,jgrid,kgrid,14) )
  allocate( r(igrid,jgrid,kgrid,4) )

!  call output_xyzgrd
  call initialize
  
 ! call cpu_time(t0)
 ! t1 = t0
  
  do n = 1, nstep

    call equilibrium(1, igrid, 1, jgrid, 1,kgrid,feq)
    call collision(1, igrid, 1, jgrid,1,kgrid)
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
!     Lattice constants for the D3Q15 lattice                              
!***********************************************************************
  subroutine define_c15
    implicit none
    integer :: l

  ! D3Q15 model (for incompressible fluid)
  !
  !   8,12   2   7,11
  !             |
  !   3 -  0,5,6 - 1
  !             |
  !   9,13   4   10,14
  !

  ! define particle velocity c
    c(1, 0) =  0.0
    c(1, 1) =  1.0
    c(1, 2) =  0.0
    c(1, 3) = -1.0
    c(1, 4) =  0.0
    c(1, 5) =  0.0
    c(1, 6) =  0.0
    c(1, 7) =  1.0
    c(1, 8) = -1.0
	c(1, 9) = -1.0
	c(1,10) =  1.0
	c(1,11) =  1.0
	c(1,12) = -1.0
	c(1,13) = -1.0
	c(1,14) =  1.0

    c(2, 0) =  0.0
    c(2, 1) =  0.0
    c(2, 2) =  1.0
    c(2, 3) =  0.0
    c(2, 4) = -1.0
    c(2, 5) =  1.0
    c(2, 6) =  1.0
    c(2, 7) =  1.0
    c(2, 8) =  1.0
	c(2, 9) = -1.0
	c(2,10) = -1.0
	c(2,11) =  1.0
	c(2,12) =  1.0
	c(2,13) = -1.0
	c(2,14) = -1.0
	
	c(3, 0) =  0.0
	c(3, 1) =  0.0
	c(3, 2) =  0.0
	c(3, 3) =  0.0
	c(3, 4) =  0.0
	c(3, 5) =  1.0
	c(3, 6) = -1.0
	c(3, 7) =  1.0
	c(3, 8) =  1.0
	c(3, 9) =  1.0
	c(3,10) =  1.0
	c(3,11) = -1.0
	c(3,12) = -1.0
	c(3,13) = -1.0
	c(3,14) = -1.0

    w(0) = 2.0d0/ 9.0d0
    do l = 1, 6
      w(l) = 1.0d0/ 9.0d0
    enddo
    do l = 7, 14
      w(l) = 1.0d0/72.0d0
    enddo

    lb(1) = 3
    lb(2) = 4
    lb(3) = 1
    lb(4) = 2
    lb(5) = 6
    lb(6) = 5
    lb(7) = 13
    lb(8) = 14
	lb(9) = 11
	lb(10) = 12
	lb(11) = 9
	lb(12) = 10
	lb(13) =  7
	lb(14) =  8

  end subroutine

!***********************************************************************
!    Initialize the simulation at an equilibrium distribution                              
!***********************************************************************
  subroutine initialize
    implicit none
    integer :: i, j, k

    do j = 1, jgrid
    do i = 1, igrid
      x(i, j) = real(i-1)/real(ngrid-1)
      y(i, j) = real(j-1)/real(ngrid-1)
    enddo; enddo

  !++++ cavity flow ++++++++++++++++++++++++++++++++
	do k = 1, kgrid
    do j = 1, jgrid
    do i = 1, igrid
       q(i, j, k, 1) = 1.0
       q(i, j, k, 2) = 0.0
       q(i, j, k, 3) = 0.0
       q(i, j, k, 4) = 0.0
	   q(i, l, k, 5) = 1.0/3.0
    enddo; enddo;enddo
	
    do k= 1, kgrid
	do i = 1, igrid
       q(i, jgrid, k,1) = 1.0
       q(i, jgrid, k,2) = uc
       q(i, jgrid, k,3) = 0.0
	   q(i, jgrid, k,4) = 0.0
       q(i, jgrid, k,5) = 0.33
    enddo;enddo

  ! initialize distribution function
    call equilibrium(1, igrid, 1, jgrid, f)

  end subroutine

!***********************************************************************
!    Compute equilibrium distribution                           
!***********************************************************************
  subroutine equilibrium(ist, ien, jst, jen, kst, ken, f)
    implicit none
    integer :: ist, ien, jst, jen, kst, ken
    real(lbkind), pointer :: f(:,:,:)
    integer      :: i, j, k, l
    real(lbkind) :: uu, cu

	do k = kst, ken
    do j = ist, ien
    do i = jst, jen
      uu = q(i, j, 2) * q(i, j, 2) + q(i, j, 3) * q(i, j, 3) + q(i, j, 4) * q(i, j, 4)
      do l = 0, 14
        cu = c(1, l) * q(i, j, 2) + c(2, l) * q(i, j, 3) + c(3, l) * q(i, j, 4)
        f(i, j, l) = w(l) * q(i, j, 1) * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu)
      enddo

    enddo; enddo;enddo

  end subroutine

!***********************************************************************
!    Collision step
!***********************************************************************
  subroutine collision(ist, ien, jst, jen, kst, ken)
    implicit none
    integer :: ist, ien, jst, jen, kst, ken
    integer :: i, j, k, l

	do k = kst, ken
    do i = ist, ien
    do j = jst, jen
    do l = 0, 14
      f(i, j, k, l) = f(i, j, k, l) - 1.0/tau * (f(i, j, k, l) - feq(i, j, k, l))
    enddo; enddo; enddo; enddo

  end subroutine

!***********************************************************************
!    Translation step
!***********************************************************************
  subroutine translation
    implicit none
    integer :: i, j, k, l, ii, jj, kk, id, jd, kd
    real(lbkind), allocatable :: ft(:,:,:)

    allocate(ft(igrid,jgrid,kgrid))

    do l = 0, 14

      id = c(1,l)
      jd = c(2,l)
	  kd = c(3,l)

	  do k = 1, kgrid
      do j = 1, jgrid
      do i = 1, igrid
        ii = max(min(i-id,igrid),1)
        jj = max(min(j-jd,jgrid),1)
		kk = max(min(k-kd,kgrid),1)
        ft(i,j,k) = f(ii,jj,kk,l)
      enddo; enddo;enddo

	  do k= 1,kgrid
      do j = 1, jgrid
      do i = 1, igrid
        f(i,j,k,l) = ft(i,j,k)
      enddo; enddo;enddo

    enddo

    deallocate(ft)

  end subroutine
  
!***********************************************************************
!    Print out simulation parameters to screen
!***********************************************************************
  subroutine macro(n)
    implicit none
    integer      :: i, j, k, l, n
    real(lbkind) :: rho, u, v, w
	real(lbkind) :: resi, qref

    resi = 0.0d0

	do k = 1, kgrid
    do j = 1, jgrid
    do i = 1, igrid

      rho = 0.0d0
      u   = 0.0d0
      v   = 0.0d0
	! Calculation of density and momentum
      do l = 0, 14
        rho = rho + f(i, j, k,l)
        u   = u   + f(i, j, k,l) * c(1, l)
        v   = v   + f(i, j, k,l) * c(2, l)
		w  = w  + f(i,j,k,l) * c(3,l)
      enddo

      resi = resi + abs(q(i, j, k,2) - u / rho)**2.0

      q(i, j, 1) = rho
      q(i, j, 2) = u / rho  !> x - velocity
      q(i, j, 3) = v / rho  !> y - velocity
	  q(i, j, 4) = w / rho !> z - velocity
      q(i, j, 5) = rho / 3.0 !> pressure

      qref = qref + abs(q(i, j, k,2))**2.0

    enddo; enddo

    resi = sqrt( resi / real(igrid * jgrid*kgrid) )
    qref = sqrt( qref / real(igrid * jgrid*kgrid) )
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
    integer :: i, j,k, l
    real(lbkind) :: fsmach, alpha, time
	real(lbkind),allocatable :: qq(:,:,:,:)
    character(len=7) :: 

	allocate(qq(igrid,jgrid,kgrid,5))
    open(1, file='res.xyz', form='unformatted', status='unknown')
    write(1) igrid, jgrid, kgrid
    write(1) (((real(x(i, j, k)), i = 1, igrid), j = 1, jgrid),k=1,kgrid), &
            (((real(y(i, j, k)), i = 1, igrid), j = 1, jgrid),k=1,kgrid), &
			(((real(z(i, j, k)), i = 1, igrid), j = 1, jgrid),k=1,kgrid)
    close(1)

	do k = 1, kgrid
	do j = 1, jgrid
    do i = 1, igrid
	  qq(i,j,k,1) = q(i,j,k,1)
      qq(i, j, k,2) = q(i, j, k,1) * q(i, j, k,2)
      qq(i, j, k,3) = q(i, j, k,1) * q(i, j, k,3)
	  qq(i,j,k,4) = q(i,j,k,1) * q(i,j,k,4)
	  qq(i,j,k,5) = q(i,j,k,5)
    enddo; enddo;enddo

    write(cn,'(i7.7)') n
    open(2, file='res'//cn//'.q', form='unformatted', status='unknown')
    write(2) igrid, jgrid,kgrid
    write(2) fsmach, alpha, re, time
    write(2) ((((real(qq(i, j,k, l)), i = 1, igrid), j = 1, jgrid), k=1,kgrid), l = 1, 5)
    close(2)


!	do j = 1, jgrid
!    do i = 1, igrid
!      q(i, j, 2) = q(i, j, 2) / q(i, j, 1)
!      q(i, j, 3) = q(i, j, 3) / q(i, j, 1)
!    enddo; enddo

	deallocate(qq)
  end subroutine

  include'lbm_param.f90'
	
 end
