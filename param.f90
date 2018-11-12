!********************************************************************
  module param
!********************************************************************
  use precision_module
  use D3Q15model
  implicit none

  integer      :: nstep,ostep,lcont
  real(lbkind) :: fsmach,alpha,time
  real(lbkind) :: uc,re,tau

  integer :: g_imin,g_imax,g_jmin,g_jmax,g_kmin,g_kmax
!  integer :: imin,imax,jmin,jmax,kmin,kmax

  contains

!*******************************************************************
  subroutine condition
!*******************************************************************
  character(len=20) :: cdmy

    open(1,file = 'condition.dat', status = 'unknown')
    !> free-stream velocity (normalized by particle velocity)
    read(1,*) cdmy,uc
    !> Reynolds number
    read(1,*) cdmy,re
    !> i direction number of grid
    read(1,*) cdmy,g_imax
    !> j direction number of grid
    read(1,*) cdmy,g_jmax
    !> k direction number of grid
    read(1,*) cdmy,g_kmax
    !> number of time step
    read(1,*) cdmy,nstep
    !> output interval step
    read(1,*) cdmy,ostep
    !> output file number to resume calculation
    read(1,*) cdmy,lcont
    close(1)

    g_imin = 1
    g_jmin = 1
    g_kmin = 1

!    imin = g_imin
!    imax = g_imax
!    jmin = g_jmin
!    jmax = g_jmax
!    kmin = g_kmin
!    kmax = g_kmax

    !> calculate the relaxation parameter
    tau = 0.5d0 * (6.0d0 * uc * real(g_jmax-1)/re + 1.0d0)

    write(6,'(a,f15.7)') '  reynolds number      = ', re
    write(6,'(a,f15.7)') '  relaxation time      = ', tau

  end subroutine condition

  end module
!************************************************************************
