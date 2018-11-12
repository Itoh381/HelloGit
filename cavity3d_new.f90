!********************************************************************
  program main
!********************************************************************
  use lbm
  implicit none
  integer :: n,i,j,k, ith, itm, its
  real(4) :: t0, t1, t2
  real(4) :: resi
!********************************************************************

  call condition

  call init_D3Q15

  call output_xyzgrd

  call initialize

  open(10,file='./resi.dat',status='unknown')

  call cpu_time(t0)

  t1 = t0

!> start loop
  do  n = 1, nstep

    call equilibrium

    call collision

!    call boundary_b
    call translation

    call boundary_a
!    call boundary

    call macro(resi)

    if( mod(n,1) == 0 ) then
      write(10,*) n, resi
      call cpu_time(t2)
      write(6,'(a,i6,a,e14.7,a,f10.5,a)') '  step = ', n, ',  residual = ', resi, ',  cpu time = ', t2-t1, '[s]'
      t1 = t2
    endif

    if( mod(n,ostep) == 0 )then
      call output_xyzq(n)
    endif

  enddo
!> end loop

  call cpu_time(t1)

  ith = aint((t1-t0)/3600.0)
  itm = aint(amod((t1-t0),3600.0)/60.0)
  its = amod(amod((t1-t0),3600.0),60.0)

  write(6,'(a,3(i3,a))') 'total cpu time=  ',ith,'h',itm,'m',its,'s'

  close(10)

  end program
!********************************************************************
