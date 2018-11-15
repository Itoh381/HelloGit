!*********************************************************************
  module lbm
!*********************************************************************
  use precision_module
  use param
  use D3Q15model
  implicit none

  real(lbkind), allocatable :: f(:,:,:,:)
  real(lbkind), allocatable :: q(:,:,:,:)
  real(lbkind), allocatable :: ft(:,:,:,:)
  real(lbkind), allocatable :: fe(:,:,:,:)

  public :: initialize
  public :: output_xyzgrd
  public :: output_xyzq

  public :: boundary_a
  public :: collision
  public :: translation
  public :: macro
  public :: equilibrium

  contains

!*********************************************************************
  subroutine collision
!*********************************************************************
    integer      :: i,j,k,l
    real(lbkind) :: tau_i

    tau_i = 1.0d0 / tau
    do l = 0, pdm
    do k = g_kmin, g_kmax
    do j = g_jmin, g_jmax
    do i = g_imin, g_imax
      f(i,j,k,l) = f(i,j,k,l) - tau_i * (f(i,j,k,l) - fe(i,j,k,l))
    enddo; enddo; enddo; enddo

  end subroutine collision

!********************************************************************
  subroutine translation
!********************************************************************
    integer :: i,j,k,l,ii,jj,kk,id,jd,kd
    real(lbkind), allocatable :: ft(:,:,:)

    allocate(ft(g_imax,g_jmax,g_kmax))

    do l = 0 ,pdm

      id = c(1,l)
      jd = c(2,l)
      kd = c(3,l)

      do k = g_kmin, g_kmax
      do j = g_jmin, g_jmax
      do i = g_imin, g_imax
        ii = max(min(i-id,g_imax),g_imin)
        jj = max(min(j-jd,g_jmax),g_jmin)
        kk = max(min(k-kd,g_kmax),g_kmin)
        ft(i,j,k) = f(ii,jj,kk,l)
      enddo; enddo; enddo

      do k = g_kmin, g_kmax
      do j = g_jmin, g_jmax
      do i = g_imin, g_imax
        f(i,j,k,l) = ft(i,j,k)
      enddo; enddo; enddo

    enddo

    deallocate(ft)

  end subroutine translation

!*******************************************************************
  subroutine macro(resi)
!*******************************************************************
    integer :: i,j,k!,l
    !real(lblind) :: rho,u,v,w
    real(lbkind) :: ut!,ft(0:pdm)
    real(4) :: resi

    resi = 0.0d0

    do k = g_kmin+1, g_kmax-1
    do j = g_jmin+1, g_jmax-1
    do i = g_imin+1, g_imax-1
    !do k = g_kmin, g_kmax
    !do j = g_jmin, g_jmax
    !do i = g_imin, g_imax
      ut = q(i,j,k,2)
      !ft(0:pdm) = f(i,j,k,0:pdm)
      !q(i,j,k,1:5) = macro_(ft)
      q(i,j,k,1:5) = macro_(f(i,j,k,0:pdm))
      resi = resi + abs( ut - q(i,j,k,2) )/uc

      !rho = 0.0d0
      !u   = 0.0d0
      !v   = 0.0d0
      !w   = 0.0d0

      !!> calculation of density and momentum
      !do l = 0, pdm
      !   rho = rho + f(i,j,k,l)        !> density
      !   u   = u   + f(i,j,k,l)*c(1,l) !> x - momentum
      !   v   = v   + f(i,j,k,l)*c(2,l) !> y - momentum
      !   w   = w   + f(i,j,k,l)*c(3,l) !> z - momentum
      !enddo

      !resi = resi + abs(q(i,j,k,2) - u / rho) / uc

      !q(i,j,k,1) = rho
      !q(i,j,k,2) = u / rho         !> x - velocity
      !q(i,j,k,3) = v / rho         !> y - velocity
      !q(i,j,k,4) = w / rho         !> z - velocity
      !q(i,j,k,5) = rho / 3.d0      !> pressure

    enddo; enddo; enddo

    resi = resi/real(g_imax*g_jmax*g_kmax)

  end subroutine macro

!*******************************************************************
  subroutine equilibrium
!*******************************************************************
    integer      :: i,j,k!,l
    !real(lbkind) :: qt(5), ft(0:pdm)
    !real(lbkind) :: uu, cu

    do k = g_kmin, g_kmax
    do j = g_jmin, g_jmax
    do i = g_imin, g_imax
      !qt = q(i,j,k,1:5)
      !ft = equi_(qt)
      !fe(i,j,k,0:pdm) = ft(0:pdm)
      fe(i,j,k,0:pdm) = equi_(q(i,j,k,1:5))

      !uu = q(i,j,k,2)*q(i,j,k,2) + q(i,j,k,3)*q(i,j,k,3) + q(i,j,k,4)*q(i,j,k,4)
      !do l = 0, pdm
      !  cu = c(1,l)*q(i,j,k,2) + c(2,l)*q(i,j,k,3) + c(3,l)*q(i,j,k,4)
      !  feq(i,j,k,l) = q(i,j,k,1)*wf(l)*( 1.d0 + 3.d0*cu + 4.5d0*cu*cu - 1.5d0*uu )
      !enddo

    enddo; enddo; enddo

  end subroutine equilibrium

!*********************************************************************
  subroutine boundary_a
!*********************************************************************
    implicit none
    integer :: i,j,k,l,m,id,jd,kd
    real(4) :: resi

    do l=1,pdm

      id = c(1,l)
      jd = c(2,l)
      kd = c(3,l)
      m = bounce_back(l)

      i = g_imin
      do k = g_kmin, g_kmax
      do j = g_jmin, g_jmax
        if( (i-id)<g_imin .or. (i-id)>g_imax .or. &
            (j-jd)<g_jmin .or. (j-jd)>g_jmax .or. &
            (k-kd)<g_kmin .or. (k-kd)>g_kmax )then
          f(i,j,k,l) = f(i,j,k,m)
        endif
      enddo; enddo
      i = g_imax
      do k = g_kmin, g_kmax
      do j = g_jmin, g_jmax
        if( (i-id)<g_imin .or. (i-id)>g_imax .or. &
            (j-jd)<g_jmin .or. (j-jd)>g_jmax .or. &
            (k-kd)<g_kmin .or. (k-kd)>g_kmax )then
          f(i,j,k,l) = f(i,j,k,m)
        endif
      enddo; enddo

      j = g_jmin
      do k = g_kmin, g_kmax
      do i = g_imin, g_imax
        if( (i-id)<g_imin .or. (i-id)>g_imax .or. &
            (j-jd)<g_jmin .or. (j-jd)>g_jmax .or. &
            (k-kd)<g_kmin .or. (k-kd)>g_kmax )then
          f(i,j,k,l) = f(i,j,k,m)
        endif
      enddo; enddo

      k = g_kmin
      do j = g_jmin, g_jmax
      do i = g_imin, g_imax
        if( (i-id)<g_imin .or. (i-id)>g_imax .or. &
            (j-jd)<g_jmin .or. (j-jd)>g_jmax .or. &
            (k-kd)<g_kmin .or. (k-kd)>g_kmax )then
          f(i,j,k,l) = f(i,j,k,m)
        endif
      enddo; enddo
      k = g_kmax
      do j = g_jmin, g_jmax
      do i = g_imin, g_imax
        if( (i-id)<g_imin .or. (i-id)>g_imax .or. &
            (j-jd)<g_jmin .or. (j-jd)>g_jmax .or. &
            (k-kd)<g_kmin .or. (k-kd)>g_kmax )then
          f(i,j,k,l) = f(i,j,k,m)
        endif
      enddo; enddo

    enddo

    i = g_imin
    do k = g_kmin, g_kmax
    do j = g_jmin, g_jmax
       q(i,j,k,1:5) = macro_(f(i,j,k,0:pdm))
       q(i,j,k,2) = 0.0d0
       q(i,j,k,3) = 0.0d0
       q(i,j,k,4) = 0.0d0
    enddo; enddo
    i = g_imax
    do k = g_kmin, g_kmax
    do j = g_jmin, g_jmax
       q(i,j,k,1:5) = macro_(f(i,j,k,0:pdm))
       q(i,j,k,2) = 0.0d0
       q(i,j,k,3) = 0.0d0
       q(i,j,k,4) = 0.0d0
    enddo; enddo

    j = g_jmin
    do k = g_kmin, g_kmax
    do i = g_imin, g_imax
       q(i,j,k,1:5) = macro_(f(i,j,k,0:pdm))
       q(i,j,k,2) = 0.0d0
       q(i,j,k,3) = 0.0d0
       q(i,j,k,4) = 0.0d0
    enddo; enddo

    k = g_kmin
    do j = g_jmin, g_jmax
    do i = g_imin, g_imax
       q(i,j,k,1:5) = macro_(f(i,j,k,0:pdm))
       q(i,j,k,2) = 0.0d0
       q(i,j,k,3) = 0.0d0
       q(i,j,k,4) = 0.0d0
    enddo; enddo
    k = g_kmax
    do j = g_jmin, g_jmax
    do i = g_imin, g_imax
       q(i,j,k,1:5) = macro_(f(i,j,k,0:pdm))
       q(i,j,k,2) = 0.0d0
       q(i,j,k,3) = 0.0d0
       q(i,j,k,4) = 0.0d0
    enddo; enddo

    j = g_jmax
    do k = g_kmin, g_kmax
    do i = g_imin, g_imax
!       q(i,j,k,1) = 1.0d0
       q(i,j,k,1) = q(i,j-1,k,1)
       q(i,j,k,2) = uc
       q(i,j,k,3) = 0.0d0
       q(i,j,k,4) = 0.0d0
       f(i,j,k,0:pdm) = equi_(q(i,j,k,1:5))
    enddo; enddo

  end subroutine boundary_a
!***********************************************************
  subroutine initialize
!******************************************************************
   integer          :: i,j,k,l
   character(len=4) :: chr2

   allocate(  f(g_imax,g_jmax,g_kmax,0:pdm) )
   allocate( ft(g_imax,g_jmax,g_kmax,0:pdm) )
   allocate( fe(g_imax,g_jmax,g_kmax,0:pdm) )
   allocate(  q(g_imax,g_jmax,g_kmax,    5) )

   if( lcont == 0 )then

     do k = g_kmin, g_kmax
     do j = g_jmin, g_jmax
     do i = g_imin, g_imax
        q(i,j,k,1) = 1.0d0
        q(i,j,k,2) = 0.0d0
        q(i,j,k,3) = 0.0d0
        q(i,j,k,4) = 0.0d0
        q(i,j,k,5) = 1.0d0/3.0d0
     enddo; enddo; enddo

     j = g_jmax
     do k = g_kmin, g_kmax
     do i = g_imin, g_imax
        q(i,j,k,1) = 1.0d0
        q(i,j,k,2) = uc
        q(i,j,k,3) = 0.0d0
        q(i,j,k,4) = 0.0d0
        q(i,j,k,5) = 1.0d0/3.0d0
     enddo; enddo

     do k = g_kmin, g_kmax
     do j = g_jmin, g_jmax
     do i = g_imin, g_imax
       f(i,j,k,0:pdm) = equi_(q(i,j,k,1:5))
     enddo; enddo; enddo

   else

     write(chr2,'(i4)') lcont

     open(1,file = './output'//trim(adjustl(chr2))//'.dat', form='unformatted', status='unknown')
     do l= 0, pdm
       read(1) (((f(i,j,k,l), i = g_imin, g_imax), j = g_jmin, g_jmax), k = g_kmin, g_kmax)
     enddo
     close(1)

     do k = g_kmin, g_kmax
     do j = g_jmin, g_jmax
     do i = g_imin, g_imax
       q(i,j,k,1:5) = macro_(f(i,j,k,0:pdm))
     enddo; enddo; enddo

   endif

  end subroutine initialize

!*******************************************************************
  subroutine output_xyzgrd
!*******************************************************************
    integer :: i,j,k
    real(lbkind), allocatable :: x(:,:,:),y(:,:,:),z(:,:,:)

    allocate( x(g_imin:g_imax,g_jmin:g_jmax,g_kmin:g_kmax) )
    allocate( y(g_imin:g_imax,g_jmin:g_jmax,g_kmin:g_kmax) )
    allocate( z(g_imin:g_imax,g_jmin:g_jmax,g_kmin:g_kmax) )

    do k = g_kmin, g_kmax
    do j = g_jmin, g_jmax
    do i = g_imin, g_imax
      x(i,j,k) = real(i-1)/real(g_imax-1)
      y(i,j,k) = real(j-1)/real(g_jmax-1)
      z(i,j,k) = real(k-1)/real(g_kmax-1)
    enddo; enddo; enddo

    open(4,file ='./0.xyz',form = 'unformatted',status = 'unknown')
    write(4) g_imax,g_jmax,g_kmax
    write(4) (((x(i,j,k),i=g_imin,g_imax),j=g_jmin,g_jmax),k=g_kmin,g_kmax), &
             (((y(i,j,k),i=g_imin,g_imax),j=g_jmin,g_jmax),k=g_kmin,g_kmax), &
             (((z(i,j,k),i=g_imin,g_imax),j=g_jmin,g_jmax),k=g_kmin,g_kmax)

    close(4)

  end subroutine output_xyzgrd

!*******************************************************************
  subroutine output_xyzq(n)
!*******************************************************************
    integer                   :: i,j,k,l,n,m
    real(lbkind), allocatable :: qq(:,:,:,:)
    character(len=4)          :: chr2

    m = n/ostep + lcont

    write(chr2,'(i4)') m

    allocate(qq(g_imax,g_jmax,g_kmax,5))

    do k = g_kmin, g_kmax
    do j = g_jmin, g_jmax
    do i = g_imin, g_imax
      qq(i,j,k,1) = q(i,j,k,1)
      qq(i,j,k,2) = q(i,j,k,1) * q(i,j,k,2)
      qq(i,j,k,3) = q(i,j,k,1) * q(i,j,k,3)
      qq(i,j,k,4) = q(i,j,k,1) * q(i,j,k,4)
      qq(i,j,k,5) = q(i,j,k,5)
    enddo; enddo; enddo

    open(2,file ='./'//trim(adjustl(chr2))//'.q',form = 'unformatted',status = 'unknown')
    write(2) g_imax,g_jmax,g_kmax
    write(2) fsmach,alpha,re,time
    write(2) ((((qq(i,j,k,l),i=g_imin,g_imax),j=g_jmin,g_jmax),k=g_kmin,g_kmax),l=1,5)
    close(2)

!    open(3,file ='./output.dat',form = 'unformatted',status = 'unknown')
    open(3,file ='./output'//trim(adjustl(chr2))//'.dat',form = 'unformatted',status = 'unknown')
    do l= 0, pdm
      write(3) (((f(i,j,k,l),i=g_imin,g_imax),j=g_jmin,g_jmax),k=g_kmin,g_kmax)
    enddo
    close(3)

    !open(4,file ='./u'//trim(adjustl(chr2))//'.dat',status = 'unknown')
    !i = 1+(g_imax-1)/2
    !k = 1+(g_kmax-1)/2
    !do j = g_jmin, g_jmax
    !  write(4,*) y(i,j,k), q(i,j,k,2)
    !enddo
    !close(4)

    !open(4,file ='./v'//trim(adjustl(chr2))//'.dat',status = 'unknown')
    !j = 1+(g_jmax-1)/2
    !k = 1+(g_kmax-1)/2
    !do i = g_imin, g_imax
    !  write(4,*) x(i,j,k), q(i,j,k,3)
    !enddo
    !close(4)

    deallocate(qq)

  end subroutine output_xyzq

end module lbm
!***********************************************************
