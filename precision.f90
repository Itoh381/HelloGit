module precision_module
  implicit none

  integer, parameter :: sp = kind(0.0)
  integer, parameter :: dp = kind(0.0d0)

  private :: sp, dp

 ! integer, parameter :: lbkind = sp
   integer, parameter :: lbkind = dp

  public :: lbkind

end module precision_module
