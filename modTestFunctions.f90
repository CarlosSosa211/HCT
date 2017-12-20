module testFunctions
  implicit none
  public :: f, dxf, dyf
  !
contains
  double precision function f(x, y)
    implicit none
    double precision :: x, y
    !
    !f = exp(x + y)
    f = y*y*y - 2*x*y*y - 5*x*x*y + 10*x*y + 1
  end function f
  !
  !
  double precision function dxf(x, y)
    implicit none
    double precision :: x, y
    !
    !dxf = exp(x + y)
    dxf = -2*y*y -10*x*y + 10*y
  end function dxf
  !
  !
  double precision function dyf(x, y)
    implicit none
    double precision :: x, y
    !
    !dyf = exp(x + y)
    dyf = 3*y*y - 4*x*y - 5*x*x + 10*x
  end function dyf
end module testFunctions
