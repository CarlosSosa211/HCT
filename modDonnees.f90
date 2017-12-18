module donnees
  implicit none
  public :: constrDonnees, free
  integer :: n, ntri
  double precision, allocatable :: coord(:, :)
  double precision, allocatable :: fonc(:), derivx(:), derivy(:)
  integer, allocatable :: tri(:,:)
  !
contains
  subroutine calcfDeriv()
    use testFunctions
    implicit none
    integer :: i
    !
    allocate(fonc(n), derivx(n), derivy(n))
    do i = 1, n
       fonc(i) = f(coord(i, 1), coord(i, 2))
       derivx(i) = dxf(coord(i, 1), coord(i, 2))
       derivy(i) = dyf(coord(i, 1), coord(i, 2))
    end do
  end subroutine calcfDeriv
  !
  !
  subroutine constrDonnees()
    call lecFiPoints()
    call lecFiTri()
    call calcfDeriv()
    call ecrFiDonnees()
  end subroutine constrDonnees
  !
  !
  subroutine ecrFiDonnees()
    implicit none
    integer :: i, nunit
    !
    nunit = 7
    open(unit = nunit, file = "hct.res")
    do i = 1, ntri
       write(nunit, *) i, tri(1, i), tri(2, i), tri(3, i)
    end do
    !
    do i = 1, n
       write(nunit, *) i, coord(i, 1), coord(i, 2), fonc(i), &
            derivx(i), derivy(i)
    end do
    close(nunit)
 
  end subroutine ecrFiDonnees
  !
  !
  subroutine free()
    deallocate(coord, fonc, derivx, derivy)
  end subroutine free
  !
  !
  subroutine lecFiPoints()
    implicit none
    integer :: i, j, nunit
    !
    nunit = 7
    open(unit = nunit, file = "hct.pts")
    read(nunit, *) n
    allocate(coord(n, 2))
    read(nunit, *) ((coord(i, j), j = 1, 2), i = 1, n)
    close(nunit)
  end subroutine lecFiPoints
  !
  !
  subroutine lecFiTri()
    implicit none
    integer :: i, j, nunit
    !
    nunit = 7
    open(unit = nunit, file = "hct.tri")
    read(nunit, *) ntri
    allocate(tri(3, ntri))
    read(nunit, *) ((tri(i, j), i = 1, 3), j = 1, ntri)
  end subroutine lecFiTri
end module donnees
