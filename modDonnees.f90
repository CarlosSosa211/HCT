module donnees
  implicit none
  public :: constrDonnees, freeDonnees
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
       write(nunit, '(i3, 5e15.7)') i, coord(i, 1), coord(i, 2), &
            fonc(i), derivx(i), derivy(i)
    end do
    close(nunit)

  end subroutine ecrFiDonnees
  !
  !
  subroutine freeDonnees()
    deallocate(coord, fonc, derivx, derivy)
  end subroutine freeDonnees
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
    do i = 1, n
       read(nunit, *) coord(i, 1), coord(i, 1), coord(i,2)
    end do
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
    do j = 1, ntri
       read(nunit, *) tri(1, j), tri(1, j), tri(2, j), tri(3, j)
    end do
  end subroutine lecFiTri
end module donnees
