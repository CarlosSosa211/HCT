module testPoints
  implicit none
  integer :: ntest
  double precision, allocatable :: err(:), S(:), testPts(:,:)
  public :: constrGrille, freeTestPts, lecFiTestPts
  private :: lecFiGrille
  !
contains
  subroutine constrGrille(fiGrille)
    implicit none
    character(len = *), intent(in) :: fiGrille
    integer :: ntestx, ntesty
    double precision :: alpha, beta, gamma, delta, pasx, pasy
    integer :: i, j
    !
    call lecFiGrille(fiGrille, ntestx, ntesty, alpha, beta, gamma, &
         delta)
    pasx = (beta - alpha) / (ntestx - 1)
    pasy = (delta - gamma) / (ntesty - 1)
    ntest = ntestx * ntesty
    allocate(err(ntest), S(ntest), testPts(ntest, 2))
    do i = 0, ntesty - 1
       do j = 0, ntestx - 1
          testPts(i*ntestx + 1 + j, 1) = alpha + j*pasx
          testPts(i*ntestx + 1 + j, 2) = gamma + i*pasy
       end do
    end do
  end subroutine constrGrille
  !
  !
  subroutine freeTestPts()
    deallocate(err, S, testPts)
  end subroutine freeTestPts
  !
  !
  subroutine lecFiGrille(fiGrille, ntestx, ntesty, alpha, beta, &
       gamma, delta)
    implicit none
    character(len = *), intent(in) :: fiGrille
    integer, intent (out):: ntestx, ntesty
    double precision, intent(out) :: alpha, beta, gamma, delta
    integer :: nunit
    !
    nunit = 7
    open(unit = nunit, file = fiGrille)
    read(nunit, *) ntestx, ntesty
    read(nunit, *) alpha, beta, gamma, delta
    close(nunit)
  end subroutine lecFiGrille
  !
  !
  subroutine lecFiTestPts(fiTestPts)
    implicit none
    character(len = *), intent(in) :: fiTestPts
    integer :: i, nunit
    !
    nunit = 7
    open(unit = nunit, file = fiTestPts)
    read(nunit, *) ntest
    allocate(err(ntest), S(ntest), testPts(ntest, 2))
    do i = 1, ntest
       read(nunit, *) testPts(i, 1), testPts(i, 2)
    end do
    close(nunit)
  end subroutine lecFiTestPts
end module testPoints
