program test_hct
  use donnees
  use hct
  use testFunctions
  implicit none
  integer :: ntest, ntestx, ntesty
  double precision :: alpha, beta, gamma, delta, pasx, pasy
  double precision, allocatable :: err(:), S(:), testPts(:,:)
  integer :: i, j, nunit
  !
  call constrDonnees()
  call ecrFiDonnees()
  !
  nunit = 7
  open(unit = nunit, file = "testPoints.don")
  read(nunit, *) ntest
  allocate(S(ntest), testPts(ntest, 2))
  do i = 1, ntest
     read(nunit, *) testPts(i, 1), testPts(i, 2)
  end do
  close(nunit)
  !
  call interp(testPts, S, ntest)
  !
  open(unit = nunit, file = "hct.res", position = "append")
  do i = 1, ntest
     write(nunit, '("S(",f11.7,",",f11.7,") = ", f11.7)') &
          testPts(i, 1), testPts(i, 2), S(i)
  end do
  close(nunit)
  !
  deallocate(S, testPts)
  !
  call lecFiGrille(ntestx, ntesty, alpha, beta, gamma, delta)
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
  !
  call interp(testPts, S, ntest)
  !
  open(unit = nunit, file = "S.res")
  do i = 1, ntest
     err(i) = abs(f(testPts(i, 1), testPts(i, 2)) - S(i))
     write(nunit, '(2e15.7)') S(i), err(i)
  end do
  !
  open(unit = nunit, file = "hct.res", position = "append")
  write(nunit, '("Error_max = ", e15.7)') maxval(err)
  write(nunit, '("Error_min = ", e15.7)') minval(err)
  close(nunit)
  deallocate(S, testPts)
  call freeDonnees()
end program test_hct
