program test_hct
  use donnees
  use hct
  use testFunctions
  use testPoints
  implicit none
  character ::  fiPts*7, fiTri*7
  character :: fiRes*7, fiSRes*5
  character :: fiGrille*10, fiTestPts*14
  integer :: i, nunit
  !
  fiPts = "hct.pts"
  fiTri = "hct.tri"
  fiRes = "hct.res"
  fiSRes = "S.res"
  fiGrille = "grille.don"
  fiTestPts = "testPoints.don"
  !
  call constrDonnees(fiPts, fiTri)
  call ecrFiDonnees(fiRes)
  call lecFiTestPts(fiTestPts)
  call interp(testPts, S, ntest)
  nunit = 7
  open(unit = nunit, file = fiRes, position = "append")
  do i = 1, ntest
     write(nunit, '("S(",f11.7,",",f11.7,") = ", f11.7)') &
          testPts(i, 1), testPts(i, 2), S(i)
  end do
  close(nunit)
  call freeTestPts()
  !
  call constrGrille(fiGrille)
  call interp(testPts, S, ntest)
  open(unit = nunit, file = fiSRes)
  do i = 1, ntest
     err(i) = abs(f(testPts(i, 1), testPts(i, 2)) - S(i))
     write(nunit, '(2e15.7)') S(i), err(i)
  end do
  open(unit = nunit, file = fiRes, position = "append")
  write(nunit, '("Error_max = ", e15.7)') maxval(err)
  write(nunit, '("Error_min = ", e15.7)') minval(err)
  close(nunit)
  call freeTestPts()
  call freeDonnees()
end program test_hct
!
!
