program hct
  use donnees
  implicit none
  integer :: ntestx, ntesty
  double precision ::  alpha, beta, gamma, delta, pasx, pasy
  double precision, allocatable :: S(:), testPts(:,:)
  double precision :: coordT(3, 2), coordOmega(3)
  double precision :: foncT(3), gradT(3, 2)
  double precision :: a(3), b(3), c(3), d(3), e(3), g(3), omega
  double precision :: p(3), q(3), u(3)
  integer :: i, j, nunit
  !
  call constrDonnees()
  !
  ntestx = 30
  ntesty = 30
  alpha = 0.0d0
  beta = 3.0d0
  gamma = 0.0d0
  delta = 4.0d0
  pasx = (beta - alpha) / (ntestx - 1)
  pasy = (delta - gamma) / (ntesty - 1)
  allocate(S(ntestx * ntesty), testPts(ntestx * ntesty, 2))
  do i = 0, ntestx - 1
     do j = 0, ntesty - 1
        testPts(i*ntesty + 1 + j, 1) = alpha + i*pasx
        testPts(i*ntesty + 1 + j, 2) = gamma + j*pasy
     end do
  end do
  !
  !
  do i = 1, ntri
     call calcCoordT(coord, n, tri(:, i), coordT)
     call calcFoncT(fonc, n, tri(:, i), foncT)
     call calcGradT(derivx, derivy, n, tri(:, i), gradT)
     call calcpq(coordT, gradT, p, q)
     call calcu(coordT, coordOmega, u)
     call calcCoeff(foncT, p, q, u, a, b, c, d, e, g, omega)
     !S() = a(k)*lambda(3)**3 + c(k)*3*lambda(2)*lambda(3)*lambda(3)
     !+ b(j)*3*lambda(2)*lambda(2)*lambda(3) + a(j)*lambda(2)**3
     !+ d(k)*3*lambda(1)*lambda(3)*lambda(3)
     !+ g(i)*6*lambda(1)*lambda(2)*lambda(3)
     !+ d(j)*3*lambda(1)*lambda(2)*lambda(2)
     !+ e(k)*3*lambda(1)*lambda(1)*lambda(3)
     !+ e(j)*3*lambda(1)*lambda(1)*lambda(2) + omega*lambda(1)**3   
  end do
  deallocate(S, testPts)
  call free()
end program hct
!
!
subroutine calcCoeff(foncT, p, q, u, a, b, c, d, e, g, omega)
  implicit none
  double precision, intent(in) :: foncT(3), p(3), q(3), u(3)
  double precision, intent(out) :: a(3), b(3), c(3), d(3)
  double precision, intent(out) :: e(3), g(3), omega
  integer :: i, j, k
  !
  do i = 1, 3
     a(i) = foncT(i)
     b(i) = a(i) + p(i)/3
     c(i) = a(i) + q(i)/3
     d(i) = (a(i) + b(i) + c(i)) / 3
  end do
  !
  do i = 1, 3
     j = mod(i, 3) + 1
     k = mod(j, 3) + 1
     g(i) = 0.25 * ((2*(d(k) + d(j)) + (4 - 3*u(i))*c(k) + &
          (u(i) - 2)*a(k) + (3*u(i) - 2)*b(j) - u(i)*a(j)))
  end do
  !
  do k = 1, 3
     i = mod(k + 1, 3)
     j = mod(k + 2, 3)
     e(k) = (d(k) + g(i) + g(j)) / 3
  end do
  !
  omega = (e(1) + e(2) + e(3)) / 3
end subroutine calcCoeff
!
!
subroutine calcFoncT(fonc, n, trii, foncT)
  implicit none
  integer, intent(in) :: n
  integer, intent(in) :: trii(3)
  double precision, intent(in) :: fonc(n)
  double precision, intent(out) :: foncT(3)
  integer :: j
  !
  do j  = 1, 3
     foncT(j) = fonc(trii(j))
  end do
end subroutine calcFoncT
!
!
subroutine calcCoordT(coord, n, trii, coordT)
  implicit none
  integer, intent(in) :: n
  integer, intent(in) :: trii(3)
  double precision, intent(in) :: coord(n, 2)
  double precision, intent(out) :: coordT(3, 2)
  integer :: j
  !
  do j = 1, 3
     coordT(j, 1) = coord(trii(j), 1)
     coordT(j, 2) = coord(trii(j), 2)
  end do
end subroutine calcCoordT
!
!
subroutine calcGradT(derivx, derivy, n, trii, gradT)
  implicit none
  integer, intent(in) :: n
  integer, intent(in) :: trii(3)
  double precision, intent(in) :: derivx(n), derivy(n)
  double precision, intent(out) :: gradT(3, 2)
  integer :: j
  !
  do j = 1, 3
     gradT(j, 1) = derivx(trii(j))
     gradT(j, 2) = derivy(trii(j))
  end do
end subroutine calcGradT
!
!
subroutine calcpq(A, gradA, p, q)
  implicit none
  double precision, intent(in) :: A(3, 2), gradA(3, 2)
  double precision, intent(out) :: p(3), q(3)
  integer :: i, j
  !
  do i = 1, 3
     j = mod(i, 3) + 1
     p(i) = gradA(i, 1) * (A(j, 1) - A(i, 1)) + &
          gradA(i, 2) * (A(j, 2) - A(i, 2))
     q(i) = gradA(i, 1) * (A(i, 1) - A(j, 1)) + &
          gradA(i, 2) * (A(i, 2) - A(j, 2))
  end do
end subroutine calcpq
!
!
subroutine calcu(A, Omega, u)
  implicit none
  double precision, intent(in) :: A(3, 2), Omega(2)
  double precision, intent(out) :: u(3)
  double precision :: AjAk(2), OmegaAk(2)
  integer :: i, j, k
  !
  do i = 1,3
     j = mod(i, 3) + 1
     k = mod(j, 3) + 1
     !
     AjAk(1) = A(k, 1) - A(j, 1)
     AjAk(2) = A(k, 2) - A(j, 2)
     ! 
     OmegaAk(1) = A(k, 1) - Omega(1)
     OmegaAk(2) = A(k, 2) - Omega(2)
     !
     u(i) = 2 * (AjAk(1) * OmegaAk(1) + AjAk(2) * OmegaAk(2)) / &
          (AjAk(1) * AjAk(1) + AjAk(2) * AjAk(2))
  end do
end subroutine calcu
