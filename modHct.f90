module hct
  implicit none
  private
  public :: interp
  !
contains
  subroutine calcBaryc(A, M, lambda, dansT)
    implicit none
    logical, intent(out) :: dansT
    double precision, intent(in) :: A(3, 2), M(2)
    double precision, intent(out) :: lambda(3)
    double precision :: det
    double precision, parameter :: tol = 1.0d-15
    integer :: i
    !
    det = (A(2, 2) - A(3, 2)) * (A(1, 1) - A(3, 1)) + &
         (A(3, 1) - A(2, 1)) * (A(1, 2) - A(3, 2))
    lambda(1) = ((A(2, 2) - A(3, 2)) * (M(1) - A(3, 1)) + &
         (A(3, 1) - A(2, 1)) * (M(2) - A(3, 2))) / det
    lambda(2) = ((A(3, 2) - A(1, 2)) * (M(1) - A(3, 1)) + &
         (A(1, 1) - A(3, 1)) * (M(2) - A(3, 2))) / det
    lambda(3) = 1.0d0 - lambda(1) - lambda(2)
    !
    dansT = .true.
    do i = 1, 3
       if(abs(lambda(i)) < tol) then
          lambda(i) = 0.0d0
       end if
       if(lambda(i) < 0.0d0 .or. lambda(i) > 1.0d0) then
          dansT = .false.
       end if
    end do
  end subroutine calcBaryc
  !
  !
  subroutine calcCoeff(foncT, p, q, u, a, b, c, d, e, g, omega)
    implicit none
    double precision, intent(in) :: foncT(3), p(3), q(3), u(3)
    double precision, intent(out) :: a(3), b(3), c(3), d(3)
    double precision, intent(out) :: e(3), g(3), omega
    integer :: i, j, k
    !
    do j = 1, 3
       a(j) = foncT(j)
       b(j) = a(j) + p(j) / 3.0d0
    end do
    do k = 1, 3
       c(k) = a(k) + q(k) / 3.0d0
       d(k) = (a(k) + b(k) + c(k)) / 3.0d0
    end do
    !
    do i = 1, 3
       j = mod(i, 3) + 1
       k = mod(j, 3) + 1
       g(i) = 2.5d-1 * ((2.0d0 * (d(k) + d(j)) + &
            (4.0d0 - 3.0d0*u(i)) * c(k) + (u(i) - 2.0d0) * a(k) + &
            (3.0d0*u(i) - 2.0d0) * b(j) - u(i) * a(j)))
    end do
    !
    do i = 1, 3
       j = mod(i, 3) + 1
       k = mod(j, 3) + 1
       e(k) = (d(k) + g(i) + g(j)) / 3.0d0
    end do
    !
    omega = (e(1) + e(2) + e(3)) / 3.0d0
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
  subroutine calcCoordTi(coordT, coordOmega, i, coordTi)
    implicit none
    integer, intent(in) :: i
    double precision, intent(in) :: coordT(3, 2), coordOmega(2)
    double precision, intent(out) :: coordTi(3, 2)
    integer :: j, k
    !
    j = mod(i, 3) + 1
    k = mod(j, 3) + 1
    coordTi(1, 1) = coordOmega(1)
    coordTi(1, 2) = coordOmega(2)
    coordTi(2, 1) = coordT(j, 1)
    coordTi(2, 2) = coordT(j, 2)
    coordTi(3, 1) = coordT(k, 1)
    coordTi(3, 2) = coordT(k, 2)
  end subroutine calcCoordTi
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
  subroutine calcOmega(A, Omega)
    implicit none
    double precision, intent(in) :: A(3, 2)
    double precision, intent(out) :: Omega(2)
    !
    Omega(1) = (A(1, 1) + A(2, 1) + A(3, 1)) / 3.0d0
    Omega(2) = (A(1, 2) + A(2, 2) + A(3, 2)) / 3.0d0
  end subroutine calcOmega
  !
  !
  subroutine calcpq(A, gradA, p, q)
    implicit none
    double precision, intent(in) :: A(3, 2), gradA(3, 2)
    double precision, intent(out) :: p(3), q(3)
    integer :: j, k
    !
    do j = 1, 3
       k = mod(j, 3) + 1
       p(j) = gradA(j, 1) * (A(k, 1) - A(j, 1)) + &
            gradA(j, 2) * (A(k, 2) - A(j, 2))
       q(k) = gradA(k, 1) * (A(j, 1) - A(k, 1)) + &
            gradA(k, 2) * (A(j, 2) - A(k, 2))
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
    do i = 1, 3
       j = mod(i, 3) + 1
       k = mod(j, 3) + 1
       !
       AjAk(1) = A(k, 1) - A(j, 1)
       AjAk(2) = A(k, 2) - A(j, 2)
       ! 
       OmegaAk(1) = A(k, 1) - Omega(1)
       OmegaAk(2) = A(k, 2) - Omega(2)
       !
       u(i) = 2.0d0 * (AjAk(1) * OmegaAk(1) + AjAk(2) * OmegaAk(2)) / &
            (AjAk(1) * AjAk(1) + AjAk(2) * AjAk(2))
    end do
  end subroutine calcu
  !
  !
  double precision function calcS(a, b, c, d, e, g, omega, lambda, i)
    double precision :: a(3), b(3), c(3), d(3), e(3), g(3), omega
    double precision :: lambda(3)
    integer :: i, j, k
    !
    j = mod(i, 3) + 1
    k = mod(j, 3) + 1
    calcS = a(k) * lambda(3)**3 + &
         c(k) * 3.0d0*lambda(2) * lambda(3)*lambda(3) + &
         b(j) * 3.0d0*lambda(2) * lambda(2) * lambda(3) + &
         a(j) * lambda(2)**3 + &
         d(k) * 3.0d0*lambda(1) * lambda(3) * lambda(3) + &
         g(i) * 6.0d0*lambda(1) * lambda(2) * lambda(3) + &
         d(j) * 3.0d0*lambda(1) * lambda(2)*lambda(2) + &
         e(k) * 3.0d0*lambda(1)*lambda(1) * lambda(3) + &
         e(j) * 3.0d0*lambda(1)*lambda(1) * lambda(2) + &
         omega * lambda(1)**3
  end function calcS
  !
  !
  subroutine interp(testPts, S, ntest)
    use donnees
    implicit none
    integer, intent(in) :: ntest
    double precision, intent(in) :: testPts(ntest, 2)
    double precision, intent(out) :: S(ntest)
    double precision :: M(2)
    double precision :: coordOmega(3), coordT(3, 2), coordTi(3, 2)
    double precision :: lambda(3)
    double precision :: foncT(3), gradT(3, 2)
    double precision :: a(3), b(3), c(3), d(3), e(3), g(3), omega
    double precision :: p(3), q(3), u(3)
    logical :: dansT, dansTi
    integer :: i, j, k, l, niter
    !
    M(1) = testPts(1, 1)
    M(2) = testPts(1, 2)
    !
    l = 0
    dansT = .false.
    do while (.not. dansT .and. l /= ntri)
       l = l + 1
       call calcCoordT(coord, n, tri(:, l), coordT)
       call calcBaryc(coordT, M, lambda, dansT)
    end do
    !
    if(dansT) then
       call calcFoncT(fonc, n, tri(:, l), foncT)
       call calcGradT(derivx, derivy, n, tri(:, l), gradT)
       call calcpq(coordT, gradT, p, q)
       call calcOmega(coordT, coordOmega)
       call calcu(coordT, coordOmega, u)
       call calcCoeff(foncT, p, q, u, a, b, c, d, e, g, omega)
       !
       i = 0
       dansTi = .false.
       do while (.not. dansTi)
          i = mod(i, 3) + 1
          call calcCoordTi(coordT, coordOmega, i, coordTi)
          call calcBaryc(coordTi, M, lambda, dansTi)
       end do
       !
       S(1) = calcS(a, b, c, d, e, g, omega, lambda, i)
       !
    else 
       write(6, '("Le point (", f15.7, ",", f15.7, ")")') &
            M(1), M(2)
       write(6, *) "n est pas dans la triangulation"
    end if
    !
    if (ntest > 1) then
       do k = 2, ntest
          M(1) = testPts(k, 1)
          M(2) = testPts(k, 2)
          !
          l = mod(l + ntri - 2, ntri) + 1
          niter = 0
          dansT = .false.
          do while (.not. dansT .and. niter /= ntri)
             l = mod(l, ntri) + 1
             call calcCoordT(coord, n, tri(:, l), coordT)
             call calcBaryc(coordT, M, lambda, dansT)
             niter = niter + 1
          end do
          !
          if(dansT) then
             if(niter /= 1) then
                call calcFoncT(fonc, n, tri(:, l), foncT)
                call calcGradT(derivx, derivy, n, tri(:, l), gradT)
                call calcpq(coordT, gradT, p, q)
                call calcOmega(coordT, coordOmega)
                call calcu(coordT, coordOmega, u)
                call calcCoeff(foncT, p, q, u, a, b, c, d, e, g, omega)
             end if
             dansTi = .false.
             i = mod(i + 1, 3) + 1
             do while (.not. dansTi)
                i = mod(i, 3) + 1
                call calcCoordTi(coordT, coordOmega, i, coordTi)
                call calcBaryc(coordTi, M, lambda, dansTi)
             end do
             S(k) = calcS(a, b, c, d, e, g, omega, lambda, i)
             !
          else 
             write(6, '("Le point (",f15.7, ",", f15.7, ")")') &
                  M(1), M(2)
             write(6, *) "n est pas dans la triangulation"
          end if
       end do
    end if
  end subroutine interp
end module hct
