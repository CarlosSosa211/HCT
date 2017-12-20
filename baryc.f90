subroutine baryc(A, r, lambda)
  implicit none
  double precision :: A(3,2), lambda(3), r(2)
  double precision :: det

  det = (A(2,2)-A(3,2))*(A(1,1)-A(3,1))+(A(3,1)-A(2,1))*(A(1,2)-A(3,2))
  lambda(1) = ((A(2,2) - A(3,2)) * (r(1) - A(3,1)) + (A(3,1) - A(2,1)) * (r(2) - A(3,2))) / det
  lambda(2) = ((A(3,2) - A(1,2)) * (r(1) - A(3,2)) + (A(1,1) - A(3,1)) * (r(2) - A(3,2))) / det
  lambda(3) = 1 - lambda(1) - lambda(2)
  
end subroutine baryc
