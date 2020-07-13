program test_choose
 
  implicit none
 
  write (*, '(i0)') choose (14_8, 3_8)
 
contains
 
  function factorial (n) result (res)
 
    implicit none
    integer(8), intent (in) :: n
    integer(8) :: res
    integer(8) :: i
 
    res = product ((/(i, i = 1, n)/))
 
  end function factorial
 
  function choose (n, k) result (res)
 
    implicit none
    integer(8), intent (in) :: n
    integer(8), intent (in) :: k
    integer(8) :: res
 
    res = factorial (n) / (factorial (k) * factorial (n - k))
 
  end function choose
 
end program test_choose
