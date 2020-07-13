program test_choose
 
  implicit none
 
  write (*, '(i0)') choose (140_16, 3_16)
 
contains
 
  function factorial (n) result (res)
 
    implicit none
    integer(16), intent (in) :: n
    integer(16) :: res
    integer(16) :: i
 
    res = product ((/(i, i = 1, n)/))
 
  end function factorial
 
  function choose (n, k) result (res)
 
    implicit none
    integer(16), intent (in) :: n
    integer(16), intent (in) :: k
    integer(16) :: res
 
    res = factorial (n) / (factorial (k) * factorial (n - k))
 
  end function choose
 
end program test_choose
