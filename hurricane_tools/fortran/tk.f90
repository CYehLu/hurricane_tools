subroutine calc_tk(n, pres, theta, temp)
    ! calculate temperature (1d array)
    
    !f2py threadsafe
    implicit none
    
    integer, intent(in) :: n
    real(kind=8), dimension(n), intent(in) :: pres, theta
    real(kind=8), dimension(n), intent(out) :: temp
    
    integer :: i
    real(kind=8), parameter :: Rd = 287
    real(kind=8), parameter :: Cp = 1004.5
    real(kind=8), parameter :: Pbase = 100000
    
    !$omp parallel do schedule(runtime)
    do i = 1, n
        temp(i) = theta(i) * (pres(i) / Pbase) ** (Rd/Cp)
    end do
    !$omp end parallel do
    
    return
end subroutine calc_tk