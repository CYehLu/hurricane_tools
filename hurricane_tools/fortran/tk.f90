subroutine calc_tk(nx, ny, nz, pres, theta, temp)
    ! from pressure and theta (potential temperature), calculate temperature
    ! (unit: K)
    
    implicit none
    integer,  intent(in) :: nx, ny, nz
    real(kind=8), dimension(nx, ny, nz), intent(in) :: pres, theta 
    real(kind=8), dimension(nx, ny, nz), intent(out) :: temp
    real(kind=8), parameter :: Rd = 287
    real(kind=8), parameter :: Cp = 1004.5
    real(kind=8), parameter :: Pbase = 100000

    temp = theta * (pres / Pbase) ** (Rd/Cp)
 
    return
end subroutine calc_tk


subroutine calc_tk_nd(nx, ny, nz, nt, pres, theta, temp)
    ! four dimension version of `calc_tk`
    
    use omp_lib
    !$ use omp_lib
    !f2py threadsafe
    
    implicit none
    
    integer, intent(in) :: nx, ny, nz, nt
    real(kind=8), dimension(nx, ny, nz, nt), intent(in) :: pres, theta
    real(kind=8), dimension(nx, ny, nz, nt), intent(out) :: temp
    
    integer :: it
    
    call omp_set_num_threads(16)
    !$omp parallel do
    do it = 1, nt
        call calc_tk(nx, ny, nz, pres(:,:,:,it), theta(:,:,:,it), temp(:,:,:,it))
    end do
    !$omp end parallel do
    
    return
end subroutine calc_tk_nd


