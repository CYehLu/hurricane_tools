    !subroutine find_level_n...(nx, ny, nz, zdata, nlev, levels, out)
    !    implicit none
    !
    !    ! arguments
    !    integer, intent(in) :: nx, ny, nz
    !    real(...), dimension(nx, ny, nz), intent(in) :: zdata
    !    integer, intent(in) :: nlev
    !    real(...), dimension(nlev), intent(in) :: levels
    !    integer, dimension(nx, ny, nlev), intent(out) :: out
    !
    !    ! local variables
    !    integer :: ilev
    !
    !    ! ----------------------------------------------------
        
        !$omp parallel do
        do ilev = 1, nlev
            call find_level_1(nx, ny, nz, zdata, levels(ilev), out(:,:,ilev))
        end do
        !$omp end parallel do

    !    return
    !end subroutine find_level_n...