    !subroutine calc_weights_n...(nx, ny, nz, zdata, nlev, levels, lev_idx, w1)
    !    implicit none
    !
    !    ! arguments
    !    integer, intent(in) :: nx, ny, nz, nlev
    !    real(...), dimension(nx, ny, nz), intent(in) :: zdata
    !    real(...), dimension(nlev), intent(in) :: levels
    !    integer, dimension(nx, ny, nlev), intent(in) :: lev_idx
    !    real(...), dimension(nx, ny, nlev), intent(out) :: w1
    !
    !    ! local variables
    !    integer :: ilev
    !
    !    ! ---------------------------------------------------------

        !$omp parallel do
        do ilev = 1, nlev
            call calc_weights_1(nx, ny, nz, zdata, levels(ilev), lev_idx(:,:,ilev), w1(:,:,ilev))
        end do
        !$omp end parallel do

    !    return
    !end subroutine calc_weights_n....