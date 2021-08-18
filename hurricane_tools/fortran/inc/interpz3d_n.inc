    !subroutine interpz3d_n...(nx, ny, nz, var, nlev, lev_idx, w1, var_interp)
    !    implicit none
    !    !f2py threadsafe
    !
    !    ! arguments
    !    integer, intent(in) :: nx, ny, nz, nlev
    !    real(...), dimension(nx, ny, nz), intent(in) :: var
    !    integer, dimension(nx, ny, nlev), intent(in) :: lev_idx
    !    real(...), dimension(nx, ny, nlev), intent(in) :: w1
    !    real(...), dimension(nx, ny, nlev), intent(out) :: var_interp
    !
    !    ! local variables
    !    integer :: ilev
    !
    !    ! ------------------------------------------------------------

        !$omp parallel do
        do ilev = 1, nlev
            call interpz3d_1(nx, ny, nz, var, lev_idx(:,:,ilev), w1(:,:,ilev), var_interp(:,:,ilev))
        end do
        !$omp end parallel do

    !    return
    !end subroutine interpz3d_n...