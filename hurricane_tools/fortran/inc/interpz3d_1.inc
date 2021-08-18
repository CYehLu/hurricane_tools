    !subroutine interpz3d_1...(nx, ny, nz, var, lev_idx, w1, var_interp)
    !    implicit none
    !
    !    ! arguments
    !    integer, intent(in) :: nx, ny, nz
    !    real(...), dimension(nx, ny, nz), intent(in) :: var
    !    integer, dimension(nx, ny), intent(in) :: lev_idx
    !    real(...), dimension(nx, ny), intent(in) :: w1
    !    real(...), dimension(nx, ny), intent(out) :: var_interp
    !
    !    ! local variables
    !    integer :: i, j
    !    real(sp), parameter :: missing_val = -99999999.e0
    !
    !    ! ------------------------------------------------------

        do j = 1, ny
            do i = 1, nx
                if (w1(i,j) >= 0) then
                    var_interp(i,j) = w1(i,j) * var(i,j,lev_idx(i,j)) + (1-w1(i,j)) * var(i,j,lev_idx(i,j)+1)
                else
                    var_interp(i,j) = missing_val
                end if
            end do
        end do

    !    return
    !end subroutine interpz3d_1...