    !subroutine calc_weights_1...(nx, ny, nz, zdata, level, lev_idx, w1)
    !    implicit none
    !
    !    ! arguments
    !    integer, intent(in) :: nx, ny, nz
    !    real(...), dimension(nx, ny, nz), intent(in) :: zdata
    !    real(...), intent(in) :: level
    !    integer, dimension(nx, ny), intent(in) :: lev_idx
    !    real(...), dimension(nx, ny), intent(out) :: w1
    !
    !    ! local variables
    !    integer :: i, j, idx
    !    real(...), parameter :: missing_val = -99999999.e0
    !    ! -------------------------------------------------
    
        do j = 1, ny
            do i = 1, nx
                if ((lev_idx(i,j) .ne. 0) .and. (lev_idx(i,j) .ne. nz)) then
                    idx = lev_idx(i,j)
                    w1(i,j) = (level - zdata(i,j,idx+1)) / (zdata(i,j,idx) - zdata(i,j,idx+1))
                else
                    w1(i,j) = missing_val
                end if
            end do
        end do

    !    return
    !end subroutine calc_weights_1....