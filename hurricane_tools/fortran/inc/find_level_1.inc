    !subroutine find_level_1...(nx, ny, nz, zdata, level, out)
    !    implicit none
    !
    !    ! arguments
    !    integer, intent(in) :: nx, ny, nz
    !    real(...), dimension(nx, ny, nz), intent(in) :: zdata
    !    real(...), intent(in) :: level
    !    integer, dimension(nx, ny), intent(out) :: out
    !
    !    ! local variables
    !    integer :: i, j, k, prev_k
    !
    !    ! ----------------------------------------------------
    
        if (zdata(1,1,1) > zdata(1,1,nz)) then
            ! for zdata is descent order, like pressure
            ! find k such that : zdata(i,j,k) >= level > zdata(i,j,k+1)

            do i = 1, nx

                ! the first grid : upward search
                j = 1
                k = 1
                do while ((k <= nz) .and. (zdata(i,j,k) >= level))
                    k = k + 1
                end do
                out(i,j) = k - 1
                prev_k = k - 1

                do j = 2, ny  

                    ! the result of previous grid is fail (because of `level` out of range)
                    ! upward search
                    if ((prev_k == 0) .or. (prev_k == nz)) then
                        k = 1
                        do while ((k <= nz) .and. (zdata(i,j,k) >= level))
                            k = k + 1
                        end do
                        out(i,j) = k - 1
                        prev_k = k - 1

                    ! previous result can derectly apply to current grid
                    else if ((zdata(i,j,prev_k) >= level) .and. (zdata(i,j,prev_k+1) < level)) then 
                        out(i,j) = prev_k

                    ! downward adjustment (downward search)
                    else if (zdata(i,j,prev_k) < level) then
                        k = prev_k
                        do while ((k > 0) .and. (zdata(i,j,k) < level))
                            k = k - 1
                        end do
                        out(i,j) = k
                        prev_k = k

                    ! upward adjustment (upward search)
                    else if (zdata(i,j,prev_k+1) >= level) then
                        k = prev_k + 1
                        do while ((k <= nz) .and. (zdata(i,j,k) >= level))
                            k = k + 1
                        end do
                        out(i,j) = k - 1
                        prev_k = k - 1
                    end if

                end do
            end do

        else
            ! for zdata is ascent order, like height
            ! find `k` such that : zdata(i,j,k) <= level < zdata(i,j,k+1)
            ! I change the search direction (decent case : from bottom to top. here: from top to bottom),
            ! to be consistent with wrf-python

            do i = 1, nx

                ! the first grid : downward search
                j = 1
                k = nz
                do while ((k > 0) .and. (zdata(i,j,k) > level))
                    k = k - 1
                end do
                out(i,j) = k
                prev_k = k

                do j = 2, ny

                    ! the result of previous grid is fail (because of `level` out of range)
                    ! downward search
                    if ((prev_k == 0) .or. (prev_k == nz)) then
                        k = nz
                        do while ((k > 0) .and. (zdata(i,j,k) > level))
                            k = k - 1
                        end do
                        out(i,j) = k 
                        prev_k = k

                    ! previous result can derectly apply to current grid
                    else if ((zdata(i,j,prev_k) <= level) .and. (zdata(i,j,prev_k+1) > level)) then 
                        out(i,j) = prev_k

                    ! downward adjustment (downward search)
                    else if (zdata(i,j,prev_k) > level) then
                        k = prev_k
                        do while ((k > 0) .and. (zdata(i,j,k) > level))
                            k = k - 1
                        end do
                        out(i,j) = k
                        prev_k = k

                    ! upward adjustment (upward search)
                    else if (zdata(i,j,prev_k+1) <= level) then
                        k = prev_k + 1
                        do while ((k <= nz) .and. (zdata(i,j,k) <= level))
                            k = k + 1
                        end do
                        out(i,j) = k - 1
                        prev_k = k - 1
                    end if

                end do
            end do

        end if

    !    return 
    !end subroutine find_level_1...