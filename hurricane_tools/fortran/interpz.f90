!!!
!!! v.1 2020/08/09, Chun-Yeh Lu, NTU, window120909120909@gmail.com
!!! v.2 2020/10/14, Chun-Yeh Lu, NTU, window120909120909@gmail.com
!!!

!!! These subroutines are called from `hurricane_tools.getvar.Interpz3d`, and it would
!!! much faster than `wrf.interpz3d` in wrf-python.

!!! subroutine list
!!! ---------------
!!!     find_level_1(nx, ny, nz, zdata, level, out)
!!!         find level index for interpolating variable on a vertical level
!!!     find_level_n(nx, ny, nz, zdata, nlev, levels, out)
!!!         find level index for interpolating variable on multiple vertical levels
!!!
!!!     calc_weights_1(nx, ny, nz, zdata, level, lev_idx, w1)
!!!         calculate weight for interpolating variable on a vertical level
!!!     calc_weights_n(nx, ny, nz, zdata, nlev, levels, lev_idx, w1)
!!!         calculate weight for interpolating variable on multiple vertical levels
!!!
!!!     interpz3d_1(nx, ny, nz, var, lev_idx, w1, var_interp)
!!!         interpolate variable on a vertical level
!!!     interpz3d_n(nx, ny, nz, var, nlev, lev_idx, w1, var_interp)
!!!         interpolate variable on multiple vertical levels

!!! In `find_level_1`, there are some similar paragraphs apear several times (upward search
!!! and downward search). I did not wrap them into subroutines or functions because of 
!!! efficiency.
!!! Although -O3 flag turns on -finline-functions flag, I still found that if I wraped
!!! these similar paragraphs into subroutines, execution speed would be slower than this
!!! version.
!!! So for efficiency reasons, I sacrificed the code readability.

!!! NOTE:
!!! When `find_level_1` and `find_level_n` used in python directly, the result should
!!! minus 1 because python index is start from 0 and fortran is from 1.



subroutine find_level_1(nx, ny, nz, zdata, level, out)
    !! find level index for interpolating variable on a vertical level
    !!
    !! intput
    !! ------
    !! nx, ny, nz : int
    !!     spatial dimension size
    !! zdata(nx, ny, nz) : real*8
    !!     vertical coordinate, e.g pressure or height
    !! level : real*8
    !!     interpolating level
    !!
    !! output
    !! ------
    !! out(nx, ny) : int
    !!     the level information. For example, if out(i, j) = 5, it means 
    !!     that :
    !!         zdata(i, j, 5) <= level < zdata(i, j, 6)
    !!     for descent order `zdata`, or
    !!         zdata(i, j, 5) >= level > zdata(i, j, 6)
    !!     for ascent order `zdata`.
    !!     if out(i,j) = 0 or nz, it means that the level is out of the range of zdata(i,j,:)
    
    implicit none
    
    !f2py threadsafe
    
    ! arguments
    integer, intent(in) :: nx, ny, nz
    real(kind=8), dimension(nx, ny, nz), intent(in) :: zdata
    real(kind=8), intent(in) :: level
    integer, dimension(nx, ny), intent(out) :: out
    
    ! local variables
    integer :: i, j, k, prev_k
    
    
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
    
    return 
end subroutine find_level_1


subroutine find_level_n(nx, ny, nz, zdata, nlev, levels, out)
    !! find level index for interpolating variable on multiple vertical levels
    !! It is similar to `find_level_1`, but is multiple levels instead
    !! of one level.
    !!
    !! intput
    !! ------
    !! nx, ny, nz : int
    !!     spatial dimension size
    !! zdata(nx, ny, nz) : real*8
    !!     vertical coordinate, e.g pressure or height
    !! nlev : int
    !!     number of interpolating levels
    !! levels(nlev) : real*8
    !!     interpolating levels
    !!
    !! output
    !! ------
    !! out(nx, ny, nlev) : int
    !!     the level information. For example, if out(i, j, ilev) = 5, 
    !!     it means that :
    !!         zdata(i, j, 5) <= levels(ilev) < zdata(i, j, 6)
    !!     for descent order `zdata`, or
    !!         zdata(i, j, 5) >= levels(ilev) > zdata(i, j, 6)
    !!     for ascent order `zdata`.   
    !!     if out(i,j,ilev) = 0 or nz, it means that the levels(ilev) is out of the range of zdata(i,j,:)
    
    implicit none
    !f2py threadsafe
    
    ! arguments
    integer, intent(in) :: nx, ny, nz
    real(kind=8), dimension(nx, ny, nz), intent(in) :: zdata
    integer, intent(in) :: nlev
    real(kind=8), dimension(nlev), intent(in) :: levels
    integer, dimension(nx, ny, nlev), intent(out) :: out
    
    ! local variables
    integer :: ilev
    
    !$omp parallel do
    do ilev = 1, nlev
        call find_level_1(nx, ny, nz, zdata, levels(ilev), out(:,:,ilev))
    end do
    !$omp end parallel do
    
    return
end subroutine find_level_n


subroutine calc_weights_1(nx, ny, nz, zdata, level, lev_idx, w1)
    !! calculate weight for interpolating variable on a vertical level
    !!
    !! input
    !! -----
    !! nx, ny, nz : int
    !!     spatial dimension size
    !! zdata(nx, ny, nz) : real*8
    !!     vertical coordinate, e.g pressure or height
    !! level : real*8
    !!     interpolating level
    !! lev_idx(nx, ny) : int
    !!     The level index, which satisfies :
    !!         zdata(i, j, lev_idx(i,j)) >= level > zdata(i, j, lev_idx(i,j)+1)
    !!     for descent order `zdata` (like pressure), or
    !!         zdata(i, j, lev_idx(i,j)) <= level < zdata(i, j, lev_idx(i,j)+1)
    !!     for ascent order `zdata` (like height).
    !!     This is the output variable of `find_level_1`
    !! 
    !! output
    !! ------
    !! w1(nx, ny) : real*8
    !!     weights for interpolation
    !!     if we want to interpolate variable at zdata=500 : 
    !!
    !!              zdata            variable value
    !!                530  --------  10
    !!                500  --------  ?
    !!                450  --------  20
    !!
    !!     so that the w1 = (500 - 530) / (450 - 530) = 3 / 8 = 0.375,
    !!     and the interpolated value is : ? = w1 * 20 + (1 - w1) * 10 = 13.75
    !!     w1 would between 0 ~ 1
    !!     if w1(i,j) = -99999999, indicates that `level` exceeds the range of `zdata(i,j,:)`
    
    implicit none
    
    ! arguments
    integer, intent(in) :: nx, ny, nz
    real(kind=8), dimension(nx, ny, nz), intent(in) :: zdata
    real(kind=8), intent(in) :: level
    integer, dimension(nx, ny), intent(in) :: lev_idx
    real(kind=8), dimension(nx, ny), intent(out) :: w1
    
    ! local variables
    integer :: i, j, idx

    do j = 1, ny
        do i = 1, nx
            if ((lev_idx(i,j) .ne. 0) .and. (lev_idx(i,j) .ne. nz)) then
                idx = lev_idx(i,j)
                w1(i,j) = (level - zdata(i,j,idx+1)) / (zdata(i,j,idx) - zdata(i,j,idx+1))
            else
                w1(i,j) = -99999999.
            end if
        end do
    end do
    
    return
end subroutine calc_weights_1


subroutine calc_weights_n(nx, ny, nz, zdata, nlev, levels, lev_idx, w1)
    !! calculate weight for interpolating variable on multiple vertical levels
    !! 
    !! input parameters are almost identical to `calc_weights_1`, but with additional
    !! dimension `nlev` for `levels`, `lev_idx` and `w1`.
    !! see `calc_weights_1`
    
    implicit none
    !f2py threadsafe
    
    ! arguments
    integer, intent(in) :: nx, ny, nz, nlev
    real(kind=8), dimension(nx, ny, nz), intent(in) :: zdata
    real(kind=8), dimension(nlev), intent(in) :: levels
    integer, dimension(nx, ny, nlev), intent(in) :: lev_idx
    real(kind=8), dimension(nx, ny, nlev), intent(out) :: w1
    
    ! local variables
    integer :: ilev
    
    !$omp parallel do
    do ilev = 1, nlev
        call calc_weights_1(nx, ny, nz, zdata, levels(ilev), lev_idx(:,:,ilev), w1(:,:,ilev))
    end do
    !$omp end parallel do

    return
end subroutine calc_weights_n


subroutine interpz3d_1(nx, ny, nz, var, lev_idx, w1, var_interp)
    !! interpolate variable on a vertical level
    !! 
    !! input
    !! -----
    !! nx, ny, nz : int
    !!     spatial dimension size
    !! var(nx, ny, nz) : real*8
    !!     the variable which would be interpolated on the specified level
    !! lev_idx(nx, ny) : int
    !!     level index. this is the output of `find_level_1`
    !! w1(nx, ny) : real*8
    !!     weight for interpolation. this is the output of `calc_weight_1`
    !!
    !! output
    !! ------
    !! var_interp(nx, ny) : real*8
    !!     interpolated variable.
    
    implicit none
    
    ! arguments
    integer, intent(in) :: nx, ny, nz
    real(kind=8), dimension(nx, ny, nz), intent(in) :: var
    integer, dimension(nx, ny), intent(in) :: lev_idx
    real(kind=8), dimension(nx, ny), intent(in) :: w1
    real(kind=8), dimension(nx, ny), intent(out) :: var_interp
    
    ! local variables
    integer :: i, j
    
    do j = 1, ny
        do i = 1, nx
            if (w1(i,j) >= 0) then
                var_interp(i,j) = w1(i,j) * var(i,j,lev_idx(i,j)) + (1-w1(i,j)) * var(i,j,lev_idx(i,j)+1)
            else
                var_interp(i,j) = -99999999.
            end if
        end do
    end do
    
end subroutine interpz3d_1
    
    
subroutine interpz3d_n(nx, ny, nz, var, nlev, lev_idx, w1, var_interp)
    !! interpolate variable on multiple vertical levels
    !! 
    !! input
    !! -----
    !! nx, ny, nz : int
    !!     spatial dimension size
    !! var(nx, ny, nz) : real*8
    !!     the variable which would be interpolated on specified levels
    !! nlev : int
    !!     dimension size of specified vertical levels
    !! lev_idx(nx, ny, nlev) : int
    !!     level index. this is the output of `find_level_n`
    !! w1(nx, ny, nlev) : real*8
    !!     weight for interpolation. this is the output of `calc_weight_n`
    !!
    !! output
    !! ------
    !! var_interp(nx, ny, nlev) : real*8
    !!     interpolated variable.
    
    implicit none
    !f2py threadsafe
    
    ! arguments
    integer, intent(in) :: nx, ny, nz, nlev
    real(kind=8), dimension(nx, ny, nz), intent(in) :: var
    integer, dimension(nx, ny, nlev), intent(in) :: lev_idx
    real(kind=8), dimension(nx, ny, nlev), intent(in) :: w1
    real(kind=8), dimension(nx, ny, nlev), intent(out) :: var_interp
    
    ! local variables
    integer :: ilev
    
    !$omp parallel do
    do ilev = 1, nlev
        call interpz3d_1(nx, ny, nz, var, lev_idx(:,:,ilev), w1(:,:,ilev), var_interp(:,:,ilev))
    end do
    !$omp end parallel do
    
    return
end subroutine interpz3d_n
    