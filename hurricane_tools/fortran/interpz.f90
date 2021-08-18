!!!
!!! v.1 2020/08/09, Chun-Yeh Lu, NTU, window120909120909@gmail.com
!!! v.2 2020/10/14, Chun-Yeh Lu, NTU, window120909120909@gmail.com
!!! v.3 2021/08/18, Chun-Yeh Lu, NTU, window120909120909@gmail.com
!!!

!!! These subroutines are called from `hurricane_tools.getvar.Interpz3d`, and it would
!!! be faster than `wrf.interpz3d` in wrf-python.

! name convention:
! ---------------
!     subroutine_name_<1/n>_r<32/64>
!         <1/n>: 
!             `1` indicates that there is only one interpolated vertical level
!             `n` indicates that there are multiple interpolated vertical levels
!         r<32/64>: 
!             `r32` indicates that the input float arguments are single precision (sp)
!             `r64` indicates that the input float arguments are double precision (dp)

! subroutine list (ignore `r32` or `r64` postfix):
! -----------------------------------------------
!     find_level_1(nx, ny, nz, zdata, level, out)
!         find level index for interpolating variable on a vertical level
!     find_level_n(nx, ny, nz, zdata, nlev, levels, out)
!         find level index for interpolating variable on multiple vertical levels
!
!     calc_weights_1(nx, ny, nz, zdata, level, lev_idx, w1)
!         calculate weight for interpolating variable on a vertical level
!     calc_weights_n(nx, ny, nz, zdata, nlev, levels, lev_idx, w1)
!         calculate weight for interpolating variable on multiple vertical levels
!
!     interpz3d_1(nx, ny, nz, var, lev_idx, w1, var_interp)
!         interpolate variable on a vertical level
!     interpz3d_n(nx, ny, nz, var, nlev, lev_idx, w1, var_interp)
!         interpolate variable on multiple vertical levels

!!! NOTE:
!!! ----
!!! (1)
!!! In `find_level_1`, there are some similar paragraphs apear several times (upward search
!!! and downward search). I did not wrap them into subroutines or functions because of 
!!! efficiency.
!!! Although -O3 flag turns on -finline-functions flag, I still found that if I wraped
!!! these similar paragraphs into subroutines, execution speed would be slower than this
!!! version.
!!! So for efficiency reasons, I sacrificed the code readability.
!!!
!!! (2)
!!! When `find_level_1` and `find_level_n` used in python directly, the result should
!!! minus 1 because python index is start from 0 and fortran is from 1.


module mod_interpz
    implicit none
    integer, parameter :: sp = kind(0.e0)
    integer, parameter :: dp = kind(0.d0)   
    
    ! it seems that f2py doesn't support function overloading...
    ! the interface here is useless at the python-end...
    
    interface find_level_1
        module procedure find_level_1_r32
        module procedure find_level_1_r64
    end interface
    
    interface find_level_n
        module procedure find_level_n_r32
        module procedure find_level_n_r64
    end interface
    
    interface calc_weights_1
        module procedure calc_weights_1_r32
        module procedure calc_weights_1_r64
    end interface
    
    interface calc_weights_n
        module procedure calc_weights_n_r32
        module procedure calc_weights_n_r64
    end interface
    
    interface interpz3d_1
        module procedure interpz3d_1_r32
        module procedure interpz3d_1_r64
    end interface
    
    interface interpz3d_n
        module procedure interpz3d_n_r32
        module procedure interpz3d_n_r64
    end interface
    
contains

    ! ====================================================
    ! ====================================================
    ! ------ subroutines for r32 (single precision) ------
    ! ====================================================
    ! ====================================================

    subroutine find_level_1_r32(nx, ny, nz, zdata, level, out)
        !! find level index for interpolating variable on a vertical level
        !!
        !! intput
        !! ------
        !! nx, ny, nz : int
        !!     spatial dimension size
        !! zdata(nx, ny, nz) : real*4
        !!     vertical coordinate, e.g pressure or height
        !! level : real*4
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
        real(sp), dimension(nx, ny, nz), intent(in) :: zdata
        real(sp), intent(in) :: level
        integer, dimension(nx, ny), intent(out) :: out

        ! local variables
        integer :: i, j, k, prev_k
        
        ! ---------------------------------------------------
        
        include "./inc/find_level_1.inc"

        return 
    end subroutine find_level_1_r32


    subroutine find_level_n_r32(nx, ny, nz, zdata, nlev, levels, out)
        !! find level index for interpolating variable on multiple vertical levels
        !! It is similar to `find_level_1`, but is multiple levels instead
        !! of one level.
        !!
        !! intput
        !! ------
        !! nx, ny, nz : int
        !!     spatial dimension size
        !! zdata(nx, ny, nz) : real*4
        !!     vertical coordinate, e.g pressure or height
        !! nlev : int
        !!     number of interpolating levels
        !! levels(nlev) : real*4
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
        real(sp), dimension(nx, ny, nz), intent(in) :: zdata
        integer, intent(in) :: nlev
        real(sp), dimension(nlev), intent(in) :: levels
        integer, dimension(nx, ny, nlev), intent(out) :: out

        ! local variables
        integer :: ilev
        
        ! ---------------------------------------------------
        
        include "./inc/find_level_n.inc"

        return
    end subroutine find_level_n_r32


    subroutine calc_weights_1_r32(nx, ny, nz, zdata, level, lev_idx, w1)
        !! calculate weight for interpolating variable on a vertical level
        !!
        !! input
        !! -----
        !! nx, ny, nz : int
        !!     spatial dimension size
        !! zdata(nx, ny, nz) : real*4
        !!     vertical coordinate, e.g pressure or height
        !! level : real*4
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
        !! w1(nx, ny) : real*4
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
        real(sp), dimension(nx, ny, nz), intent(in) :: zdata
        real(sp), intent(in) :: level
        integer, dimension(nx, ny), intent(in) :: lev_idx
        real(sp), dimension(nx, ny), intent(out) :: w1

        ! local variables
        integer :: i, j, idx
        real(sp), parameter :: missing_val = -99999999.e0

        ! --------------------------------------------------
        
        include "./inc/calc_weights_1.inc"

        return
    end subroutine calc_weights_1_r32


    subroutine calc_weights_n_r32(nx, ny, nz, zdata, nlev, levels, lev_idx, w1)
        !! calculate weight for interpolating variable on multiple vertical levels
        !! 
        !! input parameters are almost identical to `calc_weights_1`, but with additional
        !! dimension `nlev` for `levels`, `lev_idx` and `w1`.
        !! see `calc_weights_1`

        implicit none
        !f2py threadsafe

        ! arguments
        integer, intent(in) :: nx, ny, nz, nlev
        real(sp), dimension(nx, ny, nz), intent(in) :: zdata
        real(sp), dimension(nlev), intent(in) :: levels
        integer, dimension(nx, ny, nlev), intent(in) :: lev_idx
        real(sp), dimension(nx, ny, nlev), intent(out) :: w1

        ! local variables
        integer :: ilev
        
        ! -------------------------------------------------------
        
        include "./inc/calc_weights_n.inc"

        return
    end subroutine calc_weights_n_r32


    subroutine interpz3d_1_r32(nx, ny, nz, var, lev_idx, w1, var_interp)
        !! interpolate variable on a vertical level
        !! 
        !! input
        !! -----
        !! nx, ny, nz : int
        !!     spatial dimension size
        !! var(nx, ny, nz) : real*4
        !!     the variable which would be interpolated on the specified level
        !! lev_idx(nx, ny) : int
        !!     level index. this is the output of `find_level_1`
        !! w1(nx, ny) : real*4
        !!     weight for interpolation. this is the output of `calc_weight_1`
        !!
        !! output
        !! ------
        !! var_interp(nx, ny) : real*4
        !!     interpolated variable.

        implicit none

        ! arguments
        integer, intent(in) :: nx, ny, nz
        real(sp), dimension(nx, ny, nz), intent(in) :: var
        integer, dimension(nx, ny), intent(in) :: lev_idx
        real(sp), dimension(nx, ny), intent(in) :: w1
        real(sp), dimension(nx, ny), intent(out) :: var_interp

        ! local variables
        integer :: i, j
        real(sp), parameter :: missing_val = -99999999.e0
        
        ! ---------------------------------------------------
        
        include "./inc/interpz3d_1.inc"

        return
    end subroutine interpz3d_1_r32


    subroutine interpz3d_n_r32(nx, ny, nz, var, nlev, lev_idx, w1, var_interp)
        !! interpolate variable on multiple vertical levels
        !! 
        !! input
        !! -----
        !! nx, ny, nz : int
        !!     spatial dimension size
        !! var(nx, ny, nz) : real*4
        !!     the variable which would be interpolated on specified levels
        !! nlev : int
        !!     dimension size of specified vertical levels
        !! lev_idx(nx, ny, nlev) : int
        !!     level index. this is the output of `find_level_n`
        !! w1(nx, ny, nlev) : real*4
        !!     weight for interpolation. this is the output of `calc_weight_n`
        !!
        !! output
        !! ------
        !! var_interp(nx, ny, nlev) : real*4
        !!     interpolated variable.

        implicit none
        !f2py threadsafe

        ! arguments
        integer, intent(in) :: nx, ny, nz, nlev
        real(sp), dimension(nx, ny, nz), intent(in) :: var
        integer, dimension(nx, ny, nlev), intent(in) :: lev_idx
        real(sp), dimension(nx, ny, nlev), intent(in) :: w1
        real(sp), dimension(nx, ny, nlev), intent(out) :: var_interp

        ! local variables
        integer :: ilev
        
        ! ------------------------------------------------------
        
        include "./inc/interpz3d_n.inc"

        return
    end subroutine interpz3d_n_r32
    
    
    ! ====================================================
    ! ====================================================
    ! ------ subroutines for r64 (double precision) ------
    ! ====================================================
    ! ====================================================
    

    subroutine find_level_1_r64(nx, ny, nz, zdata, level, out)
        implicit none
        !f2py threadsafe

        ! arguments
        integer, intent(in) :: nx, ny, nz
        real(dp), dimension(nx, ny, nz), intent(in) :: zdata
        real(dp), intent(in) :: level
        integer, dimension(nx, ny), intent(out) :: out

        ! local variables
        integer :: i, j, k, prev_k
        
        ! ---------------------------------------------------
        
        include "./inc/find_level_1.inc"

        return 
    end subroutine find_level_1_r64


    subroutine find_level_n_r64(nx, ny, nz, zdata, nlev, levels, out)
        implicit none
        !f2py threadsafe

        ! arguments
        integer, intent(in) :: nx, ny, nz
        real(dp), dimension(nx, ny, nz), intent(in) :: zdata
        integer, intent(in) :: nlev
        real(dp), dimension(nlev), intent(in) :: levels
        integer, dimension(nx, ny, nlev), intent(out) :: out

        ! local variables
        integer :: ilev
        
        ! -------------------------------------------------------
        
        include "./inc/find_level_n.inc"

        return
    end subroutine find_level_n_r64


    subroutine calc_weights_1_r64(nx, ny, nz, zdata, level, lev_idx, w1)
        implicit none
        
        ! arguments
        integer, intent(in) :: nx, ny, nz
        real(dp), dimension(nx, ny, nz), intent(in) :: zdata
        real(dp), intent(in) :: level
        integer, dimension(nx, ny), intent(in) :: lev_idx
        real(dp), dimension(nx, ny), intent(out) :: w1

        ! local variables
        integer :: i, j, idx
        real(dp), parameter :: missing_val = -99999999.d0
        
        ! ---------------------------------------------------
        
        include "./inc/calc_weights_1.inc"

        return
    end subroutine calc_weights_1_r64


    subroutine calc_weights_n_r64(nx, ny, nz, zdata, nlev, levels, lev_idx, w1)
        implicit none
        !f2py threadsafe

        ! arguments
        integer, intent(in) :: nx, ny, nz, nlev
        real(dp), dimension(nx, ny, nz), intent(in) :: zdata
        real(dp), dimension(nlev), intent(in) :: levels
        integer, dimension(nx, ny, nlev), intent(in) :: lev_idx
        real(dp), dimension(nx, ny, nlev), intent(out) :: w1

        ! local variables
        integer :: ilev
        
        ! -----------------------------------------------------
        
        include "./inc/calc_weights_n.inc"
        
        return
    end subroutine calc_weights_n_r64


    subroutine interpz3d_1_r64(nx, ny, nz, var, lev_idx, w1, var_interp)
        implicit none

        ! arguments
        integer, intent(in) :: nx, ny, nz
        real(dp), dimension(nx, ny, nz), intent(in) :: var
        integer, dimension(nx, ny), intent(in) :: lev_idx
        real(dp), dimension(nx, ny), intent(in) :: w1
        real(dp), dimension(nx, ny), intent(out) :: var_interp

        ! local variables
        integer :: i, j
        real(dp), parameter :: missing_val = -99999999.d0
        
        ! ---------------------------------------------------
        
        include "./inc/interpz3d_1.inc"
        
        return
    end subroutine interpz3d_1_r64


    subroutine interpz3d_n_r64(nx, ny, nz, var, nlev, lev_idx, w1, var_interp)
        implicit none
        !f2py threadsafe

        ! arguments
        integer, intent(in) :: nx, ny, nz, nlev
        real(dp), dimension(nx, ny, nz), intent(in) :: var
        integer, dimension(nx, ny, nlev), intent(in) :: lev_idx
        real(dp), dimension(nx, ny, nlev), intent(in) :: w1
        real(dp), dimension(nx, ny, nlev), intent(out) :: var_interp

        ! local variables
        integer :: ilev
        
        ! ------------------------------------------------------
        
        include "./inc/interpz3d_n.inc"
        
        return
    end subroutine interpz3d_n_r64
    
end module mod_interpz