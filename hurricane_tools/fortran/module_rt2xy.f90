module rt2xy
    implicit none
    real*4, parameter :: PI = acos(-1.)
    real*4 :: UNDEF = -9999.
    
contains

    real*4 function great_circle_distance(lon1, lat1, lon2, lat2)
        ! Calculate great-circle distance between two longitude-latitude points
        ! by Haversine formula.
        ! Return the distance (km) between (lon1, lat1) and (lon2, lat2).
        ! reference:
        ! https://stackoverflow.com/questions/27928/
        implicit none
        real*4, intent(in) :: lon1, lat1, lon2, lat2
        real*4 :: rlon1, rlat1, rlon2, rlat2
        real*4 :: a, c
        real*4, parameter :: R = 6373.0
        
        ! degree to radian
        rlon1 = lon1 * PI / 180
        rlat1 = lat1 * PI / 180
        rlon2 = lon2 * PI / 180
        rlat2 = lat2 * PI / 180
        
        a = sin((rlat2-rlat1)/2)**2 + cos(rlat1) * cos(rlat2) * sin((rlon2-rlon1)/2)**2
        c = 2 * atan2(sqrt(a), sqrt(1-a))
        great_circle_distance = R * c
        return
    end function great_circle_distance
    
    
    integer function find_array_idx(n, arr, x)
        ! find `idx`:
        !     arr(idx) <= x < arr(idx+1)
        ! `arr` is ascending order.
        ! it assums that it can always find the result
        implicit none
        integer, intent(in) :: n
        real*4, intent(in) :: arr(n)
        real*4, intent(in) :: x
        integer :: idx
        
        do idx = 1, n-1
            if ((arr(idx) <= x) .and. (arr(idx+1) > x)) then
                exit
            end if
        end do
        
        find_array_idx = idx    ! if not found, idx = n
        
        return
    end function find_array_idx
    
    
    subroutine get_radtheidx_radthe(nx, ny, nrad, nthe, lon, lat, clon, clat, radius, thetas, radthe_idx, radthe)
        implicit none
        
        ! parameters
        integer, intent(in) :: nx, ny, nrad, nthe
        real*4, intent(in) :: lon(ny,nx), lat(ny,nx)
        real*4, intent(in) :: clon, clat
        real*4, intent(in) :: radius(nrad), thetas(nthe)  
        integer, intent(out) :: radthe_idx(2,ny,nx)
        real*4, intent(out) :: radthe(2,ny,nx)
        
        ! local variables
        real*4 :: thetas_new(nthe)
        integer :: i, j
        real*4 :: tmp_x, tmp_y
        real*4 :: rad, the
        
        ! ---------------------------------------------------------
        
        ! make it be ascent order. without it, `find_array_idx` will return wrong result.
        ! e.g. (unit in degree, just for convenience)
        ! thetas = (/ 50, 90, 130, 170, 210, 250, 290, 330, 10 /)
        ! thetas_new = (/ 50, 90, 130, 170, 210, 250, 290, 330, 370 /)
        do i = 1, nthe
            if (thetas(i) < thetas(1)) then
                thetas_new(i) = thetas(i) + 2*PI
            else
                thetas_new(i) = thetas(i)
            end if
        end do

        do i = 1, nx
            do j = 1, ny
                tmp_x = great_circle_distance(clon, clat, lon(j,i), clat) * sign(1., lon(j,i)-clon)
                tmp_y = great_circle_distance(clon, clat, clon, lat(j,i)) * sign(1., lat(j,i)-clat)
                the = atan2(tmp_y, tmp_x)
                rad = great_circle_distance(clon, clat, lon(j,i), lat(j,i))
                
                if (the < 0) then
                    the = the + 2*PI    ! range : -pi ~ pi => 0 ~ 2*pi
                end if
                
                if (the < thetas_new(1)) then
                    the = the + 2*PI    ! make it falls in the range of `thetas_new`
                end if
                                
                radthe_idx(1,j,i) = find_array_idx(nrad, radius, rad)
                radthe_idx(2,j,i) = find_array_idx(nthe, thetas_new, the)  
                
                radthe(1,j,i) = rad
                radthe(2,j,i) = the
            end do
        end do
        
        return
    end subroutine get_radtheidx_radthe
        
    
    subroutine interp_rt2xy(nx, ny, nrad, nthe, radius, thetas, radthe_idx, radthe, var_rt, var_xy)
        implicit none
        
        ! parameters
        integer, intent(in) :: nx, ny, nrad, nthe
        real*4, intent(in) :: radius(nrad), thetas(nthe)
        integer, intent(in) :: radthe_idx(2,ny,nx)
        real*4, intent(in) :: radthe(2,ny,nx)
        real*4, intent(in) :: var_rt(nrad,nthe)
        real*4, intent(out) :: var_xy(ny,nx)
        
        ! local variables
        integer :: i, j
        real*4 :: thetas_new(nthe)
        real*4 :: maxrad
        real*4 :: rad, the
        integer :: irad, ithe
        real*4 :: r1, r2, t1, t2
        real*4 :: v11, v21, v12, v22
        real*4 :: vr1, vr2
        
        ! --------------------------------------------
        
        do i = 1, nthe
            if (thetas(i) < thetas(1)) then
                thetas_new(i) = thetas(i) + 2*PI
            else
                thetas_new(i) = thetas(i)
            end if
        end do
        
        maxrad = maxval(radius)
        
        do i = 1, nx
            do j = 1, ny
            
                irad = radthe_idx(1,j,i)
                ithe = radthe_idx(2,j,i)
                
                rad = radthe(1,j,i)
                the = radthe(2,j,i)
                
                if (rad > maxrad) then
                    var_xy(j,i) = UNDEF
                    cycle
                end if
                
                r1 = radius(irad)
                r2 = radius(irad+1)
                t1 = thetas_new(ithe)
                !t2 = thetas_new(mod(ithe+1, nthe))
                t2 = thetas_new(mod(ithe, nthe)+1)
                
                if (ithe == nthe) then
                    t2 = t2 + 2*PI
                end if
                
                v11 = var_rt(irad, ithe)
                v21 = var_rt(irad+1, ithe)
                v12 = var_rt(irad, mod(ithe, nthe)+1)
                v22 = var_rt(irad+1, mod(ithe, nthe)+1)
                
                ! interpolate along radius-axis
                vr1 = (r2-rad)/(r2-r1) * v11 + (rad-r1)/(r2-r1) * v21
                vr2 = (r2-rad)/(r2-r1) * v12 + (rad-r1)/(r2-r1) * v22
                
                ! interpolate along azimuth-axis
                var_xy(j,i) = (t2-the)/(t2-t1) * vr1 + (the-t1)/(t2-t1) * vr2      
                
            end do
        end do
                
        return
    end subroutine interp_rt2xy
    
end module rt2xy