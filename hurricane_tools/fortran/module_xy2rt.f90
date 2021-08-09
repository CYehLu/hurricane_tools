module xy2rt
    implicit none
    real*4 :: UNDEF = -9999.
        
contains   

    subroutine find_grid(nx, ny, lon, lat, clon, clat, idx_x, idx_y)
        ! clon/clat may not located at the grid point.
        ! find the grid point (idx_x, idx_y), which satisfies:
        !     lon(idx_y, idx_x) <= clon < lon(idx_y, idx_x+1)
        !     lat(idx_y, idx_x) <= clat < lat(idx_y+1, idx_x)
        implicit none
        integer, intent(in) :: nx, ny
        real*4, intent(in) :: lon(ny,nx), lat(ny,nx)
        real*4, intent(in) :: clon, clat
        integer, intent(out) :: idx_x, idx_y

        integer :: i, j
        logical :: isfind
        ! -------------------------------------------------------

        do i = 1, nx
            do j = 1, ny
                if ((lon(j,i) <= clon) .and. (clon < lon(j,i+1)) .and. (lat(j,i) <= clat) .and. (clat < lat(j+1,i))) then
                    idx_x = i
                    idx_y = j
                    isfind = .true.
                    exit
                end if 
            end do
        end do

        if (.not. isfind) then
            write(*,*) "[Warning] Failed in `find_grid`. Set idx_x=-999, idx_y=-999."
            idx_x = -999
            idx_y = -999
        end if

        return
    end subroutine find_grid
    
    
    subroutine get_grididx_and_weights(nx, ny, nrad, nthe, lon, lat, clon, clat, radius, thetas, dx, dy, grididx, weights)
        implicit none
        integer, intent(in) :: nx, ny, nrad, nthe
        real*4, intent(in) :: lon(ny,nx), lat(ny,nx)
        real*4, intent(in) :: clon, clat
        real*4, intent(in) :: radius(nrad), thetas(nthe)
        real*4, intent(in) :: dx, dy
        integer, intent(out) :: grididx(2,nrad,nthe)
        real*4, intent(out) :: weights(4,nrad,nthe)
        
        integer :: irad, ithe
        integer :: cidx_x, cidx_y
        real*4 :: grid_ratio_x, grid_ratio_y
        real*4 :: tmp
        integer :: di, dj
        real*4 :: iw, jw
        ! ------------------------------------------------------------------
        
        ! find (cidx_x, cidx_y):
        !     lon(cidx_y, cidx_x) <= clon < lon(cidx_y, cidx_x+1)
        !     lat(cidx_y, cidx_x) <= clat < lat(cidx_y+1, cidx_x)
        call find_grid(nx, ny, lon, lat, clon, clat, cidx_x, cidx_y)
        
        ! grid_ratio: 0 ~ 1 
        ! if grid_ratio_x is closer to 0, indicates clon is closer to lon(cidx_y,cidx_x)
        grid_ratio_x = (clon - lon(cidx_y,cidx_x)) / (lon(cidx_y,cidx_x+1) - lon(cidx_y,cidx_x))
        grid_ratio_y = (clat - lat(cidx_y,cidx_x)) / (lat(cidx_y+1,cidx_x) - lat(cidx_y,cidx_x))
        
        ! find grididx and weights
        do ithe = 1, nthe
            do irad = 1, nrad
                ! x direction: di & iw
                tmp = radius(irad) * cos(thetas(ithe)) / dx + grid_ratio_x
                if (tmp >= 0.) then
                    di = int(tmp)
                    iw = tmp - di
                else
                    di = int(tmp) - 1
                    iw = tmp - di
                end if
                
                ! y direction: dj & jw
                tmp = radius(irad) * sin(thetas(ithe)) / dy + grid_ratio_y
                if (tmp >= 0.) then
                    dj = int(tmp)
                    jw = tmp - dj
                else
                    dj = int(tmp) - 1
                    jw = tmp - dj
                end if
                
                grididx(1,irad,ithe) = cidx_y + dj
                grididx(2,irad,ithe) = cidx_x + di
                
                weights(1,irad,ithe) = (1-iw) * (1-jw)
                weights(2,irad,ithe) = iw * (1-jw)
                weights(3,irad,ithe) = iw * jw
                weights(4,irad,ithe) = (1-iw) * jw
                
            end do
        end do
        
        return
    end subroutine get_grididx_and_weights
    
    
    subroutine interp_xy2rt(nx, ny, nrad, nthe, grididx, weights, var, var_rt)
        implicit none
        integer, intent(in) :: nx, ny, nrad, nthe
        integer, intent(in) :: grididx(2,nrad,nthe)
        real*4, intent(in) :: weights(4,nrad,nthe)
        real*4, intent(in) :: var(ny,nx)
        real*4, intent(out) :: var_rt(nrad,nthe)
        
        integer :: irad, ithe
        integer :: idx_y, idx_x
        logical :: isout
        real*4 :: v1, v2, v3, v4
        real*4 :: w1, w2, w3, w4
        ! --------------------------------------------------------
        
        isout = .false.
        
        do ithe = 1, nthe
            do irad = 1, nrad
                idx_y = grididx(1,irad,ithe)
                idx_x = grididx(2,irad,ithe)
                
                if ((idx_y >= ny) .or. (idx_y <= 0) .or. (idx_x >= nx) .or. (idx_x <= 0)) then
                    isout = .true.
                    var_rt(irad,ithe) = UNDEF
                    cycle
                end if
                
                v1 = var(idx_y, idx_x)
                v2 = var(idx_y, idx_x+1)
                v3 = var(idx_y+1, idx_x+1)
                v4 = var(idx_y+1, idx_x)
                
                w1 = weights(1,irad,ithe)
                w2 = weights(2,irad,ithe)
                w3 = weights(3,irad,ithe)
                w4 = weights(4,irad,ithe)
                
                var_rt(irad,ithe) = w1*v1 + w2*v2 + w3*v3 + w4*v4
            
            end do
        end do
        
        if (isout) then
            write(*,*) "[Warning] Some values of `var_rt` are `UNDEF` because of extrapolation. Suggest to use smaller `radius`."
        end if
        
        return
    end subroutine interp_xy2rt
    

end module xy2rt
