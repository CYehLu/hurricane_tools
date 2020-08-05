!!
!! This fortran code is modified from wrf-python:
!! https://github.com/NCAR/wrf-python/blob/develop/fortran/wrf_user.f90#L355
!!
!! Chun-Yeh Lu, NTU, window120909120909@gmail.com
!! 2020-08-03
!! This file is compiled by:
!!     f2py -c --opt='-O3' --f90flags='-fopenmp' -lgomp slp.f90 -m f90slp
!!

subroutine dcomputeseaprs(nx, ny, nz, z, t, p, q, sea_level_pressure)
    use omp_lib
    !$ use omp_lib
    !f2py threadsafe

    !!
    !! calculate sea level pressure.
    !!
    !! nx, ny, nz: dimensional size
    !! z: potential height (m), t: temperature (K)
    !! p: pressure (Pa), q: water vapor mixing ratio (kg/kg)
    !! sea_level_pressure: sea level pressure (hPa)
    !!
    !! When sea_level_pressure is all -999 or -998, means that
    !! there is an error. 
    !! -999 indicates that there is at least one grid point that
    !! fail to find the zeta level where is 100 hPa above the
    !! surface.
    !! -998 indicates that `klo == khi`.
    !!
    
    implicit none

    ! some constants
    real(kind=8), parameter :: RD = 287.D0
    real(kind=8), parameter :: G = 9.81D0
    real(kind=8), parameter :: USSALR = 0.0065D0
    real(kind=8), parameter :: TC=273.16D0+17.5D0, PCONST=10000
    
    ! arguments
    integer, intent(in) :: nx, ny, nz
    real(kind=8), dimension(nx, ny, nz), intent(in) :: z, t, p, q
    real(kind=8), dimension(nx, ny), intent(out) :: sea_level_pressure
    
    ! local variables (in wrf-python version, they are also arguments)
    integer, dimension(nx, ny) :: level
    real(kind=8), dimension(nx, ny) :: t_surf, t_sea_level

    ! local variables
    integer :: i, j, k
    integer :: klo, khi
    integer :: errcnt
    real(kind=8) :: plo, phi, tlo, thi, zlo, zhi
    real(kind=8) :: p_at_pconst, t_at_pconst, z_at_pconst
    logical :: l1, l2, l3, found

    logical, parameter :: ridiculous_mm5_test=.True.


    ! Find least zeta level that is PCONST Pa above the surface.  We
    ! later use this level to extrapolate a surface pressure and
    ! temperature, which is supposed to reduce the effect of the
    ! diurnal heating cycle in the pressure field.

    !$omp parallel do
    do j = 1, ny
        do i = 1, nx
            level(i,j) = -1

            k = 1
            found = .FALSE.
            do while ((.not. found) .and. (k <= nz))
                if (p(i,j,k) < p(i,j,1) - PCONST) then
                    level(i,j) = k
                    found = .TRUE.
                end if
                k = k + 1
            end do
            
        end do
    end do
    !$omp end parallel do

    ! === ERROR occur ===
    ! If there is any grid point that fail to find the zeta level
    ! where is PCONST Pa above the surface.
    if (any(level == -1)) then
        sea_level_pressure = -999.
        return
    end if


    ! Get temperature PCONST Pa above surface. Use this to extrapolate
    ! the temperature at the surface and down to sea level.
    
    errcnt = 0
    
    !$omp parallel do
    do j = 1, ny
        do i = 1, nx
            klo = max(level(i,j)-1, 1)
            khi = min(klo+1, nz-1)

            ! === ERROR occur ===
            if (klo == khi) then
                errcnt = errcnt + 1
            end if

            plo = p(i,j,klo)
            phi = p(i,j,khi)
            tlo = t(i,j,klo) * (1.D0 + 0.608D0 * q(i,j,klo))
            thi = t(i,j,khi) * (1.D0 + 0.608D0 * q(i,j,khi))
            zlo = z(i,j,klo)
            zhi = z(i,j,khi)
            
            p_at_pconst = p(i,j,1) - PCONST
            t_at_pconst = thi - (thi-tlo) * log(p_at_pconst/phi) * log(plo/phi)
            z_at_pconst = zhi - (zhi-zlo) * log(p_at_pconst/phi) * log(plo/phi)
            t_surf(i,j) = t_at_pconst * (p(i,j,1)/p_at_pconst) ** (USSALR*RD/G)
            t_sea_level(i,j) = t_at_pconst + USSALR * z_at_pconst         
        end do
    end do
    !$omp end parallel do
    
    ! === ERROR occur ===
    if (errcnt > 0) then
        sea_level_pressure = -998.
        return
    end if


    ! Correction to the sea level temperature if both the surface and 
    ! sea level temperatures are 'too' hot.
    ! This scope is originally controled by `if (ridiculous_mm5_test)`,
    ! and the default option is `ridiculous_mm5_test=.true.`
    ! So here I remove the if statement to speed up  (although the benefit 
    ! is small).

    !if (ridiculous_mm5_test) then
    !$omp parallel do
    do j = 1, ny
        do i = 1, nx
            l1 = t_sea_level(i,j) < TC
            l2 = t_surf(i,j) <= TC
            l3 = .not. l1
            if (l2 .and. l3) then
                t_sea_level(i,j) = TC
            else
                t_sea_level(i,j) = TC - 0.005D0 * (t_surf(i,j)-TC)**2
            end if
        end do
    end do
    !$omp end parallel do
    !end if


    ! Finally: compute sea level pressure
    !$omp parallel do
    do j = 1, ny
        do i = 1, nx
            sea_level_pressure(i,j) = 0.01 * p(i,j,1) * exp( 2*G*z(i,j,1) /        &
                                          (RD * (t_sea_level(i,j) + t_surf(i,j))) )
        end do
    end do
    !$omp end parallel do
    
    ! sea_level_pressure = 0.01 * p(:,:,1) * exp( 2*G*z(:,:,1) / (RD * (t_sea_level + t_surf)) ) 

    return
end subroutine dcomputeseaprs


subroutine dcomputeseaprs_nt(nx, ny, nz, nt, z, t, p, q, sea_level_pressure)
    ! Three dimensional (x, y, t) version of `dcomputeseaprs`
    
    use omp_lib
    !$ use omp_lib
    !f2py threadsafe
    
    implicit none
    
    ! arguments
    integer, intent(in) :: nx, ny, nz, nt
    real(kind=8), dimension(nx, ny, nz, nt), intent(in) :: z, t, p, q
    real(kind=8), dimension(nx, ny, nt), intent(out) :: sea_level_pressure
    
    ! local variables
    integer :: it
    
    call omp_set_num_threads(16)
    !$omp parallel do
    do it = 1, nt
        call dcomputeseaprs(    &
            nx, ny, nz,         & 
            z(:,:,:,it), t(:,:,:,it), p(:,:,:,it), q(:,:,:,it),     & 
            sea_level_pressure(:,:,it)    &
        )
    end do
    !$omp end parallel do
    
    return
end subroutine dcomputeseaprs_nt
    