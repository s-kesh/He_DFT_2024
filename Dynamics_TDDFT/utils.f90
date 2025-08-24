module utils
    implicit none
    private
    public dscal, zscal, ddot, zdot, zdotc

    contains

        subroutine dscal(nx, ny, nz, scaler, vector)
            integer, intent(in) :: nx, ny, nz
            double precision, intent(in) :: scaler
            double precision, intent(inout) :: vector(nx, ny, nz)
            integer :: ix, iy, iz

            !$omp parallel do default(shared) private(ix,iy,iz) collapse(2)
            do iz=1, nz; do iy=1, ny
                do ix=1, nx
                    vector(ix,iy,iz) = vector(ix,iy,iz) * scaler
                end do
            end do; end do
            !$omp end parallel do
        end subroutine

        subroutine zscal(nx, ny, nz, scaler, vector)
            integer, intent(in) :: nx, ny, nz
            complex, intent(in) :: scaler
            complex, intent(inout) :: vector(nx, ny, nz)
            integer :: ix, iy, iz

            !$omp parallel do default(shared) private(ix,iy,iz) collapse(2)
            do iz=1, nz; do iy=1, ny
                do ix=1, nx
                    vector(ix,iy,iz) = vector(ix,iy,iz) * scaler
                end do
            end do; end do
            !$omp end parallel do
        end subroutine

        function ddot(nx, ny, nz, vector1, vector2) result(prod)
            integer, intent(in) :: nx, ny, nz
            double precision, intent(in) :: vector1(nx,ny,nz)
            double precision, intent(in) :: vector2(nx,ny,nz)
            double precision :: prod
            integer :: ix, iy, iz

            prod = 0.d0
            !$omp parallel do default(shared) private(ix,iy,iz) collapse(2) reduction(+:prod)
            do iz=1, nz; do iy=1, ny
                do ix=1, nx
                    prod = prod + vector1(ix,iy,iz)*vector2(ix,iy,iz)
                end do
            end do; end do
            !$omp end parallel do
        end function ddot

        function zdot(nx, ny, nz, vector1, vector2) result(prod)
            integer, intent(in) :: nx, ny, nz
            complex, intent(in) :: vector1(nx,ny,nz)
            complex, intent(in) :: vector2(nx,ny,nz)
            complex :: prod
            integer :: ix, iy, iz

            prod = (0.d0, 0.d0)
            !$omp parallel do default(shared) private(ix,iy,iz) collapse(2) reduction(+:prod)
            do iz=1, nz; do iy=1, ny
                do ix=1, nx
                    prod = prod + vector1(ix,iy,iz)*vector2(ix,iy,iz)
                end do
            end do; end do
            !$omp end parallel do
        end function zdot

        function zdotc(nx, ny, nz, vector1, vector2) result(prod)
            integer, intent(in) :: nx, ny, nz
            complex, intent(in) :: vector1(nx,ny,nz)
            complex, intent(in) :: vector2(nx,ny,nz)
            complex :: prod
            integer :: ix, iy, iz

            prod = (0.d0, 0.d0)
            !$omp parallel do default(shared) private(ix,iy,iz) collapse(2) reduction(+:prod)
            do iz=1, nz; do iy=1, ny
                do ix=1, nx
                    prod = prod + conjg(vector1(ix,iy,iz))*vector2(ix,iy,iz)
                end do
            end do; end do
            !$omp end parallel do
        end function zdotc

end module utils
