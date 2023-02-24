module analysis_m
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use            :: defined_types_m, only: datablock_t
    use            :: therm_m,         only: pbc
    implicit none

    public :: g_r

contains

    subroutine g_r(gr_mat, pos, parambox, switch_case)
        ! Notes
        ! gr_mat(1,:) -> valors de distancia
        ! gr_mat(2,:) -> numero de partícules a aquesta distancia
        implicit none
        ! In/Out variables
        real(kind=dp), dimension(:,:), intent(inout) :: gr_mat
        real(kind=dp), dimension(:,:), intent(in)    :: pos
        type(datablock_t), intent(in)                :: parambox
        integer(kind=i64), intent(in)                :: switch_case
        ! Internal variables
        integer(kind=i64), save                      :: i_ax, index_mat, j_ax, n_p, n_gdr
        real(kind=dp), save                          :: dr, dist, dv, ndg, dens
        real(kind=dp), parameter                     :: PI = 4.0_dp * atan(1.0_dp)
        real(kind=dp), dimension(3)                  :: rij

        select case (switch_case)
            case (1_i64)
                
                ! SWITCH = 1 => definim la memoria de la funció
                n_p = size(pos, dim=2, kind=i64)
                dr = parambox%rdf_max_dist / real(parambox%rdf_num_bins, kind=dp)
                dens = real(n_p, kind=dp) / (parambox%box_side ** 3)

                gr_mat(1,:) = [(real(i_ax, kind=dp)*dr, i_ax=1, parambox%rdf_num_bins)]
                gr_mat(2,:) = 0.0_dp

                n_gdr = 0_i64
            
            case (2_i64)
                
                ! SWITCH = 2 => Calculem n(r)
                n_gdr = n_gdr + 1_i64
                do j_ax = 1, n_p - 1
                    do i_ax = j_ax + 1, n_p
                        ! Calculem rij
                        rij = pos(:, j_ax) - pos(:, i_ax)
                        call pbc(x=rij, l_=parambox%box_side)
                        dist = norm2(rij)
                        
                        ! Apliquem el cutoff de maxima distancia
                        if (dist < parambox%rdf_max_dist) then
                            index_mat = int(dist/dr, kind=i64) + 1_i64
                            gr_mat(2, index_mat) = gr_mat(2, index_mat) + 2.0_dp
                        end if
                    end do
                end do
            
            case (3_i64)
        
                ! SWITCH = 3 => Calculem g(r)
                do i_ax = 1, parambox%rdf_num_bins
                    associate(gdr => gr_mat(2, i_ax))
                        dv = (((real(i_ax, kind=dp) + 1.0_dp) ** 3) - (real(i_ax, kind=dp) ** 3)) * (dr ** 3)
                        ndg = (4.0_dp / 3.0_dp) * pi * dv * dens
                        gdr = gdr / (real(n_p, kind=dp) * ndg * real(n_gdr, kind=dp))
                    end associate
                end do
            
            end select

    end subroutine g_r

end module analysis_m
