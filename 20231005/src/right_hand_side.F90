module mdl_rhs
  use mdl_input_manager
  use mdl_decompo
  use mdl_mpisub_sbsp
  use mdl_param
  use mdl_refs
  use mdl_block
  use mdl_halo
  use mdl_grad
  use mdl_flux_pw
  use mdl_flux_riemann_roe
  use mdl_flux_riemann_fvs_real_gas_extension
  use mdl_flux_riemann_fds_real_gas_extension
  use mdl_flux_viscous
  use mdl_sgs_wale
  use mdl_tbl_sa
  use mdl_schs_source
  implicit none

contains

  subroutine right_hand_side(BLK)
    type(block), intent(inout) :: BLK

    real(8), allocatable, dimension(:,:,:,:) :: RHSI
    real(8), allocatable, dimension(:,:,:,:) :: RHSV

    allocate(RHSI(BLK%imax_wh,BLK%jmax_wh,BLK%kmax_wh,neqns))
    allocate(RHSV(BLK%imax_wh,BLK%jmax_wh,BLK%kmax_wh,neqns))

    RHSI(:,:,:,:) = 0.0d0
    RHSV(:,:,:,:) = 0.0d0

    call halo_exchange_primitive(BLK)

    call gradient(BLK)

    call sgs_wale(BLK)

    if (iact_tbl_sa>0) call tbl_sa_result(BLK)

    call halo_exchange_gradient_and_tbl_coefficient(BLK)

    !---------------------------------------------------------------------------

!    call flux_pw(BLK,RHSI)
!    call flux_riemann_roe(BLK,RHSI)
    call flux_riemann_fvs_real_gas_extension(BLK,RHSI)
!    call flux_riemann_fds_real_gas_extension(BLK,RHSI)

    !---------------------------------------------------------------------------

    call flux_viscous(BLK,RHSV)

    if (iact_tbl_sa>0) call tbl_sa_source(BLK,RHSV)
    if (iact_passive_scalar>0) call schs_source(BLK,RHSV)

    BLK%RHSS(:,:,:,:) = RHSI(:,:,:,:)+RHSV(:,:,:,:)

    deallocate(RHSI)
    deallocate(RHSV)
  end subroutine right_hand_side
end module mdl_rhs
