module mdl_bound
  use mdl_decompo
  use mdl_mpisub_sbsp
  use mdl_param
  use mdl_refs
  use mdl_block
  use mdl_bound_wall
  use mdl_bound_inflow_outflow_tpipe
  implicit none

contains

  subroutine bound(BLK)
    type(block), intent(inout) :: BLK

    call bound_inflow_outflow_tpipe(BLK)

    call bound_wall(BLK)

  end subroutine bound

end module mdl_bound
