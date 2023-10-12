module mdl_restart
  use mdl_input_manager
  use mdl_decompo
  use mdl_mpisub_sbsp
  use mdl_param
  use mdl_refs
  use mdl_block
  use mdl_steps
  implicit none

  character(13), private, parameter :: prefix = "restart/rst_n"
  character(4 ), private, parameter :: suffix = ".dat"

  !input------------------------------------------------------------------------
  integer       , private, save :: step
  integer       , private, save :: mode !<0>:2d->2d or 3d->3d <1>:2d->3d
  !-----------------------------------------------------------------------------

contains

  subroutine restart(read_or_write,BLK)
    integer    , intent(in   ) :: read_or_write !<0>:read <1>:write
    type(block), intent(inout) :: BLK

    call input_restart(step,mode)

    if (read_or_write==0) then
      if (step>0.and.mode==0) then
        call read_restart(step,BLK)
      else if (step>0.and.mode==1) then
        call read_restart_from_2D_to_3D(step,BLK)
      endif
    else if (read_or_write==1) then
      call output_restart(BLK)
    endif

  end subroutine restart
  !-----------------------------------------------------------------------------
  subroutine output_restart(BLK)
    type(block), intent(in) :: BLK

    integer :: myrank
    integer :: id_block
    integer :: imax,jmax,kmax,lmax
    integer :: i,j,k,l,ii,jj,kk

    real(8), allocatable, dimension(:,:,:,:) :: QB

    integer :: nstep

    lmax = neqns

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

    allocate(QB(imax,jmax,kmax,lmax))

    do l=1,neqns
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            QB(ii,jj,kk,l) = BLK%QCS1(i,j,k,l)
          enddo
        enddo
      enddo
    enddo

    call steps_get_step_now(nstep)

    call mpisub_sbsp_out_restart(prefix,suffix,nstep,QB,imax,jmax,kmax,lmax)

    deallocate(QB)

  end subroutine output_restart
!-------------------------------------------------------------------------------
  subroutine read_restart(nstep,BLK)
    integer, intent(in) :: nstep
    type(block), intent(inout) :: BLK

    integer :: myrank
    integer :: id_block
    integer :: imax,jmax,kmax,lmax
    integer :: i,j,k,l,ii,jj,kk

    real(8), allocatable, dimension(:,:,:,:) :: QB

    lmax = neqns

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

    allocate(QB(imax,jmax,kmax,lmax))

    call mpisub_sbsp_read_restart(prefix,suffix,nstep,QB,imax,jmax,kmax,lmax)

    do l=1,neqns
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            BLK%QCS1(i,j,k,l) = QB(ii,jj,kk,l)
            BLK%QCS2(i,j,k,l) = BLK%QCS1(i,j,k,l)
            BLK%QCS3(i,j,k,l) = BLK%QCS1(i,j,k,l)
          enddo
        enddo
      enddo
    enddo

    deallocate(QB)

  end subroutine read_restart
!-------------------------------------------------------------------------------
  subroutine read_restart_from_2D_to_3D(nstep,BLK)
    integer, intent(in) :: nstep
    type(block), intent(inout) :: BLK

    integer :: myrank
    integer :: id_block
    integer :: imax,jmax,kmax,lmax
    integer :: i,j,k,l,ii,jj,kk

    real(8), allocatable, dimension(:,:,:,:) :: QB

    lmax = neqns

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

    allocate(QB(imax,jmax,kmax,lmax))

    call mpisub_sbsp_read_restart_from_2D_to_3D(prefix,suffix,nstep,QB,imax,jmax,kmax,lmax)

    do l=1,neqns
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            ii = i-BLK%ista+1
            jj = j-BLK%jsta+1
            kk = k-BLK%ksta+1
            BLK%QCS1(i,j,k,l) = QB(ii,jj,kk,l)
            BLK%QCS2(i,j,k,l) = BLK%QCS1(i,j,k,l)
            BLK%QCS3(i,j,k,l) = BLK%QCS1(i,j,k,l)
          enddo
        enddo
      enddo
    enddo

    deallocate(QB)

  end subroutine read_restart_from_2D_to_3D

end module mdl_restart
