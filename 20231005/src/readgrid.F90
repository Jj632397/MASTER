module mdl_readgrid
  use mdl_input_manager
  use mdl_decompo
  use mdl_mpisub_sbsp
  use mdl_param
  use mdl_refs
  use mdl_block
  implicit none

  !input------------------------------------------------------------------------
  character(100), private, save :: fname
  real(8)       , private, save :: scale
  !-----------------------------------------------------------------------------

contains

  subroutine readgrid(BLK)
    type(block), intent(inout) :: BLK

    integer :: myrank
    integer :: id_block
    integer :: imax,jmax,kmax
    integer :: i,j,k,ii,jj,kk

    real(8), allocatable, dimension(:,:,:,:) :: QB

    call input_readgrid(fname,scale)

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

    allocate(QB(imax+1,jmax+1,kmax+1,3))

    call mpisub_sbsp_readgrid_plot3d(fname,QB,imax+1,jmax+1,kmax+1,3,scale)

    do k=BLK%ksta,BLK%kend+1
      do j=BLK%jsta,BLK%jend+1
        do i=BLK%ista,BLK%iend+1
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1
          BLK%XCRD_NODE(i,j,k) = QB(ii,jj,kk,1)
          BLK%YCRD_NODE(i,j,k) = QB(ii,jj,kk,2)
          BLK%ZCRD_NODE(i,j,k) = QB(ii,jj,kk,3)
        enddo
      enddo
    enddo

    deallocate(QB)

  end subroutine readgrid

end module mdl_readgrid
