module mdl_metric
  use mdl_input_manager
  use mdl_decompo
  use mdl_mpisub_sbsp
  use mdl_param
  use mdl_refs
  use mdl_block
  use mdl_scmm
  use mdl_exchange_domain_halo
  implicit none

contains

  subroutine metric(BLK)
    type(block), intent(inout) :: BLK

    integer :: myrank
    integer :: id_block

    integer :: imax,jmax,kmax
    integer :: i,j,k,ii,jj,kk

    real(8), allocatable, dimension(:,:,:,:) :: XCRD_NODE_TEMP
    real(8), allocatable, dimension(:,:,:,:) :: XCRD_FACE_ETAZTA_TEMP
    real(8), allocatable, dimension(:,:,:,:) :: XCRD_FACE_ZTAXSI_TEMP
    real(8), allocatable, dimension(:,:,:,:) :: XCRD_FACE_XSIETA_TEMP

    real(8), allocatable, dimension(:,:,:,:) :: DXSX_FACE_TEMP
    real(8), allocatable, dimension(:,:,:,:) :: DETX_FACE_TEMP
    real(8), allocatable, dimension(:,:,:,:) :: DZTX_FACE_TEMP

    real(8), allocatable, dimension(:,:,:,:) :: XCRD_TEMP

    real(8), allocatable, dimension(:,:,:) :: RJCB_TEMP
    real(8), allocatable, dimension(:,:,:) :: GRSC_TEMP

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)

    allocate(XCRD_NODE_TEMP(imax+1,jmax+1,kmax+1,3))

    allocate(XCRD_FACE_ETAZTA_TEMP(imax+1,jmax,kmax,3))
    allocate(XCRD_FACE_ZTAXSI_TEMP(imax,jmax+1,kmax,3))
    allocate(XCRD_FACE_XSIETA_TEMP(imax,jmax,kmax+1,3))

    allocate(DXSX_FACE_TEMP(imax+1,jmax,kmax,3))
    allocate(DETX_FACE_TEMP(imax,jmax+1,kmax,3))
    allocate(DZTX_FACE_TEMP(imax,jmax,kmax+1,3))

    allocate(XCRD_TEMP(imax,jmax,kmax,3))

    allocate(RJCB_TEMP(imax,jmax,kmax))
    allocate(GRSC_TEMP(imax,jmax,kmax))

    XCRD_NODE_TEMP(:,:,:,:) = 0.0d0

    XCRD_FACE_ETAZTA_TEMP(:,:,:,:) = 0.0d0
    XCRD_FACE_ZTAXSI_TEMP(:,:,:,:) = 0.0d0
    XCRD_FACE_XSIETA_TEMP(:,:,:,:) = 0.0d0

    DXSX_FACE_TEMP(:,:,:,:) = 0.0d0
    DETX_FACE_TEMP(:,:,:,:) = 0.0d0
    DZTX_FACE_TEMP(:,:,:,:) = 0.0d0

    XCRD_TEMP(:,:,:,:) = 0.0d0

    RJCB_TEMP(:,:,:) = 0.0d0
    GRSC_TEMP(:,:,:) = 0.0d0

    do k=BLK%ksta,BLK%kend+1
      do j=BLK%jsta,BLK%jend+1
        do i=BLK%ista,BLK%iend+1
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1
          XCRD_NODE_TEMP(ii,jj,kk,1) = BLK%XCRD_NODE(i,j,k)
          XCRD_NODE_TEMP(ii,jj,kk,2) = BLK%YCRD_NODE(i,j,k)
          XCRD_NODE_TEMP(ii,jj,kk,3) = BLK%ZCRD_NODE(i,j,k)
        enddo
      enddo
    enddo

    call node_to_cell_face(XCRD_NODE_TEMP(:,:,:,1),imax+1,jmax+1,kmax+1, &
                           XCRD_FACE_ETAZTA_TEMP(:,:,:,1), &
                           XCRD_FACE_ZTAXSI_TEMP(:,:,:,1), &
                           XCRD_FACE_XSIETA_TEMP(:,:,:,1))

    call node_to_cell_face(XCRD_NODE_TEMP(:,:,:,2),imax+1,jmax+1,kmax+1, &
                           XCRD_FACE_ETAZTA_TEMP(:,:,:,2), &
                           XCRD_FACE_ZTAXSI_TEMP(:,:,:,2), &
                           XCRD_FACE_XSIETA_TEMP(:,:,:,2))

    call node_to_cell_face(XCRD_NODE_TEMP(:,:,:,3),imax+1,jmax+1,kmax+1, &
                           XCRD_FACE_ETAZTA_TEMP(:,:,:,3), &
                           XCRD_FACE_ZTAXSI_TEMP(:,:,:,3), &
                           XCRD_FACE_XSIETA_TEMP(:,:,:,3))

    call node_to_cell_center(XCRD_NODE_TEMP(:,:,:,1),imax+1,jmax+1,kmax+1, &
                             XCRD_TEMP(:,:,:,1))

    call node_to_cell_center(XCRD_NODE_TEMP(:,:,:,2),imax+1,jmax+1,kmax+1, &
                             XCRD_TEMP(:,:,:,2))

    call node_to_cell_center(XCRD_NODE_TEMP(:,:,:,3),imax+1,jmax+1,kmax+1, &
                             XCRD_TEMP(:,:,:,3))

    call scmm(XCRD_NODE_TEMP(:,:,:,1),XCRD_NODE_TEMP(:,:,:,2),XCRD_NODE_TEMP(:,:,:,3), &
              imax+1,jmax+1,kmax+1, &
              DXSX_FACE_TEMP(:,:,:,1),DXSX_FACE_TEMP(:,:,:,2),DXSX_FACE_TEMP(:,:,:,3), &
              DETX_FACE_TEMP(:,:,:,1),DETX_FACE_TEMP(:,:,:,2),DETX_FACE_TEMP(:,:,:,3), &
              DZTX_FACE_TEMP(:,:,:,1),DZTX_FACE_TEMP(:,:,:,2),DZTX_FACE_TEMP(:,:,:,3), &
              RJCB_TEMP)

    call grid_scale(XCRD_NODE_TEMP(:,:,:,1),XCRD_NODE_TEMP(:,:,:,2),XCRD_NODE_TEMP(:,:,:,3), &
                    imax+1,jmax+1,kmax+1,GRSC_TEMP)

    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1
          BLK%XCRD(i,j,k) = XCRD_TEMP(ii,jj,kk,1)
          BLK%YCRD(i,j,k) = XCRD_TEMP(ii,jj,kk,2)
          BLK%ZCRD(i,j,k) = XCRD_TEMP(ii,jj,kk,3)
          BLK%RJCB(i,j,k) = RJCB_TEMP(ii,jj,kk)
          BLK%GRSC(i,j,k) = GRSC_TEMP(ii,jj,kk)
        enddo
      enddo
    enddo

    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend+1
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1
          BLK%XCRD_FACE_ETAZTA(i,j,k) = XCRD_FACE_ETAZTA_TEMP(ii,jj,kk,1)
          BLK%YCRD_FACE_ETAZTA(i,j,k) = XCRD_FACE_ETAZTA_TEMP(ii,jj,kk,2)
          BLK%ZCRD_FACE_ETAZTA(i,j,k) = XCRD_FACE_ETAZTA_TEMP(ii,jj,kk,3)
          BLK%DXSX_FACE(i,j,k) = DXSX_FACE_TEMP(ii,jj,kk,1)
          BLK%DXSY_FACE(i,j,k) = DXSX_FACE_TEMP(ii,jj,kk,2)
          BLK%DXSZ_FACE(i,j,k) = DXSX_FACE_TEMP(ii,jj,kk,3)
          if (i2dimension>0) then
            BLK%DXSZ_FACE(i,j,k) = 0.0d0
          endif
        enddo
      enddo
    enddo

    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend+1
        do i=BLK%ista,BLK%iend
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1
          BLK%XCRD_FACE_ZTAXSI(i,j,k) = XCRD_FACE_ZTAXSI_TEMP(ii,jj,kk,1)
          BLK%YCRD_FACE_ZTAXSI(i,j,k) = XCRD_FACE_ZTAXSI_TEMP(ii,jj,kk,2)
          BLK%ZCRD_FACE_ZTAXSI(i,j,k) = XCRD_FACE_ZTAXSI_TEMP(ii,jj,kk,3)
          BLK%DETX_FACE(i,j,k) = DETX_FACE_TEMP(ii,jj,kk,1)
          BLK%DETY_FACE(i,j,k) = DETX_FACE_TEMP(ii,jj,kk,2)
          BLK%DETZ_FACE(i,j,k) = DETX_FACE_TEMP(ii,jj,kk,3)
          if (i2dimension>0) then
            BLK%DETZ_FACE(i,j,k) = 0.0d0
          endif
        enddo
      enddo
    enddo

    do k=BLK%ksta,BLK%kend+1
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1
          BLK%XCRD_FACE_XSIETA(i,j,k) = XCRD_FACE_XSIETA_TEMP(ii,jj,kk,1)
          BLK%YCRD_FACE_XSIETA(i,j,k) = XCRD_FACE_XSIETA_TEMP(ii,jj,kk,2)
          BLK%ZCRD_FACE_XSIETA(i,j,k) = XCRD_FACE_XSIETA_TEMP(ii,jj,kk,3)
          BLK%DZTX_FACE(i,j,k) = DZTX_FACE_TEMP(ii,jj,kk,1)
          BLK%DZTY_FACE(i,j,k) = DZTX_FACE_TEMP(ii,jj,kk,2)
          BLK%DZTZ_FACE(i,j,k) = DZTX_FACE_TEMP(ii,jj,kk,3)
          if (i2dimension>0) then
            BLK%DZTX_FACE(i,j,k) = 0.0d0
            BLK%DZTY_FACE(i,j,k) = 0.0d0
            BLK%DZTZ_FACE(i,j,k) = 0.0d0
          endif
        enddo
      enddo
    enddo

    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend
          BLK%DXSX(i,j,k) = 0.5d0*(BLK%DXSX_FACE(i+1,j,k)+BLK%DXSX_FACE(i,j,k))
          BLK%DXSY(i,j,k) = 0.5d0*(BLK%DXSY_FACE(i+1,j,k)+BLK%DXSY_FACE(i,j,k))
          BLK%DXSZ(i,j,k) = 0.5d0*(BLK%DXSZ_FACE(i+1,j,k)+BLK%DXSZ_FACE(i,j,k))
          !--------------------------------------------------------------------
          BLK%DETX(i,j,k) = 0.5d0*(BLK%DETX_FACE(i,j+1,k)+BLK%DETX_FACE(i,j,k))
          BLK%DETY(i,j,k) = 0.5d0*(BLK%DETY_FACE(i,j+1,k)+BLK%DETY_FACE(i,j,k))
          BLK%DETZ(i,j,k) = 0.5d0*(BLK%DETZ_FACE(i,j+1,k)+BLK%DETZ_FACE(i,j,k))
          !--------------------------------------------------------------------
          BLK%DZTX(i,j,k) = 0.5d0*(BLK%DZTX_FACE(i,j,k+1)+BLK%DZTX_FACE(i,j,k))
          BLK%DZTY(i,j,k) = 0.5d0*(BLK%DZTY_FACE(i,j,k+1)+BLK%DZTY_FACE(i,j,k))
          BLK%DZTZ(i,j,k) = 0.5d0*(BLK%DZTZ_FACE(i,j,k+1)+BLK%DZTZ_FACE(i,j,k))
        enddo
      enddo
    enddo

    deallocate(XCRD_NODE_TEMP)
    deallocate(XCRD_FACE_ETAZTA_TEMP)
    deallocate(XCRD_FACE_ZTAXSI_TEMP)
    deallocate(XCRD_FACE_XSIETA_TEMP)
    deallocate(DXSX_FACE_TEMP)
    deallocate(DETX_FACE_TEMP)
    deallocate(DZTX_FACE_TEMP)
    deallocate(XCRD_TEMP)
    deallocate(RJCB_TEMP)
    deallocate(GRSC_TEMP)
  end subroutine metric

end module mdl_metric
