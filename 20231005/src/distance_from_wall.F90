module mdl_dfw
#ifdef _USE_MPI
  use mpi
#endif
  use mdl_decompo
  use mdl_mpisub_sbsp
  use mdl_param
  use mdl_refs
  use mdl_block
  implicit none

  integer, private, parameter :: nbound = 7
  integer, private, parameter, dimension(nbound) :: id_dms_wall = (/1,1,2,3,3,4,4/) !DomainID
  integer, private, parameter, dimension(nbound) :: id_bds_wall = (/3,4,4,3,4,1,2/) !FaceID

contains

  subroutine distance_from_wall(BLK)
    type(block), intent(inout) :: BLK

    integer :: myrank
    integer :: id_block,id_domain,id_subdomain,id_dm
    integer :: i,j,k,ii,jj,kk
    integer :: imax_b,jmax_b,kmax_b
    integer :: imax_d,jmax_d,kmax_d

    integer :: ierror

    integer :: npwall,npcount
    real(8), allocatable, dimension(:,:) :: XYZ_WALL

    integer :: m,n
    real(8) :: x,y,z,xw,yw,zw
    real(8) :: dl,dl_min

    real(8), allocatable, dimension(:,:,:,:) :: QB
    real(8), allocatable, dimension(:,:,:,:) :: QD

    call mpisub_sbsp_get_myrank_world(myrank)

    call Get_BlockID_From_Rank(myrank,id_block)

    call Get_DomainID_And_SubDomainID_From_BlockID(id_block,id_domain,id_subdomain)


    !Set QB---------------------------------------------------------------------
    call Get_BlockIJKNum_From_BlockID(id_block,imax_b,jmax_b,kmax_b)
    allocate(QB(imax_b,jmax_b,kmax_b,3))

    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend
          ii = i-BLK%ista+1
          jj = j-BLK%jsta+1
          kk = k-BLK%ksta+1

          QB(ii,jj,kk,1) = BLK%XCRD(i,j,k)
          QB(ii,jj,kk,2) = BLK%YCRD(i,j,k)
          QB(ii,jj,kk,3) = BLK%ZCRD(i,j,k)
        enddo
      enddo
    enddo

    !Gather wall coordinate-----------------------------------------------------
    npwall = 0
    do n=1,nbound
      id_dm = id_dms_wall(n)
      call Get_DomainIJKNum_From_DomainID(id_dm,imax_d,jmax_d,kmax_d)
      if (id_bds_wall(n)==1.or.id_bds_wall(n)==2) then
        npwall = npwall+jmax_d*kmax_d
      else if (id_bds_wall(n)==3.or.id_bds_wall(n)==4) then
        npwall = npwall+imax_d*kmax_d
      else if (id_bds_wall(n)==5.or.id_bds_wall(n)==6) then
        npwall = npwall+imax_d*jmax_d
      endif
    enddo

    allocate(XYZ_WALL(npwall,3))

    npcount = 1
    do n=1,nbound
      id_dm = id_dms_wall(n)
      call Get_DomainIJKNum_From_DomainID(id_dm,imax_d,jmax_d,kmax_d)
      allocate(QD(imax_d,jmax_d,kmax_d,3))

      call mpisub_sbsp_send_domain_to_a_rank(id_dm,QB,0,QD,3,0)

      if (myrank==0) then
        if (id_bds_wall(n)==1.or.id_bds_wall(n)==2) then
          if (id_bds_wall(n)==1) then
            i = 1
          else if (id_bds_wall(n)==2) then
            i = imax_d
          endif
          do k=1,kmax_d
            do j=1,jmax_d
              XYZ_WALL(npcount,1) = QD(i,j,k,1)
              XYZ_WALL(npcount,2) = QD(i,j,k,2)
              XYZ_WALL(npcount,3) = QD(i,j,k,3)
              npcount = npcount+1
            enddo
          enddo
        else if (id_bds_wall(n)==3.or.id_bds_wall(n)==4) then
          if (id_bds_wall(n)==3) then
            j = 1
          else if (id_bds_wall(n)==4) then
            j = jmax_d
          endif
          do k=1,kmax_d
            do i=1,imax_d
              XYZ_WALL(npcount,1) = QD(i,j,k,1)
              XYZ_WALL(npcount,2) = QD(i,j,k,2)
              XYZ_WALL(npcount,3) = QD(i,j,k,3)
              npcount = npcount+1
            enddo
          enddo
        else if (id_bds_wall(n)==5.or.id_bds_wall(n)==6) then
          if (id_bds_wall(n)==5) then
            k = 1
          else if (id_bds_wall(n)==6) then
            k = kmax_d
          endif
          do j=1,jmax_d
            do i=1,imax_d
              XYZ_WALL(npcount,1) = QD(i,j,k,1)
              XYZ_WALL(npcount,2) = QD(i,j,k,2)
              XYZ_WALL(npcount,3) = QD(i,j,k,3)
              npcount = npcount+1
            enddo
          enddo
        endif
      endif

      deallocate(QD)
    enddo

#ifdef _USE_MPI
    call MPI_Bcast(XYZ_WALL,npwall*3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
#endif

    !Calc distance from wall----------------------------------------------------
    do k=BLK%ksta,BLK%kend
      do j=BLK%jsta,BLK%jend
        do i=BLK%ista,BLK%iend
          x = BLK%XCRD(i,j,k)
          y = BLK%YCRD(i,j,k)
          z = BLK%ZCRD(i,j,k)

          dl_min = 1.0d16

          do m=1,npwall
            xw = XYZ_WALL(m,1)
            yw = XYZ_WALL(m,2)
            zw = XYZ_WALL(m,3)
            dl = dsqrt((x-xw)*(x-xw)+(y-yw)*(y-yw)+(z-zw)*(z-zw))
            dl_min = dmin1(dl,dl_min)
          enddo

          BLK%YNKK(i,j,k) = dl_min

        enddo
      enddo
    enddo

    deallocate(XYZ_WALL)
    deallocate(QB)
  end subroutine distance_from_wall

end module mdl_dfw
