module mdl_mpisub_sbsp
#ifdef _USE_MPI
  use mpi
#endif
  use mdl_input_manager
  use mdl_decompo
  implicit none

  !For SX-Aurora

  integer, save, private :: myrank_world

  integer, save, private :: nprocs_world

  integer, save, private :: id_block
  integer, save, private :: id_domain
  integer, save, private :: id_subdomain

  integer, save, private, allocatable, dimension(:) :: BlockID_To_Rank_Domain
  integer, save, private, allocatable, dimension(:) :: BlockID_To_Rank_DomainPlane_I
  integer, save, private, allocatable, dimension(:) :: BlockID_To_Rank_DomainPlane_J
  integer, save, private, allocatable, dimension(:) :: BlockID_To_Rank_DomainPlane_K

  integer, save, private :: MPI_COMM_DOMAIN
  integer, save, private :: myrank_domain

  integer, save, private :: MPI_COMM_DOMAINPLANE_I
  integer, save, private :: myrank_domainplane_i
  !----------------------------------------------
  integer, save, private :: MPI_COMM_DOMAINPLANE_J
  integer, save, private :: myrank_domainplane_j
  !----------------------------------------------
  integer, save, private :: MPI_COMM_DOMAINPLANE_K
  integer, save, private :: myrank_domainplane_k

contains

  subroutine mpisub_sbsp_get_myrank_world(myrank)
    integer, intent(out) :: myrank
    myrank = myrank_world
  end subroutine mpisub_sbsp_get_myrank_world
!-------------------------------------------------------------------------------
  subroutine mpisub_sbsp_get_nprocs_world(nprocs)
    integer, intent(out) :: nprocs
    nprocs = nprocs_world
  end subroutine mpisub_sbsp_get_nprocs_world
!-------------------------------------------------------------------------------
  subroutine mpisub_sbsp_allreduce_min(val)
    real(8), intent(inout) :: val
    integer :: ierror
    real(8) :: tmp
#ifdef _USE_MPI
    if (nprocs_world==1) return
    call MPI_Allreduce(val,tmp,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierror)
    val = tmp
#else
    return
#endif
  end subroutine mpisub_sbsp_allreduce_min
!-------------------------------------------------------------------------------
  subroutine mpisub_sbsp_allreduce_max(val)
    real(8), intent(inout) :: val
    integer :: ierror
    real(8) :: tmp
#ifdef _USE_MPI
    if (nprocs_world==1) return
    call MPI_Allreduce(val,tmp,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierror)
    val = tmp
#else
    return
#endif
  end subroutine mpisub_sbsp_allreduce_max
!-------------------------------------------------------------------------------
  subroutine mpisub_sbsp_init
    integer :: ierror
    integer :: ndomain
    integer :: nblocks,is,js,ks
    integer :: i,id_dm,id_sbdm,id_blk
    integer :: isplt,jsplt,ksplt
    integer :: key_domainplane_i,key_domainplane_j,key_domainplane_k
    integer :: col_domainplane_i,col_domainplane_j,col_domainplane_k
    integer :: key_domain
    integer :: col_domain

#ifdef _USE_MPI
    call MPI_Comm_size(MPI_COMM_WORLD,nprocs_world,ierror)
    call MPI_Comm_rank(MPI_COMM_WORLD,myrank_world,ierror)
#else
    nprocs_world = 1
    myrank_world = 0
#endif

    call Get_DomainNum(ndomain)

    call Get_BlockNum(nblocks)

    call Get_BlockID_From_Rank(myrank_world,id_block)
    call Get_DomainID_And_SubDomainID_From_BlockID(id_block,id_domain,id_subdomain)

    allocate(BlockID_To_Rank_Domain(nblocks))
    allocate(BlockID_To_Rank_DomainPlane_I(nblocks))
    allocate(BlockID_To_Rank_DomainPlane_J(nblocks))
    allocate(BlockID_To_Rank_DomainPlane_K(nblocks))

    do id_blk=1,nblocks
      call Get_DomainID_And_SubDomainID_From_BlockID(id_blk,id_dm,id_sbdm)
      BlockID_To_Rank_Domain(id_blk) = id_sbdm-1
    enddo

    do id_dm=1,ndomain
      call Get_DomainSplit_From_DomainID(id_dm,isplt,jsplt,ksplt)
      do ks=1,ksplt
        do js=1,jsplt
          do is=1,isplt
            call Get_SubDomainID_From_DomainID_And_SplitID(id_dm,is,js,ks,id_sbdm)
            call Get_BlockID_From_DomainID_And_SubDomainID(id_dm,id_sbdm,id_blk)
            BlockID_To_Rank_DomainPlane_I(id_blk) = js+jsplt*(ks-1)-1
            BlockID_To_Rank_DomainPlane_J(id_blk) = is+isplt*(ks-1)-1
            BlockID_To_Rank_DomainPlane_K(id_blk) = is+isplt*(js-1)-1
          enddo
        enddo
      enddo
    enddo

    col_domain = id_domain
    key_domain = BlockID_To_Rank_Domain(id_block)

#ifdef _USE_MPI
    call MPI_Comm_split(MPI_COMM_WORLD,col_domain,key_domain,MPI_COMM_DOMAIN,ierror)
    call MPI_Comm_rank(MPI_COMM_DOMAIN,myrank_domain,ierror)
#else
    MPI_COMM_DOMAIN = 0
    myrank_domain = 0
#endif

    call Get_SplitID_From_DomainID_And_SubDomainID(id_domain,id_subdomain,is,js,ks)

    col_domainplane_i = is
    col_domainplane_j = js
    col_domainplane_k = ks

    key_domainplane_i = BlockID_To_Rank_DomainPlane_I(id_block)
    key_domainplane_j = BlockID_To_Rank_DomainPlane_J(id_block)
    key_domainplane_k = BlockID_To_Rank_DomainPlane_K(id_block)

#ifdef _USE_MPI
    call MPI_Comm_split(MPI_COMM_DOMAIN,col_domainplane_i,key_domainplane_i,MPI_COMM_DOMAINPLANE_I,ierror)
    call MPI_Comm_split(MPI_COMM_DOMAIN,col_domainplane_j,key_domainplane_j,MPI_COMM_DOMAINPLANE_J,ierror)
    call MPI_Comm_split(MPI_COMM_DOMAIN,col_domainplane_k,key_domainplane_k,MPI_COMM_DOMAINPLANE_K,ierror)

    call MPI_Comm_rank(MPI_COMM_DOMAINPLANE_I,myrank_domainplane_i,ierror)
    call MPI_Comm_rank(MPI_COMM_DOMAINPLANE_J,myrank_domainplane_j,ierror)
    call MPI_Comm_rank(MPI_COMM_DOMAINPLANE_K,myrank_domainplane_k,ierror)
#else
    MPI_COMM_DOMAINPLANE_I = 0
    MPI_COMM_DOMAINPLANE_J = 0
    MPI_COMM_DOMAINPLANE_K = 0
    myrank_domainplane_i = 0
    myrank_domainplane_j = 0
    myrank_domainplane_k = 0
#endif

    if (    (myrank_domain       /=BlockID_To_Rank_Domain(id_block))        &
        .or.(myrank_domainplane_i/=BlockID_To_Rank_DomainPlane_I(id_block)) &
        .or.(myrank_domainplane_j/=BlockID_To_Rank_DomainPlane_J(id_block)) &
        .or.(myrank_domainplane_k/=BlockID_To_Rank_DomainPlane_K(id_block)) ) then
        write(6,*) 'Error:mpisub_sbsp_init:Something wrong!!'
    endif

    call mpisub_sbsp_output_info

  end subroutine mpisub_sbsp_init
  !-----------------------------------------------------------------------------
  subroutine mpisub_sbsp_final
    deallocate(BlockID_To_Rank_Domain)
    deallocate(BlockID_To_Rank_DomainPlane_I)
    deallocate(BlockID_To_Rank_DomainPlane_J)
    deallocate(BlockID_To_Rank_DomainPlane_K)
  end subroutine mpisub_sbsp_final
  !-----------------------------------------------------------------------------
  subroutine mpisub_sbsp_readgrid_plot3d(fname,QB,inum_b,jnum_b,knum_b,lnum_b,scale)
    !Read Node Coordinate
    character(*), intent(in) :: fname
    real(8), intent(out), dimension(:,:,:,:) :: QB
    integer, intent(in) :: inum_b,jnum_b,knum_b,lnum_b
    real(8), intent(in) :: scale

    integer :: id_blk,id_dm,ndomain,ndomain_tmp
    integer :: imax,jmax,kmax
    integer :: imax_tmp,jmax_tmp,kmax_tmp
    integer :: sendrecvnum,rank_world_send,rank_world_recv
    integer :: rank_domain_root
    integer :: i,j,k,l,ii,jj,kk
    integer :: istart,iend,jstart,jend,kstart,kend
    integer :: ierror
    real(8), allocatable, dimension(:,:,:,:) :: XD
#ifdef _USE_MPI
    integer, dimension(MPI_STATUS_SIZE) :: istatus
#endif

    101 format(1i9)
    103 format(3i9)
    404 format(1e22.15,3e23.15)

    call Get_DomainNum(ndomain)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)
    if (inum_b/=imax+1.or.jnum_b/=jmax+1.or.knum_b/=kmax+1.or.lnum_b/=3) then
      write(6,*) 'Error:mpisub_sbsp_readgrid_plot3d:Invalid arguments 1'
    endif

    if (myrank_world==0) then
      open(12,file=fname,form='formatted')
    endif

    if (myrank_world==0) then

      read(12,101) ndomain_tmp

      if (ndomain_tmp/=ndomain) then
        write(6,*) 'Error:mpisub_sbsp_readgrid_plot3d:Invalid arguments 2'
      endif

      do id_dm=1,ndomain
        read(12,103) imax_tmp,jmax_tmp,kmax_tmp
        call Get_DomainIJKNum_From_DomainID(id_dm,imax,jmax,kmax)
        if ((imax_tmp/=imax+1).or.(jmax_tmp/=jmax+1).or.(kmax_tmp/=kmax+1)) then
          write(6,*) 'Error:mpisub_sbsp_readgrid_plot3d:Invalid arguments 3'
        endif
      enddo

    endif


    do id_dm=1,ndomain
      call Get_DomainIJKNum_From_DomainID(id_dm,imax,jmax,kmax)

      sendrecvnum = (imax+1)*(jmax+1)*(kmax+1)*3

      if (myrank_world==0.or.id_domain==id_dm) then
        allocate(XD(imax+1,jmax+1,kmax+1,3))
      else
        allocate(XD(1     ,1     ,1     ,1))
      endif

      if (myrank_world==0) then
        read(12,404) (((XD(i,j,k,1),i=1,imax+1),j=1,jmax+1),k=1,kmax+1)
        read(12,404) (((XD(i,j,k,2),i=1,imax+1),j=1,jmax+1),k=1,kmax+1)
        read(12,404) (((XD(i,j,k,3),i=1,imax+1),j=1,jmax+1),k=1,kmax+1)
        XD = XD*scale
      endif

      rank_world_send = 0

      call Get_BlockID_From_DomainID_And_SubDomainID(id_dm,1,id_blk)
      call Get_Rank_From_BlockID(id_blk,rank_world_recv)
      rank_domain_root = BlockID_To_Rank_Domain(id_blk)

      if (rank_world_send/=rank_world_recv) then
        if (myrank_world==rank_world_send) then
#ifdef _USE_MPI
          call MPI_Send(XD,sendrecvnum,MPI_DOUBLE_PRECISION,rank_world_recv,0, &
                        MPI_COMM_WORLD,ierror)
#endif
        endif
        if (myrank_world==rank_world_recv) then
#ifdef _USE_MPI
          call MPI_Recv(XD,sendrecvnum,MPI_DOUBLE_PRECISION,rank_world_send,0, &
                        MPI_COMM_WORLD,istatus,ierror)
#endif
        endif
      endif

      if (id_dm==id_domain) then
#ifdef _USE_MPI
        call MPI_Bcast(XD,sendrecvnum,MPI_DOUBLE_PRECISION,rank_domain_root,MPI_COMM_DOMAIN,ierror)
#endif

        call Get_SubDomainIJKStartEnd_From_DomainID_And_SubDomainID(id_domain,id_subdomain, &
                                                        istart,iend,jstart,jend,kstart,kend)

        do l=1,3
          do k=kstart,kend+1
            do j=jstart,jend+1
              do i=istart,iend+1
                ii = i-istart+1
                jj = j-jstart+1
                kk = k-kstart+1
                QB(ii,jj,kk,l) = XD(i,j,k,l)
              enddo
            enddo
          enddo
        enddo

      endif

      deallocate(XD)
    enddo

    if (myrank_world==0) then
      close(12)
    endif
  end subroutine mpisub_sbsp_readgrid_plot3d
  !-----------------------------------------------------------------------------
  subroutine mpisub_sbsp_out_plot3d(fname,mode,QB,inum_b,jnum_b,knum_b,lnum_b)
    character(*), intent(in) :: fname
    integer, intent(in) :: mode !<0>:grid <1>:functoins
    real(8), intent(in), dimension(:,:,:,:) :: QB
    integer, intent(in) :: inum_b,jnum_b,knum_b,lnum_b

    integer :: i,j,k,l,ierror
    integer :: imax,jmax,kmax,lmax
    integer :: id_dm,ndomain
    real(8), allocatable, dimension(:,:,:,:) :: QD

    call Get_DomainNum(ndomain)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)
    if (inum_b/=imax.or.jnum_b/=jmax.or.knum_b/=kmax) then
      write(6,*) 'Error:mpisub_sbsp_out_plot3d:Invalid arguments 1'
    endif

    if (myrank_world==0) then
      open(21,file=fname,status='replace',form='unformatted',access='stream')
    endif

    if (myrank_world==0) then

      write(21) ndomain

      do id_dm=1,ndomain
        call Get_DomainIJKNum_From_DomainID(id_dm,imax,jmax,kmax)
        lmax = lnum_b
        write(21) imax
        write(21) jmax
        write(21) kmax
        if (mode==1) then
          write(21) lmax
        endif
      enddo

    endif

    do id_dm=1,ndomain
      call Get_DomainIJKNum_From_DomainID(id_dm,imax,jmax,kmax)
      lmax = lnum_b

      if (myrank_world==0) then
        allocate(QD(imax,jmax,kmax,lmax))
      else
        allocate(QD(1   ,1   ,1   ,1   ))
      endif

      call mpisub_sbsp_send_domain_to_a_rank(id_dm,QB,0,QD,lmax,0)

      if (myrank_world==0) then
        write(21) real(QD(1:imax,1:jmax,1:kmax,1:lmax))
      endif

      deallocate(QD)
    enddo

    if (myrank_world==0) then
      close(21)
    endif

#ifdef _USE_MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
#endif
  end subroutine mpisub_sbsp_out_plot3d
  !-----------------------------------------------------------------------------
  subroutine mpisub_sbsp_out_plot3d_planes(fname,mode,QB,inum_b,jnum_b,knum_b,lnum_b, &
                                           id_dms,id_dirs,idxs,nplanes)
    character(*), intent(in) :: fname
    integer, intent(in) :: mode !<0>:grid <1>:functions
    real(8), intent(in), dimension(:,:,:,:) :: QB
    integer, intent(in) :: inum_b,jnum_b,knum_b,lnum_b
    integer, intent(in) :: nplanes
    integer, intent(inout), dimension(nplanes) :: id_dms
    integer, intent(inout), dimension(nplanes) :: id_dirs
    integer, intent(inout), dimension(nplanes) :: idxs

    integer :: id_dm,ndomain
    integer :: imax,jmax,kmax,lmax
    integer :: id_pln,itmp
    integer :: i,j,k,l,ierror
    real(8), allocatable, dimension(:,:,:,:) :: QD

    call Get_DomainNum(ndomain)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)
    if (inum_b/=imax.or.jnum_b/=jmax.or.knum_b/=kmax) then
      write(6,*) 'Error::mpisub_sbsp_out_plot3d_planes::Invalid arguments 1'
    endif

    do id_pln=1,nplanes
      if (id_dms(id_pln)<1.or.ndomain<id_dms(id_pln)) then
        write(6,*) 'error::mpisub_sbsp_out_plot3d_planes::invalid arguments 2'
      endif
    enddo

    do id_pln=1,nplanes
      if (.not.(id_dirs(id_pln)==1.or.id_dirs(id_pln)==2.or.id_dirs(id_pln)==3)) then
        write(6,*) 'error::mpisub_sbsp_out_plot3d_planes::invalid arguments 3'
      endif
    enddo

    do id_pln=1,nplanes
      id_dm = id_dms(id_pln)
      call Get_DomainIJKNum_From_DomainID(id_dm,imax,jmax,kmax)
      if (id_dirs(id_pln)==1) then
        if (idxs(id_pln)<1.or.imax<idxs(id_pln)) then
          write(6,*) 'error::mpisub_sbsp_out_plot3d_planes::invalid arguments 4'
        endif
      else if (id_dirs(id_pln)==2) then
        if (idxs(id_pln)<1.or.jmax<idxs(id_pln)) then
          write(6,*) 'error::mpisub_sbsp_out_plot3d_planes::invalid arguments 5'
        endif
      else if (id_dirs(id_pln)==3) then
        if (idxs(id_pln)<1.or.kmax<idxs(id_pln)) then
          write(6,*) 'error::mpisub_sbsp_out_plot3d_planes::invalid arguments 6'
        endif
      endif
    enddo

    !Sort----------------------
    do j=1,nplanes
      do i=1,nplanes-j
        if (id_dms(i)>id_dms(i+1)) then
          itmp = id_dms(i)
          id_dms(i) = id_dms(i+1)
          id_dms(i+1) = itmp
          !----------------------
          itmp = id_dirs(i)
          id_dirs(i) = id_dirs(i+1)
          id_dirs(i+1) = itmp
          !----------------------
          itmp = idxs(i)
          idxs(i) = idxs(i+1)
          idxs(i+1) = itmp
        else if (id_dms(i)==id_dms(i+1)) then
          if (id_dirs(i)>id_dirs(i+1)) then
            itmp = id_dirs(i)
            id_dirs(i) = id_dirs(i+1)
            id_dirs(i+1) = itmp
            !----------------------
            itmp = idxs(i)
            idxs(i) = idxs(i+1)
            idxs(i+1) = itmp
          else if (id_dirs(i)==id_dirs(i+1)) then
            if (idxs(i)>idxs(i+1)) then
              itmp = idxs(i)
              idxs(i) = idxs(i+1)
              idxs(i+1) = itmp
            endif
          endif
        endif
      enddo
    enddo
    !--------------------------

    if (myrank_world==0) then
      open(21,file=fname,status='replace',form='unformatted',access='stream')
    endif

    if (myrank_world==0) then

      write(21) nplanes

      do id_pln=1,nplanes
        id_dm = id_dms(id_pln)
        call Get_DomainIJKNum_From_DomainID(id_dm,imax,jmax,kmax)
        lmax = lnum_b
        if (id_dirs(id_pln)==1) then
          write(21) 1
          write(21) jmax
          write(21) kmax
          if (mode==1) then
            write(21) lmax
          endif
        else if (id_dirs(id_pln)==2) then
          write(21) imax
          write(21) 1
          write(21) kmax
          if (mode==1) then
            write(21) lmax
          endif
        else if (id_dirs(id_pln)==3) then
          write(21) imax
          write(21) jmax
          write(21) 1
          if (mode==1) then
            write(21) lmax
          endif
        endif
      enddo

    endif

    id_pln = 1

    do id_dm=1,ndomain
      call Get_DomainIJKNum_From_DomainID(id_dm,imax,jmax,kmax)
      lmax = lnum_b

      if (myrank_world==0) then
        allocate(QD(imax,jmax,kmax,lmax))
      else
        allocate(QD(1   ,1   ,1   ,1   ))
      endif

      call mpisub_sbsp_send_domain_to_a_rank(id_dm,QB,0,QD,lmax,0)

      if (myrank_world==0) then

        do while (id_pln<=nplanes.and.id_dms(id_pln)==id_dm)
          if (id_dirs(id_pln)==1) then
            i = idxs(id_pln)
            write(21) real(QD(i,1:jmax,1:kmax,1:lmax))
          else if (id_dirs(id_pln)==2) then
            j = idxs(id_pln)
            write(21) real(QD(1:imax,j,1:kmax,1:lmax))
          else if (id_dirs(id_pln)==3) then
            k = idxs(id_pln)
            write(21) real(QD(1:imax,1:jmax,k,1:lmax))
          endif
          id_pln = id_pln+1
        enddo

      endif

      deallocate(QD)
    enddo

    if (myrank_world==0) then
      close(21)
    endif

#ifdef _USE_MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
#endif
  end subroutine mpisub_sbsp_out_plot3d_planes
!-------------------------------------------------------------------------------
  subroutine mpisub_sbsp_out_plot3d_block_with_halo(fname,mode,halo_width,iflag_halo, &
             QB_WH,inum_b_wh,jnum_b_wh,knum_b_wh,lnum_b)
    character(*), intent(in) :: fname
    integer, intent(in) :: mode !<0>:grid <1>:functoins
    integer, intent(in) :: halo_width
    integer, intent(in), dimension(6) :: iflag_halo
    real(8), intent(in), dimension(:,:,:,:) :: QB_WH
    integer, intent(in) :: inum_b_wh,jnum_b_wh,knum_b_wh,lnum_b

    integer :: i,j,k,ii,jj,kk,l,ierror
    integer :: nblocks,id_blk,irank
    integer :: nbo
    integer :: inum_b,jnum_b,knum_b

    real(8), allocatable, dimension(:,:,:,:) :: QB_TMP

    inum_b = inum_b_wh-2*halo_width
    jnum_b = jnum_b_wh-2*halo_width
    knum_b = knum_b_wh-2*halo_width

    call Get_BlockNum(nblocks)

    do id_blk=1,nblocks
      call Get_Rank_From_BlockID(id_blk,irank)
      if (myrank_world==irank) then

        nbo = 1
        do i=1,6
          if (iflag_halo(i)==1) nbo = nbo+1
        enddo

        open(21,file=fname,status='replace',form='unformatted',access='stream')

        write(21) nbo

        write(21) inum_b
        write(21) jnum_b
        write(21) knum_b
        if (mode==1) then
          write(21) lnum_b
        endif

        if (iflag_halo(1)==1) then
          write(21) halo_width
          write(21) jnum_b
          write(21) knum_b
          if (mode==1) then
            write(21) lnum_b
          endif
        endif

        if (iflag_halo(2)==1) then
          write(21) halo_width
          write(21) jnum_b
          write(21) knum_b
          if (mode==1) then
            write(21) lnum_b
          endif
        endif

        if (iflag_halo(3)==1) then
          write(21) inum_b
          write(21) halo_width
          write(21) knum_b
          if (mode==1) then
            write(21) lnum_b
          endif
        endif

        if (iflag_halo(4)==1) then
          write(21) inum_b
          write(21) halo_width
          write(21) knum_b
          if (mode==1) then
            write(21) lnum_b
          endif
        endif

        if (iflag_halo(5)==1) then
          write(21) inum_b
          write(21) jnum_b
          write(21) halo_width
          if (mode==1) then
            write(21) lnum_b
          endif
        endif

        if (iflag_halo(6)==1) then
          write(21) inum_b
          write(21) jnum_b
          write(21) halo_width
          if (mode==1) then
            write(21) lnum_b
          endif
        endif

        allocate(QB_TMP(inum_b,jnum_b,knum_b,lnum_b))

        do l=1,lnum_b
          do k=halo_width+1,halo_width+knum_b
            do j=halo_width+1,halo_width+jnum_b
              do i=halo_width+1,halo_width+inum_b
                ii = i-(halo_width+1)+1
                jj = j-(halo_width+1)+1
                kk = k-(halo_width+1)+1
                QB_TMP(ii,jj,kk,l) = QB_WH(i,j,k,l)
              enddo
            enddo
          enddo
        enddo

        write(21) real(QB_TMP(1:inum_b,1:jnum_b,1:knum_b,1:lnum_b))

        deallocate(QB_TMP)

        if (iflag_halo(1)==1) then

          allocate(QB_TMP(halo_width,jnum_b,knum_b,lnum_b))

          do l=1,lnum_b
            do k=halo_width+1,halo_width+knum_b
              do j=halo_width+1,halo_width+jnum_b
                do i=1,halo_width
                  ii = i
                  jj = j-(halo_width+1)+1
                  kk = k-(halo_width+1)+1
                  QB_TMP(ii,jj,kk,l) = QB_WH(i,j,k,l)
                enddo
              enddo
            enddo
          enddo

          write(21) real(QB_TMP(1:halo_width,1:jnum_b,1:knum_b,1:lnum_b))

          deallocate(QB_TMP)

        endif

        if (iflag_halo(2)==1) then

          allocate(QB_TMP(halo_width,jnum_b,knum_b,lnum_b))

          do l=1,lnum_b
            do k=halo_width+1,halo_width+knum_b
              do j=halo_width+1,halo_width+jnum_b
                do i=halo_width+inum_b+1,inum_b_wh
                  ii = i-(halo_width+inum_b+1)+1
                  jj = j-(halo_width+1)+1
                  kk = k-(halo_width+1)+1
                  QB_TMP(ii,jj,kk,l) = QB_WH(i,j,k,l)
                enddo
              enddo
            enddo
          enddo

          write(21) real(QB_TMP(1:halo_width,1:jnum_b,1:knum_b,1:lnum_b))

          deallocate(QB_TMP)

        endif

        if (iflag_halo(3)==1) then

          allocate(QB_TMP(inum_b,halo_width,knum_b,lnum_b))

          do l=1,lnum_b
            do k=halo_width+1,halo_width+knum_b
              do j=1,halo_width
                do i=halo_width+1,halo_width+inum_b
                  ii = i-(halo_width+1)+1
                  jj = j
                  kk = k-(halo_width+1)+1
                  QB_TMP(ii,jj,kk,l) = QB_WH(i,j,k,l)
                enddo
              enddo
            enddo
          enddo

          write(21) real(QB_TMP(1:inum_b,1:halo_width,1:knum_b,1:lnum_b))

          deallocate(QB_TMP)

        endif

        if (iflag_halo(4)==1) then

          allocate(QB_TMP(inum_b,halo_width,knum_b,lnum_b))

          do l=1,lnum_b
            do k=halo_width+1,halo_width+knum_b
              do j=halo_width+jnum_b+1,jnum_b_wh
                do i=halo_width+1,halo_width+inum_b
                  ii = i-(halo_width+1)+1
                  jj = j-(halo_width+jnum_b+1)+1
                  kk = k-(halo_width+1)+1
                  QB_TMP(ii,jj,kk,l) = QB_WH(i,j,k,l)
                enddo
              enddo
            enddo
          enddo

          write(21) real(QB_TMP(1:inum_b,1:halo_width,1:knum_b,1:lnum_b))

          deallocate(QB_TMP)

        endif

        if (iflag_halo(5)==1) then

          allocate(QB_TMP(inum_b,jnum_b,halo_width,lnum_b))

          do l=1,lnum_b
            do k=1,halo_width
              do j=halo_width+1,halo_width+jnum_b
                do i=halo_width+1,halo_width+inum_b
                  ii = i-(halo_width+1)+1
                  jj = j-(halo_width+1)+1
                  kk = k
                  QB_TMP(ii,jj,kk,l) = QB_WH(i,j,k,l)
                enddo
              enddo
            enddo
          enddo

          write(21) real(QB_TMP(1:inum_b,1:jnum_b,1:halo_width,1:lnum_b))

          deallocate(QB_TMP)

        endif

        if (iflag_halo(6)==1) then

          allocate(QB_TMP(inum_b,jnum_b,halo_width,lnum_b))

          do l=1,lnum_b
            do k=halo_width+knum_b+1,knum_b_wh
              do j=halo_width+1,halo_width+jnum_b
                do i=halo_width+1,halo_width+inum_b
                  ii = i-(halo_width+1)+1
                  jj = j-(halo_width+1)+1
                  kk = k-(halo_width+knum_b+1)+1
                  QB_TMP(ii,jj,kk,l) = QB_WH(i,j,k,l)
                enddo
              enddo
            enddo
          enddo

          write(21) real(QB_TMP(1:inum_b,1:jnum_b,1:halo_width,1:lnum_b))

          deallocate(QB_TMP)

        endif

        close(21)

        call flush(21)

      endif
#ifdef _USE_MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
#endif
    enddo

end subroutine mpisub_sbsp_out_plot3d_block_with_halo
  !-----------------------------------------------------------------------------
  subroutine mpisub_sbsp_out_restart(prefix,suffix,nstep, &
                                     QB,inum_b,jnum_b,knum_b,lnum_b)
    character(*), intent(in) :: prefix
    character(*), intent(in) :: suffix
    integer, intent(in) :: nstep
    real(8), intent(in), dimension(:,:,:,:) :: QB
    integer, intent(in) :: inum_b,jnum_b,knum_b,lnum_b

    integer :: ndomain,id_dm,id_blk
    integer :: rank_domain_root
    integer :: imax,jmax,kmax,lmax
    integer :: ierror
    real(8), allocatable, dimension(:,:,:,:) :: QD

    character(20 ) :: affix
    character(100) :: fname

 800 format(i6.6,'_',i4.4)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)
    if (inum_b/=imax.or.jnum_b/=jmax.or.knum_b/=kmax) then
      write(6,*) 'Error::mpisub_sbsp_out_restart::Invalid arguments 1'
    endif

    call Get_DomainNum(ndomain)

    do id_dm=1,ndomain
      if (id_domain==id_dm) then
        call Get_BlockID_From_DomainID_And_SubDomainID(id_dm,1,id_blk)
        rank_domain_root = BlockID_To_Rank_Domain(id_blk)

        call Get_DomainIJKNum_From_DomainID(id_dm,imax,jmax,kmax)
        lmax = lnum_b

        if (myrank_domain==rank_domain_root) then
          allocate(QD(imax,jmax,kmax,lmax))
        else
          allocate(QD(1   ,1   ,1   ,1   ))
        endif

        call mpisub_sbsp_gatherv_from_domain(id_dm,QB,rank_domain_root,QD,lmax,0)

        write(affix,800) nstep,id_dm

        fname = trim(prefix)//trim(affix)//trim(suffix)

        if (myrank_domain==rank_domain_root) then
          open(19,file=fname,form='unformatted')
          write(19) imax,jmax,kmax,lmax
          write(19) QD
          close(19)
        endif

        deallocate(QD)
      endif
#ifdef _USE_MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
#endif
    enddo

  end subroutine mpisub_sbsp_out_restart
  !-----------------------------------------------------------------------------
  subroutine mpisub_sbsp_read_restart(prefix,suffix,nstep, &
                                      QB,inum_b,jnum_b,knum_b,lnum_b)
    character(*), intent(in) :: prefix
    character(*), intent(in) :: suffix
    integer, intent(in) :: nstep
    real(8), intent(out), dimension(:,:,:,:) :: QB
    integer, intent(in) :: inum_b,jnum_b,knum_b,lnum_b

    integer :: ndomain,id_dm,id_blk
    integer :: rank_domain_root
    integer :: imax,jmax,kmax,lmax
    integer :: imax_tmp,jmax_tmp,kmax_tmp,lmax_tmp
    integer :: sendrecvnum
    integer :: ierror
    integer :: i,j,k,l,ii,jj,kk
    integer :: istart,iend,jstart,jend,kstart,kend
    real(8), allocatable, dimension(:,:,:,:) :: QD

    character(20 ) :: affix
    character(100) :: fname

 800 format(i6.6,'_',i4.4)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)
    if (inum_b/=imax.or.jnum_b/=jmax.or.knum_b/=kmax) then
      write(6,*) 'Error::mpisub_sbsp_read_restart::Invalid arguments 1'
    endif

    call Get_DomainNum(ndomain)

    do id_dm=1,ndomain
      if (id_domain==id_dm) then
        call Get_BlockID_From_DomainID_And_SubDomainID(id_dm,1,id_blk)
        rank_domain_root = BlockID_To_Rank_Domain(id_blk)

        call Get_DomainIJKNum_From_DomainID(id_dm,imax,jmax,kmax)
        lmax = lnum_b

        allocate(QD(imax,jmax,kmax,lmax))

        write(affix,800) nstep,id_dm

        fname = trim(prefix)//trim(affix)//trim(suffix)

        if (myrank_domain==rank_domain_root) then
          open(19,file=fname,form='unformatted')
          read(19) imax_tmp,jmax_tmp,kmax_tmp,lmax_tmp
          if (imax_tmp/=imax.or.jmax_tmp/=jmax.or.kmax_tmp/=kmax.or.lmax_tmp/=lmax) then
            write(6,*) 'Error::mpisub_sbsp_read_restart::Invalid arguments 2'
          endif
          read(19) QD
          close(19)
        endif

        sendrecvnum = imax*jmax*kmax*lmax

#ifdef _USE_MPI
        call MPI_Bcast(QD,sendrecvnum,MPI_DOUBLE_PRECISION,rank_domain_root, &
                                                        MPI_COMM_DOMAIN,ierror)
#endif

        call Get_SubDomainIJKStartEnd_From_DomainID_And_SubDomainID( &
                    id_domain,id_subdomain,istart,iend,jstart,jend,kstart,kend)

        do l=1,lmax
          do k=kstart,kend
            do j=jstart,jend
              do i=istart,iend
                ii = i-istart+1
                jj = j-jstart+1
                kk = k-kstart+1
                QB(ii,jj,kk,l) = QD(i,j,k,l)
              enddo
            enddo
          enddo
        enddo

        deallocate(QD)
      endif
#ifdef _USE_MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
#endif
    enddo

  end subroutine mpisub_sbsp_read_restart
  !-----------------------------------------------------------------------------
  subroutine mpisub_sbsp_read_restart_from_2D_to_3D(prefix,suffix,nstep, &
                                                QB,inum_b,jnum_b,knum_b,lnum_b)
    character(*), intent(in) :: prefix
    character(*), intent(in) :: suffix
    integer, intent(in) :: nstep
    real(8), intent(out), dimension(:,:,:,:) :: QB
    integer, intent(in) :: inum_b,jnum_b,knum_b,lnum_b

    integer :: ndomain,id_dm,id_blk
    integer :: rank_domain_root
    integer :: imax,jmax,kmax,lmax
    integer :: imax_tmp,jmax_tmp,kmax_tmp,lmax_tmp
    integer :: sendrecvnum
    integer :: ierror
    integer :: i,j,k,l,ii,jj,kk
    integer :: istart,iend,jstart,jend,kstart,kend
    real(8), allocatable, dimension(:,:,:,:) :: QD
    real(8), allocatable, dimension(:,:,:,:) :: QD_TMP

    character(20 ) :: affix
    character(100) :: fname

 800 format(i6.6,'_',i4.4)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)
    if (inum_b/=imax.or.jnum_b/=jmax.or.knum_b/=kmax) then
      write(6,*) 'Error::mpisub_sbsp_read_restart_from_2D_to_3D::Invalid arguments 1'
    endif

    call Get_DomainNum(ndomain)

    do id_dm=1,ndomain
      if (id_domain==id_dm) then
        call Get_BlockID_From_DomainID_And_SubDomainID(id_dm,1,id_blk)
        rank_domain_root = BlockID_To_Rank_Domain(id_blk)

        call Get_DomainIJKNum_From_DomainID(id_dm,imax,jmax,kmax)
        lmax = lnum_b

        allocate(QD(imax,jmax,kmax,lmax))

        write(affix,800) nstep,id_dm

        fname = trim(prefix)//trim(affix)//trim(suffix)

        if (myrank_domain==rank_domain_root) then
          allocate(QD_TMP(imax,jmax,1,lmax))
          open(19,file=fname,form='unformatted')
          read(19) imax_tmp,jmax_tmp,kmax_tmp,lmax_tmp
          if (imax_tmp/=imax.or.jmax_tmp/=jmax.or.kmax_tmp/=1.or.lmax_tmp/=lmax) then
            write(6,*) 'Error::mpisub_sbsp_read_restart_from_2D_to_3D::Invalid arguments 2'
          endif
          read(19) QD_TMP
          close(19)
          do l=1,lmax
            do k=1,kmax
              do j=1,jmax
                do i=1,imax
                  QD(i,j,k,l) = QD_TMP(i,j,1,l)
                enddo
              enddo
            enddo
          enddo
          deallocate(QD_TMP)
        endif

        sendrecvnum = imax*jmax*kmax*lmax

#ifdef _USE_MPI
        call MPI_Bcast(QD,sendrecvnum,MPI_DOUBLE_PRECISION,rank_domain_root, &
                                                        MPI_COMM_DOMAIN,ierror)
#endif

        call Get_SubDomainIJKStartEnd_From_DomainID_And_SubDomainID( &
                    id_domain,id_subdomain,istart,iend,jstart,jend,kstart,kend)

        do l=1,lmax
          do k=kstart,kend
            do j=jstart,jend
              do i=istart,iend
                ii = i-istart+1
                jj = j-jstart+1
                kk = k-kstart+1
                QB(ii,jj,kk,l) = QD(i,j,k,l)
              enddo
            enddo
          enddo
        enddo

        deallocate(QD)
      endif
#ifdef _USE_MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
#endif
    enddo

  end subroutine mpisub_sbsp_read_restart_from_2D_to_3D
  !-----------------------------------------------------------------------------
  subroutine mpisub_sbsp_calc_l2norm(QB,inum_b,jnum_b,knum_b,lnum_b,l2norm)
    real(8), intent(in), dimension(:,:,:,:) :: QB
    integer, intent(in) :: inum_b,jnum_b,knum_b,lnum_b
    real(8), intent(inout), dimension(:) :: l2norm

    integer :: i,j,k,l
    integer :: nblocks,id_blk
    integer :: npsum
    integer :: imax,jmax,kmax,lmax
    real(8), allocatable, dimension(:) :: l2norm_tmp
    integer :: ierror

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)
    if (inum_b/=imax.or.jnum_b/=jmax.or.knum_b/=kmax) then
        write(6,*) 'Error::mpisub_sbsp_calc_l2norm::Invalid arguments 1'
    endif

    call Get_BlockNum(nblocks)

    npsum = 0
    do id_blk=1,nblocks
      call Get_BlockIJKNum_From_BlockID(id_blk,imax,jmax,kmax)
      npsum = npsum+imax*jmax*kmax
    enddo

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)
    lmax = lnum_b

    allocate(l2norm_tmp(lmax))
    l2norm_tmp = 0.0d0

    do l=1,lmax
      do k=1,kmax
        do j=1,jmax
          do i=1,imax
            l2norm_tmp(l) = l2norm_tmp(l)+QB(i,j,k,l)*QB(i,j,k,l)/dble(npsum)
          enddo
        enddo
      enddo
    enddo

    l2norm = 0.0d0

#ifdef _USE_MPI
    call MPI_Allreduce(l2norm_tmp,l2norm,lmax,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,MPI_COMM_WORLD,ierror)
#endif

    do l=1,lmax
      l2norm(l) = dsqrt(l2norm(l))
    enddo

    deallocate(l2norm_tmp)
  end subroutine mpisub_sbsp_calc_l2norm
  !-----------------------------------------------------------------------------
  subroutine mpisub_sbsp_out_plot3d_mean_along_zta(fname,QB,inum_b,jnum_b,knum_b,lnum_b)
    character(*), intent(in) :: fname
    real(8), intent(in), dimension(:,:,:,:) :: QB
    integer, intent(in) :: inum_b,jnum_b,knum_b,lnum_b

    integer :: i,j,k,l,ierror
    integer :: imax,jmax,kmax,lmax
    integer :: id_dm,ndomain
    real(8), allocatable, dimension(:,:,:,:) :: QD
    real(8), allocatable, dimension(:,:,:)   :: QD_2D

    call Get_DomainNum(ndomain)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)
    if (inum_b/=imax.or.jnum_b/=jmax.or.knum_b/=kmax) then
      write(6,*) 'Error:mpisub_sbsp_out_plot3d:Invalid arguments 1'
    endif

    if (myrank_world==0) then
      open(21,file=fname,status='replace',form='unformatted',access='stream')
    endif

    if (myrank_world==0) then

      write(21) ndomain

      do id_dm=1,ndomain
        call Get_DomainIJKNum_From_DomainID(id_dm,imax,jmax,kmax)
        lmax = lnum_b
        write(21) imax
        write(21) jmax
        write(21) 1
        write(21) lmax
      enddo

    endif

    do id_dm=1,ndomain
      call Get_DomainIJKNum_From_DomainID(id_dm,imax,jmax,kmax)
      lmax = lnum_b

      if (myrank_world==0) then
        allocate(QD(imax,jmax,kmax,lmax))
      else
        allocate(QD(1   ,1   ,1   ,1   ))
      endif

      if (myrank_world==0) then
        allocate(QD_2D(imax,jmax,lmax))
      else
        allocate(QD_2D(1   ,1   ,1   ))
      endif

      call mpisub_sbsp_send_domain_to_a_rank(id_dm,QB,0,QD,lmax,0)

      if (myrank_world==0) then

        QD_2D = 0.0d0

        do l=1,lmax
          do k=1,kmax
            do j=1,jmax
              do i=1,imax
                QD_2D(i,j,l) = QD_2D(i,j,l) + QD(i,j,k,l)
              enddo
            enddo
          enddo
        enddo

        do l=1,lmax
          do j=1,jmax
            do i=1,imax
              QD_2D(i,j,l) = QD_2D(i,j,l)/dble(kmax)
            enddo
          enddo
        enddo

        write(21) real(QD_2D(1:imax,1:jmax,1:lmax))
      endif

      deallocate(QD)
      deallocate(QD_2D)
    enddo

    if (myrank_world==0) then
      close(21)
    endif

#ifdef _USE_MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
#endif
  end subroutine mpisub_sbsp_out_plot3d_mean_along_zta
  !-----------------------------------------------------------------------------
  subroutine mpisub_sbsp_out_restart_mean_along_zta(prefix,suffix,nstep, &
                                     QB,inum_b,jnum_b,knum_b,lnum_b)
    character(*), intent(in) :: prefix
    character(*), intent(in) :: suffix
    integer, intent(in) :: nstep
    real(8), intent(in), dimension(:,:,:,:) :: QB
    integer, intent(in) :: inum_b,jnum_b,knum_b,lnum_b

    integer :: i,j,k,l
    integer :: ndomain,id_dm,id_blk
    integer :: rank_domain_root
    integer :: imax,jmax,kmax,lmax
    integer :: ierror
    real(8), allocatable, dimension(:,:,:,:) :: QD
    real(8), allocatable, dimension(:,:,:)   :: QD_2D

    character(20 ) :: affix
    character(100) :: fname

 800 format(i6.6,'_',i4.4)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)
    if (inum_b/=imax.or.jnum_b/=jmax.or.knum_b/=kmax) then
      write(6,*) 'Error::mpisub_sbsp_out_restart::Invalid arguments 1'
    endif

    call Get_DomainNum(ndomain)

    do id_dm=1,ndomain
      if (id_domain==id_dm) then
        call Get_BlockID_From_DomainID_And_SubDomainID(id_dm,1,id_blk)
        rank_domain_root = BlockID_To_Rank_Domain(id_blk)

        call Get_DomainIJKNum_From_DomainID(id_dm,imax,jmax,kmax)
        lmax = lnum_b

        if (myrank_domain==rank_domain_root) then
          allocate(QD(imax,jmax,kmax,lmax))
        else
          allocate(QD(1   ,1   ,1   ,1   ))
        endif

        if (myrank_domain==rank_domain_root) then
          allocate(QD_2D(imax,jmax,lmax))
        else
          allocate(QD_2D(1   ,1   ,1   ))
        endif

        call mpisub_sbsp_gatherv_from_domain(id_dm,QB,rank_domain_root,QD,lmax,0)

        write(affix,800) nstep,id_dm

        fname = trim(prefix)//trim(affix)//trim(suffix)

        if (myrank_domain==rank_domain_root) then

          QD_2D = 0.0d0

          do l=1,lmax
            do k=1,kmax
              do j=1,jmax
                do i=1,imax
                  QD_2D(i,j,l) = QD_2D(i,j,l) + QD(i,j,k,l)
                enddo
              enddo
            enddo
          enddo

          do l=1,lmax
            do j=1,jmax
              do i=1,imax
                QD_2D(i,j,l) = QD_2D(i,j,l)/dble(kmax)
              enddo
            enddo
          enddo

          open(19,file=fname,form='unformatted')
          write(19) imax,jmax,1,lmax
          write(19) QD_2D
          close(19)
        endif

        deallocate(QD)
        deallocate(QD_2D)
      endif
#ifdef _USE_MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
#endif
    enddo

  end subroutine mpisub_sbsp_out_restart_mean_along_zta
  !-----------------------------------------------------------------------------
  subroutine mpisub_sbsp_out_restart_mean_along_xsi_zta(prefix,suffix,nstep, &
                                     QB,inum_b,jnum_b,knum_b,lnum_b)
    character(*), intent(in) :: prefix
    character(*), intent(in) :: suffix
    integer, intent(in) :: nstep
    real(8), intent(in), dimension(:,:,:,:) :: QB
    integer, intent(in) :: inum_b,jnum_b,knum_b,lnum_b

    integer :: i,j,k,l
    integer :: ndomain,id_dm,id_blk
    integer :: rank_domain_root
    integer :: imax,jmax,kmax,lmax
    integer :: ierror
    real(8), allocatable, dimension(:,:,:,:) :: QD
    real(8), allocatable, dimension(:,:,:)   :: QD_2D
    real(8), allocatable, dimension(:,:)     :: QD_1D

    character(20 ) :: affix
    character(100) :: fname

 800 format(i6.6,'_',i4.4)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)
    if (inum_b/=imax.or.jnum_b/=jmax.or.knum_b/=kmax) then
      write(6,*) 'Error::mpisub_sbsp_out_restart::Invalid arguments 1'
    endif

    call Get_DomainNum(ndomain)

    do id_dm=1,ndomain
      if (id_domain==id_dm) then
        call Get_BlockID_From_DomainID_And_SubDomainID(id_dm,1,id_blk)
        rank_domain_root = BlockID_To_Rank_Domain(id_blk)

        call Get_DomainIJKNum_From_DomainID(id_dm,imax,jmax,kmax)
        lmax = lnum_b

        if (myrank_domain==rank_domain_root) then
          allocate(QD(imax,jmax,kmax,lmax))
        else
          allocate(QD(1   ,1   ,1   ,1   ))
        endif

        if (myrank_domain==rank_domain_root) then
          allocate(QD_2D(imax,jmax,lmax))
          allocate(QD_1D(jmax,lmax))
        else
          allocate(QD_2D(1   ,1   ,1   ))
          allocate(QD_1D(1   ,1   ))
        endif

        call mpisub_sbsp_gatherv_from_domain(id_dm,QB,rank_domain_root,QD,lmax,0)

        write(affix,800) nstep,id_dm

        fname = trim(prefix)//trim(affix)//trim(suffix)

        if (myrank_domain==rank_domain_root) then

          QD_1D = 0.0d0
          QD_2D = 0.0d0

          do l=1,lmax
            do k=1,kmax
              do j=1,jmax
                do i=1,imax
                  QD_1D(j,l) = QD_1D(j,l) + QD(i,j,k,l)
                enddo
              enddo
            enddo
          enddo

          do l=1,lmax
            do j=1,jmax
              QD_1D(j,l) = QD_1D(j,l)/dble(imax*kmax)
            enddo
          enddo

          do l=1,lmax
            do j=1,jmax
              do i=1,imax
                QD_2D(i,j,l) = QD_1D(j,l)
              enddo
            enddo
          enddo

          open(19,file=fname,form='unformatted')
          write(19) imax,jmax,1,lmax
          write(19) QD_2D
          close(19)
        endif

        deallocate(QD)
        deallocate(QD_2D)
        deallocate(QD_1D)
      endif
#ifdef _USE_MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
#endif
    enddo

  end subroutine mpisub_sbsp_out_restart_mean_along_xsi_zta
  !-----------------------------------------------------------------------------
   subroutine mpisub_sbsp_out_mean_along_eta_zta(fname,QB,inum_b,jnum_b,knum_b,lnum_b)
    character(*), intent(in) :: fname
    real(8), intent(in), dimension(:,:,:,:) :: QB
    integer, intent(in) :: inum_b,jnum_b,knum_b,lnum_b

    integer :: i,j,k,l,ierror
    integer :: imax,jmax,kmax,lmax
    integer :: id_dm,ndomain
    real(8), allocatable, dimension(:,:,:,:) :: QD
    real(8), allocatable, dimension(:,:)     :: QD_1D

    call Get_DomainNum(ndomain)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)
    if (inum_b/=imax.or.jnum_b/=jmax.or.knum_b/=kmax) then
      write(6,*) 'Error:mpisub_sbsp_out_mean_along_eta_zta:Invalid arguments 1'
    endif

    if (ndomain/=1) then
      write(6,*) 'Error:mpisub_sbsp_out_mean_along_eta_zta:Invalid arguments 2'
    endif

    if (myrank_world==0) then
      open(21,file=fname,status='replace',form='formatted')
    endif


    do id_dm=1,ndomain
      call Get_DomainIJKNum_From_DomainID(id_dm,imax,jmax,kmax)
      lmax = lnum_b

      if (myrank_world==0) then
        allocate(QD(imax,jmax,kmax,lmax))
      else
        allocate(QD(1   ,1   ,1   ,1   ))
      endif

      if (myrank_world==0) then
        allocate(QD_1D(imax,lmax))
      else
        allocate(QD_1D(1   ,1   ))
      endif

      call mpisub_sbsp_send_domain_to_a_rank(id_dm,QB,0,QD,lmax,0)

      if (myrank_world==0) then

        QD_1D = 0.0d0

        do l=1,lmax
          do k=1,kmax
            do j=1,jmax
              do i=1,imax
                QD_1D(i,l) = QD_1D(i,l) + QD(i,j,k,l)
              enddo
            enddo
          enddo
        enddo

        do l=1,lmax
          do i=1,imax
            QD_1D(i,l) = QD_1D(i,l)/dble(jmax*kmax)
          enddo
        enddo

        do i=1,imax
          write(21,'(100e14.6)') real(QD_1D(i,1:lmax))
        enddo
      endif

      deallocate(QD)
      deallocate(QD_1D)
    enddo

    if (myrank_world==0) then
      close(21)
    endif

#ifdef _USE_MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
#endif
  end subroutine mpisub_sbsp_out_mean_along_eta_zta
  !-----------------------------------------------------------------------------
  subroutine mpisub_sbsp_send_domain_to_a_rank(id_dm,QB,rank_world_recv,QD,lnum,mode)
    integer, intent(in) :: id_dm
    integer, intent(in) :: rank_world_recv
    real(8), intent(in) , dimension(:,:,:,:) :: QB
    real(8), intent(out), dimension(:,:,:,:) :: QD
    integer, intent(in) :: lnum
    integer, intent(in) :: mode !0:cell 1:node

    integer :: rank_world_send
    integer :: rank_domain_root
    integer :: imax,jmax,kmax
    integer :: id_blk
    real(8), allocatable, dimension(:,:,:,:) :: QD_TMP
    integer :: sendcount,recvcount
    integer :: ierror
#ifdef _USE_MPI
    integer, dimension(MPI_STATUS_SIZE) :: istatus
#endif
    if (.not.((id_domain==id_dm).or.(myrank_world==rank_world_recv))) return


    call Get_BlockID_From_DomainID_And_SubDomainID(id_dm,1,id_blk)
    call Get_Rank_From_BlockID(id_blk,rank_world_send)
    rank_domain_root = BlockID_To_Rank_Domain(id_blk)

    call Get_DomainIJKNum_From_DomainID(id_dm,imax,jmax,kmax)
    if (mode==1) then
      imax = imax+1
      jmax = jmax+1
      kmax = kmax+1
    endif

    if (myrank_world==rank_world_send) then
      allocate(QD_TMP(imax,jmax,kmax,lnum))
    else
      allocate(QD_TMP(1   ,1   ,1   ,1   ))
    endif

    call mpisub_sbsp_gatherv_from_domain(id_dm,QB,rank_domain_root,QD_TMP,lnum,mode)

    sendcount = imax*jmax*kmax*lnum
    recvcount = imax*jmax*kmax*lnum

    if ((myrank_world==rank_world_send).and.(myrank_world==rank_world_recv)) then

      QD(:,:,:,:) = QD_TMP(:,:,:,:)

    else

      if (myrank_world==rank_world_send) then
#ifdef _USE_MPI
        call MPI_Send(QD_TMP,sendcount,MPI_DOUBLE_PRECISION,rank_world_recv,0, &
                      MPI_COMM_WORLD,ierror)
#endif
      endif
      if (myrank_world==rank_world_recv) then
#ifdef _USE_MPI
        call MPI_Recv(QD    ,recvcount,MPI_DOUBLE_PRECISION,rank_world_send,0, &
                      MPI_COMM_WORLD,istatus,ierror)
#endif
      endif

    endif

    deallocate(QD_TMP)
  end subroutine mpisub_sbsp_send_domain_to_a_rank
  !-----------------------------------------------------------------------------
  subroutine mpisub_sbsp_send_domainplane_to_domainplane(            &
                            id_dm_send,id_dir_send,id_split_send,QB, &
                            id_dm_recv,id_dir_recv,id_split_recv,QD,depth,lnum,mode)
    integer, intent(in) :: id_dm_send
    integer, intent(in) :: id_dir_send
    integer, intent(in) :: id_split_send
    integer, intent(in) :: id_dm_recv
    integer, intent(in) :: id_dir_recv
    integer, intent(in) :: id_split_recv
    real(8), intent(in) , dimension(:,:,:,:) :: QB
    real(8), intent(out), dimension(:,:,:,:) :: QD
    integer, intent(in) :: depth
    integer, intent(in) :: lnum
    integer, intent(in) :: mode !0:cell 1:node

    integer :: is,js,ks
    integer :: isplt,jsplt,ksplt
    integer :: imax,jmax,kmax
    integer :: is_sender,is_receiver
    integer :: id_blk,id_sbdm
    integer :: rank_world_recv
    integer :: rank_domainplane_root
    integer :: buffcount,ierror
    integer :: MPI_COMM_TEMP

    call Get_SplitID_From_DomainID_And_SubDomainID(id_domain,id_subdomain,is,js,ks)

    is_sender = 0
    is_receiver = 0

    if (id_dir_send==1) then
      if ((id_dm_send==id_domain).and.(is==id_split_send)) is_sender = 1
    else if (id_dir_send==2) then
      if ((id_dm_send==id_domain).and.(js==id_split_send)) is_sender = 1
    else if (id_dir_send==3) then
      if ((id_dm_send==id_domain).and.(ks==id_split_send)) is_sender = 1
    endif

    if (id_dir_recv==1) then
      if ((id_dm_recv==id_domain).and.(is==id_split_recv)) is_receiver = 1
    else if (id_dir_recv==2) then
      if ((id_dm_recv==id_domain).and.(js==id_split_recv)) is_receiver = 1
    else if (id_dir_recv==3) then
      if ((id_dm_recv==id_domain).and.(ks==id_split_recv)) is_receiver = 1
    endif

    if (.not.((is_sender==1).or.(is_receiver==1))) return

    call Get_DomainSplit_From_DomainID(id_dm_recv,isplt,jsplt,ksplt)
    if (id_dir_recv==1) then
      is = id_split_recv
      js = jsplt/2+mod(jsplt,2)
      ks = ksplt/2+mod(ksplt,2)
    else if (id_dir_recv==2) then
      is = isplt/2+mod(isplt,2)
      js = id_split_recv
      ks = ksplt/2+mod(ksplt,2)
    else if (id_dir_recv==3) then
      is = isplt/2+mod(isplt,2)
      js = jsplt/2+mod(jsplt,2)
      ks = id_split_recv
    endif
    call Get_SubDomainID_From_DomainID_And_SplitID(id_dm_recv,is,js,ks,id_sbdm)
    call Get_BlockID_From_DomainID_And_SubDomainID(id_dm_recv,id_sbdm,id_blk)
    call Get_Rank_From_BlockID(id_blk,rank_world_recv)
    if (id_dir_recv==1) then
      rank_domainplane_root = BlockID_To_Rank_DomainPlane_I(id_blk)
    else if (id_dir_recv==2) then
      rank_domainplane_root = BlockID_To_Rank_DomainPlane_J(id_blk)
    else if (id_dir_recv==3) then
      rank_domainplane_root = BlockID_To_Rank_DomainPlane_K(id_blk)
    endif

    if (id_dir_send==1) then
      call mpisub_sbsp_send_domainplane_i_to_a_rank(id_dm_send,id_split_send,QB, &
                                                    rank_world_recv,QD,depth,lnum,mode)
    else if (id_dir_send==2) then
      call mpisub_sbsp_send_domainplane_j_to_a_rank(id_dm_send,id_split_send,QB, &
                                                    rank_world_recv,QD,depth,lnum,mode)
    else if (id_dir_send==3) then
      call mpisub_sbsp_send_domainplane_k_to_a_rank(id_dm_send,id_split_send,QB, &
                                                    rank_world_recv,QD,depth,lnum,mode)
    endif

    call Get_DomainIJKNum_From_DomainID(id_dm_send,imax,jmax,kmax)
    if (mode==1) then
      imax = imax+1
      jmax = jmax+1
      kmax = kmax+1
    endif
    if (id_dir_send==1) then
      buffcount = depth*jmax*kmax*lnum
    else if (id_dir_send==2) then
      buffcount = imax*depth*kmax*lnum
    else if (id_dir_send==3) then
      buffcount = imax*jmax*depth*lnum
    endif

    if (id_dir_recv==1) then
      MPI_COMM_TEMP = MPI_COMM_DOMAINPLANE_I
    else if (id_dir_recv==2) then
      MPI_COMM_TEMP = MPI_COMM_DOMAINPLANE_J
    else if (id_dir_recv==3) then
      MPI_COMM_TEMP = MPI_COMM_DOMAINPLANE_K
    endif

    if (is_receiver==1) then
#ifdef _USE_MPI
      call MPI_Bcast(QD,buffcount,MPI_DOUBLE_PRECISION,rank_domainplane_root, &
                     MPI_COMM_TEMP,ierror)
#endif
    endif

  end subroutine mpisub_sbsp_send_domainplane_to_domainplane
  !-----------------------------------------------------------------------------
  subroutine mpisub_sbsp_send_domainplane_i_to_a_rank(id_dm,id_split,QB, &
                                                     rank_world_recv,QD,inum,lnum,mode)
    integer, intent(in) :: id_dm
    integer, intent(in) :: id_split
    integer, intent(in) :: rank_world_recv
    real(8), intent(in) , dimension(:,:,:,:) :: QB
    real(8), intent(out), dimension(:,:,:,:) :: QD
    integer, intent(in) :: inum,lnum
    integer, intent(in) :: mode !0:cell 1:node

    integer :: is,js,ks
    integer :: isplt,jsplt,ksplt
    integer :: imax,jmax,kmax
    integer :: rank_world_send
    integer :: rank_domainplane_root
    integer :: id_blk,id_sbdm
    real(8), allocatable, dimension(:,:,:,:) :: QD_TMP
    integer :: sendcount,recvcount
    integer :: ierror
#ifdef _USE_MPI
    integer, dimension(MPI_STATUS_SIZE) :: istatus
#endif
    call Get_SplitID_From_DomainID_And_SubDomainID(id_domain,id_subdomain,is,js,ks)

    if (.not.(((id_domain==id_dm).and.(is==id_split)).or.(myrank_world==rank_world_recv))) return


    call Get_DomainSplit_From_DomainID(id_dm,isplt,jsplt,ksplt)
    is = id_split
    js = jsplt/2+mod(jsplt,2)
    ks = ksplt/2+mod(ksplt,2)
    call Get_SubDomainID_From_DomainID_And_SplitID(id_dm,is,js,ks,id_sbdm)
    call Get_BlockID_From_DomainID_And_SubDomainID(id_dm,id_sbdm,id_blk)
    call Get_Rank_From_BlockID(id_blk,rank_world_send)
    rank_domainplane_root = BlockID_To_Rank_DomainPlane_I(id_blk)

    call Get_DomainIJKNum_From_DomainID(id_dm,imax,jmax,kmax)
    if (mode==1) then
      imax = imax+1
      jmax = jmax+1
      kmax = kmax+1
    endif

    if (myrank_world==rank_world_send) then
      allocate(QD_TMP(inum,jmax,kmax,lnum))
    else
      allocate(QD_TMP(1   ,1   ,1   ,1   ))
    endif

    call mpisub_sbsp_gatherv_from_domainplane_i(id_dm,id_split,QB    , &
                                         rank_domainplane_root,QD_TMP,inum,lnum,mode)

    sendcount = inum*jmax*kmax*lnum
    recvcount = inum*jmax*kmax*lnum

    if ((myrank_world==rank_world_send).and.(myrank_world==rank_world_recv)) then

      QD(:,:,:,:) = QD_TMP(:,:,:,:)

    else

      if (myrank_world==rank_world_send) then
#ifdef _USE_MPI
        call MPI_Send(QD_TMP,sendcount,MPI_DOUBLE_PRECISION,rank_world_recv,0, &
                      MPI_COMM_WORLD,ierror)
#endif
      endif
      if (myrank_world==rank_world_recv) then
#ifdef _USE_MPI
        call MPI_Recv(QD    ,recvcount,MPI_DOUBLE_PRECISION,rank_world_send,0, &
                      MPI_COMM_WORLD,istatus,ierror)
#endif
      endif

    endif

    deallocate(QD_TMP)
  end subroutine mpisub_sbsp_send_domainplane_i_to_a_rank
  !-----------------------------------------------------------------------------
  subroutine mpisub_sbsp_send_domainplane_j_to_a_rank(id_dm,id_split,QB, &
                                                     rank_world_recv,QD,jnum,lnum,mode)
    integer, intent(in) :: id_dm
    integer, intent(in) :: id_split
    integer, intent(in) :: rank_world_recv
    real(8), intent(in) , dimension(:,:,:,:) :: QB
    real(8), intent(out), dimension(:,:,:,:) :: QD
    integer, intent(in) :: jnum,lnum
    integer, intent(in) :: mode !0:cell 1:node

    integer :: is,js,ks
    integer :: isplt,jsplt,ksplt
    integer :: imax,jmax,kmax
    integer :: rank_world_send
    integer :: rank_domainplane_root
    integer :: id_blk,id_sbdm
    real(8), allocatable, dimension(:,:,:,:) :: QD_TMP
    integer :: sendcount,recvcount
    integer :: ierror
#ifdef _USE_MPI
    integer, dimension(MPI_STATUS_SIZE) :: istatus
#endif
    call Get_SplitID_From_DomainID_And_SubDomainID(id_domain,id_subdomain,is,js,ks)

    if (.not.(((id_domain==id_dm).and.(js==id_split)).or.(myrank_world==rank_world_recv))) return


    call Get_DomainSplit_From_DomainID(id_dm,isplt,jsplt,ksplt)
    is = isplt/2+mod(isplt,2)
    js = id_split
    ks = ksplt/2+mod(ksplt,2)
    call Get_SubDomainID_From_DomainID_And_SplitID(id_dm,is,js,ks,id_sbdm)
    call Get_BlockID_From_DomainID_And_SubDomainID(id_dm,id_sbdm,id_blk)
    call Get_Rank_From_BlockID(id_blk,rank_world_send)
    rank_domainplane_root = BlockID_To_Rank_DomainPlane_J(id_blk)

    call Get_DomainIJKNum_From_DomainID(id_dm,imax,jmax,kmax)
    if (mode==1) then
      imax = imax+1
      jmax = jmax+1
      kmax = kmax+1
    endif

    if (myrank_world==rank_world_send) then
      allocate(QD_TMP(imax,jnum,kmax,lnum))
    else
      allocate(QD_TMP(1   ,1   ,1   ,1   ))
    endif

    call mpisub_sbsp_gatherv_from_domainplane_j(id_dm,id_split,QB    , &
                                         rank_domainplane_root,QD_TMP,jnum,lnum,mode)

    sendcount = imax*jnum*kmax*lnum
    recvcount = imax*jnum*kmax*lnum

    if ((myrank_world==rank_world_send).and.(myrank_world==rank_world_recv)) then

      QD(:,:,:,:) = QD_TMP(:,:,:,:)

    else

      if (myrank_world==rank_world_send) then
#ifdef _USE_MPI
        call MPI_Send(QD_TMP,sendcount,MPI_DOUBLE_PRECISION,rank_world_recv,0, &
                      MPI_COMM_WORLD,ierror)
#endif
      endif
      if (myrank_world==rank_world_recv) then
#ifdef _USE_MPI
        call MPI_Recv(QD    ,recvcount,MPI_DOUBLE_PRECISION,rank_world_send,0, &
                      MPI_COMM_WORLD,istatus,ierror)
#endif
      endif

    endif

    deallocate(QD_TMP)
  end subroutine mpisub_sbsp_send_domainplane_j_to_a_rank
  !-----------------------------------------------------------------------------
  subroutine mpisub_sbsp_send_domainplane_k_to_a_rank(id_dm,id_split,QB, &
                                                     rank_world_recv,QD,knum,lnum,mode)
    integer, intent(in) :: id_dm
    integer, intent(in) :: id_split
    integer, intent(in) :: rank_world_recv
    real(8), intent(in ), dimension(:,:,:,:) :: QB
    real(8), intent(out), dimension(:,:,:,:) :: QD
    integer, intent(in) :: knum,lnum
    integer, intent(in) :: mode !0:cell 1:node

    integer :: is,js,ks
    integer :: isplt,jsplt,ksplt
    integer :: imax,jmax,kmax
    integer :: rank_world_send
    integer :: rank_domainplane_root
    integer :: id_blk,id_sbdm
    real(8), allocatable, dimension(:,:,:,:) :: QD_TMP
    integer :: sendcount,recvcount
    integer :: ierror
#ifdef _USE_MPI
    integer, dimension(MPI_STATUS_SIZE) :: istatus
#endif
    call Get_SplitID_From_DomainID_And_SubDomainID(id_domain,id_subdomain,is,js,ks)

    if (.not.(((id_domain==id_dm).and.(ks==id_split)).or.(myrank_world==rank_world_recv))) return


    call Get_DomainSplit_From_DomainID(id_dm,isplt,jsplt,ksplt)
    is = isplt/2+mod(isplt,2)
    js = jsplt/2+mod(jsplt,2)
    ks = id_split
    call Get_SubDomainID_From_DomainID_And_SplitID(id_dm,is,js,ks,id_sbdm)
    call Get_BlockID_From_DomainID_And_SubDomainID(id_dm,id_sbdm,id_blk)
    call Get_Rank_From_BlockID(id_blk,rank_world_send)
    rank_domainplane_root = BlockID_To_Rank_DomainPlane_K(id_blk)

    call Get_DomainIJKNum_From_DomainID(id_dm,imax,jmax,kmax)
    if (mode==1) then
      imax = imax+1
      jmax = jmax+1
      kmax = kmax+1
    endif

    if (myrank_world==rank_world_send) then
      allocate(QD_TMP(imax,jmax,knum,lnum))
    else
      allocate(QD_TMP(1   ,1   ,1   ,1   ))
    endif

    call mpisub_sbsp_gatherv_from_domainplane_k(id_dm,id_split,QB    , &
                                         rank_domainplane_root,QD_TMP,knum,lnum,mode)

    sendcount = imax*jmax*knum*lnum
    recvcount = imax*jmax*knum*lnum

    if ((myrank_world==rank_world_send).and.(myrank_world==rank_world_recv)) then

      QD(:,:,:,:) = QD_TMP(:,:,:,:)

    else

      if (myrank_world==rank_world_send) then
#ifdef _USE_MPI
        call MPI_Send(QD_TMP,sendcount,MPI_DOUBLE_PRECISION,rank_world_recv,0, &
                      MPI_COMM_WORLD,ierror)
#endif
      endif
      if (myrank_world==rank_world_recv) then
#ifdef _USE_MPI
        call MPI_Recv(QD    ,recvcount,MPI_DOUBLE_PRECISION,rank_world_send,0, &
                      MPI_COMM_WORLD,istatus,ierror)
#endif
      endif

    endif

    deallocate(QD_TMP)
  end subroutine mpisub_sbsp_send_domainplane_k_to_a_rank
  !-----------------------------------------------------------------------------
  subroutine mpisub_sbsp_gatherv_from_domain(id_dm,QB,rank_domain_root,QD,lnum,mode)
    integer, intent(in) :: id_dm
    integer, intent(in) :: rank_domain_root
    real(8), intent(in) , dimension(:,:,:,:) :: QB
    real(8), intent(out), dimension(:,:,:,:) :: QD
    integer, intent(in) :: lnum
    integer, intent(in) :: mode !0:cell 1:node

    integer :: id_blk,nsubdomain,id_sbdm,irank
    real(8), allocatable, dimension(:) :: sendbuff,recvbuff
    integer, allocatable, dimension(:) :: recvcount,displs
    integer :: sendcount,recvcountsum
    integer :: inum,jnum,knum,nhalo,is,js,ks
    integer :: istart,iend,jstart,jend,kstart,kend
    integer :: i,j,k,l,m,ii,jj,kk,moffst
    integer :: jstride,kstride,lstride
    integer :: ierror

    if (id_domain/=id_dm) return

    call Get_SubDomainNum_From_DomainID(id_dm,nsubdomain)

    if (nsubdomain==1) then
      QD(:,:,:,:) = QB(:,:,:,:)
      return
    endif

    allocate(recvcount(0:nsubdomain-1))
    allocate(displs(0:nsubdomain-1))

    do id_sbdm=1,nsubdomain
      call Get_BlockID_From_DomainID_And_SubDomainID(id_dm,id_sbdm,id_blk)
      call Get_BlockIJKNum_NoHalo_From_BlockID(id_blk,inum,jnum,knum)
      if (mode==1) then
        inum = inum+1
        jnum = jnum+1
        knum = knum+1
      endif
      irank = BlockID_To_Rank_Domain(id_blk)
      recvcount(irank) = inum*jnum*knum*lnum
    enddo

    do irank=0,nsubdomain-1
      if (irank==0) then
        displs(irank) = 0
        recvcountsum = recvcount(irank)
      else
        displs(irank) = displs(irank-1)+recvcount(irank-1)
        recvcountsum = recvcountsum+recvcount(irank)
      endif
    enddo

    call Get_HaloNum(nhalo)

    call Get_BlockIJKNum_NoHalo_From_BlockID(id_block,inum,jnum,knum)
    if (mode==1) then
      inum = inum+1
      jnum = jnum+1
      knum = knum+1
    endif

    sendcount = inum*jnum*knum*lnum

    allocate(sendbuff(sendcount))

    call Get_SplitID_From_DomainID_And_SubDomainID(id_domain,id_subdomain,is,js,ks)

    istart = 1
    if (is/=1) istart = nhalo+1
    iend = istart+inum-1
    jstart = 1
    if (js/=1) jstart = nhalo+1
    jend = jstart+jnum-1
    kstart = 1
    if (ks/=1) kstart = nhalo+1
    kend = kstart+knum-1

    jstride = inum
    kstride = inum*jnum
    lstride = inum*jnum*knum

    do l=1,lnum
      do k=kstart,kend
        do j=jstart,jend
          do i=istart,iend
            ii = i-istart+1
            jj = j-jstart+1
            kk = k-kstart+1
            m = ii+jstride*(jj-1)+kstride*(kk-1)+lstride*(l-1)
            sendbuff(m) = QB(i,j,k,l)
          enddo
        enddo
      enddo
    enddo

    if (myrank_domain==rank_domain_root) then
      allocate(recvbuff(recvcountsum))
    else
      allocate(recvbuff(1))
    endif

#ifdef _USE_MPI
    call MPI_Gatherv(sendbuff,sendcount,MPI_DOUBLE_PRECISION,        &
                     recvbuff,recvcount,displs,MPI_DOUBLE_PRECISION, &
                     rank_domain_root,MPI_COMM_DOMAIN,ierror)
#endif

    if (myrank_domain==rank_domain_root) then

      do id_sbdm=1,nsubdomain
        call Get_BlockID_From_DomainID_And_SubDomainID(id_dm,id_sbdm,id_blk)
        call Get_BlockIJKNum_NoHalo_From_BlockID(id_blk,inum,jnum,knum)
        if (mode==1) then
          inum = inum+1
          jnum = jnum+1
          knum = knum+1
        endif

        call Get_SubDomainIJKStartEnd_NoHalo_From_DomainID_And_SubDomainID(id_dm,id_sbdm, &
                          istart,iend,jstart,jend,kstart,kend)

        iend = istart+inum-1
        jend = jstart+jnum-1
        kend = kstart+knum-1

        irank = BlockID_To_Rank_Domain(id_blk)

        moffst = displs(irank)

        jstride = inum
        kstride = inum*jnum
        lstride = inum*jnum*knum

        do l=1,lnum
          do k=kstart,kend
            do j=jstart,jend
              do i=istart,iend
                ii = i-istart+1
                jj = j-jstart+1
                kk = k-kstart+1
                m  = ii+jstride*(jj-1)+kstride*(kk-1)+lstride*(l-1)+moffst
                QD(i,j,k,l) = recvbuff(m)
              enddo
            enddo
          enddo
        enddo

      enddo

    endif

    deallocate(recvbuff)
    deallocate(sendbuff)
    deallocate(recvcount)
    deallocate(displs)
  end subroutine mpisub_sbsp_gatherv_from_domain
  !-----------------------------------------------------------------------------
  subroutine mpisub_sbsp_gatherv_from_domainplane_i(id_dm,id_split,QB, &
                                             rank_domainplane_root,QD,inum,lnum,mode)
    integer, intent(in) :: id_dm
    integer, intent(in) :: id_split
    integer, intent(in) :: rank_domainplane_root
    real(8), intent(in) , dimension(:,:,:,:) :: QB
    real(8), intent(out), dimension(:,:,:,:) :: QD
    integer, intent(in) :: inum,lnum
    integer, intent(in) :: mode !0:cell 1:node

    integer :: i,j,k,l,m,ii,jj,kk,moffst
    integer :: is,js,ks
    integer :: jnum,knum,inum_tmp
    integer :: isplt,jsplt,ksplt
    real(8), allocatable, dimension(:) :: sendbuff,recvbuff
    integer, allocatable, dimension(:) :: recvcount,displs
    integer :: sendcount,recvcountsum
    integer :: id_blk,id_sbdm,irank
    integer :: istart,iend,jstart,jend,kstart,kend
    integer :: nhalo
    integer :: jstride,kstride,lstride
    integer :: ierror

    call Get_SplitID_From_DomainID_And_SubDomainID(id_domain,id_subdomain,is,js,ks)

    if (.not.(id_domain==id_dm.and.is==id_split)) return


    call Get_DomainSplit_From_DomainID(id_dm,isplt,jsplt,ksplt)

    if (jsplt*ksplt==1) then
      QD(:,:,:,:) = QB(:,:,:,:)
      return
    endif

    allocate(recvcount(0:jsplt*ksplt-1))
    allocate(displs(0:jsplt*ksplt-1))

    do ks=1,ksplt
      do js=1,jsplt
        is=id_split
        call Get_SubDomainID_From_DomainID_And_SplitID(id_dm,is,js,ks,id_sbdm)
        call Get_BlockID_From_DomainID_And_SubDomainID(id_dm,id_sbdm,id_blk)
        call Get_BlockIJKNum_NoHalo_From_BlockID(id_blk,inum_tmp,jnum,knum)
        if (mode==1) then
          jnum = jnum+1
          knum = knum+1
        endif
        irank = BlockID_To_Rank_DomainPlane_I(id_blk)
        recvcount(irank) = inum*jnum*knum*lnum
      enddo
    enddo

    do irank=0,jsplt*ksplt-1
      if (irank==0) then
        displs(irank) = 0
        recvcountsum = recvcount(irank)
      else
        displs(irank) = displs(irank-1)+recvcount(irank-1)
        recvcountsum = recvcountsum+recvcount(irank)
      endif
    enddo

    call Get_HaloNum(nhalo)

    call Get_BlockIJKNum_NoHalo_From_BlockID(id_block,inum_tmp,jnum,knum)
    if (mode==1) then
      jnum = jnum+1
      knum = knum+1
    endif

    call Get_SplitID_From_DomainID_And_SubDomainID(id_domain,id_subdomain,is,js,ks)

    istart = 1
    iend = istart+inum-1
    jstart = 1
    if (js/=1) jstart = nhalo+1
    jend = jstart+jnum-1
    kstart = 1
    if (ks/=1) kstart = nhalo+1
    kend = kstart+knum-1

    sendcount = inum*jnum*knum*lnum

    jstride = inum
    kstride = inum*jnum
    lstride = inum*jnum*knum

    allocate(sendbuff(sendcount))

    do l=1,lnum
      do k=kstart,kend
        do j=jstart,jend
          do i=istart,iend
            ii = i-istart+1
            jj = j-jstart+1
            kk = k-kstart+1
            m = ii+jstride*(jj-1)+kstride*(kk-1)+lstride*(l-1)
            sendbuff(m) = QB(i,j,k,l)
          enddo
        enddo
      enddo
    enddo

    if (myrank_domainplane_i==rank_domainplane_root) then
      allocate(recvbuff(recvcountsum))
    else
      allocate(recvbuff(1))
    endif

#ifdef _USE_MPI
    call MPI_Gatherv(sendbuff,sendcount,MPI_DOUBLE_PRECISION,        &
                     recvbuff,recvcount,displs,MPI_DOUBLE_PRECISION, &
                     rank_domainplane_root,MPI_COMM_DOMAINPLANE_I,ierror)
#endif

    if (myrank_domainplane_i==rank_domainplane_root) then

      do ks=1,ksplt
        do js=1,jsplt
          is=id_split
          call Get_SubDomainID_From_DomainID_And_SplitID(id_dm,is,js,ks,id_sbdm)
          call Get_BlockID_From_DomainID_And_SubDomainID(id_dm,id_sbdm,id_blk)
          call Get_BlockIJKNum_NoHalo_From_BlockID(id_blk,inum_tmp,jnum,knum)
          if (mode==1) then
            jnum = jnum+1
            knum = knum+1
          endif

          irank = BlockID_To_Rank_DomainPlane_I(id_blk)

          call Get_SubDomainIJKStartEnd_NoHalo_From_DomainID_And_SubDomainID(id_dm,id_sbdm, &
                istart,iend,jstart,jend,kstart,kend)
          istart = 1
          iend = istart+inum-1
          jend = jstart+jnum-1
          kend = kstart+knum-1

          moffst = displs(irank)

          jstride = inum
          kstride = inum*jnum
          lstride = inum*jnum*knum

          do l=1,lnum
            do k=kstart,kend
              do j=jstart,jend
                do i=istart,iend
                  ii = i-istart+1
                  jj = j-jstart+1
                  kk = k-kstart+1
                  m  = ii+jstride*(jj-1)+kstride*(kk-1)+lstride*(l-1)+moffst
                  QD(i,j,k,l) = recvbuff(m)
                enddo
              enddo
            enddo
          enddo

        enddo
      enddo

    endif

    deallocate(sendbuff)
    deallocate(recvbuff)
    deallocate(recvcount)
    deallocate(displs)
  end subroutine mpisub_sbsp_gatherv_from_domainplane_i
  !-----------------------------------------------------------------------------
  subroutine mpisub_sbsp_gatherv_from_domainplane_j(id_dm,id_split,QB, &
                                             rank_domainplane_root,QD,jnum,lnum,mode)
    integer, intent(in) :: id_dm
    integer, intent(in) :: id_split
    integer, intent(in) :: rank_domainplane_root
    real(8), intent(in) , dimension(:,:,:,:) :: QB
    real(8), intent(out), dimension(:,:,:,:) :: QD
    integer, intent(in) :: jnum,lnum
    integer, intent(in) :: mode !0:cell 1:node

    integer :: i,j,k,l,m,ii,jj,kk,moffst
    integer :: is,js,ks
    integer :: inum,knum,jnum_tmp
    integer :: isplt,jsplt,ksplt
    real(8), allocatable, dimension(:) :: sendbuff,recvbuff
    integer, allocatable, dimension(:) :: recvcount,displs
    integer :: sendcount,recvcountsum
    integer :: id_blk,id_sbdm,irank
    integer :: istart,iend,jstart,jend,kstart,kend
    integer :: nhalo
    integer :: jstride,kstride,lstride
    integer :: ierror

    call Get_SplitID_From_DomainID_And_SubDomainID(id_domain,id_subdomain,is,js,ks)

    if (.not.(id_domain==id_dm.and.js==id_split)) return


    call Get_DomainSplit_From_DomainID(id_dm,isplt,jsplt,ksplt)

    if (isplt*ksplt==1) then
      QD(:,:,:,:) = QB(:,:,:,:)
      return
    endif

    allocate(recvcount(0:isplt*ksplt-1))
    allocate(displs(0:isplt*ksplt-1))

    do ks=1,ksplt
      do is=1,isplt
        js=id_split
        call Get_SubDomainID_From_DomainID_And_SplitID(id_dm,is,js,ks,id_sbdm)
        call Get_BlockID_From_DomainID_And_SubDomainID(id_dm,id_sbdm,id_blk)
        call Get_BlockIJKNum_NoHalo_From_BlockID(id_blk,inum,jnum_tmp,knum)
        if (mode==1) then
          inum = inum+1
          knum = knum+1
        endif
        irank = BlockID_To_Rank_DomainPlane_J(id_blk)
        recvcount(irank) = inum*jnum*knum*lnum
      enddo
    enddo

    do irank=0,isplt*ksplt-1
      if (irank==0) then
        displs(irank) = 0
        recvcountsum = recvcount(irank)
      else
        displs(irank) = displs(irank-1)+recvcount(irank-1)
        recvcountsum = recvcountsum+recvcount(irank)
      endif
    enddo

    call Get_HaloNum(nhalo)

    call Get_BlockIJKNum_NoHalo_From_BlockID(id_block,inum,jnum_tmp,knum)
    if (mode==1) then
      inum = inum+1
      knum = knum+1
    endif

    call Get_SplitID_From_DomainID_And_SubDomainID(id_domain,id_subdomain,is,js,ks)

    istart = 1
    if (is/=1) istart = nhalo+1
    iend = istart+inum-1
    jstart = 1
    jend = jstart+jnum-1
    kstart = 1
    if (ks/=1) kstart = nhalo+1
    kend = kstart+knum-1

    sendcount = inum*jnum*knum*lnum

    jstride = inum
    kstride = inum*jnum
    lstride = inum*jnum*knum

    allocate(sendbuff(sendcount))

    do l=1,lnum
      do k=kstart,kend
        do j=jstart,jend
          do i=istart,iend
            ii = i-istart+1
            jj = j-jstart+1
            kk = k-kstart+1
            m = ii+jstride*(jj-1)+kstride*(kk-1)+lstride*(l-1)
            sendbuff(m) = QB(i,j,k,l)
          enddo
        enddo
      enddo
    enddo

    if (myrank_domainplane_j==rank_domainplane_root) then
      allocate(recvbuff(recvcountsum))
    else
      allocate(recvbuff(1))
    endif

#ifdef _USE_MPI
    call MPI_Gatherv(sendbuff,sendcount,MPI_DOUBLE_PRECISION,        &
                     recvbuff,recvcount,displs,MPI_DOUBLE_PRECISION, &
                     rank_domainplane_root,MPI_COMM_DOMAINPLANE_J,ierror)
#endif

    if (myrank_domainplane_j==rank_domainplane_root) then

      do ks=1,ksplt
        do is=1,isplt
          js=id_split
          call Get_SubDomainID_From_DomainID_And_SplitID(id_dm,is,js,ks,id_sbdm)
          call Get_BlockID_From_DomainID_And_SubDomainID(id_dm,id_sbdm,id_blk)
          call Get_BlockIJKNum_NoHalo_From_BlockID(id_blk,inum,jnum_tmp,knum)
          if (mode==1) then
            inum = inum+1
            knum = knum+1
          endif

          irank = BlockID_To_Rank_DomainPlane_J(id_blk)

          call Get_SubDomainIJKStartEnd_NoHalo_From_DomainID_And_SubDomainID(id_dm,id_sbdm, &
                istart,iend,jstart,jend,kstart,kend)
          iend = istart+inum-1
          jstart = 1
          jend = jstart+jnum-1
          kend = kstart+knum-1

          moffst = displs(irank)

          jstride = inum
          kstride = inum*jnum
          lstride = inum*jnum*knum

          do l=1,lnum
            do k=kstart,kend
              do j=jstart,jend
                do i=istart,iend
                  ii = i-istart+1
                  jj = j-jstart+1
                  kk = k-kstart+1
                  m  = ii+jstride*(jj-1)+kstride*(kk-1)+lstride*(l-1)+moffst
                  QD(i,j,k,l) = recvbuff(m)
                enddo
              enddo
            enddo
          enddo

        enddo
      enddo

    endif

    deallocate(sendbuff)
    deallocate(recvbuff)
    deallocate(recvcount)
    deallocate(displs)
  end subroutine mpisub_sbsp_gatherv_from_domainplane_j
  !-----------------------------------------------------------------------------
  subroutine mpisub_sbsp_gatherv_from_domainplane_k(id_dm,id_split,QB, &
                                             rank_domainplane_root,QD,knum,lnum,mode)
    integer, intent(in) :: id_dm
    integer, intent(in) :: id_split
    integer, intent(in) :: rank_domainplane_root
    real(8), intent(in) , dimension(:,:,:,:) :: QB
    real(8), intent(out), dimension(:,:,:,:) :: QD
    integer, intent(in) :: knum,lnum
    integer, intent(in) :: mode !0:cell 1:node

    integer :: i,j,k,l,m,ii,jj,kk,moffst
    integer :: is,js,ks
    integer :: inum,jnum,knum_tmp
    integer :: isplt,jsplt,ksplt
    real(8), allocatable, dimension(:) :: sendbuff,recvbuff
    integer, allocatable, dimension(:) :: recvcount,displs
    integer :: sendcount,recvcountsum
    integer :: id_blk,id_sbdm,irank
    integer :: istart,iend,jstart,jend,kstart,kend
    integer :: nhalo
    integer :: jstride,kstride,lstride
    integer :: ierror

    call Get_SplitID_From_DomainID_And_SubDomainID(id_domain,id_subdomain,is,js,ks)

    if (.not.(id_domain==id_dm.and.ks==id_split)) return


    call Get_DomainSplit_From_DomainID(id_dm,isplt,jsplt,ksplt)

    if (isplt*jsplt==1) then
      QD(:,:,:,:) = QB(:,:,:,:)
      return
    endif

    allocate(recvcount(0:isplt*jsplt-1))
    allocate(displs(0:isplt*jsplt-1))

    do js=1,jsplt
      do is=1,isplt
        ks=id_split
        call Get_SubDomainID_From_DomainID_And_SplitID(id_dm,is,js,ks,id_sbdm)
        call Get_BlockID_From_DomainID_And_SubDomainID(id_dm,id_sbdm,id_blk)
        call Get_BlockIJKNum_NoHalo_From_BlockID(id_blk,inum,jnum,knum_tmp)
        irank = BlockID_To_Rank_DomainPlane_K(id_blk)
        if (mode==1) then
          inum = inum+1
          jnum = jnum+1
        endif
        recvcount(irank) = inum*jnum*knum*lnum
      enddo
    enddo

    do irank=0,isplt*jsplt-1
      if (irank==0) then
        displs(irank) = 0
        recvcountsum = recvcount(irank)
      else
        displs(irank) = displs(irank-1)+recvcount(irank-1)
        recvcountsum = recvcountsum+recvcount(irank)
      endif
    enddo

    call Get_HaloNum(nhalo)

    call Get_BlockIJKNum_NoHalo_From_BlockID(id_block,inum,jnum,knum_tmp)
    if (mode==1) then
      inum = inum+1
      jnum = jnum+1
    endif

    call Get_SplitID_From_DomainID_And_SubDomainID(id_domain,id_subdomain,is,js,ks)

    istart = 1
    if (is/=1) istart = nhalo+1
    iend = istart+inum-1
    jstart = 1
    if (js/=1) jstart = nhalo+1
    jend = jstart+jnum-1
    kstart = 1
    kend = kstart+knum-1

    sendcount = inum*jnum*knum*lnum

    jstride = inum
    kstride = inum*jnum
    lstride = inum*jnum*knum

    allocate(sendbuff(sendcount))

    do l=1,lnum
      do k=kstart,kend
        do j=jstart,jend
          do i=istart,iend
            ii = i-istart+1
            jj = j-jstart+1
            kk = k-kstart+1
            m = ii+jstride*(jj-1)+kstride*(kk-1)+lstride*(l-1)
            sendbuff(m) = QB(i,j,k,l)
          enddo
        enddo
      enddo
    enddo

    if (myrank_domainplane_k==rank_domainplane_root) then
      allocate(recvbuff(recvcountsum))
    else
      allocate(recvbuff(1))
    endif

#ifdef _USE_MPI
    call MPI_Gatherv(sendbuff,sendcount,MPI_DOUBLE_PRECISION,        &
                     recvbuff,recvcount,displs,MPI_DOUBLE_PRECISION, &
                     rank_domainplane_root,MPI_COMM_DOMAINPLANE_K,ierror)
#endif

    if (myrank_domainplane_k==rank_domainplane_root) then

      do js=1,jsplt
        do is=1,isplt
          ks=id_split
          call Get_SubDomainID_From_DomainID_And_SplitID(id_dm,is,js,ks,id_sbdm)
          call Get_BlockID_From_DomainID_And_SubDomainID(id_dm,id_sbdm,id_blk)
          call Get_BlockIJKNum_NoHalo_From_BlockID(id_blk,inum,jnum,knum_tmp)
          if (mode==1) then
            inum = inum+1
            jnum = jnum+1
          endif

          irank = BlockID_To_Rank_DomainPlane_K(id_blk)

          call Get_SubDomainIJKStartEnd_NoHalo_From_DomainID_And_SubDomainID(id_dm,id_sbdm, &
                istart,iend,jstart,jend,kstart,kend)
          iend = istart+inum-1
          jend = jstart+jnum-1
          kstart = 1
          kend = kstart+knum-1

          moffst = displs(irank)

          jstride = inum
          kstride = inum*jnum
          lstride = inum*jnum*knum

          do l=1,lnum
            do k=kstart,kend
              do j=jstart,jend
                do i=istart,iend
                  ii = i-istart+1
                  jj = j-jstart+1
                  kk = k-kstart+1
                  m  = ii+jstride*(jj-1)+kstride*(kk-1)+lstride*(l-1)+moffst
                  QD(i,j,k,l) = recvbuff(m)
                enddo
              enddo
            enddo
          enddo

        enddo
      enddo

    endif

    deallocate(sendbuff)
    deallocate(recvbuff)
    deallocate(recvcount)
    deallocate(displs)
  end subroutine mpisub_sbsp_gatherv_from_domainplane_k
  !-----------------------------------------------------------------------------
  subroutine mpisub_sbsp_exchange_subdomain_halo_xsi(QB_1_M,QB_1_P, &
                                                     QB_2_M,QB_2_P,lmax,depth,iflag,mode)
    real(8), intent(inout), dimension(:,:,:,:) :: QB_1_M,QB_2_M
    real(8), intent(inout), dimension(:,:,:,:) :: QB_1_P,QB_2_P
    integer, intent(in) :: lmax
    integer, intent(in) :: depth
    integer, intent(out), dimension(:) :: iflag
    integer, intent(in) :: mode !0:cell 1:node

    integer :: n
    integer :: istart,iend
    integer :: is_blk,js_blk,ks_blk,isplt_blk,jsplt_blk,ksplt_blk
    integer :: id_sbdm1,id_sbdm2
    integer :: id_blk1,id_blk2
    integer :: irank1,irank2
    integer :: imax,jmax,kmax
    integer :: sendrecvnum
    integer :: ierror
#ifdef _USE_MPI
    integer, dimension(MPI_STATUS_SIZE) :: istatus
#endif
    call Get_DomainSplit_From_DomainID(id_domain,isplt_blk,jsplt_blk,ksplt_blk)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)
    if (mode==1) then
      imax = imax+1
      jmax = jmax+1
      kmax = kmax+1
    endif

    sendrecvnum = depth*jmax*kmax*lmax

    do n=1,2
      istart=n
      iend=isplt_blk
      do ks_blk=1,ksplt_blk
        do js_blk=1,jsplt_blk
          do is_blk=istart,iend,2

            if (is_blk+1<=isplt_blk) then
              call Get_SubDomainID_From_DomainID_And_SplitID(id_domain,is_blk  ,js_blk,ks_blk,id_sbdm1)
              call Get_SubDomainID_From_DomainID_And_SplitID(id_domain,is_blk+1,js_blk,ks_blk,id_sbdm2)
              call Get_BlockID_From_DomainID_And_SubDomainID(id_domain,id_sbdm1,id_blk1)
              call Get_BlockID_From_DomainID_And_SubDomainID(id_domain,id_sbdm2,id_blk2)

              irank1 = BlockID_To_Rank_Domain(id_blk1)
              irank2 = BlockID_To_Rank_Domain(id_blk2)

              if (irank1==myrank_domain) then
#ifdef _USE_MPI
                call MPI_Send(QB_2_M,sendrecvnum,MPI_DOUBLE_PRECISION,irank2,0, &
                              MPI_COMM_DOMAIN,ierror)
#endif
#ifdef _USE_MPI
                call MPI_Recv(QB_2_P,sendrecvnum,MPI_DOUBLE_PRECISION,irank2,0, &
                              MPI_COMM_DOMAIN,istatus,ierror)
#endif
                iflag(2) = 1

              endif

              if (irank2==myrank_domain) then
#ifdef _USE_MPI
                call MPI_Recv(QB_1_P,sendrecvnum,MPI_DOUBLE_PRECISION,irank1,0, &
                              MPI_COMM_DOMAIN,istatus,ierror)
#endif
                iflag(1) = 1
#ifdef _USE_MPI
                call MPI_Send(QB_1_M,sendrecvnum,MPI_DOUBLE_PRECISION,irank1,0, &
                              MPI_COMM_DOMAIN,ierror)
#endif
              endif

            endif

          enddo
        enddo
      enddo
    enddo

  end subroutine mpisub_sbsp_exchange_subdomain_halo_xsi
!-------------------------------------------------------------------------------
  subroutine mpisub_sbsp_exchange_subdomain_halo_eta(QB_3_M,QB_3_P, &
                                                     QB_4_M,QB_4_P,lmax,depth,iflag,mode)
    real(8), intent(inout), dimension(:,:,:,:) :: QB_3_M,QB_4_M
    real(8), intent(inout), dimension(:,:,:,:) :: QB_3_P,QB_4_P
    integer, intent(in) :: lmax
    integer, intent(in) :: depth
    integer, intent(out), dimension(:) :: iflag
    integer, intent(in) :: mode !0:cell 1:halo

    integer :: n
    integer :: jstart,jend
    integer :: is_blk,js_blk,ks_blk,isplt_blk,jsplt_blk,ksplt_blk
    integer :: id_sbdm1,id_sbdm2
    integer :: id_blk1,id_blk2
    integer :: irank1,irank2
    integer :: imax,jmax,kmax
    integer :: sendrecvnum
    integer :: ierror
#ifdef _USE_MPI
    integer, dimension(MPI_STATUS_SIZE) :: istatus
#endif
    call Get_DomainSplit_From_DomainID(id_domain,isplt_blk,jsplt_blk,ksplt_blk)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)
    if (mode==1) then
      imax = imax+1
      jmax = jmax+1
      kmax = kmax+1
    endif

    sendrecvnum = imax*depth*kmax*lmax

    do n=1,2
      jstart=n
      jend=jsplt_blk
      do ks_blk=1,ksplt_blk
        do is_blk=1,isplt_blk
          do js_blk=jstart,jend,2

            if (js_blk+1<=jsplt_blk) then
              call Get_SubDomainID_From_DomainID_And_SplitID(id_domain,is_blk,js_blk  ,ks_blk,id_sbdm1)
              call Get_SubDomainID_From_DomainID_And_SplitID(id_domain,is_blk,js_blk+1,ks_blk,id_sbdm2)
              call Get_BlockID_From_DomainID_And_SubDomainID(id_domain,id_sbdm1,id_blk1)
              call Get_BlockID_From_DomainID_And_SubDomainID(id_domain,id_sbdm2,id_blk2)

              irank1 = BlockID_To_Rank_Domain(id_blk1)
              irank2 = BlockID_To_Rank_Domain(id_blk2)

              if (irank1==myrank_domain) then
#ifdef _USE_MPI
                call MPI_Send(QB_4_M,sendrecvnum,MPI_DOUBLE_PRECISION,irank2,0, &
                              MPI_COMM_DOMAIN,ierror)
#endif
#ifdef _USE_MPI
                call MPI_Recv(QB_4_P,sendrecvnum,MPI_DOUBLE_PRECISION,irank2,0, &
                              MPI_COMM_DOMAIN,istatus,ierror)
#endif
                iflag(4) = 1

              endif

              if (irank2==myrank_domain) then
#ifdef _USE_MPI
                call MPI_Recv(QB_3_P,sendrecvnum,MPI_DOUBLE_PRECISION,irank1,0, &
                              MPI_COMM_DOMAIN,istatus,ierror)
#endif
                iflag(3) = 1
#ifdef _USE_MPI
                call MPI_Send(QB_3_M,sendrecvnum,MPI_DOUBLE_PRECISION,irank1,0, &
                              MPI_COMM_DOMAIN,ierror)
#endif
              endif

            endif

          enddo
        enddo
      enddo
    enddo

  end subroutine mpisub_sbsp_exchange_subdomain_halo_eta
!-------------------------------------------------------------------------------
  subroutine mpisub_sbsp_exchange_subdomain_halo_zta(QB_5_M,QB_5_P, &
                                                     QB_6_M,QB_6_P,lmax,depth,iflag,mode)
    real(8), intent(inout), dimension(:,:,:,:) :: QB_5_M,QB_6_M
    real(8), intent(inout), dimension(:,:,:,:) :: QB_5_P,QB_6_P
    integer, intent(in) :: lmax
    integer, intent(in) :: depth
    integer, intent(out), dimension(:) :: iflag
    integer, intent(in) :: mode !0:cell 1:node

    integer :: n
    integer :: kstart,kend
    integer :: is_blk,js_blk,ks_blk,isplt_blk,jsplt_blk,ksplt_blk
    integer :: id_sbdm1,id_sbdm2
    integer :: id_blk1,id_blk2
    integer :: irank1,irank2
    integer :: imax,jmax,kmax
    integer :: sendrecvnum
    integer :: ierror
#ifdef _USE_MPI
    integer, dimension(MPI_STATUS_SIZE) :: istatus
#endif
    call Get_DomainSplit_From_DomainID(id_domain,isplt_blk,jsplt_blk,ksplt_blk)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)
    if (mode==1) then
      imax = imax+1
      jmax = jmax+1
      kmax = kmax+1
    endif

    sendrecvnum = imax*jmax*depth*lmax

    do n=1,2
      kstart=n
      kend=ksplt_blk
      do js_blk=1,jsplt_blk
        do is_blk=1,isplt_blk
          do ks_blk=kstart,kend,2

            if (ks_blk+1<=ksplt_blk) then
              call Get_SubDomainID_From_DomainID_And_SplitID(id_domain,is_blk,js_blk,ks_blk  ,id_sbdm1)
              call Get_SubDomainID_From_DomainID_And_SplitID(id_domain,is_blk,js_blk,ks_blk+1,id_sbdm2)
              call Get_BlockID_From_DomainID_And_SubDomainID(id_domain,id_sbdm1,id_blk1)
              call Get_BlockID_From_DomainID_And_SubDomainID(id_domain,id_sbdm2,id_blk2)

              irank1 = BlockID_To_Rank_Domain(id_blk1)
              irank2 = BlockID_To_Rank_Domain(id_blk2)

              if (irank1==myrank_domain) then
#ifdef _USE_MPI
                call MPI_Send(QB_6_M,sendrecvnum,MPI_DOUBLE_PRECISION,irank2,0, &
                              MPI_COMM_DOMAIN,ierror)
#endif
#ifdef _USE_MPI
                call MPI_Recv(QB_6_P,sendrecvnum,MPI_DOUBLE_PRECISION,irank2,0, &
                              MPI_COMM_DOMAIN,istatus,ierror)
#endif
                iflag(6) = 1

              endif

              if (irank2==myrank_domain) then
#ifdef _USE_MPI
                call MPI_Recv(QB_5_P,sendrecvnum,MPI_DOUBLE_PRECISION,irank1,0, &
                              MPI_COMM_DOMAIN,istatus,ierror)
#endif
                iflag(5) = 1
#ifdef _USE_MPI
                call MPI_Send(QB_5_M,sendrecvnum,MPI_DOUBLE_PRECISION,irank1,0, &
                              MPI_COMM_DOMAIN,ierror)
#endif
              endif

            endif

          enddo
        enddo
      enddo
    enddo

  end subroutine mpisub_sbsp_exchange_subdomain_halo_zta
!-------------------------------------------------------------------------------
  subroutine mpisub_sbsp_exchange_domain_halo_xsi(QB_1_M,QB_1_P, &
                                                  QB_2_M,QB_2_P,lmax,depth,iflag,mode)
    real(8), intent(inout), dimension(:,:,:,:) :: QB_1_M,QB_2_M
    real(8), intent(inout), dimension(:,:,:,:) :: QB_1_P,QB_2_P
    integer, intent(in) :: lmax
    integer, intent(in) :: depth
    integer, intent(out), dimension(:) :: iflag
    integer, intent(in) :: mode !0:cell 1:node

    integer :: imax,jmax,kmax
    integer :: isplt_blk,jsplt_blk,ksplt_blk,js_blk,ks_blk
    integer :: id_sbdm1,id_blk1,irank1
    integer :: id_sbdm2,id_blk2,irank2
    integer :: ierror
#ifdef _USE_MPI
    integer, dimension(MPI_STATUS_SIZE) :: istatus
#endif
    integer :: sendrecvnum

    call Get_DomainSplit_From_DomainID(id_domain,isplt_blk,jsplt_blk,ksplt_blk)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)
    if (mode==1) then
      imax = imax+1
      jmax = jmax+1
      kmax = kmax+1
    endif

    sendrecvnum = depth*jmax*kmax*lmax

    !i direction
    do ks_blk=1,ksplt_blk
      do js_blk=1,jsplt_blk

        call Get_SubDomainID_From_DomainID_And_SplitID(id_domain,isplt_blk,js_blk,ks_blk,id_sbdm1)
        call Get_SubDomainID_From_DomainID_And_SplitID(id_domain,1        ,js_blk,ks_blk,id_sbdm2)
        call Get_BlockID_From_DomainID_And_SubDomainID(id_domain,id_sbdm1,id_blk1)
        call Get_BlockID_From_DomainID_And_SubDomainID(id_domain,id_sbdm2,id_blk2)

        irank1 = BlockID_To_Rank_Domain(id_blk1)
        irank2 = BlockID_To_Rank_Domain(id_blk2)

        if (.not.(irank1==myrank_domain.or.irank2==myrank_domain)) cycle

        if (irank1==irank2) then

          QB_2_P(:,:,:,:) = QB_1_M(:,:,:,:)
          iflag(2) = 1
          QB_1_P(:,:,:,:) = QB_2_M(:,:,:,:)
          iflag(1) = 1

        else

          if (irank1==myrank_domain) then
#ifdef _USE_MPI
            call MPI_Send(QB_2_M,sendrecvnum,MPI_DOUBLE_PRECISION,irank2,0, &
                          MPI_COMM_DOMAIN,ierror)
#endif
#ifdef _USE_MPI
            call MPI_Recv(QB_2_P,sendrecvnum,MPI_DOUBLE_PRECISION,irank2,0, &
                          MPI_COMM_DOMAIN,istatus,ierror)
#endif
            iflag(2) = 1

          endif

          if (irank2==myrank_domain) then
#ifdef _USE_MPI
            call MPI_Recv(QB_1_P,sendrecvnum,MPI_DOUBLE_PRECISION,irank1,0, &
                          MPI_COMM_DOMAIN,istatus,ierror)
#endif
            iflag(1) = 1
#ifdef _USE_MPI
            call MPI_Send(QB_1_M,sendrecvnum,MPI_DOUBLE_PRECISION,irank1,0, &
                          MPI_COMM_DOMAIN,ierror)
#endif
          endif

        endif

      enddo
    enddo

  end subroutine mpisub_sbsp_exchange_domain_halo_xsi
!-------------------------------------------------------------------------------
  subroutine mpisub_sbsp_exchange_domain_halo_eta(QB_3_M,QB_3_P, &
                                                  QB_4_M,QB_4_P,lmax,depth,iflag,mode)
    real(8), intent(inout), dimension(:,:,:,:) :: QB_3_M,QB_4_M
    real(8), intent(inout), dimension(:,:,:,:) :: QB_3_P,QB_4_P
    integer, intent(in) :: lmax
    integer, intent(in) :: depth
    integer, intent(out), dimension(:) :: iflag
    integer, intent(in) :: mode !0:cell 1:node

    integer :: imax,jmax,kmax
    integer :: isplt_blk,jsplt_blk,ksplt_blk,is_blk,ks_blk
    integer :: id_sbdm1,id_blk1,irank1
    integer :: id_sbdm2,id_blk2,irank2
    integer :: ierror
#ifdef _USE_MPI
    integer, dimension(MPI_STATUS_SIZE) :: istatus
#endif
    integer :: sendrecvnum

    call Get_DomainSplit_From_DomainID(id_domain,isplt_blk,jsplt_blk,ksplt_blk)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)
    if (mode==1) then
      imax = imax+1
      jmax = jmax+1
      kmax = kmax+1
    endif

    sendrecvnum = imax*depth*kmax*lmax

    !j direction
    do ks_blk=1,ksplt_blk
      do is_blk=1,isplt_blk

        call Get_SubDomainID_From_DomainID_And_SplitID(id_domain,is_blk,jsplt_blk,ks_blk,id_sbdm1)
        call Get_SubDomainID_From_DomainID_And_SplitID(id_domain,is_blk,1        ,ks_blk,id_sbdm2)
        call Get_BlockID_From_DomainID_And_SubDomainID(id_domain,id_sbdm1,id_blk1)
        call Get_BlockID_From_DomainID_And_SubDomainID(id_domain,id_sbdm2,id_blk2)

        irank1 = BlockID_To_Rank_Domain(id_blk1)
        irank2 = BlockID_To_Rank_Domain(id_blk2)

        if (.not.(irank1==myrank_domain.or.irank2==myrank_domain)) cycle

        if (irank1==irank2) then

          QB_4_P(:,:,:,:) = QB_3_M(:,:,:,:)
          iflag(4) = 1
          QB_3_P(:,:,:,:) = QB_4_M(:,:,:,:)
          iflag(3) = 1

        else

          if (irank1==myrank_domain) then
#ifdef _USE_MPI
            call MPI_Send(QB_4_M,sendrecvnum,MPI_DOUBLE_PRECISION,irank2,0, &
                          MPI_COMM_DOMAIN,ierror)
#endif
#ifdef _USE_MPI
            call MPI_Recv(QB_4_P,sendrecvnum,MPI_DOUBLE_PRECISION,irank2,0, &
                          MPI_COMM_DOMAIN,istatus,ierror)
#endif
            iflag(4) = 1

          endif

          if (irank2==myrank_domain) then
#ifdef _USE_MPI
            call MPI_Recv(QB_3_P,sendrecvnum,MPI_DOUBLE_PRECISION,irank1,0, &
                          MPI_COMM_DOMAIN,istatus,ierror)
#endif
            iflag(3) = 1
#ifdef _USE_MPI
            call MPI_Send(QB_3_M,sendrecvnum,MPI_DOUBLE_PRECISION,irank1,0, &
                          MPI_COMM_DOMAIN,ierror)
#endif
          endif

        endif

      enddo
    enddo

  end subroutine mpisub_sbsp_exchange_domain_halo_eta
!-------------------------------------------------------------------------------
  subroutine mpisub_sbsp_exchange_domain_halo_zta(QB_5_M,QB_5_P, &
                                                  QB_6_M,QB_6_P,lmax,depth,iflag,mode)
    real(8), intent(inout), dimension(:,:,:,:) :: QB_5_M,QB_6_M
    real(8), intent(inout), dimension(:,:,:,:) :: QB_5_P,QB_6_P
    integer, intent(in) :: lmax
    integer, intent(in) :: depth
    integer, intent(out), dimension(:) :: iflag
    integer, intent(in) :: mode !0:cell 1:node

    integer :: imax,jmax,kmax
    integer :: isplt_blk,jsplt_blk,ksplt_blk,is_blk,js_blk
    integer :: id_sbdm1,id_blk1,irank1
    integer :: id_sbdm2,id_blk2,irank2
    integer :: ierror
#ifdef _USE_MPI
    integer, dimension(MPI_STATUS_SIZE) :: istatus
#endif

    integer :: sendrecvnum

    call Get_DomainSplit_From_DomainID(id_domain,isplt_blk,jsplt_blk,ksplt_blk)

    call Get_BlockIJKNum_From_BlockID(id_block,imax,jmax,kmax)
    if (mode==1) then
      imax = imax+1
      jmax = jmax+1
      kmax = kmax+1
    endif

    sendrecvnum = imax*jmax*depth*lmax

    !k direction
    do js_blk=1,jsplt_blk
      do is_blk=1,isplt_blk

        call Get_SubDomainID_From_DomainID_And_SplitID(id_domain,is_blk,js_blk,ksplt_blk,id_sbdm1)
        call Get_SubDomainID_From_DomainID_And_SplitID(id_domain,is_blk,js_blk,1        ,id_sbdm2)
        call Get_BlockID_From_DomainID_And_SubDomainID(id_domain,id_sbdm1,id_blk1)
        call Get_BlockID_From_DomainID_And_SubDomainID(id_domain,id_sbdm2,id_blk2)

        irank1 = BlockID_To_Rank_Domain(id_blk1)
        irank2 = BlockID_To_Rank_Domain(id_blk2)

        if (.not.(irank1==myrank_domain.or.irank2==myrank_domain)) cycle

        if (irank1==irank2) then

          QB_6_P(:,:,:,:) = QB_5_M(:,:,:,:)
          iflag(6) = 1
          QB_5_P(:,:,:,:) = QB_6_M(:,:,:,:)
          iflag(5) = 1

        else

          if (irank1==myrank_domain) then

#ifdef _USE_MPI
            call MPI_Send(QB_6_M,sendrecvnum,MPI_DOUBLE_PRECISION,irank2,0, &
                          MPI_COMM_DOMAIN,ierror)
#endif
#ifdef _USE_MPI
            call MPI_Recv(QB_6_P,sendrecvnum,MPI_DOUBLE_PRECISION,irank2,0, &
                          MPI_COMM_DOMAIN,istatus,ierror)
#endif
            iflag(6) = 1

          endif

          if (irank2==myrank_domain) then

#ifdef _USE_MPI
            call MPI_Recv(QB_5_P,sendrecvnum,MPI_DOUBLE_PRECISION,irank1,0, &
                          MPI_COMM_DOMAIN,istatus,ierror)
#endif
            iflag(5) = 1
#ifdef _USE_MPI
            call MPI_Send(QB_5_M,sendrecvnum,MPI_DOUBLE_PRECISION,irank1,0, &
                          MPI_COMM_DOMAIN,ierror)
#endif
          endif

        endif

      enddo
    enddo

  end subroutine mpisub_sbsp_exchange_domain_halo_zta
!-------------------------------------------------------------------------------
  subroutine mpisub_sbsp_output_info
    integer :: ierror
    integer :: id_dm,id_sbdm,id_blk,nblocks
    integer :: is,js,ks
    integer :: irank_world,irank_domain
    integer :: irank_domainplane_i,irank_domainplane_j,irank_domainplane_k

    call Get_BlockNum(nblocks)

    if (myrank_world==0) then

      write(6,*) '/-----------------------------------------------------------/'
      write(6,*) '/                mpisub_sbsp_output_info                    /'
      write(6,*) '/-----------------------------------------------------------/'

      do id_blk=1,nblocks
        call Get_DomainID_And_SubDomainID_From_BlockID(id_blk,id_dm,id_sbdm)
        call Get_SplitID_From_DomainID_And_SubDomainID(id_dm,id_sbdm,is,js,ks)
        call Get_Rank_From_BlockID(id_blk,irank_world)
        irank_domain = BlockID_To_Rank_Domain(id_blk)
        irank_domainplane_i = BlockID_To_Rank_DomainPlane_I(id_blk)
        irank_domainplane_j = BlockID_To_Rank_DomainPlane_J(id_blk)
        irank_domainplane_k = BlockID_To_Rank_DomainPlane_K(id_blk)

        write(6,'(A21, I5)') 'BlockID            : ', id_blk
        write(6,'(A21, I5)') 'DomainID           : ', id_dm
        write(6,'(A21, I5)') 'SubDomainID        : ', id_sbdm
        write(6,'(A21,3I5)') 'is, js, ks         : ', is,js,ks
        write(6,'(A21, I5)') 'Rank_World         : ', irank_world
        write(6,'(A21, I5)') 'Rank_Domain        : ', irank_domain
        write(6,'(A21, I5)') 'Rank_DomainPlane_I : ', irank_domainplane_i
        write(6,'(A21, I5)') 'Rank_DomainPlane_J : ', irank_domainplane_j
        write(6,'(A21, I5)') 'Rank_DomainPlane_K : ', irank_domainplane_k
        write(6,*) '-----------------------------------------------------------'
      enddo

      call flush(6)

    endif

#ifdef _USE_MPI
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
#endif

  end subroutine mpisub_sbsp_output_info

end module mdl_mpisub_sbsp
