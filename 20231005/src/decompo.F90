module mdl_decompo
#ifdef _USE_MPI
  use mpi
#endif
  use mdl_input_manager
  implicit none

  !input------------------------------------------------------------------------
  character(100), private, save :: fname
  !-----------------------------------------------------------------------------

  integer       , private, parameter :: nhalo = 0 !Should be 0 in Halo implementation

  integer, private, save :: nblocks
  integer, private, save :: ndomain

  integer, private, save, allocatable, dimension(:)   :: Rank_To_BlockID

  integer, private, save, allocatable, dimension(:)   :: BlockID_To_Rank

  integer, private, save, allocatable, dimension(:,:) :: BlockID_To_BlockIJKNum
  integer, private, save, allocatable, dimension(:,:) :: BlockID_To_BlockIJKNum_NoHalo

  integer, private, save, allocatable, dimension(:)   :: BlockID_To_DomainID
  integer, private, save, allocatable, dimension(:)   :: BlockID_To_SubDomainID
  integer, private, save, allocatable, dimension(:,:) :: BlockID_To_SubDomainIJKStartEnd
  integer, private, save, allocatable, dimension(:,:) :: BlockID_To_SubDomainIJKStartEnd_NoHalo

  integer, private, save, allocatable, dimension(:,:) :: DomainID_To_DomainIJKNum
  integer, private, save, allocatable, dimension(:,:) :: DomainID_To_DomainSplit

  integer, private, save, allocatable, dimension(:)   :: DomainID_To_Offset_Of_SubDomainID_To_BlockID
  integer, private, save, allocatable, dimension(:)   :: SubDomainID_To_BlockID

  integer, private, save, allocatable, dimension(:,:) :: DomainID_To_DomainPeriodicFlag

contains

  subroutine Get_HaloNum(nhl)
    integer, intent(out) :: nhl
    nhl = nhalo
  end subroutine Get_HaloNum

  subroutine Get_BlockNum(nblk)
    integer, intent(out) :: nblk
    nblk = nblocks
  end subroutine Get_BlockNum

  subroutine Get_DomainNum(ndmn)
    integer, intent(out) :: ndmn
    ndmn = ndomain
  end subroutine Get_DomainNum

  subroutine Get_DomainPeriodicFlag_From_DomainID(id_dm,fxsi,feta,fzta)
    integer, intent(in) :: id_dm
    integer, intent(out):: fxsi,feta,fzta
    fxsi = DomainID_To_DomainPeriodicFlag(1,id_dm)
    feta = DomainID_To_DomainPeriodicFlag(2,id_dm)
    fzta = DomainID_To_DomainPeriodicFlag(3,id_dm)
  end subroutine Get_DomainPeriodicFlag_From_DomainID

  subroutine Get_DomainIJKNum_From_DomainID(id_dm,inum,jnum,knum)
    integer, intent(in) :: id_dm
    integer, intent(out):: inum,jnum,knum
    inum = DomainID_To_DomainIJKNum(1,id_dm)
    jnum = DomainID_To_DomainIJKNum(2,id_dm)
    knum = DomainID_To_DomainIJKNum(3,id_dm)
  end subroutine Get_DomainIJKNum_From_DomainID

  subroutine Get_DomainSplit_From_DomainID(id_dm,isplt,jsplt,ksplt)
    integer, intent(in) :: id_dm
    integer, intent(out):: isplt,jsplt,ksplt
    isplt = DomainID_To_DomainSplit(1,id_dm)
    jsplt = DomainID_To_DomainSplit(2,id_dm)
    ksplt = DomainID_To_DomainSplit(3,id_dm)
  end subroutine Get_DomainSplit_From_DomainID

  subroutine Get_SubDomainNum_From_DomainID(id_dm,nsbdm)
    integer, intent(in) :: id_dm
    integer, intent(out):: nsbdm
    integer :: isplt,jsplt,ksplt
    isplt = DomainID_To_DomainSplit(1,id_dm)
    jsplt = DomainID_To_DomainSplit(2,id_dm)
    ksplt = DomainID_To_DomainSplit(3,id_dm)
    nsbdm = isplt*jsplt*ksplt
  end subroutine Get_SubDomainNum_From_DomainID

  subroutine Get_SubDomainID_From_DomainID_And_SplitID(id_dm,is,js,ks,id_sbdm)
    integer, intent(in) :: id_dm,is,js,ks
    integer, intent(out):: id_sbdm
    integer :: isplt,jsplt,ksplt
    isplt = DomainID_To_DomainSplit(1,id_dm)
    jsplt = DomainID_To_DomainSplit(2,id_dm)
    ksplt = DomainID_To_DomainSplit(3,id_dm)
    id_sbdm = is+isplt*(js-1)+isplt*jsplt*(ks-1)
  end subroutine Get_SubDomainID_From_DomainID_And_SplitID

  subroutine Get_SplitID_From_DomainID_And_SubDomainID(id_dm,id_sbdm,is,js,ks)
    integer, intent(in)    :: id_dm,id_sbdm
    integer, intent(inout) :: is,js,ks
    integer :: isplt,jsplt,ksplt,id_tmp
    call Get_DomainSplit_From_DomainID(id_dm,isplt,jsplt,ksplt)
    ks = (id_sbdm-1)/(isplt*jsplt)+1
    id_tmp = id_sbdm-isplt*jsplt*(ks-1)
    js =  (id_tmp-1)/isplt+1
    is = id_tmp-isplt*(js-1)
  end subroutine Get_SplitID_From_DomainID_And_SubDomainID

  subroutine Get_BlockID_From_Rank(irank,id_blk)
    integer, intent(in) :: irank
    integer, intent(out):: id_blk
    id_blk = Rank_To_BlockID(irank)
  end subroutine Get_BlockID_From_Rank

  subroutine Get_BlockID_From_DomainID_And_SubDomainID(id_dm,id_sbdm,id_blk)
    integer, intent(in) :: id_dm,id_sbdm
    integer, intent(out):: id_blk
    integer :: offset
    offset = DomainID_To_Offset_Of_SubDomainID_To_BlockID(id_dm)
    id_blk = SubDomainID_To_BlockID(offset+id_sbdm)
  end subroutine Get_BlockID_From_DomainID_And_SubDomainID

  subroutine Get_BlockIJKNum_From_BlockID(id_blk,inum,jnum,knum)
    integer, intent(in) :: id_blk
    integer, intent(out):: inum,jnum,knum
    inum = BlockID_To_BlockIJKNum(1,id_blk)
    jnum = BlockID_To_BlockIJKNum(2,id_blk)
    knum = BlockID_To_BlockIJKNum(3,id_blk)
  end subroutine Get_BlockIJKNum_From_BlockID

  subroutine Get_BlockIJKNum_NoHalo_From_BlockID(id_blk,inum,jnum,knum)
    integer, intent(in) :: id_blk
    integer, intent(out):: inum,jnum,knum
    inum = BlockID_To_BlockIJKNum_NoHalo(1,id_blk)
    jnum = BlockID_To_BlockIJKNum_NoHalo(2,id_blk)
    knum = BlockID_To_BlockIJKNum_NoHalo(3,id_blk)
  end subroutine Get_BlockIJKNum_NoHalo_From_BlockID

  subroutine Get_Rank_From_BlockID(id_blk,irank)
    integer, intent(in) :: id_blk
    integer, intent(out):: irank
    irank = BlockID_To_Rank(id_blk)
  end subroutine Get_Rank_From_BlockID

  subroutine Get_DomainID_And_SubDomainID_From_BlockID(id_blk,id_dm,id_sbdm)
    integer, intent(in) :: id_blk
    integer, intent(out):: id_dm,id_sbdm
    id_dm = BlockID_To_DomainID(id_blk)
    id_sbdm = BlockID_To_SubDomainID(id_blk)
  end subroutine Get_DomainID_And_SubDomainID_From_BlockID

  subroutine Get_SubDomainIJKStartEnd_From_DomainID_And_SubDomainID(id_dm,id_sbdm, &
                                                istart,iend,jstart,jend,kstart,kend)
    integer, intent(in) :: id_dm,id_sbdm
    integer, intent(out):: istart,iend,jstart,jend,kstart,kend
    integer :: id_blk
    call Get_BlockID_From_DomainID_And_SubDomainID(id_dm,id_sbdm,id_blk)
    istart = BlockID_To_SubDomainIJKStartEnd(1,id_blk)
    jstart = BlockID_To_SubDomainIJKStartEnd(2,id_blk)
    kstart = BlockID_To_SubDomainIJKStartEnd(3,id_blk)
    iend   = BlockID_To_SubDomainIJKStartEnd(4,id_blk)
    jend   = BlockID_To_SubDomainIJKStartEnd(5,id_blk)
    kend   = BlockID_To_SubDomainIJKStartEnd(6,id_blk)
  end subroutine Get_SubDomainIJKStartEnd_From_DomainID_And_SubDomainID

  subroutine Get_SubDomainIJKStartEnd_NoHalo_From_DomainID_And_SubDomainID(id_dm,id_sbdm, &
                                                istart,iend,jstart,jend,kstart,kend)
    integer, intent(in) :: id_dm,id_sbdm
    integer, intent(out):: istart,iend,jstart,jend,kstart,kend
    integer :: id_blk
    call Get_BlockID_From_DomainID_And_SubDomainID(id_dm,id_sbdm,id_blk)
    istart = BlockID_To_SubDomainIJKStartEnd_NoHalo(1,id_blk)
    jstart = BlockID_To_SubDomainIJKStartEnd_NoHalo(2,id_blk)
    kstart = BlockID_To_SubDomainIJKStartEnd_NoHalo(3,id_blk)
    iend   = BlockID_To_SubDomainIJKStartEnd_NoHalo(4,id_blk)
    jend   = BlockID_To_SubDomainIJKStartEnd_NoHalo(5,id_blk)
    kend   = BlockID_To_SubDomainIJKStartEnd_NoHalo(6,id_blk)
  end subroutine Get_SubDomainIJKStartEnd_NoHalo_From_DomainID_And_SubDomainID

  subroutine Get_Rank_From_DomainID_And_SubdomainID(id_dm,id_sbdm,irank)
    integer, intent(in) :: id_dm,id_sbdm
    integer, intent(out):: irank
    integer :: id_blk
    call Get_BlockID_From_DomainID_And_SubDomainID(id_dm,id_sbdm,id_blk)
    call Get_Rank_From_BlockID(id_blk,irank)
  end subroutine Get_Rank_From_DomainID_And_SubdomainID

  subroutine decompo_init
    integer :: myrank,nprocs,iproc,ierror
    integer :: i,j,l,ii,jj,kk,inum,jnum,knum,icnt,imax,jmax,kmax
    integer :: id_dm
    integer :: isplt,jsplt,ksplt
    integer, dimension(:,:), allocatable :: DomainIJKNumNode
    integer, dimension(:,:), allocatable :: DomainIJKNum
    integer, dimension(:,:), allocatable :: DomainSplit
    integer, dimension(:,:), allocatable :: DomainPeriodic
    integer, dimension(:,:,:,:), allocatable :: ijk_max1,ijk_max2
    integer, dimension(:,:,:,:), allocatable :: ijk_sta1,ijk_sta2
    integer, dimension(:,:,:,:), allocatable :: ijk_end1,ijk_end2

    call input_decompo(fname)

    !Read input file---------------
    myrank = 0
    nprocs = 1

#ifdef _USE_MPI
    call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierror)
    call MPI_Comm_size(MPI_COMM_WORLD,nprocs,ierror)
#endif

    do iproc=0,nprocs-1
      if (myrank==iproc) then
        open(28,file=fname,form='formatted')
        read(28,*) ndomain
        allocate(DomainIJKNumNode(3,ndomain))
        allocate(DomainIJKNum(3,ndomain))
        allocate(DomainSplit(3,ndomain))
        allocate(DomainPeriodic(3,ndomain))
        do i=1,ndomain
          read(28,*) id_dm, &
                     DomainIJKNumNode(1,i),DomainIJKNumNode(2,i),DomainIJKNumNode(3,i), &
                     DomainSplit(1,i)     ,DomainSplit(2,i)     ,DomainSplit(3,i)     , &
                     DomainPeriodic(1,i)  ,DomainPeriodic(2,i)  ,DomainPeriodic(3,i)
          DomainIJKNum(1,i) = max0(DomainIJKNumNode(1,i)-1,1)
          DomainIJKNum(2,i) = max0(DomainIJKNumNode(2,i)-1,1)
          DomainIJKNum(3,i) = max0(DomainIJKNumNode(3,i)-1,1)
        enddo
        close(28)
      endif

#ifdef _USE_MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
#endif

    enddo
    !------------------------------

    allocate(DomainID_To_DomainIJKNum(3,ndomain))
    allocate(DomainID_To_DomainSplit(3,ndomain))
    allocate(DomainID_To_Offset_Of_SubDomainID_To_BlockID(ndomain))

    allocate(DomainID_To_DomainPeriodicFlag(3,ndomain))

    do j=1,ndomain
      do i=1,3
        DomainID_To_DomainIJKNum(i,j) = DomainIJKNum(i,j)
        DomainID_To_DomainSplit(i,j) = DomainSplit(i,j)
        DomainID_To_DomainPeriodicFlag(i,j) = DomainPeriodic(i,j)
      enddo
    enddo

    nblocks = 0
    do j=1,ndomain
      isplt = DomainID_To_DomainSplit(1,j)
      jsplt = DomainID_To_DomainSplit(2,j)
      ksplt = DomainID_To_DomainSplit(3,j)
      nblocks = nblocks + isplt*jsplt*ksplt
    enddo

    allocate(Rank_To_BlockID(0:nblocks-1))

    allocate(BlockID_To_Rank(nblocks))
    allocate(BlockID_To_BlockIJKNum(3,nblocks))
    allocate(BlockID_To_BlockIJKNum_NoHalo(3,nblocks))
    allocate(BlockID_To_DomainID(nblocks))
    allocate(BlockID_To_SubDomainID(nblocks))
    allocate(SubDomainID_To_BlockID(nblocks))
    allocate(BlockID_To_SubDomainIJKStartEnd(6,nblocks))
    allocate(BlockID_To_SubDomainIJKStartEnd_NoHalo(6,nblocks))

    icnt = 1
    do j=1,ndomain
      isplt = DomainID_To_DomainSplit(1,j)
      jsplt = DomainID_To_DomainSplit(2,j)
      ksplt = DomainID_To_DomainSplit(3,j)
      do i=1,isplt*jsplt*ksplt
        BlockID_To_DomainID(icnt) = j
        BlockID_To_SubDomainID(icnt) = i
        SubDomainID_To_BlockID(icnt) = icnt
        icnt = icnt+1
      enddo
      if (j.eq.1) then
        DomainID_To_Offset_Of_SubDomainID_To_BlockID(j) = 0
      else
        DomainID_To_Offset_Of_SubDomainID_To_BlockID(j) &
         = DomainID_To_Offset_Of_SubDomainID_To_BlockID(j-1) &
         + DomainID_To_DomainSplit(1,j-1) &
         * DomainID_To_DomainSplit(2,j-1) &
         * DomainID_To_DomainSplit(3,j-1)
      endif
    enddo

    icnt = 1
    do j=1,ndomain
      isplt = DomainID_To_DomainSplit(1,j)
      jsplt = DomainID_To_DomainSplit(2,j)
      ksplt = DomainID_To_DomainSplit(3,j)
      imax = DomainID_To_DomainIJKNum(1,j)
      jmax = DomainID_To_DomainIJKNum(2,j)
      kmax = DomainID_To_DomainIJKNum(3,j)
      inum = imax
      jnum = jmax
      knum = kmax

      allocate(ijk_max1(3,isplt,jsplt,ksplt))
      allocate(ijk_sta1(3,isplt,jsplt,ksplt))
      allocate(ijk_end1(3,isplt,jsplt,ksplt))
      allocate(ijk_max2(3,isplt,jsplt,ksplt))
      allocate(ijk_sta2(3,isplt,jsplt,ksplt))
      allocate(ijk_end2(3,isplt,jsplt,ksplt))
      ijk_max1 = 0
      ijk_sta1 = 0
      ijk_end1 = 0
      ijk_max2 = 0
      ijk_sta2 = 0
      ijk_end2 = 0

      do l=1,inum
        ii = mod(l,isplt)
        if (ii.eq.0) ii = isplt
        ijk_max1(1,ii,1,1) = ijk_max1(1,ii,1,1) + 1
      enddo

      do l=1,jnum
        jj = mod(l,jsplt)
        if (jj.eq.0) jj = jsplt
        ijk_max1(2,1,jj,1) = ijk_max1(2,1,jj,1) + 1
      enddo

      do l=1,knum
        kk = mod(l,ksplt)
        if (kk.eq.0) kk = ksplt
        ijk_max1(3,1,1,kk) = ijk_max1(3,1,1,kk) + 1
      enddo

      do kk=1,ksplt
        do jj=1,jsplt
          do ii=1,isplt
            ijk_max1(1,ii,jj,kk) = ijk_max1(1,ii,1,1)
            ijk_max1(2,ii,jj,kk) = ijk_max1(2,1,jj,1)
            ijk_max1(3,ii,jj,kk) = ijk_max1(3,1,1,kk)
          enddo
        enddo
      enddo

      do ii=1,isplt
        if (ii==1) then
          ijk_sta1(1,ii,1,1) = 1
        else
          ijk_sta1(1,ii,1,1) = ijk_end1(1,ii-1,1,1)+1
        endif
        ijk_end1(1,ii,1,1) = ijk_sta1(1,ii,1,1)+ijk_max1(1,ii,1,1)-1
      enddo
      do jj=1,jsplt
        if (jj==1) then
          ijk_sta1(2,1,jj,1) = 1
        else
          ijk_sta1(2,1,jj,1) = ijk_end1(2,1,jj-1,1)+1
        endif
        ijk_end1(2,1,jj,1) = ijk_sta1(2,1,jj,1)+ijk_max1(2,1,jj,1)-1
      enddo
      do kk=1,ksplt
        if (kk==1) then
          ijk_sta1(3,1,1,kk) = 1
        else
          ijk_sta1(3,1,1,kk) = ijk_end1(3,1,1,kk-1)+1
        endif
        ijk_end1(3,1,1,kk) = ijk_sta1(3,1,1,kk)+ijk_max1(3,1,1,kk)-1
      enddo
      do kk=1,ksplt
        do jj=1,jsplt
          do ii=1,isplt
            ijk_sta1(1,ii,jj,kk) = ijk_sta1(1,ii,1,1)
            ijk_end1(1,ii,jj,kk) = ijk_end1(1,ii,1,1)
            ijk_sta1(2,ii,jj,kk) = ijk_sta1(2,1,jj,1)
            ijk_end1(2,ii,jj,kk) = ijk_end1(2,1,jj,1)
            ijk_sta1(3,ii,jj,kk) = ijk_sta1(3,1,1,kk)
            ijk_end1(3,ii,jj,kk) = ijk_end1(3,1,1,kk)
          enddo
        enddo
      enddo

      ijk_max2 = ijk_max1
      ijk_sta2 = ijk_sta1
      ijk_end2 = ijk_end1

      do kk=1,ksplt
        do jj=1,jsplt
          do ii=1,isplt
            if (ii/=1    ) ijk_max2(1,ii,jj,kk) = ijk_max2(1,ii,jj,kk)+nhalo
            if (ii/=isplt) ijk_max2(1,ii,jj,kk) = ijk_max2(1,ii,jj,kk)+nhalo
            if (jj/=1    ) ijk_max2(2,ii,jj,kk) = ijk_max2(2,ii,jj,kk)+nhalo
            if (jj/=jsplt) ijk_max2(2,ii,jj,kk) = ijk_max2(2,ii,jj,kk)+nhalo
            if (kk/=1    ) ijk_max2(3,ii,jj,kk) = ijk_max2(3,ii,jj,kk)+nhalo
            if (kk/=ksplt) ijk_max2(3,ii,jj,kk) = ijk_max2(3,ii,jj,kk)+nhalo
          enddo
        enddo
      enddo

      do kk=1,ksplt
        do jj=1,jsplt
          do ii=1,isplt
            if (ii/=1    ) ijk_sta2(1,ii,jj,kk) = ijk_sta2(1,ii,jj,kk)-nhalo
            if (jj/=1    ) ijk_sta2(2,ii,jj,kk) = ijk_sta2(2,ii,jj,kk)-nhalo
            if (kk/=1    ) ijk_sta2(3,ii,jj,kk) = ijk_sta2(3,ii,jj,kk)-nhalo
            if (ii/=isplt) ijk_end2(1,ii,jj,kk) = ijk_end2(1,ii,jj,kk)+nhalo
            if (jj/=jsplt) ijk_end2(2,ii,jj,kk) = ijk_end2(2,ii,jj,kk)+nhalo
            if (kk/=ksplt) ijk_end2(3,ii,jj,kk) = ijk_end2(3,ii,jj,kk)+nhalo
          enddo
        enddo
      enddo

      do kk=1,ksplt
        do jj=1,jsplt
          do ii=1,isplt

            BlockID_To_BlockIJKNum_NoHalo(1,icnt) = ijk_max1(1,ii,jj,kk)
            BlockID_To_BlockIJKNum_NoHalo(2,icnt) = ijk_max1(2,ii,jj,kk)
            BlockID_To_BlockIJKNum_NoHalo(3,icnt) = ijk_max1(3,ii,jj,kk)

            BlockID_To_SubDomainIJKStartEnd_NoHalo(1,icnt) = ijk_sta1(1,ii,jj,kk)
            BlockID_To_SubDomainIJKStartEnd_NoHalo(2,icnt) = ijk_sta1(2,ii,jj,kk)
            BlockID_To_SubDomainIJKStartEnd_NoHalo(3,icnt) = ijk_sta1(3,ii,jj,kk)
            BlockID_To_SubDomainIJKStartEnd_NoHalo(4,icnt) = ijk_end1(1,ii,jj,kk)
            BlockID_To_SubDomainIJKStartEnd_NoHalo(5,icnt) = ijk_end1(2,ii,jj,kk)
            BlockID_To_SubDomainIJKStartEnd_NoHalo(6,icnt) = ijk_end1(3,ii,jj,kk)

            BlockID_To_BlockIJKNum(1,icnt) = ijk_max2(1,ii,jj,kk)
            BlockID_To_BlockIJKNum(2,icnt) = ijk_max2(2,ii,jj,kk)
            BlockID_To_BlockIJKNum(3,icnt) = ijk_max2(3,ii,jj,kk)

            BlockID_To_SubDomainIJKStartEnd(1,icnt) = ijk_sta2(1,ii,jj,kk)
            BlockID_To_SubDomainIJKStartEnd(2,icnt) = ijk_sta2(2,ii,jj,kk)
            BlockID_To_SubDomainIJKStartEnd(3,icnt) = ijk_sta2(3,ii,jj,kk)
            BlockID_To_SubDomainIJKStartEnd(4,icnt) = ijk_end2(1,ii,jj,kk)
            BlockID_To_SubDomainIJKStartEnd(5,icnt) = ijk_end2(2,ii,jj,kk)
            BlockID_To_SubDomainIJKStartEnd(6,icnt) = ijk_end2(3,ii,jj,kk)

            icnt = icnt+1
          enddo
        enddo
      enddo

      deallocate(ijk_max1)
      deallocate(ijk_sta1)
      deallocate(ijk_end1)
      deallocate(ijk_max2)
      deallocate(ijk_sta2)
      deallocate(ijk_end2)
    enddo

    do i=1,nblocks
      Rank_To_BlockID(i-1) = i
      BlockID_To_Rank(i  ) = i-1
    enddo

    do iproc=0,nprocs-1
      if (myrank==iproc) then
        if (nprocs/=nblocks) then
          write(6,*) 'The number of processes is invalid.'
          write(6,*) nblocks, ' processes are required. Good luck !!'
          call flush(6)
        endif
      endif

#ifdef _USE_MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
#endif

    enddo

    call decompo_output_info

    deallocate(DomainIJKNumNode)
    deallocate(DomainIJKNum)
    deallocate(DomainSplit)
    deallocate(DomainPeriodic)
  end subroutine decompo_init

  subroutine decompo_final
    deallocate(Rank_To_BlockID)
    deallocate(BlockID_To_Rank)
    deallocate(BlockID_To_BlockIJKNum)
    deallocate(BlockID_To_BlockIJKNum_NoHalo)
    deallocate(BlockID_To_DomainID)
    deallocate(BlockID_To_SubDomainID)
    deallocate(DomainID_To_DomainIJKNum)
    deallocate(DomainID_To_DomainSplit)
    deallocate(DomainID_To_Offset_Of_SubDomainID_To_BlockID)
    deallocate(SubDomainID_To_BlockID)
    deallocate(BlockID_To_SubDomainIJKStartEnd)
    deallocate(BlockID_To_SubDomainIJKStartEnd_NoHalo)
    deallocate(DomainID_To_DomainPeriodicFlag)
  end subroutine decompo_final

  subroutine decompo_output_info

    integer :: myrank,ierror
    integer :: id_dm,id_blk,id_sbdm
    integer :: isplt,jsplt,ksplt
    integer :: offset

    800 format(A9,6A7)
    801 format(I9,6I7)
    802 format(A9,4A7)
    803 format(I9,4I7)
    804 format(A12,A8,6A7)
    805 format(I12,I8,6I7)

    myrank = 0

#ifdef _USE_MPI
    call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierror)
#endif

    if (myrank==0) then

      write(6,*) '/-----------------------------------------------------------/'
      write(6,*) '/                  decompo_output_info                      /'
      write(6,*) '/-----------------------------------------------------------/'

      write(6,*) 'ndomain = ', ndomain
      write(6,800) 'DomainID','Inum','Jnum','Knum','Isplit','Jsplit','Ksplit'
      do id_dm=1,ndomain
        write(6,801) id_dm, &
                     DomainID_To_DomainIJKNum(1,id_dm), &
                     DomainID_To_DomainIJKNum(2,id_dm), &
                     DomainID_To_DomainIJKNum(3,id_dm), &
                     DomainID_To_DomainSplit(1,id_dm),  &
                     DomainID_To_DomainSplit(2,id_dm),  &
                     DomainID_To_DomainSplit(3,id_dm)
      enddo

      write(6,*) '-------------------------------------------------------------'

      write(6,*) 'nblocks = ', nblocks
      write(6,802) 'BlockID','Inum','Jnum','Knum','Rank'
      do id_blk = 1,nblocks
        write(6,803) id_blk, &
                     BlockID_To_BlockIJKNum(1,id_blk), &
                     BlockID_To_BlockIJKNum(2,id_blk), &
                     BlockID_To_BlockIJKNum(3,id_blk), &
                     BlockID_To_Rank(id_blk)
      enddo

      write(6,*) '-------------------------------------------------------------'

      do id_dm=1,ndomain
        isplt = DomainID_To_DomainSplit(1,id_dm)
        jsplt = DomainID_To_DomainSplit(2,id_dm)
        ksplt = DomainID_To_DomainSplit(3,id_dm)
        offset= DomainID_To_Offset_Of_SubDomainID_To_BlockID(id_dm)
        write(6,*) 'DomainID = ',id_dm
        write(6,804) 'SubDomainID','BlockID','Istart','Iend','Jstart','Jend','Kstart','Kend'
        do id_sbdm=1,isplt*jsplt*ksplt
          id_blk = SubDomainID_To_BlockID(offset+id_sbdm)
          write(6,805) id_sbdm, &
                        id_blk , &
                        BlockID_To_SubDomainIJKStartEnd(1,id_blk), &
                        BlockID_To_SubDomainIJKStartEnd(4,id_blk), &
                        BlockID_To_SubDomainIJKStartEnd(2,id_blk), &
                        BlockID_To_SubDomainIJKStartEnd(5,id_blk), &
                        BlockID_To_SubDomainIJKStartEnd(3,id_blk), &
                        BlockID_To_SubDomainIJKStartEnd(6,id_blk)
        enddo
      enddo

      call flush(6)

    endif

#ifdef _USE_MPI
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
#endif

  end subroutine decompo_output_info

end module mdl_decompo
