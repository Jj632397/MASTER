module mdl_exchange_domain_halo
#ifdef _USE_MPI
  use mpi
#endif
  use mdl_input_manager
  use mdl_decompo
  use mdl_mpisub_sbsp
  implicit none

  !input------------------------------------------------------------------------
  character(100), private, save :: fname
  !-----------------------------------------------------------------------------

  integer, private, save :: nexchange

  integer, private, save, allocatable, dimension(:) :: id_dm_1
  integer, private, save, allocatable, dimension(:) :: id_bd_1
  integer, private, save, allocatable, dimension(:) :: id_ax1_1
  integer, private, save, allocatable, dimension(:) :: id_ax2_1
  !------------------------------------------------------------
  integer, private, save, allocatable, dimension(:) :: id_dm_2
  integer, private, save, allocatable, dimension(:) :: id_bd_2
  integer, private, save, allocatable, dimension(:) :: id_ax1_2
  integer, private, save, allocatable, dimension(:) :: id_ax2_2


contains

  subroutine exchange_domain_halo_init
    integer :: myrank,nprocs,iproc
    integer :: i,ierror,index

    call input_exchange_domain_halo(fname)

#ifdef _USE_MPI
    call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierror)
    call MPI_Comm_size(MPI_COMM_WORLD,nprocs,ierror)
#else
    myrank = 0
    nprocs = 1
#endif

    do iproc=0,nprocs-1
      if (myrank==iproc) then
        open(41,file=fname,form='formatted')
        read(41,*) nexchange
        allocate(id_dm_1(nexchange))
        allocate(id_bd_1(nexchange))
        allocate(id_ax1_1(nexchange))
        allocate(id_ax2_1(nexchange))
        !----------------------------
        allocate(id_dm_2(nexchange))
        allocate(id_bd_2(nexchange))
        allocate(id_ax1_2(nexchange))
        allocate(id_ax2_2(nexchange))
        !----------------------------
        do i=1,nexchange
          read(41,*) index,id_dm_1(i),id_bd_1(i),id_ax1_1(i),id_ax2_1(i), &
                           id_dm_2(i),id_bd_2(i),id_ax1_2(i),id_ax2_2(i)
        enddo
        close(41)
      endif
#ifdef _USE_MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
#endif
    enddo

  end subroutine exchange_domain_halo_init

  subroutine exchange_domain_halo_final
    deallocate(id_dm_1)
    deallocate(id_bd_1)
    deallocate(id_ax1_1)
    deallocate(id_ax2_1)
    !------------------
    deallocate(id_dm_2)
    deallocate(id_bd_2)
    deallocate(id_ax1_2)
    deallocate(id_ax2_2)
  end subroutine exchange_domain_halo_final

  subroutine exchange_domain_halo(QB_1_M, QB_1_P, &
                                  QB_2_M, QB_2_P, &
                                  QB_3_M, QB_3_P, &
                                  QB_4_M, QB_4_P, &
                                  QB_5_M, QB_5_P, &
                                  QB_6_M, QB_6_P, lmax, depth, iflag, mode)
    real(8), intent(in) , dimension(:,:,:,:) :: QB_1_M,QB_2_M,QB_3_M
    real(8), intent(in) , dimension(:,:,:,:) :: QB_4_M,QB_5_M,QB_6_M
    real(8), intent(out), dimension(:,:,:,:) :: QB_1_P,QB_2_P,QB_3_P
    real(8), intent(out), dimension(:,:,:,:) :: QB_4_P,QB_5_P,QB_6_P
    integer, intent(in) :: lmax
    integer, intent(in) :: depth
    integer, intent(out), dimension(:) :: iflag
    integer, intent(in) :: mode !0:cell 1:node

    integer :: myrank,id_block,id_domain,id_subdomain
    integer :: n
    integer :: imax_1,jmax_1,kmax_1
    integer :: imax_2,jmax_2,kmax_2
    integer :: iimax_1,jjmax_1
    integer :: iimax_2,jjmax_2
    integer :: is,js,ks
    integer :: isplt_1,jsplt_1,ksplt_1
    integer :: isplt_2,jsplt_2,ksplt_2
    integer :: is_plane_1,is_plane_2
    integer :: id_dm_send,id_dir_send,id_split_send
    integer :: id_dm_recv,id_dir_recv,id_split_recv
    integer :: i,j,k,l,m,ii,jj,kk,mm,iii,jjj,kkk
    integer :: istart,iend,jstart,jend,kstart,kend

    real(8), allocatable, dimension(:,:,:,:) :: QD_1
    real(8), allocatable, dimension(:,:,:,:) :: QD_2

    real(8), allocatable, dimension(:,:,:,:) :: QD_1_ALN
    real(8), allocatable, dimension(:,:,:,:) :: QD_2_ALN

    real(8), allocatable, dimension(:,:,:,:) :: QD_TMP

    do n=1,nexchange!-----------------------------------------------------------

      call Get_DomainIJKNum_From_DomainID(id_dm_1(n),imax_1,jmax_1,kmax_1)
      call Get_DomainIJKNum_From_DomainID(id_dm_2(n),imax_2,jmax_2,kmax_2)
      if (mode==1) then
        imax_1 = imax_1+1
        jmax_1 = jmax_1+1
        kmax_1 = kmax_1+1
        !----------------
        imax_2 = imax_2+1
        jmax_2 = jmax_2+1
        kmax_2 = kmax_2+1
      endif

      if (iabs(id_ax1_1(n))==1) then
        iimax_1 = imax_1
      else if (iabs(id_ax1_1(n))==2) then
        iimax_1 = jmax_1
      else if (iabs(id_ax1_1(n))==3) then
        iimax_1 = kmax_1
      endif

      if (iabs(id_ax2_1(n))==1) then
        jjmax_1 = imax_1
      else if (iabs(id_ax2_1(n))==2) then
        jjmax_1 = jmax_1
      else if (iabs(id_ax2_1(n))==3) then
        jjmax_1 = kmax_1
      endif

      if (iabs(id_ax1_2(n))==1) then
        iimax_2 = imax_2
      else if (iabs(id_ax1_2(n))==2) then
        iimax_2 = jmax_2
      else if (iabs(id_ax1_2(n))==3) then
        iimax_2 = kmax_2
      endif

      if (iabs(id_ax2_2(n))==1) then
        jjmax_2 = imax_2
      else if (iabs(id_ax2_2(n))==2) then
        jjmax_2 = jmax_2
      else if (iabs(id_ax2_2(n))==3) then
        jjmax_2 = kmax_2
      endif

      if (iimax_1/=iimax_2.or.jjmax_1/=jjmax_2) then
        write(6,*) 'Error:exchange_domain_halo:Invalid arguments 1'
      endif

      call Get_DomainSplit_From_DomainID(id_dm_1(n),isplt_1,jsplt_1,ksplt_1)
      call Get_DomainSplit_From_DomainID(id_dm_2(n),isplt_2,jsplt_2,ksplt_2)

      call mpisub_sbsp_get_myrank_world(myrank)

      call Get_BlockID_From_Rank(myrank,id_block)

      call Get_DomainID_And_SubDomainID_From_BlockID(id_block,id_domain,id_subdomain)

      call Get_SplitID_From_DomainID_And_SubDomainID(id_domain,id_subdomain,is,js,ks)

      is_plane_1 = 0
      is_plane_2 = 0

      if (id_bd_1(n)==1) then
        if (id_domain==id_dm_1(n).and.is==1      ) is_plane_1 = 1
      else if (id_bd_1(n)==2) then
        if (id_domain==id_dm_1(n).and.is==isplt_1) is_plane_1 = 1
      else if (id_bd_1(n)==3) then
        if (id_domain==id_dm_1(n).and.js==1      ) is_plane_1 = 1
      else if (id_bd_1(n)==4) then
        if (id_domain==id_dm_1(n).and.js==jsplt_1) is_plane_1 = 1
      else if (id_bd_1(n)==5) then
        if (id_domain==id_dm_1(n).and.ks==1      ) is_plane_1 = 1
      else if (id_bd_1(n)==6) then
        if (id_domain==id_dm_1(n).and.ks==ksplt_1) is_plane_1 = 1
      endif

      if (id_bd_2(n)==1) then
        if (id_domain==id_dm_2(n).and.is==1      ) is_plane_2 = 1
      else if (id_bd_2(n)==2) then
        if (id_domain==id_dm_2(n).and.is==isplt_2) is_plane_2 = 1
      else if (id_bd_2(n)==3) then
        if (id_domain==id_dm_2(n).and.js==1      ) is_plane_2 = 1
      else if (id_bd_2(n)==4) then
        if (id_domain==id_dm_2(n).and.js==jsplt_2) is_plane_2 = 1
      else if (id_bd_2(n)==5) then
        if (id_domain==id_dm_2(n).and.ks==1      ) is_plane_2 = 1
      else if (id_bd_2(n)==6) then
        if (id_domain==id_dm_2(n).and.ks==ksplt_2) is_plane_2 = 1
      endif

      if (.not.(is_plane_1==1.or.is_plane_2==1)) cycle

      if (is_plane_1==1) then
        if (id_bd_2(n)==1.or.id_bd_2(n)==2) then
          allocate(QD_2(depth,jmax_2,kmax_2,lmax))
        else if (id_bd_2(n)==3.or.id_bd_2(n)==4) then
          allocate(QD_2(imax_2,depth,kmax_2,lmax))
        else if (id_bd_2(n)==5.or.id_bd_2(n)==6) then
          allocate(QD_2(imax_2,jmax_2,depth,lmax))
        endif
      endif

      if (is_plane_2==1) then
        if (id_bd_1(n)==1.or.id_bd_1(n)==2) then
          allocate(QD_1(depth,jmax_1,kmax_1,lmax))
        else if (id_bd_1(n)==3.or.id_bd_1(n)==4) then
          allocate(QD_1(imax_1,depth,kmax_1,lmax))
        else if (id_bd_1(n)==5.or.id_bd_1(n)==6) then
          allocate(QD_1(imax_1,jmax_1,depth,lmax))
        endif
      endif

      !Send plane_1 to plane_2
      id_dm_send = id_dm_1(n)
      id_dm_recv = id_dm_2(n)
      if (id_bd_2(n)==1) then
        id_dir_recv = 1
        id_split_recv = 1
      else if (id_bd_2(n)==2) then
        id_dir_recv = 1
        id_split_recv = isplt_2
      else if (id_bd_2(n)==3) then
        id_dir_recv = 2
        id_split_recv = 1
      else if (id_bd_2(n)==4) then
        id_dir_recv = 2
        id_split_recv = jsplt_2
      else if (id_bd_2(n)==5) then
        id_dir_recv = 3
        id_split_recv = 1
      else if (id_bd_2(n)==6) then
        id_dir_recv = 3
        id_split_recv = ksplt_2
      endif
      if (id_bd_1(n)==1) then
        id_dir_send = 1
        id_split_send = 1
        call mpisub_sbsp_send_domainplane_to_domainplane( &
             id_dm_send,id_dir_send,id_split_send,QB_1_M, &
             id_dm_recv,id_dir_recv,id_split_recv,QD_1  ,depth,lmax,mode)
      else if (id_bd_1(n)==2) then
        id_dir_send = 1
        id_split_send = isplt_1
        call mpisub_sbsp_send_domainplane_to_domainplane( &
             id_dm_send,id_dir_send,id_split_send,QB_2_M, &
             id_dm_recv,id_dir_recv,id_split_recv,QD_1  ,depth,lmax,mode)
      else if (id_bd_1(n)==3) then
        id_dir_send = 2
        id_split_send = 1
        call mpisub_sbsp_send_domainplane_to_domainplane( &
             id_dm_send,id_dir_send,id_split_send,QB_3_M, &
             id_dm_recv,id_dir_recv,id_split_recv,QD_1  ,depth,lmax,mode)
      else if (id_bd_1(n)==4) then
        id_dir_send = 2
        id_split_send = jsplt_1
        call mpisub_sbsp_send_domainplane_to_domainplane( &
             id_dm_send,id_dir_send,id_split_send,QB_4_M, &
             id_dm_recv,id_dir_recv,id_split_recv,QD_1  ,depth,lmax,mode)
      else if (id_bd_1(n)==5) then
        id_dir_send = 3
        id_split_send = 1
        call mpisub_sbsp_send_domainplane_to_domainplane( &
             id_dm_send,id_dir_send,id_split_send,QB_5_M, &
             id_dm_recv,id_dir_recv,id_split_recv,QD_1  ,depth,lmax,mode)
      else if (id_bd_1(n)==6) then
        id_dir_send = 3
        id_split_send = ksplt_1
        call mpisub_sbsp_send_domainplane_to_domainplane( &
             id_dm_send,id_dir_send,id_split_send,QB_6_M, &
             id_dm_recv,id_dir_recv,id_split_recv,QD_1  ,depth,lmax,mode)
      endif


      !Send plane_2 to plane_1
      id_dm_send = id_dm_2(n)
      id_dm_recv = id_dm_1(n)
      if (id_bd_1(n)==1) then
        id_dir_recv = 1
        id_split_recv = 1
      else if (id_bd_1(n)==2) then
        id_dir_recv = 1
        id_split_recv = isplt_1
      else if (id_bd_1(n)==3) then
        id_dir_recv = 2
        id_split_recv = 1
      else if (id_bd_1(n)==4) then
        id_dir_recv = 2
        id_split_recv = jsplt_1
      else if (id_bd_1(n)==5) then
        id_dir_recv = 3
        id_split_recv = 1
      else if (id_bd_1(n)==6) then
        id_dir_recv = 3
        id_split_recv = ksplt_1
      endif
      if (id_bd_2(n)==1) then
        id_dir_send = 1
        id_split_send = 1
        call mpisub_sbsp_send_domainplane_to_domainplane( &
             id_dm_send,id_dir_send,id_split_send,QB_1_M, &
             id_dm_recv,id_dir_recv,id_split_recv,QD_2  ,depth,lmax,mode)
      else if (id_bd_2(n)==2) then
        id_dir_send = 1
        id_split_send = isplt_2
        call mpisub_sbsp_send_domainplane_to_domainplane( &
             id_dm_send,id_dir_send,id_split_send,QB_2_M, &
             id_dm_recv,id_dir_recv,id_split_recv,QD_2  ,depth,lmax,mode)
      else if (id_bd_2(n)==3) then
        id_dir_send = 2
        id_split_send = 1
        call mpisub_sbsp_send_domainplane_to_domainplane( &
             id_dm_send,id_dir_send,id_split_send,QB_3_M, &
             id_dm_recv,id_dir_recv,id_split_recv,QD_2  ,depth,lmax,mode)
      else if (id_bd_2(n)==4) then
        id_dir_send = 2
        id_split_send = jsplt_2
        call mpisub_sbsp_send_domainplane_to_domainplane( &
             id_dm_send,id_dir_send,id_split_send,QB_4_M, &
             id_dm_recv,id_dir_recv,id_split_recv,QD_2  ,depth,lmax,mode)
      else if (id_bd_2(n)==5) then
        id_dir_send = 3
        id_split_send = 1
        call mpisub_sbsp_send_domainplane_to_domainplane( &
             id_dm_send,id_dir_send,id_split_send,QB_5_M, &
             id_dm_recv,id_dir_recv,id_split_recv,QD_2  ,depth,lmax,mode)
      else if (id_bd_2(n)==6) then
        id_dir_send = 3
        id_split_send = ksplt_2
        call mpisub_sbsp_send_domainplane_to_domainplane( &
             id_dm_send,id_dir_send,id_split_send,QB_6_M, &
             id_dm_recv,id_dir_recv,id_split_recv,QD_2  ,depth,lmax,mode)
      endif

      if (is_plane_1==1) then
        allocate(QD_TMP(iimax_2,jjmax_2,lmax,depth))
        if (id_bd_1(n)==1.or.id_bd_1(n)==2) then
          allocate(QD_2_ALN(depth,jmax_1,kmax_1,lmax))
        else if (id_bd_1(n)==3.or.id_bd_1(n)==4) then
          allocate(QD_2_ALN(imax_1,depth,kmax_1,lmax))
        else if (id_bd_1(n)==5.or.id_bd_1(n)==6) then
          allocate(QD_2_ALN(imax_1,jmax_1,depth,lmax))
        endif

        if (id_bd_2(n)==1.or.id_bd_2(n)==2) then
          if (iabs(id_ax1_2(n))==2.and.iabs(id_ax2_2(n))==3) then
            do m=1,depth
              do l=1,lmax
                do jj=1,jjmax_2
                  do ii=1,iimax_2
                    iii = ii
                    jjj = jj
                    if (id_ax1_2(n)<0) iii = iimax_2-ii+1
                    if (id_ax2_2(n)<0) jjj = jjmax_2-jj+1
                    QD_TMP(ii,jj,l,m) = QD_2(m,iii,jjj,l)
                  enddo
                enddo
              enddo
            enddo
          else if (iabs(id_ax1_2(n))==3.and.iabs(id_ax2_2(n))==2) then
            do m=1,depth
              do l=1,lmax
                do jj=1,jjmax_2
                  do ii=1,iimax_2
                    iii = ii
                    jjj = jj
                    if (id_ax1_2(n)<0) iii = iimax_2-ii+1
                    if (id_ax2_2(n)<0) jjj = jjmax_2-jj+1
                    QD_TMP(ii,jj,l,m) = QD_2(m,jjj,iii,l)
                  enddo
                enddo
              enddo
            enddo
          endif
        else if (id_bd_2(n)==3.or.id_bd_2(n)==4) then
          if (iabs(id_ax1_2(n))==1.and.iabs(id_ax2_2(n))==3) then
            do m=1,depth
              do l=1,lmax
                do jj=1,jjmax_2
                  do ii=1,iimax_2
                    iii = ii
                    jjj = jj
                    if (id_ax1_2(n)<0) iii = iimax_2-ii+1
                    if (id_ax2_2(n)<0) jjj = jjmax_2-jj+1
                    QD_TMP(ii,jj,l,m) = QD_2(iii,m,jjj,l)
                  enddo
                enddo
              enddo
            enddo
          else if (iabs(id_ax1_2(n))==3.and.iabs(id_ax2_2(n))==1) then
            do m=1,depth
              do l=1,lmax
                do jj=1,jjmax_2
                  do ii=1,iimax_2
                    iii = ii
                    jjj = jj
                    if (id_ax1_2(n)<0) iii = iimax_2-ii+1
                    if (id_ax2_2(n)<0) jjj = jjmax_2-jj+1
                    QD_TMP(ii,jj,l,m) = QD_2(jjj,m,iii,l)
                  enddo
                enddo
              enddo
            enddo
          endif
        else if (id_bd_2(n)==5.or.id_bd_2(n)==6) then
          if (iabs(id_ax1_2(n))==1.and.iabs(id_ax2_2(n))==2) then
            do m=1,depth
              do l=1,lmax
                do jj=1,jjmax_2
                  do ii=1,iimax_2
                    iii = ii
                    jjj = jj
                    if (id_ax1_2(n)<0) iii = iimax_2-ii+1
                    if (id_ax2_2(n)<0) jjj = jjmax_2-jj+1
                    QD_TMP(ii,jj,l,m) = QD_2(iii,jjj,m,l)
                  enddo
                enddo
              enddo
            enddo
          else if (iabs(id_ax1_2(n))==2.and.iabs(id_ax2_2(n))==1) then
            do m=1,depth
              do l=1,lmax
                do jj=1,jjmax_2
                  do ii=1,iimax_2
                    iii = ii
                    jjj = jj
                    if (id_ax1_2(n)<0) iii = iimax_2-ii+1
                    if (id_ax2_2(n)<0) jjj = jjmax_2-jj+1
                    QD_TMP(ii,jj,l,m) = QD_2(jjj,iii,m,l)
                  enddo
                enddo
              enddo
            enddo
          endif
        endif

        if (id_bd_1(n)==1.or.id_bd_1(n)==2) then
          if (iabs(id_ax1_1(n))==2.and.iabs(id_ax2_1(n))==3) then
            do m=1,depth
              do l=1,lmax
                do jj=1,jjmax_1
                  do ii=1,iimax_1
                    iii = ii
                    jjj = jj
                    if (id_ax1_1(n)<0) iii = iimax_1-ii+1
                    if (id_ax2_1(n)<0) jjj = jjmax_1-jj+1
                    QD_2_ALN(m,iii,jjj,l) = QD_TMP(ii,jj,l,m)
                  enddo
                enddo
              enddo
            enddo
          else if (iabs(id_ax1_1(n))==3.and.iabs(id_ax2_1(n))==2) then
            do m=1,depth
              do l=1,lmax
                do jj=1,jjmax_1
                  do ii=1,iimax_1
                    iii = ii
                    jjj = jj
                    if (id_ax1_1(n)<0) iii = iimax_1-ii+1
                    if (id_ax2_1(n)<0) jjj = jjmax_1-jj+1
                    QD_2_ALN(m,jjj,iii,l) = QD_TMP(ii,jj,l,m)
                  enddo
                enddo
              enddo
            enddo
          endif
        else if (id_bd_1(n)==3.or.id_bd_1(n)==4) then
          if (iabs(id_ax1_1(n))==1.and.iabs(id_ax2_1(n))==3) then
            do m=1,depth
              do l=1,lmax
                do jj=1,jjmax_1
                  do ii=1,iimax_1
                    iii = ii
                    jjj = jj
                    if (id_ax1_1(n)<0) iii = iimax_1-ii+1
                    if (id_ax2_1(n)<0) jjj = jjmax_1-jj+1
                    QD_2_ALN(iii,m,jjj,l) = QD_TMP(ii,jj,l,m)
                  enddo
                enddo
              enddo
            enddo
          else if (iabs(id_ax1_1(n))==3.and.iabs(id_ax2_1(n))==1) then
            do m=1,depth
              do l=1,lmax
                do jj=1,jjmax_1
                  do ii=1,iimax_1
                    iii = ii
                    jjj = jj
                    if (id_ax1_1(n)<0) iii = iimax_1-ii+1
                    if (id_ax2_1(n)<0) jjj = jjmax_1-jj+1
                    QD_2_ALN(jjj,m,iii,l) = QD_TMP(ii,jj,l,m)
                  enddo
                enddo
              enddo
            enddo
          endif
        else if (id_bd_1(n)==5.or.id_bd_1(n)==6) then
          if (iabs(id_ax1_1(n))==1.and.iabs(id_ax2_1(n))==2) then
            do m=1,depth
              do l=1,lmax
                do jj=1,jjmax_1
                  do ii=1,iimax_1
                    iii = ii
                    jjj = jj
                    if (id_ax1_1(n)<0) iii = iimax_1-ii+1
                    if (id_ax2_1(n)<0) jjj = jjmax_1-jj+1
                    QD_2_ALN(iii,jjj,m,l) = QD_TMP(ii,jj,l,m)
                  enddo
                enddo
              enddo
            enddo
          else if (iabs(id_ax1_1(n))==2.and.iabs(id_ax2_1(n))==1) then
            do m=1,depth
              do l=1,lmax
                do jj=1,jjmax_1
                  do ii=1,iimax_1
                    iii = ii
                    jjj = jj
                    if (id_ax1_1(n)<0) iii = iimax_1-ii+1
                    if (id_ax2_1(n)<0) jjj = jjmax_1-jj+1
                    QD_2_ALN(jjj,iii,m,l) = QD_TMP(ii,jj,l,m)
                  enddo
                enddo
              enddo
            enddo
          endif
        endif

        call Get_SubDomainIJKStartEnd_From_DomainID_And_SubDomainID(id_domain,id_subdomain, &
                          istart,iend,jstart,jend,kstart,kend)
        if (mode==1) then
          iend = iend+1
          jend = jend+1
          kend = kend+1
        endif

        if (id_bd_1(n)==1) then
          do m=1,depth
            do l=1,lmax
              do k=kstart,kend
                do j=jstart,jend
                  jj = j-jstart+1
                  kk = k-kstart+1
                  mm = m
                  if (id_bd_2(n)==1.or.id_bd_2(n)==3.or.id_bd_2(n)==5) mm = depth-m+1
                  QB_1_P(m,jj,kk,l) = QD_2_ALN(mm,j,k,l)
                enddo
              enddo
            enddo
          enddo
          iflag(1) = 1
        else if (id_bd_1(n)==2) then
          do m=1,depth
            do l=1,lmax
              do k=kstart,kend
                do j=jstart,jend
                  jj = j-jstart+1
                  kk = k-kstart+1
                  mm = m
                  if (id_bd_2(n)==2.or.id_bd_2(n)==4.or.id_bd_2(n)==6) mm = depth-m+1
                  QB_2_P(m,jj,kk,l) = QD_2_ALN(mm,j,k,l)
                enddo
              enddo
            enddo
          enddo
          iflag(2) = 1
        else if (id_bd_1(n)==3) then
          do m=1,depth
            do l=1,lmax
              do k=kstart,kend
                do i=istart,iend
                  ii = i-istart+1
                  kk = k-kstart+1
                  mm = m
                  if (id_bd_2(n)==1.or.id_bd_2(n)==3.or.id_bd_2(n)==5) mm = depth-m+1
                  QB_3_P(ii,m,kk,l) = QD_2_ALN(i,mm,k,l)
                enddo
              enddo
            enddo
          enddo
          iflag(3) = 1
        else if (id_bd_1(n)==4) then
          do m=1,depth
            do l=1,lmax
              do k=kstart,kend
                do i=istart,iend
                  ii = i-istart+1
                  kk = k-kstart+1
                  mm = m
                  if (id_bd_2(n)==2.or.id_bd_2(n)==4.or.id_bd_2(n)==6) mm = depth-m+1
                  QB_4_P(ii,m,kk,l) = QD_2_ALN(i,mm,k,l)
                enddo
              enddo
            enddo
          enddo
          iflag(4) = 1
        else if (id_bd_1(n)==5) then
          do m=1,depth
            do l=1,lmax
              do j=jstart,jend
                do i=istart,iend
                  ii = i-istart+1
                  jj = j-jstart+1
                  mm = m
                  if (id_bd_2(n)==1.or.id_bd_2(n)==3.or.id_bd_2(n)==5) mm = depth-m+1
                  QB_5_P(ii,jj,m,l) = QD_2_ALN(i,j,mm,l)
                enddo
              enddo
            enddo
          enddo
          iflag(5) = 1
        else if (id_bd_1(n)==6) then
          do m=1,depth
            do l=1,lmax
              do j=jstart,jend
                do i=istart,iend
                  ii = i-istart+1
                  jj = j-jstart+1
                  mm = m
                  if (id_bd_2(n)==2.or.id_bd_2(n)==4.or.id_bd_2(n)==6) mm = depth-m+1
                  QB_6_P(ii,jj,m,l) = QD_2_ALN(i,j,mm,l)
                enddo
              enddo
            enddo
          enddo
          iflag(6) = 1
        endif

        deallocate(QD_TMP)
        deallocate(QD_2_ALN)
      endif

      if (is_plane_2==1) then
        allocate(QD_TMP(iimax_1,jjmax_1,lmax,depth))
        if (id_bd_2(n)==1.or.id_bd_2(n)==2) then
          allocate(QD_1_ALN(depth,jmax_2,kmax_2,lmax))
        else if (id_bd_2(n)==3.or.id_bd_2(n)==4) then
          allocate(QD_1_ALN(imax_2,depth,kmax_2,lmax))
        else if (id_bd_2(n)==5.or.id_bd_2(n)==6) then
          allocate(QD_1_ALN(imax_2,jmax_2,depth,lmax))
        endif

        if (id_bd_1(n)==1.or.id_bd_1(n)==2) then
          if (iabs(id_ax1_1(n))==2.and.iabs(id_ax2_1(n))==3) then
            do m=1,depth
              do l=1,lmax
                do jj=1,jjmax_1
                  do ii=1,iimax_1
                    iii = ii
                    jjj = jj
                    if (id_ax1_1(n)<0) iii = iimax_1-ii+1
                    if (id_ax2_1(n)<0) jjj = jjmax_1-jj+1
                    QD_TMP(ii,jj,l,m) = QD_1(m,iii,jjj,l)
                  enddo
                enddo
              enddo
            enddo
          else if (iabs(id_ax1_1(n))==3.and.iabs(id_ax2_1(n))==2) then
            do m=1,depth
              do l=1,lmax
                do jj=1,jjmax_1
                  do ii=1,iimax_1
                    iii = ii
                    jjj = jj
                    if (id_ax1_1(n)<0) iii = iimax_1-ii+1
                    if (id_ax2_1(n)<0) jjj = jjmax_1-jj+1
                    QD_TMP(ii,jj,l,m) = QD_1(m,jjj,iii,l)
                  enddo
                enddo
              enddo
            enddo
          endif
        else if (id_bd_1(n)==3.or.id_bd_1(n)==4) then
          if (iabs(id_ax1_1(n))==1.and.iabs(id_ax2_1(n))==3) then
            do m=1,depth
              do l=1,lmax
                do jj=1,jjmax_1
                  do ii=1,iimax_1
                    iii = ii
                    jjj = jj
                    if (id_ax1_1(n)<0) iii = iimax_1-ii+1
                    if (id_ax2_1(n)<0) jjj = jjmax_1-jj+1
                    QD_TMP(ii,jj,l,m) = QD_1(iii,m,jjj,l)
                  enddo
                enddo
              enddo
            enddo
          else if (iabs(id_ax1_1(n))==3.and.iabs(id_ax2_1(n))==1) then
            do m=1,depth
              do l=1,lmax
                do jj=1,jjmax_1
                  do ii=1,iimax_1
                    iii = ii
                    jjj = jj
                    if (id_ax1_1(n)<0) iii = iimax_1-ii+1
                    if (id_ax2_1(n)<0) jjj = jjmax_1-jj+1
                    QD_TMP(ii,jj,l,m) = QD_1(jjj,m,iii,l)
                  enddo
                enddo
              enddo
            enddo
          endif
        else if (id_bd_1(n)==5.or.id_bd_1(n)==6) then
          if (iabs(id_ax1_1(n))==1.and.iabs(id_ax2_1(n))==2) then
            do m=1,depth
              do l=1,lmax
                do jj=1,jjmax_1
                  do ii=1,iimax_1
                    iii = ii
                    jjj = jj
                    if (id_ax1_1(n)<0) iii = iimax_1-ii+1
                    if (id_ax2_1(n)<0) jjj = jjmax_1-jj+1
                    QD_TMP(ii,jj,l,m) = QD_1(iii,jjj,m,l)
                  enddo
                enddo
              enddo
            enddo
          else if (iabs(id_ax1_1(n))==2.and.iabs(id_ax2_1(n))==1) then
            do m=1,depth
              do l=1,lmax
                do jj=1,jjmax_1
                  do ii=1,iimax_1
                    iii = ii
                    jjj = jj
                    if (id_ax1_1(n)<0) iii = iimax_1-ii+1
                    if (id_ax2_1(n)<0) jjj = jjmax_1-jj+1
                    QD_TMP(ii,jj,l,m) = QD_1(jjj,iii,m,l)
                  enddo
                enddo
              enddo
            enddo
          endif
        endif

        if (id_bd_2(n)==1.or.id_bd_2(n)==2) then
          if (iabs(id_ax1_2(n))==2.and.iabs(id_ax2_2(n))==3) then
            do m=1,depth
              do l=1,lmax
                do jj=1,jjmax_2
                  do ii=1,iimax_2
                    iii = ii
                    jjj = jj
                    if (id_ax1_2(n)<0) iii = iimax_2-ii+1
                    if (id_ax2_2(n)<0) jjj = jjmax_2-jj+1
                    QD_1_ALN(m,iii,jjj,l) = QD_TMP(ii,jj,l,m)
                  enddo
                enddo
              enddo
            enddo
          else if (iabs(id_ax1_2(n))==3.and.iabs(id_ax2_2(n))==2) then
            do m=1,depth
              do l=1,lmax
                do jj=1,jjmax_2
                  do ii=1,iimax_2
                    iii = ii
                    jjj = jj
                    if (id_ax1_2(n)<0) iii = iimax_2-ii+1
                    if (id_ax2_2(n)<0) jjj = jjmax_2-jj+1
                    QD_1_ALN(m,jjj,iii,l) = QD_TMP(ii,jj,l,m)
                  enddo
                enddo
              enddo
            enddo
          endif
        else if (id_bd_2(n)==3.or.id_bd_2(n)==4) then
          if (iabs(id_ax1_2(n))==1.and.iabs(id_ax2_2(n))==3) then
            do m=1,depth
              do l=1,lmax
                do jj=1,jjmax_2
                  do ii=1,iimax_2
                    iii = ii
                    jjj = jj
                    if (id_ax1_2(n)<0) iii = iimax_2-ii+1
                    if (id_ax2_2(n)<0) jjj = jjmax_2-jj+1
                    QD_1_ALN(iii,m,jjj,l) = QD_TMP(ii,jj,l,m)
                  enddo
                enddo
              enddo
            enddo
          else if (iabs(id_ax1_2(n))==3.and.iabs(id_ax2_2(n))==1) then
            do m=1,depth
              do l=1,lmax
                do jj=1,jjmax_2
                  do ii=1,iimax_2
                    iii = ii
                    jjj = jj
                    if (id_ax1_2(n)<0) iii = iimax_2-ii+1
                    if (id_ax2_2(n)<0) jjj = jjmax_2-jj+1
                    QD_1_ALN(jjj,m,iii,l) = QD_TMP(ii,jj,l,m)
                  enddo
                enddo
              enddo
            enddo
          endif
        else if (id_bd_2(n)==5.or.id_bd_2(n)==6) then
          if (iabs(id_ax1_2(n))==1.and.iabs(id_ax2_2(n))==2) then
            do m=1,depth
              do l=1,lmax
                do jj=1,jjmax_2
                  do ii=1,iimax_2
                    iii = ii
                    jjj = jj
                    if (id_ax1_2(n)<0) iii = iimax_2-ii+1
                    if (id_ax2_2(n)<0) jjj = jjmax_2-jj+1
                    QD_1_ALN(iii,jjj,m,l) = QD_TMP(ii,jj,l,m)
                  enddo
                enddo
              enddo
            enddo
          else if (iabs(id_ax1_2(n))==2.and.iabs(id_ax2_2(n))==1) then
            do m=1,depth
              do l=1,lmax
                do jj=1,jjmax_2
                  do ii=1,iimax_2
                    iii = ii
                    jjj = jj
                    if (id_ax1_2(n)<0) iii = iimax_2-ii+1
                    if (id_ax2_2(n)<0) jjj = jjmax_2-jj+1
                    QD_1_ALN(jjj,iii,m,l) = QD_TMP(ii,jj,l,m)
                  enddo
                enddo
              enddo
            enddo
          endif
        endif

        call Get_SubDomainIJKStartEnd_From_DomainID_And_SubDomainID(id_domain,id_subdomain, &
                          istart,iend,jstart,jend,kstart,kend)
        if (mode==1) then
          iend = iend+1
          jend = jend+1
          kend = kend+1
        endif

        if (id_bd_2(n)==1) then
          do m=1,depth
            do l=1,lmax
              do k=kstart,kend
                do j=jstart,jend
                  jj = j-jstart+1
                  kk = k-kstart+1
                  mm = m
                  if (id_bd_1(n)==1.or.id_bd_1(n)==3.or.id_bd_1(n)==5) mm = depth-m+1
                  QB_1_P(m,jj,kk,l) = QD_1_ALN(mm,j,k,l)
                enddo
              enddo
            enddo
          enddo
          iflag(1) = 1
        else if (id_bd_2(n)==2) then
          do m=1,depth
            do l=1,lmax
              do k=kstart,kend
                do j=jstart,jend
                  jj = j-jstart+1
                  kk = k-kstart+1
                  mm = m
                  if (id_bd_1(n)==2.or.id_bd_1(n)==4.or.id_bd_1(n)==6) mm = depth-m+1
                  QB_2_P(m,jj,kk,l) = QD_1_ALN(mm,j,k,l)
                enddo
              enddo
            enddo
          enddo
          iflag(2) = 1
        else if (id_bd_2(n)==3) then
          do m=1,depth
            do l=1,lmax
              do k=kstart,kend
                do i=istart,iend
                  ii = i-istart+1
                  kk = k-kstart+1
                  mm = m
                  if (id_bd_1(n)==1.or.id_bd_1(n)==3.or.id_bd_1(n)==5) mm = depth-m+1
                  QB_3_P(ii,m,kk,l) = QD_1_ALN(i,mm,k,l)
                enddo
              enddo
            enddo
          enddo
          iflag(3) = 1
        else if (id_bd_2(n)==4) then
          do m=1,depth
            do l=1,lmax
              do k=kstart,kend
                do i=istart,iend
                  ii = i-istart+1
                  kk = k-kstart+1
                  mm = m
                  if (id_bd_1(n)==2.or.id_bd_1(n)==4.or.id_bd_1(n)==6) mm = depth-m+1
                  QB_4_P(ii,m,kk,l) = QD_1_ALN(i,mm,k,l)
                enddo
              enddo
            enddo
          enddo
          iflag(4) = 1
        else if (id_bd_2(n)==5) then
          do m=1,depth
            do l=1,lmax
              do j=jstart,jend
                do i=istart,iend
                  ii = i-istart+1
                  jj = j-jstart+1
                  mm = m
                  if (id_bd_1(n)==1.or.id_bd_1(n)==3.or.id_bd_1(n)==5) mm = depth-m+1
                  QB_5_P(ii,jj,m,l) = QD_1_ALN(i,j,mm,l)
                enddo
              enddo
            enddo
          enddo
          iflag(5) = 1
        else if (id_bd_2(n)==6) then
          do m=1,depth
            do l=1,lmax
              do j=jstart,jend
                do i=istart,iend
                  ii = i-istart+1
                  jj = j-jstart+1
                  mm = m
                  if (id_bd_1(n)==2.or.id_bd_1(n)==4.or.id_bd_1(n)==6) mm = depth-m+1
                  QB_6_P(ii,jj,m,l) = QD_1_ALN(i,j,mm,l)
                enddo
              enddo
            enddo
          enddo
          iflag(6) = 1
        endif

        deallocate(QD_TMP)
        deallocate(QD_1_ALN)
      endif

      if (is_plane_1==1) then
        deallocate(QD_2)
      endif

      if (is_plane_2==1) then
        deallocate(QD_1)
      endif

    enddo!----------------------------------------------------------------------

  end subroutine exchange_domain_halo

end module mdl_exchange_domain_halo
