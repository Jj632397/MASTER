module mdl_table
#ifdef _USE_MPI
  use mpi
#endif
  use mdl_look_up_table
  implicit none

  type(look_up_table), public :: TBL_TP_to_DNST
  type(look_up_table), public :: TBL_TP_to_HTPS
  type(look_up_table), public :: TBL_TP_to_DRTM
  type(look_up_table), public :: TBL_TP_to_DRPR
  type(look_up_table), public :: TBL_TP_to_SHCP
  type(look_up_table), public :: TBL_TP_to_DHPR
  type(look_up_table), public :: TBL_TP_to_VSSV
  type(look_up_table), public :: TBL_TP_to_CKPV
  type(look_up_table), public :: TBL_TP_to_DIELEC

contains

  subroutine table_init
    integer :: myrank,nprocs,iproc
    integer :: ierror

#ifdef _USE_MPI
    call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierror)
    call MPI_Comm_size(MPI_COMM_WORLD,nprocs,ierror)
#else
    myrank = 0
    nprocs = 1
#endif

    do iproc=0,nprocs-1
      if (myrank==iproc) then
        call look_up_table_read("./data/table_water/TBL_TP_to_DNST.tbl",TBL_TP_to_DNST)
        call look_up_table_read("./data/table_water/TBL_TP_to_HTPS.tbl",TBL_TP_to_HTPS)
        call look_up_table_read("./data/table_water/TBL_TP_to_DRTM.tbl",TBL_TP_to_DRTM)
        call look_up_table_read("./data/table_water/TBL_TP_to_DRPR.tbl",TBL_TP_to_DRPR)
        call look_up_table_read("./data/table_water/TBL_TP_to_SHCP.tbl",TBL_TP_to_SHCP)
        call look_up_table_read("./data/table_water/TBL_TP_to_DHPR.tbl",TBL_TP_to_DHPR)
        call look_up_table_read("./data/table_water/TBL_TP_to_VSSV.tbl",TBL_TP_to_VSSV)
        call look_up_table_read("./data/table_water/TBL_TP_to_CKPV.tbl",TBL_TP_to_CKPV)
        call look_up_table_read("./data/table_water/TBL_TP_to_DIELEC.tbl",TBL_TP_to_DIELEC)
      endif
#ifdef _USE_MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
#endif
    enddo
  end subroutine table_init
  !-----------------------------------------------------------------------------
  subroutine table_final
    call look_up_table_deallocate(TBL_TP_to_DNST)
    call look_up_table_deallocate(TBL_TP_to_HTPS)
    call look_up_table_deallocate(TBL_TP_to_DRTM)
    call look_up_table_deallocate(TBL_TP_to_DRPR)
    call look_up_table_deallocate(TBL_TP_to_SHCP)
    call look_up_table_deallocate(TBL_TP_to_DHPR)
    call look_up_table_deallocate(TBL_TP_to_VSSV)
    call look_up_table_deallocate(TBL_TP_to_CKPV)
    call look_up_table_deallocate(TBL_TP_to_DIELEC)
  end subroutine table_final

end module mdl_table
