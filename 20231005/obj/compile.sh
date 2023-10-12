#!/bin/sh

FC="mpinfort -D_USE_MPI -O3 -ftrace -fdiag-vector=2 -report-format"

#FC="nfort -O3 -ftrace -fdiag-vector=2 -report-format"

EXE="../LetsDoThis"

OBJ1="input_manager.o decompo.o mpisub_sbsp.o steps.o"
OBJ2="band.o fdm_1st.o lp_filter.o muscl.o weno.o scmm.o look_up_table.o param.o refs.o"
OBJ3="table.o initmp_tpipe.o exchange_domain_halo.o block.o halo.o metric.o time_delta.o"
OBJ4="sgs_wale.o tbl_sa.o schs_source.o"
OBJ5="gradient.o flux_pw.o flux_riemann_roe.o flux_riemann_fvs_real_gas_extension.o flux_riemann_fds_real_gas_extension.o flux_viscous.o result_PT.o"
OBJ6="readgrid.o output.o restart.o"
OBJ7="right_hand_side.o lu_sgs_real_gas_extension.o"
OBJ8="bound_wall.o bound_inflow_outflow_tpipe.o bound.o filter.o explicit.o implicit.o"
OBJ9="distance_from_wall.o init_tpipe.o means.o output_tpipe.o main.o"

rm $OBJ1 $OBJ2 $OBJ3 $OBJ4 $OBJ5 $OBJ6 $OBJ7 $OBJ8 $OBJ9

$FC -c ../src/input_manager.F90 -o input_manager.o

$FC -c ../src/decompo.F90 -o decompo.o

$FC -c ../src/mpisub_sbsp.F90 -o mpisub_sbsp.o

$FC -c ../src/steps.F90 -o steps.o

$FC -c ../src/band.F90 -o band.o

$FC -c ../src/fdm_1st.F90 -o fdm_1st.o

$FC -c ../src/lp_filter.F90 -o lp_filter.o

$FC -c ../src/muscl.F90 -o muscl.o

$FC -c ../src/weno.F90 -o weno.o

$FC -c ../src/scmm.F90 -o scmm.o

$FC -c ../src/look_up_table.F90 -o look_up_table.o

$FC -c ../src/param.F90 -o param.o

$FC -c ../src/refs.F90 -o refs.o

$FC -c ../src/table.F90 -o table.o

$FC -c ../src/initmp_tpipe.F90 -o initmp_tpipe.o

$FC -c ../src/exchange_domain_halo.F90 -o exchange_domain_halo.o

$FC -c ../src/block.F90 -o block.o

$FC -c ../src/halo.F90 -o halo.o

$FC -c ../src/metric.F90 -o metric.o

$FC -c ../src/time_delta.F90 -o time_delta.o

$FC -c ../src/sgs_wale.F90 -o sgs_wale.o

$FC -c ../src/tbl_sa.F90 -o tbl_sa.o

$FC -c ../src/schs_source.F90 -o schs_source.o

$FC -c ../src/gradient.F90 -o gradient.o

$FC -c ../src/flux_pw.F90 -o flux_pw.o

$FC -c ../src/flux_riemann_roe.F90 -o flux_riemann_roe.o

$FC -c ../src/flux_riemann_fvs_real_gas_extension.F90 -o flux_riemann_fvs_real_gas_extension.o

$FC -c ../src/flux_riemann_fds_real_gas_extension.F90 -o flux_riemann_fds_real_gas_extension.o

$FC -c ../src/flux_viscous.F90 -o flux_viscous.o

$FC -c ../src/result_PT.F90 -o result_PT.o

$FC -c ../src/readgrid.F90 -o readgrid.o

$FC -c ../src/output.F90 -o output.o

$FC -c ../src/restart.F90 -o restart.o

$FC -c ../src/right_hand_side.F90 -o right_hand_side.o

$FC -c ../src/lu_sgs_real_gas_extension.F90 -o lu_sgs_real_gas_extension.o

$FC -c ../src/bound_wall.F90 -o bound_wall.o

$FC -c ../src/bound_inflow_outflow_tpipe.F90 -o bound_inflow_outflow_tpipe.o

$FC -c ../src/bound.F90 -o bound.o

$FC -c ../src/filter.F90 -o filter.o

$FC -c ../src/explicit.F90 -o explicit.o

$FC -c ../src/implicit.F90 -o implicit.o

$FC -c ../src/distance_from_wall.F90 -o distance_from_wall.o

$FC -c ../src/init_tpipe.F90 -o init_tpipe.o

$FC -c ../src/means.F90 -o means.o

$FC -c ../src/output_tpipe.F90 -o output_tpipe.o

$FC -c ../src/main.F90 -o main.o

$FC $OBJ1 $OBJ2 $OBJ3 $OBJ4 $OBJ5 $OBJ6 $OBJ7 $OBJ8 $OBJ9 -o $EXE
