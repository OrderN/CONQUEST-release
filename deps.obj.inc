main.o: DiagModule.o ScalapackFormat.o atomic_density.o atoms.module.o \
	blip_module.o block_module.o control.o cover_module.o \
	datatypes.module.o density_module.o dimens.module.o ewald.o \
	fft_module.o force_module.o functions_on_grid_module.o \
	generic_comms.o global_module.o grid_index.module.o group_module.o \
	initialisation_module.o logicals_module.o matrix_data_module.o \
	memory_module.o minimise.o mult_module.o numbers.module.o \
	ol_int_datatypes.module.o pao_format.o potential_module.o \
	primary_module.o pseudo_tm_info.o pseudopotential.module.o \
	pseudopotential_common.o species_module.o support_spec_format.o \
	timer_stdclocks_module.o
numbers.module.o: datatypes.module.o
atoms.module.o: datatypes.module.o generic_comms.o global_module.o \
	group_module.o memory_module.o primary_module.o timer_module.o \
	timer_stdclocks_module.o
mult_module.o: basic_types.o comms_module.o cover_module.o datatypes.module.o \
	generic_blas.o generic_comms.o global_module.o matrix_comms_module.o \
	matrix_data_module.o matrix_elements_module.o matrix_module.o \
	maxima_module.o memory_module.o mpi.o mult_init_module.o \
	multiply_module.o numbers.module.o primary_module.o timer_module.o \
	timer_stdclocks_module.o
multiply_module.o: basic_types.o comms_module.o datatypes.module.o \
	generic_comms.o global_module.o matrix_comms_module.o matrix_module.o \
	maxima_module.o mpi.o multiply_kernel_default.o numbers.module.o \
	primary_module.o timer_module.o timer_stdclocks_module.o \
	trans_module.o
multiply_kernel_default.o: basic_types.o datatypes.module.o matrix_module.o \
	primary_module.o
io_module.o: basic_types.o block_module.o construct_module.o cover_module.o \
	datatypes.module.o dimens.module.o functions_on_grid_module.o \
	generic_comms.o global_module.o group_module.o input_module.o \
	matrix_data_module.o matrix_module.o maxima_module.o memory_module.o \
	mpi.o mult_module.o numbers.module.o primary_module.o \
	sfc_partitions_module.o species_module.o timer_module.o \
	timer_stdclocks_module.o units.module.o
sfc_partitions_module.o: Hilbert3D.o basic_types.o construct_module.o \
	datatypes.module.o dimens.module.o generic_comms.o global_module.o \
	maxima_module.o memory_module.o numbers.module.o species_module.o \
	timer_module.o timer_stdclocks_module.o
Hilbert3D.o: generic_comms.o global_module.o memory_module.o
matrix_elements_module.o: basic_types.o datatypes.module.o generic_comms.o \
	global_module.o group_module.o matrix_module.o maxima_module.o mpi.o \
	species_module.o trans_module.o
comms_module.o: datatypes.module.o generic_comms.o global_module.o \
	matrix_comms_module.o matrix_data_module.o matrix_module.o mpi.o \
	support_spec_format.o
matrix_comms_module.o: basic_types.o datatypes.module.o generic_comms.o \
	global_module.o matrix_module.o maxima_module.o mpi.o
basic_types.o: datatypes.module.o
group_module.o: basic_types.o datatypes.module.o generic_comms.o \
	global_module.o memory_module.o timer_module.o \
	timer_stdclocks_module.o
mult_init_module.o: basic_types.o datatypes.module.o generic_comms.o \
	matrix_comms_module.o matrix_module.o
trans_module.o: basic_types.o datatypes.module.o generic_comms.o \
	matrix_comms_module.o matrix_module.o mpi.o
matrix_module.o: basic_types.o datatypes.module.o generic_comms.o \
	global_module.o timer_module.o timer_stdclocks_module.o
construct_module.o: basic_types.o cover_module.o group_module.o \
	primary_module.o timer_module.o
primary_module.o: basic_types.o datatypes.module.o generic_comms.o \
	global_module.o timer_module.o
cover_module.o: basic_types.o datatypes.module.o generic_comms.o \
	global_module.o memory_module.o timer_module.o \
	timer_stdclocks_module.o
maxima_module.o: matrix_module.o
generic_blas.o: datatypes.module.o global_module.o mpi.o numbers.module.o
matrix_data_module.o: datatypes.module.o matrix_module.o
atom_dispenser_module.o: cover_module.o datatypes.module.o dimens.module.o \
	generic_comms.o global_module.o group_module.o input_module.o \
	io_module.o numbers.module.o
io_module2.o: basic_types.o cover_module.o datatypes.module.o generic_comms.o \
	global_module.o group_module.o input_module.o io_module.o \
	matrix_data_module.o mult_module.o numbers.module.o primary_module.o \
	species_module.o
UpdateInfo_module.o: atom_dispenser_module.o cover_module.o \
	datatypes.module.o generic_comms.o global_module.o group_module.o \
	input_module.o io_module.o io_module2.o matrix_data_module.o \
	matrix_module.o mpi.o mult_module.o numbers.module.o primary_module.o
UpdateMember_module.o: atom_dispenser_module.o basic_types.o cover_module.o \
	datatypes.module.o density_module.o ewald.o generic_comms.o \
	global_module.o group_module.o input_module.o io_module.o \
	maxima_module.o primary_module.o species_module.o timer_module.o \
	timer_stdclocks_module.o
XLBOMD_module.o: UpdateInfo_module.o basic_types.o datatypes.module.o \
	generic_comms.o global_module.o input_module.o io_module.o \
	io_module2.o matrix_data_module.o matrix_module.o mult_init_module.o \
	mult_module.o numbers.module.o
generic_comms.o: datatypes.module.o global_module.o mpi.o timer_module.o
global_module.o: datatypes.module.o numbers.module.o
hartree.o: datatypes.module.o dimens.module.o fft_module.o generic_comms.o \
	global_module.o memory_module.o numbers.module.o
units.module.o: datatypes.module.o numbers.module.o
input_module.o: datatypes.module.o generic_comms.o timer_module.o
memory_module.o: datatypes.module.o generic_comms.o global_module.o \
	units.module.o
timer_module.o: datatypes.module.o global_module.o mpi.o
timer_stdclocks_module.o: timer_module.o
auxiliary_types.o: datatypes.module.o
fft_module.o: datatypes.module.o dimens.module.o fft_ffte.o generic_comms.o \
	global_module.o grid_index.module.o maxima_module.o memory_module.o \
	numbers.module.o timer_module.o timer_stdclocks_module.o
fft_ffte.o: datatypes.module.o
blip_module.o: datatypes.module.o dimens.module.o fft_module.o generic_blas.o \
	generic_comms.o global_module.o group_module.o matrix_data_module.o \
	memory_module.o mult_module.o numbers.module.o primary_module.o \
	species_module.o support_spec_format.o timer_module.o \
	timer_stdclocks_module.o
dimens.module.o: block_module.o datatypes.module.o generic_comms.o \
	global_module.o input_module.o matrix_data_module.o numbers.module.o \
	pseudopotential_common.o species_module.o timer_module.o
species_module.o: datatypes.module.o
spline_module.o: datatypes.module.o generic_comms.o global_module.o \
	memory_module.o numbers.module.o
H_matrix_module.o: PAO_grid_transform_module.o XC_module.o \
	blip_grid_transform_module.o blip_module.o block_module.o \
	build_PAO_matrices.module.o calc_matrix_elements_module.o \
	cdft_data_module.o cover_module.o datatypes.module.o density_module.o \
	dimens.module.o energy.module.o exx_kernel_default.o exx_module.o \
	exx_types.o fft_module.o functions_on_grid_module.o generic_blas.o \
	generic_comms.o global_module.o group_module.o hartree.o io_module.o \
	matrix_data_module.o maxima_module.o memory_module.o mult_module.o \
	nlpf2blip.o numbers.module.o potential_module.o primary_module.o \
	pseudo_tm_info.o pseudopotential.module.o pseudopotential_common.o \
	set_blipgrid_module.o set_bucket_module.o species_module.o \
	support_spec_format.o timer_module.o timer_stdclocks_module.o
SelfCon_module.o: DMMinModule.o EarlySC_module.o H_matrix_module.o \
	PosTan_module.o Pulay.o atoms.module.o datatypes.module.o \
	density_module.o dimens.module.o energy.module.o generic_blas.o \
	generic_comms.o global_module.o hartree.o input_module.o io_module.o \
	matrix_data_module.o maxima_module.o memory_module.o mult_module.o \
	numbers.module.o primary_module.o timer_module.o \
	timer_stdclocks_module.o
EarlySC_module.o: DMMinModule.o DiagModule.o H_matrix_module.o cdft_module.o \
	datatypes.module.o density_module.o dimens.module.o energy.module.o \
	functions_on_grid_module.o generic_blas.o generic_comms.o \
	global_module.o logicals_module.o memory_module.o mult_module.o \
	numbers.module.o timer_module.o timer_stdclocks_module.o
PosTan_module.o: datatypes.module.o generic_comms.o global_module.o \
	memory_module.o numbers.module.o
DMMinModule.o: DiagModule.o H_matrix_module.o McWeeny.o PosTan_module.o \
	Pulay.o XLBOMD_module.o datatypes.module.o density_module.o \
	dimens.module.o energy.module.o exx_kernel_default.o \
	functions_on_grid_module.o generic_comms.o global_module.o \
	io_module.o io_module2.o logicals_module.o matrix_data_module.o \
	maxima_module.o memory_module.o mult_module.o numbers.module.o \
	primary_module.o timer_module.o timer_stdclocks_module.o
DiagModule.o: ScalapackFormat.o cdft_data_module.o cover_module.o \
	datatypes.module.o energy.module.o generic_blas.o generic_comms.o \
	global_module.o group_module.o matrix_data_module.o matrix_module.o \
	maxima_module.o memory_module.o mpi.o mult_module.o numbers.module.o \
	primary_module.o species_module.o timer_module.o \
	timer_stdclocks_module.o units.module.o
ScalapackFormat.o: generic_comms.o global_module.o group_module.o \
	maxima_module.o species_module.o timer_module.o \
	timer_stdclocks_module.o
blip_minimisation.module.o: DiagModule.o PosTan_module.o S_matrix_module.o \
	SelfCon_module.o blip_gradient.module.o blip_module.o \
	datatypes.module.o dimens.module.o functions_on_grid_module.o \
	generic_blas.o generic_comms.o global_module.o io_module.o \
	logicals_module.o memory_module.o mult_module.o numbers.module.o \
	primary_module.o support_spec_format.o timer_module.o \
	timer_stdclocks_module.o units.module.o
blip_gradient.module.o: blip_grid_transform_module.o blip_module.o \
	calc_matrix_elements_module.o cover_module.o datatypes.module.o \
	dimens.module.o functions_on_grid_module.o generic_blas.o \
	generic_comms.o global_module.o group_module.o matrix_data_module.o \
	memory_module.o mult_module.o nlpf2blip.o numbers.module.o \
	primary_module.o pseudopotential.module.o set_bucket_module.o \
	species_module.o support_spec_format.o timer_module.o \
	timer_stdclocks_module.o
pseudopotential.module.o: block_module.o cover_module.o datatypes.module.o \
	dimens.module.o functions_on_grid_module.o generic_blas.o \
	generic_comms.o global_module.o group_module.o hartree.o \
	input_module.o maxima_module.o numbers.module.o primary_module.o \
	pseudopotential_common.o set_blipgrid_module.o species_module.o \
	spline_module.o timer_module.o timer_stdclocks_module.o
potential_module.o: datatypes.module.o
ewald.o: atoms.module.o construct_module.o cover_module.o datatypes.module.o \
	dimens.module.o generic_comms.o global_module.o group_module.o \
	maxima_module.o memory_module.o numbers.module.o primary_module.o \
	species_module.o timer_module.o units.module.o
minimise.o: DMMinModule.o SelfCon_module.o blip_minimisation.module.o \
	datatypes.module.o density_module.o energy.module.o force_module.o \
	generic_comms.o global_module.o input_module.o \
	pao_minimisation.module.o timer_module.o timer_stdclocks_module.o \
	units.module.o vdWDFT_module.o
Pulay.o: datatypes.module.o generic_blas.o generic_comms.o global_module.o \
	numbers.module.o
S_matrix_module.o: DiagModule.o PAO_grid_transform_module.o \
	UpdateInfo_module.o blip_grid_transform_module.o blip_module.o \
	build_PAO_matrices.module.o calc_matrix_elements_module.o \
	cdft_module.o comms_module.o cover_module.o datatypes.module.o \
	functions_on_grid_module.o generic_blas.o generic_comms.o \
	global_module.o group_module.o input_module.o io_module.o \
	io_module2.o matrix_data_module.o mpi.o mult_module.o nlpf2blip.o \
	numbers.module.o primary_module.o set_bucket_module.o \
	species_module.o support_spec_format.o timer_module.o \
	timer_stdclocks_module.o
McWeeny.o: cover_module.o datatypes.module.o generic_blas.o generic_comms.o \
	global_module.o matrix_data_module.o matrix_module.o maxima_module.o \
	memory_module.o mult_module.o numbers.module.o primary_module.o \
	species_module.o timer_module.o timer_stdclocks_module.o
control.o: DiagModule.o Integrators_module.o UpdateInfo_module.o \
	XLBOMD_module.o constraint_module.o datatypes.module.o \
	dimens.module.o ewald.o force_module.o generic_blas.o generic_comms.o \
	global_module.o group_module.o input_module.o io_module.o \
	io_module2.o matrix_data_module.o memory_module.o minimise.o \
	move_atoms.module.o mult_module.o numbers.module.o primary_module.o \
	pseudopotential.module.o timer_module.o units.module.o
density_module.o: atomic_density.o block_module.o \
	calc_matrix_elements_module.o cdft_data_module.o cover_module.o \
	datatypes.module.o dimens.module.o functions_on_grid_module.o \
	generic_blas.o generic_comms.o global_module.o group_module.o \
	mult_module.o numbers.module.o primary_module.o pseudo_tm_info.o \
	set_blipgrid_module.o set_bucket_module.o spline_module.o \
	timer_module.o timer_stdclocks_module.o units.module.o
initialisation_module.o: DFT_D2_module.o DMMinModule.o DiagModule.o \
	H_matrix_module.o S_matrix_module.o SelfCon_module.o \
	UpdateInfo_module.o XLBOMD_module.o atoms.module.o \
	blip_grid_transform_module.o blip_module.o block_module.o \
	calc_matrix_elements_module.o cdft_module.o construct_module.o \
	cover_module.o datatypes.module.o density_module.o dimens.module.o \
	energy.module.o ewald.o fft_module.o functions_on_grid_module.o \
	generic_comms.o global_module.o grid_index.module.o group_module.o \
	initial_read.module.o input_module.o io_module.o io_module2.o \
	ionic_data.o logicals_module.o matrix_data_module.o maxima_module.o \
	memory_module.o minimise.o mult_module.o nlpf2blip.o numbers.module.o \
	ol_ang_coeff_subs.o ol_rad_table_subs.o pao2blip.o potential_module.o \
	primary_module.o pseudo_tm_module.o pseudopotential.module.o \
	pseudopotential_common.o read_support_spec.o set_blipgrid_module.o \
	set_bucket_module.o species_module.o support_spec_format.o \
	test_force_module.o timer_module.o timer_stdclocks_module.o
initial_read.module.o: DMMinModule.o DiagModule.o H_matrix_module.o \
	S_matrix_module.o ScalapackFormat.o SelfCon_module.o XLBOMD_module.o \
	atomic_density.o auxiliary_types.o blip_module.o block_module.o \
	cdft_data_module.o constraint_module.o construct_module.o control.o \
	datatypes.module.o datestamp.o dimens.module.o energy.module.o \
	ewald.o exx_types.o force_module.o generic_comms.o global_module.o \
	group_module.o input_module.o io_module.o maxima_module.o \
	memory_module.o minimise.o numbers.module.o pao_format.o \
	pao_minimisation.module.o primary_module.o pseudo_tm_module.o \
	pseudopotential.module.o pseudopotential_common.o read_pao_info.o \
	read_support_spec.o sfc_partitions_module.o species_module.o \
	support_spec_format.o test_force_module.o timer_module.o \
	timer_stdclocks_module.o units.module.o
energy.module.o: DFT_D2_module.o datatypes.module.o density_module.o ewald.o \
	generic_comms.o global_module.o mult_module.o numbers.module.o \
	pseudopotential_common.o timer_module.o units.module.o
test_force_module.o: DMMinModule.o H_matrix_module.o \
	PAO_grid_transform_module.o S_matrix_module.o SelfCon_module.o \
	blip_grid_transform_module.o cdft_module.o cover_module.o \
	datatypes.module.o density_module.o energy.module.o ewald.o \
	force_module.o functions_on_grid_module.o generic_comms.o \
	global_module.o group_module.o io_module.o matrix_data_module.o \
	maxima_module.o memory_module.o move_atoms.module.o mult_module.o \
	numbers.module.o primary_module.o pseudo_tm_module.o \
	pseudopotential.module.o pseudopotential_common.o timer_module.o \
	timer_stdclocks_module.o
functions_on_grid_module.o: datatypes.module.o generic_comms.o \
	global_module.o memory_module.o numbers.module.o timer_module.o \
	timer_stdclocks_module.o
cdft_data_module.o: datatypes.module.o
cdft_module.o: DMMinModule.o DiagModule.o cdft_data_module.o \
	datatypes.module.o density_module.o energy.module.o generic_comms.o \
	global_module.o io_module.o logicals_module.o matrix_data_module.o \
	maxima_module.o mult_module.o numbers.module.o
DFT_D2_module.o: cover_module.o datatypes.module.o dimens.module.o \
	generic_comms.o global_module.o group_module.o input_module.o \
	io_module.o numbers.module.o primary_module.o species_module.o \
	units.module.o
XC_module.o: block_module.o datatypes.module.o dimens.module.o fft_module.o \
	generic_comms.o global_module.o group_module.o memory_module.o \
	numbers.module.o
vdWMesh_module.o: datatypes.module.o generic_comms.o global_module.o \
	memory_module.o numbers.module.o spline_module.o
vdWDFT_module.o: XC_module.o datatypes.module.o dimens.module.o fft_module.o \
	generic_comms.o global_module.o grid_index.module.o input_module.o \
	maxima_module.o memory_module.o numbers.module.o ol_fft_procs.o \
	spline_module.o timer_module.o vdWMesh_module.o
atomic_density.o: datatypes.module.o generic_comms.o global_module.o \
	input_module.o memory_module.o numbers.module.o pao_format.o \
	spline_module.o timer_module.o timer_stdclocks_module.o
ionic_data.o: atomic_density.o blip_module.o datatypes.module.o \
	generic_comms.o global_module.o input_module.o ol_rad_table_subs.o \
	pseudopotential_common.o read_pao_info.o species_module.o \
	timer_module.o
pao_format.o: datatypes.module.o
read_pao_info.o: datatypes.module.o generic_comms.o global_module.o \
	input_module.o memory_module.o numbers.module.o pao_format.o \
	spline_module.o timer_module.o timer_stdclocks_module.o
force_module.o: DFT_D2_module.o DiagModule.o H_matrix_module.o \
	PAO_grid_transform_module.o S_matrix_module.o XC_module.o \
	atomic_density.o blip_gradient.module.o blip_grid_transform_module.o \
	blip_module.o block_module.o build_PAO_matrices.module.o \
	calc_matrix_elements_module.o comms_module.o cover_module.o \
	datatypes.module.o density_module.o dimens.module.o ewald.o \
	functions_on_grid_module.o generic_blas.o generic_comms.o \
	global_module.o grid_index.module.o group_module.o hartree.o \
	logicals_module.o matrix_data_module.o matrix_module.o \
	maxima_module.o memory_module.o mpi.o mult_module.o nlpf2blip.o \
	numbers.module.o potential_module.o primary_module.o pseudo_tm_info.o \
	pseudo_tm_module.o pseudopotential.module.o pseudopotential_common.o \
	set_blipgrid_module.o set_bucket_module.o species_module.o \
	spline_module.o support_spec_format.o timer_module.o \
	timer_stdclocks_module.o units.module.o
move_atoms.module.o: DFT_D2_module.o DiagModule.o H_matrix_module.o Pulay.o \
	S_matrix_module.o SelfCon_module.o UpdateInfo_module.o \
	UpdateMember_module.o XLBOMD_module.o atoms.module.o basic_types.o \
	cover_module.o datatypes.module.o density_module.o dimens.module.o \
	ewald.o force_module.o functions_on_grid_module.o generic_blas.o \
	generic_comms.o global_module.o group_module.o input_module.o \
	io_module.o io_module2.o logicals_module.o matrix_data_module.o \
	matrix_module.o maxima_module.o minimise.o mult_module.o \
	numbers.module.o primary_module.o pseudo_tm_module.o \
	pseudopotential.module.o pseudopotential_common.o \
	set_blipgrid_module.o set_bucket_module.o species_module.o \
	timer_module.o timer_stdclocks_module.o units.module.o
Integrators_module.o: datatypes.module.o global_module.o move_atoms.module.o \
	numbers.module.o species_module.o
constraint_module.o: auxiliary_types.o datatypes.module.o generic_comms.o \
	global_module.o numbers.module.o species_module.o
block_module.o: construct_module.o datatypes.module.o generic_comms.o \
	global_module.o grid_index.module.o group_module.o input_module.o \
	maxima_module.o numbers.module.o primary_module.o timer_module.o \
	timer_stdclocks_module.o
comm_array_module.o: datatypes.module.o generic_comms.o global_module.o \
	numbers.module.o
naba_blk_module.o: datatypes.module.o generic_comms.o global_module.o \
	memory_module.o timer_module.o timer_stdclocks_module.o
set_blipgrid_module.o: atomic_density.o atoms.module.o blip_module.o \
	block_module.o comm_array_module.o cover_module.o datatypes.module.o \
	dimens.module.o functions_on_grid_module.o generic_comms.o \
	global_module.o group_module.o maxima_module.o mpi.o \
	naba_blk_module.o numbers.module.o primary_module.o species_module.o \
	timer_module.o timer_stdclocks_module.o
bucket_module.o: datatypes.module.o generic_comms.o global_module.o \
	memory_module.o naba_blk_module.o set_blipgrid_module.o \
	timer_module.o timer_stdclocks_module.o
set_bucket_module.o: bucket_module.o comm_array_module.o cover_module.o \
	datatypes.module.o generic_comms.o global_module.o group_module.o \
	matrix_data_module.o matrix_module.o mpi.o naba_blk_module.o \
	numbers.module.o primary_module.o set_blipgrid_module.o \
	species_module.o timer_module.o timer_stdclocks_module.o
calc_matrix_elements_module.o: block_module.o bucket_module.o \
	comm_array_module.o datatypes.module.o dimens.module.o \
	functions_on_grid_module.o generic_blas.o generic_comms.o \
	global_module.o matrix_data_module.o matrix_module.o memory_module.o \
	mpi.o mult_module.o naba_blk_module.o numbers.module.o \
	primary_module.o set_blipgrid_module.o timer_module.o \
	timer_stdclocks_module.o
blip_grid_transform_module.o: atoms.module.o blip_module.o block_module.o \
	comm_array_module.o cover_module.o datatypes.module.o fft_module.o \
	functions_on_grid_module.o generic_comms.o global_module.o \
	group_module.o memory_module.o mpi.o numbers.module.o \
	primary_module.o set_blipgrid_module.o species_module.o \
	support_spec_format.o timer_module.o timer_stdclocks_module.o
blip_pao_values.o: datatypes.module.o generic_comms.o global_module.o \
	memory_module.o numbers.module.o pao_format.o timer_module.o \
	timer_stdclocks_module.o
linear_equation_solver.o: datatypes.module.o generic_comms.o numbers.module.o
ol_ang_coeff_subs.o: datatypes.module.o generic_comms.o global_module.o \
	numbers.module.o ol_bessel_transforms.o ol_int_datatypes.module.o \
	ol_rad_table_subs.o ol_spline_subs.o pao_format.o spline_module.o \
	timer_module.o timer_stdclocks_module.o
ol_bessel_transforms.o: datatypes.module.o generic_comms.o global_module.o \
	memory_module.o numbers.module.o ol_fft_procs.o
ol_fft_procs.o: datatypes.module.o numbers.module.o
ol_int_datatypes.module.o: datatypes.module.o
ol_rad_table_subs.o: datatypes.module.o generic_comms.o global_module.o \
	input_module.o memory_module.o numbers.module.o \
	ol_bessel_transforms.o ol_int_datatypes.module.o ol_spline_subs.o \
	pao_format.o primary_module.o pseudo_tm_info.o species_module.o \
	support_spec_format.o timer_module.o timer_stdclocks_module.o
ol_spline_subs.o: datatypes.module.o generic_comms.o global_module.o \
	numbers.module.o
pao2blip.o: blip_module.o blip_pao_values.o datatypes.module.o \
	dimens.module.o generic_comms.o global_module.o \
	linear_equation_solver.o maxima_module.o memory_module.o \
	numbers.module.o pao_format.o primary_module.o read_support_spec.o \
	species_module.o support_spec_format.o symmetry.o timer_module.o \
	timer_stdclocks_module.o
pao_minimisation.module.o: DMMinModule.o DiagModule.o H_matrix_module.o \
	PAO_grid_transform_module.o PosTan_module.o Pulay.o S_matrix_module.o \
	SelfCon_module.o cover_module.o datatypes.module.o density_module.o \
	energy.module.o functions_on_grid_module.o generic_blas.o \
	generic_comms.o global_module.o group_module.o logicals_module.o \
	matrix_data_module.o matrix_module.o maxima_module.o memory_module.o \
	mult_module.o numbers.module.o ol_rad_table_subs.o primary_module.o \
	species_module.o support_spec_format.o timer_module.o \
	timer_stdclocks_module.o
PAO_grid_transform_module.o: blip_grid_transform_module.o block_module.o \
	comm_array_module.o cover_module.o datatypes.module.o dimens.module.o \
	functions_on_grid_module.o generic_comms.o global_module.o \
	group_module.o numbers.module.o ol_ang_coeff_subs.o primary_module.o \
	set_blipgrid_module.o species_module.o support_spec_format.o \
	timer_module.o timer_stdclocks_module.o
build_PAO_matrices.module.o: cover_module.o datatypes.module.o \
	generic_comms.o global_module.o group_module.o matrix_data_module.o \
	matrix_module.o mult_module.o numbers.module.o ol_ang_coeff_subs.o \
	primary_module.o pseudo_tm_info.o species_module.o \
	support_spec_format.o timer_module.o timer_stdclocks_module.o
read_support_spec.o: datatypes.module.o generic_comms.o global_module.o \
	input_module.o maxima_module.o numbers.module.o species_module.o
support_spec_format.o: datatypes.module.o generic_comms.o global_module.o \
	memory_module.o timer_module.o timer_stdclocks_module.o
symmetry.o: datatypes.module.o generic_comms.o numbers.module.o
pseudo_tm_module.o: block_module.o cover_module.o datatypes.module.o \
	dimens.module.o ewald.o functions_on_grid_module.o generic_blas.o \
	generic_comms.o global_module.o group_module.o hartree.o \
	maxima_module.o numbers.module.o ol_ang_coeff_subs.o primary_module.o \
	pseudo_tm_info.o pseudopotential_common.o set_blipgrid_module.o \
	species_module.o timer_module.o timer_stdclocks_module.o
pseudo_tm_info.o: datatypes.module.o dimens.module.o ewald.o generic_comms.o \
	global_module.o input_module.o memory_module.o numbers.module.o \
	pao_format.o pseudopotential_common.o species_module.o \
	spline_module.o timer_module.o timer_stdclocks_module.o
nlpf2blip.o: blip_grid_transform_module.o blip_module.o datatypes.module.o \
	dimens.module.o generic_blas.o generic_comms.o global_module.o \
	io_module.o mult_module.o numbers.module.o pseudo_tm_info.o \
	species_module.o support_spec_format.o timer_module.o \
	timer_stdclocks_module.o
pseudopotential_common.o: datatypes.module.o
exx_module.o: DiagModule.o PAO_grid_transform_module.o atomic_density.o \
	cover_module.o datatypes.module.o dimens.module.o exx_io.o \
	exx_types.o fft_ffte.o generic_comms.o global_module.o group_module.o \
	matrix_data_module.o numbers.module.o primary_module.o \
	species_module.o support_spec_format.o units.module.o
exx_kernel_default.o: PAO_grid_transform_module.o basic_types.o \
	comms_module.o cover_module.o datatypes.module.o energy.module.o \
	exx_evalpao.o exx_io.o exx_memory.o exx_module.o exx_types.o \
	generic_comms.o global_module.o group_module.o matrix_comms_module.o \
	matrix_data_module.o matrix_module.o memory_module.o mpi.o \
	mult_module.o multiply_module.o numbers.module.o primary_module.o \
	species_module.o timer_module.o timer_stdclocks_module.o
exx_types.o: datatypes.module.o timer_module.o timer_stdclocks_module.o
exx_memory.o: datatypes.module.o exx_types.o fft_ffte.o generic_comms.o \
	global_module.o memory_module.o numbers.module.o timer_module.o
exx_io.o: atomic_density.o cover_module.o datatypes.module.o exx_types.o \
	generic_comms.o global_module.o group_module.o io_module.o \
	matrix_data_module.o numbers.module.o primary_module.o \
	species_module.o support_spec_format.o timer_module.o units.module.o
exx_evalpao.o: datatypes.module.o dimens.module.o exx_types.o \
	functions_on_grid_module.o generic_comms.o numbers.module.o \
	ol_ang_coeff_subs.o pao_format.o support_spec_format.o timer_module.o
