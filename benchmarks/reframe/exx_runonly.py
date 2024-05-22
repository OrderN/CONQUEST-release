# Import modules from reframe and excalibur-tests
import reframe as rfm
import reframe.utility.sanity as sn

class ConquestBaseBenchmark(rfm.RunOnlyRegressionTest):

    # Run configuration
    ## Mandatory ReFrame setup
    valid_systems = ['-gpu']
    valid_prog_environs = ['default']

    ## Executable
    executable_opts = ['']

    ## Scheduler options
    time_limit = '30m'

    conquest_base_dir = variable(str, value='foo')

    @run_after('setup')
    def setup_variables(self):
        self.executable = f"{self.conquest_base_dir}/bin/Conquest"
        self.num_tasks = self.num_tasks_param
        self.num_cpus_per_task = self.num_cpus_per_task_param
        self.env_vars['OMP_NUM_THREADS'] = f'{self.num_cpus_per_task}'
        self.env_vars['OMP_STACKSIZE'] = '100M'

        if self.current_partition.scheduler.registered_name == 'sge':
            # `self.num_tasks` or `self.num_cpus_per_task` may be `None`, here
            # we default to `1` if not set.
            num_tasks = self.num_tasks or 1
            num_cpus_per_task = self.num_cpus_per_task or 1
            # Set the total number of CPUs to be requested for the SGE scheduler.
            # Set to a full node size to reduce runtime variance.
            self.extra_resources['mpi'] = {'num_slots': self.current_partition.processor.num_cpus} #num_tasks * num_cpus_per_task}

    @sanity_function
    def validate_solution(self):
        return sn.assert_found(r'Total run time was:.*', self.stdout)

    @performance_function('s', perf_key='Runtime')
    def extract_runtime_perf(self):
        return sn.extractsingle(r'Total run time was:\s+(\S+)\s+seconds', self.stdout, 1, float)
    
    @performance_function('s', perf_key='exx_setup_runtime')
    def extract_exx_setup_runtime_perf(self):
        return sn.extractsingle(r'Time spent in exx_setup\s+time:\s+=\s+(\S+)\ss', self.stdout, 1, float)
    
    @performance_function('s', perf_key='exx_fetch_runtime')
    def extract_exx_fetch_runtime_perf(self):
        return sn.extractsingle(r'Time spent in exx_fetch\s+time:\s+=\s+(\S+)\ss', self.stdout, 1, float)
    
    @performance_function('s', perf_key='exx_evalpao_runtime')
    def extract_exx_evalpao_runtime_perf(self):
        return sn.extractsingle(r'Time spent in exx_evalpao\s+time:\s+=\s+(\S+)\ss', self.stdout, 1, float)
    
    @performance_function('s', perf_key='exx_allocat_runtime')
    def extract_exx_allocat_runtime_perf(self):
        return sn.extractsingle(r'Time spent in exx_allocat\s+time:\s+=\s+(\S+)\ss', self.stdout, 1, float)
    
    @performance_function('s', perf_key='exx_dealloc_runtime')
    def extract_exx_dealloc_runtime_perf(self):
        return sn.extractsingle(r'Time spent in exx_dealloc\s+time:\s+=\s+(\S+)\ss', self.stdout, 1, float)
    
    @performance_function('s', perf_key='exx_kernel_runtime')
    def extract_exx_kernel_runtime_perf(self):
        return sn.extractsingle(r'Time spent in exx_kernel\s+time:\s+=\s+(\S+)\ss', self.stdout, 1, float)
    
    @performance_function('s', perf_key='exx_total_runtime')
    def extract_exx_total_runtime_perf(self):
        return sn.extractsingle(r'Time spent in exx_total\s+time:\s+=\s+(\S+)\ss', self.stdout, 1, float)

    @performance_function('MB', perf_key='Memory')
    def extract_memory_perf(self):
        return sn.extractsingle(r'Max total mem use is\s+(\S+)\s+MB', self.stdout, 1, float)

@rfm.simple_test
class test_EXX_isol_C2H4_4proc_PBE0CRI_fullDZP_0_2(ConquestBaseBenchmark):

    tags = {"test_EXX_isol_C2H4_4proc_PBE0CRI_fullDZP_0.2"}
    num_tasks_param = parameter([2,4,8])
    num_cpus_per_task_param = parameter([1,2,4,8])

    @run_before('run')
    def get_input(self):
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0CRI_fullDZP_0.2/C_PBE_DZP_CQ.ion .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0CRI_fullDZP_0.2/Conquest_coord .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0CRI_fullDZP_0.2/Conquest_input .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0CRI_fullDZP_0.2/H_PBE_DZP_CQ.ion .")

@rfm.simple_test
class test_EXX_isol_C2H4_4proc_PBE0CRI_fullDZP_0_4(ConquestBaseBenchmark):

    tags = {"test_EXX_isol_C2H4_4proc_PBE0CRI_fullDZP_0.4"}
    num_tasks_param = parameter([2,4,8])
    num_cpus_per_task_param = parameter([1,2,4,8])

    @run_before('run')
    def get_input(self):
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0CRI_fullDZP_0.4/C_PBE_DZP_CQ.ion .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0CRI_fullDZP_0.4/Conquest_coord .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0CRI_fullDZP_0.4/Conquest_input .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0CRI_fullDZP_0.4/H_PBE_DZP_CQ.ion .")

@rfm.simple_test
class test_EXX_isol_C2H4_4proc_PBE0CRI_fullDZP_0_6(ConquestBaseBenchmark):

    tags = {"test_EXX_isol_C2H4_4proc_PBE0CRI_fullDZP_0.6"}
    num_tasks_param = parameter([2,4,8])
    num_cpus_per_task_param = parameter([1,2,4,8])

    @run_before('run')
    def get_input(self):
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0CRI_fullDZP_0.6/C_PBE_DZP_CQ.ion .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0CRI_fullDZP_0.6/Conquest_coord .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0CRI_fullDZP_0.6/Conquest_input .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0CRI_fullDZP_0.6/H_PBE_DZP_CQ.ion .")

@rfm.simple_test
class test_EXX_isol_C2H4_4proc_PBE0CRI_fullSZP_0_4_SCF(ConquestBaseBenchmark):

    tags = {"test_EXX_isol_C2H4_4proc_PBE0CRI_fullSZP_0.4_SCF"}
    num_tasks_param = parameter([2,4,8])
    num_cpus_per_task_param = parameter([1,2,4,8])

    @run_before('run')
    def get_input(self):
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0CRI_fullSZP_0.4_SCF/C_PBE_SZP_CQ.ion .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0CRI_fullSZP_0.4_SCF/Conquest_coord .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0CRI_fullSZP_0.4_SCF/Conquest_input .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0CRI_fullSZP_0.4_SCF/H_PBE_SZP_CQ.ion .")

@rfm.simple_test
class test_EXX_isol_C2H4_4proc_PBE0CRI_fullTZTP_0_6(ConquestBaseBenchmark):

    tags = {"test_EXX_isol_C2H4_4proc_PBE0CRI_fullTZTP_0.6"}
    num_tasks_param = parameter([2,4,8])
    num_cpus_per_task_param = parameter([1,2,4,8])

    @run_before('run')
    def get_input(self):
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0CRI_fullTZTP_0.6/C_PBE_TZTP_CQ.ion .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0CRI_fullTZTP_0.6/Conquest_coord .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0CRI_fullTZTP_0.6/Conquest_input .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0CRI_fullTZTP_0.6/H_PBE_TZTP_CQ.ion .")

@rfm.simple_test
class test_EXX_isol_C2H4_4proc_PBE0ERI_fullSZP_0_4_SCF(ConquestBaseBenchmark):

    tags = {"test_EXX_isol_C2H4_4proc_PBE0ERI_fullSZP_0.4_SCF"}
    num_tasks_param = parameter([2,4,8])
    num_cpus_per_task_param = parameter([1,2,4,8])

    @run_before('run')
    def get_input(self):
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0ERI_fullSZP_0.4_SCF/C_PBE_SZP_CQ.ion .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0ERI_fullSZP_0.4_SCF/Conquest_coord .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0ERI_fullSZP_0.4_SCF/Conquest_input .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0ERI_fullSZP_0.4_SCF/H_PBE_SZP_CQ.ion .")

@rfm.simple_test
class test_EXX_isol_C2H4_4proc_PBE0ERI_fullSZP_GTO_SCF(ConquestBaseBenchmark):

    tags = {"test_EXX_isol_C2H4_4proc_PBE0ERI_fullSZP_GTO_SCF"}
    num_tasks_param = parameter([2,4,8])
    num_cpus_per_task_param = parameter([1,2,4,8])

    @run_before('run')
    def get_input(self):
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0ERI_fullSZP_GTO_SCF/C_PBE_SZP_CQ.ion .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0ERI_fullSZP_GTO_SCF/C.gto .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0ERI_fullSZP_GTO_SCF/Conquest_coord .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0ERI_fullSZP_GTO_SCF/Conquest_input .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0ERI_fullSZP_GTO_SCF/H_PBE_SZP_CQ.ion .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C2H4_4proc_PBE0ERI_fullSZP_GTO_SCF/H.gto .")

@rfm.simple_test
class test_EXX_isol_C6H6_4proc_PBE0CRI_fullDZP_0_6(ConquestBaseBenchmark):

    tags = {"test_EXX_isol_C6H6_4proc_PBE0CRI_fullDZP_0.6"}
    num_tasks_param = parameter([2,4,8])
    num_cpus_per_task_param = parameter([1,2,4,8])

    @run_before('run')
    def get_input(self):
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C6H6_4proc_PBE0CRI_fullDZP_0.6/C_PBE_DZP_CQ.ion .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C6H6_4proc_PBE0CRI_fullDZP_0.6/Conquest_coord .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C6H6_4proc_PBE0CRI_fullDZP_0.6/Conquest_input .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C6H6_4proc_PBE0CRI_fullDZP_0.6/H_PBE_DZP_CQ.ion .")

@rfm.simple_test
class test_EXX_isol_C6H6_4proc_PBE0CRI_fullSZP_0_6(ConquestBaseBenchmark):

    tags = {"test_EXX_isol_C6H6_4proc_PBE0CRI_fullSZP_0.6"}
    num_tasks_param = parameter([2,4,8])
    num_cpus_per_task_param = parameter([1,2,4,8])

    @run_before('run')
    def get_input(self):
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C6H6_4proc_PBE0CRI_fullSZP_0.6/C_PBE_SZP_CQ.ion .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C6H6_4proc_PBE0CRI_fullSZP_0.6/Conquest_coord .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C6H6_4proc_PBE0CRI_fullSZP_0.6/Conquest_input .")
        self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/test_EXX_isol_C6H6_4proc_PBE0CRI_fullSZP_0.6/H_PBE_SZP_CQ.ion .")


# @rfm.simple_test
# class Water64(ConquestBaseBenchmark):

#     tags = {"water64"}
#     num_tasks_param = 8
#     num_cpus_per_task = 4

#     @run_before('run')
#     def get_input(self):
#         self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/water_64mols/Conquest_input .")
#         self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/water_64mols/H2O_coord.in .")
#         self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/water_64mols/H_SZ.ion .")
#         self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/water_64mols/H_SZP.ion .")
#         self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/water_64mols/O_SZ.ion .")
#         self.prerun_cmds.append(f"cp {self.conquest_base_dir}/benchmarks/water_64mols/O_SZP.ion .")
