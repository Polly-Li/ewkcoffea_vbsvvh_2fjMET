#scripts
working_dir = '/home/users/pyli/analysis/vvh/'
coffea_scripts_dir = '/home/users/pyli/analysis/vvh/coffea_scripts/'
cutflow_yamls_dir = '/home/users/pyli/analysis/vvh/extra/configs/cutflows/'

#files
default_output_dir = '/data/userdata/pyli/projects/VVHjj/outputs/'

#var stuffs
basic_cutflow = f'{cutflow_yamls_dir}/MET_basic.yaml'
from extra.configs.config_handling import get_cutflow
objsel_cf = get_cutflow(basic_cutflow,'MET_objsel')
default_cutflow_yaml = f'{cutflow_yamls_dir}/all.yaml'