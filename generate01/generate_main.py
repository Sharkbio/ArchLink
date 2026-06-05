from .generate_folder import *
from .generate_frequency import *


def generate_init(args):
    binning_dir = args.output_path
    num_threads = args.num_threads
    output_file = args.output_path + '/binning/frequency/frequency.pkl'
    # generate folder
    process_selected_binning_methods(binning_dir)
    # generate frequency
    result_files_pattern = os.path.join(binning_dir+'/binning/frequency/', "*", "result")
    input_files = glob.glob(result_files_pattern)
    # print(result_files_pattern)

    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    process_files_and_generate_dict(input_files, output_file, num_threads)
    

if __name__ == "__main__":
    generate_init()