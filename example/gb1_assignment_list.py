import os
import nmrassign.analysis

# Load the data
data_dir = '/Users/kjf/git/nmrassign/tests/data/gb1_1'
control_file = os.path.join(data_dir, 'control_nsga.txt')
with open(control_file, 'r') as control_fid:
    params = nmrassign.fileio.read_control(control_fid, parse_input=True)

output_file = os.path.join(data_dir, 'outdata.txt')
with open(output_file, 'r') as output_fid:
    spectra_assignments = nmrassign.fileio.read_outdata(params, output_fid)


spectrum_assign_list = nmrassign.analysis.assignment_list(spectra_assignments)


for res_assignment in zip(*spectrum_assign_list):
    print(res_assignment)
