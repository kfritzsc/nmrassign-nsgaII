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


# Analyze the data
scores = nmrassign.fileio.read_outtab_scores(params)
top_paretos = nmrassign.analysis.keep_pareto_order(scores)
pareto_assigns = nmrassign.analysis.pareto_filter(
    spectra_assignments, top_paretos)

spectrum_assign_list = nmrassign.analysis.assignment_list(pareto_assigns)

with open('gb1_output_list.txt', 'w') as fid:

    for n, res_assignment in enumerate(zip(*spectrum_assign_list)):

        nca = ', '.join(map(str, res_assignment[0][0]))
        nco = ', '.join(map(str, res_assignment[1][0]))

        score_nca = ', '.join(map(str, res_assignment[0][1]))
        score_nco = ', '.join(map(str, res_assignment[1][1]))

        line = '\t'.join(map(str, [n+1, nca, nco, score_nca, score_nco]))
        line +='\n'

        fid.write(line)



