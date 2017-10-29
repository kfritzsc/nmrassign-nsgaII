import os
import numpy as np
import matplotlib.pyplot as plt

import nmrassign.fileio
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
mode_assignments, mode_probabilities = nmrassign.analysis.assignment_stats(
    pareto_assigns)
mode_probabilities = np.array(mode_probabilities)

# Plot the probabilities
fig = plt.figure()
for k in range(len(mode_probabilities)):
    ax = plt.subplot(1, 2, k + 1)
    plt.plot(mode_probabilities[k], 'o--', color='green')

    title = os.path.basename(params['spectra_files'][k])
    title = os.path.splitext(title)[0]
    plt.title(title)
    plt.ylabel('Probability (%)')
    plt.xlabel('Residue number')

    plt.xlim(0, len(mode_probabilities[k]))
    plt.ylim([-0.1, 1.1])

plt.show()
