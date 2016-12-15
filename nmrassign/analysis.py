"""
Tools for analysing the results of GA/MC Assignment programs.
"""


def most_likely_assignment(params, assignments):
    pass

if __name__ == '__main__':

    import os
    from nmrassign.fileio import read_control, read_outdata

    data_dir = '/Users/kjf/git/optassign/tests/data/gb1_1/'
    control_file = os.path.join(data_dir, 'control_nsga.txt')
    with open(control_file, 'r') as control_fid:
        params = read_control(control_fid, parse_input=True)

    output_file = os.path.join(data_dir, 'outdata.txt')
    with open(output_file, 'r') as output_fid:
        spectra_assignments = read_outdata(params, output_fid)

    print(spectra_assignments)