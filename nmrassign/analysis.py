"""
Tools for analysing the results of GA/MC Assignment programs.
"""


def keep_pareto_order(scores, order=1):
    """
    Returns a bool array of len == group size. The value is true if
    the the pareto_order is less then order

    :param scores: list of `base.scores`
    :param order: The highest pareto order to be kept
    :return:
    """
    top_paretos = []
    for score in scores:
        if score.pareto <= order:
            top_paretos.append(True)
        else:
            top_paretos.append(False)

    return top_paretos


def pareto_filter(assignmets, top_paretos):
    """
    Remove assignment groups that are dominated.

    :param assignmets: n_specta x group_size x len(seq) array
    :param top_paretos: bool array of len group_size
    """
    new_assignments = []
    for spectra_assigns in assignmets:
        new_spectrum_assigns = []
        for assigns, top_pareto in zip(spectra_assigns, top_paretos):
            if top_pareto:
                new_spectrum_assigns.append(assigns)
        new_assignments.append(new_spectrum_assigns)

    return new_assignments





if __name__ == '__main__':

    import os
    from nmrassign.fileio import read_control, read_outdata

    data_dir = '/Users/kjf/git/nmrassign/tests/data/gb1_2/'
    control_file = os.path.join(data_dir, 'control_nsga.txt')
    with open(control_file, 'r') as control_fid:
        params = read_control(control_fid, parse_input=True)

    output_file = os.path.join(data_dir, 'outdata.txt')
    with open(output_file, 'r') as output_fid:
        spectra_assignments = read_outdata(params, output_fid)

    print(spectra_assignments)
