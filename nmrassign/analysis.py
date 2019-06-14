"""
Tools for analysing the results of GA/MC Assignment programs.
"""
from collections import Counter


def keep_pareto_order(scores, order=1):
    """
    Returns a bool array of len group size. The value is True if
    the the Pareto order is less then order

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


def pareto_filter(assignments, top_paretos):
    """
    Remove assignment groups that are dominated.

    :param assignments: n_specta x group_size x len(seq) array
    :param top_paretos: bool array of len group_size
    """
    new_assignments = []
    for spectra_assigns in assignments:
        new_spectrum_assigns = []
        for assigns, top_pareto in zip(spectra_assigns, top_paretos):
            if top_pareto:
                new_spectrum_assigns.append(assigns)
        new_assignments.append(new_spectrum_assigns)

    return new_assignments


def assignment_stats(assignments):
    """
    Find the most probably assignment and the probability
    :param assignments:
    """

    mode_assignments = []
    mode_probs = []

    for spectrum_assign in assignments:
        mode_spectrum_assign = []
        mode_spectrum_prob = []
        for res_assign in zip(*spectrum_assign):
            assignment, assignment_count = \
                Counter(res_assign).most_common(1)[0]

            mode_spectrum_assign.append(assignment)
            mode_spectrum_prob.append(assignment_count)

        n = float(len(res_assign))
        mode_spectrum_prob = [x/n for x in mode_spectrum_prob]

        mode_assignments.append(mode_spectrum_assign)
        mode_probs.append(mode_spectrum_prob)

    return mode_assignments, mode_probs


def assignment_list(assignments):
    spectrum_assign_list = []

    for spectrum_assign in assignments:
        assignments_list = []
        scores_list = []
        for res_assign in zip(*spectrum_assign):
            count = Counter(res_assign)

            # Grab the key and counts
            assignments = tuple(count.keys())
            scores = tuple(count[x] for x in assignments)

            assignments_list.append(assignments)
            scores_list.append(scores)

        spectrum_assign_list.append(zip(assignments_list, scores_list))
    return spectrum_assign_list