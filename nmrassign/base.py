"""
Common NMR assignment objects.
"""
from collections import namedtuple


Connections = namedtuple('Connections', ['spectra', 'ids', 'residues'])


Peak = namedtuple('Peak', ['poss', 'width'])


class PeakGroup(object):
    """
    A group of NMR peaks with ambiguity and list of possible assignments.

    :param peaks: list of peaks in group
    :param num: number of possible assignments, normally 1
    :param assignments: list of possible assignments
    """
    def __init__(self, peaks, assignments, degeneracy=1, num_used=0):
        self.peaks = peaks
        self.degeneracy = degeneracy
        self.assignments = assignments
        self.num_used = num_used

    def __str__(self):
        return '{}, {}, {}'.format(self.peaks, self.assignments, self.num)

    def spectrum_line_str(self, delimiter=' '):
        line = [x.poss for x in self.peaks]
        line += [x.width for x in self.peaks]
        line += [self.degeneracy, ''.join(self.assignments)]
        return delimiter.join(map(str, line))


# class Assign(object):
#     """
#     Assign spectrum tables with connection to a sequence.
#
#     :param data: lists of list of PeakGroups
#     :param connections: list of Connections
#     :param sequence: protein sequence
#     :param current: current assignment
#     :param random_seed: set the random_seed
#     """
#     def __init__(self, data, connections, sequence, current=None,
#                  random_seed=None):
#         self.seq = sequence
#         self.n = len(sequence)
#
#         self.data = data
#         self.connect = connections
#         if current is None:
#             self.current = self.random_guess()
#         self.scores = self.scores()
#
#         if random_seed is not None:
#             np.random.seed(random_seed)
#
#     def random_guess(self):
#         raise NotImplementedError
#
#     def scores(self):
#         raise NotImplementedError
