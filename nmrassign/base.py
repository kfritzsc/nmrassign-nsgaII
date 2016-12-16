"""
Common NMR assignment objects.
"""
from collections import namedtuple

Connections = namedtuple('Connections',
                         ['spectra', 'ids', 'residues'])


Peak = namedtuple('Peak', ['poss', 'width'])


Score = namedtuple('Score', ['ng', 'nb', 'ne', 'nu', 'pareto'])


class PeakGroup(object):
    """
    A group of NMR peaks with ambiguity and list of possible
    assignments.

    :param peaks: list of peaks in group
    :param assignments: list of possible assignments
    :param degeneracy: number of possible assignments, normally 1
    :param num_used: number of of times the peak group cab be used
    """
    def __init__(self, peaks, assignments, degeneracy=1, num_used=0):
        self.peaks = peaks
        self.res_poss = assignments
        self.degeneracy = degeneracy
        self.num_used = num_used

    def __str__(self):
        return '{}, {}, {}'.format(self.peaks, self.res_poss,
                                   self.num)

    def spectrum_line_str(self, delimiter=' '):
        line = [x.poss for x in self.peaks]
        line += [x.width for x in self.peaks]
        line += [self.degeneracy, ''.join(self.res_poss)]
        return delimiter.join(map(str, line))
