"""
Functions for reading and writing files for the NMR GA-Assign program.
"""
import re
import nmrassign.base


def read_sequence(seq_file):
    """
    Read sequence file.

    :param seq_file: sequence file name with path. File should have a list
        of 1-letter amino acid codes. FASTA type files are also supported but
        if more then one sequence is in the file an error will be thrown.
    :return: residue_sequence
    :rtype: list
    """
    sequence_count = 0
    residue_sequence = []
    for line in seq_file:
        line = list(line.strip())

        # if the line is empty
        if not line:
            continue

        # if the line is a fast header line.
        if line[0] == '>':
            sequence_count += 1

            # if there has been more then one fasta header file
            if sequence_count > 1:
                mesg = '{} contains more then one seq.'.format(seq_file)
                raise ValueError(mesg)
            # there is not sequence info on this line skip it.
            continue

        residue_sequence.extend(line)
    return residue_sequence


def write_sequence(seq, seq_file):
    """
    Write sequence to file.

    :seq: iterable of 1-letter amino acids codes.
    :param seq_file: file object
    """
    if seq_file.name[-6:] == '.fasta':
        seq_file.write('>\n')
    seq_file.write(''.join(seq))


def read_spectrum(spectrum_file):
    """
    Read spectrum file.

    :param spectrum_file: file name with path.
    """
    # Read the header line to get the number of peak and frequencies
    line = next(spectrum_file)
    num_peak, num_freq = map(int, line.split())

    # Read the rest of the lines.
    peak_groups = []

    # regex expression to match (and remove) things in parentheses
    parentheses = re.compile('\(.+?\)')

    for line in spectrum_file:
        line = line.split()

        # parse the shifts and lines.
        shifts = map(float, line[0:num_freq])
        widths = map(float, line[num_freq: num_freq * 2])

        # group the shifts and line and parse the rest of the line
        peaks = [(s, w) for (s, w) in zip(shifts, widths)]
        peaks = [nmrassign.base.Peak(*x) for x in peaks]
        num = line[-2]
        assignments = parentheses.sub('', line[-1])

        peak_group = nmrassign.base.PeakGroup(peaks, assignments, num)
        peak_groups.append(peak_group)
    return peak_groups


def write_spectrum(peak_groups, spectrum_file, deliminator=' '):
    """
    Write spectrum file.

    :param peak_groups: list of base.PeakGroup
    :param spectrum_file: file name with path.
    :param deliminator: delimiter default
    """
    # Read the header line to get the number of peak and frequencies
    num_group = len(peak_groups)
    num_freq = len(peak_groups[0].peaks)

    spectrum_file.write('{} {}\n'.format(num_group, num_freq))

    for peak_group in peak_groups:
        spectrum_file.write(peak_group.spectrum_line_str(deliminator))
        spectrum_file.write('\n')


def read_connection(connection_file):
    """
    Read connection file.

    File Description:
    The first line is count of atom lines. Each of the other lines has the
    connection information for one atom. There are three piece of information
    in each line.. The length of the line divided by 3.0 is the number of
    spectra (n) that the atoms is in. The line[:n] is a list of spectra that
    the atom is in. line[n:n * 2] is the atom index in the respective spectra.
    line[-n:] says how the atoms connect to the surrounding residues. If 0 the
    the convention is intra-residue.

    :param connection_file: file name with path.
    """
    # Read the header line to get the number of peak and frequencies

    next(connection_file)
    # line = next(fid)
    # number of peaks = int(line.strip())

    atoms = []
    # Read the rest of the lines.
    for line in connection_file:
        line = list(map(int, line.split()))

        # The number of spectra in the atoms line.
        n = int(len(line) / 3.0)
        spectra = [x - 1 for x in line[:n]]
        ids = [x - 1 for x in line[n:n * 2]]
        residue = line[-n:]
        atoms.append(nmrassign.base.Connections(spectra, ids, residue))
    return atoms


def write_connection(atom_connections, connections_file, delimiter=' '):
    """
    Write connections file.

    :param atom_connections: list of base.PeakGroup
    :param connections_file: file name with path.
    """
    # Read the header line to get the number of peak and frequencies
    num_atoms = len(atom_connections)

    connections_file.write('{}\n'.format(num_atoms))

    for connects in atom_connections:
        # convert back to non-zero indexing
        line = [x + 1 for x in connects.spectra]
        line += [x + 1 for x in connects.ids]
        line += connects.residues
        connections_file.write(delimiter.join(map(str, line)))
        connections_file.write('\n')
