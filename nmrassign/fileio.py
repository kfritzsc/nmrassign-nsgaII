"""
Functions for reading and writing files for the NMR GA/MC Assignment
program.
"""
import os
import glob
import re
import itertools

import nmrassign.base


def grouper(n, iterable, fillvalue=None):
    """
    Grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx

    :param n:
    :param iterable:
    :param fillvalue:
    """
    args = [iter(iterable)] * n
    return itertools.zip_longest(fillvalue=fillvalue, *args)


def read_sequence(seq_file):
    """
    Read sequence file.

    :param seq_file: sequence file name with path. File should have a
        list of 1-letter amino acid codes. FASTA type files are also
        supported but if more then one sequence is in the file an
        error will be thrown.
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
                mesg = '{} has more then one seq.'.format(seq_file)
                raise ValueError(mesg)
            # there is not sequence info on this line skip it.
            continue

        residue_sequence.extend(line)
    return residue_sequence


def write_sequence(seq, seq_file):
    """
    Write sequence to file.

    :param seq: iterable of 1-letter amino acids codes.
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
    remove_parentheses = re.compile(r'\(.+?\)')

    for line in spectrum_file:
        line = line.split()

        # parse the shifts and lines.
        shifts = map(float, line[0:num_freq])
        widths = map(float, line[num_freq: num_freq * 2])

        # group the shifts and line and parse the rest of the line
        peaks = [(s, w) for (s, w) in zip(shifts, widths)]
        peaks = [nmrassign.base.Peak(*x) for x in peaks]
        num = line[-2]
        assignments = remove_parentheses.sub('', line[-1])

        peak_group = nmrassign.base.PeakGroup(peaks, assignments,
                                              num)
        peak_groups.append(peak_group)
    return peak_groups


def write_spectrum(peak_groups, spectrum_file, deliminator=' '):
    """
    Write spectrum file.

    :param peak_groups: list of `base.PeakGroup`
    :param spectrum_file: file name with path.
    :param deliminator: delimiter default
    """
    # Read the header line to get the number of peak and frequencies
    num_group = len(peak_groups)
    num_freq = len(peak_groups[0].peaks)

    spectrum_file.write('{} {}\n'.format(num_group, num_freq))

    for peak_group in peak_groups:
        spectrum_file.write(
            peak_group.spectrum_line_str(deliminator))
        spectrum_file.write('\n')


def read_connection(connection_file):
    """
    Read connection file.

    File Description:
    The first line is count of atom lines. Each of the other lines
    has the connection information for one atom. There are three
    piece of information in each line. The length of the line divided
    by 3 is the number of spectra (n) that the atom is in. The
    line[:n] is a list of spectra that the atom is in. line[n:n * 2]
    is the atom index in the respective spectra. line[-n:] says how
    the atoms connect to the surrounding residues. If 0 the
    connection is intraresidue.

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
        atoms.append(nmrassign.base.Connections(spectra, ids,
                                                residue))
    return atoms


def write_connection(atom_connections, connections_file,
                     delimiter=' '):
    """
    Write connections file.

    :param atom_connections: list of base.PeakGroup
    :param connections_file: file name with path.
    :param delimiter: type of delimiter
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


def read_control(control_file, parse_input=True):
    """
    Read control file

    :param control_file: file name with path.
    :param parse_input: Bool, if True parse input file.
    """

    params = dict()

    params['method'] = next(control_file).split()[0]

    next(control_file)
    params['seq_file'] = next(control_file).split()[0]

    if parse_input:
        with open(params['seq_file'], 'r') as fid:
            params['seq'] = read_sequence(fid)

    params['n_spectra'] = int(next(control_file).split()[0])

    # There are n_spectra spectrum files.
    spectra_files = []
    for _ in range(params['n_spectra']):
        spectra_files.append(next(control_file).split()[0])
    params['spectra_files'] = spectra_files

    if parse_input:
        spectra_peak_groups = []
        for spectrum_file in params['spectra_files']:
            with open(spectrum_file, 'r') as fid:
                peak_groups = read_spectrum(fid)
                spectra_peak_groups.append(peak_groups)
        params['spectra_peak_groups'] = spectra_peak_groups

    params['connection_file'] = next(control_file).split()[0]
    if parse_input:
        with open(params['connection_file'], 'r') as fid:
            params['connections'] = read_connection(fid)

    params['initial_file'] = next(control_file).split()[0]
    if params['initial_file'].lower() in {'null', 'none'}:
        params['initial_file'] = None
    # TODO: Parse initial input file.

    params['output_file'] = next(control_file).split()[0]
    params['outtab_folder'] = next(control_file).split()[0]

    # read all the algorithm parameters
    next(control_file)

    if params['method'].lower() == 'nsga2_assign':
        params['group_size'] = int(next(control_file).split()[0])
        params['pool_size'] = int(next(control_file).split()[0])
        params['steps'] = int(next(control_file).split()[0])
        params['nsga_attempts'] = int(next(control_file).split()[0])
        params['free_steps'] = int(next(control_file).split()[0])
        params['mutation_rate'] = float(
            next(control_file).split()[0])
        params['additional_rate'] = float(
            next(control_file).split()[0])
        params['crossover_rate'] = float(
            next(control_file).split()[0])
        params['null_prob'] = float(next(control_file).split()[0])
        line = next(control_file).split()[0]
        if line[0] == '!':
            params['random_seed'] = None
        else:
            params['random_seed'] = int(line)

    elif params['method'].lower() == 'nsga2_mc_assign':
        params['group_size'] = int(next(control_file).split()[0])
        params['pool_size'] = int(next(control_file).split()[0])
        params['steps'] = int(next(control_file).split()[0])
        params['nsga_attempts'] = int(next(control_file).split()[0])
        params['mc_attempts'] = int(next(control_file).split()[0])
        params['free_steps'] = int(next(control_file).split()[0])
        params['mutation_rate'] = float(
            next(control_file).split()[0])
        params['additional_rate'] = float(
            next(control_file).split()[0])
        params['crossover_rate'] = float(
            next(control_file).split()[0])
        params['null_prob'] = float(next(control_file).split()[0])
        params['w1_min'] = float(next(control_file).split()[0])
        params['w2_min'] = float(next(control_file).split()[0])
        params['w3_min'] = float(next(control_file).split()[0])
        params['w4_min'] = float(next(control_file).split()[0])
        params['w1_max'] = float(next(control_file).split()[0])
        params['w2_max'] = float(next(control_file).split()[0])
        params['w3_max'] = float(next(control_file).split()[0])
        params['w4_max'] = float(next(control_file).split()[0])
        line = next(control_file).split()[0]
        if line != '!':
            params['random_seed'] = line
        else:
            params['random_seed'] = None
    else:
        msg = 'Method {} not NSGA2_ASSIGN or NSGA2_MC_ASSIGN'.format(
            params['method'])
        raise ValueError(msg)
    return params


def write_control(params, control_file, write_input=False):
    """
    Write the control file from the params dictionary. Optionally,
    write all input files.

    :param control_file: file name with path.
    :param params: dictionary with all info to write to file(s)
    :param write_input: default False, write all out put file from
        param data. The directory will be the basename of the
        control_file.
    """

    control_file.write(params['method'] + ' input\n')

    # write file parameters
    control_file.write('!' + '*' * 20 + '\n')
    control_file.write(params['seq_file'] + '\t! seq_file\n')
    if write_input:
        with open(params['seq_file'], 'w') as fid:
            write_sequence(params['seq'], fid)

    control_file.write(str(params['n_spectra']) + '\t! seq_file\n')

    for spectrum_file in params['spectra_files']:
        control_file.write(spectrum_file + '\t! spectrum_file\n')

    if write_input:
        for n in params['n_spectra']:
            with open(params['spectra_files'][n], 'w') as fid:
                write_spectrum(
                    params['spectra_peak_groups'][n], fid)

    control_file.write(
        params['connection_file'] + '\t! connection file\n')
    if write_input:
        with open(params['connection_file'], 'w') as fid:
            write_connection(params['connections'], fid)

    if not params['initial_file']:
        control_file.write('NULL' + '\t! initial file\n')
    else:
        control_file.write(
            params['initial_file'] + '\t! initial state file\n')
    # TODO: Output initial input state file.

    control_file.write(params['output_file'] + '\t! output file\n')
    control_file.write(
        params['outtab_folder'] + '\t! outtab folder\n')

    # Write the rest of the algorithm parameters
    control_file.write('!' + '*' * 20 + '\n')
    if params['method'].lower() == 'nsga2_assign':
        control_file.write(
            str(params['group_size']) + '\t! group size\n')
        control_file.write(
            str(params['pool_size']) + '\t! pool size\n')
        control_file.write(
            str(params['steps']) + '\t! steps\n')
        control_file.write(
            str(params['nsga_attempts']) + '\t! steps\n')
        control_file.write(
            str(params['free_steps']) + '\t! free steps\n')
        control_file.write(
            str(params['mutation_rate']) + '\t! mutation_rate\n')
        control_file.write(
            str(params['additional_rate']) + '\t! additional rate\n')
        control_file.write(
            str(params['crossover_rate']) + '\t! crossover rate\n')
        control_file.write(
            str(params['null_prob']) + '\t! null prob\n')

        if params['random_seed']:
            control_file.write(
                str(params['random_seed']) + '\t! random seed\n')
        control_file.write('!' + '*' * 20 + '\n')

    elif params['method'].lower() == 'nsga2_mc_assign':
        control_file.write(
            str(params['group_size']) + '\t! group size\n')
        control_file.write(
            str(params['pool_size']) + '\t! pool size\n')
        control_file.write(
            str(params['steps']) + '\t! steps\n')
        control_file.write(
            str(params['nsga_attempts']) + '\t! steps\n')
        control_file.write(
            str(params['mc_attempts']) + '\t! steps\n')
        control_file.write(
            str(params['free_steps']) + '\t! free steps\n')
        control_file.write(
            str(params['mutation_rate']) + '\t! mutation_rate\n')
        control_file.write(
            str(params['additional_rate']) + '\t! additional rate\n')
        control_file.write(
            str(params['crossover_rate']) + '\t! crossover rate\n')
        control_file.write(
            str(params['null_prob']) + '\t! null prob\n')
        control_file.write(
            str(params['w1_min']) + '\t! w1_min\n')
        control_file.write(
            str(params['w2_min']) + '\t! w2_min\n')
        control_file.write(
            str(params['w3_min']) + '\t! w3_min\n')
        control_file.write(
            str(params['w4_min']) + '\t! w4_min\n')
        control_file.write(
            str(params['w1_max']) + '\t! w1_max\n')
        control_file.write(
            str(params['w2_max']) + '\t! w2_max\n')
        control_file.write(
            str(params['w3_max']) + '\t! w3_max\n')
        control_file.write(
            str(params['w4_max']) + '\t! w4_max\n')

        if params['random_seed']:
            control_file.write(
                str(params['random_seed']) + '\t! random seed\n')
        control_file.write('!' + '*' * 20 + '\n')


def read_outdata(params, outdata_file):
    """
    Read the output data
    :param params: algorithm parameters from control file
    :param outdata_file: open file handel
    :return spectra assignment:  list of length params['n_spectra']
        each of length params['group_size'] each of length
        params['seq'].
    """
    next(outdata_file)

    n = params['n_spectra']
    spectra_assignments = [[] for _ in range(n)]

    for lines in grouper(n, outdata_file):

        for k, line in enumerate(lines):
            values = [int(x)-1 for x in line.split()]
            values = [x if x >= 0 else None for x in values]
            spectra_assignments[k].append(values)

    return spectra_assignments


def write_outdata(spectra_assignments, outdata_file, deliminator=' '):
    """
    Write the output data
    :param spectra_assignments:  list of length params['n_spectra']
        each of length params['group_size'] each of length
        params['seq'].
    :param outdata_file: open file handel
    :param deliminator: deliminator use ' ', to avoid F90 errors.
    """
    outdata_file.write('Output_data\n')

    for assigns in zip(*spectra_assignments):
        for assign in assigns:
            line = [0 if x is None else x+1 for x in assign]
            line = deliminator.join(map(str, line)) + '\n'
            outdata_file.write(line)


def read_outtab_score(outtab_file):
    """
    Read the score from a single output table file.

    :param outtab_file: file name with path.
    """
    header = next(outtab_file)
    values = map(int, header.split()[-5:])
    score = nmrassign.base.Score(*values)

    return score


def read_outtab_scores(params):
    """
    Read the scores from all the output table files.

    :param params: dict of control params
    """
    try:
        outab_file_stub = params['outtab_folder']
    except KeyError:
        msg = 'outtab_folder must be a key in params'
        raise KeyError(msg)

    outtab_stub = os.path.basename(outab_file_stub)
    outtab_folder = os.path.abspath(os.path.dirname(outab_file_stub))

    match = os.path.join(outtab_folder, outtab_stub+'*')
    outtab_files = glob.glob(match)

    scores = []
    for outtab_file in outtab_files:
        with open(outtab_file, 'r') as outtab_fid:
            score = read_outtab_score(outtab_fid)

        scores.append(score)

    return scores
