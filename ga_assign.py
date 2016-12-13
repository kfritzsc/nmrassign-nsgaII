"""
Python script to mac hypothetical NMR data and run a genetic
algorithm that assigns NMR data.
"""

from nmrassign.fileio import read_sequence, read_spectrum, read_connection

if __name__ == "__main__":
    import os

    # File locations
    file_dir = os.path.abspath(__file__)
    base_dir = os.path.dirname(file_dir)
    data_dir = os.path.join(base_dir, 'test', 'data', 'gb1_1')

    file_seq = os.path.join(data_dir, 'seq.txt')
    file_connection = os.path.join(data_dir, 'connection.txt')
    file_spectra = [os.path.join(data_dir, 'NCACX.txt'),
                    os.path.join(data_dir, 'NCOCX.txt')]

    # Read the files.
    seq = read_sequence(file_seq)
    connection = read_connection(file_connection)

    spectra = []
    for file_spectrum in file_spectra:
        spectra.append(read_spectrum(file_spectrum))
