"""
Unit tests for `nmrassign.fileio.`
"""

import os
import unittest
import shutil
import tempfile

import nmrassign.fileio as fileio


class FileIOTest(unittest.TestCase):
    """
    Test `nmrassign.fileio.` read and write functions.

    Use nearly ideal GB1 data to test reading and writing files for
    GA assign.
    """
    def setUp(self):

        # Put any temporary writing files into this temporary dir.
        self.test_dir = tempfile.mkdtemp()

        # The test data files are here.
        test_data_dir = os.path.dirname(__file__)

        # sequence files
        seq_file = os.path.join('data', 'gb1_1', 'seq.txt')
        self.gb1_seq = os.path.join(test_data_dir, seq_file)

        seq_file = os.path.join('data', 'gb1_1', 'seq.fasta')
        self.gb1_seq_fasta = os.path.join(test_data_dir, seq_file)

        seq_file = os.path.join('data', 'gb1_1', 'bad_seq.fasta')
        self.gb1_bad_seq_fasta = os.path.join(test_data_dir,
                                              seq_file)

        # other files.
        gb1_spec_file_1 = os.path.join('data', 'gb1_1', 'NCACX.txt')
        self.gb1_spec_1 = os.path.join(test_data_dir,
                                       gb1_spec_file_1)
        gb1_con_file = os.path.join('data', 'gb1_1',
                                    'connection.txt')
        self.gb1_con = os.path.join(test_data_dir, gb1_con_file)

        # control file
        control_nsga_file = os.path.join('data', 'gb1_1',
                                         'control_nsga.txt')
        self.gb1_control_nsga_file = os.path.join(test_data_dir,
                                                  control_nsga_file)

        # output files
        gb1_tab_file = os.path.join('data', 'gb1_1', 'output',
                                    'outtab001.txt')
        self.gb1_tab_file = os.path.join(test_data_dir, gb1_tab_file)

    def tearDown(self):
        # Remove the directory after the test.
        shutil.rmtree(self.test_dir)

    def test_read_sequence(self):
        """
        Read a regular correct sequence file.
        """
        with open(self.gb1_seq, 'r') as fid:
            sequence = fileio.read_sequence(fid)
        self.assertEqual(sequence[:5], list('MQYKL'))

    def test_read_sequence_fasta(self):
        """
        Read a fasta like sequence.
        """
        with open(self.gb1_seq_fasta, 'r') as fid:
            sequence = fileio.read_sequence(fid)
            self.assertEqual(sequence[:5], list('MQYKL'))

    def test_read_sequence_bad_fasta(self):
        """
        Read_sequence for a fasta file with multiple sequences.

        In this case it is not known which sequence should be used
        so a ValueError should be thrown.
        """
        with open(self.gb1_bad_seq_fasta, 'r') as fid:
            self.assertRaises(ValueError, fileio.read_sequence, fid)

    def test_sequence_round_trip(self):
        """
        Read sequence, write sequence, read sequence, compare.
        """
        sequence = fileio.read_sequence(self.gb1_seq)

        with open(os.path.join(self.test_dir, 'test_seq.fasta'),
                  'w') as fid:
            fileio.write_sequence(sequence, fid)

        with open(os.path.join(self.test_dir, 'test_seq.fasta'),
                  'r') as fid:
            read_seq = fileio.read_sequence(fid)

        self.assertEqual(sequence, read_seq)

    def test_read_groups_assignments(self):
        """
        Read a list of assignment groups from an ideal file.
        """
        with open(self.gb1_spec_1, 'r') as fid:
            groups = fileio.read_spectrum(fid)

        self.assertEqual(groups[0].assignments, 'KM')

    def test_groups_assignments_round_trip(self):
        """
        Read, write then read group assignments and compare.
        """
        with open(self.gb1_spec_1, 'r') as fid:
            groups = fileio.read_spectrum(fid)

        with open(os.path.join(self.test_dir, 'test_spec_1.txt'),
                  'w') as fid:
            fileio.write_spectrum(groups, fid)

        with open(os.path.join(self.test_dir, 'test_spec_1.txt'),
                  'r') as fid:
            new_groups = fileio.read_spectrum(fid)

        self.assertEqual(groups[10].assignments,
                         new_groups[10].assignments)

    def test_read_connection(self):
        """
        Read a list of connection form an ideal file.
        """
        with open(self.gb1_con, 'r') as fid:
            connections = fileio.read_connection(fid)
            connection = connections[0]
        # The id number are converted to zero based convention.
        self.assertEqual(connection.ids[0], 0)

    def test_connections_round_trip(self):
        """
        Read, write then read group assignments and compare.
        """
        with open(self.gb1_con, 'r') as fid:
            connections = fileio.read_connection(fid)

        with open(os.path.join(self.test_dir, 'test_conn.txt'),
                  'w') as fid:
            fileio.write_connection(connections, fid)

        with open(os.path.join(self.test_dir, 'test_conn.txt'),
                  'r') as fid:
            new_connections = fileio.read_connection(fid)

        self.assertEqual(connections, new_connections)

    def test_read_control_spot1(self):
        """
        Read control file.
        """
        with open(self.gb1_control_nsga_file, 'r') as fid:
            params = fileio.read_control(fid)
            self.assertEqual(params['steps'], 20)

    def test_read_control_spot2(self):
        """
        Read control file.
        """
        with open(self.gb1_control_nsga_file, 'r') as fid:
            params = fileio.read_control(fid)
            self.assertEqual(params['random_seed'], 1)
            self.assertEqual(params['random_seed'], 1)

    def test_control_round_trip(self):
        """
        Read control file , write control, read control compare.
        """
        with open(self.gb1_control_nsga_file, 'r') as fid:
            params = fileio.read_control(fid, parse_input=False)

        with open(os.path.join(self.test_dir, 'test_control.txt'),
                  'w') as fid:
            fileio.write_control(params, fid, write_input=False)

        with open(os.path.join(self.test_dir, 'test_control.txt'),
                  'r') as fid:
            read_params = fileio.read_control(fid, parse_input=False)

        self.assertEqual(params, read_params)

if __name__ == "__main__":
    unittest.main()
