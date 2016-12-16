"""
Unit tests for `nmrassign.analysis.`
"""
import os
import unittest
import shutil
import tempfile

import nmrassign.analysis as analysis
import nmrassign.fileio as fileio

class AnalysisTest(unittest.TestCase):
    """
    Test `nmrassign.analysis` functions.

    Use nearly ideal GB1 data to test reading and writing files for
    GA assign.
    """
    def setUp(self):

        # Put any temporary writing files into this temporary dir.
        self.test_dir = tempfile.mkdtemp()

        # The test data files are here.
        test_data_dir = os.path.dirname(__file__)

        # gb1_1 data
        control_file_path = os.path.join('data', 'gb1_1',
                                         'control_nsga.txt')
        self.control_file = os.path.join(test_data_dir,
                                         control_file_path)

        with open(self.control_file) as control_fid:
            self.params = fileio.read_control(control_fid,
                                              parse_input=True)
        self.scores = fileio.read_outtab_scores(self.params)

        output_file = os.path.join('data', 'gb1_1', 'outdata.txt')
        with open(output_file, 'r') as output_fid:
            self.assigns = fileio.read_outdata(self.params,
                                               output_fid)

        # gb1_2 data
        control_file_path2 = os.path.join('data', 'gb1_2',
                                         'control_nsga.txt')
        self.control_file2 = os.path.join(test_data_dir,
                                         control_file_path2)

        with open(self.control_file2) as control_fid2:
            self.params2 = fileio.read_control(control_fid2,
                                              parse_input=True)
        self.scores2 = fileio.read_outtab_scores(self.params2)

        output_file2 = os.path.join('data', 'gb1_2', 'outdata.txt')
        with open(output_file2, 'r') as output_fid2:
            self.assigns2 = fileio.read_outdata(self.params2,
                                               output_fid2)

    def tearDown(self):
        # Remove the directory after the test.
        shutil.rmtree(self.test_dir)

    def test_keep_pareto_order(self):
        groups = analysis.keep_pareto_order(self.scores)
        self.assertEqual(sum([1 for x in groups if x]), 1)

    def test_keep_pareto_order2(self):
        top_paretos = analysis.keep_pareto_order(self.scores2)
        self.assertEqual(sum([1 for x in top_paretos if x]), 3)

    def test_pareto_filter(self):
        top_paretos = analysis.keep_pareto_order(self.scores2)
        new_assigns = analysis.pareto_filter(self.assigns2,
                                    top_paretos=top_paretos)

        n_top = sum([1 for x in top_paretos if x])
        self.assertEqual(len(new_assigns[0]), n_top)


if __name__ == '__main__':
    unittest.main()
