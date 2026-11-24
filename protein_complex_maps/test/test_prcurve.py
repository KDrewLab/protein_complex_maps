#!/usr/bin/env python3

#claude: NOT EVALUATED BY A HUMAN CODER - written by claude, review before relying on it
#claude: tests for evaluation/plots/prcurve.py: (1) the local _binary_clf_curve copy, via precision_recall_curve_weighted
#claude: with weight_pos=1.0, gives the same curve as sklearn's precision_recall_curve; (2) smoke test that main() runs
#claude: end to end on tiny input files and writes a plot

import os
import shutil
import sys
import tempfile
import unittest

import numpy as np
from sklearn.metrics import precision_recall_curve

import protein_complex_maps.evaluation.plots.prcurve as pr


#claude: NOT EVALUATED BY A HUMAN CODER - written by claude, review before relying on it
class _QuietStdout(object):
    #claude: prcurve prints a lot of intermediate values; hide them (works on Python 2 and 3)
    def __enter__(self):
        self.saved = sys.stdout
        sys.stdout = open(os.devnull, 'w')
    def __exit__(self, *exc):
        sys.stdout.close()
        sys.stdout = self.saved


#claude: NOT EVALUATED BY A HUMAN CODER - written by claude, review before relying on it
class PRCurveTest(unittest.TestCase):

    def setUp(self,):
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self,):
        pr.plt.close('all')
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def write(self, filename, lines):
        path = os.path.join(self.temp_dir, filename)
        with open(path, 'w') as f:
            f.write(''.join(line + '\n' for line in lines))
        return path

    def testWeightedMatchesSklearn(self,):
        #claude: labels as prcurve.main builds them (1 positive, -1 negative), with tied scores. The lowest-scored item is
        #claude: a positive: sklearn >= 1.1 no longer truncates the curve once full recall is reached, while
        #claude: precision_recall_curve_weighted still does (as sklearn 0.24 did), so with trailing negatives the curves
        #claude: differ by extra points at recall 1.0
        y_true = [1, -1, 1, 1, -1, -1, 1, -1, 1, -1, -1, 1]
        y_score = [0.9, 0.8, 0.8, 0.7, 0.6, 0.6, 0.5, 0.4, 0.3, 0.3, 0.2, 0.1]

        with _QuietStdout():
            precision, recall, thresholds = pr.precision_recall_curve_weighted(y_true, y_score, weight_pos=1.0)
        sk_precision, sk_recall, sk_thresholds = precision_recall_curve(y_true, y_score)

        np.testing.assert_allclose(precision, sk_precision)
        np.testing.assert_allclose(recall, sk_recall)
        np.testing.assert_allclose(thresholds, sk_thresholds)
        #claude: hand-checked endpoints: at threshold 0.1 all 12 predictions are called, 6 positive
        self.assertAlmostEqual(precision[0], 0.5)
        self.assertAlmostEqual(recall[0], 1.0)

    def run_main(self, extra_args):
        results = self.write('results.pairsWprob', [
            'A B 0.95', 'A C 0.9', 'B C 0.85', 'C D 0.7', 'D E 0.65', 'A D 0.6',
            'E F 0.5', 'B E 0.4', 'D F 0.3', 'A F 0.2',
        ])
        positives = self.write('positives.txt', ['A B', 'A C', 'B C', 'D E', 'E F', 'D F'])
        output_file = os.path.join(self.temp_dir, 'prcurve.png')
        argv = ['prcurve.py', '--results_wprob', results, '--input_positives', positives,
                '--output_file', output_file] + extra_args(self.write)
        saved_argv = sys.argv
        sys.argv = argv
        try:
            with _QuietStdout():
                pr.main()
        finally:
            sys.argv = saved_argv
        return output_file

    def assertPNG(self, output_file):
        self.assertTrue(os.path.exists(output_file))
        with open(output_file, 'rb') as f:
            self.assertEqual(f.read(8), b'\x89PNG\r\n\x1a\n')

    def testMainSmoke(self,):
        #claude: default path: sklearn precision_recall_curve, negatives generated from the positives' proteins,
        #claude: threshold markers plotted
        output_file = self.run_main(lambda write: ['--labels', 'toy', '--plot_thresholds'])
        self.assertPNG(output_file)

    def testMainSmokeWeighted(self,):
        #claude: weighted path (local _binary_clf_curve), negatives read from a file, unscored benchmark pairs added
        output_file = self.run_main(lambda write: [
            '--input_negatives', write('negatives.txt', ['A D', 'C D', 'B E', 'A F', 'C F']),
            '--weight_positives', '0.5', '--complete_benchmark'])
        self.assertPNG(output_file)


if __name__ == "__main__":
    unittest.main()
