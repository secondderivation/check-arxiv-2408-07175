"""Negative controls: the replay comparison must reject the failures it exists to catch."""
import unittest
from pathlib import Path
from compare_results import text_compare, posterior_compare, grids_compare, compare_run

HERE = Path(__file__).parent


class ValidationTests(unittest.TestCase):
    def test_changed_verdict_rejected(self):
        self.assertTrue(text_compare('level 0.4293 (0.57 sigma) ', 'level 0.4293 (0.57 sigma) outside prior'))

    def test_changed_value_rejected(self):
        self.assertTrue(text_compare('differential -0.0360', 'differential -0.0361'))

    def test_printed_decimal_boundary(self):
        for old, new in [('-7.995', '-7.996'), ('-5.271', '-5.272'),
                         ('-0.406', '-0.405'), ('0', '1e-3')]:
            with self.subTest(old=old, new=new):
                self.assertFalse(text_compare(old, new, atol=1e-3, rtol=0))
                self.assertFalse(text_compare(new, old, atol=1e-3, rtol=0))
        self.assertTrue(text_compare('-7.995', '-7.997', atol=1e-3, rtol=0))
        self.assertTrue(text_compare('0', '0.001000000000000001', atol=1e-3, rtol=0))

    def test_relative_decimal_boundary(self):
        self.assertFalse(text_compare('999', '1000', atol=0, rtol=1e-3))
        self.assertTrue(text_compare('998.9999', '1000', atol=0, rtol=1e-3))

    def test_nonfinite_numeric_values_rejected(self):
        for value in ['nan', 'inf', '-inf', '1e999']:
            with self.subTest(value=value):
                self.assertTrue(text_compare('1', value))

    def test_fit_boundary_through_compare_run(self):
        run = {'out': 'out/fits.out'}
        old = (HERE / 'checks' / run['out']).read_text()
        self.assertIn('-7.995)', old)
        self.assertFalse(compare_run(run, HERE / '_no_such_work_copy',
                                     old.replace('-7.995)', '-7.996)', 1)))
        self.assertTrue(compare_run(run, HERE / '_no_such_work_copy',
                                    old.replace('-7.995)', '-7.997)', 1)))

    def test_posterior_tolerance(self):
        old = (HERE / 'checks/out/posteriors.out').read_text()
        self.assertFalse(posterior_compare(old, old))
        self.assertFalse(posterior_compare(old, old.replace('level 0.4293', 'level 0.4294', 1)))
        self.assertTrue(posterior_compare(old, old.replace('level 0.4293', 'level 0.4393', 1)))
        self.assertFalse(posterior_compare(old, old.replace('(8.29 sigma)', '(8.13 sigma)', 1)))
        self.assertTrue(posterior_compare(old, old.replace('(0.57 sigma)', '(0.67 sigma)', 1)))

    def test_grids_ignore_progress_and_time(self):
        old = (HERE / 'checks/out/grids.out').read_text()
        new = '120/24311 pairs, 0.3 min\n' + old + 'run time 71.4 min\n'
        self.assertFalse(grids_compare(old, new))
        self.assertTrue(grids_compare(old, new.replace('1.82e-05', '3.00e-03', 1)))
        self.assertTrue(grids_compare(old, new.replace('14 arrays', '13 arrays')))
        self.assertTrue(grids_compare(old, ''))

    def test_reference_is_the_published_output(self):
        # End to end through compare_run: the reference must come from the package, not from a work copy
        # the check may have rewritten, so a small numeric change in stdout is still caught.
        run = {'out': 'out/crossmatch.out'}
        old = (HERE / 'checks' / run['out']).read_text()
        bogus = HERE / '_no_such_work_copy'
        self.assertFalse(compare_run(run, bogus, old))
        self.assertTrue(compare_run(run, bogus, old.replace('-0.0360', '-0.0361', 1)))


if __name__ == '__main__':
    unittest.main()
