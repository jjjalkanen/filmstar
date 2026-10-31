"""The plotting command exits and its timeout does not leave a child behind."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time
import unittest


ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(shutil.which('timeout'), 'GNU timeout is required')
class ShutdownTests(unittest.TestCase):
    def assert_group_gone(self, process):
        try:
            os.killpg(process.pid, 0)
        except ProcessLookupError:
            return
        os.killpg(process.pid, 9)
        self.fail('Plot command left a child process running')

    def run_plot(self, limit):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)/'plot.png'
            env = dict(os.environ, SPLINE_RUN_TIMEOUT=limit)
            process = subprocess.Popen(
                [str(ROOT/'run.sh'), '--mesh', str(ROOT/'examples/meshes/tetrahedron.json'),
                 '--deg', '1', '--smooth', '0', '--save', str(output)],
                cwd=ROOT, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                text=True, start_new_session=True)
            try:
                stdout, stderr = process.communicate(timeout=30)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, 9)
                process.communicate()
                self.fail('Plot command did not exit within 30 seconds')
            self.assert_group_gone(process)
            return process.returncode, output.is_file(), stdout, stderr

    def test_saved_plot_exits(self):
        code, exists, stdout, stderr = self.run_plot('30s')
        self.assertEqual(code, 0, stderr)
        self.assertTrue(exists)
        self.assertIn('Saved', stdout)

    def test_timeout_exits_without_child(self):
        code, exists, _stdout, stderr = self.run_plot('0.1s')
        self.assertEqual(code, 124, stderr)
        self.assertFalse(exists)

    def test_spatial_timeout_exits_without_child(self):
        with tempfile.TemporaryDirectory() as directory:
            env = dict(os.environ, SURFACE_RUN_TIMEOUT='0.1s')
            process = subprocess.Popen(
                [str(ROOT/'run_surface.sh'), '--mesh',
                 str(ROOT/'examples/surface_meshes/dodecahedron.json'), '--deg', '6',
                 '--smooth', '1', '--save', str(Path(directory)/'plot.png')],
                cwd=ROOT, env=env, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, text=True, start_new_session=True)
            try:
                _stdout, stderr = process.communicate(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, 9)
                process.communicate()
                self.fail('Spatial plot ignored its timeout')
            self.assertEqual(process.returncode, 124, stderr)
            self.assert_group_gone(process)

    def test_external_termination_exits_without_child(self):
        with tempfile.TemporaryDirectory() as directory:
            process = subprocess.Popen(
                [str(ROOT/'run.sh'), '--mesh', str(ROOT/'examples/meshes/dodecahedron.json'),
                 '--deg', '1', '--smooth', '0', '--save',
                 str(Path(directory)/'plot.png')],
                cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                start_new_session=True)
            time.sleep(0.5)
            self.assertIsNone(process.poll(), 'Plot finished before termination test')
            process.terminate()
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, 9)
                self.fail('Plot command ignored termination')
            self.assert_group_gone(process)


if __name__ == '__main__':
    unittest.main()
