import contextlib
import os
import subprocess
from tempfile import TemporaryDirectory


class ProcessRunner(contextlib.AbstractContextManager):
    def __init__(self, env=None):
        self._env = env
        self.out = None
        self.err = None
        self.succeed = False
        self.tmpdir = None
        self.original_dir = os.getcwd()
        self._stack = contextlib.ExitStack()

    def __enter__(self):
        self.tmpdir = self._stack.enter_context(
            TemporaryDirectory(prefix='cspr_'))
        os.chdir(self.tmpdir)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        try:
            os.chdir(self.original_dir)
        finally:
            self._stack.close()

    def run(self, cmd, timeout=None):
        r = subprocess.run(cmd, capture_output=True,
                           timeout=timeout, env=self._env, text=True)
        self.succeed = r.returncode == 0
        self.out = r.stdout
        self.err = r.stderr
        return self.succeed
