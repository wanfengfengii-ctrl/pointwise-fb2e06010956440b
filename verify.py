import sys
sys.dont_write_bytecode = True
from pathlib import Path
import subprocess
import unittest


def main():
    root = Path(__file__).resolve().parent
    suite = unittest.defaultTestLoader.discover(
        str(root/"tests"), pattern="test_*.py", top_level_dir=str(root))
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():
        return 1
    # The unittest suite covers the full V1 and connectivity chains; also
    # execute the provided in-memory connectivity acceptance as delivered.
    acceptance = subprocess.run(
        [sys.executable, "-B", str(root/"connectivity_cases.py"), "--memory"],
        cwd=str(root))
    return 0 if acceptance.returncode == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
