#!/usr/bin/env python3

#claude: NOT EVALUATED BY A HUMAN CODER - written by claude, review before relying on it
#claude: smoke test that every module in the package compiles (catches Python 3 syntax regressions, e.g. print statements,
#claude: in modules no other test imports). Compiles in memory, so no .pyc files are written. The Flask websites are
#claude: not part of the package and are skipped.

import io
import os
import unittest
import warnings

import protein_complex_maps

SKIP_DIRS = ('plant_map_website', 'complex_map_website')

#claude: NOT EVALUATED BY A HUMAN CODER - written by claude, review before relying on it
def package_python_files():
    package_dir = os.path.dirname(os.path.abspath(protein_complex_maps.__file__))
    for dirpath, dirnames, filenames in os.walk(package_dir):
        dirnames[:] = sorted(d for d in dirnames if d not in SKIP_DIRS and d != '__pycache__')
        for filename in sorted(filenames):
            if filename.endswith('.py'):
                yield os.path.join(dirpath, filename)


#claude: NOT EVALUATED BY A HUMAN CODER - written by claude, review before relying on it
class CompileTest(unittest.TestCase):

    def testAllModulesCompile(self,):
        filenames = list(package_python_files())
        #claude: guard against the walk silently finding nothing (e.g. wrong package directory)
        self.assertGreater(len(filenames), 50)

        failures = []
        for filename in filenames:
            with io.open(filename, 'rb') as f:
                source = f.read()
            #claude: warnings as errors, so e.g. invalid escape sequences (SyntaxWarning on Python 3.12+) also fail
            with warnings.catch_warnings():
                warnings.simplefilter('error')
                try:
                    compile(source, filename, 'exec')
                except Exception as e:
                    failures.append('%s: %s: %s' % (filename, type(e).__name__, e))

        self.assertEqual(failures, [], '\n'.join(failures))


if __name__ == "__main__":
    unittest.main()
