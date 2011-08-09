#!/usr/bin/python

# this script counts and outputs the number of
# lines of source code for chemkit

import os

if __name__ == '__main__':
    os.system("sloccount ../src/apps ../src/chemkit ../src/io ../src/graphics ../src/plugins ../src/web ../src/widgets ../examples ../tests")
