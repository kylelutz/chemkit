#!/usr/bin/python

# this script runs all the chemkit auto tests

import os
import sys
import subprocess

TESTS_PASSED = 0
TESTS_RUN = 0

def runtest(file):
    filePath = os.path.join(os.getcwd(), file)
    os.putenv("LD_LIBRARY_PATH", os.path.join(os.getcwd(), "../lib/"))
    os.putenv("CHEMKIT_PLUGIN_PATH", os.path.join(os.getcwd(), "../share/chemkit/plugins/"))

    print "  %-32s" % (os.path.sep.join(file.split("/")[-3:-1])),

    sys.stdout.flush()

    test = subprocess.Popen([filePath, "-silent"], stdout=subprocess.PIPE, cwd=os.path.dirname(filePath))

    ret = test.wait()
   
    # e.g. "Totals: 12 passed, 0 failed, 0 skipped"
    output = test.stdout.read()
    lastline = output.split("\n")[-2];

    passed = int(lastline.split(' ')[1])
    failed = int(lastline.split(' ')[3])
    skipped = int(lastline.split(' ')[5])
    
    global TESTS_PASSED
    global TESTS_RUN
    if not failed:
        TESTS_PASSED += 1
    TESTS_RUN += 1

    print " -- %s" % (("pass","FAIL")[failed>0])

def runtests(directory):
    directory = os.path.join("../tests/auto/", directory)

    for dir in os.listdir(directory):
        dirpath = os.path.join(directory, dir)

        if not os.path.isdir(dirpath):
            continue

        testfile = os.path.join(dirpath, dir)

        if not os.path.exists(testfile):
            continue
        
        runtest(testfile)

if __name__ == '__main__':
    print "Running Tests:"
    runtests("chemkit")
    runtests("graphics") 
    runtests("widgets")
    runtests("plugins")
    runtests("apps")
    #runtests("../benchmarks")

    print "Results:"
    print "  %i passed, %i failed, %i total" % (TESTS_PASSED, (TESTS_RUN-TESTS_PASSED), TESTS_RUN)
