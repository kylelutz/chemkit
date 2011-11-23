#!/usr/bin/python

import os

os.system("git log --format='%aN' | sort -u")
