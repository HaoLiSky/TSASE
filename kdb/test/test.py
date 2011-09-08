#!/usr/bin/env python

import os
import sys
import commands

os.system("rm -rf kdb kdbmatches")

out1 = commands.getoutput("../kdbinsert.py reactant.con saddle.con product.con mode.dat")

if out1 != "good":
	print "fail"
	sys.exit()

out2 = commands.getoutput("../kdbinsert.py reactant.con saddle.con product.con mode.dat")

if out2 != "duplicate of ./kdb/Al/0":
	print "fail"
	sys.exit()

commands.getoutput("../kdbquery.py reactant.con")

out3 = commands.getoutput("diff SADDLE_0 kdbmatches/SADDLE_0")

if out3 != "":
	print "fail"
	sys.exit()

commands.getoutput("../kdbquery.py product.con")

out3 = commands.getoutput("diff SADDLE_1 kdbmatches/SADDLE_1")

if out3 != "":
	print "fail"
	sys.exit()

print "success"

os.system("rm -rf kdb kdbmatches")
