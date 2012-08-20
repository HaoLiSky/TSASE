#!/usr/bin/env python

import os
import commands
import tempfile
import json
import glob
from optparse import OptionParser
import bottle as b
from kdb import MOBILE_ATOM_CUTOFF, NEIGHBOR_FUDGE, DISTANCE_CUTOFF


QUERY_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'kdbquery.py')
INSERT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'kdbinsert.py')
KDB_PATH = os.path.join(os.path.abspath(os.getcwd()), 'kdb')


@b.post('/insert')
def insert():
    tmpdir = tempfile.mkdtemp()
    f = open(os.path.join(tmpdir, 'reactant.con'), 'w')
    f.write(b.request.forms.get('reactant'))
    f.close()
    f = open(os.path.join(tmpdir, 'saddle.con'), 'w')
    f.write(b.request.forms.get('saddle'))
    f.close()
    f = open(os.path.join(tmpdir, 'product.con'), 'w')
    f.write(b.request.forms.get('product'))
    f.close()
    f = open(os.path.join(tmpdir, 'mode.dat'), 'w')
    f.write(b.request.forms.get('mode'))
    f.close()
    output = commands.getoutput('%s %s %s %s %s --nf=%f --dc=%f --mac=%f --kdbdir=%s' % (INSERT_PATH,
                      os.path.join(tmpdir, 'reactant.con'),
                      os.path.join(tmpdir, 'saddle.con'),
                      os.path.join(tmpdir, 'product.con'),
                      os.path.join(tmpdir, 'mode.dat'),
                      NEIGHBOR_FUDGE, DISTANCE_CUTOFF, MOBILE_ATOM_CUTOFF, KDB_PATH)) 
    return output
    
@b.post('/query')
def query():
    tmpdir = tempfile.mkdtemp()
    f = open(os.path.join(tmpdir, 'reactant.con'), 'w')
    f.write(b.request.forms.get('reactant'))
    f.close()
    output = commands.getoutput('%s %s --nf=%f --dc=%f --kdbdir=%s' % 
                                (QUERY_PATH, os.path.join(tmpdir, 'reactant.con'), NEIGHBOR_FUDGE, DISTANCE_CUTOFF, KDB_PATH)) 
    filelist = glob.glob(os.path.join('kdbmatches', '*'))
    ret = {}
    ret['stdout'] = output
    ret['files'] = {}
    for filename in filelist:
        ret['files'][os.path.basename(filename)] = ''.join(open(filename, 'r').readlines())
    return json.dumps(ret)
    

@b.route('/')
def index():
    return '<b>... under construction ...</b>'

if __name__ == "__main__":
    # Parse command line options.
    usage = """%prog [options]"""
    parser = OptionParser(usage = usage)
    parser.add_option("-n", "--nf", dest = "nf", action="store", type="float", 
                      help = "neighbor fudge parameter",
                      default = NEIGHBOR_FUDGE)
    parser.add_option("-c", "--dc", dest = "dc", action="store", type="float", 
                      help = "distance cutoff parameter",
                      default = DISTANCE_CUTOFF)
    parser.add_option("-m", "--mac", dest = "mac", action="store", type="float", 
                      help = "mobile atom cutoff parameter",
                      default = MOBILE_ATOM_CUTOFF)
    parser.add_option("--host", dest = "host", help = "the hostname to use",
                      default = "localhost")
    parser.add_option("--port", dest = "port", action="store", type="int", 
                      help = "the port to use", default = 8080)
    options, args = parser.parse_args()
    NEIGHBOR_FUDGE = options.nf
    MOBILE_ATOM_CUTOFF = options.mac
    DISTANCE_CUTOFF = options.dc
    b.run(host=options.host, port=options.port, reloader=True)