#!/usr/bin/env python

import os
import commands
import tempfile
import json
import glob
from optparse import OptionParser
import bottle as b
from tsase.data import elements, num_elements
import kdb
from kdb import MOBILE_ATOM_CUTOFF, NEIGHBOR_FUDGE, DISTANCE_CUTOFF

QUERY_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'kdbquery.py')
INSERT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'kdbinsert.py')
KDB_PATH = os.path.join(os.path.abspath(os.getcwd()), 'kdb')
WWW_PATH = os.path.join(os.path.abspath(os.getcwd()), '..', 'www')
PATH = os.path.dirname(os.path.abspath(__file__))

b.TEMPLATE_PATH.append(os.path.join(PATH, 'templates'))

@b.post('/insert')
def insert():
    tmpdir = tempfile.mkdtemp()
    f = open(os.path.join(tmpdir, 'reactant'), 'w')
    f.write(b.request.forms.get('reactant'))
    f.close()
    f = open(os.path.join(tmpdir, 'saddle'), 'w')
    f.write(b.request.forms.get('saddle'))
    f.close()
    f = open(os.path.join(tmpdir, 'product'), 'w')
    f.write(b.request.forms.get('product'))
    f.close()
    mode = b.request.forms.get('mode')
    if mode:
      f = open(os.path.join(tmpdir, 'mode.dat'), 'w')
      f.write(mode)
      f.close()
    command = '%s %s %s %s --nf=%f --dc=%f --mac=%f --kdbdir=%s' % (INSERT_PATH,
                      os.path.join(tmpdir, 'reactant'),
                      os.path.join(tmpdir, 'saddle'),
                      os.path.join(tmpdir, 'product'),
                      NEIGHBOR_FUDGE, DISTANCE_CUTOFF, MOBILE_ATOM_CUTOFF, KDB_PATH)
    if mode:
        command += ' --mode=%s' % os.path.join(tmpdir, 'mode.dat')
    output = commands.getoutput(command) 
    return output


@b.post('/query')
def query():
    tmpdir = tempfile.mkdtemp()
    f = open(os.path.join(tmpdir, 'reactant'), 'w')
    f.write(b.request.forms.get('reactant'))
    f.close()
    output = commands.getoutput('%s %s --nf=%f --dc=%f --kdbdir=%s' % 
                                (QUERY_PATH, os.path.join(tmpdir, 'reactant'), NEIGHBOR_FUDGE, DISTANCE_CUTOFF, KDB_PATH)) 
    filelist = glob.glob(os.path.join('kdbmatches', '*'))
    ret = {}
    ret['stdout'] = output
    ret['files'] = {}
    for filename in filelist:
        ret['files'][os.path.basename(filename)] = ''.join(open(filename, 'r').readlines())
    return json.dumps(ret)
    

@b.route('/')
def index():
    return '... under construction ...'


def toSymbolList(items):
    symbols = []
    if items is None:
        return symbols
    for item in items:
        print item
        newSymbol = None
        try:
            newSymbol = elements[int(item)]['symbol']
        except:
            pass
        try:
            newSymbol = elements[item.lower().capitalize()]['symbol']
        except:
            pass
        if newSymbol is None:
            for i in range(1, num_elements):
                if item.lower() == elements[i]['name']:
                    newSymbol = elements[i]['symbol']
                    break
        if newSymbol is not None:
            symbols.append(newSymbol)
    return symbols


@b.route('/browse')
def browse():
    filter = toSymbolList(b.request.query.filter.split())
    results = [os.path.relpath(r) for r in kdb.query_has_all(KDB_PATH, filter)]
    return b.template('browse', filter=filter, results=results)

@b.route('/static/<filename>')
def static(filename):
    return b.static_file(filename, os.path.join(PATH, 'static'))    

@b.route('/www/<filename>')
def www(filename):
    return b.static_file(filename, WWW_PATH)    

@b.route('/kdb/<combo>/<number>/<filename>')
def kdbpath(combo, number, filename):
    return b.static_file(os.path.join(combo, number, filename), KDB_PATH)    


    

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
