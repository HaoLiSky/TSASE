#!/usr/bin/env python

import os
import commands
import tempfile
import json
import glob
import bottle as b

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
    nf = float(b.request.forms.get('nf'))
    dc = float(b.request.forms.get('dc'))
    mac = float(b.request.forms.get('mac'))
    output = commands.getoutput('%s %s %s %s %s --nf=%f --dc=%f --mac=%f --kdbdir=%s' % (INSERT_PATH,
                      os.path.join(tmpdir, 'reactant.con'),
                      os.path.join(tmpdir, 'saddle.con'),
                      os.path.join(tmpdir, 'product.con'),
                      os.path.join(tmpdir, 'mode.dat'),
                      nf, dc, mac, KDB_PATH)) 
    return output
    
@b.post('/query')
def query():
    tmpdir = tempfile.mkdtemp()
    f = open(os.path.join(tmpdir, 'reactant.con'), 'w')
    f.write(b.request.forms.get('reactant'))
    f.close()
    nf = float(b.request.forms.get('nf'))
    dc = float(b.request.forms.get('dc'))
    output = commands.getoutput('%s %s --nf=%f --dc=%f --kdbdir=%s' % 
                                (QUERY_PATH, os.path.join(tmpdir, 'reactant.con'), nf, dc, KDB_PATH)) 
    filelist = glob.glob(os.path.join('kdbmatches', '*'))
    ret = {}
    ret['stdout'] = output
    ret['files'] = {}
    for filename in filelist:
        ret['files'][os.path.basename(filename)] = ''.join(open(filename, 'r').readlines())
    return json.dumps(ret)
    

@b.route('/')
def index():
    return '<b>Hello  </b>'

b.run(host='localhost', port=8080, reloader=True)
