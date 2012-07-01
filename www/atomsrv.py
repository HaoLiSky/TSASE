#!/usr/bin/env python

import os
import sys
import json
import glob
import bottle as b

import tsase
import ase.io
from ase.io import vasp

def no_escape(string):
    return string

def get_pathlist(path):
    if os.path.isfile(path):
        path = os.path.dirname(path)
    pathlist = []
    pathlist.insert(0, '<b><a href="%s" class="directory">%s</a></b>' % (os.path.dirname(path), '..'))
    ls = glob.glob(os.path.join(path, '*'))    
    for l in ls:
        if os.path.isdir(l):
            pathlist.insert(1, '<a href="%s" class="directory">%s</a>' % (l, os.path.basename(l)))
        else:
            pathlist.append('<a href="%s" class="file">%s</a>' % (l, os.path.basename(l)))
    return pathlist
    
def atoms_to_jsatoms(atoms):
    return json.dumps([{'symbol': a.symbol, 'position': list(a.position)} for a in atoms])

@b.route('/atomsrv/<fname>')
def static(fname):
    whitelist = ['xyz.js']
    if fname in whitelist:
        return ''.join(open(fname).readlines())
    b.abort(401)

def read(filename):
    try:
        data = tsase.io.read_xdatcar(filename)
        if len(data) < 1:
            raise
    except:
        try:
            f = open(filename, 'r')
            data = []
            while True:
                try:
                    data.append(vasp.read_vasp(f))
                except:
                    f.close()
                    break
            if len(data) < 1:
                raise
        except:
            try:
                data = ase.io.read(filename + "@:")
            except:
                try:
                    data = tsase.io.read_con(filename)
                    if len(data) < 1:
                        raise
                except:
                    try:
                        data = tsase.io.read_bopfox(filename)
                        if len(data) < 1:
                            raise
                    except:
                        b.abort(401, 'Could not read %s' % filename)
    if type(data) is list:
        return data[0]
    return data
    

@b.route('/<path:path>')
def index(path):
    path = '/' + path
    pathlist = get_pathlist(path)
    atoms = None
    if os.path.isfile(path):
        atoms = read(path)
        atoms = atoms_to_jsatoms(atoms)
    tpl = b.SimpleTemplate(''.join(open('atomsrv.html','r').readlines()), escape_func=no_escape)
    return tpl.render(pathlist=pathlist, atoms=atoms)
            

b.run(host=sys.argv[1], port=int(sys.argv[2]))
