#!/usr/bin/env python

import os
import sys
import httplib
import urllib
import json
import shutil
from optparse import OptionParser


def server_insert(args, options):
    params = {}
    params['reactant'] = ''.join(open(args[1], 'r').readlines())
    params['saddle']   = ''.join(open(args[2], 'r').readlines())
    params['product']  = ''.join(open(args[3], 'r').readlines())
    params['mode'] = ''
    if options.mode:
        params['mode'] = ''.join(open(options.mode, 'r').readlines())
    params  = urllib.urlencode(params)
    headers = {'Content-type': 'application/x-www-form-urlencoded', 'Accept': 'text/plain'}
    conn = httplib.HTTPConnection(host=options.host, port=options.port)
    conn.request('POST', '/insert', params, headers)
    response = conn.getresponse()
    print response.status, response.reason
    data = response.read()
    print data    


def server_query(args, options):
    params = {}
    params['reactant'] = ''.join(open(args[1], 'r').readlines())
    params  = urllib.urlencode(params)
    headers = {'Content-type': 'application/x-www-form-urlencoded', 'Accept': 'text/plain'}
    conn = httplib.HTTPConnection(host=options.host, port=options.port)
    conn.request('POST', '/query', params, headers)
    response = conn.getresponse()
    print response.status, response.reason
    data = json.loads(response.read())
    print data['stdout']
    if os.path.isdir('kdbmatches'):
        shutil.rmtree('kdbmatches')
    os.mkdir('kdbmatches')
    for filename in data['files']:
        f = open(os.path.join('kdbmatches', filename), 'w')
        f.write(data['files'][filename])
        f.close()


if __name__ == "__main__":

    # Parse command line options.
    usage = """%prog insert [options] reactant saddle product
        - or - 
       %prog query [options] reactant"""
    parser = OptionParser(usage = usage)
    parser.add_option("-o", "--mode", dest = "mode", 
                      help = "optional mode file",
                      default = None)
    parser.add_option("--host", dest = "host", help = "the hostname of a kdbserver",
                      default = "theory.cm.utexas.edu")
    parser.add_option("--port", dest = "port", action="store", type="int", 
                      help = "the port of a kdbserver", default = 8080)
    options, args = parser.parse_args()

    # Are we inserting or querying?
    if args[0] == 'insert':
        # Make sure we get the reactant, saddle, and product files.
        if len(args) < 4:
            parser.print_help()
            sys.exit()
        # Perform the insert.
        server_insert(args, options)

    elif args[0] == 'query':
        # Make sure we get the reactant file name.
        if len(args) < 2:
            parser.print_help()
            sys.exit()
        # Perform the query.
        server_query(args, options)

    else:
        parser.print_help()
        sys.exit()
