#!/usr/bin/env python

import os
import sys
import httplib
import urllib
import json
from optparse import OptionParser


def server_insert(args, options):
    params = {}
    params['reactant'] = ''.join(open(args[1], 'r').readlines())
    params['saddle']   = ''.join(open(args[2], 'r').readlines())
    params['product']  = ''.join(open(args[3], 'r').readlines())
    params['mode']     = ''.join(open(args[4], 'r').readlines())
    params['nf']  = options.nf
    params['dc']  = options.dc
    params['mac'] = options.mac
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
    params['nf']  = options.nf
    params['dc']  = options.dc
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
    usage = """%prog insert [options] reactant saddle product mode
        - or - 
       %prog query [options] reactant"""
    parser = OptionParser(usage = usage)
    parser.add_option("-n", "--nf", dest = "nf", action="store", type="float", 
                      help = "neighbor fudge parameter",
                      default = 0.2)
    parser.add_option("-c", "--dc", dest = "dc", action="store", type="float", 
                      help = "distance cutoff parameter",
                      default = 0.3)
    parser.add_option("-m", "--mac", dest = "mac", action="store", type="float", 
                      help = "mobile atom cutoff parameter",
                      default = 0.7)
    parser.add_option("--host", dest = "host", help = "the hostname of a kdbserver",
                      default = "localhost")
    parser.add_option("--port", dest = "port", action="store", type="int", 
                      help = "the port of a kdbserver", default = 8080)
    options, args = parser.parse_args()

    # Are we inserting or querying?
    if args[0] == 'insert':
        # Make sure we get the reactant, saddle, product, and mode files.
        if len(args) < 5:
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
        sys.exit();

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        


