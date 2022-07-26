#! /usr/bin/env python

import os
# from ROOT import *
# import ROOT as R
# import numpy as np
# import re
import sys
import glob
# import math as m
# import datetime as dt
# from array import array
import subprocess as sp
import time


def fwk(fname, eventsstr):
    # sp.call(['cmd3fwk_xr', '-i', '"%s"' % fname, 'fwkcfg_online.opts'])
    cmd = ['cmd3fwk_xr', '-i', fname, eventsstr, 'jobopts/fwkcfg_online.opts']
    n = False
    for ia, a in enumerate( sys.argv ):
        if a == '-n':
            n = int( sys.argv[ ia + 1 ] )
    if n:
        cmd.append( '-n' )
        cmd.append( str( n ) )
    print '*'*80
    print '*'*80
    print ' '.join(cmd)
    print '*'*80
    print '*'*80
    # sp.call(['cmd3fwk_xr', '-i', fname, 'fwkcfg_online.opts', '-n 300'])
    # sp.call(['cmd3fwk_xr', '-i', fname, 'fwkcfg_online.opts'])
    sp.call( cmd )


def doit(fname, eventsstr):
    fnam = fname.split('/')[-1].split('.bz2')[0]+'.root'
    print fnam
    # print ' '.join(['root', '-b', '-q', """'doit.C("%s")'""" % fnam])
    # sp.call(['root', '-b', '-q', """'doit.C("%s")'""" % fnam])
    command = """nice -5 root -b -q 'doit.C("%s")'""" % fnam
    print '*'*80
    print '*'*80
    print command
    print '*'*80
    print '*'*80
    os.system(command)
    sp.call(['rm', '-rf', 'tmp/'+fnam])


def getFile(run):
    print run
    out = sp.check_output(['lfc-ls', '-r', r'%raw%'+str(run)+r'%'])
    out = out.split('\n')[2]
    print out
    if 'cmd' in out and 'online' in out:
        # out = out.replace( r'//cmd//', r'//sl10cmd//' )# tmp solution
        out = out.replace( r'//sl10cmd//', r'//cmd//' )# tmp solution
        out = out.replace( r'//slcmd//', r'//cmd//' )# tmp solution
        return True, out
    else:
        names = glob.glob( '/storeA/daqdata/online/incoming/run%s.mid -m "' % str(run) )
        if len(names) == 1:
            print '\33[33mOne file has been found :  {}\33[0m'.format( names[0] )
            return True, names[0]
        elif len(names) > 1:
            print '\33[31m{} files have been found :  {}\33[0m'.format( len(names), names )
            return False, out
        print 'The bad run number :  {}.'.format(run)
        return False, out

def getFiles( args ):
    for ia, a in enumerate( args ):
        if a == 'range':
            try:
                run0 = int( args[ia+1] )
                run1 = int( args[ia+2] )
                files = []
                for i in xrange( run0, run1 + 1 ):
                    good, f = getFile(i)
                    if good:
                        files.append( f )
                return files
            except:
                raise SystemError
        elif a == 'runs':
            files = []
            ja = ia + 1
            while True:
                try:
                    run = int( args[ ja ] )
                except:
                    break
                good, f = getFile( run )
                if good:
                    files.append(f)
                ja += 1
            print 'Obtained files:', files
            return files

    '''
    if runnums[0] == 'range':
        files = []
        for i in xrange(int(runnums[1]), int(runnums[2])+1):
            good, f = getFile(i)
            if good:
                files.append(f)
            print i
            out = sp.check_output(['lfc-ls', '-r', r'%raw%'+str(i)+r'%'])
            out = out.split('\n')[2]
            print out
            if 'cmd' in out and 'online' in out:
                files.append(out)
            else:
                print 'The bad run number :  {}.'.format(i)
        print 'Obtained files:', files
        return files
    elif runnums[0] == 'runs':
        files = []
        for run in runnums[1:]:
            good, f = getFile(int(run))
            if good:
                files.append(f)
        print 'Obtained files:', files
        return files
    else:
        usage()
        return
    '''



def usage():
    print 'Usage:'
    print '\t{} range runnum_initial runnum_final'.format(sys.argv[0])
    print '\t{} runs runnum1'.format(sys.argv[0])
    print '\t{} runs runnum1 runnum2 runnum3 ... runnumN'.format(sys.argv[0])




def main( runnum, eventsstr ):
    b, fil = getFile( runnum )
    #for fil in files:
    t0 = time.time()
    print '\33[1m\33[33mLocal time :', time.localtime( t0 ), '\33[0m'
    fwk( fil, eventsstr )
    print '\33[1m\33[33mFWK time :', time.time()-t0, 'seconds\33[0m'
    t1 = time.time()
    print '\33[1m\33[33mLocal time :', time.localtime( t1 ), '\33[0m'
    #doit( fil )
    print '\33[1m\33[33mFWK time :', time.time()-t1, 'seconds\33[0m'
    print '\33[1m\33[33mLocal time :', time.localtime( t0 ), '\33[0m'


#if __name__ == '__main__':
#    if len(sys.argv) > 1:
#        main(sys.argv[1:])
#    else:
#        usage()

