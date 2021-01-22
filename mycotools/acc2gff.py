#! /usr/bin/env python3

import pandas as pd, multiprocessing as mp
import os, sys, argparse, re
from mycotools.lib.fastatools import gff2dict, dict2gff, grabGffAcc
from mycotools.lib.dbtools import db2df, masterDB
from mycotools.lib.kontools import formatPath


def main( df, column = None, db = None, gff = None, cpu = 1 ):

    acc_comp = re.compile( r'(.*?)\[(\d+)\-(\d+)\]$' )

    gff_strs = {}
    if gff: 
        if type(df) is not str:
            accessions = list(df[column])
        else:
            accession = [df]
        for acc in accessions:
            gff_strs[acc] = grabGffAcc( gff2dict(gff), acc )
    else:
        db = db.set_index( 'internal_ome' )
        if type(df) is not str:
            omes_prep = re.findall( r'^.*?_', '\n'.join(list(df[ column ])), re.M )
            omes = set( x[:-1] for x in omes_prep )
        else:
            omes_prep = re.search( r'(^.*?)_', df )
            if omes_prep is not None:
                omes = [ omes_prep[1] ]
        for ome in omes:
            gff_str = ''
            if type(df) is not str:
                accessions = [ x for x in df[ column ] if x.startswith( ome + '_' ) ]
            else:
                accessions = [ df ]
#                acc_search = acc_comp.search( accessions_prep[0] )
 #               if acc_search:
  #                  accessions = {}
   #                 for accession in accessions_prep:
    #                    accessions[ acc_search[1] ] = [ 
     #                       int(acc_search[2]),
      #                      int(acc_search[3]),
       #                     acc_search[1]
        #                    ]
         #       else:
          #          accessions = accessions_prep
            if not pd.isnull( db['gff3'][ome] ):
                gff = os.environ['MYCOGFF3'] + '/' + db['gff3'][ome]
                if not os.path.isfile( gff ):
                    print( '\t' + ome + ' invalid gff3' )
                    continue
                gff_list = gff2dict( gff )
                if gff_list[0]['source'] == 'funannotate':
                    tag = 'ID'
                else:
                    tag = 'Alias'
                mp_cmds = []
                for acc in accessions:
                    mp_cmds.append( [gff_list, acc] )
                with mp.get_context('spawn').Pool( processes = cpu ) as pool:
                    mp_strs = pool.starmap( grabGffAcc, mp_cmds )
#                    gff_str += dict2gff( grabGffAcc( gff_list, acc, tag = tag ) ) + '\n'
                for mp_str in mp_strs:
                    gff_str += dict2gff( mp_str ) + '\n'
                gff_strs[ome] =  gff_str.rstrip()
            else:
                print( '\t' + ome + ' no gff3' )

    return gff_strs


if __name__ == '__main__':

    mp.set_start_method('spawn')
    parser = argparse.ArgumentParser( 
        description = 'Inputs gff or database and new line delimitted file of accessions.'
        )
    parser.add_argument( '-i', '--input', help = 'Input file with accessions' )
    parser.add_argument( '-a', '--accession', help = 'Input accession' )
    parser.add_argument( '-g', '--gff', help = '`gff3` reference' )
    parser.add_argument( '-c', '--column', default = 0, help = 'Column to abstract from. Must be specified if there ' \
        + 'are headers. Otherwise, the default is the first tab separated column.' )
    parser.add_argument( '-o', '--ome', action = 'store_true', help = 'Output files by ome code' )
    parser.add_argument( '-d', '--database', default = masterDB(), help = 'mycodb DEFAULT: master' )
    parser.add_argument( '--cpu', type = int, default = mp.cpu_count(), \
        help = 'CPUs to parallelize accession searches. DEFAULT: all' )
    args = parser.parse_args()

    if args.input:
        input_file = formatPath( args.input )
        if args.column == 0:
            df = pd.read_csv(input_file, sep = '\t', header = None)
        else:
            df = pd.read_csv(input_file, sep = '\t')
    elif not args.accession:
        print('\nERROR: need input file or accession')
        sys.exit( 1 )

    if args.cpu < mp.cpu_count():
        cpu = args.cpu
    else:
        cpu = mp.cpu_count()

    db_path = formatPath( args.database )
    if not args.gff:
        db = db2df( formatPath(args.database) )
        if args.accession:
            gff_strs = main( args.accession, db = db, cpu = 1 )
        else:
            gff_strs = main( df, args.column, db = db, cpu = args.cpu )
    else:
        gff_path = formatPath( args.gff )
        if args.accession:
            gff_strs = main( args.accession, gff = gff_path, cpu = 1 )
        else:
            gff_strs = main( df, args.column, gff = gff_path, cpu = args.cpu )

    if args.accession:
        print( gff_strs[list(gff_strs.keys())[0]].rstrip() )
    elif args.ome:
        if not os.path.isdir( 'header2gff' ):
            os.mkdir( 'header2gff' )
        for ome in gff_strs:
            with open( 'header2gff/' + ome + '.gff3', 'w' ) as out:
                out.write( gff_strs[ome] )
    else:
        out_str = ''
        for ome in gff_strs:
            out_str += gff_strs[ome]
            with open( input_file + '.gff3', 'w' ) as out:
                out.write( out_str )

    sys.exit( 0 )