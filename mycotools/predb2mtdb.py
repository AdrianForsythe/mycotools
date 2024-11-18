#!/usr/bin/env python3

# NEED source to reference the annotation source
# NEED to error check FAA generation simply by file size

import os
import re
import sys
import copy
import shutil
import multiprocessing as mp
import mmap
from tqdm import tqdm
from collections import Counter, defaultdict
from mycotools.lib.kontools import gunzip, mkOutput, format_path, eprint, vprint
from mycotools.lib.biotools import gff2list, list2gff, fa2dict, dict2fa, \
    gff3Comps, gff2Comps, gtfComps
from mycotools.lib.dbtools import mtdb, primaryDB, loginCheck
from mycotools.utils.gtf2gff3 import main as gtf2gff3
from mycotools.utils.curGFF3 import main as curGFF3
from mycotools.utils.gff2gff3 import main as gff2gff3
from mycotools.utils.curGFF3 import rename_and_organize as rename_and_organize
from mycotools.gff2seq import aamain as gff2seq
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed

predb_headers = [
    'assembly_accession', 'previous_ome', 
    'genus', 'species', 'strain', 'version', 'biosample',
    'assemblyPath', 'gffPath', 'genomeSource (ncbi/jgi/new)', 
    'useRestriction (yes/no)', 'published',
    'has_gff', 'has_faa'
]


def acq_forbid_omes(file_path):
    """Parse a file with forbidden ome accessions - ome codes that have been
    used before and are no longer valid"""
    if not os.path.isfile(file_path):
        return set()
    with open(file_path, 'r') as raw:
        relics = set([x.rstrip() for x in raw])
    return relics

def prep_output(base_dir):
    out_dir = mkOutput(base_dir, 'predb2mtdb')
    wrk_dir = out_dir + 'working/'
    dirs = [
        out_dir, wrk_dir, wrk_dir + 'gff3/', wrk_dir + 'fna/', wrk_dir + 'faa/'
        ]
    for dir_ in dirs:
        if not os.path.isdir(dir_):
            os.mkdir(dir_)
    return dirs[:2]

def copy_file(old_path, new_path):
    """Copy a file with better error handling"""
    try:
        shutil.copy(old_path, new_path)
        return True
    except (IOError, OSError) as e:
        eprint(f"\nERROR: Failed to copy file from {old_path} to {new_path}")
        eprint(f"Error details: {str(e)}")
        raise IOError(f"Failed to copy {old_path} to {new_path}: {str(e)}")

def move_biofile(old_path, ome, typ, wrk_dir, suffix = ''):
    """Move biological file with better error handling"""
    old_path = format_path(old_path)
    if not os.path.exists(old_path):
        raise IOError(f"Input file does not exist: {old_path}")
        
    if old_path.endswith('.gz'):
        if not os.path.isfile(old_path[:-3]):
            temp_path = gunzip(old_path)
            new_path = wrk_dir + ome + '.' + typ + suffix
        else:
            new_path = wrk_dir + ome + '.' + typ + suffix
            temp_path = old_path[:-3]
    else:
        new_path = wrk_dir + ome + '.' + typ + suffix
        temp_path = old_path

    # Ensure the destination directory exists
    os.makedirs(os.path.dirname(new_path), exist_ok=True)

    try:
        copy_file(temp_path, new_path)
    except IOError as e:
        eprint(f"\nERROR: Failed to copy {temp_path} to {new_path}")
        eprint(f"Error details: {str(e)}")
        raise

    return new_path


def gen_predb():
    example = [
        'Fibsp1', 'fibpsy1', 'Fibularhizoctonia', 'psychrophila', 'CBS',
        '1.0', 'n', '<PATH/TO/ASSEMBLY>', '<PATH/TO/GFF3>', 'jgi',
        'no', '2018'
        ]
    eprint('INSTRUCTIONS: fill in each column with the relevant information and \
        separate each column by a tab. The predb can be filled in \
        via spreadsheet software and exported as a tab delimited `.tsv`. \
        ASSEMBLY ACCESSIONS and PREVIOUS_OME fields must be unique to the \
        genome; otherwise predb2mtdb will update the corresponding database entry. \
        Novel data must be filled in as "new" for the genomeSource column.', flush = True)
    outputStr = '#' + '\t'.join(predb_headers)
    outputStr += '\n#' + '\t'.join(example) + '\n'

    return outputStr

def read_predb(predb_path, spacer = '\t'):
    """Modified to handle optional GFF/FAA"""
    required_headers = {
        'assembly_accession', 'genus', 'assemblyPath',
        'genomeSource (ncbi/jgi/new)'
    }  # Remove gffPath from required headers
    
    allowed_headers = {
        'previous_ome', 'assembly_acc', 'assembly_accession',
        'genus', 'species', 'strain', 'version', 'biosample',
        'assemblyPath', 'gffPath', 'genomeSource (ncbi/jgi/new)', 
        'useRestriction (yes/no)', 'published', 'restriction',
        'source', 'fna_path', 'gff3_path'
    }
    allowed2required = {'assembly_acc': 'assembly_accession',
                        'gff3_path': 'gffPath', 'fna_path': 'assemblyPath',
                        'source': 'genomeSource (ncbi/jgi/new)'}

#    predb, headers = {}, None
    predb = defaultdict(list)
    i2header = {}
    with open(predb_path, 'r') as raw:
        for i, line in enumerate(raw):
            # flexibly identify headers from the first column
            if line.startswith('#') and not i2header and i == 0:
                d = line.split('\t')
                for i0, head_p in enumerate(d):
                    head = head_p.rstrip()
                    if head in allowed_headers:
                        i2header[i0] = head
                    elif head.replace('#', '') in allowed_headers:
                        i2header[i0] = head.replace('#','')

                for head in i2header.values():
                    if head in allowed2required:
                        required_headers.remove(allowed2required[head])
#                    elif head in required_headers:
 #                       required_headers.remove(head)
                missing_headers = \
                    required_headers.difference(set(i2header.values()))
                if missing_headers:
                    eprint(f'{spacer}ERROR: Required columns missing: ' \
                         + f'{missing_headers}', flush = True)
                    sys.exit(4)
 #               if not headers:
  #                  predb = {x: [] for x in line.rstrip()[1:].split('\t')}
   #                 headers = list(predb.keys())
            elif not line.startswith('#') and line.rstrip():
                entry = line.split('\t')
                # proceed with default header organization scheme
                if not i2header:
                    if len(entry) != len(predb_headers):
                        eprint(spacer + 'ERROR: Incorrect columns, line ' + str(i),
                               flush = True)
                        eprint(predb_headers, '\n', entry, flush = True)
                        sys.exit(3)
                    for i1, v in enumerate(entry):
                        predb[predb_headers[i1]].append(v.rstrip())
                # allow for flexible header identification
                else:
                    used = []
                    for i1, v in enumerate(entry):
                        if i1 in i2header:
                            head = i2header[i1]
                            predb[head].append(v.rstrip()) 
                            used.append(i1)
                    for mi in set(i2header.keys()).difference(set(used)):
                        predb[i2header[mi]].append('')

    if not 'assembly_acc' in predb and not 'assembly_accession' in predb:
        raise KeyError('unique assembly accessions are required')
    elif 'assembly_accession' in predb:
        predb['assembly_acc'] = predb['assembly_accession']
        del predb['assembly_accession']

    if len(predb['assembly_acc']) != len(set(predb['assembly_acc'])):
        raise KeyError('unique assembly accessions are required')

    if not 'source' in predb:
        if 'genomeSource (ncbi/jgi/new)' in predb:
            predb['source'] = predb['genomeSource (ncbi/jgi/new)']
            del predb['genomeSource (ncbi/jgi/new)']
        elif 'genomeSource' in predb:
            predb['source'] = predb['genomeSource']
            del predb['genomeSource']
        else:
            raise KeyError
    if not 'assemblyPath' in predb:
        if 'fna_path' in predb:
            predb['assemblyPath'] = predb['fna_path']
            del predb['fna_path']
        else:
            raise KeyError
    if not 'gffPath' in predb:
        if 'gff3_path' in predb:
            predb['gffPath'] = predb['gff3_path']
            del predb['gff3_path']
        else:
            raise KeyError
    try:
        if 'restriction' not in predb \
            and 'useRestriction (yes/no)' in predb:
            predb['restriction'] = predb['useRestriction (yes/no)']
            del predb['useRestriction (yes/no)']
        if any(x.lower() not in {'y', 'n', 'yes', 'no', '', 'true', 'false'} \
               for x in predb['restriction']):
            eprint(spacer + 'ERROR: useRestriction entries must be in {y, n, yes, no}', flush = True)
            sys.exit(4)
    except KeyError:
        if not 'restriction' in predb and 'published' not in predb:
            raise KeyError('restriction/published columns required')
        elif 'published' in predb and 'restriction' not in predb:
            predb['restriction'] = [bool(x) for x in predb['published']]

    if any(x.lower() not in {'jgi', 'ncbi', 'new'} \
        for x in predb['source']):
        eprint([predb['assembly_acc'][i] for i, v in enumerate(predb['source']) \
                if v not in {'jgi', 'ncbi', 'new'}])
        eprint(spacer + 'ERROR: genomeSource entries must be in {jgi, ncbi, new}', flush = True)
        sys.exit(5)

    missing_from_predb = list(set(predb_headers).difference(set(predb.keys())))
    for i, path in enumerate(predb['assemblyPath']):
        predb['assemblyPath'][i] = format_path(predb['assemblyPath'][i])
        predb['gffPath'][i] = format_path(predb['gffPath'][i])
        if predb['restriction'][i].lower() in {'y', 'yes', 'true'}:
            predb['restriction'][i] = True
        elif predb['restriction'][i].lower() in {'n', 'no', 'false'}:
            predb['restriction'][i] = ''
        if not predb['species'][i]:
            predb['species'][i] = 'sp.'
        # add a blank entry for each missing column from the entire predb
        for missing_header in missing_from_predb:
            predb[missing_header].append('')

    # Add validation for optional files
    if not 'has_gff' in predb:
        predb['has_gff'] = ['yes' if os.path.exists(p) else 'no' 
                           for p in predb.get('gffPath', [])]
    if not 'has_faa' in predb:
        predb['has_faa'] = ['no'] * len(predb['assembly_acc'])  # Default to no

    # Make gffPath optional
    if not 'gffPath' in predb:
        predb['gffPath'] = [''] * len(predb['assembly_acc'])

    return dict(predb)

def sub_disallowed(data, disallowed = r"""[^\w\d]"""):
    if data:
        return re.sub(disallowed, '', data)
    else:
        return data

def predb2mtdb(predb):
    infdb = mtdb()
    if not 'assemblyPath' in predb and 'fna' in predb:
        predb['assemblyPath'] = predb['fna']
    if not 'gffPath' in predb and 'gff3' in predb:
        predb['gffPath'] = predb['gff3']
    if not 'previous_ome' in predb and 'ome' in predb:
        predb['previous_ome'] = predb['ome']
    elif not 'previous_ome' in predb:
        predb['previous_ome'] = \
            ['' for x in predb[list(predb.keys())[0]]]
    for i, code in enumerate(predb['genus']):
        if not code:
            raise ValueError(f'no genus detected for genus {i}')
        toAdd = {
            'assembly_acc': predb['assembly_acc'][i],
            'ome': predb['previous_ome'][i].lower(),
            'genus': sub_disallowed(predb['genus'][i]),
            'species': sub_disallowed(predb['species'][i]),
            'strain': re.sub(r'[^a-zA-Z0-9]', '', predb['strain'][i]),
            'version': sub_disallowed(predb['version'][i]),
            'biosample': sub_disallowed(predb['biosample'][i]),
            'fna': predb['assemblyPath'][i],
            'gff3': predb['gffPath'][i],
            'source': predb['source'][i]
            }
        if predb['published'][i]:
            toAdd['published'] = predb['published'][i]
        else:
            toAdd['published'] = not predb['restriction'][i]
#        if predb['published'][i]:
 #           toAdd['published'] = predb['published'][i]
        infdb = infdb.append(toAdd)
    return infdb

def gen_omes(
    newdb, refdb = None, ome_col = 'ome', forbidden = set(),
    spacer = '\t'
    ):

    t_failed = []
    ref_ome_check = [k for k,v in Counter(refdb['ome']).items() if v > 1]
#    new_ome_check = [k for k,v in Counter(newdb['ome']).items() if v > 1 and k]
    if ref_ome_check:
        raise ValueError('corrupted reference database with non-unique omes: '
                        + str(ref_ome_check))
 #   elif new_ome_check:
  #      raise ValueError('corrupted pre-database with non-unique omes: '
   #                     + str(new_ome_check))
    tax_list = list(set(refdb['ome']).union(forbidden))
    tax_count = {}
    for tax in tax_list:
        abb = tax[:6]
        try:
            num = int(tax[6:])
        except ValueError: # version included in number
            num = int(re.search(r'(^\d+)', tax[6:])[1])
        if abb in tax_count:
            if tax_count[abb] < num:
                tax_count[abb] = num
        else:
            tax_count[abb] = num
    tax_count = Counter(tax_count)

    todel = []
    refdb_aas = refdb.set_index('assembly_acc')
    refdb_accs = set([x for x in list(refdb['assembly_acc']) if x])
    refdb_omes = set([x for x in list(refdb['ome'])])
    refdb_nover = {re.sub(r'(^.{6}\d+)\.\d+', r'\1', x): x \
                   for x in list(refdb['ome'])}
    for i, ome in enumerate(newdb['ome']):
        if not ome: # if there isn't an ome for this entry yet
            if newdb['assembly_acc'][i] in refdb_accs: # if this is an
            # established assembly accession; maybe make sure this works for
            # changed MycoCosm or NCBI assembly accs
                ome = refdb_aas[newdb['assembly_acc'][i]]['ome']
                v_search = re.search(r'\.(\d+)$', ome[6:])
                if v_search:
                    v = int(v_search[1]) + 1 # new version
                    new_ome = re.sub(r'\.\d+$', '.' + str(v), ome)
                else:
                    new_ome = ome + '.1' # first modified version
                eprint(spacer + ome + ' update -> ' + new_ome, flush = True)
                newdb['ome'][i] = new_ome
                continue
               
            try:
                name = newdb['genus'][i][:3].lower() + newdb['species'][i][:3].lower()
            except TypeError:
                todel.append(i)
                if not isinstance(newdb['assembly_acc'][i], float):
                    eprint(spacer + newdb['assembly_acc'][i] + ' no metadata - ' \
                         + 'failed', flush = True)
                elif 'index' in newdb: # for updateDB
                    if not isinstance(newdb, float):
                        eprint(spacer + newdb['index'][i] + ' no metadata - ' \
                            + 'failed', flush = True)
                    continue 
                else: # no use appending failed when there's no identifiable
                # info
                    continue
                row = {key: newdb[key][i] for key in mtdb.columns \
                       if key != 'ome'}
                t_failed.append(add2failed(row))
                continue
            name = re.sub(r'\(|\)|\[|\]|\$|\#|\@| |\||\+|\=|\%|\^|\&|\*|\'|\"|\!|\~|\`|\,|\<|\>|\?|\;|\:|\\|\{|\}', '', name)
            name.replace('(', '')
            name.replace(')', '')
            while len(name) < 6:
                name += '.'
            tax_count[name] += 1
            new_ome = name + str(tax_count[name])
            newdb['ome'][i] = new_ome
        elif ome:
            format_search = re.search(r'^[^\d_,\'";:\\\|\[\]\{\}\=\+\!@#\$\%\^' \
                                    + r'&\*\(\)]{6}\d+[^-_+=\\\|\{\[\}\]\:;' \
                                    + r'\'\"\,\<\>\?/\`\~\!\@\#\$\%\^\&\*\(' \
                                    + r'\)\w\W]*\.{0,1}\d*$', ome) # crude format check
            if not format_search:
                raise ValueError('invalid ome ' + ome)
            if ome in refdb_omes: # it's an update
                v_search = re.search(r'\.(\d+)$', ome[6:])
                if v_search:
                    v = int(v_search[1]) + 1 # new version
                    new_ome = re.sub(r'\.\d+$', '.' + str(v), ome)
                else:
                    new_ome = ome + '.1' # first modified version
                eprint(spacer + ome + ' update -> ' + new_ome, flush = True)
                newdb['ome'][i] = new_ome
            elif ome in refdb_nover: # has a version, wasn't given in predb
                version_ome = refdb_nover[ome]
                v_search = re.search(r'\.(\d+)$', version_ome[6:])
                if v_search:
                    v = int(v_search[1]) + 1 # new version
                    new_ome = ome + '.' + str(v)
                else:
                    raise TypeError('unknown error ' + ome)
                eprint(spacer + ome + ' version added -> ' + new_ome, flush = True)
                newdb['ome'][i] = new_ome
                    

    for i in reversed(todel):
        for key in mtdb.columns:
            del newdb[key][i]

    return newdb, t_failed

def cur_fna(cur_raw_fna_path, uncur_raw_fna_path, ome):
    ome_ver = re.search(r'(.{6}\d+).(\d+)$', ome)
    if ome_ver:
        less_ome = ome_ver[1]
        ver_num = ome_ver[2]
    else:
        less_ome = ome
        ver_num = 0
    with open(cur_raw_fna_path + '.tmp', 'w') as out:
        with open(uncur_raw_fna_path, 'r') as in_:
            for line in in_:
                if line.startswith('>'):
                    if not line.startswith('>' + ome + '_'):
                        if re.search(r'^>' + less_ome + '_', line) is not None:
                            new_line = line.replace('>' + less_ome + '_', '>' + ome + '_')
                        elif re.search(r'^>' + less_ome + r'\.\d+_', line) is not None:
                            new_line = re.sub(r'^>' + less_ome + r'\.\d+_', '>' + ome + '_',
                                              line)
                        else:
                            new_line = '>' + ome + '_' + line[1:]
                        out.write(new_line)
                    else: # already curated
                        out.write(line)
                else:
                    out.write(line)
    shutil.move(cur_raw_fna_path + '.tmp', cur_raw_fna_path)

def mmap_file_read(filename):
    """Read a file using memory mapping"""
    try:
        with open(filename, 'rb') as f:
            # Create memory map of file
            with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
                return mm.read().decode('utf-8')
    except (IOError, OSError) as e:
        raise IOError(f"Failed to memory map {filename}: {str(e)}")

def cur_mngr(ome, raw_fna_path, raw_gff_path, wrk_dir, 
            source, assembly_accession, exit=False,
            remove=False, spacer='\t\t\t', verbose=False,
            has_gff='no'):  # Add flag for GFF presence
    """Process individual genome files with optional GFF"""
    # Ensure working directory has trailing slash
    if not wrk_dir.endswith('/'):
        wrk_dir += '/'

    predb_dir = os.path.basename(os.path.dirname(wrk_dir[:-1])) + '/working/'

    # Construct paths with proper directory joining
    uncur_fna_path = os.path.join(wrk_dir, 'fna', f'{ome}.fna.uncur')
    cur_fna_path = os.path.join(wrk_dir, 'fna', f'{ome}.fna')
    uncur_gff_path = os.path.join(wrk_dir, 'gff3', f'{ome}.gff3.uncur')
    cur_gff_path = os.path.join(wrk_dir, 'gff3', f'{ome}.gff3')
    faa_path = os.path.join(wrk_dir, 'faa', f'{ome}.faa')
        
    # Process FNA files
    if not os.path.isfile(cur_fna_path):
        try:
            if not os.path.exists(raw_fna_path):
                eprint(f"\nERROR: Input FNA file does not exist: {raw_fna_path}")
                return ome, False, 'fna'
            
            # Use memory mapping for large FNA files
            if os.path.getsize(raw_fna_path) > 10_000_000:  # 10MB threshold
                fna_content = mmap_file_read(raw_fna_path)
                with open(uncur_fna_path, 'w') as f:
                    f.write(fna_content)
            else:
                uncur_fna_path = move_biofile(raw_fna_path, ome, 'fa', 
                                            wrk_dir + 'fna/', suffix = '.uncur')
        except (IOError, OSError) as ie:
            eprint(f"{spacer}{ome}|{assembly_accession} failed FNA parsing: {str(ie)}", 
                  flush=True)
            if exit:
                raise ie from None
            return ome, False, 'fna'

    # Only process GFF if it exists
    if has_gff == 'yes':
        if not os.path.isfile(cur_gff_path):
            try:
                # Use memory mapping for large GFF files
                if os.path.getsize(raw_gff_path) > 5_000_000:  # 5MB threshold
                    gff_content = mmap_file_read(raw_gff_path)
                    with open(uncur_gff_path, 'w') as f:
                        f.write(gff_content)
                    gff = gff2list(uncur_gff_path)
                else:
                    uncur_gff_path = move_biofile(raw_gff_path, ome, 'gff3', 
                                                wrk_dir + 'gff3/', suffix = '.uncur')
                    gff = gff2list(uncur_gff_path)
            except IOError as ie:
                eprint(f"{spacer}{ome}|{assembly_accession} failed GFF3 parsing: {str(ie)}", 
                      flush=True)
                if exit:
                    raise ie from None
                return ome, False, 'gff3'
            except IndexError:  # malformatted GFF
                return ome, False, 'gff3'
            except Exception as e:  # catch all other errors
                eprint(f"{spacer}{ome}|{assembly_accession} failed GFF3 curation: {str(e)}", 
                      flush=True)
                if exit:
                    raise e from None
                return ome, False, 'gff3'
    else:
        # Set empty paths for missing files
        cur_gff_path = ''
        faa_path = ''
        
    # Clean up temporary files if requested
    if remove:
        for path in [uncur_gff_path, raw_gff_path, uncur_fna_path, raw_fna_path]:
            if os.path.isfile(path):
                os.remove(path)
            # Also remove uncompressed versions
            uncompressed = re.sub(r'\.gz$', '', path)
            if os.path.isfile(uncompressed):
                os.remove(uncompressed)
    return ome, True, cur_fna_path, cur_gff_path, faa_path


def gff_mngr(ome, gff, cur_path, source, assembly_accession):

    gffVer, alias = None, False
    for entry in gff:
        if re.search(gff3Comps()['id'], entry['attributes']):
            gffVer = 3
            alias = re.search(gff3Comps()['Alias'], entry['attributes'])
            if alias is not None:
#                if entry['seqid'].startswith(ome + '_'):
                alias = True
                #else:
                 #   alias = False # remove the old aliases
                  #  entry['attributes'] = re.sub(r';?Alias=[^;]+', '',
                  #                              entry['attributes'])
            else:
                break
        elif re.search(gtfComps()['id'], entry['attributes']):
            gffVer = 2.5
            break
        elif re.search(gff2Comps()['id'], entry['attributes']):
            gffVer = 2
            break

    if gffVer == 3:
#        if source == 'new':
        if alias: #already curated
            try:
                new_gff = copy.deepcopy(gff)
                old_ome_p = re.search(gff3Comps()['Alias'], new_gff[0]['attributes'])[1]
                old_ome = old_ome_p[:old_ome_p.find('_')]
                for entry in new_gff:
#                    alias0 = re.search(gff3Comps()['Alias'], entry['attributes'])[1]
 #                   alias_num = alias0[alias0.find('_') + 1:]
#                    new_alias = ome + '_' + alias_num
 #                   entry['attributes'] = re.sub(
  #                      gff3Comps()['Alias'], 'Alias='+ new_alias,
   #                     entry['attributes']
    #                    )
                    entry['attributes'] = entry['attributes'].replace(old_ome, ome)
                new_gff = rename_and_organize(new_gff)
                gff = new_gff
            except:
                gff = curGFF3(gff, ome, cur_seqids = True)
#        else:
 #           gff = curGFF3(gff, ome)
        else:
            gff = curGFF3(gff, ome, cur_seqids = True)
    elif gffVer == 2.5:
        gff, trans_str, failed, flagged = gtf2gff3(gff, ome)
    else:
        gff, errors = gff2gff3(gff, ome, assembly_accession, verbose = False)

    ver_search = re.search(r'(.{6}\d+)\.(\d+)', ome)
    if ver_search is not None:
        less_ome = ver_search[1]
        ome_ver = ver_search[2]
    else:
        less_ome = ome
    for line in gff:
        seqid = line['seqid']
        if seqid.startswith(less_ome + '_'):
            line['seqid'] = re.sub(r'^' + less_ome + '_', 
                                 ome + '_', seqid)
        elif re.search(r'^' + less_ome + r'\.\d+_', seqid):
            line['seqid'] = re.sub(r'^' + less_ome + r'[^_]+_', 
                                 ome + '_', seqid)
        else:
            line['seqid'] = ome + '_' + seqid

    with open(cur_path + '.tmp', 'w') as out:
        out.write(list2gff(gff))
    shutil.move(cur_path + '.tmp', cur_path)

def add2failed(row):
    if isinstance(row['assembly_acc'], float) or not row['assembly_acc']:
        return False
    else:
        return [row['assembly_acc'], row['version']]

def batch_process_genomes(cur_cmds, max_cpus=None):
    """Process genomes in parallel with improved I/O handling and memory management"""
    n_cpus = min(max_cpus or mp.cpu_count(), mp.cpu_count())
    
    # Pre-validate all input files before processing
    valid_cmds = []
    failed_files = defaultdict(list)  # Track which files are missing for each genome
    
    for cmd in cur_cmds:
        ome, raw_fna, raw_gff, wrk_dir, *_ = cmd
        raw_fna = format_path(raw_fna)
        raw_gff = format_path(raw_gff)
        
        is_valid = True
        if not os.path.exists(raw_fna):
            failed_files[ome].append(('FNA', raw_fna))
            is_valid = False
        if not os.path.exists(raw_gff):
            failed_files[ome].append(('GFF', raw_gff))
            is_valid = False
            
        if is_valid:
            valid_cmds.append(cmd)
    
    # Report validation results
    if failed_files:
        eprint("\nValidation Failures:")
        for ome, failures in failed_files.items():
            eprint(f"\n{ome}:")
            for file_type, path in failures:
                eprint(f"  Missing {file_type}: {path}")
    
    if not valid_cmds:
        raise FileNotFoundError("No valid genomes to process - all input files missing")
    
    eprint(f"\nProcessing {len(valid_cmds)} valid genomes out of {len(cur_cmds)} total")
    
    # Use ProcessPoolExecutor for better exception handling
    results = []
    failed_processing = []
    with ProcessPoolExecutor(max_workers=n_cpus) as executor:
        futures = {executor.submit(cur_mngr, *cmd): cmd for cmd in valid_cmds}
        for f in tqdm(as_completed(futures), 
                     total=len(valid_cmds), 
                     desc='Processing genomes'):
            try:
                result = f.result()
                if result[1]:  # Check if processing was successful
                    results.append(result)
                else:
                    cmd = futures[f]
                    ome = cmd[0]
                    failed_processing.append((ome, result[2]))  # Store failure reason
            except Exception as e:
                cmd = futures[f]
                ome = cmd[0]
                failed_processing.append((ome, str(e)))
    
    # Report processing failures
    if failed_processing:
        eprint("\nProcessing Failures:")
        for ome, error in failed_processing:
            eprint(f"\n{ome}: {error}")
    
    eprint(f"\nSuccessfully processed {len(results)} out of {len(valid_cmds)} valid genomes")
    
    return results

def main(predb, refdb, wrk_dir, verbose=False, spacer='\t\t\t', 
         forbidden=set(), cpus=1, exit=False, remove=False):
    """Process predb files into MycotoolsDB format with optional GFF/FAA"""
    # Ensure working directory has trailing slash and is absolute
    wrk_dir = os.path.abspath(wrk_dir)
    if not wrk_dir.endswith('/'):
        wrk_dir += '/'

    # Create subdirectories with proper path joining
    for subdir in ['fna', 'gff3', 'faa']:
        dir_path = os.path.join(wrk_dir, subdir)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
            
    # Read the predb file and reference database
    predb_data = read_predb(predb)
    ref_db = mtdb(refdb) if refdb else mtdb(primaryDB())
    
    # Now pass the parsed predb data
    infdb = predb2mtdb(predb_data)
    
    vprint('\nGenerating omes', v=verbose, flush=True)
    omedb, failed = gen_omes(infdb, ref_db, ome_col='ome', 
                            forbidden=forbidden, spacer=spacer)    
    
    cur_cmds = []
    omedb = omedb.set_index('ome')
    for ome, row in omedb.items():
        cur_cmds.append([
            ome, row['fna'], row['gff3'], 
            wrk_dir, row['source'], row['assembly_acc'],
            exit, remove, spacer, verbose,
            row.get('has_gff', 'no')  # Add GFF flag to commands
        ])

    vprint('\nCurating data', v=verbose, flush=True)
    
    if cpus > 1:
        cur_data = batch_process_genomes(
            cur_cmds, 
            max_cpus=cpus
        )
    else:
        cur_data = []
        for cur_cmd in tqdm(cur_cmds, total=len(cur_cmds)):
            cur_data.append(cur_mngr(*cur_cmd))
            
    for data in cur_data:
        if not data[1]:
            failed.append(add2failed(omedb[data[0]]))
            del omedb[data[0]]
        else:
            ome, fna_path, gff3_path, faa_path = data
            omedb[ome]['fna'] = fna_path
            omedb[ome]['gff3'] = gff3_path if gff3_path else ''
            omedb[ome]['faa'] = faa_path if faa_path else ''
            omedb[ome]['has_gff'] = 'yes' if gff3_path else 'no'
            omedb[ome]['has_faa'] = 'yes' if faa_path else 'no'

    return omedb.reset_index(), failed


def cli():
    usage = 'Generate a predb file:\npredb2mtdb\n\nCreate a mycotoolsdb ' + \
    'from a predb file:\npredb2mtdb <PREDBFILE>\n\nCreate a mycotoolsdb ' + \
    'referencing an alternative master database:\npredb2mtdb <PREDBFILE> ' + \
    '<REFERENCEDB>\nSkip failing genomes:\npredb2mtdb <PREDBFILE> -s\n\n' + \
    'Control CPU usage:\npredb2mtdb <PREDBFILE> --cpus <NUMBER>\n' + \
    'Default behavior uses all available CPUs.'

    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('predb', nargs='?', help='Path to predb file')
    parser.add_argument('refdb', nargs='?', help='Reference database')
    parser.add_argument('-s', '--skip', action='store_true',
                       help='Skip failing genomes')
    parser.add_argument('--cpus', type=int, default=mp.cpu_count(),
                       help='Number of CPUs to use (default: all available)')
    args = parser.parse_args()

    if args.predb is None:
        print(gen_predb())
        sys.exit(0)

    wrk_dir = os.path.abspath('./predb2mtdb_working')
    if not wrk_dir.endswith('/'):
        wrk_dir += '/'
    
    main(args.predb, args.refdb, wrk_dir=wrk_dir, 
         cpus=args.cpus, exit=not args.skip)
    
if __name__ == '__main__':
    cli()
