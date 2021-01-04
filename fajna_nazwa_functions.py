from datetime import datetime
import os
from glob import glob as glob
import re
import yaml
import fajna_nazwa_classes as classes
from shutil import rmtree
from Bio import SeqIO


def parse_paths():
    dirname = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(dirname, 'fajna_nazwa_paths.txt')
    with open(path) as file:
        parsed = yaml.safe_load(file)
    return parsed['Software paths/calls'], parsed['Databases']


def parse_targetp_tabular(file):
    list_of_hits = []
    for line in file:
        line = line.split() # Name, Len, mTP, SP, other, Loc, RC
        try:
            length, mTP, SP, other, RC = int(line[1]), int(line[2]), int(line[3]), int(line[4]), int(line[6])
            name, localisation = line[0], line[5]
            list_of_hits.append([])
        except ValueError:
            continue
    return list_of_hits


def check_working_directory(path, test_flag, continue_flag=False):
    if not path:
        currentDT = datetime.now()
        helper_list = [str(currentDT.year), str(currentDT.month), str(currentDT.day), str(currentDT.hour),
                       str(currentDT.minute)]
        for i in range(len(helper_list)):
            while len(helper_list[i]) < 2:
                helper_list[i] = '0' + helper_list[i]
        path = 'working_directory_' + '_'.join(helper_list)
    if os.path.exists(path) and os.path.isdir(path):
        if not os.listdir(path):
            pass
        else:
            if continue_flag:  # TODO continue option
                subdirs = os.listdir(path)
            if test_flag:
                rmtree(path)
                os.mkdir(path)
            else:
                raise IOError('working directory not empty')
    else:
        os.mkdir(path)
    return path


def check_if_dir_or_file(list_of_paths, pattern=None):
    list_of_dewildcarded_paths = []
    for path in list_of_paths:
        if '*' in path:
            list_of_paths.extend(glob(path))
        else:
            list_of_dewildcarded_paths.append(path)
    list_of_files = []
    for path in list_of_dewildcarded_paths:
        if os.path.isdir(path):
            for file in os.listdir(path):
                if pattern:
                    if re.search(pattern, file):
                        list_of_files.append(os.path.join(path, file))
                else:
                    list_of_files.append(os.path.join(path, file))
        elif os.path.isfile(path):
            list_of_files.append(path)
    if not list_of_files:
        raise Exception('No proper files found in {}'.format(list_of_paths))
    return list_of_files


def fix_fasta(options):
    new_fasta_path = os.path.join(options.working_directory, options.working_fasta)
    dict_path = os.path.join(options.working_directory, get_basename(new_fasta_path) + '_dict.tsv')
    new2old = {}
    with open(options.input) as _input, open(new_fasta_path, 'w') as _output,\
            open(dict_path, 'w') as _dict:
        cnt = 1
        for line in _input:
            if line[0] == '>':
                old_name = '_'.join(line.strip()[1:].split())
                new_name = 'gene_' + str(cnt)
                _output.write('>' + new_name + '\n')
                _dict.write('\t'.join([old_name, new_name]) + '\n')
                if new_name in new2old.keys():
                    raise KeyError('At least two proteins have the same name!')
                else:
                    new2old[new_name] = old_name
                cnt += 1
            else:
                _output.write(line)
    options.working_fasta = new_fasta_path
    return new2old


def add_annotation(protein_name, hmmer_name, dict_of_annotations):
    if protein_name in dict_of_annotations:
        dict_of_annotations[protein_name].append(hmmer_name)
    else:
        dict_of_annotations[protein_name] = [hmmer_name]


def parse_blast_tabular(handle):
    """Parses tabular blast file, returns list of hits with evalue <= best_evalue for each query."""
    hits_dictionary = {}
    with open(handle) as file:
        working_qname, working_hits, best_evalue = '', [], 1.0
        for line in file:
            qname, sname, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.split()
            evalue = float(evalue)
            if qname == working_qname:
                if evalue < 10 * best_evalue:
                    working_hits.append(sname)
                    if evalue < best_evalue:
                        best_evalue = evalue
            else:
                hits_dictionary[working_qname] = working_hits
                working_qname = qname
                best_evalue = evalue
                working_hits = [sname]
        else:
            hits_dictionary[working_qname] = working_hits
    return hits_dictionary


def search_for_bbhs(q2db_dict, db2q_dict):
    """Checks if query is repriprocal hit of its own hits."""
    final_dict = classes.HomologyResults()
    q2db_keys, db2q_keys = q2db_dict.keys(), db2q_dict.keys()
    for key in q2db_keys:
        hits = q2db_dict[key]
        for hit in hits:
            if hit in db2q_keys:
                if key in db2q_dict[hit]:
                    final_dict.add_hit(key, hit)
    return final_dict


def fasta_for_targetp(fasta, working_directory, max_chars = 200000):
    results = []
    out_base = os.path.join(working_directory, get_basename(fasta))
    records = parse_fasta_for_targetp(fasta)  # [(header, seq), ...]
    count, out_num, print_buffer = 0, 1, []
    for i in range(len(records)):
        header, seq = records[i]
        print_buffer.append(header)
        seq = seq[:120]
        count += len(seq)
        print_buffer.append(seq)
        if count > max_chars - 120:
            output_name = out_base + '_' + str(out_num) + '.fasta'
            print_fasta_for_targetp(output_name, print_buffer)
            results.append(output_name)
            out_num += 1
            print_buffer = []
            count = 0
    if print_buffer:
        output_name = out_base + '_' + str(out_num) + '.fasta'
        print_fasta_for_targetp(output_name, print_buffer)
        results.append(output_name)
    return results


def parse_fasta_for_targetp(file):
    seqs = []
    header, seq = '', ''
    with open(file) as _input:
        for line in _input:
            if line[0] == ">":
                seqs.append((header, seq))
                header, seq = line.strip(), ''
            else:
                seq += line.strip()
        if header:
            seqs.append((header, seq))
    return seqs


def print_fasta_for_targetp(output, print_buffer):
    with open(output, 'w') as _output:
        for line in print_buffer:
            _output.write(line + '\n')


def get_basename(path):
    return os.path.splitext(os.path.basename(path))[0]


def get_output_name(working_directory, input_path, extension):
    return os.path.join(working_directory, get_basename(input_path) + extension)


# TMHMM IMPLEMENTATION
def calculate_tmd_gravy(seq):
    score = 0
    length = len(seq)
    for i in seq:
        score += hydrophobicity(i)
    return score/length


def calculate_charge(seq):
    score = 0
    for i in seq:
        score += charge(i)
    return score


def hydrophobicity(amino_acid):  # Kyte (1982) A simple method for displaying the hydropathic character of a protein
    hydrophobicity_dict = {'I': 4.5,
                           'V': 4.2,
                           'L': 3.8,
                           'F': 2.8,
                           'C': 2.5,
                           'M': 1.9,
                           'A': 1.8,
                           'X': 0,
                           '*': 0,
                           'G': -0.4,
                           'T': -0.7,
                           'S': -0.8,
                           'W': -0.9,
                           'Y': -1.3,
                           'P': -1.6,
                           'H': -3.2,
                           'E': -3.5,
                           'Q': -3.5,
                           'D': -3.5,
                           'N': -3.5,
                           'K': -3.9,
                           'R': -4.5}
    return hydrophobicity_dict[amino_acid]


def charge(amino_acid):
    charge_dict = {'H': 1,
                   'K': 1,
                   'R': 1,
                   'D': -1,
                   'E': -1}
    if amino_acid in charge_dict.keys():
        return charge_dict[amino_acid]
    else:
        return 0


def re_search(pattern, seq):
    found = None
    for found in re.finditer(pattern, str(seq)):
        pass
    if found:
        return max(0, found.start() - 300), found.end()
    else:
        return 0, 0


def choose_pattern(pattern):
    help_dict = {}
    Po = 'KRHSTNQ'
    hy = 'ACVLIMFYW'
    Hy = 'VLIMFYW'
    help_dict['0'] = '[X{0}].G..[X{1}].[X{2}]'.format(Po, hy, Hy)  # PoxGxxHyxHy
    help_dict['1'] = '[X{0}][X{2}]G[X{1}].[X{2}].[X{2}]'.format(Po, hy, Hy)  # PoHyGhyxHyxHy
    help_dict['2'] = '[X{0}][X{2}]G[X{1}][^{2}][X{2}].[X{2}]'.format(Po, hy, Hy)  # PoHyGhyH^yHyxHy
    help_dict['3'] = '[^{2}][X{2}]G[X{1}].[X{2}].[X{2}]'.format(Po, hy, Hy)  # H^yHyGhyxHyxHy
    help_dict['4'] = '[^{2}][X{2}]G[X{1}][^{2}][X{2}].[X{2}]'.format(Po, hy, Hy)  # H^yHyGhyH^yHyxHy
    return help_dict[pattern]
