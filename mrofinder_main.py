import argparse
import subprocess
import os
import re
import math
from Bio import SeqIO, SearchIO
# from glob import glob as glob

import mrofinder_classes as classes
import mrofinder_helpers as helpers


def main():
    options = parser()
    proteome = fasta2proteome(options)
    proteome.add_homology(manage_hmmers(options), 'hmmer')
    proteome.add_homology(manage_mitominer(options), 'mitominer')
    manage_targeting(proteome, options)
    manage_structure(proteome, options)
    manage_interproscan(proteome, options)
    manage_blast_nr(proteome, options)
    proteome.annotate()
    save_table(options, proteome)


def parser():
    """ Parse arguments"""
    parser = argparse.ArgumentParser()
    # general options
    parser.add_argument('-i', '--input', required=True, help='Analysed protein fasta')
    parser.add_argument('-o', '--output', required=True, help='Output file (table with all informations)')
    parser.add_argument('-w', '--working_directory', help='Working directory to be used during analyses.')
    parser.add_argument('--working_fasta', default='working.fasta',
                        help='Specify name of working fasta file (used during analyses to ensure no naming issues.')
    parser.add_argument('-c', '--continue', action='store_true',
                        help='Use existing working directory with partial results.')
    parser.add_argument('-t', '--threads', type=int, help='Number of threads to be used (only in several subprocesses)',
                        default=1)
    parser.add_argument('--test', action='store_true')
    # homology options
    hmmer = parser.add_mutually_exclusive_group(required=True)
    hmmer.add_argument('--hmmer_results', nargs='+',
                       help='Directory containing results of conserved proteins hmmer search.')
    hmmer.add_argument('--hmmer_queries', nargs='+',
                       help='Directory with profiles/profiles of conserved proteins to be used in hmmer search.')
    parser.add_argument('--q2db_blast', nargs=1, help='Proteome to database blast file')
    parser.add_argument('--db2q_blast', nargs=1, help='Database to proteome blast file')
    parser.add_argument('--mitominer_fasta', help='Mitominer+ fasta after cdhit')
    parser.add_argument('--mitominer_db', help='Mitominer+ blast database')
    parser.add_argument('--interproscan', nargs=1, help='Interproscan output (run externally).')
    blast = parser.add_mutually_exclusive_group()
    blast.add_argument('--blast_nr', nargs=1, help='Blast results in .xml format. (run externally)')
    blast.add_argument('--ncbi_nr_db', nargs=1, help='NCBI nr database', default='/mnt/databases/NCBI/nr/nr')
    # alpha/beta options
    parser.add_argument('--tmhmm', nargs='+', help='File(s) with tmhmm results.')
    parser.add_argument('--psipred_results', nargs=1, help='Directory with psipred results, names need to be \
                        PROTEINNAME.fasta.ss2 (psipred\'s defaults).')
    parser.add_argument('--pattern', default='1', choices=['0', '1', '2', '3', '4'],
                        help='0 = PoxGxxHyxHy, 1 = PoHyGhyxHyxHy (default), 2 = PoHyGhyH^yHyxHy, 3 = H^yHyGhyxHyxHy, \
                        4 = H^yHyGhyH^yHyxHy')
    # targeting options
    parser.add_argument('--targetp', nargs='+', help='Result(s) of targetp.')
    parser.add_argument('--mitofates', nargs=1, help='Results of mitofates')
    parser.add_argument('--nommpred', nargs=1, help='Results of nommpred')
    # parse and check arguments
    parsed = parser.parse_args()
    check = 0
    if parsed.q2db_blast and parsed.db2q_blast:
        check += 1
    if parsed.mitominer_fasta and parsed.mitominer_db:
        check += 1
    if check != 1:
        parser.print_help()
        raise Exception('error: (--q2db_blast and --db2q_blast) xor (--mitominer_fasta and --mitominer_db) required')
    parsed.working_directory = helpers.check_working_directory(parsed.working_directory, parsed.test)
    return parsed


# LOAD INPUT
def fasta2proteome(options):
    proteome = classes.Proteome()
    old2new_name_dict = helpers.fix_fasta(options)
    # TODO suppport for external files using 'normal' fasta
    seqs = SeqIO.parse(options.working_fasta, 'fasta')
    for seq in seqs:
        proteome.add_protein(seq.id, classes.Protein(old2new_name_dict[seq.id], seq.id, str(seq.seq)))
    proteome.add_names_dict(old2new_name_dict)
    proteome.get_helpers()
    return proteome


# HOMOLOGY
def manage_hmmers(options):
    print('Managing hmmers ...')
    if options.hmmer_results:
        list_of_results = helpers.check_if_dir_or_file(options.hmmer_results, '\.out$')
        results = parse_hmmers(list_of_results)
    else:
        list_of_queries = helpers.check_if_dir_or_file(options.hmmer_queries)
        list_of_results = run_hmmers(list_of_queries, options)
        results = parse_hmmers(list_of_results)
    return results


def parse_hmmers(list_of_files):
    results = classes.HomologyResults()
    for handle in list_of_files:
        best_pevalue = 0
        hmmers = SearchIO.parse(handle, 'hmmer3-text')
        for queryresult in hmmers:
            for hit in queryresult:
                hmmer_name = queryresult.id
                protein_name = hit.id
                pevalue = -math.log(float(hit.evalue))
                if not best_pevalue:
                    best_pevalue = pevalue
                    results.add_hit(protein_name, hmmer_name)
                elif best_pevalue < pevalue:
                    best_pevalue = pevalue
                    results.add_hit(protein_name, hmmer_name)
                elif best_pevalue <= 2 * pevalue:
                    results.add_hit(protein_name, hmmer_name)
    return results


def run_hmmers(list_of_queries, options):
    blast_wd = os.path.join(options.working_directory, 'hmmers')
    os.mkdir(blast_wd)
    list_of_results = []
    for handle in list_of_queries:
        name = os.path.basename(handle).split('.')[0]
        fasta_name = helpers.get_basename(options.input)
        output_name = os.path.join(blast_wd, fasta_name + '_' + name + '.out')
        with open(output_name, 'w') as output_file:
            _p = subprocess.run([paths['hmmsearch'], handle, options.working_fasta], stdout=output_file)
        list_of_results.append(output_name)
    return list_of_results


def manage_mitominer(options):
    if options.q2db_blast and options.db2q_blast:
        q2db_blast = options.q2db_blast[0]
        db2q_blast = options.db2q_blast[0]
    else:
        q2db_blast, db2q_blast = run_mitominer_blast(options)
    q2db_dict = helpers.parse_blast_tabular(q2db_blast)
    db2q_dict = helpers.parse_blast_tabular(db2q_blast)
    return helpers.search_for_bbhs(q2db_dict, db2q_dict)


def run_mitominer_blast(options):
    blast_wd = os.path.join(options.working_directory, 'mitominer_blast')
    os.mkdir(blast_wd)
    q2db_blast = os.path.join(blast_wd, 'q2db_blast.tsv')
    db2q_blast = os.path.join(blast_wd, 'db2q_blast.tsv')
    db_path = check_if_db(options.working_fasta)
    if not db_path:
        db_path = os.path.join(blast_wd, os.path.splitext(os.path.basename(options.working_fasta))[0])
        _s = subprocess.run([os.path.join(paths['blast'], 'makeblastdb'), '-dbtype', 'prot', '-in',
                             options.working_fasta, '-out', db_path])
    _p = subprocess.run([os.path.join(paths['blast'], 'blastp'), '-num_threads', str(options.threads), '-query',
                         options.working_fasta, '-db', options.mitominer_db, '-outfmt', '6', '-evalue', '0.001', '-out',
                         q2db_blast])
    _q = subprocess.run([os.path.join(paths['blast'], 'blastp'), '-num_threads', str(options.threads), '-query',
                         options.mitominer_fasta, '-db', db_path, '-outfmt', '6', '-evalue', '0.001', '-out',
                         db2q_blast])
    return q2db_blast, db2q_blast


def check_if_db(path):
    directory, file = os.path.split(path)
    rootname = os.path.splitext(file)[0]
    try1 = os.path.join(directory, rootname)
    if os.path.exists(try1 + '.phr') and os.path.exists(try1 + '.psq') and os.path.exists(try1 + '.pin'):
        return try1
    try2 = try1 + '.fasta'
    if os.path.exists(try2 + '.phr') and os.path.exists(try2 + '.psq') and os.path.exists(try2 + '.pin'):
        return try2
    return False


# TARGETING
def manage_targeting(proteome, options):
    list_of_targetp_results, list_of_mitofates_results = [], []
    # targetp
    if options.targetp:
        list_of_targetp_results = helpers.check_if_dir_or_file(options.targetp, '.targetp')
    elif paths['targetp']:
        list_of_targetp_results = run_targetp(options)
    else:
        raise Exception('Targetp is required')
    proteome.add_targeting(parse_targetp(list_of_targetp_results), 'targetp')
    # mitofates
    if options.mitofates:
        list_of_mitofates_results = helpers.check_if_dir_or_file(options.mitofates)
    elif paths['mitofates']:
        list_of_mitofates_results = run_mitofates(options)
    if list_of_mitofates_results:
        proteome.add_targeting(parse_mitofates(list_of_mitofates_results), 'mitofates')
    # nommpred
    if options.nommpred:
        proteome.add_targeting(parse_nommperd(options.nommpred[0]), 'nommpred')


def run_targetp(options):
    print('Running targetp...')
    targetp_wd = os.path.join(options.working_directory, 'targetp')
    os.mkdir(targetp_wd)
    smaller_fasta_files = helpers.fasta_for_targetp(options.working_fasta, targetp_wd)
    results = []
    for small_fasta in smaller_fasta_files:
        small_fasta_basename = helpers.get_basename(small_fasta)
        output_name = helpers.get_output_name(targetp_wd, small_fasta_basename, '.targetp')
        with open(output_name, 'w') as output_file:
            _p = subprocess.run([paths['targetp'], '-N', '-t', '0.50', small_fasta], stdout=output_file)
        results.append(output_name)
    return results


def parse_targetp(targetp_file_list):
    results = classes.SingularResults()
    for targetp_file in targetp_file_list:
        with open(targetp_file) as file:
            for i in range(8):
                _header = next(file)
            for line in file:
                if not line[0] == '-' and not line[-2] == '-':
                    line = line.split()
                    try:
                        name, probability, localisation = line[0], line[2], line[5]
                    except Exception as exc:
                        print(line)
                        raise exc
                    results.add_hit(name, (localisation, probability))
                else:
                    break
    return results


def run_mitofates(options):
    print('Running mitofates...')
    mitofates_wd = os.path.join(options.working_directory, 'mitofates')
    os.mkdir(mitofates_wd)
    output_name = helpers.get_output_name(mitofates_wd, options.input, '.mitofates')
    with open(output_name, 'w') as output_file:
        _p = subprocess.run([paths['mitofates'], options.working_fasta, 'fungi'], stdout=output_file)
    return [output_name]


def parse_mitofates(list_of_files):
    results = classes.SingularResults()
    for mitofates_file in list_of_files:
        with open(mitofates_file) as file:
            _header = next(file)
            for line in file:
                line = line.split('\t')
                try:
                    name, probability, localisation = line[0], line[1], ''
                except Exception as bla:
                    print(line)
                    print(list_of_files)
                    raise bla
                if line[2] == 'Possessing mitochondrial presequence':
                    localisation = 'M'
                elif line[2] == 'No mitochondrial presequence':
                    localisation = '_'
                else:
                    raise ValueError('Mitofates prediction incorrect')
                results.add_hit(name, (localisation, probability))
    return results


def parse_nommperd(nommpred_file):
    results = classes.SingularResults()
    with open(nommpred_file) as file:
        _header = next(file)
        for line in file:
            name, raw_localisation = line.split()
            if raw_localisation == 'Other':
                localisation = '_'
            elif raw_localisation == 'MRO':
                localisation = 'M'
            else:
                raise ValueError('Nommpred prediction incorrect')
            results.add_hit(name, (localisation, 'N/A'))
    return results


# ALPHA/BETA STRUCTURE
def manage_structure(proteome, options):
    manage_tmhmm(proteome, options)
    manage_beta_structure(proteome, options)


# ALPHA
def manage_tmhmm(proteome, options):
    if options.tmhmm:
        list_of_tmhmm_results = helpers.check_if_dir_or_file(options.tmhmm)
    else:
        list_of_tmhmm_results = run_tmhmm(options)
    proteome.add_tmhmm(parse_tmhmm(list_of_tmhmm_results))


def parse_tmhmm(list_of_files):
    results = classes.SingularResults()
    for tmhmm_file in list_of_files:
        with open(tmhmm_file) as file:
            for line in file:
                name, length, expaa, first60, predehel, topology = line.split()
                helices = []
                pre_topology = topology.split('=')[1]
                numbers = re.findall('(\d+)', pre_topology)
                while numbers:
                    beginning, end = int(numbers.pop(0)) - 1, int(numbers.pop(0)) - 1
                    helices.append(classes.TransMembraneHelix(beginning, end))
                results.add_hit(name, helices)
    return results


def run_tmhmm(options):
    print('Running tmhmm...')
    tmhmm_wd = os.path.join(options.working_directory, 'tmhmm')
    os.mkdir(tmhmm_wd)
    output_name = helpers.get_output_name(tmhmm_wd, options.input, '.tmhmm')
    with open(output_name, 'w') as output_file:
        _p = subprocess.run([paths['tmhmm'], '--short', options.working_fasta], stdout=output_file)
    return [output_name]


# BETA
def manage_beta_structure(proteome, options):
    print('Managing psipred ...')
    psipred_wd = os.path.join(options.working_directory, 'psipred')
    os.mkdir(psipred_wd)
    candidate_list = search_for_betasignal(proteome, helpers.choose_pattern(options.pattern))
    if options.psipred_results:
        list_of_psipred_files = helpers.check_if_dir_or_file(options.psipred_results, '\.ss2$')
    else:
        list_of_psipred_files = []
    beta_proteins_lists = run_and_analyse_psipred(list_of_psipred_files, candidate_list, psipred_wd)
    proteome.add_beta_signal(beta_proteins_lists)


def search_for_betasignal(proteome, pattern):
    hit_list = []
    for key, protein in proteome.proteins.items():
        start, end = helpers.re_search(pattern, protein.seq)
        if end - start >= 90 and not protein.helices:
            hit_list.append(classes.BetaSignalHit(protein.work_id, protein.seq, start, end))
    return hit_list


def run_and_analyse_psipred(list_of_psipred_files, candidates_list, working_directory):
    list_of_beta_proteins = []
    if list_of_psipred_files:
        prefix = os.path.dirname(list_of_psipred_files[0])
    else:
        prefix = ''
    for hit in candidates_list:
        ss2_file = os.path.join(prefix, hit.name + '.fasta.ss2')
        if ss2_file not in list_of_psipred_files:
            subseq_fasta_path = os.path.join(working_directory, hit.name + '.fasta')
            ss2_file = subseq_fasta_path + '.ss2'
            with open(subseq_fasta_path, 'w') as subseq_fasta_file:
                subseq_fasta_file.write('>' + hit.name + '\n')
                if len(hit.seq) <= 1500:
                    subseq_fasta_file.write(hit.seq + '\n')
                else:
                    subseq_fasta_file.write(hit.seq[max([0, hit.end - 1500]):hit.end])
            _p = subprocess.run([paths['psipred'], subseq_fasta_path, '-d', databases['psipred'], '-o',
                                 working_directory])
        try:
            if parse_psipred(ss2_file, hit.start, hit.end):
                list_of_beta_proteins.append(hit.name)
        except Exception as exc:
            print(working_directory)
            print(prefix)
            print(ss2_file)
            print(list_of_psipred_files)
            raise exc
    return list_of_beta_proteins


def parse_psipred(ss2_file, start, end):
    length = end - start
    with open(ss2_file) as file:
        _header = next(file)
        _blank_line = next(file)
        structure = []
        counter = {'C': 0, 'H': 0, 'E': 0}
        cnt = 0
        for line in file:
            if start <= cnt < end:
                letter = line.split()[2]
                structure.append(letter)
                counter[letter] += 1
            cnt += 1
        signal_counter = 0
        for i in structure[-8:]:
            if i == 'H':
                signal_counter += 1
        if counter['H'] / length <= 0.1 \
                and counter['E'] / length >= 0.25 \
                and signal_counter < 5:
            return True
        else:
            return False
        # as for 12.07.19 I think removal is stupid, but maybe I'll change mind in future
        # else:
        #     for file in glob(hit.name + '*'):
        #         os.remove(file)


# INTERPROSCAN

def manage_interproscan(proteome, options):
    print('Managing interproscan...')
    if not options.interproscan:
        options.interproscan = run_interproscan(options)
    go_dictionary = helpers.get_go_dictionary(paths['go_basic'])
    proteome.add_go_categories(parse_interproscan(options.interproscan), go_dictionary)


def run_interproscan(options):
    interproscan_wd = os.path.join(options.working_directory, 'interproscan')
    os.mkdir(interproscan_wd)
    output_path = helpers.get_output_name(interproscan_wd, options.input, '.interproscan')
    # TODO finish proper running
    _p = subprocess.run([paths['interproscan'], '-i', options.input, '-o', output_path, '--goterms', '-f', 'tsv',
                         '--cpu', str(options.threads)])
    return output_path


def parse_interproscan(interproscan_file, ):
    protein2go = {}
    with open(interproscan_file) as file:
        for line in file:
            line = line.strip().split('\t')
            if len(line) == 14:
                name, gos = line[0], line[13].split('|')
                if name not in protein2go.keys():
                    protein2go[name] = []
                for go in gos:
                    if go not in protein2go[name]:
                        protein2go[name].append(go)
    return protein2go


# blast_nr

def manage_blast_nr(proteome, options):
    print('Managing blast nr...')
    if options.blast_nr:
        blast_nr_file = options.blast_nr
    elif options.ncbi_nr_db:
        blast_nr_file = run_blast_nr(options)
    else:
        return
    blast_nr_dict = parse_blast_nr(blast_nr_file)
    proteome.add_blast_nr_results(blast_nr_dict)


def run_blast_nr(options):
    blast_nr_wd = os.path.join(options.working_directory, 'blast_nr')
    os.mkdir(blast_nr_wd)
    output_path = helpers.get_output_name(blast_nr_wd, options.input, '2nr.xml')
    _p = subprocess.run([os.path.join(paths['blast'], 'blastp'), '-num_threads', str(options.threads), '-query',
                         options.working_fasta, '-db', options.ncbi_nr_db, '-outfmt', '5', '-evalue', '0.001', '-out',
                         output_path])
    return output_path


def parse_blast_nr(blast_nr_file):
    blast_dict = {}
    results = SearchIO.parse(blast_nr_file, 'blast-xml')
    for result in results:
        hit = helpers.get_hit(result)
        if hit == 'flag':
            blast_dict[result.id] = ['-', '_just_hypothetical', '-']
        elif hit:
            blast_dict[result.id] = [hit.id, hit.description, helpers.get_evalue(hit)]
        else:
            blast_dict[result.id] = ['-', '_no_hit', '-']
    return blast_dict


# saving results

def save_table(options, proteome):
    with open(options.output, 'w') as file:
        # TODO naming columns
        for key, protein in proteome.proteins.items():
            file.write(str(protein) + '\n')


if __name__ == '__main__':
    paths, databases = helpers.parse_paths()
    main()
