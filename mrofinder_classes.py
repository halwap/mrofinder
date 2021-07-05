import mrofinder_helpers as helpers
import warnings


class Proteome:
    def __init__(self):
        self.name = ''
        self.proteins = {}
        self.old2new_names_dict = {}
        self.new2old_names_dict = {}
        self.proteins_keys = []
        # TODO dodac sprawdzanie co bylo zannotowane, zeby inteligentnie printowac tabelke

    def add_protein(self, working_name, protein):
        self.proteins[working_name] = protein

    # def add_uniprot(self, ?):
    #     pass

    # def add_pfam(self):
    #     pass

    def add_homology(self, results, software):
        """expected arguments
        results: SearchResult object
        type: 'hmmer', 'mitominer'
        """
        if software == 'hmmer':
            for protein_name in results.keys():
                if protein_name not in self.proteins.keys():
                    protein_name = self.old2new_names_dict[protein_name]
                self.proteins[protein_name].hmmer = results.hits(protein_name)
        elif software == 'mitominer':
            for protein_name in results.keys():
                if protein_name not in self.proteins.keys():
                    protein_name = self.old2new_names_dict[protein_name]
                self.proteins[protein_name].mitominer_bbd = results.hits(protein_name)
        else:
            raise ValueError('Accepted types of homology: hmmer, mitominer')

    def add_tmhmm(self, parsed_tmhmm):
        """takes SingularResults with list of TransMembraneHelix(es) and adds to proteins"""
        for protein_name in parsed_tmhmm.keys():
            if protein_name not in self.proteins.keys():
                protein_name = self.old2new_names_dict[protein_name]
            self.proteins[protein_name].add_helices(parsed_tmhmm.hits(protein_name))

    def add_targeting(self, parsed_targeting, software):
        """
        parsed_targeting: LocalisationResults
        software in [targetp, mitofates, nommpred]
        """
        if software in ['targetp', 'mitofates', 'nommpred']:
            for hit in parsed_targeting.keys():  # Name, Len, mTP, SP, other, Loc, RC
                self.proteins[hit].set_targeting(parsed_targeting.hits(hit), software)
        else:
            raise ValueError("Software not in [targetp, mitofates, nommpred]")

    def add_beta_signal(self, list_of_positives):
        for protein_name in list_of_positives:
            if protein_name not in self.proteins.keys():
                protein_name = self.old2new_names_dict[protein_name]
            protein = self.proteins[protein_name]
            if protein.length >= 90:
                protein.beta_signal = True

    def add_go_categories(self, protein2go, go_dictionary):
        print('add_go_categories')
        for protein_name in self.proteins_keys:
            go_ids, go_names = [], []
            if protein_name in protein2go.keys():
                print(protein_name, protein2go[protein_name])
                go_ids = protein2go[protein_name]
                go_names = []
                for go_id in go_ids:
                    print(go_id)
                    if go_id in go_dictionary.keys():
                        go_names.append(go_dictionary[go_id])
            self.proteins[protein_name].add_go(go_ids, go_names)

    def add_blast_nr_results(self, blast_dictionary):
        for protein_name in self.proteins_keys:
            keys = blast_dictionary.keys()
            if protein_name in keys:
                blast_hit = blast_dictionary[protein_name]
            else:
                blast_hit = ['-', '_no_hit', '-']
            self.proteins[protein_name].add_blast_nr(blast_hit)

    def annotate(self):
        for key, protein in self.proteins.items():
            protein.annotate()

    def add_names_dict(self, old2new_names_dict):
        self.old2new_names_dict = old2new_names_dict
        self.new2old_names_dict = {v:[k for k in old2new_names_dict if old2new_names_dict[k] == v]
                                   for v in old2new_names_dict.values()}

    def get_helpers(self):
        self.proteins_keys = self.proteins.keys()


class Protein:
    def __init__(self, name: str, working_name: str, seq: str):
        self.id = name
        self.work_id = working_name
        self.seq = seq
        self.length = len(seq)
        # self.uniprot_hit = ''
        # self.pfam_hit = ''
        self.hmmer = []
        self.mitominer_bbd = []
        self.targetp = None
        self.mitofates = None
        self.nommpred = None
        self.helices = None
        self.tmd_gravy = 0
        self.end_charge = 0
        self.beginning_charge = 0
        self.beta_signal = False
        self.go_categories = []
        self.go_descriptions = []
        self.ncbi_nr_hit = []
        # classification
        self.tail_anchored = False
        self.mbomp = False
        self.targeted_by = False
        self.homology = False
        self.interesting = False

    def __str__(self):
        list_to_print = [self.id]
        if self.hmmer:
            list_to_print.append(';'.join(self.hmmer))
        else:
            list_to_print.append('_no_hmmer_hits')
        if self.mitominer_bbd:
            list_to_print.append(';'.join(self.mitominer_bbd))
        else:
            list_to_print.append('_no_mitominer_bbds')
        if self.go_categories:
            print(type(go_categories[0]), type(go_descriptions[0]))
            list_to_print.append('|'.join(self.go_categories))
            list_to_print.append('|'.join(self.go_descriptions))
        else:
            list_to_print.append('_no_go_hits')
            list_to_print.append('_no_go_hits')
        list_to_print.extend(self.ncbi_nr_hit)
        list_to_print.append(str(self.targeted_by))
        list_to_print.extend([self.targetp.localisation, self.targetp.probability, self.mitofates.localisation,
                             self.mitofates.probability])
        if self.nommpred:
            list_to_print.append(self.nommpred.localisation)
        if self.tail_anchored:
            list_to_print.extend(['_TA', str(self.helices[0].length), str(self.tmd_gravy), str(self.beginning_charge),
                                  str(self.end_charge)])
        else:
            list_to_print.extend(['_not_TA', 'N/A', 'N/A', 'N/A', 'N/A'])
        if self.mbomp:
            list_to_print.append('_MBOMP')
        else:
            list_to_print.append('_not_MBOMP')
        return '\t'.join(list_to_print)

    def set_seq(self, seq: str):
        self.seq = seq

    # def set_uniprot_hit(self, uniprot_hit):
    #     self.uniprot_hit = uniprot_hit

    def set_targeting(self, targeting: tuple, software: str):
        """
        Targeting: (localisation, probability)
        software: targetp, mitofates, nommpred (others to be implemented soonTM)
        """
        if software == 'targetp':
            self.targetp = Targeting(targeting[0], targeting[1])
        elif software == 'mitofates':
            self.mitofates = Targeting(targeting[0], targeting[1])
        elif software == 'nommpred':
            self.nommpred = Targeting(targeting[0], targeting[1])
        else:
            raise ValueError("Software not in [targetp, mitofates, nommpred]")
        # if targeting[0] == 'M':  # ???
        #     pass

    def add_helices(self, helices):
        self.helices = helices
        if len(helices) == 1:
            tmd_subseq = self.seq[helices[0].beginning:helices[0].end]
            end_subseq = self.seq[helices[0].end:]
            beginning_subseq = self.seq[helices[0].beginning - len(end_subseq): helices[0].beginning]
            try:
                self.tmd_gravy = helpers.calculate_tmd_gravy(tmd_subseq)
            except ZeroDivisionError as exc:
                print(self.seq, tmd_subseq)
                raise exc
            self.beginning_charge = helpers.calculate_charge(beginning_subseq)
            self.end_charge = helpers.calculate_charge(end_subseq)

    def add_go(self, go_categories, go_descriptions):
        self.go_categories = go_categories
        self.go_descriptions = go_descriptions

    def add_blast_nr(self, ncbi_hit_description):
        self.ncbi_nr_hit = ncbi_hit_description

    def annotate(self):
        self.targeted_by = self.check_if_targeted()
        self.mbomp = self.check_if_mbomp()
        self.tail_anchored = self.check_if_ta()
        if self.targeted_by or \
            self.mbomp or \
            self.tail_anchored or \
            self.hmmer:
            self.interesting = True

    def check_if_targeted(self):
        cnt = 0
        if self.targetp.check():
            cnt += 1
        if self.mitofates.check():
            cnt += 1
        if self.nommpred:
            if self.nommpred.localisation == 'M':
                cnt += 1
        return cnt

    def check_if_mbomp(self):
        if self.helices:
            return False
        if self.targeted_by:
            return False
        if not self.beta_signal:
            return False
        else:
            return True

    def check_if_ta(self):
        if len(self.helices) == 1:
            if self.helices[0].beginning > self.length - 31:
                return True
        else:
            return False


class Targeting:
    def __init__(self, localisation, probability):
        self.localisation = localisation
        self.probability = probability

    def check(self):
        if self.localisation == 'M' and float(self.probability) >= 0.5:
            return True
        else:
            return False

    def __str__(self):
        return '({0},{1})'.format(self.localisation, self.probability)


class TransMembraneHelix:
    def __init__(self, beginning, end):
        if beginning < end:
            self.beginning = beginning
            self.end = end
            self.length = end - beginning
        else:
            raise ValueError("The beginning should be before the end.")

    def __str__(self):
        return '[{0},{1},{2}]'.format(self.beginning, self.end, self.length)


class BetaSignalHit:
    def __init__(self, name, seq, start, end):
        self.name = name
        self.seq = seq
        self.start = start
        self.end = end
        self.length = end - start

    def __str__(self):
        return [self.name, self.seq, self.subseq]


class Results:
    def __init__(self):
        self.dict_of_proteins = {}

    def add_hit(self, name_of_protein, name_of_hit):
        raise Exception("Not implemented")

    def keys(self):
        return self.dict_of_proteins.keys()

    def hits(self, name_of_protein):
        return self.dict_of_proteins[name_of_protein]


class HomologyResults(Results):
    def add_hit(self, name_of_protein, name_of_hit):
        if name_of_protein in self.dict_of_proteins.keys():
            self.dict_of_proteins[name_of_protein].append(name_of_hit)
        else:
            self.dict_of_proteins[name_of_protein] = [name_of_hit]


class SingularResults(Results):
    def add_hit(self, name_of_protein, hit):
        if name_of_protein in self.dict_of_proteins.keys():
            print(len(self.keys()))
            print(self.dict_of_proteins[name_of_protein])
            raise ValueError("Protein {} annotated twice.".format(name_of_protein))
        else:
            self.dict_of_proteins[name_of_protein] = hit
