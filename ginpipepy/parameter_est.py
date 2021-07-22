#parameter_est.py

"""
Infer global parameters and filter sites for origins countings.

Created on Mon Nov 30 12:50:41 2020

@author: Maria Trofimova
"""


class parameterEstimation:
    """
    parameterEstimation class.

    Infer global parameters and filter sites for origins countings

    """

    def __init__(self, seqSets, reffile, freqCutoff):
        """
        Filter sequence fingerprints based on desired base frequency.

        :param seqSets: list of sequence fingerprints in format ('id','MUT_POS_1>MUT_BASE-MUT_POS_2>MUT_BASE-...')  
        :type seqSets: list
        :param reffile: path to reference FASTA file
        :type reffile: str
        :param freqCutoff: frequency cutoff (as count) for excluding mutants from fingerprints
        :type freqCutoff: int
        """
        self.seqSets = seqSets
        self.freqCutoff = freqCutoff
        self.reffile = reffile

    def _get_alt_pos_counts(self):
        """
        Get number of mutant bases in reference by sequence fingerprints.

        :returns positions: dictionary (position,number_of_mutants) 
        :rtype: dict
        """
        # Aggregate all mutant positions in a dictionary
        positions = dict()

        for t, seqSet in enumerate(self.seqSets):
            if len(seqSet)!=0:
                for i, seq in enumerate(seqSet):
                    if seq[1]!='':
                        split_fp = seq[1].split("-")
                        for pos in split_fp:
                            spl = pos.split(">")
                            pos = int(spl[0])
                            if pos in positions:
                                positions[pos] += 1
                            else:
                                positions[pos] = 1

        return positions

    def _count_mutants(self,seqSets):
        """
        Count mutant sequences in sample.

        :param seqSets: set of bins of sequence fingerprints
        :type seqSets: list
        :returns mut_count: list of mutant sequences counts per bin
        :rtype: int
        :returns num_seqs: list of number of sequences per bin
        :rtype: int
        """
        mut_count = 0
        num_seqs = 0
        # Enumerate bins
        for t, seqSet in enumerate(seqSets):
            if len(seqSet)!=0:
                # Enumerate sequences
                for i, seq in enumerate(seqSet):
                    #print(seq)
                    if seq[1]!='':
                        mut_count += 1
                        num_seqs += 1
                    else:
                        num_seqs += 1
        return mut_count, num_seqs

    def _filter_singletons(self):
        """
        Filter our mutant positions that only occur few times in entire sample.

        frequency cutoff defined by user via config
        :returns filteredSeqSets: sequence fingerprints without singletons
        :rtype: list
        """
        filteredSeqSets = []
        # First identify positions that pop up below a threshold

        positions = self._get_alt_pos_counts()

        # Identify not so frequent mutant positions
        minor_positions = []
        for (key,value) in positions.items():
            if value<=self.freqCutoff:
                #print("Singleton position: ",key)
                #print("Count: ",value)
                minor_positions.append(key)
        #print("Positions below cutoff:")

        # Filter sequence fingerprint based on that
        for t, seqSet in enumerate(self.seqSets):
            filteredSets = []
            for i, seq in enumerate(seqSet):
                if seq[1]!='':
                    split_fp = seq[1].split("-")
                    to_keep = []
                    for pos in split_fp:
                        split_pos = pos.split(">")
                        posit = int(split_pos[0])
                        if not posit in minor_positions:
                            to_keep.append(pos)
                    merger = ''
                    if len(to_keep)!=0:
                        merger = "-".join(to_keep)
                    filteredSets.append((seq[0],merger))
                else:
                    filteredSets.append((seq[0],seq[1]))
            filteredSeqSets.append(filteredSets)

        return filteredSeqSets


    def run(self):
        """
        Run filter - remove variants that are below predefined cutoff.

        :returns mut_count_post/num_seqs_post: mutant sequences proportion
        :rtype: float
        :returns filteredSets1: sequences set filtered based on mutant positions
            frequency cutoff
        :rtype: list
        """
        # In one of the trajectories - count mutants and number of sequences
        mut_count_pre, num_seqs_pre = self._count_mutants(self.seqSets)
        print("           Before filtering: ")
        print("              Number of sequences: ",num_seqs_pre)
        print("              Number of mutant sequences: ",mut_count_pre)

        filteredSets1 = self._filter_singletons()
        # Replace sequence set
        self.seqSets = filteredSets1

        mut_count_post, num_seqs_post = self._count_mutants(filteredSets1)
        print("           After filtering: ")
        print("              Number of sequences: ",num_seqs_post)
        print("              Number of mutant sequences: ",mut_count_post)

        return mut_count_post/num_seqs_post, filteredSets1
