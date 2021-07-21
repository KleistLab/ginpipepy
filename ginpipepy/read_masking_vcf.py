#read_masking_vcf.py


"""
Created on Mon Jul 13 15:15:41 2021

@author: Maria Trofimova
"""

import vcf

class VCFreader:

	def __init__(self, filename, reffile, fp_dict):
		self.fp_dict = fp_dict
		self.filename = filename
		self.reffile = reffile

	def _getMaskingFromVCF(self):
		"""
        	Gets a list of bases from aa Variant Calling File that have a FILTER tag
		'mask'
        	:return list: posiions of bases that should be masked
        	"""
		vcfile = vcf.Reader(open(self.filename),'r')
		masking_list = []
		for record in vcfile:
			operation = record.FILTER
			if 'mask' in operation:
				masking_list.append(int(record.POS))
		return masking_list

	def maskBasesInFP(self):
		"""
        	Removes mutated poositions from fingerprints based on Variant Calling
		File, where FILTER tag has value 'mask' and returns a filtered fingerprint
		list
        	:return list: set of filtered fingerprints in format (id_string,fp_string)
        	"""
		masking_list = self._getMaskingFromVCF()
		newSeqSets = []
		for t, seqSet in enumerate(self.fp_dict):
			if len(seqSet)!=0:
				newSeqSet = []
				for i, seq in enumerate(seqSet):
					if seq[1]!='':
						new_str = []
						split_fp = seq[1].split("-")
						for pos in split_fp:
							spl = pos.split(">")
							posit = int(spl[0])
							if not posit in masking_list:
								new_str.append(pos)
						new_str_j = '-'.join(new_str)
						newSeqSet.append((seq[0],new_str_j))
					else:
						newSeqSet.append(seq)
				newSeqSets.append(newSeqSet)
		return newSeqSets
