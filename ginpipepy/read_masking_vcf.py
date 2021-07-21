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
		vcfile = vcf.Reader(open(self.filename),'r')
		print(vcfile)
		masking_list = []
		for record in vcfile:
			operation = record.INFO
			if str(operation)=='mask':
				masking_list.append(int(record.POS))
		return masking_list

	def maskBasesInFP(self):
		masking_list = self._getMaskingFromVCF()
		print(masking_list)
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
							else:
								print("Found masking position:")
								print(posit)
						new_str_j = '-'.join(new_str)
						print("Old CIGAR:")
						print(seq[1])
						print("New CIGAR:")
						print(new_str_j)
						print("\n")
						newSeqSet.append((seq[0],new_str_j))
					else:
						newSeqSet.append(seq)
				newSeqSets.append(newSeqSet)
		return newSeqSets
