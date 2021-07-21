#test_read_masking_vcf.py


"""
Created on Mon Jul 13 15:15:41 2021

@author: Maria Trofimova
"""

import vcf

class TestVCFreader:

	def __init__(self, filename):
		self.filename = filename

	def _getMaskingFromVCF(self):
		vcfile = vcf.Reader(open(self.filename),'r')
		masking_list = []
		for record in vcfile:
			operation = record.FILTER
			if 'mask' in operation:
				masking_list.append(int(record.POS))
		return masking_list

	def testVCFReader(self):
		masking_list = self._getMaskingFromVCF()
		before = "8>A-12>G-13>A-17>C-19>G-20>G-25>C-33>C-45>A-46>G-51>A-53>C-56>A-70>C-72>C-78>A-80>G-85>G-90>G-92>A-93>C-108>T-115>G-124>A-151>T-157>A-161>T-164>A-168>G-170>G-171>C-191>G-194>A"
		after = "8>A-12>G-13>A-17>C-19>G-20>G-25>C-33>C-45>A-46>G-51>A-53>C-56>A-72>C-78>A-80>G-85>G-90>G-92>A-93>C-115>G-124>A-151>T-157>A-161>T-164>A-168>G-170>G-171>C-191>G-194>A"
		new_str = []
		split_fp = before.split("-")
		for pos in split_fp:
			spl = pos.split(">")
			posit = int(spl[0])
			if not posit in masking_list:
				new_str.append(pos)
		new_str_j = '-'.join(new_str)
		if new_str_j==after:
			print("Removed masked bases successfully.")
		else:
			print("Could not remove masked bases correctly.")
