#test_read_masking_vcf.py


"""
Created on Mon Jul 13 15:15:41 2021

@author: Maria Trofimova
"""

import vcf
from ginpipepy.read_masking_vcf import VCFreader

class TestVCFreader:

	def __init__(self, filename, reffile):
		self.filename = filename
		self.fp_dict = [("ID:1","8>A-12>G-13>A-17>C-19>G-20>G-25>C-33>C-45>A-46>G-51>A-53>C-56>A-70>C-72>C-78>A-80>G-85>G-90>G-92>A-93>C-108>T-115>G-124>A-151>T-157>A-161>T-164>A-168>G-170>G-171>C-191>G-194>A")]
		self.reffile = reffile

	def testVCFReader(self):
		mask = VCFreader(self.filename, self.reffile, self.fp_dict)
        filtered_seqset = mask.maskBasesInFP()
		after_masking = "8>A-12>G-13>A-17>C-19>G-20>G-25>C-33>C-45>A-46>G-51>A-53>C-56>A-72>C-78>A-80>G-85>G-90>G-92>A-93>C-115>G-124>A-151>T-157>A-161>T-164>A-168>G-170>G-171>C-191>G-194>A"

		if filtered_seqset[0][1]==after:
			print("Removed masked bases successfully.")
		else:
			print("Could not remove masked bases correctly.")
