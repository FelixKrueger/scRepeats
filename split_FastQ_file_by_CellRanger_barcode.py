#!/usr/bin/env python
import sys, gzip, os
import getopt
from glob import glob
from time import sleep
import re

# This script splits out FastQ files by the (corrected) cell-level barcode (from Seurat)
# Script last modified 27 April 2021, Felix Krueger

fhs = {}           # storing the filehandles for all output files

def submain():
	
	print (f"Python version: {sys.version}.")
	allfiles = glob("*.fastq.gz")
	allfiles.sort() # required as glob doesn't necessarily store files in alphabetical order


	for filename in allfiles:
		print (f"Reading in FastQ file:\t >> {filename} <<\n")
		main(filename)
		fhs.clear() # resetting filehandles
		
def main(filename):

	expected_count   = 0
	unexpected_count = 0   
	count = 0              # total sequence count
	
	unique_count = 0
	dups_count = 0
	
	observed = {}  # storing all cell-level barcodes

	print (f"Reading file: >{filename}<")
	
	with gzip.open(filename) as cf:
	
		while True:
			readID  = cf.readline().decode().strip()
			seq     = cf.readline().decode().strip()
			line3   = cf.readline().decode().strip()
			qual    = cf.readline().decode().strip()
			
			if not qual:
				break
			
			count += 1

			if count%500000 == 0:
				print (f"Processed {count} reads so far")

			# print (f"{readID}\n{seq}\n{line3}\n{qual}\n")
		
			elments   = readID.split(":")
			cell_barcode = elments[-1]
			# print (cell_barcode)
			
			if not cell_barcode in observed.keys():
				observed[cell_barcode] = 0

			observed[cell_barcode] += 1	

			if cell_barcode in fhs.keys():
				pass
			else:
				# need to make new filehandle
				make_out_filehandle(cell_barcode,filename)

			
			# readID = readID.replace(" ","_") # this is required for e.g. Bowtie2 to retain the last part of the read ID (= the UMI sequence)
			fhs[cell_barcode].write (("\n".join([readID, seq, line3, qual]) + "\n").encode())

	close_filehandles()

	print (f"Sequences processed: {count}\nUnique sequences: {unique_count}\nDuplicate sequences (removed): {dups_count}\n")

	
def make_out_filehandle(sample_name,filename):
	
	print (f"Making new filehandle for sample name: {sample_name}, and file {filename}")
	# extracting useful parts from filename

	pattern = '(.*)\.(fastq.gz)'
	p = re.compile(pattern)
	# print (filename)
	m = p.findall(filename)
	sample = m[0][0]
	ending = m[0][1	]
	new_filename = f"{sample}_{sample_name}.{ending}"
	# print (new_filename)
	
	fhs[sample_name] = gzip.open (new_filename,mode='wb',compresslevel=3)
	return

def close_filehandles():
	for name in fhs.keys():
		fhs[name].close()

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)	


if __name__ == "__main__":
	submain()
else:
	print ("Just getting imported")
