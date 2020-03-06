### Santalucia_NN_Tm; a part of Thermo-Align tool for the design of template specific hybridization and priming oligonucleotides
### Updated: version-1.05 06/29/2017
### Property of Wisser Lab at University of Delaware
### Author: Felix Francis (felixfrancier@gmail.com)
#    Copyright (C) 2016-2030 by
#    Felix Francis <felixfrancier@gmail.com>
#    All rights reserved.

############################################################
#### IMPORT FUNCTIONS
############################################################
import math

############################################################
#### FUNCTIONS AND TABLES
############################################################

# SantaLucia & Hicks (2004), Annu. Rev. Biophys. Biomol. Struct 33: 415-440
# delta H (Enthalpy)(kcal/mol) and delta S (Entropy)(eu) coefficients
DNA_NN_table = { 
	'init': (0.2, -5.7), 'init_A/T': (2.2, 6.9), 'init_G/C': (0, 0), 
	'init_oneG/C': (0, 0), 'init_allA/T': (0, 0), 'init_5T/A': (0, 0),
	'sym': (0, -1.4),
	'AA/TT': (-7.6, -21.3), 'AT/TA': (-7.2, -20.4), 'TA/AT': (-7.2, -20.4), 
	'CA/GT': (-8.5, -22.7), 'GT/CA': (-8.4, -22.4), 'CT/GA': (-7.8, -21.0), 
	'GA/CT': (-8.2, -22.2), 'CG/GC': (-10.6, -27.2), 'GC/CG': (-9.8, -24.4), 
	'GG/CC': (-8.0, -19.0)}

# Internal mismatch and inosine table (DNA) 
# Allawi & SantaLucia (1997), Biochemistry 36: 10581-10594 
# Allawi & SantaLucia (1998), Biochemistry 37: 9435-9444 
# Allawi & SantaLucia (1998), Biochemistry 37: 2170-2179 
# Allawi & SantaLucia (1998), Nucl Acids Res 26: 2694-2701 
# Peyret et al. (1999), Biochemistry 38: 3468-3477 					

DNA_IMM_table = { 
	'AG/TT': (1.0, 0.9), 'AT/TG': (-2.5, -8.3), 'CG/GT': (-4.1, -11.7), 				# Allawi & SantaLucia (1997), Biochemistry 36: 10581-10594  http://pubs.acs.org/doi/pdf/10.1021/bi962590c
	'CT/GG': (-2.8, -8.0), 'GG/CT': (3.3, 10.4), 'GG/TT': (5.8, 16.3), 
	'GT/CG': (-4.4, -12.3), 'GT/TG': (4.1, 9.5), 'TG/AT': (-0.1, -1.7), 
	'TG/GT': (-1.4, -6.2), 'TT/AG': (-1.3, -5.3),
		
	'AA/TG': (-0.6, -2.3), 
	'AG/TA': (-0.7, -2.3), 'CA/GG': (-0.7, -2.3), 'CG/GA': (-4.0, -13.2), 
	'GA/CG': (-0.6, -1.0), 'GG/CA': (0.5, 3.2), 'TA/AG': (0.7, 0.7), 
	'TG/AA': (3.0, 7.4),
		
	'AC/TT': (0.7, 0.2), 'AT/TC': (-1.2, -6.2), 'CC/GT': (-0.8, -4.5), 					# Allawi & SantaLucia (1998), Nucl Acids Res 26: 2694-2701 
	'CT/GC': (-1.5, -6.1), 'GC/CT': (2.3, 5.4), 'GT/CC': (5.2, 13.5), 
	'TC/AT': (1.2, 0.7), 'TT/AC': (1.0, 0.7),
		
	'AA/TC': (2.3, 4.6), 'AC/TA': (5.3, 14.6), 'CA/GC': (1.9, 3.7), 					# Allawi & SantaLucia (1998), Biochemistry 37: 9435-9444  http://pubs.acs.org/doi/pdf/10.1021/bi9803729 
	'CC/GA': (0.6, -0.6), 'GA/CC': (5.2, 14.2), 'GC/CA': (-0.7, -3.8), 
	'TA/AC': (3.4, 8.0), 'TC/AA': (7.6, 20.2),
		
	'AA/TA': (1.2, 1.7), 'CA/GA': (-0.9, -4.2), 'GA/CA': (-2.9, -9.8), 					# Peyret et al. (1999), Biochemistry 38: 3468-3477 
	'TA/AA': (4.7, 12.9), 'AC/TC': (0.0, -4.4), 'CC/GC': (-1.5, -7.2), 
	'GC/CC': (3.6, 8.9), 'TC/AC': (6.1, 16.4), 'AG/TG': (-3.1, -9.5), 
	'CG/GG': (-4.9, -15.3), 'GG/CG': (-6.0, -15.8), 'TG/AG': (1.6, 3.6), 
	'AT/TT': (-2.7, -10.8), 'CT/GT': (-5.0, -15.8), 'GT/CT': (-2.2, -8.4), 
	'TT/AT': (0.2, -1.5)} 

# Terminal mismatch table (DNA) 
# SantaLucia & Peyret (2001) Patent Application WO 01/94611 
DNA_TMM_table = { 
	'AA/TA': (-3.1, -7.8),	'TA/AA': (-2.5, -6.3),
	'CA/GA': (-4.3, -10.7),	'GA/CA': (-8.0, -22.5), 
	'AC/TC': (-0.1, 0.5),	'TC/AC': (-0.7, -1.3),
	'CC/GC': (-2.1, -5.1),	'GC/CC': (-3.9, -10.6), 
	'AG/TG': (-1.1, -2.1),	'TG/AG': (-1.1, -2.7),
	'CG/GG': (-3.8, -9.5),	'GG/CG': (-0.7, -19.2), 
	'AT/TT': (-2.4, -6.5), 'TT/AT': (-3.2, -8.9),
	'CT/GT': (-6.1, -16.9),	'GT/CT': (-7.4, -21.2), 
	'AA/TC': (-1.6, -4.0), 'AC/TA': (-1.8, -3.8), 'CA/GC': (-2.6, -5.9), 
	'CC/GA': (-2.7, -6.0), 'GA/CC': (-5.0, -13.8), 'GC/CA': (-3.2, -7.1), 
	'TA/AC': (-2.3, -5.9), 'TC/AA': (-2.7, -7.0), 
	'AC/TT': (-0.9, -1.7), 'AT/TC': (-2.3, -6.3), 'CC/GT': (-3.2, -8.0), 
	'CT/GC': (-3.9, -10.6), 'GC/CT': (-4.9, -13.5), 'GT/CC': (-3.0, -7.8), 
	'TC/AT': (-2.5, -6.3), 'TT/AC': (-0.7, -1.2), 
	'AA/TG': (-1.9, -4.4), 'AG/TA': (-2.5, -5.9), 'CA/GG': (-3.9, -9.6), 
	'CG/GA': (-6.0, -15.5), 'GA/CG': (-4.3, -11.1), 'GG/CA': (-4.6, -11.4), 
	'TA/AG': (-2.0, -4.7), 'TG/AA': (-2.4, -5.8), 
	'AG/TT': (-3.2, -8.7), 'AT/TG': (-3.5, -9.4), 'CG/GT': (-3.8, -9.0), 
	'CT/GG': (-6.6, -18.7), 'GG/CT': (-5.7, -15.9), 'GT/CG': (-5.9, -16.1), 
	'TG/AT': (-3.9, -10.5), 'TT/AG': (-3.6, -9.8)} 

# Dangling ends table (DNA) 
# Bommarito et al. (2000), Nucl Acids Res 28: 1929-1934 
DNA_DE_table = { 
	'AA/.T': (0.2, 2.3), 'AC/.G': (-6.3, -17.1), 'AG/.C': (-3.7, -10.0), 
	'AT/.A': (-2.9, -7.6), 'CA/.T': (0.6, 3.3), 'CC/.G': (-4.4, -12.6), 
	'CG/.C': (-4.0, -11.9), 'CT/.A': (-4.1, -13.0), 'GA/.T': (-1.1, -1.6), 
	'GC/.G': (-5.1, -14.0), 'GG/.C': (-3.9, -10.9), 'GT/.A': (-4.2, -15.0), 
	'TA/.T': (-6.9, -20.0), 'TC/.G': (-4.0, -10.9), 'TG/.C': (-4.9, -13.8), 
	'TT/.A': (-0.2, -0.5), 
	'.A/AT': (-0.7, -0.8), '.C/AG': (-2.1, -3.9), '.G/AC': (-5.9, -16.5), 
	'.T/AA': (-0.5, -1.1), '.A/CT': (4.4, 14.9), '.C/CG': (-0.2, -0.1), 
	'.G/CC': (-2.6, -7.4), '.T/CA': (4.7, 14.2), '.A/GT': (-1.6, -3.6), 
	'.C/GG': (-3.9, -11.2), '.G/GC': (-3.2, -10.4), '.T/GA': (-4.1, -13.1), 
	'.A/TT': (2.9, 10.4), '.C/TG': (-4.4, -13.1), '.G/TC': (-5.2, -15.0), 
	'.T/TA': (-3.8, -12.6)} 

### DELTA S ION CORRECTION FUNCTION
def ion_correction(Na=0, K=0, Tris=0, Mg=0, dNTPs=0, seq_len=None):
	'''
	Correction for deltaS(Entropy correction) : 0.368 x (N-1) x ln[Na+]
	Reference: (SantaLucia (1998), Proc Natl Acad Sci USA 95: 1460-1465)	http://www.pnas.org/content/95/4/1460.full.pdf+html
	Provide Millimolar(mmol/L) concentrations for Na, K, Tris, Mg and dNTPs
	Von Ahsen et al. (2001, Clin Chem 47: 1956-1961)	http://www.clinchem.org/content/47/11/1956.full
	[Na_eq] = [Na+] + [K+] + [Tris]/2 + 120*([Mg2+] - [dNTPs])^0.5 
	If [dNTPs] >= [Mg2+]: [Na_eq] = [Na+] + [K+] + [Tris]/2
	effect on entropy by salt correction; von Ahsen et al 1999
		s+=0.368 * (strlen(c)-1)* log(salt_effect)
	'''
	correction_factor		=	0
	Na=float(Na)
	K=float(K)
	Tris=float(Tris)
	Mg=float(Mg)
	dNTPs=float(dNTPs)
	Monovalent_cations_mmol	=	Na + K + (Tris/2)
	if dNTPs >= Mg:
		Na_eq_mmol		=	Monovalent_cations_mmol 
	else:
		Na_eq_mmol		=	Monovalent_cations_mmol + (120 * math.sqrt(Mg-dNTPs))
	Na_eq_mol	=	Na_eq_mmol/1000.0												# Convert from millimolar to molar
	correction_factor		=	0.368 * (seq_len - 1) * math.log(Na_eq_mol)
	return correction_factor

### FUNCTION TO CHECK WHETHER TWO SEQUENCES ARE REVERSE COMPLEMENTS OF EACH OTHER
def transform_degenerate(seq):
	seq = seq.upper()
	basesub = {'R': 'A',
		'Y': 'C',
		'M': 'A',
		'K': 'G',
		'S': 'G',
		'W': 'A',
		'H': 'A',
		'B': 'G',
		'V': 'G',
		'D': 'G'}
	for base in basesub.keys():
		seq = seq.replace(base, basesub[base])
	return seq

def rev_complement(seq):
	seq = seq.upper()
	basecomplement = {'A':'T', 
					  'C':'G', 
					  'G':'C', 
					  'T':'A', 
					  '-':'-', 
					  'N':'N'}
	letters = list(seq)
	newletters = []
	for base in letters:
		if base in basecomplement:
			newletters.append(basecomplement[base])
		else:
			newletters.append('N')
	complement = (''.join(newletters))                                              #gives the complement of the bases in list letters
	return complement[::-1]                                                         #gives the reverse of the compliment   

def complement(seq):
	seq = seq.upper()
	basecomplement = {'A':'T', 
					  'C':'G', 
					  'G':'C', 
					  'T':'A', 
					  '-':'-', 
					  'N':'N'}
	letters = list(seq)
	newletters = []
	for base in letters:
		if base in basecomplement:
			newletters.append(basecomplement[base])
		else:
			newletters.append('N')
	complement = (''.join(newletters))                                                 #gives the complement of the bases in list letters
	return complement 

### FUNCTION TO QC INPUT OLIGO SEQUENCES. (covert the input seq to uppercase, remove any whitespace and characters other than 'A','C','G','T','N','-')
def seq_qc(seq):
	seq = ''.join(seq.split()).upper()
	allowed_chars =('A','C','G','T','N','-')													# Reconsider including other allowed characters '-"???
	seq = ''.join([base for base in seq if base in allowed_chars])
	return seq

### FUNCTION TO CALCULATE millimolar (mM) CONCENTRATION OF MONOVALENT AND DIVALENT CATIONS (OUTPUT IS IN mM)
def mM_monovalent(Na=0, K=0, Tris=0, Mg=0, dNTPs=0):
	Na=float(Na)
	K=float(K)
	Tris=float(Tris)
	Mg=float(Mg)
	dNTPs=float(dNTPs)
	Monovalent_cations_mmol	=	Na + K + (Tris/2)
	if dNTPs >= Mg:
		Na_eq_mmol		=	Monovalent_cations_mmol 
	else:
		Na_eq_mmol		=	Monovalent_cations_mmol + (120 * math.sqrt(Mg-dNTPs))
	return Na_eq_mmol

### FUNCTION TO CALCULATE NN TERMODYNAMICS BASED MELTING TEMPERATURE(TM)

def NN_Tm(seq=None, compl_seq=None, primer_conc=100, Na=0, K=50, Tris=10, Mg=1.5, dNTPs=0.2, ion_corr=True):
	""" Calculates the Tm using the nearest neighbor thermodynamics
	Arguments:
		-seq		: The primer sequence 			(5'->3' direction)
		-compl_seq	: The complementary sequence	(3'->5' direction!!!!!) 
		-primer_conc	: Concentration of the primer [nM]. Template strand which concentration is typically very low and may be ignored and so not included in this function.
		-Na, K, Tris, Mg, dNTPs	: Concentration of the respective ions [mM]. If any of K, Tris, Mg and dNTPS is non-zero, a 'sodium-equivalent' concentration
								  is calculated and used for salt correction (von Ahsen et al., 2001).
		-ion_corr	: See method 'Tm_GC'. Default = True. (0 means no salt correction).	
	"""
	if not seq:
		raise ValueError('Please provide an input sequence!')
	seq = seq_qc(seq)
	if not compl_seq:											# If no complementary seq is provided, the complement of the given seq is caclulated and used
			compl_seq = complement(seq)
	compl_seq = seq_qc(compl_seq)
	
	temp_seq = str(seq)
	temp_compl_seq = str(compl_seq)
	dH			= 0											# dH stands for delta H
	dS			= 0											# dS stands for delta S
	dH_index	= 0											# dH_index stands for deta H_index in a table (here position 0)
	dS_index	= 1											# dS_index stands for delta S_index in a table (here position 1)
	left_tmm	= temp_compl_seq[:2][::-1] + '/' + temp_seq[:2][::-1]
	if left_tmm in DNA_TMM_table:
		dH += DNA_TMM_table[left_tmm][dH_index]
		dS += DNA_TMM_table[left_tmm][dS_index]
		temp_compl_seq	= temp_compl_seq[1:]
		temp_seq		= temp_seq[1:]
	# right_tmm	= temp_compl_seq[-2:] + '/' + temp_seq[-2:]
	right_tmm	= temp_seq[-2:] + '/' + temp_compl_seq[-2:] 
	if right_tmm in DNA_TMM_table:
		dH += DNA_TMM_table[right_tmm][dH_index]
		dS += DNA_TMM_table[right_tmm][dS_index]	
		temp_compl_seq	= temp_compl_seq[:-1]	
		temp_seq		= temp_seq[:-1]
	# General initiation value
	dH += DNA_NN_table['init'][dH_index]
	dS += DNA_NN_table['init'][dS_index]
	# delta H and delta S correction for  Duplex with no (allA/T) or at least one (oneG/C) GC pair coefficients need not be considered while using SantaLucia & Hicks (2004) model
	# delta H and delta S correction for 5' end being T need not be considered while using SantaLucia & Hicks (2004) model 	
	# delta H and delta S correction for A/T terminal basepairs (for the original seq and not for the end trimmed temp_seq)
	### CONSIDERS TERMINAL MISMATCHES WHILE COUNTING TERMINAL A/Ts 
	terminals	=	seq[0]+'/'+compl_seq[0] + seq[-1] +'/'+compl_seq[-1]
	count_AT	=	terminals.count('A/T')+terminals.count('T/A')
	dH += DNA_NN_table['init_A/T'][dH_index] * count_AT
	dS += DNA_NN_table['init_A/T'][dS_index] * count_AT	
	# delta H and delta S calculations for matches and mimatches
	for base_index in range(len(temp_seq)-1):
		NN	=	temp_seq[base_index:base_index+2] + '/' + temp_compl_seq[base_index:base_index+2]
		if NN in DNA_IMM_table:
			dH += DNA_IMM_table[NN][dH_index]
			dS += DNA_IMM_table[NN][dS_index]
		elif NN[::-1] in DNA_IMM_table:
			dH += DNA_IMM_table[NN[::-1]][dH_index]
			dS += DNA_IMM_table[NN[::-1]][dS_index]
		elif NN in DNA_NN_table:
			dH += DNA_NN_table[NN][dH_index]
			dS += DNA_NN_table[NN][dS_index]
		elif NN [::-1] in DNA_NN_table:
			dH += DNA_NN_table[NN[::-1]][dH_index]
			dS += DNA_NN_table[NN[::-1]][dS_index]		
	if ion_corr:
		seq_len	=	len(seq)
		correction_factor = ion_correction(Na, K, Tris, Mg, dNTPs, seq_len)
		dS 	+=	correction_factor
	x	=	4															#	x = 4 if not self complementary; x = 1 if self complementary
	R	=	1.9872  													#	Universal gas constant in Cal/degrees C*Mol
	primer_conc = primer_conc/1000000000.0								#	To convert nM into molar; do not multiply with 2 if the other strand is genomic DNA template, which can be negligible = 
	Tm	= 	(1000* dH) / (dS + (R * (math.log(primer_conc/x)))) - 273.15 
	return ("{0:.2f}".format(round(Tm,2)))								# Round to 2 decimals


if __name__ == "__main__":
    print('CTTCTGCAATGCCAAGTCCAG')
    print('full match: ', NN_Tm(seq='CTTCTGCAATGCCAAGTCCAG', compl_seq='GAAGACGTTACGGTTCAGGTC'))
    print('One mismatch', NN_Tm(seq='CTTCTGCAATGCCAAGTCCAG', compl_seq='GAAGACGTCACGGTTCAGGTC'))


### Updated: 07/12/2019
# rename "xrange" to "range" in order to fit python3
# Add __main__ to debug
# setting default parameters to: primer_conc=100, Na=0, K=50, Tris=10, Mg=1.5, dNTPs=0.2, ion_corr=True in order to debug

### Updated: version-1.05 06/29/2017
# "right_tmm" definition modified
