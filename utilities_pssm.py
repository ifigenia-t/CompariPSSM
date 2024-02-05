import math

import utilities_stats

###################################################################
### PSSM
###################################################################

aas = list("CPQNTSGAVILMFYWHKRDE")

def read_pssm(pssm_file_path,delimiter="\t",read_range=[]):
	pssm_file_str = open(pssm_file_path).read()
	return process_pssm_file_string(pssm_file_str,delimiter,read_range)
	
def process_pssm_file_string(pssm_file_string,delimiter="\t",read_range=[]):
	lines = pssm_file_string.split("\n")
	
	pssm = {}
	for line in lines[1:]:
		if len(line.strip()) == 0: continue
		line_bits = line.split(delimiter)
		aa = line_bits[0]

		if len(read_range) == 0:
			pssm[aa] = [float(x) for x in line_bits[1:]]
		elif len(read_range) == 2:
			pssm[aa] = [float(x) for x in line_bits[1:]][read_range[0]:read_range[1]]
		else:
			pssm[aa] = [float(x) for x in line_bits[1:]]

	return pssm

def make_empty_pssm(length,aas=[]):
	if len(aas) == 0: aas = list("CPQNTSGAVILMFYWHKRDE")

	pssm = {}
	for aa in aas:
		pssm[aa] = []

	for i in range(0,length):	
		for aa in aas:
			pssm[aa].append(0)

	return pssm

def make_peptide_pssm(peptide):
	pssm =  make_empty_pssm(len(peptide))
	
	for i in range(0,len(peptide)):
		aa = peptide[i]
		pssm[aa][i] = 1

	return pssm

def print_pssm(pssm,precision=2,scientific=False,log=False,normalise=False,added_gini=False,aas_order=None,aas=[],print_range=[]):
	if len(aas) == 0: aas = list("CPQNTSGAVILMFYWHKRDE")
	if aas_order != None:
		aas = aas_order

	if normalise:
		pssm = normalise_pssm(pssm)

	if len(print_range) > 0:
		print_range = print_range
	else:
		print_range = range(0,len(pssm['A']))

	rows = []
	rows.append("\t".join(["aa"] + [str(x+1) for x in print_range]))
	for aa in aas:
		row = [aa]
		for i in print_range:
			try:
				if scientific:
					row.append("%1.3g"%pssm[aa][i])
				elif log:
					row.append(str(round(math.log10(pssm[aa][i]),precision)))
				else:
					row.append(str(round(pssm[aa][i],precision)))
			except:
				row.append(str(round(float(pssm[aa][i]),precision)))
				
		rows.append("\t".join(row))

	if added_gini:		
		gini = get_pssm_gini(pssm)
		row = ["gini"]
		for i in print_range:
			row.append("%1.3f"%gini[i])
		rows.append("\t".join(row))

	pssm_txt = "\n".join(rows)
	print(pssm_txt)
	return pssm_txt

def sliding_window(peptide,pssm):
	pssm_length = len(pssm['A'])
	peptide_length = len(peptide)
	peptides = []
	for seq_iter in range(-pssm_length  + 1, len(peptide)):
		if seq_iter < 0:
			padding_length = abs(seq_iter)
			tmp_peptide  = "-"*padding_length + peptide[max(0,seq_iter): max(0,seq_iter) + pssm_length + seq_iter]

		elif seq_iter > (peptide_length - pssm_length):
			padding_length = seq_iter - (peptide_length - pssm_length )
			tmp_peptide = peptide[seq_iter: seq_iter + pssm_length  - (seq_iter - (peptide_length- pssm_length ))] + "-"*padding_length
		else:
			tmp_peptide  = peptide[seq_iter: seq_iter + pssm_length ]

		if len(tmp_peptide) < pssm_length :
			tmp_peptide = tmp_peptide + "-"*(pssm_length  - len(tmp_peptide))

		peptides.append(tmp_peptide)

	return peptides

def score_peptide(peptide,pssm,sliding_window=False,gap_score=-3):
	pssm_length = len(pssm['A'])
	if sliding_window:
		return score_peptides(sliding_window(peptide,pssm),all=False)
	else:
		summer_scores = []
		for i in range(0,pssm_length):
			if peptide[i] in pssm:
				summer_scores.append(pssm[peptide[i]][i])
			elif peptide[i] == "-":
				summer_scores.append(gap_score)

		score = sum(summer_scores)
		return {"score":score,"scores":summer_scores}

def score_peptides(peptides,pssm,all=False):
	peptide_scores = {"best_hit":None,"best_score":None,"peptides":{}}

	for peptide in peptides:
		peptide_scores["peptides"][peptide] = score_peptide(peptide,pssm)

		if peptide_scores["best_hit"] == None or peptide_scores["peptides"][peptide]["score"] > peptide_scores["best_score"]:
			peptide_scores["best_hit"] = peptide
			peptide_scores["best_score"] = peptide_scores["peptides"][peptide]["score"]

	if all:
		return peptide_scores
	else:
		json = peptide_scores["peptides"][peptide_scores["best_hit"]]
		json['peptide'] = peptide_scores["best_hit"]

		return json


def range_values_pssm(pssm):
	min_value = None
	max_value = None
	for aa in aas:
		row = [aa]
		for i in range(0,len(pssm['A'])):
			if min_value == None or pssm[aa][i] < min_value:
				min_value = pssm[aa][i] 
			if max_value == None or pssm[aa][i] > max_value:
				max_value = pssm[aa][i] 
	
	return [min_value,max_value]

def normalise_pssm(pssm):
	return normalise_pssm_min_max(pssm)

def get_pssm_sparcity(pssm,no_data_value="-",aas=[]):
	if len(aas) == 0: aas = list("CPQNTSGAVILMFYWHKRDE")
	sparcity_counter = []
	
	for aa in aas:
		for i in range(0,len(pssm['A'])):
			sparcity_counter.append(pssm[aa][i] == no_data_value)

	return sparcity_counter.count(False)/len(sparcity_counter)

def get_pssm_sum(pssm,aas=[]):
	if len(aas) == 0: aas = list("CPQNTSGAVILMFYWHKRDE")
	summer = []
	
	for aa in aas:
		for i in range(0,len(pssm['A'])):
			summer.append(pssm[aa][i])

	return sum(summer)

def get_pssm_gini(pssm):
	pssm_columns = get_columns(pssm)
	gini_scores = {}
	for column in pssm_columns:
		gini_scores[column] = utilities_stats.gini_coefficient(pssm_columns[column])
	
	return gini_scores

def normalise_pssm_gini(pssm,aas=[]):
	if len(aas) == 0: aas = list("CPQNTSGAVILMFYWHKRDE")
	gini_scores = get_pssm_gini(pssm)
	normalised_pssm = {}

	for aa in aas:
		normalised_pssm[aa] = []
		for i in range(0,len(pssm['A'])):
			normalised_pssm[aa].append(pssm[aa][i]*(gini_scores[i]))

	return normalised_pssm

def pssm_normalise_column_sum_to_one(pssm,aas=[]):
	if len(aas) == 0: aas = list("CPQNTSGAVILMFYWHKRDE")
	normalised_pssm = normalise_pssm_min_max(pssm)
	column_sum_to_one_normalised_pssm = {}
	
	for aa in aas:
		column_sum_to_one_normalised_pssm[aa] = []

	for i in range(0,len(normalised_pssm['A'])):
		column_sum = 0
		for aa in aas:
			column_sum += normalised_pssm[aa][i]
	
		if column_sum == 0:
			for aa in aas:
				column_sum_to_one_normalised_pssm[aa].append(1/20)
		else:
			for aa in aas:
				column_sum_to_one_normalised_pssm[aa].append(normalised_pssm[aa][i]/column_sum)

	return column_sum_to_one_normalised_pssm

def pssm_normalise_column_sum_to_one_fold_change(pssm,logged=True,aas=[]):
	if len(aas) == 0: aas = list("CPQNTSGAVILMFYWHKRDE")
	
	normalised_pssm = normalise_pssm_min_max(pssm,aas=aas)
	column_sum_to_one_normalised_pssm = {}
	
	for aa in aas:
		column_sum_to_one_normalised_pssm[aa] = []

	for i in range(0,len(normalised_pssm['A'])):
		column_sum = 0
		for aa in aas:
			column_sum += normalised_pssm[aa][i]
	
		if column_sum == 0:
			for aa in aas:
				score = 1/20
				score = math.log2(score/0.05)
				column_sum_to_one_normalised_pssm[aa].append(score)
		else:
			for aa in aas:
				if normalised_pssm[aa][i] == 0:
					if logged:
						score = -math.log2(1/0.05)
					else:
						score = 0

					column_sum_to_one_normalised_pssm[aa].append(score)
				else:
					if logged:
						score = math.log2((normalised_pssm[aa][i]/column_sum)/0.05)
					else:
						score = ((normalised_pssm[aa][i]/column_sum))/0.05

					column_sum_to_one_normalised_pssm[aa].append(score)
	
	return column_sum_to_one_normalised_pssm

def extract_values_pssm_as_list(pssm,aas=[]):
	if len(aas) == 0: aas = list("CPQNTSGAVILMFYWHKRDE")

	values_list = []
	for aa in aas:
		for i in range(0,len(pssm['A'])):
			values_list.append(pssm[aa][i])

	return values_list

def normalise_subtraction(pssm,subtraction_value,aas=[]):
	if len(aas) == 0: aas = list("CPQNTSGAVILMFYWHKRDE")

	normalised_pssm = {}
	
	for aa in aas:
		normalised_pssm[aa] = []
		for i in range(0,len(pssm['A'])):
			normalised_pssm[aa].append(pssm[aa][i] - subtraction_value)

	return normalised_pssm

def normalise_pssm_min_max(pssm,aas=[]):
	if len(aas) == 0: aas = list("CPQNTSGAVILMFYWHKRDE")

	normalised_pssm = {}
	[min_value,max_value] = range_values_pssm(pssm)
	
	for aa in aas:
		normalised_pssm[aa] = []
		for i in range(0,len(pssm['A'])):
			if max_value-min_value != 0:
				normalised_pssm[aa].append((pssm[aa][i] - min_value)/(max_value-min_value))
			else:
				normalised_pssm[aa].append(0)

	return normalised_pssm

def range_values_pssm_column(pssm,aas=[]):
	range_values = {}
	for i in range(0,len(pssm['A'])):
		range_values[i] = {"column":[]}
		for aa in aas:
			range_values[i]["column"].append(pssm[aa][i])
		
		range_values[i]["min_value"] = min(range_values[i]["column"])
		range_values[i]["max_value"] = max(range_values[i]["column"])

	return range_values

def normalise_pssm_column_min_max(pssm,aas=[]):
	if len(aas) == 0: aas = list("CPQNTSGAVILMFYWHKRDE")

	normalised_pssm = {}
	range_values = range_values_pssm_column(pssm,aas)
	
	for aa in aas:
		normalised_pssm[aa] = []
		for i in range(0,len(pssm['A'])):
			if range_values[i]["max_value"] - range_values[i]["min_value"] != 0:
				normalised_pssm[aa].append((pssm[aa][i] - range_values[i]["min_value"] )/(range_values[i]["max_value"] -range_values[i]["min_value"] ))
			else:
				normalised_pssm[aa].append(0)

	return normalised_pssm

def positivise_pssm(pssm,aas=[]):
	if len(aas) == 0: aas = list("CPQNTSGAVILMFYWHKRDE")
	positivised_pssm = {}
	
	for aa in aas:
		positivised_pssm[aa] = []
		for i in range(0,len(pssm['A'])):
			if pssm[aa][i] > 0:
				positivised_pssm[aa].append(pssm[aa][i])
			else:
				positivised_pssm[aa].append(0)

	return positivised_pssm

def cut_pssm(pssm,start,end,aas=[]):
	if len(aas) == 0: aas = list("CPQNTSGAVILMFYWHKRDE")
	cut_pssm = {}
	for aa in aas:
		cut_pssm[aa] = []
		for i in range(start,min(end,len(pssm[aa]))):
			if i <0:
				cut_pssm[aa].append(0.0)
			else:
				cut_pssm[aa].append(pssm[aa][i])

	return cut_pssm

def pearson_correlation(column_a, column_b):
	try:
		products_mean = sum([column_a[i] * column_b[i] for i in range(0, len(column_a))])/len(column_a)
		covariance = products_mean - (sum(column_a)/len(column_a) * sum(column_b)/len(column_b))
		column_a_standard_deviation =  standard_deviation(column_a)
		column_b_standard_deviation = standard_deviation(column_b)
		pearson_correlation = covariance / (column_a_standard_deviation * column_b_standard_deviation)
		return pearson_correlation
	except:
		#logger.error([column_a, column_b])
		return 0

def standard_deviation(column):
	column_squared = [i**2 for i in column]
	squares_mean = sum(column_squared)/len(column_squared)
	column_mean = sum(column)/len(column)
	column_mean_squared = column_mean**2
	variance = squares_mean - column_mean_squared
	std_dev = math.sqrt(variance)
	return std_dev 


def compare_pssms(pssm_a,pssm_b,aas=[]):	
	if len(pssm_a['A']) != len(pssm_b['A']):
		return {"status":"Error","error_type":"PSSMs differ in size"}

	pssm_a_columns = get_columns(pssm_a)
	pssm_b_columns = get_columns(pssm_b)
	
	mae_list = []
	for i in range(0,len(pssm_a_columns)):
		column_a = pssm_a_columns[i]
		column_b = pssm_b_columns[i]
		
		mae = pearson_correlation(column_a,column_b) #sum([abs(column_a[i] - column_b[i]) for i in range(0,len(column_a))])/(sum(column_a)+sum(column_b)) #Mean absolute error
		mae_list.append(mae)

	return mae_list

def substract_pssms(pssm_a,pssm_b,aas=[]):
	if len(aas) == 0: aas = list("CPQNTSGAVILMFYWHKRDE")
	substracted_pssm = {}

	pssm_a_column_sums = get_column_sums(pssm_a)
	pssm_b_column_sums = get_column_sums(pssm_b)

	if len(pssm_a['A']) != len(pssm_b['A']):
		return {"status":"Error","error_type":"PSSMs differ in size"}

	for aa in aas:
		substracted_pssm[aa] = []
		for i in range(0,len(pssm_a['A'])):
			substracted_pssm[aa].append(pssm_a[aa][i]/pssm_a_column_sums[i] - pssm_b[aa][i]/pssm_b_column_sums[i])

	return substracted_pssm



def get_columns(pssm,aas=[]):
	if len(aas) == 0: aas = list("CPQNTSGAVILMFYWHKRDE")
	pssm_columns = {}

	for i in range(0,len(pssm['A'])):
		pssm_columns[i] = []
		for aa in aas:
			pssm_columns[i].append(pssm[aa][i])
			
	return pssm_columns

def get_column_sums(pssm,aas=[]):
	if len(aas) == 0: aas = list("CPQNTSGAVILMFYWHKRDE")
	pssm_column_sums = {}

	for i in range(0,len(pssm['A'])):
		pssm_column_sums[i] = 0 
		for aa in aas:
			pssm_column_sums[i] += pssm[aa][i]
			
	return pssm_column_sums

def get_pssm_motif(pssm,cut_off=0.5,aas=[]):
	if len(aas) == 0: aas = list("CPQNTSGAVILMFYWHKRDE")
	normalised_pssm = positivise_pssm(pssm)
	normalised_pssm = pssm_normalise_column_sum_to_one(normalised_pssm)

	motif = []
	for i in range(0,len(normalised_pssm['A'])):
		motif_aas = []
		for aa in aas:
			if normalised_pssm[aa][i] >= cut_off:
				motif_aas.append(aa)

		if len(motif_aas) == 0:
			motif.append(".")
		elif len(motif_aas) == 1:
			motif.append(motif_aas[0])
		else:
			motif.append("[" + "".join(motif_aas) + "]")

	return motif


def check_pssm(pssm):
	aaStandard = ['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y']

	aa_pssms = pssm.keys()
	if set(aa_pssms).symmetric_difference(set(aaStandard)):
		return {"status":"error", "error_type":"Incorrect keys, should be standard amino acids: " + ",".join(aaStandard)}
	
	for aa in aaStandard:
		for idx,val in enumerate(pssm[aa]):
			try:
				float(val)
			except:
				return {"status":"error", "error_type":"not float value for PSSM["+str(aa)+"]["+str(idx)+"]"}
				
	pssm_len = list(set([len(x) for x in list(pssm.values())]))
	if len(pssm_len) > 1:
		return {"status":"error", "error_type":"PSSM - incorrect number of columns for aas. All should be of equal length. Range:" + ",".join(map(str,pssm_len)) }
		
	pssm_len = pssm_len[0]
	if pssm_len < 3:
		return {"status":"error", "error_type":"PSSM is too short, PSSM length should be > 2, not " + str(pssm_len) + ". Key: pssm."}
		

	if pssm_len > 40:
		return {"status":"error", "error_type":"PSSM is too long, PSSM length should be < 40, not " + str(pssm_len) + ". Key: pssm."}
		
	return {"status":"success"}
		