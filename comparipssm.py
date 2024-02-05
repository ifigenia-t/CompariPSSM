import random,pprint,json,sys,os

#-------#-------#-------#-------#-------#-------#-------#-------#-------#-------#-------

import logging

from bisect import bisect
import functools
import operator

import utilities_options
import utilities_stats
import utilities_pssm
import utilities_error

logger = logging.getLogger("compariPSSM")

##---------------------------------------------------------------------##

class CompariPSSM:
	"""
	Pipeline to process the peptides and make oligonucleotides.
	Attributes:
	"""

	def __init__(self):
		self.options = {
			"verbose":False,
			"debug":False,
			"query_pssm":{},
			"compare_pssm":{},
			"query_pssm_file":"",
			"compare_pssm_file":"",
			"output_file":"",
			'ignore_negative_pssm_values':False, 
			'significance_cutoff':0.0001,
			'specificity_scoring_type':"gini_normalised",
			"accession":False,
			"column_scores_samples":1000000,
			"compare_accession":False,
			"compare_pssm_score_scheme":"frequency",
			"compare_region_end":False,
			"compare_region_start":False,
			"compare_regions":False,
			"comparison_dataset":"momap_alignment_pssm",
			"detailed_results":False,
			"include_pssms":False,
			"include_query_pssm":False,
			"min_overlap":0,
			"mask_positions":[],
			"mask_score":0,
			"dissimilarity_score_focus":"query",
			"pssm_normalise_count_sum_to_one":False,
			"pssm_normalise_min_max":True,
			"query_accession":False,
			"query_dataset":"momap_alignment_pssm",
			"query_pssm_score_scheme":"frequency",
			"query_region_end":False,
			"query_region_start":False,
			"query_specificity_acc":False,
			"random_pssm_length":10,
			"random_pssm_samples":0,
			"region_end":False,
			"region_start":False,
			"significant_positions_count_cutoff":2,
			"significant_positions_cutoff":0.05,
			"sliding_window_comparison":True
		}	

		self.options.update(utilities_options.load_commandline_options(self.options,self.options,purge_commandline=False))

		if self.options['verbose']:
			logging.basicConfig()
			logger.setLevel(logging.INFO)

		if self.options['debug']:
			logging.basicConfig()
			logger.setLevel(logging.DEBUG)

		self.aas = "ACDEFGHIKLMNPQRSTVWY"
		
		self.pssm_data = {}
		self.data = {}
		self.results = {}

		self.options["query_pssm_names"] = []
		self.options["compare_pssm_names"] = []
		
		self.score_type_options = ["psi_blast", "binomial log", "log odds", "log_pseudocount", "ratio_sum"]
		self.specificity_type_options = ["no_specificity_scoring","gini", "gini_normalised", "gini_range"]

	##---------------------------------##

	def calculate_dissimilarity_score(self,column_dissimilarity,column_gini_coefficients_query,column_gini_coefficients_compare):
		if self.options['dissimilarity_score_focus'] == "query":
			dissimilarity_score = (column_dissimilarity*column_gini_coefficients_query)
		elif self.options['dissimilarity_score_focus'] == "compare":
			dissimilarity_score = (column_dissimilarity*column_gini_coefficients_compare)
	
		return dissimilarity_score

	def calculate_similarity_score(self,column_pearson_corr,column_gini_coefficients_query,column_gini_coefficients_compare):
		return (column_pearson_corr)*column_gini_coefficients_query*column_gini_coefficients_compare

	##---------------------------------##

	def normalise_pssm(self,pssm):
		normalised_pssm = pssm

		if self.options["ignore_negative_pssm_values"]:
			normalised_pssm = utilities_pssm.positivise_pssm(normalised_pssm)
		
		if self.options["pssm_normalise_count_sum_to_one"]:
			normalised_pssm = utilities_pssm.pssm_normalise_column_sum_to_one(normalised_pssm)
		elif self.options["pssm_normalise_min_max"]:
			normalised_pssm = utilities_pssm.normalise_pssm_min_max(normalised_pssm)

		return normalised_pssm

	##---------------------------------##

	def process_pssm(self,pssm):
		normalised_columns = {}

		pssm_length = len(pssm[self.aas[0]])

		normalised_pssm = self.normalise_pssm(pssm)

		columns = {}
		normalised_columns = {}
		for i in range(0,pssm_length):
			columns[i] = []
			normalised_columns[i] = []
			for aa in self.aas:
				if aa == "-": continue
				columns[i].append(pssm[aa][i])
				normalised_columns[i].append(normalised_pssm[aa][i])
			
		return {
			"gini_coefficients":self.process_pssm_gini_coefficients(normalised_columns),
			"normalised_columns":normalised_columns,
			"columns":columns
		}
	
	##---------------------------------##

	def process_pssm_gini_coefficients(self,normalised_columns):
		gini_coefficients = {}
		response_gini_coefficients = {}

		for i in range(0,len(normalised_columns)):
			gini_coefficients[i] = utilities_stats.gini_coefficient(normalised_columns[i])
			if gini_coefficients[i] > 1:
				print(gini_coefficients[i] ,normalised_columns[i])

		max_gini_coefficients = max(gini_coefficients.values())
		min_gini_coefficients = min(gini_coefficients.values())

		if self.options['specificity_scoring_type'] == "no_specificity_scoring":
			for i in range(0,len(normalised_columns)):
				response_gini_coefficients[i] = 1
	
		if self.options['specificity_scoring_type'] == "gini":
			response_gini_coefficients = gini_coefficients

		if  self.options['specificity_scoring_type'] == "gini_normalised":
			for i in range(0,len(normalised_columns)):
				try:
					response_gini_coefficients[i] = gini_coefficients[i]/max_gini_coefficients
				except:
					response_gini_coefficients[i] = 0

		if  self.options['specificity_scoring_type'] == "gini_range":
			for i in range(0,len(normalised_columns)):
				try:
					response_gini_coefficients[i] = (gini_coefficients[i]-min_gini_coefficients)/(max_gini_coefficients-min_gini_coefficients)
				except:
					response_gini_coefficients[i] = 0
		
		return response_gini_coefficients
	
	##---------------------------------##

	def sample_column_scores(self):
		logger.debug("Calculating sampled column comparison scores")

		self.data["sampled_column_comparison_scores"] = []
		
		random.seed(1) 

		for i in range(0,self.options["column_scores_samples"]):	
			
			pssm_1_random = random.choice(self.options["query_pssm_names"])
			column_1_list = list(self.pssm_data[pssm_1_random]["normalised_columns"].keys())
			column_1_random = random.choice(column_1_list)

			pssm_2_random = random.choice(self.options["compare_pssm_names"])
			column_2_list = list(self.pssm_data[pssm_2_random]["normalised_columns"].keys())
			column_2_random = random.choice(column_2_list)
			
			pssm_3_random = random.choice(list(self.options["query_pssm_names"]) + list(self.options["compare_pssm_names"]))
			column_3_list = list(self.pssm_data[pssm_3_random]["normalised_columns"].keys())
			column_3_random = random.choice(column_3_list)

			pssm_4_random = random.choice(list(self.options["query_pssm_names"]) + list(self.options["compare_pssm_names"]))
			column_4_list = list(self.pssm_data[pssm_4_random]["normalised_columns"].keys())
			column_4_random = random.choice(column_4_list)

			column_a = self.pssm_data[pssm_1_random]["normalised_columns"][column_1_random]
			column_b = self.pssm_data[pssm_2_random]["normalised_columns"][column_2_random]

		
			column_pearson_corr = utilities_stats.pearson_correlation(column_a, column_b)
			column_gini_coefficients_query = self.pssm_data[pssm_1_random]["gini_coefficients"][column_1_random]
			column_gini_coefficients_compare = self.pssm_data[pssm_2_random]["gini_coefficients"][column_2_random]
				
			similarity_score = self.calculate_similarity_score(column_pearson_corr,column_gini_coefficients_query,column_gini_coefficients_compare) 
			self.data["sampled_column_comparison_scores"].append(similarity_score)		

		###---###
			
		self.data["sampled_column_comparison_scores"].sort()
		
		###---###

		logger.debug("Completed sampled column comparison scores calculation")

	##---------------------------------##

	def sampled_column_comparison_scores_p_value(self,score):
		list_bisect = bisect(self.data["sampled_column_comparison_scores"], score)
		if list_bisect == len(self.data["sampled_column_comparison_scores"]):
			list_bisect = len(self.data["sampled_column_comparison_scores"]) - 1 
		
		p_value = 1 - list_bisect/len(self.data["sampled_column_comparison_scores"])
		return max(p_value,(1/len(self.data["sampled_column_comparison_scores"]))*100) 

	##---------------------------------##

	def randomised_pssms(self):
		for i in range(0,self.options["random_pssm_samples"]):
			self.pssm_data["random_" + str(i)] = {
			"normalised_columns":{}, 
			"gini_coefficients":{}, 
			"motifs":[], 
			"pssm_type":"random"
			}

			for ii in range(0,self.options["random_pssm_length"]):
				motif_list = list(self.pssm_data.keys())
				motif_random = random.choice(motif_list[:-1])
				
				column_list = list(self.pssm_data[motif_random]["normalised_columns"].keys())
				column_random = random.choice(column_list)

				self.pssm_data["random_" + str(i)]["normalised_columns"][ii] = self.pssm_data[motif_random]["normalised_columns"][column_random]
				self.pssm_data["random_" + str(i)]["gini_coefficients"][ii] = self.pssm_data[motif_random]["gini_coefficients"][column_random]
				
				self.pssm_data["random_" + str(i)]["motifs"].append(self.pssm_data[motif_random]["motifs"][column_random])
	
	##---------------------------------##

	def get_peptide_gap_proportion(self,peptides,query_offset_start,query_offset_end):
		residues = 0
		gaps = 0
		for peptide in peptides:
			peptide_segment = peptide[query_offset_start-1:query_offset_end]
			residues +=  len(peptide_segment)
			gaps += peptide_segment.count('-')
		
		return 1 - (gaps/residues)
	
	##---------------------------------##
	
	def compare_pssm_region(self,pssm_1,pssm_2,motif_shorter,motif_longer,start_offset,pssm_offset_pairs,multiple_testing_correction=1):
		
		raw_scores = []
		raw_dissimilarity_scores = []
		raw_column_dissimilarity_scores = []
		p_scores = []
		columns_gini_coefficients_query = []
		columns_gini_coefficients_compare = []
		columns_pearson_corr = []
		motifs = [[],[]]

		pssm_comparisons = {"score":False}

		shorter_start = False

		for pssm_offset_pair in pssm_offset_pairs:
			motif_shorter_i =  pssm_offset_pair[0]
			motif_longer_i = pssm_offset_pair[1]

			if motif_shorter_i >= 0 and motif_longer_i >= 0 and shorter_start == False:
				shorter_start = pssm_offset_pair

			if motif_shorter_i in self.pssm_data[motif_shorter]["normalised_columns"] and motif_longer_i in self.pssm_data[motif_longer]["normalised_columns"]:
				
				column_pearson_corr = utilities_stats.pearson_correlation(self.pssm_data[motif_shorter]["normalised_columns"][motif_shorter_i],self.pssm_data[motif_longer]["normalised_columns"][motif_longer_i])
				column_gini_coefficients_query = self.pssm_data[motif_shorter]["gini_coefficients"][motif_shorter_i]
				column_gini_coefficients_compare = self.pssm_data[motif_longer]["gini_coefficients"][motif_longer_i]

				column_dissimilarity_score = utilities_stats.mean_absolute_error(self.pssm_data[motif_shorter]["normalised_columns"][motif_shorter_i],self.pssm_data[motif_longer]["normalised_columns"][motif_longer_i])
				dissimilarity_score = self.calculate_dissimilarity_score(column_dissimilarity_score,column_gini_coefficients_query,column_gini_coefficients_compare) 

				raw_score = self.calculate_similarity_score(column_pearson_corr,column_gini_coefficients_query,column_gini_coefficients_compare) 
				raw_scores.append(raw_score)
				p_scores.append(self.sampled_column_comparison_scores_p_value(raw_score))

				raw_column_dissimilarity_scores.append(column_dissimilarity_score)
				raw_dissimilarity_scores.append(dissimilarity_score)
				
				columns_gini_coefficients_query.append(column_gini_coefficients_query)
				columns_gini_coefficients_compare.append(column_gini_coefficients_compare)
				columns_pearson_corr.append(column_pearson_corr)

				try:
					motifs[0].append(self.pssm_data[motif_shorter]["motifs"][motif_shorter_i])
				except:
					motifs[0].append("?")

				try:
					motifs[1].append(self.pssm_data[motif_longer]["motifs"][motif_longer_i])
				except:
					motifs[1].append("?")
			
		if len(raw_scores) == 0: return

		#---------------------------#
		
		raw_score_combined = sum(raw_scores)/len(raw_scores)			
		p_score_combined = utilities_stats.cum_uniform_product(functools.reduce(operator.mul, p_scores),len(p_scores))
		p_score_combined_adj = min(p_score_combined*multiple_testing_correction,1)
		score = p_score_combined_adj
		overlap_proportion = len(p_scores)/len(pssm_offset_pairs)
		
		#---------------------------#

		core_gini = []
		flank_gini = []

		focus_pssm = pssm_1
		if self.options['dissimilarity_score_focus'] == 'compare':
			focus_pssm = pssm_2

		if focus_pssm == motif_shorter:
			for i in range(0,len(self.pssm_data[focus_pssm]["gini_coefficients"])):
				if i in range(shorter_start[0],shorter_start[0]+len(raw_scores)):
					core_gini.append(self.pssm_data[focus_pssm]["gini_coefficients"][i])
				else:
					raw_dissimilarity_scores.append(self.pssm_data[focus_pssm]["gini_coefficients"][i])
					flank_gini.append(self.pssm_data[focus_pssm]["gini_coefficients"][i])

		if focus_pssm == motif_longer:
			for i in range(0,len(self.pssm_data[focus_pssm]["gini_coefficients"])):
				if i in range(shorter_start[1],shorter_start[1]+len(raw_scores)):
					core_gini.append(self.pssm_data[focus_pssm]["gini_coefficients"][i])
				else:
					raw_dissimilarity_scores.append(self.pssm_data[focus_pssm]["gini_coefficients"][i])
					flank_gini.append(self.pssm_data[focus_pssm]["gini_coefficients"][i])

		#---------------------------#
			
		significant_positions_check = False
		significant_positions = [p_score < self.options['significant_positions_cutoff'] for p_score in p_scores]
		
		if significant_positions.count(True) >= self.options['significant_positions_count_cutoff']:
			significant_positions_check = True	
		
		try:		
			query_accession = self.pssm_data[pssm_1]["accession"]
		except:
			query_accession = pssm_1
		try:
			compare_accession = self.pssm_data[pssm_2]["accession"]
		except:
			compare_accession = pssm_2

		if self.options['detailed_results']:
			if pssm_2 not in self.results[pssm_1]["other"]:
				self.results[pssm_1]["other"][pssm_2] = {
					"raw_score":raw_score_combined,
					"p_score":p_score_combined,
					"p_score_adj":p_score_combined_adj,
					"length":len(p_scores),
					"overlap_proportion":overlap_proportion,
					"dissimilarity_scores_max":max(raw_dissimilarity_scores),
					"significant_positions_count":significant_positions.count(True),
					"query_name":self.pssm_data[pssm_1]["name"],
					"compare_name":self.pssm_data[pssm_2]["name"],
					"query_accession":query_accession,
					"compare_accession":compare_accession
					}
				
				if "region_peptide" in self.pssm_data[pssm_1] :
					self.results[pssm_1]["other"][pssm_2]['query_region_peptide'] = self.pssm_data[pssm_1]['region_peptide']
				
				if "region_peptide" in self.pssm_data[pssm_2] :
					self.results[pssm_1]["other"][pssm_2]['compare_region_peptide'] = self.pssm_data[pssm_2]['region_peptide']

				if "region_start" in self.pssm_data[pssm_1] and "region_end" in self.pssm_data[pssm_1]:
					self.results[pssm_1]["other"][pssm_2]['query_region_start'] = self.pssm_data[pssm_1]['region_start']
					self.results[pssm_1]["other"][pssm_2]['query_region_end'] = self.pssm_data[pssm_1]['region_end']

				if "region_start" in self.pssm_data[pssm_2] and "region_end" in self.pssm_data[pssm_2]:
					self.results[pssm_1]["other"][pssm_2]['compare_region_start'] = self.pssm_data[pssm_2]['region_start']
					self.results[pssm_1]["other"][pssm_2]['compare_region_end'] = self.pssm_data[pssm_2]['region_end']

			elif p_score_combined < self.results[pssm_1]["other"][pssm_2]['p_score']:
				self.results[pssm_1]["other"][pssm_2].update({
					"raw_score":raw_score_combined,
					"p_score":p_score_combined,
					"p_score_adj":p_score_combined_adj,
					"length":len(p_scores),
					"overlap_proportion":overlap_proportion,
					"dissimilarity_scores_max":max(raw_dissimilarity_scores),
					"significant_positions_count":significant_positions.count(True),
					"query_name":self.pssm_data[pssm_1]["name"],
					"compare_name":self.pssm_data[pssm_2]["name"],
					"query_accession":query_accession,
					"compare_accession":compare_accession
					})

		if pssm_comparisons["score"] == False or score < pssm_comparisons["score"]:
			pssm_comparisons = {
				"score":score,
				"significant_positions_check":significant_positions_check,
				"significant_positions":significant_positions,
				"significant_positions_count":significant_positions.count(True),
				"raw_score":raw_score_combined,
				"p_score":p_score_combined,
				"p_score_adj":p_score_combined_adj,
				"p_scores":p_scores,
				"dissimilarity_raw_scores":raw_dissimilarity_scores,
				"raw_column_dissimilarity_scores":raw_column_dissimilarity_scores,
				"dissimilarity_scores_max":max(raw_dissimilarity_scores),
				"columns_gini_coefficients_compare":columns_gini_coefficients_compare,
				"columns_gini_coefficients_query":columns_gini_coefficients_query,
				"columns_pearson_corr":columns_pearson_corr,
				"raw_scores":raw_scores,
			}

			if pssm_1 == motif_shorter:
				pssm_comparisons["query_motif_re"] = "".join(motifs[0])
				pssm_comparisons["compare_motif_re"] = "".join(motifs[1])
				pssm_comparisons["query_offset_start"] = shorter_start[0]
				pssm_comparisons["query_offset_end"] = shorter_start[0] + len(raw_scores)
				pssm_comparisons["compare_offset_start"] = shorter_start[1]
				pssm_comparisons["compare_offset_end"] = shorter_start[1] + len(raw_scores)
			else:
				pssm_comparisons["query_motif_re"] = "".join(motifs[1])
				pssm_comparisons["compare_motif_re"] = "".join(motifs[0])
				pssm_comparisons["query_offset_start"] = shorter_start[1]
				pssm_comparisons["query_offset_end"] =  shorter_start[1] + len(raw_scores)
				pssm_comparisons["compare_offset_start"] = shorter_start[0]
				pssm_comparisons["compare_offset_end"] = shorter_start[0] + len(raw_scores)
			
	
			if (self.results[pssm_1]["best"]["score"] == False or score < self.results[pssm_1]["best"]["score"]) and significant_positions_check:
				self.results[pssm_1]["best"] = pssm_comparisons
				self.results[pssm_1]["best"]["core_gini"] = core_gini
				self.results[pssm_1]["best"]["flank_gini"] = flank_gini
				self.results[pssm_1]["best"]["query_motif"] = pssm_1
				self.results[pssm_1]["best"]["compare_motif"] = pssm_2
				self.results[pssm_1]["best"]["query_name"] = self.pssm_data[pssm_1]["name"]
				self.results[pssm_1]["best"]["compare_name"] = self.pssm_data[pssm_2]["name"]
				self.results[pssm_1]["best"]["compare_pssm_score_scheme"] = self.options["compare_pssm_score_scheme"]

				if self.options["query_dataset"] not in ["proppd_dms","folder_pssms","pssm_file"]:
					self.results[pssm_1]["best"]["query_pssm_score_scheme"] = self.options["query_pssm_score_scheme"]
			
			#####------------------–#####

			if score < self.options['significance_cutoff'] and significant_positions_check:
				if self.options['detailed_results']:
					self.results[pssm_1]["significant"][pssm_2] = pssm_comparisons
					self.results[pssm_1]["significant"][pssm_2]["query_motif"] = pssm_1
					self.results[pssm_1]["significant"][pssm_2]["compare_motif"] = pssm_2
					self.results[pssm_1]["significant"][pssm_2]["query_name"] = self.pssm_data[pssm_1]["name"]
					self.results[pssm_1]["significant"][pssm_2]["compare_name"] = self.pssm_data[pssm_2]["name"]
				else:
					self.results[pssm_1]["significant"][pssm_2] = {
						"score":score,
						"p_score":p_score_combined,
						"p_score_adj":p_score_combined_adj,
						"dissimilarity_scores_max":pssm_comparisons["dissimilarity_scores_max"],
						"query_name":self.pssm_data[pssm_1]["name"],
						"compare_name":self.pssm_data[pssm_2]["name"],
						"compare_motif_re":pssm_comparisons["compare_motif_re"],
						"query_motif_re":pssm_comparisons["query_motif_re"]
					}

    ##---------------------------------##
									
	def compare_pssms(self):
		logger.debug("Comparing PSSMs: " + str(len(self.options["query_pssm_names"])) + " -> " + str(len(self.options["compare_pssm_names"])))		
		rows = ["\t".join(["sim_p","dissim","q_re","c_re","q_name","c_name"])]

		if len(self.options["query_pssm_names"]) == 0 or len(self.options["compare_pssm_names"]) == 0:
			logger.error('Nothing to compare')
			return

		for pssm_1 in self.options["query_pssm_names"]:
			
			self.results[pssm_1] = {
					"best":{"score":False},
					"significant":{}
				}
			
			if self.options['detailed_results']:
				self.results[pssm_1]["other"] = {}

			for pssm_2 in self.options["compare_pssm_names"]:
				if len(self.pssm_data[pssm_2]["normalised_columns"]) > len(self.pssm_data[pssm_1]["normalised_columns"]):
					motif_shorter = pssm_1
					motif_longer = pssm_2
				else:
					motif_shorter = pssm_2
					motif_longer = pssm_1
					
				if self.options['sliding_window_comparison']:
					multiple_testing_correction = len(self.pssm_data[motif_longer]["normalised_columns"])
				
					for start_offset in range(-(len(self.pssm_data[motif_shorter]["normalised_columns"])-self.options['min_overlap']) + 1,len(self.pssm_data[motif_longer]["normalised_columns"])-self.options['min_overlap'] + 1):
						pssm_offset_pairs = []
						for pssm_iter in range(0,len(self.pssm_data[motif_shorter]["normalised_columns"])):
							motif_shorter_i =  pssm_iter
							motif_longer_i = start_offset + pssm_iter
					
							pssm_offset_pairs.append([motif_shorter_i,motif_longer_i])
					
						self.compare_pssm_region(pssm_1,pssm_2,motif_shorter,motif_longer,start_offset,pssm_offset_pairs,multiple_testing_correction)
				else:
					pssm_offset_pairs = []
					
					for pssm_iter in range(0,len(self.pssm_data[motif_shorter]["normalised_columns"])):
						motif_shorter_i =  pssm_iter
						motif_longer_i = pssm_iter
						pssm_offset_pairs.append([motif_shorter_i,motif_longer_i])

					self.compare_pssm_region(pssm_1,pssm_2,motif_shorter,motif_longer,0,pssm_offset_pairs,1)
			
			logger.debug("Compared " + pssm_1 + " " + str(len(self.results)) + "/" + str(len(self.options["query_pssm_names"])) + " - significant hits:" + str(len(self.results[pssm_1]["significant"])))
			
			##############################################################################
			
			try:
				if self.options['include_query_pssm']:
					self.results[pssm_1]["query_pssm"] = self.pssm_data[pssm_1]["pssm"]
					
				if self.options['include_pssms']:
					self.results[pssm_1]["best"]["query_pssm"] = utilities_pssm.cut_pssm(self.pssm_data[pssm_1]["pssm"],self.results[pssm_1]["best"]["query_offset_start"],self.results[pssm_1]["best"]["query_offset_end"])
					
					self.results[pssm_1]["best"]["compare_pssm"] = utilities_pssm.cut_pssm(self.pssm_data[self.results[pssm_1]["best"]["compare_motif"]]["pssm"],self.results[pssm_1]["best"]["compare_offset_start"],self.results[pssm_1]["best"]["compare_offset_end"])

					self.results[pssm_1]["best"]["query_pssm_normalised"] = self.normalise_pssm(self.results[pssm_1]["best"]["query_pssm"] )
					self.results[pssm_1]["best"]["compare_pssm_normalised"] = self.normalise_pssm(self.results[pssm_1]["best"]["compare_pssm"] )
					
				##############################################################################			

				row = "\t".join([
					"%1.3g"%self.results[pssm_1]["best"]["score"],
					"%1.3f"%self.results[pssm_1]["best"]["dissimilarity_scores_max"],
					self.results[pssm_1]["best"]["query_motif_re"],
					self.results[pssm_1]["best"]["compare_motif_re"],
					self.pssm_data[pssm_1]["name"],
					self.pssm_data[self.results[pssm_1]["best"]["compare_motif"]]["name"],
				])


				
				rows.append(row)
			except:
				logger.error("Can't find best hit: " + pssm_1)
				utilities_error.printError()

			
			if self.options['verbose']:
				print("\n".join(rows))
	
	##---------------------------------####-------------------------####-------------------------##
	##---------------------------------####-------------------------####-------------------------##

	def setup_pssm_json(self,pssm_type):
		if pssm_type == "query":
		
			pssms = self.options['query_pssm']
			
		if pssm_type == "compare":
			pssms = self.options['compare_pssm']
			
		for pssm_id in pssms:
			self.pssm_data[pssm_type + "_" + pssm_id] = {
				"pssm_type":pssm_type,
				"name":pssm_id,
				"motifs":utilities_pssm.get_pssm_motif(pssms[pssm_id],cut_off=0.2),
				"pssm": pssms[pssm_id]
			}


	##---------------------------------####-------------------------####-------------------------##
	##---------------------------------####-------------------------####-------------------------##

	def setup_pssms(self):
		logger.debug("Setting up query data - " + self.options["query_dataset"])
		
		self.setup_pssm_json("query")
		self.setup_pssm_json("compare")
		
		self.mask_pssms()

	##---------------------------------##
	
	def mask_pssms(self):
		if len(self.options['mask_positions']) > 0:
			for pssm_id in self.pssm_data:
				for mask_position in self.options['mask_positions']:
					mask_position = int(mask_position)
					for aa in self.pssm_data[pssm_id]["pssm"]:
						self.pssm_data[pssm_id]["pssm"][aa][mask_position] = self.options['mask_score']

	##---------------------------------##

	def process_pssms(self):
		for pssm_id in self.pssm_data:
			logger.debug("Processing " + pssm_id)
			processed_pssm = self.process_pssm( self.pssm_data[pssm_id]["pssm"])
			self.pssm_data[pssm_id]["gini_coefficients"] = processed_pssm["gini_coefficients"]
			self.pssm_data[pssm_id]["normalised_columns"] = processed_pssm["normalised_columns"]

			if self.pssm_data[pssm_id]["pssm_type"] == "query":
				self.options["query_pssm_names"].append(pssm_id)
			elif self.pssm_data[pssm_id]["pssm_type"] == "compare":
				self.options["compare_pssm_names"].append(pssm_id)
			else:
				self.options["query_pssm_names"].append(pssm_id)

		self.options["query_pssm_names"].sort()
		self.options["compare_pssm_names"].sort()

		if len(self.options["query_pssm_names"]) == 0 or self.options["compare_pssm_names"] == 0:
			logger.error("No PSSMs to compare")
			logger.error(self.options["query_pssm_names"])
			logger.error(self.options["compare_pssm_names"])
			return

		self.randomised_pssms()
		self.sample_column_scores()

	##---------------------------------##

	def run_pssm_comparison_analysis(self):
		try:
			if (len(self.options["query_pssm"]) == 0 or len(self.options["compare_pssm"]) == 0) and (len(self.options["query_pssm_file"]) == 0 or len(self.options["compare_pssm_file"]) == 0):
				logger.error("Input PSSMs not set. Add using --query_pssm_file and --compare_pssm_file or (in function query_pssm and compare_pssm)")
				sys.exit()
			
			if (len(self.options["query_pssm_file"]) != 0 and len(self.options["compare_pssm_file"]) != 0):
				if not os.path.exists(self.options["query_pssm_file"]):
					logger.error(self.options["query_pssm_file"] + " output file does not exist ")
					sys.exit()
				
				if not os.path.exists(self.options["query_pssm_file"]):
					logger.error(self.options["query_pssm_file"] + " output file does not exist ")
					sys.exit()
				
				self.options['query_pssm'] = json.loads(open(self.options['query_pssm_file']).read())
				self.options['compare_pssm'] = json.loads(open(self.options['compare_pssm_file']).read())

				for pssm_type in ['query_pssm',"compare_pssm"]:
					for pssm_id in self.options[pssm_type]:
						pssm_status = utilities_pssm.check_pssm(self.options[pssm_type][pssm_id])
						if pssm_status['status'] == "error":
							logger.error(pssm_type + ":" + pssm_id + " " + pssm_status['error_type'])
							sys.exit()
						else:
							logger.debug(pssm_type + ":" + pssm_id + " PSSM check successful")
							
			if (len(self.options["output_file"]) == 0) and (len(self.options["query_pssm_file"]) != 0 and len(self.options["compare_pssm_file"]) != 0):
				logger.error("Output file not set. Add using --output_file ")
				sys.exit()
			
			self.setup_pssms()
			self.process_pssms()
			self.compare_pssms()
		except:
			logger.error("Error running compariPSSM")
			self.results = utilities_error.getError()
			raise

		if len(self.options["output_file"]) > 0:
			with open(self.options["output_file"], 'w') as f:
				json.dump( self.results, f)

		return self.results

	##---------------------------------##

if __name__ == "__main__":
	comparipssm_runner = CompariPSSM()
	pssm_comparison_response = comparipssm_runner.run_pssm_comparison_analysis()
	