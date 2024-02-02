import json, pprint

import comparipssm

##---------------------------------------------------------------------##

if __name__ == "__main__":
	comparipssm_runner = comparipssm.CompariPSSM()
	
	query_pssm = json.loads(open('./pssm_sets/query_pssm.json').read())
	compare_pssm = json.loads(open('./pssm_sets/elm_pssm.json').read())

	comparipssm_runner.options["query_pssm"] = query_pssm
	comparipssm_runner.options["compare_pssm"] = compare_pssm
	comparipssm_runner.options["detailed_results"] = False
	comparipssm_runner.options["include_pssms"] = False
	comparipssm_runner.options["significance_cutoff"] = 0.0001

	pssm_comparison_response = comparipssm_runner.run_pssm_comparison_analysis()
	pprint.pprint(pssm_comparison_response)
