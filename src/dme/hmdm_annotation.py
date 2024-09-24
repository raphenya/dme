import os, sys, json, csv, argparse
"""
This script it used to create annotations and fasta
"""
def main(args):
	working_directory = os.getcwd()
	"""
	reading hmdm.json
	"""
	with open(os.path.join(args.input), 'r') as jfile:
		data = json.load(jfile)

	# get version
	try:
		version = data["_version"]
	except Exception as e:
		print("Error: missing version number")
		exit()

	annotations = []
	"""
	write card reference fasta (FASTA format)
	"""
	with open(os.path.join(working_directory, "hmdm_database_v{}.fasta".format(version)), 'w') as fout:
		for i in data:
			if i.isdigit():
				# use homolog models
				if data[i]['model_type_id'] in ['54']:

					drug_class = []
					mechanism = []
					group = []

					if "HMDM_category" in data[i]:
						for c in data[i]["HMDM_category"]:
							if "category_aro_class_name" in data[i]["HMDM_category"][c]: 
								if data[i]["HMDM_category"][c]["category_aro_class_name"] in ["Drug Class"]:
									drug_class.append(("{}".format(data[i]["HMDM_category"][c]["category_aro_name"])))
								if data[i]["HMDM_category"][c]["category_aro_class_name"] in ["Resistance Mechanism"]:
									mechanism.append(("{}".format(data[i]["HMDM_category"][c]["category_aro_name"])))
								if data[i]["HMDM_category"][c]["category_aro_class_name"] in ["Gene Family"]:
									group.append(("{}".format(data[i]["HMDM_category"][c]["category_aro_name"])))
					try:
						for seq in data[i]['model_sequences']['sequence']:
							if args.ncbi == True:
								# header used to be able to validate CARD sequences with genbank sequences
								header = ("gb|{ncbi}|HMDM:{HMDM_accession}|ID:{model_id}|Name:{HMDM_name}".format(
									ARO_accession=data[i]['ARO_accession'],
									model_id=data[i]['model_id'],
									ARO_name=(data[i]['HMDM_name']).replace(" ", "_"),
									ncbi=data[i]['model_sequences']['sequence'][seq]["dna_sequence"]["accession"]
									))
							else:
								header = ("HMDM:{}|ID:{}|Name:{}|NCBI:{}".format(
									data[i]['HMDM_accession'],
									data[i]['model_id'],
									(data[i]['HMDM_name']).replace(" ", "_"),
									data[i]['model_sequences']['sequence'][seq]["dna_sequence"]["accession"]
									))

							sequence = data[i]['model_sequences']['sequence'][seq]["dna_sequence"]["sequence"]
							fout.write(">{}\n".format(header))
							fout.write("{}\n".format(sequence))
							annotations.append([header, "; ".join(drug_class), "; ".join(mechanism), "; ".join(group)])

					except Exception as e:
						print("No model sequences for model ({}, {}). Omitting this model and keep running.".format(data[i]['model_id'], data[i]['model_name']))

def create_parser():
    parser = argparse.ArgumentParser(prog="rgi hmdm_annotation",description='Creates HMDM annotations for DME BWT from hmdm.json')
    parser.add_argument('-i', '--input', dest="input", required=True, help="hmdm.json file")
    parser.add_argument('--ncbi', dest="ncbi", action="store_true", help="adds ncbi accession to FASTA headers")
    return parser

def run():
    parser = create_parser()
    args = parser.parse_args()
    main(args)

if __name__ == '__main__':
    run()

