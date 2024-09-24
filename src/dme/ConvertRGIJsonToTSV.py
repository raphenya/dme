import csv
from dme.settings import *
from operator import itemgetter, attrgetter
from collections import OrderedDict
import Bio
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import subprocess


class ConvertJsonToTSV(object):

    def __init__(self, filepath, homolog_file=None, variant_file=None, overexpression_file=None, rrna_file=None,
                 local_database=False):
        f_path, f_name = os.path.split(filepath)
        name, ext = os.path.splitext(f_name)
        self.filepath = os.path.join(f_path, "{}.json".format(f_name))
        if ext.lower() == ".json":
            self.filepath = os.path.join(f_path, "{}{}".format(name, ext))
        self.homolog_file = homolog_file
        self.variant_file = variant_file
        self.overexpression_file = overexpression_file
        self.rrna_file = rrna_file
        self.db = path
        if local_database:
            self.db = LOCAL_DATABASE

    def __repr__(self):
        """Returns ConvertJsonToTSV class full object."""
        return "ConvertJsonToTSV({}".format(self.__dict__)

    def parse_jsons(self, part, data):
        if part is not None:
            if data.keys():
                for k1, v1 in part.items():
                    if k1 in data:
                        data[k1].update(v1)
                    elif k1 not in data:
                        data[k1] = v1
            else:
                data.update(part)
        return data

    def combine_jsons(self):
        if self.homolog_file is not None and self.variant_file is not None and self.overexpression_file is not None and self.rrna_file is not None:
            data = {}
            if os.path.isfile(self.homolog_file):
                with open(self.homolog_file) as hf:
                    homolog = json.load(hf)
                    data = self.parse_jsons(homolog, data)
            if os.path.isfile(self.variant_file):
                with open(self.variant_file) as vf:
                    variant = json.load(vf)
                    data = self.parse_jsons(variant, data)

            if os.path.isfile(self.overexpression_file):
                with open(self.overexpression_file) as of:
                    overexpression = json.load(of)
                    data = self.parse_jsons(overexpression, data)

            if os.path.isfile(self.rrna_file):
                with open(self.rrna_file) as rf:
                    rrna = json.load(rf)
                    # if rrna is not None:
                    data = self.parse_jsons(rrna, data)

            with open(self.filepath, "w") as outfile:
                json.dump(data, outfile)

    def make_blast_database_1(self, dbtype="strains"):
        """Build BLAST database from a FASTA file."""
        if os.path.isfile(os.path.join(self.db, "{}db.fsa".format(dbtype))) == True and os.path.exists(
                os.path.join(self.db, "{}db.fsa".format(dbtype))) == True \
                and os.path.exists(os.path.join(self.db, "{}.db.nin".format(dbtype))) == True and os.path.exists(
            os.path.join(self.db, "{}.db.nhr".format(dbtype))) == True \
                and os.path.exists(os.path.join(self.db, "{}.db.nsq".format(dbtype))) == True:
            logger.info("blast DB (strains) exists.")
        else:
            logger.info("create blast DB (strains).")
            type = "nucl"
            os.system(
                'makeblastdb -in {} -dbtype {type} -out {}'.format(os.path.join(self.db, "{}db.fsa".format(dbtype)),
                                                                   os.path.join(
                                                                       self.db, "{}.db".format(dbtype)),
                                                                   type=type))

    def get_strain(self, sequence, orf=""):
        if sequence == "":
            return "N/A"
        else:
            # check if strain db exists if not, create it
            self.make_blast_database_1("strains")
            f_path, f_name = os.path.split(self.filepath)
            output = os.path.join(f_path, "{}.temp.kmers.xml".format(
                os.path.splitext(f_name)[0]))
            query_string = '>' + orf + '\n' + str(sequence)
            # blast sequence against custom strains database
            # evalue=0.001)
            blast_command = Bio.Blast.Applications.NcbiblastnCommandline(cmd='blastn',
                                                                         out=output,
                                                                         outfmt=5,
                                                                         db=os.path.join(self.db,
                                                                                         "{}.db".format("strains")),
                                                                         max_target_seqs=1)
            stdout, stderr = blast_command(stdin=query_string)
            # subprocess.run(blast_command)
            # print("stdout: ", stdout)
            # print("stderr: ", stderr)
            strain = ""
            percent_identity = ""
            # read results and return 1 top hit to strains database
            with open(output, 'r') as result_handle:
                blast_records = NCBIXML.read(result_handle)
                for alignment in blast_records.alignments:
                    # strain = alignment.hit_def
                    for hsp in alignment.hsps:
                        # print(alignment.hit_def, " => ", hsp.identities)
                        percent_identity = float(
                            format(float(hsp.identities * 100) / len(hsp.query), '.2f'))
                        strain = alignment.hit_def
            return [strain, percent_identity]

    def run(self):
        if os.path.isfile(self.filepath):
            f_path, f_name = os.path.split(self.filepath)
            with open(os.path.join(f_path, "{}.txt".format(os.path.splitext(f_name)[0])), "w") as af:
                writer = csv.writer(af, delimiter='\t', dialect='excel')
                writer.writerow(["ORF_ID",
                                 "Contig",
                                 "Start",
                                 "Stop",
                                 "Orientation",
                                 "Cut_Off",
                                 "Pass_Bitscore",
                                 "Best_Hit_Bitscore",
                                 "Best_Hit_HMDM",
                                 "Best_Identities",
                                 "HMDM",
                                 "Model_type",
                                 "SNPs_in_Best_Hit_HMDM",
                                 "Other_SNPs",
                                 "Drug Class",
                                 "Mechanism",
                                 "Gene Family",
                                 "Predicted_DNA",
                                 "Predicted_Protein",
                                 "HMDM_Protein_Sequence",
                                 "Percentage Length of Reference Sequence",
                                 "ID",
                                 "Model_ID",
                                 "Nudged",
                                 "Note",
                                 "Percent_Positives",
                                 "Drug",
                                 "Strain",
                                 "Strain_Indentites"])

                if os.path.isfile(self.filepath):
                    with open(self.filepath) as rgi_file:
                        rgi_data = json.load(rgi_file)
                    try:
                        del rgi_data["_metadata"]
                    except:
                        pass

                    for hsp in rgi_data:
                        order_perfect = []
                        order_loose = []
                        order_strict = []
                        dna = 0
                        cgList = []
                        hitID = []
                        temp2 = []
                        temp3 = []
                        best_snps = ""
                        other_snps = ""

                        nudged = ""
                        note = ""
                        orf_start_possible = ""
                        orf_end_possible = ""
                        orf_dna_sequence_possible = ""
                        orf_prot_sequence_possible = ""

                        for hit in rgi_data[hsp]:
                            if rgi_data[hsp][hit]["type_match"] == "Perfect":
                                order_perfect.append((
                                    hit, rgi_data[hsp][hit]["bit_score"], rgi_data[hsp][hit]["perc_identity"]
                                ))

                            if rgi_data[hsp][hit]["type_match"] == "Strict":
                                order_strict.append((
                                    hit, rgi_data[hsp][hit]["bit_score"], rgi_data[hsp][hit]["perc_identity"]
                                ))

                            if rgi_data[hsp][hit]["type_match"] == "Loose":
                                order_loose.append((
                                    hit, rgi_data[hsp][hit]["bit_score"], rgi_data[hsp][hit]["perc_identity"]
                                ))

                        ordered = []
                        # sort by bitscore and percent identity
                        if len(order_perfect) > 0:
                            ordered = sorted(order_perfect, key=itemgetter(
                                1, 2), reverse=True).pop(0)
                        elif len(order_strict) > 0:
                            ordered = sorted(order_strict, key=itemgetter(
                                1, 2), reverse=True).pop(0)
                        else:
                            ordered = sorted(order_loose, key=itemgetter(
                                1, 2), reverse=True).pop(0)

                        ordered = [i for i in ordered]

                        if "orf_dna_sequence" in rgi_data[hsp][ordered[0]]:
                            dna = 1
                        if "HMDM_category" in rgi_data[hsp][ordered[0]]:
                            for aroctkey in rgi_data[hsp][hit]["HMDM_category"]:
                                cgList.append(str(
                                    rgi_data[hsp][hit]["HMDM_category"][aroctkey]["category_hmdm_name"].encode('ascii',
                                                                                                               'replace').decode(
                                        "utf-8")))
                        if "hsp_num:" in rgi_data[hsp][ordered[0]]:
                            hitID.append(rgi_data[hsp][ordered[0]])

                        match_dict = {}

                        if "nudged" in rgi_data[hsp][ordered[0]].keys():
                            nudged = rgi_data[hsp][ordered[0]]["nudged"]

                        if "orf_start_possible" in rgi_data[hsp][ordered[0]].keys():
                            orf_start_possible = rgi_data[hsp][ordered[0]
                                                               ]["orf_start_possible"]

                        if "orf_end_possible" in rgi_data[hsp][ordered[0]].keys():
                            orf_end_possible = rgi_data[hsp][ordered[0]
                                                             ]["orf_end_possible"]

                        if "note" in rgi_data[hsp][ordered[0]].keys():
                            note = rgi_data[hsp][ordered[0]]["note"]

                        if "orf_dna_sequence_possible" in rgi_data[hsp][ordered[0]].keys():
                            orf_dna_sequence_possible = rgi_data[hsp][ordered[0]
                                                                      ]["orf_dna_sequence_possible"]

                        if "orf_prot_sequence_possible" in rgi_data[hsp][ordered[0]].keys():
                            orf_prot_sequence_possible = rgi_data[hsp][ordered[0]
                                                                       ]["orf_prot_sequence_possible"]

                        if dna == 1:
                            if nudged == True and rgi_data[hsp][ordered[0]]["type_match"] == "Perfect":
                                orf_start = orf_start_possible
                                orf_end = orf_end_possible
                                orf_dna = orf_dna_sequence_possible
                                orf_prot = orf_prot_sequence_possible
                            else:
                                orf_start = rgi_data[hsp][ordered[0]
                                                          ]["orf_start"]
                                orf_end = rgi_data[hsp][ordered[0]]["orf_end"]
                                orf_dna = rgi_data[hsp][ordered[0]
                                                        ]["orf_dna_sequence"]
                                orf_prot = rgi_data[hsp][ordered[0]
                                                         ]["orf_prot_sequence"]

                            if len(rgi_data[hsp]) != 0:
                                if rgi_data[hsp][hit]["model_type_id"] == 41091:
                                    if "snp" in rgi_data[hsp][ordered[0]]:
                                        for x in rgi_data[hsp].values():
                                            if "snp" in x.keys():
                                                if x['model_id'] == rgi_data[hsp][ordered[0]]['model_id']:
                                                    temp2.append(
                                                        x["snp"]["original"] + str(x["snp"]["position"]) + x["snp"][
                                                            "change"])
                                                    best_snps = ', '.join(
                                                        temp2)
                                                else:
                                                    temp3.append(
                                                        x["snp"]["original"] + str(x["snp"]["position"]) + x["snp"][
                                                            "change"] + ":" + x['model_id'])
                                                    other_snps = ', '.join(
                                                        temp3)
                                    else:
                                        best_snps = "n/a"
                                        other_snps = "n/a"
                                elif rgi_data[hsp][hit]["model_type_id"] in [40293, 40295]:
                                    if "snp" in rgi_data[hsp][ordered[0]]:
                                        for x in rgi_data[hsp].values():
                                            if "snp" in x.keys():
                                                if x['model_id'] == rgi_data[hsp][ordered[0]]['model_id']:
                                                    temp2.append(
                                                        x["snp"]["original"] + str(x["snp"]["position"]) + x["snp"][
                                                            "change"])
                                                    best_snps = ', '.join(
                                                        temp2)
                                                else:
                                                    temp3.append(
                                                        x["snp"]["original"] + str(x["snp"]["position"]) + x["snp"][
                                                            "change"] + ":" + x['model_id'])
                                                    other_snps = ', '.join(
                                                        temp3)
                                            # add unique snps
                                            temp2 = list(
                                                OrderedDict.fromkeys(temp2))
                                            best_snps = ', '.join(temp2)
                                            temp3 = list(
                                                OrderedDict.fromkeys(temp3))
                                            other_snps = ', '.join(temp3)
                                    else:
                                        best_snps = "n/a"
                                        other_snps = "n/a"
                                elif rgi_data[hsp][hit]["model_type_id"] == 54:
                                    best_snps = "n/a"
                                    other_snps = "n/a"
                                if not other_snps:
                                    other_snps = "n/a"

                                if rgi_data[hsp][hit]["model_type_id"] in [40295]:
                                    percentage_length_reference_sequence = format(((orf_end - orf_start) /
                                                                                   len(rgi_data[hsp][ordered[0]][
                                                                                       "dna_sequence_from_hmdm"])) * 100,
                                                                                  '.2f')
                                else:
                                    percentage_length_reference_sequence = format(
                                        (len(rgi_data[hsp][ordered[0]]["orf_prot_sequence"]) /
                                         len(rgi_data[hsp][ordered[0]]["sequence_from_hmdm"])) * 100, '.2f')

                                strain = ""
                                strain_percent_identity = ""
                                # only populate strain information for Perfect and Strict hits
                                if rgi_data[hsp][ordered[0]]["type_match"] in ["Perfect", "Strict"]:
                                    strain, strain_percent_identity = self.get_strain(
                                        orf_dna, rgi_data[hsp][ordered[0]]["orf_from"])

                                match_dict[hsp] = [hsp,
                                                   rgi_data[hsp][ordered[0]
                                                                 ]["orf_from"],
                                                   orf_start,
                                                   orf_end,
                                                   rgi_data[hsp][ordered[0]
                                                                 ]["orf_strand"],
                                                   rgi_data[hsp][ordered[0]
                                                                 ]["type_match"],
                                                   rgi_data[hsp][ordered[0]
                                                                 ]["pass_bitscore"],
                                                   rgi_data[hsp][ordered[0]
                                                                 ]["bit_score"],
                                                   rgi_data[hsp][ordered[0]
                                                                 ]["HMDM_name"],
                                                   rgi_data[hsp][ordered[0]
                                                                 ]["perc_identity"],
                                                   rgi_data[hsp][ordered[0]
                                                                 ]["HMDM_accession"],
                                                   rgi_data[hsp][ordered[0]
                                                                 ]["model_type"],
                                                   best_snps,
                                                   other_snps,
                                                   "; ".join(rgi_data[hsp][ordered[0]]["HMDM_category"][x][
                                                       "category_hmdm_name"] for x in
                                                       rgi_data[hsp][ordered[0]
                                                                     ]["HMDM_category"]
                                                       if rgi_data[hsp][ordered[0]]["HMDM_category"][x][
                                                       "category_hmdm_class_name"] == 'Drug Class'),
                                                   "; ".join(rgi_data[hsp][ordered[0]]["HMDM_category"][x][
                                                       "category_hmdm_name"] for x in
                                                       rgi_data[hsp][ordered[0]
                                                                     ]["HMDM_category"]
                                                       if rgi_data[hsp][ordered[0]]["HMDM_category"][x][
                                                       "category_hmdm_class_name"] == 'Mechanism'),
                                                   "; ".join(rgi_data[hsp][ordered[0]]["HMDM_category"][x][
                                                       "category_hmdm_name"] for x in
                                                       rgi_data[hsp][ordered[0]
                                                                     ]["HMDM_category"]
                                                       if rgi_data[hsp][ordered[0]]["HMDM_category"][x][
                                                       "category_hmdm_class_name"] == 'Gene Family'),
                                                   orf_dna,
                                                   orf_prot,
                                                   rgi_data[hsp][ordered[0]
                                                                 ]["sequence_from_hmdm"],
                                                   # length of hsps / length reference
                                                   percentage_length_reference_sequence,
                                                   ordered[0],
                                                   rgi_data[hsp][ordered[0]
                                                                 ]["model_id"],
                                                   nudged,
                                                   note,
                                                   format((rgi_data[hsp][ordered[0]]["positives"] / len(
                                                       rgi_data[hsp][ordered[0]]["sequence_from_hmdm"])) * 100, '.2f'),
                                                   "; ".join(rgi_data[hsp][ordered[0]]["HMDM_category"][x][
                                                       "category_hmdm_name"] for x in
                                                       rgi_data[hsp][ordered[0]]["HMDM_category"] \
                                                       if rgi_data[hsp][ordered[0]]["HMDM_category"][x][
                                                       "category_hmdm_class_name"] == 'Drug'),
                                                   strain,
                                                   strain_percent_identity
                                                   ]
                            for key, value in match_dict.items():
                                writer.writerow(value)

                        else:
                            if len(rgi_data[hsp]) != 0:
                                if rgi_data[hsp][hit]["model_type_id"] == 41091:
                                    if "snp" in rgi_data[hsp][ordered[0]]:
                                        for x in rgi_data[hsp].values():
                                            if "snp" in x.keys():
                                                if x['model_id'] == rgi_data[hsp][ordered[0]]['model_id']:
                                                    temp2.append(
                                                        x["snp"]["original"] + str(x["snp"]["position"]) + x["snp"][
                                                            "change"])
                                                    best_snps = ', '.join(
                                                        temp2)
                                                else:
                                                    temp3.append(
                                                        x["snp"]["original"] + str(x["snp"]["position"]) + x["snp"][
                                                            "change"] + ":" + x['model_id'])
                                                    other_snps = ', '.join(
                                                        temp3)
                                    else:
                                        best_snps = "n/a"
                                        other_snps = "n/a"
                                elif rgi_data[hsp][hit]["model_type_id"] == 40293:
                                    if "snp" in rgi_data[hsp][ordered[0]]:
                                        for x in rgi_data[hsp].values():
                                            if "snp" in x.keys():
                                                if x['model_id'] == rgi_data[hsp][ordered[0]]['model_id']:
                                                    temp2.append(
                                                        x["snp"]["original"] + str(x["snp"]["position"]) + x["snp"][
                                                            "change"])
                                                    best_snps = ', '.join(
                                                        temp2)
                                                else:
                                                    temp3.append(
                                                        x["snp"]["original"] + str(x["snp"]["position"]) + x["snp"][
                                                            "change"] + ":" + x['model_id'])
                                                    other_snps = ', '.join(
                                                        temp3)
                                    else:
                                        best_snps = "n/a"
                                        other_snps = "n/a"
                                elif rgi_data[hsp][hit]["model_type_id"] == 54:
                                    best_snps = "n/a"
                                    other_snps = "n/a"

                                match_dict[hsp] = [hsp, "", "", "", "",
                                                   rgi_data[hsp][ordered[0]
                                                                 ]["type_match"],
                                                   rgi_data[hsp][ordered[0]
                                                                 ]["pass_bitscore"],
                                                   rgi_data[hsp][ordered[0]
                                                                 ]["bit_score"],
                                                   rgi_data[hsp][ordered[0]
                                                                 ]["HMDM_name"],
                                                   rgi_data[hsp][ordered[0]
                                                                 ]["perc_identity"],
                                                   rgi_data[hsp][ordered[0]
                                                                 ]["HMDM_accession"],
                                                   rgi_data[hsp][ordered[0]
                                                                 ]["model_type"],
                                                   best_snps,
                                                   other_snps,
                                                   "; ".join(rgi_data[hsp][ordered[0]]["HMDM_category"][x][
                                                       "category_hmdm_name"] for x in
                                                       rgi_data[hsp][ordered[0]
                                                                     ]["HMDM_category"]
                                                       if rgi_data[hsp][ordered[0]]["HMDM_category"][x][
                                                       "category_hmdm_class_name"] == 'Drug Class'),
                                                   "; ".join(rgi_data[hsp][ordered[0]]["HMDM_category"][x][
                                                       "category_hmdm_name"] for x in
                                                       rgi_data[hsp][ordered[0]
                                                                     ]["HMDM_category"]
                                                       if rgi_data[hsp][ordered[0]]["HMDM_category"][x][
                                                       "category_hmdm_class_name"] == 'Mechanism'),
                                                   "; ".join(rgi_data[hsp][ordered[0]]["HMDM_category"][x][
                                                       "category_hmdm_name"] for x in
                                                       rgi_data[hsp][ordered[0]
                                                                     ]["HMDM_category"]
                                                       if rgi_data[hsp][ordered[0]]["HMDM_category"][x][
                                                       "category_hmdm_class_name"] == 'Gene Family'),
                                                   "",
                                                   rgi_data[hsp][ordered[0]
                                                                 ]["orf_prot_sequence"],
                                                   rgi_data[hsp][ordered[0]
                                                                 ]["sequence_from_hmdm"],
                                                   format((len(rgi_data[hsp][ordered[0]]["orf_prot_sequence"]) / len(
                                                       rgi_data[hsp][ordered[0]]["sequence_from_hmdm"])) * 100, '.2f'),
                                                   ordered[0],
                                                   rgi_data[hsp][ordered[0]
                                                                 ]["model_id"],
                                                   nudged,
                                                   note,
                                                   format((rgi_data[hsp][ordered[0]]["positives"] / len(
                                                       rgi_data[hsp][ordered[0]]["sequence_from_hmdm"])) * 100, '.2f'),
                                                   "; ".join(rgi_data[hsp][ordered[0]]["HMDM_category"][x][
                                                       "category_hmdm_name"] for x in
                                                       rgi_data[hsp][ordered[0]
                                                                     ]["HMDM_category"]
                                                       if rgi_data[hsp][ordered[0]]["HMDM_category"][x][
                                                       "category_hmdm_class_name"] == 'Drug'),
                                                   "",
                                                   ""
                                                   ]

                            for key, value in match_dict.items():
                                writer.writerow(value)

    def manual():
        h = {}
        h["ORF_ID"] = "Open Reading Frame identifier (internal to DME)"
        h["Contig"] = "Source Sequence"
        h["Start"] = "Start co-ordinate of ORF"
        h["Stop"] = "End co-ordinate of ORF"
        h["Orientation"] = "Strand of ORF"
        h["Cut_Off"] = "DME Detection Paradigm"
        h["Pass_Bitscore"] = "STRICT detection model bitscore value cut-off"
        h["Best_Hit_Bitscore"] = "bitscore value of match to top hit in HMDM"
        h["Best_Hit_HMDM"] = "HMDM term of top hit in HMDM"
        h["Best_Identities"] = "Percent identity of match to top hit in HMDM"
        h["HMDM"] = "HMDM accession of top hit in HMDM"
        h["Model_type"] = "HMDM detection model type"
        h["SNPs_in_Best_Hit_HMDM"] = "Mutations observed in the HMDM term of top hit in HMDM (if applicable)"
        h["Other_SNPs"] = "Mutations observed in HMDM terms of other hits indicated by model id (if applicable)"
        h["Drug Class"] = "HMDM Categorization"
        h["Mechanism"] = "HMDM Categorization"
        h["Gene Family"] = "HMDM Categorization"
        h["Predicted_DNA"] = "ORF predicted nucleotide sequence"
        h["Predicted_Protein"] = "ORF predicted protein sequence"
        h["HMDM_Protein_Sequence"] = "Protein sequence of top hit in HMDM"
        h["Percentage Length of Reference Sequence"] = "Percentage Length of Reference Sequence"
        h["ID"] = "HSP identifier (internal to DME)"
        h["Model_ID"] = "HMDM detection model id"
        h["Percent_Positives"] = "Percent positives value of match to top hit in HMDM"
        h["Drug"] = "Drug(s) metabolized by enzyme"
        h["Strain"] = "Bacterial strain carrying drug-metabolizing enzyme"
        h["Strain_Indentites"] = "Percent identify for top hit to bacterial strain"
        print("\n")
        print("COLUMN", "\t\t\t", "HELP_MESSAGE")
        for i in h:
            print(i, "\t\t\t", h[i])
        print("\n")
