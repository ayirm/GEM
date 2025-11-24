import sys
import json
import os
import time

from Bio import SeqIO
from bioservices import QuickGO
from tqdm import tqdm

class UniprotGO:
    def __init__(self, annPath, prkName="prokkaAnnotes", cacheFile="go_cache.json"):
        self.annPath = annPath
        self.prkName = prkName
        self.cacheFile = cacheFile

        self.cachePath = os.path.join(self.annPath, self.cacheFile)
        self.go_data = {}
        self.go = QuickGO()

    def __save_cache(self):
        try:
            with open(self.cachePath, "w") as cache:
                json.dump(self.go_data, cache, indent=2)
            print(f"Saved the cache in {self.cachePath}")
        except Exception as e:
            print(f"An error happened \n {e}", file=sys.stderr)

    def __check_cache(self):
        self.cache_exists = False
        if os.path.exists(self.cachePath):
            print(f"Found cache at {self.cachePath}")
            try:
                with open(self.cachePath, "r") as cache:
                    self.go_data = json.load(cache)
                print("\n Loaded from cache")
                self.cache_exists = True
                return self.cache_exists
            except Exception as e:
                print(f"Something happened: \n {e}")
                self.go_data = {}
                return self.cache_exists
        else:
            print("\n No cache found")
            self.go_data = {}
            return self.cache_exists


    def __parse_gb_file(self, gbFile):
        """
        Gets the following parts from the provided gbf or gbk file
        - gene="thrA"
        - EC_number="2.7.1.39"
        - inference="similar to AA sequence:UniProtKB:P00561"
        - product="Bifunctional aspartokinase/homoserine dehydrogenase 1"

        Only the cds records that have inference id's are saved since they are mostly enzymes
        """

        cds_values = []
        total_cds = 0

        for record in SeqIO.parse(gbFile, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    total_cds += 1
                    
                    qualifiers = feature.qualifiers
                    gene = qualifiers.get("gene", [None])[0]
                    protein = qualifiers.get("product", [None])[0]
                    ec_number = qualifiers.get("EC_number", [None])[0]

                    uni_inference = None
                    for inf in qualifiers.get("inference"):
                        if "UniProtKB" in inf:
                            uni_inference = inf.split("UniProtKB:")[-1]
                            break

                    if uni_inference:
                        cds_values.append({
                            "gene":gene,
                            "UniProt_ID":uni_inference,
                            "Protein": protein,
                            "EC_number": ec_number
                        })
        sys.stdout.write(f"Parsing {len(cds_values)}/{total_cds} valid cds records")

        return cds_values

    def __search_GO(self, cds_values):
        """
        Strips the UniProt_ID's out of cds_values then checks if there are new id's to write
        This check is done by cache. After the first run it saves a cache of the GO results
        """
        uniprotID = [entry["UniProt_ID"] for entry in cds_values if entry["UniProt_ID"]]

        new_ids = [uid for uid in uniprotID if uid not in self.go_data]
        if not new_ids:
            print("---Cache: No new id's to check")
            return

        
        for uid in tqdm(new_ids, desc="Fetching GO terms"):
            try:
                response = self.go.Annotation(
                    geneProductId=uid, 
                    includeFields="goName"
                )
                self.go_data[uid] = response.get("results", [])
            except Exception as e:
                print(f"Something happened for {uid}: \n {e}")
                # self.go_data[uid] = []
            time.sleep(0.1)

        sys.stdout.write(f"Searching {len(uid)}/{len(new_ids)} for the GO terms")
        sys.stdout.flush()

        self.__save_cache()

    def organize_GO(self):
        """
        Checks the genbank file(.gbf or .gbk) then gets the uniprot id's from the file
        Which is then checked against the GO site with QuickGO
        """
        gb_file = os.path.join(self.annPath, f"{self.prkName}.gbf")
        if not os.path.exists(gb_file):
            gb_file = os.path.join(self.annPath, f"{self.prkName}.gbk")

        cds_values = self.__parse_gb_file(gbFile=gb_file)
        
        self.__check_cache()
        if self.cache_exists:
            pass
        else:
            self.__search_GO(cds_values=cds_values)

        for entry in cds_values:
            uid = entry["UniProt_ID"]
            go_terms = self.go_data.get(uid, [])
            #print(f"DEBUG: {uid} GO data: {go_terms}")
            entry["GO_Terms"] = "\n".join([f"{g.get('goID', 'N/A')}: {g["goName"]}" for g in go_terms]) if go_terms else "N/A"

        print("GO Term search is complete")
        return cds_values