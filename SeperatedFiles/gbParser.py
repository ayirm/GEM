import os
import json
import sys

from Bio import SeqIO

class gbParser:
    def __init__(self,annotedFile, outLoc, prkName="prokkaAnnotates", cacheName="parser_cache.json"):
        self.annotedFile = annotedFile
        self.prkName = prkName
        
        self.cachePath = os.path.join(outLoc, cacheName)
        self.parserData = []

    def __cache_saving(self):
        "Saves the cache into the outLoc/parser Any error would stop the steps that come after this(unless cache already exists)"
        try:
            with open(self.cachePath, "w") as file:
                json.dump(self.parserData, file, indent=2)
                # indent allows for more a pretty print.
            print(f"---CacheSuccess: Cache saved at {self.cachePath} \n")
        except Exception as e:
            print(f"---CacheError: Something went wrong with cache saving \n Error: {e} \n", file=sys.stderr)

    def __cache_checking(self):
        self.cacheExists = False
        if os.path.exists(self.cachePath):
            print(f"---CacheSuccess: Already exists at {self.cachePath}, skipping parsing \n")
            try:
                with open(self.cachePath, "r") as cache:
                    self.parserData = json.load(cache)
                self.cacheExists = True
            except Exception as e:
                print(f"---CacheError: Something is wrong with cache at {self.cachePath} \n Error: {e} \n", file=sys.stderr)
                
        return self.cacheExists
            
    # Parsing job

    def __parse_gb_file(self, gbFile):
        """
        Gets the following parts from the provided gbf or gbk file
        - gene="thrA"
        - EC_number="2.7.1.39"
        - inference="similar to AA sequence:UniProtKB:P00561"
        - product="Bifunctional aspartokinase/homoserine dehydrogenase 1"

        Only the cds records that have UniProt id's are saved. Other areas can be empty but UniProtID will always be filled
        """
        totalCds = 0

        for record in SeqIO.parse(gbFile, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    totalCds += 1
                    qualifiers = feature.qualifiers

                    gene = qualifiers.get("gene", [None])[0]
                    protein = qualifiers.get("product", [None])[0]
                    ecNumber = qualifiers.get("EC_number", [None])[0]

                    uniProtID = None
                    for ids in qualifiers.get("inference"):
                        if "UniProtKB" in ids:
                            uniProtID = ids.split("UniProtKB:")[-1]
                            break

                    if uniProtID:
                        self.parserData.append({
                            "Gene":gene,
                            "Protein":protein,
                            "EC_Number":ecNumber,
                            "UniProt_ID":uniProtID
                        })
        sys.stdout.write(f"Parsing {len(self.parserData)}th value out of {totalCds} \n")

    def gbFileParser(self):
        gbFile = os.path.join(self.annotedFile, f"{self.prkName}.gbf")
        if not os.path.exists(gbFile):
            gbFile = os.path.join(self.annotedFile, f"{self.prkName}.gbk")

        self.__cache_checking()
        if not self.cacheExists:
            self.__parse_gb_file(gbFile)
            self.__cache_saving()

        return self.parserData