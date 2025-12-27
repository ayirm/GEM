import os
import json
import sys
import re
from tqdm import tqdm

from .keggParser import keggParser
from .keggRequests import keggRequests

class keggOrganizer:
    def __init__(self, cdsValues, outLoc, cacheFile="KEGG_cache.json"):
        self.cdsValues = cdsValues
        # Ensure outLoc exists and use an absolute path for the cache file so
        # the cache is written to a predictable location regardless of CWD.
        os.makedirs(outLoc, exist_ok=True)
        self.cachePath = os.path.abspath(os.path.join(outLoc, cacheFile))
        self.keggResults = {}

    def __cache_saving(self):
        "Saves the cache into the outLoc/kegg \n Any error would stop the steps that come after this(unless cache already exists)"
        try:
            with open(self.cachePath, "w") as file:
                json.dump(self.keggResults, file, indent=2)
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
                    self.keggResults = json.load(cache)
                self.cacheExists = True
            except Exception as e:
                print(f"---CacheError: Something is wrong with cache at {self.cachePath} \n Error: {e} \n", file=sys.stderr)
                
        return self.cacheExists
    
    def __kegg_dict_making(self, keggID, webRequests, pageParsing):
        results = {
            "geneID": keggID,
            "KO_Terms": [],
            "pathways": [],
            "BRITE": "N/A",
            "BRITE_id": [],
            "reactions": {},
            "compounds": {}
        }

        idPage = webRequests.get_id_page(keggID=keggID)
        if not idPage:
            if getattr(webRequests, 'debug', False):
                print(f"[KEGG DEBUG] No id page returned for {keggID}", file=sys.stderr)
            return results
        else:
            # quick sanity check for expected sections
            if getattr(webRequests, 'debug', False):
                has_orth = 'ORTHOLOGY' in idPage
                has_path = 'PATHWAY' in idPage
                has_brite = 'BRITE' in idPage
                snippet = (idPage[:400] + '...') if len(idPage) > 400 else idPage
                print(f"[KEGG DEBUG] idPage for {keggID}: ORTHOLOGY={has_orth}, PATHWAY={has_path}, BRITE={has_brite}\n{snippet}", file=sys.stderr)
            
        results["KO_Terms"] = pageParsing.find_ko_numbers(idPage)
        results["BRITE"] = pageParsing.find_full_brite(idPage)
        results["BRITE_id"] = pageParsing.find_brite_numbers(idPage)
        
        pathways = pageParsing.find_pathways(idPage)
        results["pathways"] = pathways

        allReactions = set()
        processedMaps = {}

        for pID in pathways:
            mapMatch = re.search(r'(\d{5})$', pID)
            if not mapMatch:
                continue
            mapName = f"map{mapMatch.group(1)}"
            # Turns the path:ecj00260 into map00260

            rID = []
            if mapName in processedMaps:
                rID = processedMaps[mapName]
            else:
                rID = webRequests.get_linked_page(targetDB="reaction", searchID=mapName)
                processedMaps[mapName] = rID

            if rID:
                results["reactions"][pID] = rID
                allReactions.update(rID)

        pbar = tqdm(allReactions, desc="Processing reactions", unit="rxn")
        for rnID in pbar:
            # rnID is like 'rn:R00220' — fetch the reaction entry page, not links
            reactionPage = webRequests.get_id_page(keggID=rnID)
            equations = pageParsing.find_compounds(reactionPage)
            results["compounds"][rnID] = equations
            pbar.set_postfix({"rnID": rnID})
        
        return results
        
    def keggResultsOrganizer(self):
        self.__cache_checking()
        
        if not self.cacheExists:
            unique_kegg_ids = {entry["KEGG_ID"] for entry in self.cdsValues if entry.get("KEGG_ID")}
            
            webRequests = keggRequests()
            pageParsing = keggParser()
            
            for kID in unique_kegg_ids:
                results = self.__kegg_dict_making(kID, webRequests, pageParsing)
                self.keggResults[kID] = results
            
            for entry in self.cdsValues:
                kID = entry.get("KEGG_ID")
                if kID and kID in self.keggResults:
                    kegg_data = self.keggResults[kID]
                    entry.update({
                        "KO_Terms": kegg_data.get("KO_Terms", []),
                        "BRITE": kegg_data.get("BRITE", "N/A"),
                        "BRITE_id": kegg_data.get("BRITE_id", []),
                        "reactions": kegg_data.get("reactions", {}),
                        "compounds": kegg_data.get("compounds", {})
                    })
            
            self.__cache_saving()

        return self.cdsValues
