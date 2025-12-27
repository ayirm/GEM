import os
import json
import sys
import tqdm

from bioservices import QuickGO

class goTermFinder():
    def __init__(self, parsedFile, outLoc, cacheFile="go_cache.json"):
        self.parsedFile = parsedFile

        self.cachePath = os.path.join(outLoc, cacheFile)
        self.goData = {}
        self.go = QuickGO()

    def __cache_saving(self):
        "Saves the cache into the outLoc/goTerms \n Any error would stop the steps that come after this(unless cache already exists)"
        try:
            with open(self.cachePath, "w") as file:
                json.dump(self.goData, file, indent=2)
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
                    self.goData = json.load(cache)
                self.cacheExists = True
            except Exception as e:
                print(f"---CacheError: Something is wrong with cache at {self.cachePath} \n Error: {e} \n", file=sys.stderr)
                
        return self.cacheExists
    
    def __find_go_terms(self, parsedFile):
        """
        First gets the UniProt_ID inside the parsed file then run QuickGO from bioservices to get all of the terms
        As the result, modifies the original parsed file to add the go terms into it and returns the modified file
        """
        # !TODO: go_cache also contains the ascepts like (BP, MF and CC) as "goAspect", so it might be possible to seperate them based on that. I would need to read cache then append the info like this goId, goName, goAspect

        uniProtID = [ids["UniProt_ID"] for ids in parsedFile if ids["UniProt_ID"]]

        for uids in tqdm.tqdm(uniProtID, desc="Fetching GO Terms"):
            try:
                response = self.go.Annotation(
                    geneProductId=uids,
                    includeFields="goName"
                )
                self.goData[uids] = response.get("results", [])
            except Exception as e:
                print(f"---GoError: Error for {uids}: \n {e} \n")

        for entry in parsedFile:
            uid = entry["UniProt_ID"]
            goTerm = self.goData.get(uid, [])
            if goTerm:
                go_terms_list = []
                for g in goTerm:
                    go_id = g.get('goID', g.get('id', 'N/A'))
                    go_name = g.get('goName', g.get('name', 'N/A'))
                    go_aspect = g.get('goAspect', g.get('aspect', 'N/A'))
                    go_terms_list.append(f"{go_id}: {go_name}: {go_aspect}")
                entry["GO_Terms"] = "\n".join(go_terms_list)
            else:
                entry["GO_Terms"] = "N/A"
        
    def goFinder(self):
        self.__cache_checking()
        if not self.cacheExists:
            self.__find_go_terms(self.parsedFile)
            self.__cache_saving()

        return self.parsedFile
