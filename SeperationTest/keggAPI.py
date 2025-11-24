import os
import json
import time
import requests
import sys
import re
from tqdm import tqdm

class KeggAPI:
    def __init__(self, cdsValues, cacheFile="Kegg_Cache.json", outDir="results"):
        self.cacheFile = cacheFile
        self.rateLimit = 0.35
        self.lastRequestTime = 0
        self.cdsValues = cdsValues

        self.cachePath = os.path.join(outDir, self.cacheFile)
        self.keggResults = self.__load_cache()

    def __load_cache(self):
        if os.path.exists(self.cachePath):
            try:
                with open(self.cachePath, "r") as f:
                    data = json.load(f)
                print(f"[---Cache: Loaded {len(data)} pathways from {self.cachePath}")
                return data
            except json.JSONDecodeError:
                print(f"---Cache: {self.cachePath} is broken. Starting new")
        return {}

    def __save_cache(self):
        with open(self.cachePath, "w") as f:
            json.dump(self.keggResults, f, indent=2)
            print("---Cache: Saved the cache")

    BASE_URL = "https://rest.kegg.jp"

    def __make_request(self, endUrl):
        elapsed = time.time() - self.lastRequestTime
        if elapsed < self.rateLimit:
            sleep_time = self.rateLimit - elapsed
            time.sleep(sleep_time)
        # Basically checks if the time since the last request is smaller than 0.35
        # if it is smaller, sleeps for that amount
        # ex: time elapsed is 0,25 then it sleeps for 0,10. 0,20 for 0,15 etc.

        reqUrl = f"{self.BASE_URL}{endUrl}"
        try:
            response = requests.get(reqUrl)
            self.lastRequestTime = time.time()
            response.raise_for_status()
            return response.text
        except requests.exceptions.RequestException as e:
            print(f"\n[KEGGClient Error] Error querying URL '{endUrl}': {e}", file=sys.stderr)
            return None

    def __parse_page_entries(self, entryLocation, sectionName):
        """
        Searches the entry inside the provided page, returns the part of the searched Entry
        """
        if not entryLocation:
            return []
        
        lines = []
        entryExists = False
        # All KEGG entries align section names. Content starts 12 spaces in.
        content_indent = "            " # 12 spaces

        for line in entryLocation.split("\n"):
            if line.startswith(sectionName):
                entryExists = True
                lines.append(line[len(sectionName):].strip())
            elif entryExists and line.startswith(content_indent):
                lines.append(line.strip())
            elif entryExists and not line.startswith(" "):
                break
        
        return lines
    
    def __get_id_page(self, keggID):
        endURL = f"/get/{keggID}"
        return self.__make_request(endUrl=endURL)

    def __get_linked_page(self, targetDB, searchID):
        endURL = f"/link/{targetDB}/{searchID}"
        results = self.__make_request(endUrl=endURL)

        if not results:
            return []
        
        links = []
        for line in results.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) == 2:
                links.append(parts[1])

        return links

    # --- Parsing(↓) 

    def __find_ko_numbers(self, pageResult):
        """
        ORTHOLOGY   K12524  bifunctional aspartokinase / homoserine dehydrogenase 1 [EC:2.7.2.4 1.1.1.3]
        ^ Finds this in the page and gets the K12524 part.
        """
        koLines = self.__parse_page_entries(pageResult, "ORTHOLOGY")
        koID = []

        # Regex to find one or more K numbers
        koRegex = re.compile(r'(K\d{5,})')

        for line in koLines:
            matches = koRegex.findall(line)
            for ko in matches:
                koID.append(f"KO:{ko}")

        return koID

    def __find_full_brite(self, pageResult):
        """
        BRITE       KEGG Orthology (KO) [BR:ecj00001]
             09100 Metabolism
              09105 Amino acid metabolism
               00260 Glycine, serine and threonine metabolism
                JW0001 (thrA)
               00270 Cysteine and methionine metabolism
                JW0001 (thrA)
               00300 Lysine biosynthesis
                JW0001 (thrA)
              09110 Biosynthesis of other secondary metabolites
               00261 Monobactam biosynthesis
                JW0001 (thrA)
            Enzymes [BR:ecj01000]
             1. Oxidoreductases
              1.1  Acting on the CH-OH group of donors
               1.1.1  With NAD+ or NADP+ as acceptor
                1.1.1.3  homoserine dehydrogenase
                 JW0001 (thrA)
             2. Transferases
              2.7  Transferring phosphorus-containing groups
               2.7.2  Phosphotransferases with a carboxy group as acceptor
                2.7.2.4  aspartate kinase
                 JW0001 (thrA)
        ^ Finds this then just returns it
        """
        briteResults = self.__parse_page_entries(pageResult, "BRITE")
        if not briteResults:
            return "No BRITE section"
        
        return "\n".join(briteResults)

    def __find_brite_numbers(self, pageResult):
        briteResults = self.__parse_page_entries(pageResult, "BRITE")
        if not briteResults:
            return []
        
        briteID = []
        briteRegex = re.compile(r"\[(BR:[a-z0-9]+)\]")

        for line in briteResults:
            matches = briteRegex.findall(line)
            for brID in matches:
                briteID.append(brID)

        return briteID

    def __find_pathways(self, pageResult):
        """
        PATHWAY     ecj00260  Glycine, serine and threonine metabolism
            ecj00261  Monobactam biosynthesis
            ecj00270  Cysteine and methionine metabolism
            ecj00300  Lysine biosynthesis
            ecj01100  Metabolic pathways
            ecj01110  Biosynthesis of secondary metabolites
            ecj01120  Microbial metabolism in diverse environments
            ecj01230  Biosynthesis of amino acids
        ^ Finds this and then seperates the ecj... parts
        
        Then returns the ecj part. No convertion to map happens here
        """
        pathwayResults = self.__parse_page_entries(pageResult, "PATHWAY")
        pathwayID = []
        pathwayRegex = re.compile(r'([a-z]{3,5}\d{5})')

        for line in pathwayResults:
            # Find the first ID on the line
            match = pathwayRegex.match(line)
            if match:
                pathwayID.append(f"path:{match.group(1)}") 
        # Add 'path:' prefix for consistency, normally just map00260 works but other ones have br:, ko: etc.
        # It was purely an aesthetic choice

        return pathwayID

    def __find_compounds(self, reactionID):
        """
        Makes a get/ request for the reaction page then gets
        DEFINITION  ATP + L-Aspartate <=> ADP + 4-Phospho-L-aspartate
        EQUATION    C00002 + C00049 <=> C00008 + C03082
        parts

        Then sends them back as equation(DEFINITION) and cmpd_id(EQUATION)
        """

        reactionPage = self.__get_id_page(reactionID)
        if not reactionPage:
            return {
                "equation" : "N/A(ERROR)",
                "cmpdID" : "N/A(ERROR)"
            }
        
        readableEQ = " ".join(self.__parse_page_entries(reactionPage, "DEFINITION"))
        cmpdEQ = " ".join(self.__parse_page_entries(reactionPage, "EQUATION"))

        return {
            "equation" : readableEQ if readableEQ else "N/A",
            "cmpdID" : cmpdEQ if cmpdEQ else "N/A"
        }

    def organizeKeggResults(self, keggID):
        if keggID in self.keggResults:
            print(f"---Cache: Using cached KEGG data for {keggID}")
            return self.keggResults[keggID]

        print(f"---Cache: Fetching KEGG data for {keggID}")
        idPage = self.__get_id_page(keggID=keggID)
        if not idPage:
            print(f"No page for {keggID}. Cancelling Kegg requests")
            return None
        
        resultsDict = {
            "geneID" : keggID,
            "KO_Terms" : [],
            "pathways" : [],
            "BRITE" : "N/A",
            "BRITE_id" : [],
            "reactions" : {}, # map_id : list of reactions
            "compunds" : {} # reaction_id : {equation, cmpdID}
        }

        resultsDict["KO_Terms"] = self.__find_ko_numbers(idPage)
        # print("---Debugging: Finished KO Terms")

        resultsDict["BRITE"] = self.__find_full_brite(idPage)
        # print("---Debuggging: Finished BRITE hieacrhy")

        resultsDict["BRITE_id"] = self.__find_brite_numbers(idPage)
        # print("---Debugging: Found brite id's")

        pathways = self.__find_pathways(idPage)
        resultsDict["pathways"] = pathways
        # print("---Debugging: Found the pathways")

        allReactions = set()
        processedMaps = {}

        for pID in pathways:
            mapMatch = re.search(r'(\d{5})$', pID)
            mapID = f"map{mapMatch.group(1)}"
            # Turns the path:ecj00260 into map00260

            reactionID = []
            if mapID in processedMaps:
                reactionID = processedMaps[mapID]
            else:
                reactionID = self.__get_linked_page("reaction", mapID)
                processedMaps[mapID] = reactionID

            if reactionID:
                resultsDict["reactions"][pID] = reactionID
                allReactions.update(reactionID)
        # print("---Debugging: Reactions per mapID found")

        #for i, rnID in enumerate(allReactions, 1):
        #    sys.stdout.write(f"Processing reaction {i}/{len(allReactions)} ({rnID})")
        #    sys.stdout.flush()
        #    equations = self.__find_compounds(rnID)
        #    resultsDict["compunds"] = equations

        pbar = tqdm(allReactions, desc="Processing reactions", unit="rxn")
        for rnID in pbar:
            equations = self.__find_compounds(rnID)
            resultsDict["compunds"] = equations
            pbar.set_postfix({"rnID": rnID})
        print("---Debugging: Finished compound info")

        self.__save_cache()

        return resultsDict