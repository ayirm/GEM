import json
import requests
import os
import sys
import time

class idMapping:
    def __init__(self, cdsValues, outLoc, chunkSize=500, maxRetries=3, cacheFile="idMapping.json"):
        self.chunkSize = chunkSize
        self.maxRetries = maxRetries
        self.cdsValues = cdsValues

        self.cachePath = os.path.join(outLoc, cacheFile)
        self.mappingData = {}

    BASE_URL = "https://rest.uniprot.org/idmapping/"

    def __cache_saving(self):
        "Saves the cache into the outLoc/idMaps \n Any error would stop the steps that come after this(unless cache already exists)"
        try:
            with open(self.cachePath, "w") as file:
                json.dump(self.mappingData, file, indent=2)
                # indent allows for more a pretty print.
            print(f"---CacheSuccess: Cache saved at {self.cachePath}")
        except Exception as e:
            print(f"---CacheError: Something went wrong with cache saving \n Error: {e} \n", file=sys.stderr)

    def __cache_checking(self):
        self.cacheExists = False
        if os.path.exists(self.cachePath):
            print(f"---CacheSuccess: Already exists at {self.cachePath}, skipping parsing")
            try:
                with open(self.cachePath, "r") as cache:
                    self.mappingData = json.load(cache)
                self.cacheExists = True
            except Exception as e:
                print(f"---CacheError: Something is wrong with cache at {self.cachePath} \n Error: {e} \n", file=sys.stderr)
                
        return self.cacheExists
    
    def __submitting_job(self, chunks, fromDB="UniProtKB_AC-ID", toDB="KEGG"):
        # curl -X POST "https://rest.uniprot.org/idmapping/run" -d "from=UniProtKB_AC-ID" -d "to=KEGG" -d "ids=P00561"
        # {"jobId":"6QC3nJs0lv"}

        bundle = {
            "from" : fromDB,
            "to" : toDB,
            "ids" : ",".join(chunks)
        }

        response = requests.post(f"{self.BASE_URL}run", data=bundle)
        response.raise_for_status()
        jobID = response.json()["jobId"]

        return jobID
    
    def __checking_job(self, jobID):
        # curl "https://rest.uniprot.org/idmapping/status/6QC3nJs0lv"
        # {"jobStatus":"FINISHED"}
        jobStatus = False

        while not jobStatus:
            statusUrl = f"{self.BASE_URL}status/{jobID}"
            response = requests.get(statusUrl)

            if response.status_code == 404:
                sys.stdout.write("Polling job status... (job not ready, waiting)")
                sys.stdout.flush()
                time.sleep(5)
                continue

            response.raise_for_status()
            statusData = response.json()

            if "results" in statusData or statusData.get("jobStatus") == "FINISHED":
                jobStatus = True
                return jobStatus
            elif statusData.get("jobStatus") in ["RUNNING", "QUEUED"]:
                sys.stdout.write("Polling job status...")
                sys.stdout.flush()
                time.sleep(2)
            else:
                raise RuntimeError(f"Error in jobStatus: \n {json.dumps(statusData, indent=2)}")
    
    def __getting_results(self, jobID):
        resultsUrl = f"{self.BASE_URL}results/{jobID}"
        pageCount = 1

        while resultsUrl:
            response = requests.get(resultsUrl)
            response.raise_for_status()
            
            for item in response.json().get("results", []):
                uid = item["from"]
                if uid not in self.mappingData:
                    self.mappingData[uid] = item["to"]

            # Check if there's a next page
            if response.links and "next" in response.links:
                resultsUrl = response.links["next"]["url"]
                pageCount += 1
                sys.stdout.write(f"\rFetching page {pageCount}...")
                sys.stdout.flush()
            else:
                resultsUrl = None
                sys.stdout.write("\n")
                sys.stdout.flush()

    def __run_with_retries(self, bundle):
    # !TODO: It repepats the whole thing 3 times even it is finished
        for attempt in range(self.maxRetries):
            try:
                jobID = self.__submitting_job(chunks=bundle)
                jobStatus = self.__checking_job(jobID)
                if jobStatus:
                    self.__getting_results(jobID)
                    self.__cache_saving()
            except (requests.exceptions.RequestException, RuntimeError) as e:
                print(f"Attempt {attempt+1}/{self.maxRetries} failed: {e}", file=sys.stderr)
                if attempt +1 == self.maxRetries:
                    print("---MappingFailure: Giving up on chunk")
                else:
                    delay = 2 ** attempt
                    print(f"Retry in {delay} seconds")
                    time.sleep(delay)

    def idMapper(self):
        uniprotID = [entry.get("UniProt_ID") for entry in self.cdsValues if entry["UniProt_ID"]]

        self.__cache_checking()
        if not self.cacheExists:
            for i in range(0, len(uniprotID), self.chunkSize):
                bundle = uniprotID[i:i+self.chunkSize]
                # Running of the other functions happen in __run_with_retries, bcs if i add it here it becomes too much nested
                self.__run_with_retries(bundle)

                for entry in self.cdsValues:
                    uid = entry.get("UniProt_ID")
                    entry["KEGG_ID"] = self.mappingData.get(uid)
        
        return self.cdsValues


            
