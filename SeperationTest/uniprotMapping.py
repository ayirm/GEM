import requests, json, sys, time, os

class UniprotMapping:
    
    BASE_URL = "https://rest.uniprot.org/idmapping/"
    
    def __init__(self, chunkSize=500, maxRetries=3, cacheFile="UniToKegg.json"):
        self.chunkSize = chunkSize
        self.maxRetries = maxRetries
        self.mappingDict = self.__load_cache()
        self.cacheFile = cacheFile

    def __load_cache(self):
        if os.path.exists(self.cacheFile):
            try:
                with open(self.cacheFile, "r") as f:
                    data = json.load(f)
                print(f"[---Cache: Loaded {len(data)} mappings from {self.cacheFile}")
                return data
            except json.JSONDecodeError:
                print(f"---Cache: cache file {self.cacheFile} is broken. Starting new")
        return {}

    def __save_cache(self):
        with open(self.cacheFile, "w") as f:
            json.dump(self.mappingDict, f, indent=2)
        print(f"[Cache] Saved {len(self.mappingDict)} total mappings to {self.cacheFile}")

    def __submit_chunks(self, chunks, from_db="UniProtKB_AC-ID", to_db="KEGG"):
        
        # curl -X POST "https://rest.uniprot.org/idmapping/run" -d "from=UniProtKB_AC-ID" -d "to=KEGG" -d "ids=P00561"
        bundle = {
            "from" : from_db,
            "to" : to_db,
            "ids" : ",".join(chunks)
        }
        # {"jobId":"6QC3nJs0lv"}
        response = requests.post(f"{self.BASE_URL}/run", data=bundle)
        response.raise_for_status()
        jobID = response.json()["jobID"]

        return jobID
    
    def __check_job_status(self, jobID):
        
        spinner_idx = 0
        job_status = False

        while not job_status: # It would work if the job_status is False

            # curl "https://rest.uniprot.org/idmapping/status/6QC3nJs0lv"
            status_url = f"{self.BASE_URL}/status/{jobID}"
            response = requests.get(status_url)

            if response.status_code == 404:
                sys.stdout.write(f"Polling job status... (job not ready, waiting) {'/-|'[spinner_idx % 4]}")
                sys.stdout.flush()
                spinner_idx += 1
                time.sleep(5)
                continue

            response.raise_for_status()
            status_data = response.json()

            # {"jobStatus":"FINISHED"}
            if "results" in status_data or status_data.get("jobStatus") == "FINISHED":
                job_status = True
                return job_status
            elif status_data.get("jobStatus") in ["RUNNING", "QUEUED"]:
                sys.stdout.write(f"\rPolling job status... {'/-|'[spinner_idx % 4]}")
                sys.stdout.flush()
                spinner_idx += 1
                time.sleep(2)
            else:
                raise RuntimeError(f"Unexpected job status: \n {json.dumps(status_data, indent=2)}")
            
    def __read_the_results(self, jobID):

        results_url = f"{self.BASE_URL}/results/{jobID}"
        page_count = 1

        while results_url:
            response = requests.get(results_url)

            for item in response.json().get("results", []):
                u_id = item["from"]
                if u_id not in self.mappingDict:
                    self.mappingDict[u_id] = item["to"]

                if "next" in response.links:
                    results_url = response.links["next"]["url"]
                    sys.stdout.write(f"\rFetching page {page_count + 1}...")
                    sys.stdout.flush()
                    page_count += 1
                else:
                    results_url = None

    def organize_Mapping(self, uniProt_ids):
        """
        
        """
        new_ids = [uid for uid in uniProt_ids if uid not in self.mappingDict]

        if not new_ids:
            print("---Cache: No new id's to check")
            return self.mappingDict

        for i in range(0,len(new_ids), self.chunkSize):
            bundle = new_ids[i:i+self.chunkSize]
            if not bundle:
                continue

            print(f"Submitting bundle {i//self.chunk_size + 1}/{(len(new_ids) - 1 )//self.chunk_size} for mapping")

            for attempt in range(self.maxRetries):
                try:
                    job_id = self.__submit_chunks(chunks=bundle)
                    job_status = self.__check_job_status(jobID=job_id)
                    if job_status:
                        self.__read_the_results(jobID=job_id)
                        self.__save_cache()
                    else:
                        print(f"Job {job_id} failed or returned unexpected status. Skipping this chunk.", file=sys.stderr)
                    break
                except (requests.exceptions.RequestException, RuntimeError) as e:
                    print(f"\nAttempt {attempt + 1}/{self.max_retries} failed: {e}", file=sys.stderr)
                    if attempt + 1 == self.max_retries:
                        print(f"--- Giving up on this chunk. ---", file=sys.stderr)
                    else:
                        delay = 2 ** attempt
                        print(f"--- Retrying in {delay} seconds... ---")
                        time.sleep(delay)
            
            return self.mappingDict