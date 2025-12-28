import json
import requests
import os
import sys
import time
import argparse

class idMapping:
    def __init__(self, gbk_json, chunkSize=500, out_loc= ".", maxRetries=3, mapping_json="KEGG_mapping.json"):
        self.chunkSize = chunkSize
        self.maxRetries = maxRetries
        self.gbk_json = gbk_json

        self.mapping_json = os.path.join(out_loc, mapping_json)
        self.mapping_data = {}

    BASE_URL = "https://rest.uniprot.org/idmapping/"

    def __mapping_save(self):
        """
        Saves the mapping as a json, this is then passed to other nextflow processes \n
        Structure of json: \n
        {
        "P00561": "ecj:JW0001",
        "P0A8I3": "ecj:JW0005",
        "Q45068": "bsu:BSU18120",
        }
        """
        try:
            with open(self.mapping_json, "w") as file:
                json.dump(self.mapping_data, file, indent=2)
            print(f"---Json(Mapping) Success: File in {self.mapping_json}, total mappings: {len(self.mapping_data)}")
        except Exception as e:
            print(f"---Json(Mapping) Error: {e}")
    
    def __submitting_job(self, chunks, fromDB="UniProtKB_AC-ID", toDB="KEGG"):
        # curl -X POST "https://rest.uniprot.org/idmapping/run" -d "from=UniProtKB_AC-ID" -d "to=KEGG" -d "ids=P00561"
        # {"jobId":"6QC3nJs0lv"}
        payload = {
            "from": fromDB,
            "to": toDB,
            "ids": ",".join(chunks)
        }

        response = requests.post(f"{self.BASE_URL}run", data=payload, timeout=30)
        response.raise_for_status()
        return response.json()["jobId"]
    
    def __checking_job(self, jobID):
        # curl "https://rest.uniprot.org/idmapping/status/6QC3nJs0lv"
        # {"jobStatus":"FINISHED"}
        while True:
            response = requests.get(f"{self.BASE_URL}status/{jobID}", timeout=30)

            if response.status_code == 404:
                time.sleep(5)
                continue

            response.raise_for_status()
            status = response.json()

            if status.get("jobStatus") == "FINISHED" or "results" in status:
                return
            elif status.get("jobStatus") in ("RUNNING", "QUEUED"):
                time.sleep(2)
            else:
                raise RuntimeError(f"Unexpected job status: {status}")
    
    def __getting_results(self, jobID):
        resultsUrl = f"{self.BASE_URL}results/{jobID}"
        page = 1

        while resultsUrl:
            response = requests.get(resultsUrl,timeout=30)
            response.raise_for_status()

            data = response.json()

            for item in data.get("results", []):
                self.mapping_data[item["from"]] = item["to"]

            # Check if there's a next page
            if "next" in response.links:
                resultsUrl = response.links["next"]["url"]
                page += 1
                sys.stdout.write(f"\rFetching results page {page}...")
                sys.stdout.flush()
            else:
                resultsUrl = None
                sys.stdout.write("\n")
                sys.stdout.flush()

    def __run_with_retries(self, bundle):
        for attempt in range(1, self.maxRetries + 1):
            try:
                jobID = self.__submitting_job(bundle)
                self.__checking_job(jobID)
                self.__getting_results(jobID)
                return 
            except (requests.exceptions.RequestException, RuntimeError) as e:
                if attempt == self.maxRetries:
                    raise RuntimeError(f"Mapping failed after {self.maxRetries} attempts: {e}")
                sleep_time = 2 ** (attempt - 1)
                print(f"Retry {attempt}/{self.maxRetries} in {sleep_time}s", file=sys.stderr)
                time.sleep(sleep_time)

    def idMapper(self):
        uniprot_ids = sorted({
            entry["UniProt_ID"]
            for entry in self.gbk_json
            if entry.get("UniProt_ID")
        })

        if not uniprot_ids:
            raise RuntimeError("No UniProt IDs found in input JSON")

        for i in range(0, len(uniprot_ids), self.chunkSize):
            self.__run_with_retries(uniprot_ids[i:i + self.chunkSize])

        self.__mapping_save()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--parse_json", required=True)
    parser.add_argument("--mapping_json", required=True)
    args = parser.parse_args()

    with open(args.parse_json) as f:
        gbk_json = json.load(f)

    mapper = idMapping(
        gbk_json=gbk_json,
        out_loc=".",
        chunkSize=500,
        maxRetries=3,
        mapping_json=args.mapping_json
    )

    mapper.idMapper()

    with open("versions.yml", "w") as f:
        f.write("uniprot_mapping: 1.0\n")
