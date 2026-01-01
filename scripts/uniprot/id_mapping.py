#!/usr/bin/env python3

import json
import requests
import sys
import time
import argparse

from datetime import datetime, timezone

BASE_URL = "https://rest.uniprot.org/idmapping/"

def save_mapping(mapping_data, out_path):
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
        with open(out_path, "w") as fh:
            json.dump(mapping_data, fh, indent=2)
        print(f"---Json(Mapping) Success: File in {out_path}, total mappings: {len(mapping_data)}")
    except Exception as e:
        raise RuntimeError(f"Failed to write mapping JSON: {e}")

def submit_job(chunks, from_db="UniProtKB_AC-ID", to_db="KEGG"):
    # curl -X POST "https://rest.uniprot.org/idmapping/run" -d "from=UniProtKB_AC-ID" -d "to=KEGG" -d "ids=P00561"
    # {"jobId":"6QC3nJs0lv"}
    payload = {
        "from": from_db,
        "to": to_db,
        "ids": ",".join(chunks),
    }

    response = requests.post(f"{BASE_URL}run", data=payload, timeout=30)
    response.raise_for_status()
    return response.json()["jobId"]

def wait_for_job(job_id):
    # curl "https://rest.uniprot.org/idmapping/status/6QC3nJs0lv"
    # {"jobStatus":"FINISHED"}
    while True:
        response = requests.get(f"{BASE_URL}status/{job_id}", timeout=30)

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


def fetch_results(job_id, mapping_data):
    results_url = f"{BASE_URL}results/{job_id}"
    page = 1

    while results_url:
        response = requests.get(results_url, timeout=30)
        response.raise_for_status()
        data = response.json()

        for item in data.get("results", []):
            mapping_data[item["from"]] = item["to"]

        # Check if there's a next page
        if "next" in response.links:
            results_url = response.links["next"]["url"]
            page += 1
            sys.stdout.write(f"\rFetching results page {page}...")
            sys.stdout.flush()
        else:
            results_url = None
            sys.stdout.write("\n")
            sys.stdout.flush()


def run_with_retries(bundle, mapping_data, max_retries):
    for attempt in range(1, max_retries + 1):
        try:
            job_id = submit_job(bundle)
            wait_for_job(job_id)
            fetch_results(job_id, mapping_data)
            return
        except (requests.exceptions.RequestException, RuntimeError) as e:
            if attempt == max_retries:
                raise RuntimeError(
                    f"Mapping failed after {max_retries} attempts: {e}"
                )
            sleep_time = 2 ** (attempt - 1)
            print(
                f"Retry {attempt}/{max_retries} in {sleep_time}s",
                file=sys.stderr,
            )
            time.sleep(sleep_time)


def run_id_mapping(gbk_json, out_path, chunk_size=500, max_retries=3):
    uniprot_ids = sorted({
        entry["UniProt_ID"]
        for entry in gbk_json
        if entry.get("UniProt_ID")
    })

    if not uniprot_ids:
        raise RuntimeError("No UniProt IDs found in input JSON")

    mapping_data = {}

    for i in range(0, len(uniprot_ids), chunk_size):
        run_with_retries(
            uniprot_ids[i:i + chunk_size],
            mapping_data,
            max_retries,
        )

    save_mapping(mapping_data, out_path)

def write_versions():
    versions = {
        "json_merging": {
            "python": sys.version.split()[0],
            "timestamp": datetime.now(timezone.utc).isoformat()
        }
    }

    with open("versions.yml", "w") as f:
        json.dump(versions, f, indent=2)

def main():
    parser = argparse.ArgumentParser(
        description="Map UniProt IDs to KEGG IDs using UniProt REST API"
    )
    parser.add_argument("--parse_json", required=True)
    parser.add_argument("--mapping_json", required=True)
    parser.add_argument("--chunk_size", type=int, default=500)
    parser.add_argument("--max_retries",type=int, default=3)
    args = parser.parse_args()

    with open(args.parse_json) as fh:
        gbk_json = json.load(fh)

    run_id_mapping(
        gbk_json=gbk_json,
        out_path=args.mapping_json,
        chunk_size=args.chunk_size,
        max_retries=args.max_retries,
    )

    write_versions()


if __name__ == "__main__":
    main()