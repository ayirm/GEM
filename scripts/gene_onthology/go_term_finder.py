import os
import json
import tqdm
import argparse

from bioservices import QuickGO

class goTermFinder():
    def __init__(self, gbk_json, out_loc=".", go_json="go_terms.json"):
        self.gbk_json = gbk_json

        self.go_json = os.path.join(out_loc, go_json)
        self.go_data = {}
        self.go = QuickGO()

    def __term_save(self):
        """
        Saves the term finding results as a json. Structure is this: \n
        {
            "P00561": [
                {
                    "id": "UniProtKB:P00561!612814768",
                    "geneProductId": "UniProtKB:P00561",
                    "qualifier": "acts_upstream_of_or_within",
                    "goId": "GO:0009089",
                    "goName": "L-lysine biosynthetic process via diaminopimelate",
                    "goEvidence": "IDA",
                    "goAspect": "biological_process",
                    "evidenceCode": "ECO:0000314",
                    "reference": "PMID:11368768",
                    "withFrom": null,
                    "taxonId": 83333,
                    "taxonName": null,
                    "assignedBy": "EcoCyc",
                    "extensions": null,
                    "targetSets": null,
                    "symbol": "thrA",
                    "date": "20070802",
                    "synonyms": null,
                    "name": null
                },
        """
        try:
            with open(self.go_json, "w") as file:
                json.dump(self.go_data, file, indent=2)
            print("---Json(GoTerms) Success")
        except Exception as e:
            print(f"---Json(GoTerms) Error: \n {e}")

    def __find_go_terms(self, json_file):
        uniProtIDs = {
            entry["UniProt_ID"]
            for entry in json_file
            if entry.get("UniProt_ID")
        }

        for uid in tqdm.tqdm(uniProtIDs, desc="Fetching GO Terms"):
            try:
                response = self.go.Annotation(
                    geneProductId=uid,
                    includeFields="goName,goAspect"
                )
                self.go_data[uid] = response.get("results", [])
            except Exception as e:
                print(f"---GoError: {uid}\n{e}")
                self.go_data[uid] = []

    def goFinder(self):
        self.__find_go_terms(self.gbk_json)
        self.__term_save()
        return self.gbk_json

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_json", required=True)
    parser.add_argument("--go_json", required=True)

    args = parser.parse_args()

    with open(args.input_json) as f:
        gbk_json = json.load(f)

    termFinder = goTermFinder(
        gbk_json=gbk_json,
        out_loc=".",
        go_json=args.go_json
    )

    termFinder.goFinder()

    with open("versions.yml", "w") as f:
        f.write("go_term_finding: 1.0\n")
