import requests
import time
import sys

class keggRequests:
    def __init__(self, rateLimit=0.35, debug=False):
        self.rateLimit = rateLimit
        self.debug = debug

        self.lastRequestTime = 0

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
            text = response.text
            if self.debug:
                snippet = (text[:500] + '...') if text and len(text) > 500 else text
                print(f"[KEGG DEBUG] GET {reqUrl} -> len={len(text) if text else 0}\n{snippet}", file=sys.stderr)
            return text
        except requests.exceptions.RequestException as e:
            print(f"---KeggRequestError: Error at '{endUrl}': {e}", file=sys.stderr)
            return None
        
    def get_id_page(self, keggID):
        endURL = f"/get/{keggID}"
        return self.__make_request(endURL)
    
    def get_linked_page(self, targetDB, searchID):
        endURL = f"/link/{targetDB}/{searchID}"
        results = self.__make_request(endURL)

        if not results:
            return []
        
        links = []
        for line in results.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) == 2:
                links.append(parts[1])

        return links
