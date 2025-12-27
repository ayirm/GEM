import os

from gbParser import gbParser
from goTermFinder import goTermFinder
from idMapping import idMapping
from keggFunctions.keggOrganizer import keggOrganizer

def run_pipeline(outDir="SeperatedFiles/results",excelName="annotationResults.xlsx"):
    os.makedirs(outDir, exist_ok=True)

    parserResults = os.path.join(outDir, "parser")
    goResults = os.path.join(outDir, "goTerms")
    idResults = os.path.join(outDir, "idMaps")
    keggResults = os.path.join(outDir, "kegg")

    os.makedirs(parserResults, exist_ok=True)
    os.makedirs(goResults, exist_ok=True)
    os.makedirs(idResults, exist_ok=True)
    os.makedirs(keggResults, exist_ok=True)

    # annotedFile is the output location of prokka annotation
    # prkName is the prefix that will name the prokka results
    annotedFile = "results/prk"
    prkName = "prokkaAnnotes"
    print("---: Parsing of gb file starts \n")
    gbParsing = gbParser(annotedFile, parserResults, prkName)
    cdsValues = gbParsing.gbFileParser()
    print("---: Parsing of gb file ended \n")

    print("---: GO Term finding started \n")
    goFinding = goTermFinder(cdsValues, goResults)
    goValues = goFinding.goFinder()
    print("---: GO Term finding ended \n")

    print("---: ID mapping started \n")
    idMap = idMapping(goValues, idResults)
    idMaps = idMap.idMapper()
    print("---: ID mapping ended \n")

    print("---: KEGG search started \n")
    keggSearch = keggOrganizer(idMaps, keggResults)
    fullResults = keggSearch.keggResultsOrganizer()
    print("---: KEGG search ended \n")

    print(f"---: Excel file created at {outDir}")

if __name__ == "__main__":
    run_pipeline()