import os

from Bio import Blast, Entrez, SeqIO
from Bio.Blast import NCBIXML # Old BLAST but couldn't find any other way to parse the file (in short time) !NOTE: Come back to check this
# from Bio.KEGG import REST, Enzyme #only parsers and writers for compound, enzyme, and map are currently implemented
from bioservices import *

Blast.email = "dogayasemen.testere@ogr.iu.edu.tr"
Entrez.email = "dogayasemen.testere@ogr.iu.edu.tr"

nuc_sequence = input("Sequence to be searched:")
loc_path = input("Directory to save results: \n")


#  GGTGGTCTGCCTCGCATAAAGCGGTATGAAAATGGATTGAAGCCCGGGCCGTGGATTCTACTCAACTTTCGTCTTTCGAGAAAGACTCCGGGATCCTGAGTATTAAAAAGAAGATCTTTATATAGAGATCTGTTCTATTGTGATCTCTTATTAGGATCGACTCTCTTTGTGGATAAGTCGGATCCGCGAAGTAAGATCAAAAGCTTAAGAAGGATCACTATCTGTGAATGATCGGTGATC
# example sequence

def run_blast_search(nuc_sequence:str, loc_path:str):
    """
    Runs a blastx search for finding the protein corresponding to the sequence
    """

    print("Starting BLASTx search")
    try:
        result_stream = Blast.qblast("blastx","nr", nuc_sequence)
        blast_record = result_stream.read()
        print("Downloading! Bytes:", len(blast_record))

        if loc_path:
            if not os.path.exists(loc_path):
                os.makedirs(loc_path, exist_ok=True)

        save_path = os.path.join(loc_path, 'Blastresult.xml')
        with open(save_path, "wb") as out_handle:
            out_handle.write(blast_record)
        result_stream.close()
        print(f"BLAST record saved to {save_path}")
        return save_path


    except Exception as e:
        print(f"An error occured during search \n {e}")
        return None

def parse_blast_and_fetch_details(blast_xml_path: str):
    """
    Parses the BLAST XML result, gets the top hit, and fetches its details.
    """
    try:
        with open(blast_xml_path) as result_handle:
            blast_record = NCBIXML.read(result_handle)

        if not blast_record.alignments:
            print("No hits found in BLAST result.")
            return

        # Get the top hit
        top_alignment = blast_record.alignments[0]
        top_hit_title = top_alignment.title
        print(f"Top hit found: {top_hit_title}")

        # Extract the protein ID from a common format, e.g., >gi|5524211|gb|AAD44166.1|
        try:
            protein_id = top_hit_title.split('|')[3]
        except IndexError:
            print("Could not parse standard protein ID from title. Using accession.")
            protein_id = top_alignment.accession

        print(f"Fetching details for Protein ID: {protein_id}")

        # Fetch details from NCBI Protein database
        handle = Entrez.efetch(db="protein", id=protein_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        protein_name = record.description
        gene_name = "N/A"
        # Look for the gene name in the record's features
        for feature in record.features:
            if feature.type == "gene" and "gene" in feature.qualifiers:
                gene_name = feature.qualifiers['gene'][0]
                break  # Found it

        print("\n--- Protein Details ---")
        print(f"ID: {protein_id}")
        print(f"Protein Name: {protein_name}")
        print(f"Gene Name: {gene_name}")
        print("-----------------------\n")

        return protein_id


    except FileNotFoundError:
        print(f"Error: Could not find Blastresult.xml at the path '{blast_xml_path}'")
    except Exception as e:
        print(f"An error occurred while parsing or fetching details: {e}")

def run_kegg_search_for_reactions(accs_id:str):
    """
    Checks UniPortDB for the protein with the accession id. !Important: Does not check any other database like UniRef or UniParc
    """
    u = UniProt(verbose=False)
    result = u.search(accs_id)

    print(result)


# --- Main script execution ---

blast_xml_path = os.path.join(loc_path, 'Blastresult.xml')

# Check if BLAST result already exists
if os.path.exists(blast_xml_path):
    print(f"Found existing BLAST result at {blast_xml_path}. Skipping search.")
else:
    # If it doesn't exist, run the search
    print("BLAST result not found. Running search...")
    run_blast_search(nuc_sequence, loc_path)

# Proceed with parsing if the file exists (either pre-existing or newly created)
if os.path.exists(blast_xml_path):
    protein_id = parse_blast_and_fetch_details(blast_xml_path)
    if protein_id:
        run_kegg_search_for_reactions(protein_id)
else:
    print("Could not find or create BLAST result file. Exiting.")
