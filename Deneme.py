# bioinformatics_pipeline.py
# Description: This script takes a nucleotide sequence, identifies the corresponding gene and protein,
# finds associated metabolic reactions, and outputs the information to an Excel file.

import pandas as pd
from Bio.Blast import NCBIXML, NCBIWWW
from Bio import Entrez
import requests
import io
import time

# --- CONFIGURATION ---
# NCBI requires you to identify yourself. Please replace with your own email.
Entrez.email = "dogayasemen.testere@ogr.iu.edu.tr"
# Input nucleotide sequence (Example: E. coli bglA gene)
NUCLEOTIDE_SEQUENCE = "ATGACACAAAAACAAAAAACGTTTTAACGTTAACTTTGACGTTTTATTTATCGGTTTTACGCCCAGCGGCGCAGAAATGCTGTTGATATTTGGCGCGGTGCGATGGGTAACGGTAATATTGTCGGTGGTCGTTCTGCGTTTGTTGCTGTTGAACAGGGCGATCAGTTAGCGTTGGCGCATTCCGTTGTCGTACTGGAAGCGGATGGCGGTAAAGCGGGTGGCGCGGTAGGTGAAGTTGCTGTTGTTGAAGGTGAAGGTAAAGGTCTGGGTGCGGCGGCAGGTGCGGAAGTTGGTGCGATTGGTGCGAAAGGTGCGGGTGCGAAAGGTCTGGCAGGTGCGAAAGGTCTGGGTGCGGGTGAAGGTGAAGTTGCTGTTGTTGAAGGTGAAGGTAAAGGTCTGGGTGCGGCGGCAGGTGCGGAAGTTGGTGCGATTGGTGCGAAAGGTGCGGGTGCGAAAGGTCTGGCAGGTGCGAAAGGTCTGGGTGCGGGTGAAGGTGAAGTTGCTGTTGTTGAAGGTGAAGGTAAAGGTCTGGGTGCGGCGGCAGGTGCGGAAGTTGGTGCGATTGGTGCGAAAGGTGCGGGTGCGAAAGGTCTGGCAGGTGCGAAAGGTCTGGGTGCGGGTGAAGTGAA"
OUTPUT_FILENAME = "biochemical_analysis_results.xlsx"


def run_blast_search(sequence: str):
    """
    Performs a BLAST search for a given nucleotide sequence against the 'nt' database.
    Returns the top BLAST record.
    """
    print("Performing BLAST search... (This may take a moment)")
    try:
        # Use qblast to perform the search
        result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
        blast_records = NCBIXML.parse(result_handle)

        # Get the first record (highest significance)
        first_record = next(blast_records, None)
        return first_record

    except Exception as e:
        print(f"An error occurred during BLAST search: {e}")
        return None


def get_gene_and_protein_info(blast_record):
    """
    Parses a BLAST record to extract gene and protein information using Entrez.
    """
    if not blast_record or not blast_record.alignments:
        print("No significant BLAST hits found.")
        return None, None

    # Get the accession number of the top hit
    top_alignment = blast_record.alignments[0]
    accession_id = top_alignment.accession
    print(f"Top hit found: {accession_id}")

    print("Fetching gene and protein details from Entrez...")
    try:
        # Fetch the GenBank record for the accession ID
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text")
        record_text = handle.read()
        handle.close()

        gene_name = None
        protein_name = None
        protein_id = None

        # Parse the GenBank record manually to find relevant features
        for line in record_text.split('\n'):
            if "/gene=" in line:
                gene_name = line.split('"')[1]
            if "/product=" in line:
                protein_name = line.split('"')[1]
            if "/protein_id=" in line:
                protein_id = line.split('"')[1]
        
        print(f"Found Gene: {gene_name}, Protein: {protein_name}, Protein ID: {protein_id}")
        return gene_name, protein_name

    except Exception as e:
        print(f"An error occurred fetching data from Entrez: {e}")
        return None, None


def get_kegg_reactions(protein_name: str):
    """
    Searches the KEGG database for a protein name to find associated reactions.
    """
    if not protein_name:
        return []

    print(f"Searching KEGG for protein: '{protein_name}'...")
    reactions_data = []
    
    # --- Find KEGG gene entry for the protein ---
    # We search for the protein name in the E. coli (eco) organism database
    search_url = f"http://rest.kegg.jp/find/genes/{protein_name}"
    response = requests.get(search_url)

    if response.status_code != 200 or not response.text:
        print(f"Could not find KEGG entry for protein '{protein_name}'.")
        return []
        
    # Get the first gene entry ID
    gene_id = response.text.split('\t')[0]
    print(f"Found KEGG Gene ID: {gene_id}")
    
    # --- Get detailed info for the gene to find associated pathways/reactions ---
    get_url = f"http://rest.kegg.jp/get/{gene_id}"
    response = requests.get(get_url)
    
    if response.status_code != 200:
        return []

    # Find the EC number linked to the gene
    ec_number = None
    for line in response.text.split('\n'):
        if "ORTHOLOGY" in line and "EC:" in line:
            # Extract EC number, e.g., from "EC:3.2.1.21"
            ec_number = line.split("EC:")[1].split("]")[0].strip()
            print(f"Found EC Number: {ec_number}")
            break
            
    if not ec_number:
        print("Could not determine EC number for the gene.")
        return []
        
    # --- Find reactions associated with the EC number ---
    print(f"Finding reactions for EC:{ec_number}...")
    link_url = f"http://rest.kegg.jp/link/reaction/ec:{ec_number}"
    response = requests.get(link_url)
    
    if response.status_code != 200 or not response.text:
        print(f"No reactions found for EC:{ec_number}")
        return []

    reaction_ids = [line.split('\t')[1] for line in response.text.strip().split('\n')]
    
    # --- Get details for each reaction ---
    for reaction_id in reaction_ids:
        print(f"Fetching details for reaction: {reaction_id}...")
        get_reaction_url = f"http://rest.kegg.jp/get/{reaction_id}"
        reaction_response = requests.get(get_reaction_url)
        
        if reaction_response.status_code == 200:
            name = ""
            equation = ""
            subsystem = ""
            for line in reaction_response.text.split('\n'):
                if line.startswith("NAME"):
                    name = line.replace("NAME", "").strip()
                elif line.startswith("EQUATION"):
                    equation = line.replace("EQUATION", "").strip()
                elif line.startswith("PATHWAY"):
                    subsystem = line.replace("PATHWAY", "").strip().split("  ")[1]
            
            reactions_data.append({
                "Reaction ID": reaction_id,
                "Reaction Name": name,
                "Stoichiometry": equation,
                "Subsystem": subsystem
            })
        # Be nice to the API
        time.sleep(0.1)

    return reactions_data


def main():
    """
    Main function to run the bioinformatics pipeline.
    """
    print("--- Starting Bioinformatics Pipeline ---")
    
    # 1. Recognize nucleotide string and search database
    blast_result = run_blast_search(NUCLEOTIDE_SEQUENCE)
    
    if not blast_result:
        print("Pipeline stopped due to BLAST search failure.")
        return

    # 2 & 3. Find gene and associated protein
    gene_name, protein_name = get_gene_and_protein_info(blast_result)

    if not protein_name:
        print("Could not identify protein. Cannot proceed.")
        return

    # 4 & 5. Find reactions, stoichiometry, and subsystems
    reactions = get_kegg_reactions(protein_name)
    
    if not reactions:
        print("No associated reactions found in KEGG.")
    
    # 6. Write information to an Excel-like format
    final_data = []
    if reactions:
        for reaction in reactions:
            final_data.append({
                "Input Sequence": NUCLEOTIDE_SEQUENCE[:30] + "...", # Truncate for display
                "Identified Gene": gene_name,
                "Associated Protein": protein_name,
                "Reaction ID": reaction["Reaction ID"],
                "Reaction Name": reaction["Reaction Name"],
                "Reaction Stoichiometry": reaction["Stoichiometry"],
                "Reaction Subsystem": reaction["Subsystem"],
            })
    else:
        # Add basic info even if no reactions are found
        final_data.append({
            "Input Sequence": NUCLEOTIDE_SEQUENCE[:30] + "...",
            "Identified Gene": gene_name,
            "Associated Protein": protein_name,
            "Reaction ID": "N/A",
            "Reaction Name": "N/A",
            "Reaction Stoichiometry": "N/A",
            "Reaction Subsystem": "N/A",
        })

    df = pd.DataFrame(final_data)

    try:
        df.to_excel(OUTPUT_FILENAME, index=False, engine='openpyxl')
        print(f"\n✅ Success! Data has been written to '{OUTPUT_FILENAME}'")
    except Exception as e:
        print(f"\n❌ Error writing to Excel file: {e}")

    print("--- Pipeline Finished ---")


if __name__ == "__main__":
    main()
