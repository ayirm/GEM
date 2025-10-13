import os
import re
import pandas as pd

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

# TTTTTTTTATACCTTCCAGAGCAATCTCACGTCTTGCAAAAACAGCCTGCGTTTTCATCAGTAATAGTTGGAATTTTGTAAATCTCCCGTTACCCTGATAGCGGACTTCCCTTCTGTAACCATAATGGAACCTCGTCATGTTTGAGAACATTACCGCCGCTCCTGCCGACCCGATTCTGGGCCTGGCCGATCTGTTTCGTGCCGATGAACGTCCCGGCAAAATTAACCTCGGGATTGGTGTCTATAAAGATGAGACGGGCAAAACCCCGGTACTGACCAGCGTGAAAAAGGCTGAACAGTATCTGCTCGAAAATGAAACCACCAAAAATTACCTCGGCATTGACGGCATCCCTGAATTTGGTCGCTGCACTCAGGAACTGCTGTTTGGTAAAGGTAGCGCCCTGATCAATGACAAACGTGCTCGCACGGCACAGACTCCGGGGGGCACTGGCGCACTACGCGTGGCTGCCGATTTCCTGGCAAAAAATACCAGCGTTAAGCGTGTGTGGGTGAGCAACCCAAGCTGGCCGAACCATAAGAGCGTCTTTAACTCTGCAGGTCTGGAAGTTCGTGAATACGCTTATTATGATGCGGAAAATCACACTCTTGACTTCGATGCACTGATTAACAGCCTGAATGAAGCTCAGGCTGGCGACGTAGTGCTGTTCCATGGCTGCTGCCATAACCCAACCGGTATCGACCCTACGCTGGAACAATGGCAAACACTGGCACAACTCTCCGTTGAGAAAGGCTGGTTACCGCTGTTTGACTTCGCTTACCAGGGTTTTGCCCGTGGTCTGGAAGAAGATGCTGAAGGACTGCGCGCTTTCGCGGCTATGCATAAAGAGCTGATTGTTGCCAGTTCCTACTCTAAAAACTTTGGCCTGTACAACGAGCGTGTTGGCGCTTGTACTCTGGTTGCTGCCGACAGTGAAACCGTTGATCGCGCATTCAGCCAAATGAAAGCGGCGATTCGCGCTAACTACTCTAACCCACCAGCACACGGCGCTTCTGTTGTTGCCACCATCCTGAGCAACGATGCGTTACGTGCGATTTGGGAACAAGAGCTGACTGATATGCGCCAGCGTATTCAGCGTATGCGTCAGTTGTTCGTCAATACGCTGCAGGAAAAAGGCGCAAACCGCGACTTCAGCTTTATCATCAAACAGAACGGCATGTTCTCCTTCAGTGGCCTGACAAAAGAACAAGTGCTGCGTCTGCGCGAAGAGTTTGGCGTATATGCGGTTGCTTCTGGTCGCGTAAATGTGGCCGGGATGACACCAGATAACATGGCTCCGCTGTGCGAAGCGATTGTGGCAGTGCTGTAAGCATTAAAAACAATGAAGCCCGCTGAAAAGCGGGCTGAGACTGATGACAAACGCAACATTGCCTGATGCGCTACGCTTATCAGGCCT
# Real sequence in aspC gene for aspartate aminotransferase

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
        # print(f"Top hit found: {top_hit_title}")

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

        # print("\n--- Protein Details ---")
        # print(f"ID: {protein_id}")
        # print(f"Protein Name: {protein_name}")
        # print(f"Gene Name: {gene_name}")
        # print("-----------------------\n")

        return protein_id


    except FileNotFoundError:
        print(f"Error: Could not find Blastresult.xml at the path '{blast_xml_path}'")
    except Exception as e:
        print(f"An error occurred while parsing or fetching details: {e}")

def run_kegg_search_for_reactions(accs_id:str):
    """
    Checks UniPortDB for the protein with the accession id. !Important: Does not check any other database like UniRef or UniParc
    Afterwards it finds the protein pathways on kegg
    """
    u = UniProt(verbose=False)
    k = KEGG(verbose=False)

    # result = u.search(accs_id) # !TODO: Check why columns attribute doesn't work with more than one option
    # print("Full results: \n", result)

    # entry = result.split("\n")[1].split("\t")
    
    entry = u.search(accs_id, columns="id")
    just_entry = entry.split("\n")[1].split("\t")
    print(f"Entry in full: \n {entry} \n Just entry name: \n {just_entry}")

    kegg_equ = u.mapping("UniProtKB_AC-ID","KEGG",just_entry)
    print("KEGG conversion Results:\n", kegg_equ)

    to_list = [r.get('to') for r in kegg_equ.get('results', []) if r.get('to')]
    eco_hits = [t for t in to_list if t.startswith('eco:')]
    # print(f"Kegg Conversion: {to_list} \n ECO value hopefully: {eco_hits}")

    eco_val = eco_hits[0]
    organism_code, gene_id = eco_val.split(':', 1)
    # print(f"Organism code: {organism_code} \n Gene Value: {gene_id}")

    kegg_result = k.get_pathway_by_gene(gene=gene_id,organism=organism_code)
    print("Kegg result:",kegg_result)

    return kegg_result, entry, eco_val

def parse_uniprot_table(tbl_text):
    """
    Expect header line + one data line.
    Uses a conservative split to keep 'Protein names' intact.
    Returns a dict of columns.
    """
    lines = [l for l in tbl_text.splitlines() if l.strip()]
    if not lines:
        return {}
    headers = lines[0].split()
    # maxsplit = number of columns - 1 to preserve protein names which can contain spaces
    maxsplit = max(0, len(headers) - 1)
    row_parts = lines[1].split(None, maxsplit)
    # if columns mismatch, pad
    while len(row_parts) < len(headers):
        row_parts.append('')
    return dict(zip(headers, row_parts))

def strip_ec_from_protein_name(protein_name):
    """
    Remove the '(EC x.x.x.x)' substring if present and trim whitespace.
    Leaves other parentheses (like synonyms) but removes the EC bracketed term.
    """
    # remove the substring "(EC ...)" including surrounding whitespace
    cleaned = re.sub(r'\s*\(EC\s*\d+\.\d+\.\d+\.\d+\)\s*', ' ', protein_name, flags=re.IGNORECASE)
    # collapse multiple spaces and strip
    cleaned = re.sub(r'\s+', ' ', cleaned).strip()
    return cleaned




# --- Main script execution ---

blast_xml_path = os.path.join(loc_path, 'Blastresult.xml')

# Check if BLAST result already exists
if os.path.exists(blast_xml_path):
    resp = input(f"Found existing BLAST XML at {blast_xml_path}. Skip BLAST and use existing file? [Y/n]: ")
    if resp.strip().lower() in ['', 'y', 'yes']:
        print("Using existing BLAST XML file.")
        blast_result_file = blast_xml_path
    else:
        blast_result_file = run_blast_search(nuc_sequence, loc_path)
else:
    # If it doesn't exist, run the search
    print("BLAST result not found. Running search...")
    run_blast_search(nuc_sequence, loc_path)

# Proceed with parsing if the file exists (either pre-existing or newly created)
if os.path.exists(blast_xml_path):
    protein_id = parse_blast_and_fetch_details(blast_xml_path)
    if protein_id:
        kegg_result, entry, eco_val = run_kegg_search_for_reactions(protein_id)
        if kegg_result and entry and eco_val:
            uniprot_info = parse_uniprot_table(entry)
            protein_name_raw = uniprot_info.get('Protein', '') or uniprot_info.get('Protein', '') or uniprot_info.get('Protein names', '')
            if not protein_name_raw:
                # try alternative header names if split produced different header tokens
                # join everything between the known fields: fallback to the second element if present
                protein_name_raw = list(uniprot_info.values())[3] if len(uniprot_info) >= 4 else ''
            protein_name_clean = strip_ec_from_protein_name(protein_name_raw)
            ec_matches = re.findall(r'\d+\.\d+\.\d+\.\d+', protein_name_raw)
            protein_name_cell = protein_name_clean + ('; EC=' + ';'.join(ec_matches) if ec_matches else '')

            genes_cell = ';'.join(eco_val) if isinstance(eco_val, (list, tuple)) else str(eco_val)
            pathways_cell = '; '.join([f"{k}: {v}" for k, v in kegg_result.items()])

            row = {
                'Nucleotide sequence': nuc_sequence,
                'Protein Name': protein_name_cell,
                'Genes': genes_cell,
                'KEGG Pathways': pathways_cell
            }

            df = pd.DataFrame([row])

            # Write to Excel (single sheet)
            out_file = "combined_results_single_sheet.xlsx"
            df.to_excel(out_file, index=False, sheet_name='Results')
            print("Wrote", out_file)

else:
    print("Could not find or create BLAST result file. Exiting.")
