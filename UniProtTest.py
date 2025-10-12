import os
import time
import pandas as pd
from bioservices import *
from xml.etree import ElementTree as ET
from Bio.Blast import NCBIXML
import io
from Bio import Entrez, SeqIO
import requests

Entrez.email = "dogayasemen.testere@ogr.iu.edu.tr"

# GGTGGTCTGCCTCGCATAAAGCGGTATGAAAATGGATTGAAGCCCGGGCCGTGGATTCTACTCAACTTTCGTCTTTCGAGAAAGACTCCGGGATCCTGAGTATTAAAAAGAAGATCTTTATATAGAGATCTGTTCTATTGTGATCTCTTATTAGGATCGACTCTCTTTGTGGATAAGTCGGATCCGCGAAGTAAGATCAAAAGCTTAAGAAGGATCACTATCTGTGAATGATCGGTGATC
# Test String

def run_ncbi_blast_uniprotkb(sequence: str, loc_path: str):
    """
    Submits a BLASTx search to NCBI, targeting the UniProtKB database,
    and saves the result to a file.
    """
    print("Submitting BLASTx search via NCBI to UniProtKB database...")
    s = NCBIblast(verbose=False)
    
    # Use keyword arguments for clarity and correctness
    job_id = s.run(
        program='blastx',
        sequence=sequence,
        stype='dna',
        database='uniprotkb',
        email='dogayasemen.testere@ogr.iu.edu.tr'  # NCBI BLAST requires an email address
    )
    print(f"Job submitted with ID: {job_id}")

    # Poll for job completion instead of using a blocking wait
    while True:
        status = s.get_status(job_id)
        print(f"Current job status: {status}")
        if status == 'FINISHED':
            break
        elif status in ['FAILED', 'UNKNOWN']:
            print(f"Job failed or has an unknown status: {status}")
            return None
        # Wait for 30 seconds before checking again
        print("Waiting for 30 seconds before the next check...")
        time.sleep(30)

    print("BLAST search finished. Retrieving results...")
    # Retrieve results in XML format
    xml_results = s.get_result(job_id, 'xml')

    # Save the XML results to a file before parsing
    save_path = os.path.join(loc_path, 'ncbi_blast_result.xml')
    try:
        with open(save_path, "w") as out_handle:
            out_handle.write(xml_results)
        print(f"Successfully saved BLAST XML results to {save_path}")
        return save_path
    except Exception as e:
        print(f"Error saving XML file: {e}")
        return None

def parse_ncbi_blast(xml_file_path: str, skip_accessions=None):
    """
    Parses the NCBI BLAST XML result from a file to get the top hit's UniProt accession.
    This function handles two cases:
      - Direct NCBI BLAST XML (parsed by Biopython's NCBIXML)
      - EBIApplicationResult wrapper produced by the EBI tools (extract embedded XML from <output> elements)
    """
    try:
        # Read raw file first to inspect header
        with open(xml_file_path, 'r', encoding='utf-8', errors='replace') as result_handle:
            raw = result_handle.read()

        # Quick check for EBI wrapper
        if raw.lstrip().startswith("<?xml") and "EBIApplicationResult" in raw.split('\n', 2)[1]:
            # Parse wrapper
            try:
                root = ET.fromstring(raw)
            except Exception as e:
                print(f"EBI wrapper XML parsed failed: {e}")
                return None

            # Search for any <output> elements and look for embedded BLAST XML
            found_inner = None
            for elem in root.iter():
                # use local-name match to be robust to namespaces
                tag = elem.tag
                if isinstance(tag, str) and tag.lower().endswith('output'):
                    # text may contain CDATA with BLAST XML
                    text = (elem.text or '').strip()
                    if '<BlastOutput' in text or '<Blast' in text or '<Iteration>' in text:
                        found_inner = text
                        break
                    # sometimes child contains the xml
                    for child in elem:
                        inner_text = (child.text or '').strip()
                        if '<BlastOutput' in inner_text or '<Blast' in inner_text or '<Iteration>' in inner_text:
                            found_inner = inner_text
                            break
                    if found_inner:
                        break

            if not found_inner:
                print("EBI wrapper found but no embedded BLAST XML was detected in <output> elements.")
                # Try to extract top hit accession directly from the EBI wrapper's SequenceSimilaritySearchResult
                ns = {'ebi': 'http://www.ebi.ac.uk/schema'}
                hits = root.findall('.//ebi:hit', ns)
                if hits:
                    # iterate and respect skip_accessions
                    for top in hits:
                        ac = top.get('ac') or top.get('id')
                        if ac:
                            if skip_accessions and ac in skip_accessions:
                                continue
                            print(f"Extracted accession from EBI wrapper top hit: {ac}")
                            return ac
                # Fallback: try regex extraction of accession attributes from the raw XML text
                import re as _re
                ac_matches = _re.findall(r'ac="([A-Za-z0-9_\-\.]+)"', raw)
                if ac_matches:
                    for ac in ac_matches:
                        if skip_accessions and ac in skip_accessions:
                            continue
                        print(f"Extracted accession from EBI wrapper (regex): {ac}")
                        return ac
                print("Available <output> element tags and their sizes:")
                for elem in root.iter():
                    tag = elem.tag
                    if isinstance(tag, str) and tag.lower().endswith('output'):
                        text = (elem.text or '')
                        print(f" - {tag} (length {len(text)})")
                return None

            # We have inner BLAST XML as string in found_inner
            xml_str = found_inner
            # Save inner xml for debugging
            inner_save = xml_file_path.replace('.xml', '.inner.xml')
            with open(inner_save, 'w', encoding='utf-8') as fh:
                fh.write(xml_str)
            print(f"Extracted embedded BLAST XML saved to {inner_save}")

            # Parse using Biopython
            try:
                result_handle = io.StringIO(xml_str)
                blast_record = NCBIXML.read(result_handle)
            except Exception as e:
                print(f"Failed to parse extracted BLAST XML: {e}")
                return None

        else:
            # Not an EBI wrapper: parse directly
            with open(xml_file_path, 'r', encoding='utf-8', errors='replace') as result_handle:
                blast_record = NCBIXML.read(result_handle)

                if not blast_record.alignments:
                    print("No hits found in BLAST result.")
                    return None

                # For uniprotkb searches, the accession is the UniProt ID
                for aln in blast_record.alignments:
                    acc = aln.accession
                    if skip_accessions and acc in skip_accessions:
                        continue
                    print(f"Found top hit UniProt ID: {acc}")
                    return acc
                print("No non-skipped hits found in BLAST record.")
                return None

    except Exception as e:
        print(f"An error occurred while parsing BLAST XML: {e}")
        print("--- Displaying the first 10 lines of the problematic XML file ---")
        try:
            with open(xml_file_path, 'r', encoding='utf-8', errors='replace') as f:
                for i, line in enumerate(f):
                    if i >= 10:
                        break
                    print(line.strip())
        except Exception as read_e:
            print(f"Could not read the file to display content: {read_e}")
        print("--- End of file content ---")
        return None


def fetch_uniprot_json(uniprot_id: str):
    """Return UniProt JSON dict for accession or None on failure."""
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            return r.json()
    except Exception:
        return None
    return None


def extract_ecs_from_uniprot_json(data: dict):
    """Return list of EC numbers found in UniProt JSON structure."""
    import re

    def collect_ec(obj):
        ecs = []
        if isinstance(obj, dict):
            for k, v in obj.items():
                if k.lower().startswith('ec') and isinstance(v, (str, dict)):
                    if isinstance(v, str):
                        ecs += re.findall(r'\d+\.\d+\.\d+\.\d+', v)
                    elif isinstance(v, dict):
                        val = v.get('value')
                        if val:
                            ecs += re.findall(r'\d+\.\d+\.\d+\.\d+', val)
                else:
                    ecs += collect_ec(v)
        elif isinstance(obj, list):
            for item in obj:
                ecs += collect_ec(item)
        return ecs

    ecs = collect_ec(data.get('proteinDescription', {}))
    ecs += collect_ec(data.get('comments', []))
    # dedupe
    ecs = [e for i, e in enumerate(ecs) if e and e not in ecs[:i]]
    return ecs


def extract_accessions_from_blast_xml(xml_file_path: str):
    """Return a list of accession-like strings found in a BLAST XML/wrapper file.
    Uses XML parsing + regex fallback to collect ac/id attributes.
    """
    import re
    accs = []
    try:
        with open(xml_file_path, 'r', encoding='utf-8', errors='replace') as fh:
            raw = fh.read()

        # regex for ac="..." and id="..."
        acs = re.findall(r'ac\s*=\s*"([A-Za-z0-9_\-\.]+)"', raw)
        ids = re.findall(r'id\s*=\s*"([A-Za-z0-9_\-\.]+)"', raw)
        accs = acs + ids

        # also try to parse any 'hit' elements under SequenceSimilaritySearchResult if present
        try:
            root = ET.fromstring(raw)
            for hit in root.findall('.//{http://www.ebi.ac.uk/schema}hit'):
                a = hit.get('ac') or hit.get('id')
                if a:
                    accs.append(a)
        except Exception:
            pass

        # dedupe while preserving order
        seen = set()
        out = []
        for a in accs:
            if a not in seen:
                seen.add(a)
                out.append(a)
        return out
    except Exception as e:
        print(f"Could not extract accessions from BLAST XML: {e}")
        return []

def get_details_and_reactions(uniprot_id: str, loc_path: str):
    """
    Gets protein details from UniProt, finds reactions from KEGG, and exports to Excel.
    """
    try:
        # 1. Get protein details from UniProt
        u = UniProt(verbose=False)
        print(f"Fetching details for UniProt ID: {uniprot_id}")

        # First try direct UniProt REST accession JSON (more reliable than search)
        result_str = None
        try:
            url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
            r = requests.get(url, timeout=10)
            if r.status_code == 200:
                data = r.json()
                # extract protein name
                seq_from_uniprot = None
                protein_name = None
                pd = data.get('proteinDescription', {})
                if 'recommendedName' in pd and pd['recommendedName']:
                    rn = pd['recommendedName']
                    fn = rn.get('fullName')
                    if isinstance(fn, dict):
                        protein_name = fn.get('value')
                    elif isinstance(fn, str):
                        protein_name = fn
                if not protein_name and 'submissionNames' in pd and pd['submissionNames']:
                    sub = pd['submissionNames'][0]
                    fn = sub.get('fullName')
                    if isinstance(fn, dict):
                        protein_name = fn.get('value')
                    elif isinstance(fn, str):
                        protein_name = fn

                # extract gene name (prefer geneName then orfNames)
                gene_name = 'N/A'
                genes = data.get('genes', [])
                if genes:
                    g0 = genes[0]
                    if 'geneName' in g0 and g0['geneName']:
                        gene_name = g0['geneName'].get('value', 'N/A')
                    elif 'orfNames' in g0 and g0['orfNames']:
                        gene_name = g0['orfNames'][0].get('value', 'N/A')

                # extract organism scientific name for KEGG mapping
                organism_name = None
                try:
                    organism_name = data.get('organism', {}).get('scientificName')
                except Exception:
                    organism_name = None

                # extract EC numbers by searching common locations
                ec_numbers = []
                # check recommendedName.ecNumber or submissionNames.ecNumber
                def collect_ec(obj):
                    ecs = []
                    if isinstance(obj, dict):
                        for k, v in obj.items():
                            if k.lower().startswith('ec') and isinstance(v, (str, dict)):
                                if isinstance(v, str):
                                    ecs += re.findall(r'\d+\.\d+\.\d+\.\d+', v)
                                elif isinstance(v, dict):
                                    val = v.get('value')
                                    if val:
                                        ecs += re.findall(r'\d+\.\d+\.\d+\.\d+', val)
                            else:
                                ecs += collect_ec(v)
                    elif isinstance(obj, list):
                        for item in obj:
                            ecs += collect_ec(item)
                    return ecs

                import re
                ec_numbers = collect_ec(data.get('proteinDescription', {}))
                # capture sequence from UniProt JSON for use in BLASTp fallback
                try:
                    seq_from_uniprot = data.get('sequence', {}).get('value')
                    if seq_from_uniprot:
                        seq_from_uniprot = seq_from_uniprot.replace('\n', '').strip()
                except Exception:
                    seq_from_uniprot = None
                # also search comments and crossrefs
                ec_numbers += collect_ec(data.get('comments', []))
                ec_numbers = [e for i, e in enumerate(ec_numbers) if e and e not in ec_numbers[:i]]
                ec_numbers_str = ';'.join(ec_numbers)

                header = ['Entry', 'Protein names', 'Gene Names', 'EC number']
                data_line = [uniprot_id, protein_name or 'N/A', gene_name or 'N/A', ec_numbers_str]
                result_str = '\t'.join(header) + '\n' + '\t'.join(data_line)
                print('Fetched UniProt JSON successfully')
        except Exception as e:
            result_str = None

        # If direct JSON fetch failed, try several query styles with bioservices
        if not result_str:
            # Try several query styles to be robust against accession vs entry name issues
            query_variants = [f"accession:{uniprot_id}", f"id:{uniprot_id}", uniprot_id]
            for q in query_variants:
                try:
                    result_str = u.search(q, columns="id,protein_name,genes,ec")
                except Exception:
                    result_str = None
                if result_str:
                    print(f"UniProt search successful with query: {q}")
                    break

        # Fallback: try to retrieve the entry directly
        if not result_str:
            try:
                print("Search returned no results — attempting to retrieve the UniProt entry directly with u.get()...")
                entry_text = u.get(uniprot_id)
                if entry_text:
                    # crude parsing: look for lines starting with 'ID', 'DE', 'GN', 'DR   EC;'
                    protein_name = 'N/A'
                    gene_name = 'N/A'
                    ec_numbers = []
                    for line in entry_text.splitlines():
                        if line.startswith('DE   RecName: Full=') and protein_name == 'N/A':
                            protein_name = line.replace('DE   RecName: Full=', '').strip().rstrip(';')
                        if line.startswith('GN   Name=') and gene_name == 'N/A':
                            gene_name = line.replace('GN   Name=', '').split(';')[0].strip()
                        if 'EC=' in line:
                            parts = line.split('EC=')
                            for p in parts[1:]:
                                ec = p.split(';')[0].strip()
                                if ec:
                                    ec_numbers.append(ec)
                    if not ec_numbers:
                        ec_numbers_str = ''
                    else:
                        ec_numbers_str = ';'.join(ec_numbers)
                    # construct a pseudo result_str to reuse later logic
                    header = ['Entry', 'Protein names', 'Gene Names', 'EC number']
                    data = [uniprot_id, protein_name, gene_name, ec_numbers_str]
                    result_str = '\t'.join(header) + '\n' + '\t'.join(data)
            except Exception:
                result_str = None

        if not result_str:
            print(f"Could not find details for UniProt ID {uniprot_id}")
            return

        # Parse the tab-separated result string
        lines = result_str.strip().split('\n')
        header = lines[0].split('\t')
        data = lines[1].split('\t')
        protein_info = dict(zip(header, data))

        protein_name = protein_info.get('Protein names', 'N/A')
        gene_name = protein_info.get('Gene Names', 'N/A')
        ec_numbers_str = protein_info.get('EC number', '')

        if not ec_numbers_str:
            print(f"No EC number found for UniProt ID {uniprot_id}.")

            # --- Fallback 0: Try KEGG gene lookup directly using UniProt gene/orf name and organism ---
            try:
                if gene_name and gene_name != 'N/A' and organism_name:
                    print(f"Attempting direct KEGG gene lookup using organism '{organism_name}' and gene '{gene_name}'...")
                    ktmp = KEGG(verbose=False)
                    # try to find the KEGG organism code by scanning the organism list for a substring match
                    org_list = ktmp.list('organism')
                    org_code = None
                    for line in org_list.splitlines():
                        parts = line.split('\t')
                        if len(parts) >= 3:
                            code = parts[1]
                            name = parts[2]
                            if organism_name.lower() in name.lower():
                                org_code = code
                                break
                    if not org_code:
                        # fallback: try first two words of organism name
                        short = ' '.join(organism_name.split()[:2])
                        res = ktmp.find('organism', short)
                        if res:
                            org_code = res.strip().split('\n')[0].split('\t')[1]

                    if org_code:
                        # try both alias code (parts[1]) and the T... code (parts[0]) if available
                        possible_codes = [org_code]
                        # attempt to detect a T code from the same line by searching the list again
                        # some KEGG list lines look like: 'T03934\tcama\tCitrobacter amalonaticus ...'
                        for line in org_list.splitlines():
                            if organism_name.lower() in line.lower():
                                parts = line.split('\t')
                                if parts and parts[0] not in possible_codes:
                                    possible_codes.insert(0, parts[0])
                                break

                        entry = None
                        for code in possible_codes:
                            kegg_id = f"{code}:{gene_name}"
                            print(f"Trying KEGG ID: {kegg_id}")
                            entry = ktmp.get(kegg_id)
                            if entry:
                                break
                        if entry:
                            parsed = ktmp.parse(entry)
                            ec_from_gene = None
                            if 'DBLINKS' in parsed and 'EC' in parsed['DBLINKS']:
                                ec_from_gene = parsed['DBLINKS']['EC'].split()[0]
                            if not ec_from_gene and 'ENZYME' in parsed:
                                ec_from_gene = parsed['ENZYME'][0]
                            if ec_from_gene:
                                ec_numbers = [ec_from_gene]
                                ec_numbers_str = ';'.join(ec_numbers)
                                print(f"Found EC from KEGG gene entry: {ec_from_gene}")
                            else:
                                print(f"No EC found in KEGG gene entry for {kegg_id}")
                        else:
                            print(f"No KEGG entry for {kegg_id}")
                    else:
                        print(f"Could not determine KEGG organism code for '{organism_name}'")
            except Exception as e:
                print(f"Direct KEGG gene lookup failed: {e}")

            # --- Fallback 1: Use NCBI protein record to get locus_tag + organism, then query KEGG ---
            try:
                print("Attempting NCBI fallback: fetch protein record to extract locus_tag and organism...")
                # Try to find a matching NCBI protein by accession
                handle = Entrez.esearch(db='protein', term=f'{uniprot_id}[Accession]')
                rec = Entrez.read(handle)
                handle.close()
                idlist = rec.get('IdList', [])
                protein_rec = None
                if idlist:
                    prot_id = idlist[0]
                    print(f"Found NCBI protein ID: {prot_id}, fetching GenBank record...")
                    h = Entrez.efetch(db='protein', id=prot_id, rettype='gb', retmode='text')
                    protein_rec = SeqIO.read(h, 'genbank')
                    h.close()
                else:
                    print("No NCBI protein found by accession. Trying search by accession as text...")
                    h2 = Entrez.esearch(db='protein', term=uniprot_id)
                    r2 = Entrez.read(h2)
                    h2.close()
                    if r2.get('IdList'):
                        prot_id = r2['IdList'][0]
                        h3 = Entrez.efetch(db='protein', id=prot_id, rettype='gb', retmode='text')
                        protein_rec = SeqIO.read(h3, 'genbank')
                        h3.close()

                if protein_rec:
                    locus_tag = None
                    organism = protein_rec.annotations.get('organism')
                    for f in protein_rec.features:
                        if f.type == 'CDS' and 'locus_tag' in f.qualifiers:
                            locus_tag = f.qualifiers['locus_tag'][0]
                            break

                    print(f"NCBI record organism: {organism}, locus_tag: {locus_tag}")

                    if locus_tag and organism:
                        k = KEGG(verbose=False)
                        search_term = ' '.join(organism.split()[:2])
                        org_search_result = k.find('organism', search_term)
                        if org_search_result:
                            org_code = org_search_result.strip().split('\n')[0].split('\t')[1]
                            kegg_id = f"{org_code}:{locus_tag}"
                            print(f"Constructed KEGG ID: {kegg_id}. Fetching entry...")
                            gene_entry_str = k.get(kegg_id)
                            if gene_entry_str:
                                gene_entry = k.parse(gene_entry_str)
                                ec = None
                                if 'DBLINKS' in gene_entry and 'EC' in gene_entry['DBLINKS']:
                                    ec = gene_entry['DBLINKS']['EC'].split()[0]
                                # some entries list EC under 'ENZYME' or 'ORTHOLOGY' or 'COMMENT'
                                if not ec and 'ENZYME' in gene_entry:
                                    ec = gene_entry['ENZYME'][0]
                                if ec:
                                    print(f"Found EC from KEGG gene entry: {ec}")
                                    ec_numbers = [ec]
                                    # proceed to collect reactions below
                                else:
                                    print("No EC found in KEGG gene entry for locus tag.")
                            else:
                                print(f"Could not fetch KEGG entry for {kegg_id}")
                        else:
                            print(f"Could not map organism '{organism}' to a KEGG code.")
                else:
                    print("NCBI protein record not found during fallback.")
            except Exception as e:
                print(f"NCBI fallback failed: {e}")

            # If after NCBI fallback we still don't have EC numbers, try BLASTp fallback (homology)
            if not ec_numbers_str and ('ec_numbers' not in locals() or not ec_numbers):
                try:
                    print("Attempting BLASTp fallback: fetch protein sequence and BLASTp against UniProtKB to find homologs...")
                    # obtain protein sequence from NCBI if we have protein_rec
                    seq_to_blast = None
                    if 'protein_rec' in locals() and protein_rec:
                        seq_to_blast = str(protein_rec.seq)
                    else:
                        # try to fetch sequence from UniProt entry text if possible
                        try:
                            u = UniProt(verbose=False)
                            fasta = u.get(uniprot_id, frmt='fasta') if hasattr(u, 'get') else None
                        except Exception:
                            fasta = None
                        if fasta and fasta.startswith('>'):
                            # extract sequence
                            seq_to_blast = ''.join(fasta.split('\n')[1:])
                    # fallback: if we fetched UniProt JSON earlier, use its sequence
                        if not seq_to_blast and 'seq_from_uniprot' in locals() and seq_from_uniprot:
                            seq_to_blast = seq_from_uniprot

                        # debug: report which sequence source will be used (if any)
                        if seq_to_blast:
                            print(f"BLASTp will use a sequence of length {len(seq_to_blast)}")
                        else:
                            print("No protein sequence available from NCBI, UniProt FASTA, or UniProt JSON to run BLASTp fallback.")

                    if seq_to_blast:
                        # before submitting a new BLASTp job, check if a previous blastp_result.xml exists and reuse it
                        tmp = os.path.join(loc_path, 'blastp_result.xml')
                        if os.path.exists(tmp):
                            print(f"Existing BLASTp result found at {tmp}, parsing it before submitting a new job...")
                            top_acc = parse_ncbi_blast(tmp, skip_accessions=[uniprot_id])
                            if top_acc and top_acc != uniprot_id:
                                print(f"Found homolog in existing BLASTp: {top_acc}")
                                return get_details_and_reactions(top_acc, loc_path)
                            else:
                                print("No useful homolog found in existing BLASTp result, will submit a new BLASTp job against Swiss-Prot (reviewed).")

                        print("Protein sequence obtained, submitting BLASTp to Swiss-Prot (curated UniProt) to prioritize annotated homologs...")
                        s = NCBIblast(verbose=False)
                        # use 'swissprot' database to favor reviewed entries with EC annotations
                        # attempt to submit BLASTp, preferring swissprot but falling back to uniprotkb
                        job = None
                        for db_choice in ('swissprot', 'uniprotkb'):
                            try:
                                job = s.run(program='blastp', sequence=seq_to_blast, stype='protein', database=db_choice, email=Entrez.email)
                                print(f"BLASTp job submitted to {db_choice}: {job}")
                                # if the returned job looks like an HTTP error code or empty, treat as failure
                                if not job or str(job).isdigit() or str(job).upper().startswith('400') or str(job).upper().startswith('ERROR'):
                                    print(f"Submission to {db_choice} returned an unexpected response: {job}")
                                    job = None
                                    continue
                                break
                            except Exception as _e:
                                print(f"Submission to {db_choice} failed: {_e}")
                                job = None
                        if not job:
                            print("Could not submit BLASTp job to any supported database.")
                        else:
                            # poll status with timeout/backoff
                            max_wait = 600  # seconds
                            waited = 0
                            backoff = 5
                            while waited < max_wait:
                                st = s.get_status(job)
                                print(f"BLASTp job status: {st}")
                                if st == 'FINISHED':
                                    break
                                if st in ['FAILED', 'UNKNOWN', 'NOT_FOUND']:
                                    print(f"BLASTp job failed or not available: {st}")
                                    break
                                time.sleep(backoff)
                                waited += backoff
                                # gentle exponential backoff up to 30s
                                backoff = min(backoff * 1.5, 30)

                            if st == 'FINISHED':
                                try:
                                    xml = s.get_result(job, 'xml')
                                    # save and parse similarly to previous
                                    tmp = os.path.join(loc_path, 'blastp_result.xml')
                                    with open(tmp, 'w', encoding='utf-8') as fh:
                                        fh.write(xml)
                                    print(f"Saved BLASTp result to {tmp}")
                                    # Instead of only taking the top hit, scan the saved XML for all accession attributes
                                    import re
                                    raw = xml
                                    accs = re.findall(r'ac="([A-Za-z0-9_\-\.]+)"', raw)
                                    accs = [a for a in accs if a and a != uniprot_id]
                                    found_ec = None
                                    for acc in accs:
                                        print(f"Checking homolog accession from BLASTp: {acc}")
                                        json_h = fetch_uniprot_json(acc)
                                        if not json_h:
                                            print(f" - Could not fetch UniProt JSON for {acc}")
                                            continue
                                        ecs = extract_ecs_from_uniprot_json(json_h)
                                        if ecs:
                                            print(f" - Found EC(s) {ecs} on homolog {acc}")
                                            found_ec = (acc, ecs)
                                            break
                                    if found_ec:
                                        homolog_acc, homolog_ecs = found_ec
                                        # attach ec numbers and proceed
                                        ec_numbers = homolog_ecs
                                        ec_numbers_str = ';'.join(ec_numbers)
                                        print(f"Proceeding with EC(s) from homolog {homolog_acc}: {ec_numbers_str}")
                                        return get_details_and_reactions(homolog_acc, loc_path)
                                    else:
                                        print("No homologs with EC annotations found in BLASTp results.")
                                except Exception as _e:
                                    print(f"Failed to retrieve/parse BLASTp result: {_e}")
                            else:
                                print("BLASTp did not finish successfully within timeout or failed; no result to parse.")
                    else:
                        print("Could not obtain a protein sequence to run BLASTp fallback.")
                except Exception as e:
                    print(f"BLASTp fallback failed: {e}")

            # After all fallbacks, if ec_numbers still not found, give up
            if ('ec_numbers' not in locals() or not ec_numbers):
                print("All fallbacks failed — no EC numbers found. Cannot continue to KEGG reactions.")
                return

            # If we obtained ec_numbers via fallback, continue to reaction fetching below
            if 'ec_numbers' in locals() and ec_numbers:
                print(f"Proceeding with EC numbers from fallback: {ec_numbers}")
                # reuse existing KEGG reaction collection logic by assigning ec_numbers_str
                ec_numbers_str = ';'.join(ec_numbers)

        # An entry can have multiple EC numbers, separated by ';'
        ec_numbers = [ec.strip() for ec in ec_numbers_str.split(';') if ec.strip()]
        print(f"Found EC Number(s): {', '.join(ec_numbers)}")

        # 2. Get reactions from KEGG for each EC number
        k = KEGG(verbose=False)
        all_reaction_data = []

        for ec_number in ec_numbers:
            print(f"Fetching reactions for EC: {ec_number}")
            enzyme_entry_str = k.get(f"ec:{ec_number}")
            if not enzyme_entry_str:
                print(f"Could not get entry for EC {ec_number}")
                continue

            enzyme_entry = k.parse(enzyme_entry_str)

            if 'REACTION' not in enzyme_entry:
                print(f"No reactions found for EC number {ec_number}.")
                continue

            reaction_ids = enzyme_entry['REACTION']
            print(f"Found {len(reaction_ids)} reactions for {ec_number}: {', '.join(reaction_ids)}")

            for rxn_id in reaction_ids:
                reaction_entry_str = k.get(f"rn:{rxn_id}")
                reaction_entry = k.parse(reaction_entry_str)

                reaction_info = {
                    "Gene Name": gene_name,
                    "Protein Name": protein_name,
                    "EC Number": ec_number,
                    "Reaction ID": rxn_id,
                    "Reaction Name": reaction_entry.get('NAME', ['N/A'])[0],
                    "Stoichiometry": reaction_entry.get('EQUATION', 'N/A'),
                    "Reversibility": 'Reversible' if '<=>' in reaction_entry.get('EQUATION', '') else 'Irreversible'
                }
                all_reaction_data.append(reaction_info)

        # 3. Export to Excel
        if not all_reaction_data:
            print("No reaction data to export.")
            return

        df = pd.DataFrame(all_reaction_data)
        excel_path = os.path.join(loc_path, 'uniprot_kegg_reactions.xlsx')
        df.to_excel(excel_path, index=False)
        print(f"\nSuccessfully exported reaction data to {excel_path}")

    except Exception as e:
        print(f"An error occurred: {e}")


# --- Main script execution ---
if __name__ == "__main__":
    nuc_sequence = input("Sequence to be searched: ")
    loc_path = input("Directory to save results: \n")

    if loc_path and not os.path.exists(loc_path):
        os.makedirs(loc_path, exist_ok=True)

    # Path where the BLAST XML will be saved
    save_path = os.path.join(loc_path, 'ncbi_blast_result.xml')

    # If a saved XML exists, offer the option to skip running BLAST
    if os.path.exists(save_path):
        resp = input(f"Found existing BLAST XML at {save_path}. Skip BLAST and use existing file? [Y/n]: ")
        if resp.strip().lower() in ['', 'y', 'yes']:
            print("Using existing BLAST XML file.")
            blast_result_file = save_path
        else:
            blast_result_file = run_ncbi_blast_uniprotkb(nuc_sequence, loc_path)
    else:
        # No existing file, run BLAST and save
        blast_result_file = run_ncbi_blast_uniprotkb(nuc_sequence, loc_path)

    if blast_result_file:
        # First, extract all accessions from the BLAST XML and try UniProt JSON for each
        accessions = extract_accessions_from_blast_xml(blast_result_file)
        print(f"Accessions extracted from BLAST XML: {accessions}")
        found = False
        for acc in accessions:
            print(f"Probing UniProt for accession {acc}...")
            jsond = fetch_uniprot_json(acc)
            if not jsond:
                print(f" - No UniProt JSON available for {acc}")
                continue
            ecs = extract_ecs_from_uniprot_json(jsond)
            if ecs:
                print(f" - Found EC(s) {ecs} on accession {acc}. Proceeding to KEGG mapping.")
                get_details_and_reactions(acc, loc_path)
                found = True
                break
            else:
                print(f" - No ECs on UniProt entry {acc}.")

        if not found:
            # fallback: use the original parse behavior to get the top hit and run the existing pipeline
            top_uniprot_id = parse_ncbi_blast(blast_result_file)
            if top_uniprot_id:
                get_details_and_reactions(top_uniprot_id, loc_path)
