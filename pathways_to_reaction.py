import requests
from bs4 import BeautifulSoup
from bioservices import KEGG
import pandas as pd
import time
import re

# ==============================================================================
# 1. GİRDİ TÜRÜNÜ VE İLGİLİ ID'LERİ BELİRLEYEN FONKSİYONLAR
# ==============================================================================

def is_protein_sequence(text):
    """Girdinin bir protein sekansı olup olmadığını anlamak için basit bir kontrol."""
    valid_chars = "ACDEFGHIKLMNPQRSTVWY"
    if len(text) > 20 and sum(c.upper() in valid_chars for c in text) / len(text) > 0.9:
        return True
    return False

def get_gene_from_sequence(sequence):
    """Verilen protein sekansını KEGG BLAST ile arar ve en iyi E. coli gen ID'sini döndürür."""
    print("Protein sekansı KEGG BLAST ile E. coli veritabanında aranıyor...")
    print("(Bu işlem 15-30 saniye sürebilir, lütfen bekleyin...)")
    
    # --- DEĞİŞİKLİK BURADA: ESKİ URL YENİSİYLE GÜNCELLENDİ ---
    url = "https://www.genome.jp/tools/blast/blast.cgi" 
    
    payload = {
        'dbselect': 'eco', 'program': 'blastp', 'query': sequence, 'bexec': 'Bexec'
    }
    
    try:
        response = requests.post(url, data=payload, timeout=90) # Timeout süresini biraz artırdık
        response.raise_for_status()
        soup = BeautifulSoup(response.text, 'html.parser')
        pre_tag_content = soup.find('pre').text if soup.find('pre') else ''
        match = re.search(r"(eco:[b-d]\d{4})", pre_tag_content)
        
        if match:
            gene_id = match.group(1)
            print(f"✅ BLAST Eşleşmesi Bulundu: {gene_id}")
            return gene_id
        else:
            print("❌ BLAST aramasında E. coli için anlamlı bir gen eşleşmesi bulunamadı.")
            return None
    except requests.exceptions.RequestException as e:
        print(f"HATA: KEGG BLAST'a bağlanırken bir sorun oluştu: {e}")
        return None

def get_reactions_from_gene(gene_input):
    """Gen ID'si veya Gen adına göre reaksiyonları bulur."""
    k = KEGG()
    print(f"Gen Adı/ID '{gene_input}' için arama yapılıyor...")
    
    gene_id = f"eco:{gene_input}"

    print(f"'{gene_id}' ile ilişkili reaksiyonlar aranıyor...")
    try:
        reaction_links = k.link("reaction", gene_id)
        if not isinstance(reaction_links, str) or not reaction_links.strip():
             print(f"'{gene_id}' için bağlantılı reaksiyon bulunamadı.")
             return None
        
        reaction_ids = [line.split('\t')[1].split(':')[1] for line in reaction_links.strip().split('\n')]
        return reaction_ids
    except Exception:
        print(f"'{gene_id}' geçerli bir ID gibi görünmüyor veya bağlantılı reaksiyonu yok.")
        return None

def get_reactions_from_pathway(pathway_name):
    """Yolak adına göre reaksiyonları bulur."""
    k = KEGG()
    print(f"Yolak Adı '{pathway_name}' için arama yapılıyor...")
    try:
        all_paths_str = k.list("pathway", "eco")
        if not isinstance(all_paths_str, str):
             print("⚠️ KEGG'den yolak listesi alınamadı (API hatası).")
             return None
        matches = [p for p in all_paths_str.split("\n") if pathway_name.lower() in p.lower()]
        
        if not matches:
            print(f"E. coli ('eco') içinde '{pathway_name}' ile eşleşen bir yolak bulunamadı.")
            return None
            
        first_match = matches[0]
        print(f"🔎 Eşleşme bulundu: {first_match}")
        path_id = first_match.split("\t")[0].split(":")[-1]
        map_id = "map" + path_id[3:]
        
        print(f"'{map_id}' haritasıyla ilişkili reaksiyonlar aranıyor...")
        reaction_links = k.link("reaction", map_id)
        
        reaction_ids = [line.split('\t')[1].split(':')[1] for line in reaction_links.strip().split('\n')]
        return reaction_ids
    except Exception as e:
        print(f"⚠️ Yolak işlenirken hata oluştu: {e}")
        return None

# ==============================================================================
# 2. REAKSİYON DETAYLARINI ÇEKEN VE EXCEL'E YAZAN ANA FONKSİYON
# ==============================================================================

def fetch_reaction_details(reaction_ids, output_name):
    """Verilen reaksiyon ID listesi için detayları çeker ve bir DataFrame döndürür."""
    k = KEGG()
    k.TIMEOUT = 60
    compound_name_cache = {}
    
    if not reaction_ids:
        print("Detayları alınacak reaksiyon bulunamadı.")
        return pd.DataFrame()

    print(f"Toplam {len(reaction_ids)} benzersiz reaksiyon bulundu. Detaylar çekiliyor...")
    print("UYARI: Okunabilir stokiyometri oluşturmak için ek API sorguları yapılacak, bu işlem yavaş olabilir.")

    reaction_data = []
    for i, rxn_id in enumerate(reaction_ids, 1):
        print(f"İşleniyor: {i}/{len(reaction_ids)} ({rxn_id})")
        try:
            rxn_data_str = k.get(rxn_id)
            if not rxn_data_str: continue

            rxn_name, equation_ids, ec_numbers, reversibility = "", "", "", "N/A"
            readable_stoichiometry = ""

            for line in rxn_data_str.split("\n"):
                if line.startswith("NAME"):
                    rxn_name = line.replace("NAME", "").strip()
                elif line.startswith("ENZYME"):
                    found_ecs = re.findall(r"\d+\.\d+\.\d+\.(?:\d+|\-)", line)
                    ec_numbers = " ".join(found_ecs)
                elif line.startswith("EQUATION"):
                    equation_ids = line.replace("EQUATION", "").strip()
                    reversibility = "Reversible" if "<=>" in equation_ids else "Irreversible"
                    readable_stoichiometry = equation_ids
                    compound_ids = set(re.findall(r"C\d{5}", equation_ids))
                    
                    for cpd_id in compound_ids:
                        cpd_name = compound_name_cache.get(cpd_id)
                        if not cpd_name:
                            cpd_data = k.get(cpd_id)
                            cpd_name = "Unknown"
                            for cpd_line in cpd_data.split('\n'):
                                if cpd_line.startswith("NAME"):
                                    cpd_name = cpd_line.replace("NAME", "").strip().split(';')[0]
                                    break
                            compound_name_cache[cpd_id] = cpd_name
                            time.sleep(0.05)
                        readable_stoichiometry = re.sub(r'\b' + cpd_id + r'\b', cpd_name, readable_stoichiometry)

            reaction_data.append({
                "RxnID": rxn_id, "Reaction name": rxn_name, "EC Number": ec_numbers,
                "Stoichiometry (IDs)": equation_ids, "Readable Stoichiometry": readable_stoichiometry,
                "Reversibility": reversibility, "Evidence/Source": "KEGG"
            })
            time.sleep(0.1)
        except Exception as e:
            print(f"Hata: {rxn_id} - {e}")
            continue

    print("Veri çekme işlemi tamamlandı.\n")
    df = pd.DataFrame(reaction_data)
    
    if not df.empty:
        output_file = f"kegg_reactions_{output_name}.xlsx"
        df.to_excel(output_file, index=False)
        print(f"💾 Sonuçlar kaydedildi: {output_file}")
    else:
        print("⚠️ Veri toplanamadığı için Excel dosyası oluşturulmadı.")
    return df

# ==============================================================================
# 3. ANA AKIŞ (PROGRAMIN BAŞLANGIÇ NOKTASI)
# ==============================================================================

if __name__ == "__main__":
    user_input = input(
        "Lütfen bir yolak adı (örn: glycolysis), \n"
        "bir E. coli gen adı/ID'si (örn: serA veya b2926), \n"
        "VEYA bir protein sekansı girin:\n"
    ).strip()
    
    reaction_ids = None
    if user_input:
        if is_protein_sequence(user_input):
            print("\n--- Protein Sekansı Modunda Çalışılıyor ---")
            gene_id = get_gene_from_sequence(user_input)
            if gene_id:
                reaction_ids = get_reactions_from_gene(gene_id.split(':')[1])

        elif re.match(r"^[a-zA-Z]{3,4}[A-Z]?$", user_input) or re.match(r"^[b-d]\d{4}$", user_input.lower()):
            print("\n--- Gen Modunda Çalışılıyor ---")
            reaction_ids = get_reactions_from_gene(user_input)
        
        else:
            print("\n--- Yolak Modunda Çalışılıyor ---")
            reaction_ids = get_reactions_from_pathway(user_input)

        if reaction_ids:
            safe_filename = re.sub(r'[\\/*?:"<>|]', "", user_input.replace(" ", "_"))[:40]
            fetch_reaction_details(reaction_ids, safe_filename)
        else:
            print("\n⚠️ İşlem tamamlandı ancak veri bulunamadı.")
    else:
        print("Geçerli bir girdi girmediniz.")