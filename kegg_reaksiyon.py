from bioservices import KEGG
import pandas as pd
import time
import re

def get_reaction_data_from_kegg(pathway_name):
    k = KEGG()
    k.TIMEOUT = 60
    
    # Bileşik ID'lerini isimleriyle eşleştirmek için bir önbellek (cache) oluşturalım.
    # Bu, aynı bileşiği tekrar tekrar sormamızı engelleyerek hızı artırır.
    compound_name_cache = {}

    print(f"KEGG veritabanında '{pathway_name}' için arama yapılıyor...")

    path_id = None
    try:
        # E. coli için tüm yolakları listele ve içinde ara (en sağlam yöntem)
        all_paths_str = k.list("pathway", "eco")
        if not isinstance(all_paths_str, str):
             print(f"⚠️ KEGG'den yolak listesi alınamadı (API hatası).")
             return None

        matches = [p for p in all_paths_str.split("\n") if pathway_name.lower() in p.lower()]
        
        if matches:
            first_match = matches[0]
            print(f"🔎 Eşleşme bulundu: {first_match}")
            path_id = first_match.split("\t")[0].split(":")[-1]
        else:
            print(f"E. coli ('eco') içinde '{pathway_name}' ile eşleşen bir yolak bulunamadı.")
            return None
            
    except Exception as e:
        print(f"⚠️ Yolak listesi alınamadı: {e}")
        return None

    print(f"KEGG yolağı bulundu: {path_id}")

    # Reaksiyonları 'link' fonksiyonu ile çekmek en güvenilir yoldur.
    map_id = "map" + path_id[3:]
    print(f"'{map_id}' haritasıyla ilişkili reaksiyonlar aranıyor...")

    unique_reaction_ids = []
    try:
        reaction_links = k.link("reaction", map_id)
        if isinstance(reaction_links, str) and reaction_links.strip():
            for line in reaction_links.strip().split('\n'):
                parts = line.split('\t')
                if len(parts) == 2:
                    reaction_part = parts[1].split(':')[1]
                    unique_reaction_ids.append(reaction_part)
        else:
            print(f"'{map_id}' haritası için bağlantılı reaksiyon bulunamadı.")
            return None
    except Exception as e:
        print(f"⚠️ Reaksiyon bağlantıları alınırken hata oluştu: {e}")
        return None

    if not unique_reaction_ids:
        print(f"'{map_id}' haritası içinde reaksiyon ID'si bulunamadı.")
        return None

    print(f"Toplam {len(unique_reaction_ids)} benzersiz reaksiyon bulundu. Detaylar çekiliyor...")
    print("UYARI: Okunabilir stokiyometri oluşturmak için ek API sorguları yapılacak, bu işlem yavaş olabilir.")

    reaction_data = []
    for i, rxn_id in enumerate(unique_reaction_ids, 1):
        print(f"İşleniyor: {i}/{len(unique_reaction_ids)} ({rxn_id})")
        try:
            rxn_data_str = k.get(rxn_id)
            if not rxn_data_str: continue

            # Değişkenleri başlangıçta tanımla
            rxn_name, equation_ids, ec_numbers, reversibility = "", "", "", "N/A"
            readable_stoichiometry = ""

            for line in rxn_data_str.split("\n"):
                if line.startswith("NAME"):
                    rxn_name = line.replace("NAME", "").strip()
                elif line.startswith("ENZYME"):
                    # Sadece EC numaralarını bul (örn: 1.2.3.4 veya 1.2.3.-)
                    found_ecs = re.findall(r"\d+\.\d+\.\d+\.(?:\d+|\-)", line)
                    ec_numbers = " ".join(found_ecs)
                elif line.startswith("EQUATION"):
                    equation_ids = line.replace("EQUATION", "").strip()
                    reversibility = "Reversible" if "<=>" in equation_ids else "Irreversible"
                    
                    # Okunabilir stokiyometri oluşturma
                    readable_stoichiometry = equation_ids
                    compound_ids = set(re.findall(r"C\d{5}", equation_ids))
                    
                    for cpd_id in compound_ids:
                        if cpd_id in compound_name_cache:
                            cpd_name = compound_name_cache[cpd_id]
                        else:
                            # Önbellekte yoksa KEGG'den çek
                            cpd_data = k.get(cpd_id)
                            cpd_name = "Unknown"
                            for cpd_line in cpd_data.split('\n'):
                                if cpd_line.startswith("NAME"):
                                    cpd_name = cpd_line.replace("NAME", "").strip().split(';')[0]
                                    break
                            compound_name_cache[cpd_id] = cpd_name
                            time.sleep(0.05) # API'yi yormamak için küçük bekleme
                        
                        # ID'yi isimle değiştir (regex ile tam eşleşme sağlanır)
                        readable_stoichiometry = re.sub(r'\b' + cpd_id + r'\b', cpd_name, readable_stoichiometry)

            reaction_data.append({
                "RxnID": rxn_id,
                "Reaction name": rxn_name,
                "EC Number": ec_numbers,
                "Stoichiometry (IDs)": equation_ids,
                "Readable Stoichiometry": readable_stoichiometry,
                "Reversibility": reversibility,
                "Evidence/Source": "KEGG"
            })
            time.sleep(0.1)

        except Exception as e:
            print(f"Hata: {rxn_id} - {e}")
            continue

    print("Veri çekme işlemi tamamlandı.\n")

    df = pd.DataFrame(reaction_data)
    if not df.empty:
        output_file = f"kegg_reactions_{pathway_name.replace(' ', '_')}.xlsx"
        df.to_excel(output_file, index=False)
        print(f"💾 Sonuçlar kaydedildi: {output_file}")
    else:
        print("⚠️ Veri toplanamadığı için Excel dosyası oluşturulmadı.")
    return df

# Ana akış
if __name__ == "__main__":
    pathway_input = input("Lütfen E. coli'de aramak istediğiniz yolağın adını girin (örneğin, glycolysis): ").strip()
    if pathway_input:
        final_df = get_reaction_data_from_kegg(pathway_input)