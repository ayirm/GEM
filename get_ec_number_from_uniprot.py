import requests
import pandas as pd
import sys
import os
import time
import re

# --- Global Liste (Yeni Sayfa için) ---
# Bulunan tüm (KO, ReaksiyonID) çiftlerini depolamak için
g_ko_reaction_pairs = set()

# --- KEGG API Fonksiyonları (v7'den - Değişiklik yok) ---

def get_reactions_from_kegg(ec_number):
    """(YÖNTEM 1: Doğrudan EC -> Reaksiyon)"""
    url = f"https://rest.kegg.jp/link/reaction/ec:{ec_number}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        reaction_ids = [line.split('\t')[1].split(':')[1] for line in response.text.strip().split('\n') if line and len(line.split('\t')) == 2]
        return reaction_ids
    except requests.exceptions.RequestException as e:
        print(f"  -> HATA (KEGG link {ec_number}): {e}", file=sys.stderr)
        return []

def get_reaction_details(reaction_id):
    """'R' ID'si olan reaksiyonun detaylarını çeker."""
    url = f"https://rest.kegg.jp/get/{reaction_id}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        definition, equation = "", ""
        for line in response.text.split('\n'):
            if line.startswith("DEFINITION"):
                definition = line.replace("DEFINITION", "").strip()
            elif line.startswith("EQUATION"):
                equation = line.replace("EQUATION", "").strip()
            if definition and equation: break
        return definition if definition else equation if equation else "Tanım/Denklem bulunamadı."
    except requests.exceptions.RequestException as e:
        print(f"  -> HATA (KEGG get {reaction_id}): {e}", file=sys.stderr)
        return "Detay alınamadı."

def get_ko_from_kegg(ec_number):
    """(YÖNTEM 2 - Adım 1: EC -> KO)"""
    url = f"https://rest.kegg.jp/link/ko/ec:{ec_number}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        ko_ids = [line.split('\t')[1].split(':')[1] for line in response.text.strip().split('\n') if line and len(line.split('\t')) == 2]
        return list(set(ko_ids))
    except requests.exceptions.RequestException as e:
        print(f"  -> HATA (KEGG KO link {ec_number}): {e}", file=sys.stderr)
        return []

def get_reactions_from_ko(ko_id):
    """(YÖNTEM 2 - Adım 2: KO -> Reaksiyon)"""
    url = f"https://rest.kegg.jp/link/reaction/ko:{ko_id}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        reaction_ids = [line.split('\t')[1].split(':')[1] for line in response.text.strip().split('\n') if line and len(line.split('\t')) == 2]
        
        # --- YENİ SAYFA İÇİN VERİ TOPLAMA ---
        if not reaction_ids:
            # Reaksiyonu olmasa bile KO'yu listeye ekle
            g_ko_reaction_pairs.add((ko_id, "N/A"))
        else:
            for rxn_id in reaction_ids:
                g_ko_reaction_pairs.add((ko_id, rxn_id))
        # ------------------------------------
        
        return reaction_ids
    except requests.exceptions.RequestException as e:
        print(f"  -> HATA (KEGG RXN link from KO {ko_id}): {e}", file=sys.stderr)
        return []

# --- Enzim Sayfası Okuma (GÜNCELLENDİ v8 - Hata Düzeltmesi) ---

def parse_enzyme_page(ec_number):
    """
    (YÖNTEM 3 & 4 Birleşik)
    EC ana sayfasını (/get/ec:...) sorgular.
    'Transferred to' (Obsolete) ve 'REACTION'/'DEFINITION' (Fallback) arar.
    """
    url = f"https://rest.kegg.jp/get/ec:{ec_number}"
    new_ec, definition_text, reaction_text = None, "", ""
    
    try:
        response = requests.get(url)
        response.raise_for_status()
        lines = response.text.split('\n')
        
        for i, line in enumerate(lines):
            # YÖNTEM 3: Obsolete (Geçersiz) kontrolü
            if "Transferred to" in line:
                match = re.search(r"(\d+\.\d+\.\d+\.\d+)", line)
                if match:
                    new_ec = match.group(1)
                    print(f"  -> BİLGİ: EC {ec_number} geçersiz. Yeni EC {new_ec} bulundu.")
                    break # Yeni EC bulduysak, bu sayfanın tanımına gerek yok.
            
            # YÖNTEM 4: Reaksiyon/Tanım metni arama
            
            # --- HATA DÜZELTMESİ (v8): 'elif' -> 'if' ---
            # Bazı girişlerde hem REACTION hem DEFINITION olabilir.
            
            if line.startswith("REACTION(IUBMB)"):
                reaction_text = line.replace("REACTION(IUBMB)", "").strip()
                j = i + 1
                while j < len(lines) and lines[j].startswith("            "): # 12 boşluk
                    reaction_text += " " + lines[j].strip()
                    j += 1
                    
            if line.startswith("DEFINITION"):
                definition_text = line.replace("DEFINITION", "").strip()
                j = i + 1
                while j < len(lines) and lines[j].startswith("            "):
                    definition_text += " " + lines[j].strip()
                    j += 1
        
        if new_ec:
            return {'new_ec': new_ec, 'definition': None}
        
        # REACTION(IUBMB) tanımını, DEFINITION'a tercih et
        final_definition = reaction_text if reaction_text else definition_text
        
        if final_definition:
             print(f"  -> YÖNTEM 4: Enzim sayfasından tanım metni bulundu (EC: {ec_number}).")
             return {'new_ec': None, 'definition': final_definition}
        
        return {'new_ec': None, 'definition': None}

    except requests.exceptions.RequestException:
        print(f"  -> BİLGİ: KEGG'de {ec_number} için bir sayfa bulunamadı (404 veya diğer).")
        return {'new_ec': None, 'definition': None}

# --- Dosya Okuma (Değişiklik yok) ---

def read_ec_file(filepath):
    """Dosyayı okur ve geçerli EC numaralarını döndürür."""
    if not os.path.exists(filepath):
        print(f"HATA: Girdi dosyası bulunamadı: {filepath}", file=sys.stderr)
        return None
    ec_pattern = re.compile(r"^\d+\.\d+\.\d+\.\d+$")
    ec_numbers = []
    print(f"'{filepath}' dosyası okunuyor...")
    with open(filepath, 'r') as f:
        for line in f:
            cleaned_line = line.strip().replace("EC:", "").replace("ec:", "").strip()
            if not cleaned_line or cleaned_line.startswith('#'): continue
            if ec_pattern.match(cleaned_line):
                ec_numbers.append(cleaned_line)
            else:
                print(f"  -> Uyarı: Geçersiz veya kısmi EC no atlanıyor: {line.strip()}")
    return sorted(list(set(ec_numbers)))

# --- Ana İşleme Fonksiyonları (v7'den - Değişiklik yok) ---

def process_single_ec(ec, original_ec=None):
    """Tek bir EC numarasını 4 yöntemle işler."""
    
    current_ec_for_display = ec
    original_ec_for_display = original_ec if original_ec else ec
    results_list = []
    
    # --- YÖNTEM 1: Doğrudan EC -> Reaksiyon ---
    time.sleep(0.1)
    reaction_ids = get_reactions_from_kegg(ec)
    
    if reaction_ids:
        print(f"  -> YÖNTEM 1 (Doğrudan): {len(reaction_ids)} reaksiyon bulundu (EC: {ec})...")
        for rxn_id in reaction_ids:
            time.sleep(0.1)
            details = get_reaction_details(rxn_id)
            results_list.append({
                "Original_EC_Number": original_ec_for_display,
                "Searched_EC_Number": current_ec_for_display,
                "Reaction_ID": rxn_id,
                "Definition": details,
                "Find_Method": "Direct (EC->RXN)",
                "Associated_KO": "N/A"
            })
        return results_list # Başarılı, bitti

    # --- YÖNTEM 2: Dolaylı EC -> KO -> Reaksiyon ---
    print(f"  -> YÖNTEM 1 (Doğrudan) başarısız (EC: {ec}). YÖNTEM 2 (Dolaylı) deneniyor...")
    time.sleep(0.1)
    ko_ids = get_ko_from_kegg(ec)
    found_kos_str = ""
    
    if ko_ids:
        print(f"  -> YÖNTEM 2: {len(ko_ids)} KO ID bulundu (EC: {ec}). Reaksiyonlar aranıyor...")
        found_kos_str = ", ".join(ko_ids)
        all_indirect_reaction_ids = set()
        ko_reaction_map = {}
        for ko_id in ko_ids:
            time.sleep(0.1)
            # get_reactions_from_ko, global listeyi (g_ko_reaction_pairs) doldurur
            reactions_for_this_ko = get_reactions_from_ko(ko_id)
            for rxn_id in reactions_for_this_ko:
                all_indirect_reaction_ids.add(rxn_id)
                if rxn_id not in ko_reaction_map: ko_reaction_map[rxn_id] = []
                ko_reaction_map[rxn_id].append(ko_id)

        if all_indirect_reaction_ids:
            # YÖNTEM 2 BAŞARILI
            print(f"  -> YÖNTEM 2: {len(all_indirect_reaction_ids)} benzersiz reaksiyon bulundu. Detaylar alınıyor...")
            for rxn_id in sorted(list(all_indirect_reaction_ids)):
                time.sleep(0.1)
                details = get_reaction_details(rxn_id)
                associated_kos_str_map = ", ".join(ko_reaction_map.get(rxn_id, []))
                results_list.append({
                    "Original_EC_Number": original_ec_for_display,
                    "Searched_EC_Number": current_ec_for_display,
                    "Reaction_ID": rxn_id,
                    "Definition": details,
                    "Find_Method": "Indirect (EC->KO->RXN)",
                    "Associated_KO": associated_kos_str_map
                })
            return results_list # Başarılı, bitti
        else:
            print(f"  -> YÖNTEM 2: KO ID'leri bulundu (EC: {ec}) ancak ilişkili 'R' reaksiyonu bulunamadı.")
    else:
        print(f"  -> YÖNTEM 2 başarısız (EC: {ec}). KO ID bulunamadı.")


    # --- YÖNTEM 3 & 4: Enzim Sayfasını Tara ---
    print(f"  -> Yöntem 1 & 2 'R' ID'si bulamadı. Yöntem 3/4 (Enzim Sayfası Kontrolü) deneniyor...")
    time.sleep(0.1)
    page_data = parse_enzyme_page(ec)
    
    # YÖNTEM 3 BAŞARILI: Obsolete (Geçersiz)
    if page_data.get('new_ec'):
        print(f"  -> YÖNTEM 3: {ec} -> {page_data['new_ec']} yönlendirmesi bulundu. Yeni EC için işlemler yeniden başlatılıyor...")
        return process_single_ec(page_data['new_ec'], original_ec=original_ec_for_display) 
    
    # YÖNTEM 4 BAŞARILI: Yedek Tanım bulundu
    if page_data.get('definition'):
        results_list.append({
            "Original_EC_Number": original_ec_for_display,
            "Searched_EC_Number": current_ec_for_display,
            "Reaction_ID": "N/A (Metin)",
            "Definition": page_data['definition'],
            "Find_Method": "Fallback (EC->Definition)",
            "Associated_KO": found_kos_str if found_kos_str else "N/A"
        })
        return results_list

    # --- TÜM YÖNTEMLER BAŞARISIZ OLDU ---
    if found_kos_str: # KO bulduysak ama hiçbir reaksiyon/tanım bulamadıysak
        results_list.append({
            "Original_EC_Number": original_ec_for_display,
            "Searched_EC_Number": current_ec_for_display,
            "Reaction_ID": "N/A",
            "Definition": "KO(lar) bulundu ancak ilişkili reaksiyon veya tanım metni yok.",
            "Find_Method": "Indirect (EC->KO) - Failed",
            "Associated_KO": found_kos_str
        })
    else: # Hiçbir şey bulamadıysak
        print(f"  -> YÖNTEM 3/4: {ec} için yönlendirme veya tanım metni bulunamadı. Tüm yöntemler tükendi.")
        results_list.append({
            "Original_EC_Number": original_ec_for_display,
            "Searched_EC_Number": current_ec_for_display,
            "Reaction_ID": "N/A",
            "Definition": "KEGG'de reaksiyon, KO, yönlendirme veya tanım metni bulunamadı.",
            "Find_Method": "Failed",
            "Associated_KO": "N/A"
        })
    return results_list


def process_all_ec_numbers(ec_list):
    """Tüm EC listesini işler."""
    print(f"\nToplam {len(ec_list)} GEÇERLİ ve benzersiz EC numarası işlenecek...")
    final_results = []
    for i, ec in enumerate(ec_list):
        print(f"İşleniyor ({i+1}/{len(ec_list)}): EC {ec}")
        results_for_this_ec = process_single_ec(ec)
        final_results.extend(results_for_this_ec)
    return final_results

# --- Excel'e Kaydetme (GÜNCELLENDİ v8 - Yeni Sayfa Ekleme) ---

def process_ko_sheet_data(ko_pairs):
    """
    Yeni KO sayfası için veriyi işler.
    Tüm Reaksiyon ID'lerinin detaylarını çeker.
    """
    if not ko_pairs:
        return pd.DataFrame(columns=["KO_ID", "Reaction_ID", "Reaction_Definition"])
        
    print("\n--- Yeni KO Sayfası Verisi Hazırlanıyor ---")
    
    # Benzersiz Reaction ID'lerini topla (N/A hariç)
    unique_rxn_ids = set(rxn_id for ko, rxn_id in ko_pairs if rxn_id != "N/A")
    reaction_details_map = {}

    print(f"Bulunan {len(unique_rxn_ids)} benzersiz KO-reaksiyonu için detaylar alınıyor...")
    
    for i, rxn_id in enumerate(list(unique_rxn_ids)):
        if (i+1) % 10 == 0: # Her 10 istekte bir küçük ilerleme bildirimi
            print(f"  -> KO reaksiyon detayı alınıyor ({i+1}/{len(unique_rxn_ids)})...")
        time.sleep(0.1)
        reaction_details_map[rxn_id] = get_reaction_details(rxn_id)

    print("KO reaksiyon detayları alındı.")

    # Son listeyi oluştur
    ko_sheet_data = []
    for ko, rxn_id in sorted(list(ko_pairs)):
        ko_sheet_data.append({
            "KO_ID": ko,
            "Reaction_ID": rxn_id,
            "Reaction_Definition": reaction_details_map.get(rxn_id, "N/A")
        })
    
    return pd.DataFrame(ko_sheet_data)


def save_to_excel(results, output_filepath, ko_pairs_data):
    """
    Sonuçları Excel'e kaydeder.
    SAYFA 1: EC Reaksiyonları
    SAYFA 2: KO Reaksiyonları
    """
    if not results:
        print("Kaydedilecek (Sayfa 1) veri bulunamadı.")
        return

    print(f"\nSonuçlar {output_filepath} dosyasına kaydediliyor...")
    
    try:
        # Ana EC Sonuçları (Sayfa 1)
        df_results = pd.DataFrame(results)
        columns_order = [
            "Original_EC_Number", "Searched_EC_Number", "Reaction_ID", 
            "Definition", "Find_Method", "Associated_KO"
        ]
        final_columns = [col for col in columns_order if col in df_results.columns]
        df_results = df_results[final_columns]
        
        # KO Sonuçları (Sayfa 2)
        df_ko_results = process_ko_sheet_data(ko_pairs_data)

        # Excel Yazıcısı kullanarak iki sayfayı da yaz
        with pd.ExcelWriter(output_filepath, engine='xlsxwriter') as writer:
            df_results.to_excel(writer, sheet_name="Reaksiyonlar (EC-bazlı)", index=False)
            df_ko_results.to_excel(writer, sheet_name="KO-Reaksiyon-Listesi", index=False)
            
            # (Opsiyonel) Sayfa 1'deki sütunları biraz genişlet
            worksheet1 = writer.sheets["Reaksiyonlar (EC-bazlı)"]
            worksheet1.set_column('A:B', 20)
            worksheet1.set_column('C:C', 15)
            worksheet1.set_column('D:D', 70) # Definition sütunu
            worksheet1.set_column('E:E', 25)
            worksheet1.set_column('F:F', 30) # KO sütunu
            
            # (Opsiyonel) Sayfa 2'deki sütunları biraz genişlet
            worksheet2 = writer.sheets["KO-Reaksiyon-Listesi"]
            worksheet2.set_column('A:B', 15)
            worksheet2.set_column('C:C', 70) # Definition sütunu

        print(f"Başarıyla tamamlandı!")
        print(f"  -> Sayfa 1 (EC-bazlı): {len(df_results)} kayıt yazıldı.")
        print(f"  -> Sayfa 2 (KO-Listesi): {len(df_ko_results)} kayıt yazıldı.")
        
    except Exception as e:
        print(f"HATA: Excel dosyası yazılırken bir sorun oluştu: {e}", file=sys.stderr)

# --- ANA SCRIPT (GÜNCELLENDİ v8) ---
if __name__ == "__main__":
    
    print("--- KEGG Reaksiyon Bulucu v8 (Hata Düzeltmeli + KO Sayfalı) ---")
    
    INPUT_FILE = input("Lütfen EC numaralarını içeren .txt dosyasının adını girin (örn: ec_listesi.txt): ")
    OUTPUT_FILE = input("Lütfen sonuçların kaydedileceği Excel dosyasının adını girin (örn: sonuclar.xlsx): ")
    
    if not OUTPUT_FILE.endswith('.xlsx'):
        OUTPUT_FILE += '.xlsx'
        print(f"  -> Çıktı dosyası adı '{OUTPUT_FILE}' olarak ayarlandı.\n")
    
    ec_list = read_ec_file(INPUT_FILE)
    
    if ec_list is None or not ec_list:
        print("Dosyada işlenecek geçerli EC numarası bulunamadı. Script durduruldu.")
        time.sleep(5) 
        sys.exit(1)
        
    # Sayfa 1 verisini işle (EC-bazlı)
    reaction_data = process_all_ec_numbers(ec_list)
    
    # Sonuçları Excel'e kaydet (Artık 2 sayfayı da işleyecek)
    # g_ko_reaction_pairs, process_all_ec_numbers çalışırken doldu
    save_to_excel(reaction_data, OUTPUT_FILE, g_ko_reaction_pairs)
    
    print("\nTüm işlemler bitti.")
    time.sleep(5)