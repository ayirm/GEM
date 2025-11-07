import requests
import time
import pandas as pd
import re # Metin ayrıştırma (parsing) için
import sys

# KEGG API'sinin temel URL'si
KEGG_API_URL = "http://rest.kegg.jp"

def search_kegg_all(query):
    """
    Verilen bir sorgu terimi için KEGG veritabanlarında kapsamlı bir arama yapar
    ve ham metin kayıtlarını bir sözlük olarak döndürür.
    (Bu fonksiyon değişmedi)
    """
    
    print(f"KEGG'de '{query}' için kapsamlı arama başlatılıyor...")
    
    # 1. Aşama: 'find'
    databases_to_search = ['pathway', 'module', 'compound', 'drug', 'disease', 'ko', 'genes', 'reaction'] # 'reaction' eklendi
    initial_ids = set()
    print(f"Aşama 1: '{query}' terimi şu veritabanlarında aranıyor: {', '.join(databases_to_search)}")
    for db in databases_to_search:
        try:
            url = f"{KEGG_API_URL}/find/{db}/{query}"
            response = requests.get(url)
            response.raise_for_status() 
            lines = response.text.strip().split('\n')
            for line in lines:
                if line:
                    entry_id = line.split('\t')[0]
                    initial_ids.add(entry_id)
            time.sleep(0.1) 
        except requests.exceptions.RequestException as e:
            print(f"Uyarı: '{db}' veritabanı aranırken hata oluştu: {e}")

    if not initial_ids:
        print("İlk arama sonucunda hiçbir kayıt bulunamadı.")
        return {}

    print(f"Aşama 1 tamamlandı. {len(initial_ids)} adet birincil KEGG ID bulundu.")
    
    # 2. Aşama: 'link'
    all_related_ids = set(initial_ids)
    link_map = {
        'pathway': ['gene', 'compound', 'ko', 'module', 'reaction'],
        'module': ['gene', 'compound', 'ko', 'pathway', 'reaction'],
        'compound': ['pathway', 'module', 'drug', 'disease', 'reaction'],
        'drug': ['compound', 'pathway', 'disease'],
        'disease': ['pathway', 'gene', 'drug', 'compound'],
        'ko': ['pathway', 'module', 'gene', 'reaction'],
        'genes': ['pathway', 'module', 'ko', 'disease', 'reaction'],
        'reaction': ['pathway', 'module', 'ko', 'gene', 'compound']
    }
    
    print("Aşama 2: Bulunan ID'ler arası bağlantılar (linkler) aranıyor...")
    for entry_id in initial_ids:
        try:
            db, local_id = entry_id.split(':', 1)
            if db in link_map:
                for target_db in link_map[db]:
                    url = f"{KEGG_API_URL}/link/{target_db}/{entry_id}"
                    response = requests.get(url)
                    response.raise_for_status()
                    lines = response.text.strip().split('\n')
                    for line in lines:
                        if line:
                            linked_id = line.split('\t')[1]
                            all_related_ids.add(linked_id)
                    time.sleep(0.1)
        except ValueError: pass 
        except requests.exceptions.RequestException as e:
            print(f"Uyarı: '{entry_id}' linklenirken hata: {e}")

    print(f"Aşama 2 tamamlandı. Bağlantılarla birlikte toplam {len(all_related_ids)} adet KEGG ID bulundu.")

    # 3. Aşama: 'get'
    results = {}
    print(f"Aşama 3: Toplam {len(all_related_ids)} adet kaydın tam içeriği çekiliyor...")
    id_list = list(all_related_ids)
    for i in range(0, len(id_list), 10):
        chunk = id_list[i:i+10]
        query_str = "+".join(chunk)
        try:
            url = f"{KEGG_API_URL}/get/{query_str}"
            response = requests.get(url)
            response.raise_for_status()
            entries = response.text.strip().split('///\n')
            
            for entry_text in entries:
                entry_text = entry_text.strip()
                if not entry_text: continue 

                first_line = entry_text.split('\n')[0]
                current_id = None
                for an_id in chunk: 
                     id_to_check = an_id.split(':')[-1]
                     if id_to_check in first_line:
                        current_id = an_id
                        break
                
                if not current_id:
                    line_parts = first_line.split()
                    if len(line_parts) > 1:
                        current_id = line_parts[1] 
                    else:
                        print(f"Uyarı: Tanımlanamayan bir kayıt metni bulundu, atlanıyor: '{first_line}...'")
                        continue 

                results[current_id] = entry_text.strip()
            
            print(f"  {min(i + len(chunk), len(id_list))} / {len(id_list)} kayıt çekildi...")
            time.sleep(0.2)
        except requests.exceptions.RequestException as e:
            print(f"Uyarı: ID grubu çekilirken hata: {e}")

    print("Tüm aşamalar tamamlandı.")
    return results

def get_db_type(kegg_id):
    """KEGG ID'sine bakarak veritabanı türünü (Excel sayfa adını) tahmin eder."""
    if ':' in kegg_id:
        prefix = kegg_id.split(':')[0]
        db_map = {
            'map': 'Pathways', 'ko': 'KOs', 'md': 'Modules',
            'cpd': 'Compounds', 'dr': 'Drugs', 'ds': 'Diseases',
            'rn': 'Reactions', 'ec': 'Enzymes' # ec eklendi
        }
        if prefix in db_map: return db_map[prefix]
        return 'Genes' # hsa, eco, mmu vb.
    
    if re.match(r'^C\d{5}$', kegg_id): return 'Compounds'
    if re.match(r'^D\d{5}$', kegg_id): return 'Drugs'
    if re.match(r'^H\d{5}$', kegg_id): return 'Diseases'
    if re.match(r'^K\d{5}$', kegg_id): return 'KOs'
    if re.match(r'^M\d{5}$', kegg_id): return 'Modules'
    if re.match(r'^R\d{5}$', kegg_id): return 'Reactions'
    
    return 'Other'

# --- BU FONKSİYON İSTEĞİNİZE GÖRE GÜNCELLENDİ ---
def parse_kegg_record(record_text):
    """
    KEGG'den gelen tam ham metni Gelişmiş olarak ayrıştırır.
    İlgili alanları (NAME, PATHWAY, ORGANISM, ENZYME, EQUATION vb.) bulup bir sözlüğe atar.
    """
    if not isinstance(record_text, str):
        return {}

    # 1. Tüm kaydı anahtar kelimelere göre böl
    parsed_fields = {}
    current_key = None
    for line in record_text.split('\n'):
        if line.startswith('   '): # Önceki satırın devamı
            if current_key:
                parsed_fields[current_key].append(line.strip())
        else: # Yeni anahtar kelime
            parts = line.split(None, 1) # Sadece ilk boşluktan böl
            if len(parts) == 2 and parts[0].isupper():
                current_key = parts[0]
                parsed_fields[current_key] = [parts[1].strip()]
    
    # 2. Bölünen alanları temizle ve tek metin haline getir
    cleaned_data = {}
    for key, value_list in parsed_fields.items():
        cleaned_data[key] = " ".join(value_list)
        
    # 3. Temizlenen verilerden spesifik ID ve İsimleri çek (Gelişmiş)
    final_parsed = {}
    
    # Temel Bilgiler (Tüm kayıtlarda ortak olabilir)
    if 'ENTRY' in cleaned_data:
        final_parsed['KEGG_ID'] = cleaned_data['ENTRY'].split()[0]
    if 'NAME' in cleaned_data:
        # 'NAME' alanı Gen Adı, Reaksiyon Adı, Yolak Adı vb. olabilir
        final_parsed['NAME'] = cleaned_data['NAME']
    if 'DEFINITION' in cleaned_data:
        final_parsed['DEFINITION'] = cleaned_data['DEFINITION']
    
    # Gen/Organizmaya özel
    if 'ORGANISM' in cleaned_data:
        org_match = re.search(r'(\w{3,5})\s+(.*?)(?:\s\[|$)', cleaned_data['ORGANISM'])
        if org_match:
            final_parsed['Organism_Code'] = org_match.group(1)
            final_parsed['Organism_Name'] = org_match.group(2)
        else:
            final_parsed['ORGANISM_Info'] = cleaned_data['ORGANISM']

    # Reaksiyona özel (İsteğinizdeki EC, Stoichiometry, Reversibility)
    if 'EQUATION' in cleaned_data:
        eq_str = cleaned_data['EQUATION']
        final_parsed['Readable_Stoichiometry'] = eq_str
        if '<=>' in eq_str:
            final_parsed['Reversibility'] = 'Reversible'
        elif '=>' in eq_str:
            final_parsed['Reversibility'] = 'Irreversible'
        else:
            final_parsed['Reversibility'] = 'N/A'
    
    if 'ENZYME' in cleaned_data:
        # EC Numarası
        final_parsed['EC_Number'] = cleaned_data['ENZYME']

    # Bağlantı Bilgileri (Yolak, Modül, KO, Hastalık vb.)
    if 'PATHWAY' in cleaned_data:
        final_parsed['Pathway_Info'] = cleaned_data['PATHWAY']
        final_parsed['Pathway_IDs'] = ", ".join(re.findall(r'(map\d{5})', cleaned_data['PATHWAY']))

    if 'MODULE' in cleaned_data:
        final_parsed['Module_Info'] = cleaned_data['MODULE']
        final_parsed['Module_IDs'] = ", ".join(re.findall(r'(M\d{5})', cleaned_data['MODULE']))

    if 'ORTHOLOGY' in cleaned_data:
        # KEGG Orthology
        final_parsed['Orthology_Info'] = cleaned_data['ORTHOLOGY']
        final_parsed['KO_IDs'] = ", ".join(re.findall(r'(K\d{5})', cleaned_data['ORTHOLOGY']))

    if 'DISEASE' in cleaned_data:
        final_parsed['Disease_Info'] = cleaned_data['DISEASE']
        final_parsed['Disease_IDs'] = ", ".join(re.findall(r'(H\d{5})', cleaned_data['DISEASE']))

    # Component/Compound/Reaction (Modül kayıtları için)
    if 'COMPOUND' in cleaned_data:
        # Bağlı Bileşikler/Componentler
        final_parsed['Component_Compound_Info'] = cleaned_data['COMPOUND']
        final_parsed['Component_Compound_IDs'] = ", ".join(re.findall(r'(C\d{5})', cleaned_data['COMPOUND']))

    if 'REACTION' in cleaned_data:
        # Modüldeki Reaksiyonlar
        final_parsed['Component_Reaction_Info'] = cleaned_data['REACTION']
        final_parsed['Component_Reaction_IDs'] = ", ".join(re.findall(r'(R\d{5})', cleaned_data['REACTION']))
    
    # Bileşiğe (Compound) özel
    if 'FORMULA' in cleaned_data:
        final_parsed['FORMULA'] = cleaned_data['FORMULA']
    if 'EXACT_MASS' in cleaned_data:
        final_parsed['EXACT_MASS'] = cleaned_data['EXACT_MASS']
    
    return final_parsed

# --- BU FONKSİYON YENİ SÜTUNLARI İÇERECEK ŞEKİLDE GÜNCELLENDİ ---
def save_and_parse_to_excel(results_dict, filename):
    """
    Sonuçları alır, HER BİR kaydı Gelişmiş parse_kegg_record ile ayrıştırır
    ve Excel'e farklı sayfalara kaydeder.
    """
    print(f"\nSonuçlar ayrıştırılıyor ve '{filename}' dosyasına aktarılıyor...")
    
    try:
        with pd.ExcelWriter(filename, engine='openpyxl') as writer:
            
            sheets_data = {} # örn: {'Pathways': [dict, dict], 'Genes': [dict, dict]}

            for kegg_id, full_text in results_dict.items():
                
                sheet_name = get_db_type(kegg_id)
                parsed_data = parse_kegg_record(full_text)
                parsed_data['Original_Search_ID'] = kegg_id # Orijinal ID'yi de ekle
                
                if sheet_name not in sheets_data:
                    sheets_data[sheet_name] = []
                sheets_data[sheet_name].append(parsed_data)

            if not sheets_data:
                print("Excel'e yazılacak veri bulunamadı.")
                return

            print(f"Toplam {len(sheets_data)} farklı kategori (sayfa) bulundu.")
            
            for sheet_name, data_list in sheets_data.items():
                print(f"  - '{sheet_name}' sayfası yazılıyor ({len(data_list)} kayıt)...")
                
                df = pd.DataFrame(data_list)
                
                # --- YENİ SÜTUN SIRALAMASI ---
                # İstediğiniz tüm alanları içeren tercih edilen sütun sırası
                preferred_cols = [
                    'Original_Search_ID', 'KEGG_ID', 'NAME', 'DEFINITION', 
                    'EC_Number', 'Readable_Stoichiometry', 'Reversibility', # Reaksiyon
                    'Organism_Code', 'Organism_Name', # Gen
                    'FORMULA', 'EXACT_MASS', # Bileşik
                    'KO_IDs', 'Orthology_Info', # KO
                    'Pathway_IDs', 'Pathway_Info', # Yolak
                    'Module_IDs', 'Module_Info', # Modül
                    'Disease_IDs', 'Disease_Info', # Hastalık
                    'Component_Compound_IDs', 'Component_Compound_Info', # Modül Component
                    'Component_Reaction_IDs', 'Component_Reaction_Info'  # Modül Component
                ]
                
                all_cols = list(df.columns)
                ordered_cols = [col for col in preferred_cols if col in all_cols]
                other_cols = [col for col in all_cols if col not in ordered_cols]
                df = df[ordered_cols + other_cols] # Yeniden sıralanmış DF
                
                df.to_excel(writer, sheet_name=sheet_name, index=False)
        
        print("\n--- İŞLEM BAŞARILI! ---")
        print(f"Tüm sonuçlar hem çekildi, hem ayrıştırıldı hem de '{filename}' dosyasına kaydedildi.")

    except PermissionError:
        print(f"\nHATA: '{filename}' dosyasına yazılamadı. Dosya başka bir programda açık olabilir.")
    except Exception as e:
        print(f"\nExcel'e yazma sırasında beklenmedik bir hata oluştu: {e}")

# --- ANA ÇALIŞTIRMA BLOĞU (Değişmedi) ---

if __name__ == "__main__":
    try:
        user_query = input("Lütfen KEGG'de aramak istediğiniz terimi girin (örn: aspirin, glycolysis, R00200): ")
        
        if not user_query:
            print("Geçerli bir terim girmediniz. Çıkılıyor.")
            sys.exit()
            
        print(f"\n--- '{user_query}' ARAMASI BAŞLATILIYOR ---")
        
        results = search_kegg_all(user_query)
        
        if results:
            print(f"\n'{user_query}' için toplam {len(results)} ilişkili kayıt bulundu.")
            
            safe_filename = re.sub(r'[\\/*?:"<>|]', "", user_query)
            excel_filename = f"{safe_filename}_kegg_ayrintili_sonuclar.xlsx"
            
            save_and_parse_to_excel(results, excel_filename)
                
        else:
            print(f"'{user_query}' için hiçbir sonuç bulunamadı.")
                
    except KeyboardInterrupt:
        print("\nİşlem kullanıcı tarafından iptal edildi.")
    except Exception as e:
        print(f"\nAna programda beklenmedik bir hata oluştu: {e}")