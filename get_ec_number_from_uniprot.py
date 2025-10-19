import requests
import pandas as pd
import sys
import os
import time
import re  # EC numaralarını doğrulamak için 'Regular Expressions' modülü eklendi

# --- KEGG API Fonksiyonları (v2'den - Değişiklik yok) ---

def get_reactions_from_kegg(ec_number):
    """
    Verilen bir EC numarası için KEGG veritabanından reaksiyon ID'lerini çeker.
    HTTPS kullanır.
    """
    url = f"https://rest.kegg.jp/link/reaction/ec:{ec_number}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        
        reaction_ids = []
        for line in response.text.strip().split('\n'):
            if not line: continue
            parts = line.split('\t')
            if len(parts) == 2:
                reaction_id = parts[1].split(':')[1]
                reaction_ids.append(reaction_id)
        return reaction_ids
    except requests.exceptions.RequestException as e:
        print(f"HATA (KEGG link {ec_number}): {e}", file=sys.stderr)
        return []

def get_reaction_details(reaction_id):
    """
    Verilen bir KEGG Reaksiyon ID'si için reaksiyonun tanımını/denklemini çeker.
    HTTPS kullanır.
    """
    url = f"https://rest.kegg.jp/get/{reaction_id}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        
        definition = ""
        equation = ""
        for line in response.text.split('\n'):
            if line.startswith("DEFINITION"):
                definition = line.replace("DEFINITION", "").strip()
            elif line.startswith("EQUATION"):
                equation = line.replace("EQUATION", "").strip()
            if definition and equation: break
                
        return definition if definition else equation if equation else "Tanım/Denklem bulunamadı."
    except requests.exceptions.RequestException as e:
        print(f"HATA (KEGG get {reaction_id}): {e}", file=sys.stderr)
        return "Detay alınamadı."

# --- Dosya Okuma, İşleme ve Excel'e Yazma Fonksiyonları (GÜNCELLENDİ) ---

def read_ec_file(filepath):
    """
    Verilen .txt dosyasını okur ve SADECE GEÇERLİ, TAM EC numaralarının 
    bir listesini döndürür. Kısmi olanları atlar.
    """
    if not os.path.exists(filepath):
        print(f"HATA: Girdi dosyası bulunamadı: {filepath}", file=sys.stderr)
        return None
        
    # Geçerli bir EC numarası formatı: 1.2.3.4 (sadece sayılar)
    # Kısmi olanlar (örn: 1.1.1.-) bu filtreye takılacak.
    ec_pattern = re.compile(r"^\d+\.\d+\.\d+\.\d+$")
    
    ec_numbers = []
    print(f"'{filepath}' dosyası okunuyor...")
    with open(filepath, 'r') as f:
        for line in f:
            cleaned_line = line.strip().replace("EC:", "").replace("ec:", "").strip()
            
            if not cleaned_line or cleaned_line.startswith('#'):
                continue # Boş veya yorum satırlarını atla
                
            # Regex ile formatı kontrol et
            if ec_pattern.match(cleaned_line):
                ec_numbers.append(cleaned_line)
            else:
                print(f"  -> Uyarı: Geçersiz veya kısmi EC no atlanıyor: {line.strip()}")
                
    # Tekrarlanan EC numaralarını temizle
    unique_ec_numbers = sorted(list(set(ec_numbers)))
    return unique_ec_numbers

def process_all_ec_numbers(ec_list):
    """
    EC numarası listesini alır, her biri için reaksiyonları bulur.
    Excel'e yazmak için yapılandırılmış bir liste döndürür.
    """
    print(f"\nToplam {len(ec_list)} GEÇERLİ ve benzersiz EC numarası işlenecek...")
    
    results_list = []
    
    for i, ec in enumerate(ec_list):
        print(f"İşleniyor ({i+1}/{len(ec_list)}): EC {ec}")
        
        # API'ye çok yüklenmemek için her istek arası kısa bir bekleme
        time.sleep(0.1) # 100 milisaniye
        
        reaction_ids = get_reactions_from_kegg(ec)
        if not reaction_ids:
            print(f"  -> EC {ec} için reaksiyon bulunamadı.")
            results_list.append({
                "EC_Number": ec,
                "Reaction_ID": "N/A",
                "Definition": "KEGG'de reaksiyon bulunamadı."
            })
            continue
            
        print(f"  -> {len(reaction_ids)} reaksiyon bulundu. Detaylar alınıyor...")
        for rxn_id in reaction_ids:
            time.sleep(0.1) 
            details = get_reaction_details(rxn_id)
            results_list.append({
                "EC_Number": ec,
                "Reaction_ID": rxn_id,
                "Definition": details
            })

    return results_list

def save_to_excel(results, output_filepath):
    """
    Sonuç listesini bir Pandas DataFrame'e dönüştürür ve Excel'e kaydeder.
    """
    if not results:
        print("Kaydedilecek veri bulunamadı.")
        return

    print(f"\nSonuçlar {output_filepath} dosyasına kaydediliyor...")
    try:
        df = pd.DataFrame(results)
        df.to_excel(output_filepath, index=False, sheet_name="Reaksiyonlar")
        print(f"Başarıyla tamamlandı! {len(df)} kayıt Excel'e yazıldı.")
    except Exception as e:
        print(f"HATA: Excel dosyası yazılırken bir sorun oluştu: {e}", file=sys.stderr)

# --- ANA SCRIPT ---
if __name__ == "__main__":
    
    INPUT_FILE = "ec_number_list.txt"
    OUTPUT_FILE = "e.coli_prokka_reaksiyonlar_sonuc.xlsx"
    
    # 1. EC Numaralarını dosyadan oku (Filtreleyerek)
    ec_list = read_ec_file(INPUT_FILE)
    
    if ec_list is None or not ec_list:
        print("Dosyada işlenecek geçerli EC numarası bulunamadı. Script durduruldu.")
        sys.exit(1)
        
    # 2. Tüm geçerli EC numaralarını işle ve reaksiyonları al
    reaction_data = process_all_ec_numbers(ec_list)
    
    # 3. Sonuçları Excel'e kaydet
    save_to_excel(reaction_data, OUTPUT_FILE)