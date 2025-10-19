import re
import requests
import pandas as pd
from collections import defaultdict
import time

# --- KEGG API Fonksiyonları ---

KEGG_API_URL = "http://rest.kegg.jp"

def get_pathway_reactions(pathway_name: str) -> tuple[str | None, str | None, list]:
    """
    Verilen pathway adına göre KEGG'den reaksiyon ID'lerini (RNs) bulur.
    """
    print(f"'{pathway_name}' için KEGG'de pathway ID'si aranıyor...")
    try:
        # Pathway adından Pathway ID'sini bul (örn: "glycolysis" -> "map00010")
        response = requests.get(f"{KEGG_API_URL}/find/pathway/{pathway_name}")
        response.raise_for_status()

        lines = response.text.strip().split('\n')
        if not lines or not lines[0]:
            print(f"HATA: '{pathway_name}' adında bir yol bulunamadı.")
            return None, None, []

        first_match_line = lines[0]
        parts = first_match_line.split('\t')
        pathway_id = parts[0].split(':')[1]
        full_pathway_name = parts[1] if len(parts) > 1 else pathway_name
        
        print(f"Pathway bulundu: {pathway_id} ({full_pathway_name})")

        # Pathway ID'sinden reaksiyon ID'lerini al
        response = requests.get(f"{KEGG_API_URL}/link/rn/{pathway_id}")
        response.raise_for_status()

        reaction_ids = [line.split('\t')[1].split(':')[1] for line in response.text.strip().split('\n') if line]
        
        if not reaction_ids:
            print(f"HATA: {pathway_id} için reaksiyon bulunamadı.")
            return pathway_id, full_pathway_name, []

        print(f"Toplam {len(reaction_ids)} reaksiyon ID'si bulundu.")
        return pathway_id, full_pathway_name, reaction_ids

    except requests.exceptions.RequestException as e:
        print(f"KEGG API'ye bağlanırken bir ağ hatası oluştu: {e}")
        return None, None, []

def get_reaction_equation(reaction_id: str) -> str | None:
    """
    Verilen reaksiyon ID'si için KEGG'den tam denklem metnini çeker.
    """
    try:
        response = requests.get(f"{KEGG_API_URL}/get/{reaction_id}")
        response.raise_for_status()
        
        for line in response.text.split('\n'):
            if line.startswith("EQUATION"):
                return line.replace("EQUATION", "").strip()
        return None # Denklem satırı bulunamadı
    except requests.exceptions.RequestException:
        return None

# --- Stokiyometri Hesaplama Fonksiyonu ---

def calculate_stoichiometry(reaction_string: str) -> dict:
    """
    Bir reaksiyon metnini ayrıştırır ve her bileşik için net stokiyometrisini hesaplar.
    """
    if not isinstance(reaction_string, str) or ' <=> ' not in reaction_string:
        raise ValueError("Geçersiz reaksiyon formatı.")

    stoichiometry = defaultdict(int)
    compound_regex = re.compile(r'(\d*)\s*([a-zA-Z0-9_]+)')
    reactants_part, products_part = reaction_string.split(' <=> ')

    # Reaktanları işle (Tüketilenler)
    for match in compound_regex.finditer(reactants_part):
        coefficient = int(match.group(1)) if match.group(1) else 1
        stoichiometry[match.group(2)] -= coefficient

    # Ürünleri işle (Üretilenler)
    for match in compound_regex.finditer(products_part):
        coefficient = int(match.group(1)) if match.group(1) else 1
        stoichiometry[match.group(2)] += coefficient
        
    return dict(stoichiometry)

# --- Ana Program ---
if __name__ == "__main__":
    
    pathway_input = input("Lütfen analiz etmek istediğiniz pathway adını girin (örn: glycolysis, citrate cycle): ")
    if not pathway_input.strip():
        print("Geçerli bir ad girmediniz. Program sonlandırılıyor.")
        exit()

    # Pathway için tüm reaksiyon ID'lerini al
    path_id, full_name, reaction_ids = get_pathway_reactions(pathway_input)

    if not reaction_ids:
        print(f"'{pathway_input}' için analiz edilecek reaksiyon bulunamadı.")
        exit()

    all_reactions_data = []
    total_reactions = len(reaction_ids)
    
    print("\n" + "-" * 45)
    print(f"Toplam {total_reactions} reaksiyonun denklemleri alınıp analiz ediliyor...")

    # Her bir reaksiyon ID'si için döngü başlat
    for i, rn_id in enumerate(reaction_ids):
        # Kullanıcıya ilerlemeyi göster
        print(f"  İşleniyor... [{i+1}/{total_reactions}] {rn_id}", end='\r')

        # Reaksiyonun denklemini al
        equation = get_reaction_equation(rn_id)
        
        # API sunucusunu yormamak için çok kısa bir bekleme ekle
        time.sleep(0.05)
        
        if not equation:
            continue # Denklem bulunamazsa bu reaksiyonu atla

        # Denklem için stokiyometriyi hesapla
        try:
            stoichiometry_result = calculate_stoichiometry(equation)
            
            # Excel satırı için veriyi hazırla
            row_data = {'Reaction ID': rn_id, 'Reaction String': equation}
            row_data.update(stoichiometry_result)
            all_reactions_data.append(row_data)

        except ValueError as e:
            print(f"\n  [HATA] {rn_id} işlenemedi: {e}")

    # Döngü bittikten sonra, eğer veri varsa Excel'e yaz
    if not all_reactions_data:
        print("\nHiçbir reaksiyon için geçerli veri bulunamadı.")
    else:
        print("\n\n" + "-" * 45)
        print(f"Toplam {len(all_reactions_data)} reaksiyon analiz edildi. Excel dosyası oluşturuluyor...")
        
        try:
            df = pd.DataFrame(all_reactions_data).fillna(0)
            
            # ID ve denklem metni dışındaki sütunları tam sayıya çevir
            string_cols = ['Reaction ID', 'Reaction String']
            for col in df.columns:
                if col not in string_cols:
                    df[col] = df[col].astype(int)

            # Sütunları yeniden sırala (ID ve Denklem başa, diğerleri alfabetik)
            cols = sorted([col for col in df.columns if col not in string_cols])
            final_cols = string_cols + cols
            df = df[final_cols]

            # Temiz bir dosya adı oluştur
            safe_name = re.sub(r'[^\w\-_\. ]', '_', full_name.split(' / ')[0]).replace(' ', '_')
            output_filename = f"{path_id}_{safe_name}_stokiyometri.xlsx"
            
            df.to_excel(output_filename, index=False, engine='openpyxl')
            
            print(f"\n✓ Başarılı! Sonuçlar '{output_filename}' dosyasına kaydedildi.")

        except Exception as e:
            print(f"\n[HATA] Dosyaya yazarken bir hata oluştu: {e}")
