import re
import requests
import pandas as pd
import json
import argparse
import time
from collections import defaultdict
import sys

# --- KEGG API FONKSİYONLARI ---
KEGG_API_URL = "http://rest.kegg.jp"

def get_reaction_equation(reaction_id: str) -> str | None:
    """Reaction ID (RNxxxx) için KEGG'den denklemi çeker."""
    try:
        # ID başında 'rn:' yoksa ekle
        full_id = f"rn:{reaction_id}" if not reaction_id.startswith("rn:") else reaction_id
        response = requests.get(f"{KEGG_API_URL}/get/{full_id}")
        response.raise_for_status()
        
        for line in response.text.split('\n'):
            if line.startswith("EQUATION"):
                return line.replace("EQUATION", "").strip()
        return None
    except Exception:
        return None

def calculate_stoichiometry(reaction_string: str) -> dict:
    """Denklemi ayrıştırıp net stokiyometriyi hesaplar."""
    if not isinstance(reaction_string, str) or ' <=> ' not in reaction_string:
        return {}

    stoichiometry = defaultdict(int)
    # Bileşiklerin başındaki katsayıları ve ID'leri yakalar
    compound_regex = re.compile(r'(\d*)\s*([a-zA-Z0-9_]+)')
    reactants_part, products_part = reaction_string.split(' <=> ')

    # Reaktanlar (-)
    for match in compound_regex.finditer(reactants_part):
        coeff = int(match.group(1)) if match.group(1) else 1
        stoichiometry[match.group(2)] -= coeff

    # Ürünler (+)
    for match in compound_regex.finditer(products_part):
        coeff = int(match.group(1)) if match.group(1) else 1
        stoichiometry[match.group(2)] += coeff
        
    return dict(stoichiometry)

# --- ANA MANTIK ---
def main():
    parser = argparse.ArgumentParser(description="JSON çıktısından metabolik matris oluşturur.")
    parser.add_argument("--input", required=True, help="Arkadaşının oluşturduğu .json dosyasının yolu")
    parser.add_argument("--output", default="metabolic_matrix.xlsx", help="Çıktı Excel adı")
    args = parser.parse_args()

    # 1. JSON dosyasını oku ve RXN ID'leri ayıkla
    try:
        with open(args.input, 'r') as f:
            data = json.load(f)
    except Exception as e:
        print(f"HATA: JSON dosyası okunamadı: {e}")
        sys.exit(1)

    # Benzersiz reaksiyon ID'lerini topla
    unique_rxn_ids = set()
    for entry in data:
        # Arkadaşının 'reactions_by_pathway' içindeki string'i veya listeyi çekiyoruz
        rxn_info = entry.get("reactions_by_pathway", "{}")
        
        # Eğer veri string formatında bir JSON ise parse et
        if isinstance(rxn_info, str) and rxn_info.startswith("{"):
            try:
                path_dict = json.loads(rxn_info)
                for rxn_list in path_dict.values():
                    for rid in rxn_list:
                        # 'rn:R00001' -> 'R00001' temizliği
                        clean_id = rid.split(':')[-1]
                        unique_rxn_ids.add(clean_id)
            except: continue
        elif isinstance(rxn_info, dict): # Eğer doğrudan dict ise
            for rxn_list in rxn_info.values():
                for rid in rxn_list:
                    unique_rxn_ids.add(rid.split(':')[-1])

    if not unique_rxn_ids:
        print("Uyarı: JSON içinde hiç reaksiyon ID'si bulunamadı.")
        return

    print(f"Toplam {len(unique_rxn_ids)} benzersiz reaksiyon tespit edildi. Analiz başlıyor...")

    # 2. KEGG'den denklemleri çek ve stokiyometri hesapla
    all_data = []
    for i, rid in enumerate(unique_rxn_ids, 1):
        print(f" İşleniyor [{i}/{len(unique_rxn_ids)}]: {rid}", end='\r')
        eqn = get_reaction_equation(rid)
        time.sleep(0.1) # Rate limit koruması
        
        if eqn:
            stoich = calculate_stoichiometry(eqn)
            if stoich:
                row = {'Reaction_ID': rid, 'Equation': eqn}
                row.update(stoich)
                all_data.append(row)

    # 3. Excel'e aktar
    if all_data:
        df = pd.DataFrame(all_data).fillna(0)
        # Sütunları düzenle (ID ve Eq başa)
        cols = ['Reaction_ID', 'Equation'] + sorted([c for c in df.columns if c not in ['Reaction_ID', 'Equation']])
        df = df[cols]
        
        # Sayısal değerleri int yap (okunabilirlik için)
        for col in df.columns[2:]:
            df[col] = df[col].astype(int)

        df.to_excel(args.output, index=False)
        print(f"\n✓ Başarılı! {args.output} oluşturuldu.")
    else:
        print("\nVeri bulunamadı.")

if __name__ == "__main__":
    main()
