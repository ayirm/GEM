import requests
import re
import pandas as pd
from io import StringIO

# KEGG (Kyoto Encyclopedia of Genes and Genomes) API adresleri
KEGG_API_URL = "http://rest.kegg.jp"

# Redoks reaksiyonlarını, tiplerini ve elektron sayılarını tanımla
# Biyokimyasal olarak, bu kofaktörlerin değişimi 2 elektron transferi ile olur.
REDOX_REACTIONS = {
    # NAD+ / NADH (2 elektron)
    # (Giren, Çıkan): {Detaylar}
    ('C00003', 'C00004'): {'name': 'NAD+ -> NADH', 'type': 'İndirgenme', 'electrons': '+2e-'},
    ('C00004', 'C00003'): {'name': 'NADH -> NAD+', 'type': 'Yükseltgenme', 'electrons': '-2e-'},
    # NADP+ / NADPH (2 elektron)
    ('C00006', 'C00005'): {'name': 'NADP+ -> NADPH', 'type': 'İndirgenme', 'electrons': '+2e-'},
    ('C00005', 'C00006'): {'name': 'NADPH -> NADP+', 'type': 'Yükseltgenme', 'electrons': '-2e-'},
    # FAD / FADH2 (2 elektron)
    ('C00016', 'C00019'): {'name': 'FAD -> FADH2', 'type': 'İndirgenme', 'electrons': '+2e-'},
    ('C00019', 'C00016'): {'name': 'FADH2 -> FAD', 'type': 'Yükseltgenme', 'electrons': '-2e-'},
}
# Hızlı arama için tüm ilgili kofaktör ID'leri (sadece girenler yeterli)
REDOX_IDS = set(key[0] for key in REDOX_REACTIONS.keys())


def get_pathway_reactions(pathway_name):
    """
    Verilen pathway (yol) adına göre KEGG'den reaksiyon ID'lerini (RN) bulur.
    """
    print(f"'{pathway_name}' için pathway ID aranıyor...")
    
    try:
        # 1. Yol adından (örn: glycolysis) yol ID'sini (örn: map00010) bul
        response = requests.get(f"{KEGG_API_URL}/find/pathway/{pathway_name}")
        response.raise_for_status() 

        lines = response.text.strip().split('\n')
        if not lines or not lines[0]:
            print(f"Hata: '{pathway_name}' adında bir yol bulunamadı.")
            return None, None, None

        first_match_line = lines[0]
        parts = first_match_line.split('\t')
        pathway_id = parts[0].split(':')[1]
        full_pathway_name = parts[1] if len(parts) > 1 else pathway_name
        
        print(f"Pathway bulundu: {pathway_id} ({full_pathway_name})")

        # 2. Yol ID'sinden reaksiyon ID'lerini (RN) al
        response = requests.get(f"{KEGG_API_URL}/link/rn/{pathway_id}")
        response.raise_for_status()

        reaction_ids = []
        for line in response.text.strip().split('\n'):
            if line:
                rn_id = line.split('\t')[1].split(':')[1]
                reaction_ids.append(rn_id)
        
        if not reaction_ids:
            print(f"Hata: {pathway_id} için reaksiyon bulunamadı.")
            return pathway_id, full_pathway_name, []

        print(f"Toplam {len(reaction_ids)} reaksiyon bulundu.")
        return pathway_id, full_pathway_name, reaction_ids

    except requests.exceptions.RequestException as e:
        print(f"KEGG API'ye bağlanırken hata: {e}")
        return None, None, []
    except Exception as e:
        print(f"Beklenmedik bir hata oluştu: {e}")
        return None, None, []


def get_reaction_stoichiometry(reaction_id):
    """
    Verilen reaksiyon ID'sine (RN) ait denklemi (stoichiometry) çeker.
    """
    try:
        response = requests.get(f"{KEGG_API_URL}/get/{reaction_id}")
        response.raise_for_status()
        
        for line in response.text.split('\n'):
            if line.startswith("EQUATION"):
                equation_str = line.replace("EQUATION", "").strip()
                return equation_str
        return None 
    except requests.exceptions.RequestException:
        return None

def analyze_electron_exchange(equation_str):
    """
    Verilen reaksiyon denklemini (Stoichiometry IDs) analiz ederek
    elektron alışverişini arar. Bulunan her alışveriş için bir 
    sözlük listesi döndürür (Excel'e yazmak için).
    """
    if not equation_str:
        return []
    
    parts = equation_str.split(' <=> ')
    if len(parts) != 2:
        return [] # Geçerli bir denklem değil

    reactants_str, products_str = parts
    
    # ID'leri bul (örn: C00003)
    reactant_ids = set(re.findall(r'C\d{5}', reactants_str))
    product_ids = set(re.findall(r'C\d{5}', products_str))

    # Hangi redoks kofaktörlerinin reaksiyona girdiğini bul
    reactants_redox = reactant_ids.intersection(REDOX_IDS)
    # Ürünlerdeki tüm redoks ID'lerini de alalım
    products_redox_all = product_ids.intersection(set(key[1] for key in REDOX_REACTIONS.keys()))

    if not reactants_redox:
        return [] # Girenlerde redoks kofaktörü yoksa (hızlı kontrol)
            
    exchange_list = []
    
    # Girenlerdeki her redoks ID'sini kontrol et
    for r_id in reactants_redox:
        # Ürünlerdeki her redoks ID'sini kontrol et
        for p_id in products_redox_all:
            if (r_id, p_id) in REDOX_REACTIONS:
                # Eşleşen bir redoks reaksiyonu bulundu (örn: C00003 -> C00004)
                exchange_list.append(REDOX_REACTIONS[(r_id, p_id)])
                
    return exchange_list


# --- Ana Program ---
if __name__ == "__main__":
    
    # 1. KULLANICIDAN PATHWAY ADINI AL
    pathway_adı = input("Lütfen analiz etmek istediğiniz pathway adını girin (örn: glycolysis, citrate cycle): ")
    
    if not pathway_adı:
        print("Geçerli bir ad girmediniz. Çıkılıyor.")
        exit()

    path_id, full_name, reaction_list = get_pathway_reactions(pathway_adı.strip())
    
    # 3. EXCEL'E YAZMAK İÇİN SONUÇLARI BİRİKTİR
    results_list = []
    
    if path_id and reaction_list:
        print("\n" + "="*50)
        print(f"Pathway: {full_name} ({path_id})")
        print("Elektron Alışverişi (Redoks Reaksiyonları) Analizi:")
        print("="*50)

        found_redox_reactions = 0
        
        for rn_id in reaction_list:
            # Stoichiometry (IDs) bilgisini (denklemi) al
            equation = get_reaction_stoichiometry(rn_id)
            
            if not equation:
                continue
            
            # 2. ELEKTRON ALIŞVERİŞİNİ ANALİZ ET (liste döndürür)
            exchanges = analyze_electron_exchange(equation)
            
            # Eğer bir elektron alışverişi bulunduysa
            if exchanges:
                found_redox_reactions += 1
                print(f"\nReaksiyon ID: {rn_id}")
                print(f"  Denklem: {equation}")
                
                # Bir reaksiyonda birden fazla alışveriş olabilir
                for ex in exchanges:
                    print(f"  Alışveriş: {ex['name']} (Tip: {ex['type']}, Elektron: {ex['electrons']})")
                    
                    # Excel için listeye ekle
                    results_list.append({
                        'Reaksiyon ID': rn_id,
                        'Denklem (Stoichiometry)': equation,
                        'Elektron Alışverişi': ex['name'],
                        'Tip': ex['type'],
                        'Elektron Sayısı': ex['electrons']
                    })

        if found_redox_reactions == 0:
            print(f"\n'{full_name}' yolunda belirtilen kofaktörler (NAD/FAD) üzerinden")
            print("bariz bir elektron alışverişi reaksiyonu bulunamadı.")
        else:
            print(f"\nAnaliz tamamlandı. Toplam {found_redox_reactions} redoks reaksiyonu bulundu.")
            
            # 3. EXCEL DOSYASINA KAYDET
            if results_list:
                try:
                    # Dosya adını pathway ID'sinden ve adından oluştur
                    # Geçersiz dosya adı karakterlerini temizle
                    safe_name = re.sub(r'[^\w\-_\. ]', '_', full_name.split(' / ')[0])
                    excel_filename = f"{path_id}_{safe_name.replace(' ', '_')}_redoks_analizi.xlsx"
                    
                    df = pd.DataFrame(results_list)
                    # Sütun sırasını belirle
                    df = df[['Reaksiyon ID', 'Denklem (Stoichiometry)', 'Elektron Alışverişi', 'Tip', 'Elektron Sayısı']]
                    
                    df.to_excel(excel_filename, index=False, engine='openpyxl')
                    print(f"\nSonuçlar başarıyla '{excel_filename}' dosyasına kaydedildi.")
                    
                except ImportError:
                    print("\n[HATA] Excel'e kaydetmek için 'pandas' ve 'openpyxl' kütüphaneleri gerekli.")
                    print("Lütfen 'pip install pandas openpyxl' komutları ile yükleyin.")
                except Exception as e:
                    print(f"\nExcel dosyasına kaydederken bir hata oluştu: {e}")
                    
    else:
        print(f"'{pathway_adı}' için analiz yapılamadı.")