import sys
import re  # Regex kütüphanesi
import os  # Dosya yolu için eklendi

try:
    from Bio.KEGG import REST
except ImportError:
    print("Hata: Biopython kütüphanesi bulunamadı.")
    print("Lütfen 'pip install biopython' komutu ile kurun.")
    sys.exit(1)

try:
    # Excel'e kaydetmek için YENİ kütüphane
    import openpyxl
    from openpyxl.styles import Font
except ImportError:
    print("Hata: 'openpyxl' kütüphanesi bulunamadı.")
    print("Excel'e kaydetme özelliği için 'pip install openpyxl' komutu ile kurun.")
    sys.exit(1)


class KEGGClient:
    """
    KEGG veritabanıyla web servisi (API) üzerinden iletişim kurmak
    için kullanılan istemci sınıfı.
    """
    
    def __init__(self):
        """Sınıf başlatıldığında bir mesaj yazdırır."""
        print("[KEGGClient] İstemci başlatıldı. KEGG API'ye bağlanmaya hazır.")

    def find(self, database, query):
        """
        Belirli bir veritabanında anahtar kelime araması yapar ve TÜM sonuçları listeler.
        (Menü 1, 2, 3 için kullanılır)
        """
        print(f"\n[Sorgu] '{database}' veritabanında '{query}' aranıyor...")
        try:
            raw_result = REST.kegg_find(database, query).read()
            if not raw_result:
                print("[Sonuç] Sonuç bulunamadı.")
                return []
            results = []
            for line in raw_result.strip().splitlines():
                parts = line.split('\t')
                if len(parts) >= 2:
                    results.append((parts[0], parts[1]))
            print(f"[Sonuç] {len(results)} adet eşleşme bulundu.")
            return results
        except Exception as e:
            print(f"[Hata] Arama sırasında bir hata oluştu: {e}")
            return []

    def find_first_match(self, database, query):
        """
        Belirli bir veritabanında arama yapar ve bulunan İLK sonucu döndürür.
        (Menü 6 - Detaylı Analiz için kullanılır)
        """
        print(f"\n[Sorgu] '{database}' veritabanında '{query}' için ilk eşleşme aranıyor...")
        try:
            raw_result = REST.kegg_find(database, query).read()
            if not raw_result:
                print("[Sonuç] Eşleşme bulunamadı.")
                return None
            first_line = raw_result.strip().splitlines()[0]
            parts = first_line.split('\t')
            if len(parts) >= 2:
                entry_id = parts[0]
                print(f"[Sonuç] '{query}' için ilk eşleşme: {entry_id}")
                return entry_id
            else:
                print("[Sonuç] Geçerli bir eşleşme bulunamadı.")
                return None
        except Exception as e:
            print(f"[Hata] Arama sırasında bir hata oluştu: {e}")
            return None

    def get_entry(self, entry_id):
        """
        Spesifik bir KEGG ID'sine ait tüm detaylı bilgileri çeker.
        (Menü 4, 5, 6 için kullanılır)
        """
        print(f"\n[Sorgu] '{entry_id}' girdisinin detayları alınıyor...")
        try:
            raw_data = REST.kegg_get(entry_id).read()
            if not raw_data:
                print("[Sonuç] Girdi bulunamadı.")
                return None
            return raw_data
        except Exception as e:
            print(f"[Hata] Veri alınırken bir hata oluştu: {e}")
            return None

    def get_brite_hierarchy(self, brite_id):
        """
        BRITE ID'si için tam hiyerarşiyi çeker.
        (Menü 4 için kullanılır)
        """
        if not brite_id.startswith("br:"):
            brite_id = f"br:{brite_id}"
        return self.get_entry(brite_id)

    def link(self, target_db, source_id):
        """
        Bir KEGG girdisinden başka bir veritabanına olan bağlantıları bulur.
        (Menü 6 - Detaylı Analiz için kullanılır)
        """
        print(f"\n[Link Sorgusu] '{source_id}' girdisinden '{target_db}' veritabanına bağlantılar aranıyor...")
        try:
            raw_result = REST.kegg_link(target_db, source_id).read()
            if not raw_result:
                print("[Sonuç] Bağlantı bulunamadı.")
                return []
            results = []
            for line in raw_result.strip().splitlines():
                parts = line.split('\t')
                if len(parts) == 2:
                    results.append(parts[1])
            print(f"[Sonuç] {len(results)} adet bağlantı bulundu.")
            return results
        except Exception as e:
            print(f"[Hata] Link sorgusu sırasında bir hata oluştu: {e}")
            return []

# --- YARDIMCI PARSE FONKSİYONLARI (Menü 6 için) ---

def parse_pathway_for_reactions(pathway_data):
    """(Bu fonksiyon doğruydu, değiştirilmedi)"""
    if not pathway_data:
        return []
    reactions = []
    current_section = None
    for line in pathway_data.splitlines():
        if line.startswith("REACTION"):
            current_section = "REACTION"
            parts = line[12:].split()
            reactions.extend([r.replace(',', '') for r in parts if r.startswith('R')])
        elif line.startswith(" ") and current_section == "REACTION":
            parts = line.strip().split()
            reactions.extend([r.replace(',', '') for r in parts if r.startswith('R')])
        elif not line.startswith(" "):
            if current_section == "REACTION":
                break
            current_section = None
    return sorted(list(set(reactions)))

def parse_reaction(reaction_data):
    """
    DÜZELTİLMİŞ FONKSİYON (v2)
    Reaksiyon dosyasındaki 'COMPOUND' bölümünü
    araya giren 'PATHWAY', 'ENZYME' gibi bölümlere rağmen
    doğru okuyabilen versiyon.
    """
    if not reaction_data:
        return None, None, None

    equation = None
    reaction_id = None
    name_map = {}
    current_section = None  # Mevcut durumu (state) tutar

    for line in reaction_data.splitlines():
        if line.startswith("///"):
            break  # Dosya sonu
        
        # 1. Durum Belirleme (Yeni Bölüm Başlıkları)
        if line.startswith("ENTRY"):
            current_section = "ENTRY"
            try:
                reaction_id = line.split()[1]
            except IndexError:
                pass  # Hatalı satır
        
        elif line.startswith("EQUATION"):
            current_section = "EQUATION"
            # Satır "EQUATION   " (12 karakter) ile başlar
            equation = line[12:].strip()
        
        elif line.startswith("COMPOUND"):
            current_section = "COMPOUND"
            # Başlık satırının kendisi de veri içerebilir
            parts = line[12:].split(None, 1)
            if len(parts) == 2:
                name_map[parts[0]] = parts[1]
        
        elif line.startswith(" "):
            # 2. Durum Kullanma (Devam Satırları)
            # Eğer bir bölümün devamıysa, mevcut duruma göre işle
            if current_section == "EQUATION" and equation is not None:
                equation += " " + line.strip()
            elif current_section == "COMPOUND":
                parts = line.strip().split(None, 1)
                if len(parts) == 2:
                    name_map[parts[0]] = parts[1]
        
        elif not line.startswith(" ") and line.strip():
            # 3. Durum Sıfırlama (Diğer Bölüm Başlıkları, örn: PATHWAY)
            # Bu satır boşlukla başlamıyor, ama ilgilendiğimiz bir bölüm değil.
            # Durumu sıfırla ki, bir sonraki devam satırını yanlışlıkla işlemesin.
            current_section = None

    # --- Okunabilir denklemi oluştur ---
    readable_equation = equation
    if equation and name_map:
        # ID'leri en uzundan kısaya doğru sırala (C00010'u C0001'den önce işle)
        sorted_ids = sorted(name_map.keys(), key=len, reverse=True)
        
        for c_id in sorted_ids:
            name = name_map[c_id]
            # İsimde olabilecek "; açıklama" kısımlarını temizle
            # örn: "ATP ; Adenosine 5'-triphosphate" -> "ATP"
            clean_name = name.split(';')[0].strip()
            if not clean_name:  # Eğer isim "; " ile başlıyorsa
                 clean_name = name.strip()
                 
            # \b = kelime sınırı (tam eşleşme için)
            readable_equation = re.sub(r'\b' + re.escape(c_id) + r'\b', f"[{clean_name}]", readable_equation)

    return reaction_id, equation, readable_equation

# --- YENİ EXCEL KAYDETME FONKSİYONU ---

def save_to_excel(filename, gene_id, pathway_ids, enzyme_ids, reaction_details):
    """
    Analiz sonuçlarını bir Excel dosyasına kaydeder.
    İki sayfa oluşturur: 'Ozet' ve 'Reaksiyonlar'.
    """
    print(f"\n[Excel] Sonuçlar '{filename}' dosyasına kaydediliyor...")
    try:
        wb = openpyxl.Workbook()
        
        # --- Özet Sayfası ---
        ws_summary = wb.active
        ws_summary.title = "Ozet"
        
        # Başlıklar
        ws_summary["A1"] = "Analiz Parametresi"
        ws_summary["B1"] = "Değer"
        ws_summary["A1"].font = Font(bold=True)
        ws_summary["B1"].font = Font(bold=True)
        
        # Veriler
        ws_summary["A2"] = "Aranan Gen ID"
        ws_summary["B2"] = gene_id
        
        ws_summary["A3"] = "Bulunan Pathway'ler"
        # Hücreye sığması için ', ' ile birleştir
        ws_summary["B3"] = ", ".join(pathway_ids) if pathway_ids else "Yok"
        
        ws_summary["A4"] = "Bulunan Enzimler (EC)"
        ws_summary["B4"] = ", ".join(enzyme_ids) if enzyme_ids else "Yok"
        
        ws_summary["A5"] = "Toplam Unique Reaksiyon"
        ws_summary["B5"] = len(reaction_details)
        
        # Sütun genişlikleri
        ws_summary.column_dimensions['A'].width = 25
        ws_summary.column_dimensions['B'].width = 80
        
        # Reaksiyon varsa ikinci sayfayı oluştur
        if reaction_details:
            ws_reactions = wb.create_sheet(title="Reaksiyonlar")
            
            # Başlıklar
            headers = ["Reaksiyon ID", "STO Denklemi (STI)", "Okunabilir Denklem"]
            ws_reactions.append(headers)
            for col in ['A', 'B', 'C']:
                ws_reactions[f"{col}1"].font = Font(bold=True)

            # Veriler
            for rxn in reaction_details:
                ws_reactions.append([rxn['id'], rxn['sti'], rxn['readable']])
                
            ws_reactions.column_dimensions['A'].width = 15
            ws_reactions.column_dimensions['B'].width = 60
            ws_reactions.column_dimensions['C'].width = 80
        
        wb.save(filename)
        print(f"[Excel] Başarıyla kaydedildi: {os.path.abspath(filename)}")
    except PermissionError:
        print(f"[Hata] Excel dosyası kaydedilemedi. '{filename}' dosyası açık mı?")
    except Exception as e:
        print(f"[Hata] Excel dosyası kaydedilirken bir hata oluştu: {e}")


# --- ANA İŞ AKIŞI FONKSİYONLARI ---

def run_gene_analysis(kegg_istemcisi):
    """
    EXCEL'E KAYDETMEK İÇİN GÜNCELLENDİ
    Kullanıcıdan organizma ve gen alıp tam bir analiz yürüten
    ve sonuçları Excel'e kaydeden ana fonksiyon.
    """
    
    print("\n--- Detaylı Gen Analizi Başlatılıyor ---")
    org_code = input("Analiz edilecek ORGANİZMA KODUNU girin (örn: 'hsa', 'eco', 'ecj'): ").strip().lower()
    if not org_code:
        print("Organizma kodu girilmedi.")
        return

    gene_query = input(f"'{org_code}' organizmasındaki GEN ADINI/ID'SİNİ girin (örn: 'tp53', 'dnaA', 'JW0001'): ").strip()
    if not gene_query:
        print("Gen adı girilmedi.")
        return

    # --- ADIM 1: Gen ID'sini Bul ---
    gene_id = kegg_istemcisi.find_first_match(org_code, gene_query)
    if not gene_id:
        print(f"'{org_code}' içinde '{gene_query}' için bir gen bulunamadı.")
        return
    print(f"\n--- {gene_id} (ID: {gene_id}) İÇİN ANALİZ BAŞLADI ---")
    
    # Excel dosya adını belirle
    excel_filename = f"{gene_id.replace(':', '_')}_analizi.xlsx"

    # --- ADIM 2: Gene Bağlı Pathway'leri Bul (Bilgi amaçlı) ---
    print("\n[ADIM 2] Gene bağlı pathway'ler aranıyor...")
    pathway_ids = kegg_istemcisi.link("pathway", gene_id)
    if not pathway_ids:
        print(f"'{gene_id}' genine bağlı pathway bulunamadı.")
    else:
        print(f"Toplam {len(pathway_ids)} pathway bulundu:\n{pathway_ids}")

    # --- ADIM 3: Gene Bağlı Enzimleri Bul (DOĞRU MANTIK) ---
    print("\n[ADIM 3] Genin kodladığı enzimler (EC kodları) aranıyor...")
    enzyme_ids = [eid.replace("ec:", "") for eid in kegg_istemcisi.link("enzyme", gene_id)]
    if not enzyme_ids:
        print(f"\n--- '{gene_id}' genine bağlı bir Enzim (EC) kodu bulunamadı. ---")
        print("(Not: Bu gen, bir enzimi kodlamak yerine regülatör bir gen olabilir.)")
        print("\n--- ANALİZ TAMAMLANDI (Enzim/Reaksiyon Yok) ---")
        # Excel'e boş da olsa kaydet
        save_to_excel(excel_filename, gene_id, pathway_ids, enzyme_ids, [])
        return
    print(f"Toplam {len(enzyme_ids)} enzim kodu bulundu: {enzyme_ids}")

    # --- ADIM 4: Enzimlere Bağlı Reaksiyonları Bul (DOĞRU MANTIK) ---
    all_reactions = set()
    print("\n[ADIM 4] Enzimlerin katalizlediği reaksiyonlar aranıyor...")
    for i, ec_id in enumerate(enzyme_ids, 1):
        print(f"  [4.{i}] '{ec_id}' enzimi için reaksiyonlar:")
        full_ec_id = f"ec:{ec_id}"
        reaction_ids_for_this_ec = [rid.replace("rn:", "") for rid in kegg_istemcisi.link("reaction", full_ec_id)]
        if reaction_ids_for_this_ec:
            print(f"    -> {len(reaction_ids_for_this_ec)} reaksiyon bulundu: {reaction_ids_for_this_ec}")
            all_reactions.update(reaction_ids_for_this_ec)
        else:
            print(f"    -> Bu enzime bağlı reaksiyon bulunamadı.")
    
    if not all_reactions:
        print(f"\n--- Bulunan enzimlere bağlı hiçbir reaksiyon bulunamadı. ---")
        print("\n--- ANALİZ TAMAMLANDI (Reaksiyonsuz) ---")
        save_to_excel(excel_filename, gene_id, pathway_ids, enzyme_ids, [])
        return
        
    print(f"\n--- Toplam {len(all_reactions)} unique reaksiyon bulundu. Detaylar alınıyor... ---")
    
    # Excel için reaksiyon detaylarını sakla
    reaction_details_list = []

    # --- ADIM 5: Reaksiyon Detaylarını ve Denklemleri Al ---
    for i, reaction_id in enumerate(sorted(list(all_reactions)), 1):
        print(f"\n[ADIM 5.{i}] Reaksiyon '{reaction_id}' detayları:")
        reaction_data = kegg_istemcisi.get_entry(reaction_id)
        
        # DÜZELTİLMİŞ PARSER'I ÇAĞIR
        r_id, equation, readable_equation = parse_reaction(reaction_data)
        
        if equation:
            print(f"  > STO Denklemi (STI): {equation}")
            print(f"  > Okunabilir Denklem: {readable_equation}")
            # Excel listesine ekle
            reaction_details_list.append({
                'id': r_id,
                'sti': equation,
                'readable': readable_equation
            })
        else:
            print(f"  -> '{reaction_id}' için denklem verisi alınamadı.")

    print("\n--- ANALİZ TAMAMLANDI ---")
    
    # --- ADIM 6: Excel'e Kaydet ---
    save_to_excel(excel_filename, gene_id, pathway_ids, enzyme_ids, reaction_details_list)


def main_menu():
    """
    Kullanıcıya seçenekler sunan ve KEGGClient sınıfını kullanan
    ana menü fonksiyonu.
    """
    kegg_istemcisi = KEGGClient()
    while True:
        print("\n--- KEGG Veritabanı Arama Aracı ---")
        print("1. Hastalık Ara (Liste)")
        print("2. İnsan Geni Ara (Liste)")
        print("3. İlaç Ara (Liste)")
        print("4. BRITE Hiyerarşisi Al (örn: 'br08902')")
        print("5. Spesifik bir ID'nin Detayını Getir (örn: 'H00056')")
        print("6. Detaylı Gen Analizi (Pathway & Reaksiyon) [Excel'e Kaydeder]")
        print("0. Çıkış")
        
        choice = input("Seçiminiz [1-6, 0]: ")
        
        if choice == '1':
            query = input("Aranacak hastalık adı: ").strip()
            if not query: continue
            results = kegg_istemcisi.find("disease", query)
            for entry_id, description in results:
                print(f"  {entry_id} -> {description}")
                
        elif choice == '2':
            query = input("Aranacak insan geni adı (hsa): ").strip()
            if not query: continue
            results = kegg_istemcisi.find("hsa", query)
            for entry_id, description in results:
                print(f"  {entry_id} -> {description}")
                
        elif choice == '3':
            query = input("Aranacak ilaç adı: ").strip()
            if not query: continue
            results = kegg_istemcisi.find("drug", query)
            for entry_id, description in results:
                print(f"  {entry_id} -> {description}")
                
        elif choice == '4':
            query = input("Getirilecek BRITE ID'si (örn: br08902): ").strip()
            if not query: continue
            data = kegg_istemcisi.get_brite_hierarchy(query)
            if data:
                print(data[:2000] + "\n... (ve daha fazlası)")
                
        elif choice == '5':
            query = input("Detayı istenen tam KEGG ID'si: ").strip()
            if not query: continue
            data = kegg_istemcisi.get_entry(query)
            if data:
                print(data)
        
        elif choice == '6':
            run_gene_analysis(kegg_istemcisi)

        elif choice == '0':
            print("Çıkış yapılıyor...")
            break
            
        else:
            print("Geçersiz seçim. Lütfen tekrar deneyin.")

# Programı çalıştır
if __name__ == "__main__":
    try:
        main_menu()
    except KeyboardInterrupt:
        print("\nProgram kullanıcı tarafından sonlandırıldı.")
        sys.exit(0)