import pandas as pd
from bioservices import KEGG
import sys


def get_reaction_data_from_kegg(pathway_name):
    """
    KEGG veritabanından belirtilen yolağa ait reaksiyon bilgilerini çeker.
    """
    k = KEGG()

    # Zaman aşımı süresini doğrudan ayarla (bazı sürümlerde k.settings yok)
    if hasattr(k, "settings") and k.settings is not None:
        k.settings.TIMEOUT = 60
    else:
        k.TIMEOUT = 60

    print(f"KEGG veritabanında '{pathway_name}' için arama yapılıyor...")

    # E. coli (eco) için yolağı bul
    try:
        pathway_list_str = k.find("pathway", pathway_name)
        if not pathway_list_str:
            print("KEGG yanıt vermedi veya arama sonucu boş döndü.")
            return None
    except Exception as e:
        print(f"KEGG API'sine bağlanırken bir hata oluştu: {e}")
        return None

    eco_pathways = [line for line in pathway_list_str.strip().split("\n") if "eco" in line]

    if not eco_pathways:
        print(f"KEGG'de E. coli için '{pathway_name}' ile eşleşen bir yolağı bulunamadı.")
        return None

    # Genellikle ilk sonuç en alakalı olandır
    pathway_id = eco_pathways[0].split("\t")[0]
    pathway_full_name = eco_pathways[0].split("\t")[1]
    print(f"KEGG yolağı bulundu: {pathway_full_name} ({pathway_id})")

    # Reaksiyonları al
    try:
        reaction_links = k.link("reaction", pathway_id)
        if not reaction_links:
            print("Bu yolağa bağlı reaksiyon bulunamadı.")
            return None
        reaction_ids = [line.split("\t")[1].replace("rn:", "") for line in reaction_links.strip().split("\n")]
    except Exception as e:
        print(f"Reaksiyonlar alınamadı: {e}")
        return None

    print(f"Toplam {len(reaction_ids)} reaksiyon bulundu. Detaylar çekiliyor...")

    all_reactions_data = []
    for i, rxn_id in enumerate(reaction_ids):
        sys.stdout.write(f"\rİşleniyor: {i+1}/{len(reaction_ids)} ({rxn_id})")
        sys.stdout.flush()

        try:
            rxn_data_str = k.get(rxn_id)
            if not rxn_data_str:
                print(f"\n{rxn_id} için veri alınamadı. Atlanıyor.")
                continue
        except Exception as e:
            print(f"\n{rxn_id} için veri çekilemedi: {e}")
            continue

        reaction_name = ""
        ec_number = ""
        stoichiometry = ""
        reversibility = ""
        gpr = ""

        for line in rxn_data_str.strip().split("\n"):
            if line.startswith("NAME"):
                reaction_name = line.replace("NAME", "").strip()
            elif line.startswith("EQUATION"):
                stoichiometry = line.replace("EQUATION", "").strip()
                reversibility = "Reversible" if "<=>" in stoichiometry else "Irreversible"
            elif line.startswith("ENZYME"):
                ec_number = line.replace("ENZYME", "").strip()

        # GPR (Gen-Protein ilişkilendirmesi)
        try:
            gene_links = k.link("genes", f"rn:{rxn_id}")
            if gene_links:
                eco_genes = [line.split("\t")[1] for line in gene_links.strip().split("\n") if line.startswith("eco:")]
                gpr = " or ".join(sorted(list(set(eco_genes)))) if eco_genes else "N/A"
            else:
                gpr = "N/A"
        except Exception:
            gpr = "N/A"

        all_reactions_data.append({
            "RxnID": rxn_id,
            "Reaction name": reaction_name,
            "EC Number": ec_number,
            "Stoichiometry": stoichiometry,
            "Compartment": "N/A",
            "Reversibility": reversibility,
            "GPR": gpr,
            "Evidence/Source": "KEGG"
        })

    print("\nVeri çekme işlemi tamamlandı.")
    return pd.DataFrame(all_reactions_data)


if __name__ == "__main__":
    pathway_input = input("Lütfen E. coli'de aramak istediğiniz yolağın adını girin (örneğin, serine biosynthesis): ")

    final_df = get_reaction_data_from_kegg(pathway_input)

    if final_df is not None and not final_df.empty:
        output_file = f"{pathway_input.replace(' ', '_')}_reaction_data_KEGG.xlsx"
        try:
            final_df.to_excel(output_file, index=False)
            print(f"\nVeriler başarıyla '{output_file}' dosyasına kaydedildi.")
        except Exception as e:
            print(f"\nExcel dosyası kaydedilirken bir hata oluştu: {e}")
    else:
        print("\nKEGG veritabanından herhangi bir veri alınamadı.")
