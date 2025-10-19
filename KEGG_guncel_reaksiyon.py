from bioservices import KEGG
import pandas as pd
import re
import time

def safe_kegg_get(kegg, entry, retries=3, delay=5):
    """Retry mekanizması ile güvenli KEGG GET isteği"""
    for i in range(retries):
        try:
            data = kegg.get(entry)
            if data:
                return data
        except Exception as e:
            print(f"⚠️ {entry} isteğinde hata: {e}, tekrar deneniyor ({i+1}/{retries})...")
            time.sleep(delay)
    print(f"❌ {entry} alınamadı (tüm denemeler başarısız).")
    return None


def get_reaction_data_from_kegg(pathway_name):
    k = KEGG()
    k.TIMEOUT = 90  # daha uzun timeout

    print(f"🔎 KEGG veritabanında '{pathway_name}' için arama yapılıyor...")

    try:
        all_paths_str = k.list("pathway", "eco")
        matches = [p for p in all_paths_str.split("\n") if pathway_name.lower() in p.lower()]
        if not matches:
            print(f"❌ '{pathway_name}' için yolak bulunamadı.")
            return None
        first_match = matches[0]
        path_id = first_match.split("\t")[0].split(":")[-1]
        print(f"✅ Eşleşme bulundu: {first_match}")
    except Exception as e:
        print(f"⚠️ Yolak alınamadı: {e}")
        return None

    map_id = "map" + path_id[3:]
    print(f"'{map_id}' için reaksiyonlar çekiliyor...")

    try:
        reaction_links = k.link("reaction", map_id)
    except Exception as e:
        print(f"⚠️ Reaksiyon bağlantıları alınamadı: {e}")
        return None

    reaction_ids = []
    if reaction_links:
        for line in reaction_links.strip().split("\n"):
            if "\t" in line:
                rid = line.split("\t")[1].split(":")[-1]
                reaction_ids.append(rid)

    print(f"🧪 {len(reaction_ids)} reaksiyon bulundu.\n")

    data = []
    compound_cache = {}

    for idx, rxn_id in enumerate(reaction_ids, start=1):
        print(f"[{idx}/{len(reaction_ids)}] {rxn_id} inceleniyor...")

        rxn_entry = safe_kegg_get(k, rxn_id)
        if not rxn_entry:
            continue

        rxn_name = ""
        ec_numbers = []
        gene_names = []
        ko_numbers = []
        equation = ""
        readable_equation = ""
        reversibility = "N/A"
        pathways_info = []
        components = "N/A"

        # --- Reaction detaylarını çözümle ---
        for line in rxn_entry.split("\n"):
            if line.startswith("NAME"):
                rxn_name = line.replace("NAME", "").strip()
            elif line.startswith("ENZYME"):
                ec_numbers += re.findall(r"\d+\.\d+\.\d+\.(?:\d+|\-)", line)
            elif line.startswith("ORTHOLOGY"):
                ko_numbers += re.findall(r"K\d{5}", line)
            elif line.startswith("PATHWAY"):
                pathways_info.append(line.replace("PATHWAY", "").strip())
            elif line.startswith("EQUATION"):
                equation = line.replace("EQUATION", "").strip()
                reversibility = "Reversible" if "<=>" in equation else "Irreversible"
            elif line.startswith("COMPONENT") or "COMPONENT" in line:
                components = line.split("COMPONENT")[-1].strip()
            elif line.startswith("COMMENT"):
                if components == "N/A":
                    components = line.replace("COMMENT", "").strip()
            elif "REACTION COMPONENT" in line or "REACTION CLASS" in line:
                if components == "N/A":
                    components = line.split()[-1].strip()

        # --- Bileşik isimlerini çözümle ---
        if equation:
            readable_equation = equation
            compound_ids = re.findall(r"C\d{5}", equation)
            for cid in compound_ids:
                if cid not in compound_cache:
                    cpd_entry = safe_kegg_get(k, cid)
                    cpd_name = "Unknown"
                    if cpd_entry:
                        for l in cpd_entry.split("\n"):
                            if l.startswith("NAME"):
                                cpd_name = l.replace("NAME", "").strip().split(";")[0]
                                break
                    compound_cache[cid] = cpd_name
                    time.sleep(0.05)
                readable_equation = re.sub(rf"\b{cid}\b", compound_cache[cid], readable_equation)

        # --- EC üzerinden genleri bul ---
        if ec_numbers:
            for ec in ec_numbers:
                try:
                    eco_link = k.link("eco", f"ec:{ec}")
                    eco_genes = re.findall(r"eco:(\w+)", eco_link)
                    if eco_genes:
                        gene_names.extend(eco_genes)
                        time.sleep(0.05)
                except Exception:
                    continue

        data.append({
            "RxnID": rxn_id,
            "Reaction name": rxn_name or "N/A",
            "EC Number": " ".join(ec_numbers) if ec_numbers else "N/A",
            "Gene Name": " ".join(set(gene_names)) if gene_names else "N/A",
            "KEGG Orthology": " ".join(ko_numbers) if ko_numbers else "N/A",
            "Pathway Name": "; ".join(pathways_info) if pathways_info else "N/A",
            "Components": components if components else "N/A",
            "Stoichiometry (IDs)": equation or "N/A",
            "Readable Stoichiometry": readable_equation or "N/A",
            "Reversibility": reversibility,
            "Evidence/Source": "KEGG"
        })

        time.sleep(0.2)

    df = pd.DataFrame(data)
    if not df.empty:
        out_file = f"kegg_reactions_{pathway_name.replace(' ', '_')}.xlsx"
        df.to_excel(out_file, index=False)
        print(f"\n💾 Kaydedildi: {out_file}")
    else:
        print("⚠️ Veri alınamadı.")

    return df


if __name__ == "__main__":
    pathway_input = input("Lütfen E. coli'de aramak istediğiniz yolağın adını girin (örnek: glycolysis): ").strip()
    if pathway_input:
        get_reaction_data_from_kegg(pathway_input)
