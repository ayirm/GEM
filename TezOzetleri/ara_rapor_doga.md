# GEM ve Modelleme Projesi
Çalışmanın amacı *E.coli(Escherichia coli)* türündeki antibiyotik direncinin modellenmesi ve bu direnci azaltmaya yönelik yöntemleri bulmaktır.

Modellenme sürecinde ise fastq dosyalarından elde edilen okumalar(reads), *E.coli* referansları[^1] ile birlikte reference guided assembly pipeline'nından geçirelerek consesus sequence oluşturulmuş.
Bu aşamalarda:
- bowtie2
- samtools
- bcftools
- bgzip ve tabix
kullanılmıştır.

İlk olarak referans genom, bowtie ile indekslenmiş sonrasında fastq sonuçları[^2] referans genom ile eşleştirilmiştir(alignment). Bu eşleştirme sonucunda çıkan .sam dosyası önce samtools ile .bam dosyasına çevrilmiş ve ardından kromozomlara göre sıraya dizilmiştir, en son ise indekslenmiştir.

İndekslenmiş .bam dosyası bcftools ile varyansyon bilgileri korunmuş bir okuma haline(collapsed reads) getirildi ve .bcf dosyası olarak kaydedildi. Ardından bu .bcf dosyası .vcf dosya formatına çevrildi.
Oluşan .vcf dosyası önce bgzip ile sıkıştırıldı, sonrasında tabix ile indekslendi.

En son adım olarak bcftools kullanılarak referans dosyası ve okuma dosyası karşılaştırılarak ortak bir genom dosyası(consensus sequence) elde edilmiştir.

Consensus sequence oluştuktan sonra ise prokka ile annotation yapılmıştır.
```
# Çalıştırılan Kodlar

bowtie2-build --threads 4 ref_file ref_name
bowtie2 -x ref_name -1 fastq1 -2 fastq2 -S sam_file

samtools view -uS -o bam_file sam_file
samtools sort -@ 4 -T temp_file -o sorted_bam bam_file
samtools index sorted_bam

bcftools mpileup -f ref_file sorted_bam | bcftools call -mV -Ob -o bcf_file
bcftools convert -O v -o vcf_file bcf_file
bgzip -c vcf_file > gz_file
tabix -p vcf gz_file

bcftools consensus -f ref_file gz_file > consensus_fasta

prokka --outdir outputDir --prefix prokkaAnnotes --force consensus_fasta
```

## Annotation ve Sonraki Adımlar
Annotation sonrasında modelleme için gereken enzimler bulunmuş olsa bile bu enzimlerin yer aldığı reaksiyonlar ve bu reaksiyonların stoichiometric matrix(S) için gereken bilgileri elimizde bulunmamaktadır.

Bu nedenle prokka sonucundan çıkan genbank dosyalarından:
- gen adı
- enzim adı
- uniprot id
bilgileri elde edilmiştir. 
>[!info]
> Genbank dosyasının(.gbk veya .gbf) okunması `KEGG_merge.py` içerisinde `uniprot_search` fonksiyonu ile yapılmaktadır.

Sonrasında KEGG'de reaksiyonları aratabilmek için Uniprot ID'ler, Uniprot'un sitesine(`https://rest.uniprot.org/idmapping`) gönderilerek her bir Uniprot ID'ye karşılık gelen KEGG ID'ler bulunmuştur.
Elde edilen KEGG ID'ler bu sefer KEGG'in web sitesine(`https://rest.kegg.jp`) gönderilerek yolak(pathway) bilgileri elde edilmiştir. Bu işlem sonrasında her bir yolak için tekrarlanmıştır ve her bir yolakta bulunan reaksiyonlar ve bu reaksiyonlardaki moleküller(compounds) bulunmuştur.
>[!info]
>Uniprot ID'lerden KEGG ID'lere geçiş için `KEGG_merge.py` dosyasındaki `map_uniprot_to_kegg_via_rest` fonksiyonu kullanılmıştır.
### KO Terimleri
KEGG, enzimleri aratırken bu enzimlerin orthology(ortoloji) bilgilerini de barındırır. Bu bilgileri ise KEGG Orthology(KO) olarak kaydeder. 

KO terimleri aratılırsa, bu terimin:
- Yer aldığı yolaklar
- Yer aldığı reaksiyonlar
- Bu terime sahip genler
- BRITE hiyerarşisi
bulunabilir.

### BRITE Hiyerarşileri
BRITE hiyerarşileri:
- Genler ve Proteinler (KO terimleri de dahil)
- Moleküller ve Reaksiyonlar
- İlaçlar
- Hastalıklar
- Organzimalar ve Virüsler
için oluşturulmuş ve genelden özele giden bir listedir.

Bir KEGG ID'si birden fazla BRITE hiyerarşisine sahip olabilir. Modellemede sırasında en fazla KO terimleri ve Reaksiyonlar için olan BRITE hiyerarşileri elde edilmiştir.
>[!info]
>`KEGG_merge.py` içerisindeki `KEGGClient` sınıfı ve `kegg_search` fonksiyonu KO terimleri ve BRITE hiyerarşisi dışındaki bilgileri bulmaktadır.
> BRITE hiyerarşisi ve KO bilgisi `GenomeToPathway.py` dosyası ile bulunabilmektedir.

### GO Terimleri
Bunun dışında, ortak enzimleri ve işlevleri bulmak için Gene Ontology(GO) terimleri de UniProt sitesinden çekilmiştir. Bu işlem sırasında `bioservices` paketinden `QuickGO`  kullanılmıştır.

Bu terimler 3 ana başlık altında toplanır:
- *Biological Process(Biyolojik İşlev)*, gen ürünün hangi biyolojik işlemde yer aldığını gösterir. Örnekğin: Glikoliz, Oksidatif Fosforilasyon
- *Molecular Function(Moleküler İşlev)*, gen ürününü hangi biyokimyasal işlevlerde yer aldığını gösterir. Örneğin: ATP Hidrolizi, DNA Helikaz Aktivitesi
- *Cellular Component(Hücrese Konum)*, gen ürünün işlevini hücrenin hangi parçasında gerçekleştirdiğini gösterir. Örneğin: Ribozom, Hücre Zarı

Bütün bunların sonucunda ise bir gen ürününe ait bütün bilgiler sıralanmış olur ve isimlendirmedeki karışıklıklar önlenir.
>[!info]
>`KEGG_merge.py` içerisindeki `uniprot_search` fonksiyonu GO terimleri bulmasında kullanılmaktadır.

## Sonuçların Bir Araya Getirilmesi
`KEGG_merge.py` dosyasında yukarıdaki bilgilerin toplanması için kod yazılmıştır. Bunun dışında kod önceden çalışmış ve sonuç elde etmiş ise, uzun zaman harcayan adımları (annotation, KEGG ID eldesi, Pathway sorgulanması) atlayan kodlar eklenmiştir.

En son adım ise bütün bu bilgiler tek bir excel dosyası olarak kaydedilmiştir. Excel dosyasında sonuçlar şu şekilde gözükmektedir:

| gene | inference | product | EC_number | GO_Terms | kegg_id | pathways | reactions_by_pathway | compunds_by_pathway |
| ---- | --------- | ------- | --------- | -------- | ------- | -------- | -------------------- | ------------------- |
| Gen ADı | Uniprot ID'si | Enzim Adı | EC değeri | GO Terimleri | KEGG ID'si | Bulunduğu Yolaklar | Yolakların içerdiği reaksiyonlar | Reaksiyonların içeridği moleküller |

Bunun dışında S Matriksi oluşturan kod yazılmıştır ancak `kegg_merge.py` içine eklenmemiştir. 

[^1] : *Escherichia coli K-12* ve *Escherichia coli O157:H7*

[^2] : SRR35660962 ve ec1_S1_L001_R1_001.fastq.gz ile ec1_S1_L001_R2_001.fastq.gz