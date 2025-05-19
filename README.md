# Testiniai duomenys sRNA Shiny aplikacijai

Šiame kataloge pateikti testiniai duomenys, skirti naudoti su [sRNA Shiny aplikacija](ip adresas).

## 📁 Failų sąrašas

| Failas         							  | Aprašymas                                       |
|-------------------------------------------------------------------------|-------------------------------------------------|
| `test.bam`      							  | Pavyzdinis BAM failas (grupė 1)                 |
| `test2.bam`     							  | Pavyzdinis BAM failas (grupė 2)                 |
| `NC_009004.1.ptt`    							  | Genomo anotacijų failas (PTT formatas)          |
| `NC_009004.1.gff3`    						  | Genomo anotacijų failas (GFF3 formatas)         |
| `NC_009004.1.fasta`   						  | Referencinė genomo seka (FASTA formatas)        |
| `Lactococcus cremoris strain 4B0 chromosome, complete genome.fasta`     | Papildoma seka BLAST analizei                   |
| `Lactococcus cremoris strain UL033 chromosome, complete genome.fasta`   | Papildoma seka BLAST analizei                   |
| `Lactococcus lactis subsp. cremoris MG1363, complete genome.fasta`      | Papildoma seka BLAST analizei                   |

## 🚀 Naudojimo instrukcija

1. **Atidarykite aplikaciją** naršyklėje per šią nuorodą: [ip adresas]

2. **Analizės žingsniai**:
   - **APERO analizė (First step)** – jei norite pradėti nuo BAM failų analizės.
   - **APERO rezultatų analizė (Second step)** – jei turite jau paruoštus APERO rezultatus.

3. **Pasirinkite pirmą analizės žingsnį**

4. **Įkelkite atitinkamus failus**:
   - **BAM failai**: `test.bam`
   - **PTT failas**: `NC_009004.1.ptt`
   - **GFF3 failas**: `NC_009004.1.gff3`

5. **Paspauskite mygtuką „Run APERO“** ir palaukite, kol analizė bus baigta. Atsisiųskite sugeneruotą `.csv` failą.

6. **Pakartokite analizę su `test2.bam` failu. Atsisiųskite sugeneruotą `.csv` failą.**

7. **Pasirinkite antrą analizės žingsnį.**

8. **Įkelkite pirmą sugeneruotą `.csv` failą į `Upload APERO Group 1 CSVs` ir suteikite jam pavadinimą, tada įkelkite antrą sugeneruotą `.csv` failą į `Upload APERO Group 2 CSVs` ir taip pat atitinkamai jį pavadinkite. Taip pat įkelkite**:
   - **PTT failas**: `NC_009004.1.ptt`
   - **FASTA failas**: `NC_009004.1.fasta`
   - **Papildomi FASTA failai BLAST analizei**: `Lactococcus cremoris strain 4B0 chromosome, complete genome.fasta`, `Lactococcus cremoris strain UL033 chromosome, complete genome.fasta`, `Lactococcus lactis subsp. cremoris MG1363, complete genome.fasta`

4. **Pasirinkite bakterijos rūšį**: *Lactococcus lactis* arba *L. casei*- testiniai duomenys yra *Lactococcus lactis*.

5. **Paspauskite mygtuką „Process Files“** ir palaukite, kol analizė bus baigta.

6. **Peržiūrėkite rezultatus** skirtingose kortelėse:
   - **Tables** – analizės lentelės.
   - **Plots** – grafiniai atvaizdai.
   - **Genome browser** – interaktyvus genomo naršymas.
   - **Secondary structures** – antrinės struktūros.

7. **Atsisiųskite rezultatus** naudodami atitinkamus mygtukus:
   - **Excel failas** su visomis lentelėmis.
   - **ZIP archyvas** su visomis antrinių struktūrų nuotraukomis.

## 📝 Pastabos

- Visi failai yra paruošti taip, kad juos būtų galima tiesiogiai įkelti į aplikaciją be papildomo apdorojimo.
- Aplikacija veikia tik tuo metu, kai ji yra paleista ... kompiuteryje ir prieinama per pateiktą nuorodą.

## 📄 Licencija

Šie testiniai duomenys skirti tik demonstraciniams ir edukaciniams tikslams.
