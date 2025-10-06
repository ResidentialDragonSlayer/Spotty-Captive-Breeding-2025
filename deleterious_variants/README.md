## 🔹 VCF Annotation Filtering to Detect Deleterious SNPs and Autosomal SNP Extraction

This workflow is designed to process annotated Variant Call Format (VCF) files to identify potentially deleterious single nucleotide polymorphisms (SNPs) and isolate those located on autosomal chromosomes.
It consists of two main steps: filtering for high-impact variants and extracting autosomal regions

---

### 🧬 Step 1: Filter VEP Annotations

Annotated VCF file is filtered using `filter_vep` to retain only variants with **HIGH** predicted impact.
HIGH impact on gene function—such as stop-gained, frameshift, or splice-site mutations—which are more likely to be deleterious

```bash
filter_vep -i "$vcf_annotated" \
  --filter "IMPACT is HIGH" \
  --format vcf \
  -o "$vcf_high"

```

### 🧬 Step 2: Extract Autosomal Chromosomes

```
bcftools filter -Oz \
  -r NC_080782.1,NC_080783.1,NC_080784.1,NC_080785.1,NC_080786.1,NC_080787.1,NC_080788.1,NC_080789.1,NC_080790.1,NC_080791.1,NC_080792.1,NC_080793.1,NC_080794.1,NC_080795.1,NC_080796.1,NC_080797.1,NC_080798.1,NC_080799.1 \
  "$vcf_high_gz" \
  -o "$vcf_autosomes"

```
