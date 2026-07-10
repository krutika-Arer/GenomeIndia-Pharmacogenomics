import pandas as pd
import re

INPUT = "output/GenomeIndia_ClinPGx_Clinical_RSMatch.tsv"
OUTPUT = "output/GenomeIndia_ClinPGx_Clinical_SelectedPhenotypes.tsv"

# Load file
print("Loading file...")

df = pd.read_csv(
    INPUT,
    sep="\t",
    dtype=str,
    low_memory=False
).fillna("")

#############################################################
# Phenotypes to retain
#############################################################

phenotypes = [

"Severe Cutaneous Adverse Reactions",
"Stevens-Johnson Syndrome",
"Toxic Epidermal Necrolysis",
"Ovarian Neoplasms",
"Acute lymphoblastic leukemia",
"Lymphoma",
"Osteosarcoma",
"Burkitt Lymphoma",
"Colorectal Neoplasms",
"Hypertension",
"Neoplasms",
"Acute coronary syndrome",
"Angina Pectoris",
"Drug Toxicity",
"Hematopoietic stem cell transplantation",
"Psoriasis",
"Rheumatoid arthritis",
"Toxic liver disease",
"Myocardial Infarction",
"Diabetes Mellitus, Type 2",
"Agranulocytosis",
"Tuberculosis",
"Non-Small Cell Lung Carcinoma",
"Thrombocytopenia",
"HIV infectious disease",
"Drug Reaction with Eosinophilia and Systemic Symptoms",
"Acquired Immunodeficiency Syndrome",
"Inflammatory Bowel Diseases",
"Pancreatitis",
"Prostatic Neoplasms",
"Colonic Neoplasms",
"Gastrointestinal toxicity",
"Chronic myelogenous leukemia, BCR-ABL1 positive",
"Diabetes Mellitus",
"Leukemia, Myeloid, Acute",
"Mesothelioma",
"Urinary Bladder Neoplasms",
"Breast Neoplasms",
"Gastrointestinal Stromal Tumors",
"Coronary Artery Disease",
"Multiple Myeloma",
"Major Adverse Cardiac Events (MACE)",
"Hypercholesterolemia",
"Renal Cell Carcinoma",
"Pancreatic Neoplasms",
"Rectal Neoplasms",
"Brain Neoplasms",
"Essential hypertension",
"Leukopenia",
"Angioedema",
"Metastatic neoplasm",
"Hemorrhage",
"Stomach Neoplasms",
"Hyperlipidemias",
"Anemia",
"Nausea",
"Vomiting",
"Nasopharyngeal Neoplasms",
"Hypertrophy, Left Ventricular",
"Kidney Neoplasms",
"Neuroendocrine Tumors",
"Cough",
"Hyperglycemia",
"Lymphoma, Large B-Cell, Diffuse",
"Cardiovascular Disease",
"Rhabdomyolysis",
"Residual Neoplasm",
"Small cell carcinoma",
"Dermatitis",
"Myelosuppression",
"Liver cancer",
"Hepatitis C virus infection",
"Hepatocellular Carcinoma",
"Drug-induced liver injury",
"Neurotoxicity Syndromes",
"Exanthema",
"Osteonecrosis",
"Lung Neoplasms",
"Hypertriglyceridemia",
"Osteoporosis",
"Heart Failure",
"Vasculitis",
"Heart valve replacement",
"Nephrolithiasis",
"Peripheral Nervous System Diseases",
"Nephrotoxicity",
"Overdose",
"Myeloproliferative Disorder",
"Diarrhea",
"Ototoxicity",
"Hematologic Neoplasms",
"Menopause",
"Neoplasm of esophagus",
"Discontinuation",
"Mucositis"

]

#############################################################
# Build regex
#############################################################

pattern = "|".join(
    re.escape(x)
    for x in sorted(phenotypes, key=len, reverse=True)
)

#############################################################
# Filter
#############################################################

filtered = df[
    df["Clinical_Phenotypes"].str.contains(
        pattern,
        case=False,
        regex=True,
        na=False
    )
].copy()

#############################################################
# Save
#############################################################

filtered.to_csv(
    OUTPUT,
    sep="\t",
    index=False
)

print("\nFinished")
print("Original rows :", len(df))
print("Filtered rows :", len(filtered))
print("Saved to      :", OUTPUT)

#############################################################
# Summary counts
#############################################################

print("\nTop Phenotypes:")

counts = {}

for p in phenotypes:
    counts[p] = filtered["Clinical_Phenotypes"].str.contains(
        re.escape(p),
        case=False,
        regex=True,
        na=False
    ).sum()

summary = (
    pd.DataFrame({
        "Phenotype": counts.keys(),
        "Count": counts.values()
    })
    .sort_values("Count", ascending=False)
)

summary.to_csv(
    "output/SelectedPhenotype_Counts.tsv",
    sep="\t",
    index=False
)

print(summary.head(20))
