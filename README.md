import requests
import pandas as pd

# Query NCI Proteomic Data Commons API for CPTAC-3 Prostate (PDC000222)
url = "https://pdc.cancer.gov/graphql"

query = """
{
  uiGeneProteinExpression(
    pdc_study_id: "PDC000222"
    genes: ["STEAP1", "STEAP2", "KLK2", "FOLH1"]
  ) {
    gene_name
    case_submitter_id
    sample_submitter_id
    log2_ratio
  }
}
"""

response = requests.post(url, json={'query': query})
res_json = response.json()

if 'data' in res_json and res_json['data']['uiGeneProteinExpression']:
    df = pd.DataFrame(res_json['data']['uiGeneProteinExpression'])
    print("--- CPTAC PRAD (PDC000222) Proteomic Levels ---")
    print(df.head(10))
else:
    print("PDC Data Structure Response:", res_json)
