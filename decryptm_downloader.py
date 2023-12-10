from urllib.request import urlretrieve
import requests
import pandas as pd
import os

experiment_summary = pd.ExcelFile('/net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/Experiment_summary.xlsx')
experiment_summary_phospho = pd.read_excel(experiment_summary, 'PTMs')

rawfiles = experiment_summary_phospho['Raw file'][
    (experiment_summary_phospho['Drug'].isin(['Afatinib', 'Dasatinib', 'Gefitinib'])) & 
    (experiment_summary_phospho['Treatment time'] == '30min') & 
    (experiment_summary_phospho['Cell line'] == 'A431') &
    (experiment_summary_phospho["Search"].str.contains('EGFR'))
    ].tolist()

os.mkdir('/net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/all_phospho')

for rawfile in rawfiles:
    print(f'getting {rawfile}')
    url = f"https://ftp.pride.ebi.ac.uk/pride/data/archive/2023/03/PXD037285/{rawfile}"
    filename = f"/net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/all_phospho/{rawfile}"

    try:
        resp = requests.get(url).content
        with open(filename, "wb") as f:
            f.write(resp)
        print(f"{rawfile} is saved")
    except Exception as e:
        print(e)
