import requests, sys
import json

for server in ["http://rest.ensembl.org"]:#, "http://rest.ensemblgenomes.org"]:
    ext = "/info/species?"
  
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
   
    if not r.ok:
        r.raise_for_status()
        sys.exit()
            
             
    data = json.loads(r.text)
    for sp in data['species']:
        tax_id = sp[u'taxon_id']
        for name in [sp[u'name']] + sp[u'aliases']:
            print("{}\t{}".format(name, tax_id))

