#!/usr/bin/env python
import os
import sys
import pandas as pd
from gimmemotifs.motif import read_motifs
from gimmemotifs.cluster import cluster_motifs
from tempfile import NamedTemporaryFile

from Bio import Entrez
try:
    from urllib.error import HTTPError  # for Python 3
except ImportError:
    from urllib2 import HTTPError  # for Python 2


def get_species(tax_id, email, batch_size=500):
    """ Given a taxonomy id return the name of all children species

    Parameters
    ----------
    tax_id : str or int
        NCBI taxonomy id

    email : str
        e-mail address, required by Entrez API

    batch_size : int, optional
        Batch size to retrieve results from Entrez

    Yields
    ------

    str, species name
    """
    
    try:
        int(tax_id)
    except:
        raise ValueError(
            "tax_id {} does not look like a valid numerical id".format(tax_id)
            )
    
    Entrez.email = email

    # Initial query
    handle = Entrez.esearch(term="txid{}[Subtree]".format(tax_id), 
            db="taxonomy", usehistory="y")
    record = Entrez.read(handle)
    webenv = record["WebEnv"]
    query_key = record["QueryKey"]
    count = int(record["Count"])
    
    # Retrieve result
    for start in range(0, count, batch_size):
        end = min(count, start + batch_size)
        attempt = 1
        while attempt <= 3:
            try:
                fetch_handle = Entrez.efetch(db="taxonomy", rettype="summary", retmode="text",
                                         retstart=start, retmax=batch_size,
                                         webenv=webenv, query_key=query_key)
                break
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    sys.stderr.write("Received error from server %s\n" % err)
                    sys.stderr.write("Attempt %i of 3\n" % attempt)
                    attempt += 1
                    time.sleep(15)
                else:
                    raise
        
        lines = fetch_handle.read().splitlines()
        fetch_handle.close()
        for pos in range(0, len(lines), 2):
            name = lines[pos]
            clf = lines[pos + 1]
            if clf.find("species") != -1:
                yield " ".join(name.split(" ")[1:])


def cluster_pwms(pwm_dir, df_cis_bp, outfile):

    with open(outfile, "w") as f:
        motifs = []
        for family in df_cis_bp["Family_Name"].unique():
            sys.stderr.write("Family: {}\n".format(family))
            my_motifs = []
            for name in df_cis_bp[df_cis_bp["Family_Name"] == family]["Motif_ID"]:
                pwm = os.path.join(pwm_dir, name) + ".pwm"
                if os.path.exists(pwm):
                    my_motifs += read_motifs(open(pwm))
                    
            motifs += my_motifs 
            fam = family.replace("/","_").replace(" ", "_").replace(",", "_")
            
            # Write motifs to cluster to temporary file    
            tmp = NamedTemporaryFile()
            for motif in my_motifs:
                tmp.write("{}\n".format(motif.to_pfm()))
            tmp.flush()
            
    
            if len(my_motifs) > 1:
                threshold = 0.9999
                tree = cluster_motifs(tmp.name, "total", "wic", "mean", True, threshold, include_bg=True)
                clusters = tree.getResult()
    
                for cluster, members in clusters:
                    cluster.id = "{}_{}".format(fam, cluster.id)
                    f.write("# {}\t{}\n".format(cluster.id, [m.id for m in members]))
                    f.write("{}\n".format(cluster.to_pwm()))
            elif len(my_motifs) == 1:
                my_motifs[0].id = "{}_{}".format(fam, my_motifs[0].id)
                f.write("#{}\t{}\n".format(my_motifs[0].id,[my_motifs[0].id]))
                f.write("{}\n".format(my_motifs[0].to_pfm()))
    
def filter_cis_bp(df_cis_bp, tax_ids, email):
    if not isinstance(tax_ids, list):
        tax_ids = [tax_ids]
    
    species = []
    for tax_id in tax_ids:
        for name in get_species(tax_id, email):
            species.append(name.replace(" ", "_"))
    
    if len(species) == 0:
        raise ValueError(
                "no species find with tax_ids {}".format(",".join(tax_ids))
                )
    
    # return dataframe filtered by species
    return df_cis_bp[df_cis_bp["TF_Species"].isin(species)]
    
def parse_args(args):
    pass

if __name__ == "__main__":

    cis_bp = "/data/cis-bp/TF_Information_all_motifs.txt"
    tax_ids = 9605
    
    # Fungi 
    #tax_ids = [4751]

    # Plants
    #tax_ids = [33090]
    
    # Vertebrates
    #tax_ids = [7742]

    # Invertebrates (all metazoa except vertebrates)
    tax_ids = [
        10213,      # Mesozoa
        10226,      # Placozoa
        6040,       # Porifera
        7735,       # Cephalochordata 
        7712,       # Tunicata 
        10229,      # Chaetognatha 
        7586,       # Echinodermata 
        10219,      # Hemichordata 
        66780,      # Gnathostomulida
        6157,       # Platyhelminthes 
        33317,      # Protostomia
        1312402,    # Xenacoelomorpha
        6073,       # Cnidaria
        10197,      # Ctenophora
        ]

    email = "s.vanheeringen@science.ru.nl"
    pwm_dir = "/data/cis-bp/pwms"
    outfile = "invertebrate.clustered.pwm"

    sys.stderr.write("reading cis-bp table\n")
    df_cis_bp = pd.read_table(cis_bp)

    sys.stderr.write("downloading species and filtering\n")
    df_cis_bp = filter_cis_bp(df_cis_bp, tax_ids, email)
    
    df_filtered = df_cis_bp[df_cis_bp["TF_Status"] == "D"]
    sys.stderr.write("clustering motifs\n")
    cluster_pwms(pwm_dir, df_filtered, outfile)

    motifs = read_motifs(open(outfile))
    sys.stderr.write("got {} clustered motifs\n".format(len(motifs)))
