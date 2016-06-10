import pandas as pd

pwmfile = "motif_databases/gimme.vertebrate.v3.1.pwm"
df = pd.read_table("/data/cis-bp/TF_Information_all_motifs.txt")

df = df[df["TF_Species"] == "Mus_musculus"]

for line in open(pwmfile):
    if line.startswith("#"):
        name,motifs = line.strip().split("\t")
        name = name.strip("# ")
        motifs = [m.strip(" []'") for m in motifs.split(",")]
        f = df["Motif_ID"].isin(motifs)
        
        factors = set(df.loc[f, "TF_Name"])
        if len(factors) > 0:
            fam_name = list(set(df.loc[f]["Family_Name"]))[0]
            fam_name = fam_name.replace("/","_").replace(" ", "_").replace(",", "_")
            if not name.startswith(fam_name):
                name = fam_name + "_" + name
            print "{}\t{}".format(
                    name,
                    ",".join(factors)
                    )

