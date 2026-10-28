import gseapy as gp
import pandas as pd


kegg_metabolic = pd.read_csv("/io01dat03/PROJECT/zoujiping/brca/st/a_gsva/kegg_metabolic_terms.csv", index_col=0)
kegg_metab_dict = kegg_metabolic.groupby('gs_name')['gene_symbol'].apply(list).to_dict()

es = gp.gsva(data = '/io01dat03/PROJECT/zoujiping/brca/st/a_gsva/metacell_gsva_filtered_input.txt',
             gene_sets=kegg_metab_dict,
             outdir=None,
             kcdf="Gaussian",
             threads=24,
             max_size=1000,
             min_size=0,
             verbose=True)

es.res2d.pivot(index='Term', columns='Name', values='ES').to_csv('/io01dat03/PROJECT/zoujiping/brca/st/a_gsva/epi_metacell_gsva_filtered_kegg_res.csv')
