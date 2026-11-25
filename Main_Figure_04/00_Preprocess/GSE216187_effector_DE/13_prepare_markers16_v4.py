from pathlib import Path
import pandas as pd,numpy as np
r=Path(__file__).resolve().parents[1]
cr=['KLRG1','HLA-DRB5','IKZF1','GZMH','CX3CR1','ZNF683','FGFBP2','KLRD1','TRGC2']
rel=['CXCR4','DUSP2','NR4A2','TNFAIP3','FTH1','IFI6','CD27']
markers=cr+rel
old=pd.read_csv(r/'data/plot_ready/heatmap_row_zscores.tsv',sep='\t',index_col=0)
expression=pd.read_csv(r/'data/plot_ready/patient_log_expression.tsv',sep='\t',index_col=0)
assert set(markers)<=set(expression.index)
added=[g for g in markers if g not in old.index]
new_expression=expression.loc[added,old.columns]
new_z=new_expression.sub(new_expression.mean(axis=1),axis=0).div(new_expression.std(axis=1,ddof=1),axis=0).fillna(0)
z=pd.concat([old,new_z])
assert z.index.is_unique and np.allclose(z.loc[old.index],old)
with (r/'data/plot_ready/heatmap_row_zscores_markers16_v4.tsv').open('x',encoding='utf8',newline='') as f:z.to_csv(f,sep='\t',index_label='gene',float_format='%.15g')
de=pd.read_csv(r/'data/processed/CR_vs_REL_combined_effectors_DE.tsv',sep='\t').set_index('gene').loc[markers].copy()
de['requested_group']=['CR']*len(cr)+['REL']*len(rel)
de['observed_direction_matches_requested']=((de.log2FC>0)&(de.requested_group=='CR'))|((de.log2FC<0)&(de.requested_group=='REL'))
with (r/'data/processed/requested_gene_results_markers16_v4.tsv').open('x',encoding='utf8',newline='') as f:de.to_csv(f,sep='\t',index_label='gene',float_format='%.15g')
s=(r/'code/11_compact_heatmap_v3.R').read_text(encoding='utf8')
s=s.replace('heatmap_row_zscores.tsv','heatmap_row_zscores_markers16_v4.tsv')
s=s.replace("c('KLRG1','IKZF1','GZMH','CX3CR1','ZNF683','FGFBP2','KLRD1')",repr(cr).replace('[','c(').replace(']',')'))
s=s.replace("c('CXCR4','DUSP2','CD27')",repr(rel).replace('[','c(').replace(']',')'))
s=s.replace("padding=unit(.4,'mm')","padding=unit(.15,'mm')")
s=s.replace('Effector_CR_REL_patient_line_heatmap_compact_v3','Effector_CR_REL_patient_line_heatmap_markers16_compact_v4')
s=s.replace('heatmap_display_order_compact_v3.tsv','heatmap_display_order_markers16_compact_v4.tsv')
s=s.replace('COMPACT_HEATMAP_READY','MARKERS16_COMPACT_HEATMAP_READY')
with (r/'code/14_heatmap_markers16_compact_v4.R').open('x',encoding='utf8') as f:f.write(s)
direct=s[:s.index('dir.create(')]+"\ndraw(ht,heatmap_legend_side='bottom',annotation_legend_side='bottom',merge_legends=TRUE,padding=unit(c(3,4,3,4),'mm'))\n"
with (r/'code/15_direct_markers16_compact_v4.R').open('x',encoding='utf8') as f:f.write(direct)
print('Rows:',len(z),'Added rows:',added)
print(de[['log2FC','P','requested_group','observed_direction_matches_requested']].to_string())
