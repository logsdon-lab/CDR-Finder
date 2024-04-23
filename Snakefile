import pandas as pd
import numpy as np
from Bio import SeqIO

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

configfile: 'config.yaml'

TARGET_BED = config['target_bed']
METH_TSV = config['meth_tsv']
WINDOW = config['window_size']
SAMPLE = config['sample']
FASTA=config['fasta']
SLIDE=config['slide']
REPORT_THRESHOLD=config.get('report_threshold', 25000)

rule all:
	input:
		expand("results/{sample}_CDR.bed", sample=SAMPLE)

def getOverlap(alin_rec):
	a = [alin_rec['start'], alin_rec['end']]
	b =[alin_rec['bin_start'],alin_rec['bin_stop']]
	return max(0, min(a[1], b[1]) - max(a[0], b[0]))


rule subset_fasta:
	input:
		fasta = FASTA,
		bed = TARGET_BED
	output:
		rm_fasta = 'results/fasta/{sample}_subset.fasta'
	resources:
		mem=8,
		hrs=1,
	threads: 1
	shell:
		'''
		bedtools getfasta -fi {input.fasta} -bed {input.bed} > {output.rm_fasta}
		'''

rule subset_meth:
	input:
		methylation_tsv = METH_TSV,
		target_bed = TARGET_BED
	output:
		subset_bed = "results/{sample}_subset.bed"
	resources:
		mem=8,
		hrs=12,
	threads: 1
	shell:
		'''
		module load bedtools/2.29.2
		bedtools intersect -a {input.methylation_tsv} -b {input.target_bed} -wa > {output.subset_bed}
		'''
rule run_rm:
	input:
		fasta = rules.subset_fasta.output.rm_fasta
	output:
		rm_out = "results/rm/{sample}_subset.fasta.out"
	resources:
		mem=8,
		hrs=12,
	params:
		directory = "results/rm/",
	threads: 12
	shell:
		'''
		module load RepeatMasker/4.1.0
		RepeatMasker -species human -dir {params.directory} -qq -pa {threads} {input.fasta}
		'''

rule calc_windows:
	input:
		methylation_tsv = rules.subset_meth.output.subset_bed,
		target_bed = TARGET_BED
	output:
		binned_freq = "results/{sample}_binned_freq.bed"
	resources: 
		mem = 8,
		hrs = 1
	threads: 1
	params:
		window_size = WINDOW,
		slide = SLIDE
	run:
		bed_df = pd.read_csv(input.target_bed, header=None, sep="\t", names=['chrom','start','stop'])
		meth_df = pd.read_csv(input.methylation_tsv, header=None, sep="\t", names=['chrom','start','end','mod_base_name','score','strand','n1','n2','n3','coverage','freq','ncanon','nmod','nfilt','nnocall','naltmod'])
		keep_window = pd.DataFrame()

		for index in bed_df.index:
			start = bed_df.at[index,'start']
			stop = bed_df.at[index,'stop']
			chrom = bed_df.at[index,'chrom']
			chrom_meth_df = meth_df.loc[meth_df['chrom'] == chrom].copy()
			seq = np.arange(int(start), int(stop)+int(params.slide), int(params.slide))
			window_size=(int(int(params.window_size)/int(params.slide)))+1
			window_avg = pd.DataFrame(columns=['Chr','Bin','Freq'])
			
			for i in range(len(seq) - window_size + 1):
				window = seq[i: i + window_size]
				window_meth_pre = chrom_meth_df.copy()
				RM_start = window[0] - start
				RM_stop = window[-1] - start
				window_meth_pre['bin_start'] = window[0]
				window_meth_pre['bin_stop'] = window[-1]
				print(RM_start,RM_stop)
				if len(window_meth_pre) == 0:
					continue
				window_meth_pre['overlap'] = window_meth_pre.apply(getOverlap, axis=1)
				window_meth = window_meth_pre.loc[window_meth_pre['overlap'] != 0]
				if len(window_meth) == 0:
					continue
				avg_meth = window_meth['freq'].mean()
				window_avg =window_avg.append({'Chr':str(chrom), 'Bin':str(RM_start)+"-"+str(RM_stop), 'Freq':avg_meth}, ignore_index=True)

			window_avg = window_avg.round({'Freq':2})
			keep_window = keep_window.append(window_avg)

		keep_window.reset_index(inplace=True)

		with open(output.binned_freq, "w+") as outfile:
			for idx_2 in keep_window.index:
				pos1 = str(keep_window.at[idx_2,'Bin']).split("-")[0]
				pos2 = str(keep_window.at[idx_2,'Bin']).split("-")[1]
				freq = str(keep_window.at[idx_2,'Freq'])
				chrom = str(keep_window.at[idx_2,'Chr'])
				outfile.write(chrom+"\t"+pos1+"\t"+pos2+"\t"+str(freq)+"\t"+str(idx_2)+"\n")

rule format_RM:
	input:
		repeat_masker = rules.run_rm.output.rm_out,
	output:
		repeat_bed = "results/{sample}_rm.bed"
	resources:
		mem = 8,
		hrs = 1,
	threads: 1
	shell:
		'''
		awk -v OFS="\\t" '{{print $5, $6, $7, $10, $9}}' {input.repeat_masker} | grep "ALR" > {output.repeat_bed}
		'''
rule filter_RM:
	input:
		repeat_masker = rules.format_RM.output.repeat_bed
	output:
		rm_bed = "results/{sample}_rm_ALR.bed"
	resources:
		mem = 8,
		hrs = 1
	threads: 1
	run:
		RM_df = pd.read_csv(input.repeat_masker, header=None, sep="\t", names=['chr','start','stop','ALR','or'])
		RM_df['chr']=RM_df['chr'].str.split(pat=":").str[0]
		RM_df.to_csv(output.rm_bed, header=None, index=None, sep="\t")

rule intersect_RM:
	input:
		repeat_masker = rules.filter_RM.output.rm_bed,
		binned_freq = rules.calc_windows.output.binned_freq
	output:
		intersect_bed = "results/{sample}_intersect.bed",
		merged = "results/{sample}_rm_merged.bed"
	resources:
		mem = 8,
		hrs = 1
	threads: 1
	shell:
		'''
		module load bedtools/2.29.2
		bedtools merge -i {input.repeat_masker} -d 500 > {output.merged}
		bedtools intersect -a {input.binned_freq} -b {output.merged} -f 1 -wa  -u > {output.intersect_bed}
		'''
rule in_threshold:
	input:
		intersect_bed = rules.intersect_RM.output.intersect_bed,
		target_bed = TARGET_BED
	output:
		final_call = "results/{sample}_CDR.bed"
	resources:
		mem = 8,
		hrs = 1
	threads: 1
	run:
		intersect_bed = pd.read_csv(input.intersect_bed, header=None, sep="\t", names=['chrom','start','stop','freq','idx'])
		
		with open(input.target_bed) as infile, open(output.final_call, "w+") as outfile:
			for line in infile:
				window_avg = pd.DataFrame(columns=['Bin','Freq'])
				chrom, start, stop = line.split("\t")[0], line.split("\t")[1], line.split("\t")[2]
				intersect_bed_sub = intersect_bed.loc[intersect_bed['chrom'] == chrom]
				for index in intersect_bed_sub.index:
					chrom_int = intersect_bed_sub.at[index,'chrom']
					start_int = int(intersect_bed_sub.at[index,'start']) +int(start)
					stop_int = int(intersect_bed_sub.at[index,'stop']) +int(start)
					freq = intersect_bed_sub.at[index,'freq']
					idx = intersect_bed_sub.at[index,'idx']
					window_avg =window_avg.append({'Bin':str(start_int)+"-"+str(stop_int), 'Freq':float(freq), 'idx_2': int(idx)}, ignore_index=True)
				threshold = np.median(window_avg['Freq'].tolist())
				print(threshold)
				window_avg = window_avg.set_index('idx_2')
				window_avg = window_avg.loc[window_avg['Freq'] <= threshold]
				s = pd.Series(window_avg.index.tolist())
				groups = s.groupby(s.diff().ne(1).cumsum()).apply(lambda x: [x.iloc[0], x.iloc[-1]] if len(x) >= 2 else [x.iloc[0]]).tolist()
				groups = [sub for sub in groups if len(sub)>1]
				for y in groups:
					pos1 = str(window_avg.at[y[0],'Bin']).split("-")[0]
					pos2 = str(window_avg.at[y[1],'Bin']).split("-")[1]
					if int(pos2)-int(pos1) > int(REPORT_THRESHOLD):
						outfile.write(chrom+"\t"+pos1+"\t"+pos2+"\n")
