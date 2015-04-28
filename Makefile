export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files
LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'

PROJECT_HOME=~/iamp
TRIM_BEFORE_BASE=1

ER=C57C3ACXX_CV_C1_14s006573-1-1_Vesely_lane214s006573_sequence \
   C57C3ACXX_CV_C2_14s006574-1-1_Vesely_lane214s006574_sequence \
   C57C3ACXX_CV_C3_14s006575-1-1_Vesely_lane314s006575_sequence \
   C57C3ACXX_CV_C4_14s006576-1-1_Vesely_lane414s006576_sequence \
   C57C3ACXX_CV_C5_14s006577-1-1_Vesely_lane414s006577_sequence \
   C57C3ACXX_CV_C6_14s006578-1-1_Vesely_lane414s006578_sequence \
   C57C3ACXX_CV_C7_14s006579-1-1_Vesely_lane314s006579_sequence \
   C57C3ACXX_CV_C8_14s006580-1-1_Vesely_lane414s006580_sequence \
   C57C3ACXX_CV_C9_14s006581-1-1_Vesely_lane314s006581_sequence \
   C57C3ACXX_CV_C10_14s006582-1-1_Vesely_lane314s006582_sequence \
   C57C3ACXX_CV_C11_14s006583-1-1_Vesely_lane414s006583_sequence \
   C57C3ACXX_CV_C12_14s006584-1-1_Vesely_lane414s006584_sequence

IAMP=C57C3ACXX_CV_A1_14s006561-1-1_Vesely_lane114s006561_sequence \
     C57C3ACXX_CV_A2_14s006562-1-1_Vesely_lane114s006562_sequence \
     C57C3ACXX_CV_A3_14s006563-1-1_Vesely_lane114s006563_sequence \
     C57C3ACXX_CV_A4_14s006564-1-1_Vesely_lane114s006564_sequence \
     C57C3ACXX_CV_A5_14s006565-1-1_Vesely_lane114s006565_sequence \
     C57C3ACXX_CV_A7_14s006566-1-1_Vesely_lane114s006566_sequence \
     C57C3ACXX_CV_A8_14s006567-1-1_Vesely_lane114s006567_sequence \
     C57C3ACXX_CV_A9_14s006568-1-1_Vesely_lane214s006568_sequence \
     C57C3ACXX_CV_A10_14s006569-1-1_Vesely_lane214s006569_sequence \
     C57C3ACXX_CV_A11_14s006570-1-1_Vesely_lane214s006570_sequence \
     C57C3ACXX_CV_A12_14s006571-1-1_Vesely_lane214s006571_sequence \
     C57C3ACXX_CV_A13_14s006572-1-1_Vesely_lane214s006572_sequence
     
PC=C57C3ACXX_CV_D1_14s006585-1-1_Vesely_lane414s006585_sequence \
   C57C3ACXX_CV_D2_14s006586-1-1_Vesely_lane514s006586_sequence \
   C57C3ACXX_CV_D3_14s006587-1-1_Vesely_lane514s006587_sequence \
   C57C3ACXX_CV_D4_14s006588-1-1_Vesely_lane514s006588_sequence \
   C57C3ACXX_CV_D5_14s006589-1-1_Vesely_lane514s006589_sequence \
   C57C3ACXX_CV_D6_14s006590-1-1_Vesely_lane514s006590_sequence \
   C57C3ACXX_CV_D7_14s006591-1-1_Vesely_lane514s006591_sequence \
   C57C3ACXX_CV_D8_14s006592-1-1_Vesely_lane514s006592_sequence
   
S1=C57C3ACXX_CV_S1_14s006593-1-1_Vesely_lane314s006593_sequence
S2=C57C3ACXX_CV_S2_14s006594-1-1_Vesely_lane314s006594_sequence
S3=C57C3ACXX_CV_S3_14s006595-1-1_Vesely_lane314s006595_sequence

CD10=$(S1) $(S2) $(S3)
  
SAMPLES=$(ER) $(IAMP) $(PC) $(CD10)

SPACE :=
SPACE +=
COMMA := ,

#all: gsnap htseq qc blast fastqc
all: gsnap htseq qc fastqc diffexp gsea

include ~/generic/scripts/rna-seq/default.mk

#---------------
# DIFFERENTIALLY EXPRESSED GENES
#---------------

.PHONY: diffexp
diffexp: deseq/iAMP21.diff-exp-genes.deseq2.merged.pairwise.tsv

deseq/iAMP21.diff-exp-genes.deseq2.merged.pairwise.tsv: deseq/iAMP-vs-PC.tsv deseq/iAMP-vs-ER.tsv deseq/ER-vs-PC.tsv \
         										  deseq/ER-vs-immature.tsv deseq/ER-vs-preB.tsv deseq/ER-vs-mature.tsv \
                                                  deseq/PC-vs-immature.tsv deseq/PC-vs-preB.tsv deseq/PC-vs-mature.tsv \
                                                  deseq/iAMP-vs-immature.tsv deseq/iAMP-vs-preB.tsv deseq/iAMP-vs-mature.tsv 
	Rscript ~/iamp/scripts/merge-pairwise.R
	mv $@.part $@
	
deseq/iAMP-vs-PC.tsv: ~/generic/scripts/rna-seq/diff-exp.R $(foreach S, $(SAMPLES), htseq/$S.count)
	mkdir -p deseq
	Rscript ~/generic/scripts/rna-seq/diff-exp.R \
		--experiment $(subst $(SPACE),$(COMMA),$(foreach S,$(IAMP),$S.count)) \
		--control $(subst $(SPACE),$(COMMA),$(foreach S,$(PC),$S.count)) \
		--name-subst-pattern ".*CV_(.\\d+)_.*" \
		--output-tsv $@.part
	mv $@.part $@

deseq/iAMP-vs-ER.tsv: ~/generic/scripts/rna-seq/diff-exp.R $(foreach S, $(SAMPLES), htseq/$S.count)
	mkdir -p deseq
	Rscript ~/generic/scripts/rna-seq/diff-exp.R \
		--experiment $(subst $(SPACE),$(COMMA),$(foreach S,$(IAMP),$S.count)) \
		--control $(subst $(SPACE),$(COMMA),$(foreach S,$(ER),$S.count)) \
		--name-subst-pattern ".*CV_(.\\d+)_.*" \
		--output-tsv $@.part
	mv $@.part $@

deseq/ER-vs-PC.tsv: ~/generic/scripts/rna-seq/diff-exp.R $(foreach S, $(SAMPLES), htseq/$S.count)
	mkdir -p deseq
	Rscript ~/generic/scripts/rna-seq/diff-exp.R \
		--experiment $(subst $(SPACE),$(COMMA),$(foreach S,$(ER),$S.count)) \
		--control $(subst $(SPACE),$(COMMA),$(foreach S,$(PC),$S.count)) \
		--name-subst-pattern ".*CV_(.\\d+)_.*" \
		--output-tsv $@.part
	mv $@.part $@

#---------------------
# ER vs. normal B-cell
#---------------------
deseq/ER-vs-immature.tsv: ~/generic/scripts/rna-seq/diff-exp.R $(foreach S, $(SAMPLES), htseq/$S.count)
	mkdir -p deseq
	Rscript ~/generic/scripts/rna-seq/diff-exp.R \
		--experiment $(subst $(SPACE),$(COMMA),$(foreach S,$(ER),$S.count)) \
		--control $(S1).count \
		--name-subst-pattern ".*CV_(.\\d+)_.*" \
		--output-tsv $@.part
	mv $@.part $@

deseq/ER-vs-preB.tsv: ~/generic/scripts/rna-seq/diff-exp.R $(foreach S, $(SAMPLES), htseq/$S.count)
	mkdir -p deseq
	Rscript ~/generic/scripts/rna-seq/diff-exp.R \
		--experiment $(subst $(SPACE),$(COMMA),$(foreach S,$(ER),$S.count)) \
		--control $(S2).count \
		--name-subst-pattern ".*CV_(.\\d+)_.*" \
		--output-tsv $@.part
	mv $@.part $@

deseq/ER-vs-mature.tsv: ~/generic/scripts/rna-seq/diff-exp.R $(foreach S, $(SAMPLES), htseq/$S.count)
	mkdir -p deseq
	Rscript ~/generic/scripts/rna-seq/diff-exp.R \
		--experiment $(subst $(SPACE),$(COMMA),$(foreach S,$(ER),$S.count)) \
		--control $(S3).count \
		--name-subst-pattern ".*CV_(.\\d+)_.*" \
		--output-tsv $@.part
	mv $@.part $@

#---------------------
# PC vs. normal B-cell
#---------------------
deseq/PC-vs-immature.tsv: ~/generic/scripts/rna-seq/diff-exp.R $(foreach S, $(SAMPLES), htseq/$S.count)
	mkdir -p deseq
	Rscript ~/generic/scripts/rna-seq/diff-exp.R \
		--experiment $(subst $(SPACE),$(COMMA),$(foreach S,$(PC),$S.count)) \
		--control $(S1).count \
		--name-subst-pattern ".*CV_(.\\d+)_.*" \
		--output-tsv $@.part
	mv $@.part $@

deseq/PC-vs-preB.tsv: ~/generic/scripts/rna-seq/diff-exp.R $(foreach S, $(SAMPLES), htseq/$S.count)
	mkdir -p deseq
	Rscript ~/generic/scripts/rna-seq/diff-exp.R \
		--experiment $(subst $(SPACE),$(COMMA),$(foreach S,$(PC),$S.count)) \
		--control $(S2).count \
		--name-subst-pattern ".*CV_(.\\d+)_.*" \
		--output-tsv $@.part
	mv $@.part $@

deseq/PC-vs-mature.tsv: ~/generic/scripts/rna-seq/diff-exp.R $(foreach S, $(SAMPLES), htseq/$S.count)
	mkdir -p deseq
	Rscript ~/generic/scripts/rna-seq/diff-exp.R \
		--experiment $(subst $(SPACE),$(COMMA),$(foreach S,$(PC),$S.count)) \
		--control $(S3).count \
		--name-subst-pattern ".*CV_(.\\d+)_.*" \
		--output-tsv $@.part
	mv $@.part $@

#---------------------
# iAMP21 vs. normal B-cell
#---------------------
deseq/iAMP-vs-immature.tsv: ~/generic/scripts/rna-seq/diff-exp.R $(foreach S, $(SAMPLES), htseq/$S.count)
	mkdir -p deseq
	Rscript ~/generic/scripts/rna-seq/diff-exp.R \
		--experiment $(subst $(SPACE),$(COMMA),$(foreach S,$(IAMP),$S.count)) \
		--control $(S1).count \
		--name-subst-pattern ".*CV_(.\\d+)_.*" \
		--output-tsv $@.part
	mv $@.part $@

deseq/iAMP-vs-preB.tsv: ~/generic/scripts/rna-seq/diff-exp.R $(foreach S, $(SAMPLES), htseq/$S.count)
	mkdir -p deseq
	Rscript ~/generic/scripts/rna-seq/diff-exp.R \
		--experiment $(subst $(SPACE),$(COMMA),$(foreach S,$(IAMP),$S.count)) \
		--control $(S2).count \
		--name-subst-pattern ".*CV_(.\\d+)_.*" \
		--output-tsv $@.part
	mv $@.part $@

deseq/iAMP-vs-mature.tsv: ~/generic/scripts/rna-seq/diff-exp.R $(foreach S, $(SAMPLES), htseq/$S.count)
	mkdir -p deseq
	Rscript ~/generic/scripts/rna-seq/diff-exp.R \
		--experiment $(subst $(SPACE),$(COMMA),$(foreach S,$(IAMP),$S.count)) \
		--control $(S3).count \
		--name-subst-pattern ".*CV_(.\\d+)_.*" \
		--output-tsv $@.part
	mv $@.part $@

#---------------
# GSEA
#---------------

# NOTE: GSEA does not understand dashes in its filenames... so we use underscores instead
.PHONY: gsea
gsea: gsea-heatmap.chr.pdf

gsea-heatmap.chr.pdf: gsea/iAMP_vs_PC.gsea gsea/iAMP_vs_ER.gsea gsea/ER_vs_PC.gsea gsea/ER_vs_immature.gsea gsea/ER_vs_preB.gsea gsea/ER_vs_mature.gsea gsea/PC_vs_immature.gsea gsea/PC_vs_preB.gsea gsea/PC_vs_mature.gsea gsea/iAMP_vs_immature.gsea gsea/iAMP_vs_preB.gsea gsea/iAMP_vs_mature.gsea ~/iamp/scripts/gsea-heatmap.R
	Rscript ~/iamp/scripts/gsea-heatmap.R

gsea/%.gsea: gsea/%.rnk 
	mkdir -p gsea
	rm -rf $@*
	java -cp /home/STANNANET/christian.frech/tools/gsea-2.0.13/gsea2-2.0.13.jar -Xmx3048m xtools.gsea.GseaPreranked \
		-rpt_label $(notdir $@) \
		-rnk $< \
		-gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/c1.all.v4.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets/c2.all.v4.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets/c3.all.v4.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets/c4.all.v4.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets/c5.all.v4.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets/c6.all.v4.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets/c7.all.v4.0.symbols.gmt \
		-out gsea \
		-collapse false -mode Max_probe -norm meandiv -nperm 100 -scoring_scheme weighted -include_only_symbols true -make_sets true \
		-rnd_seed 149 \
		-plot_top_x 300 \
		-set_max 500 \
		-set_min 15 \
		-zip_report false \
		-gui false 
	rename 's/\.GseaPreranked\.\d+$$//' $@.*
	
.SECONDEXPANSION:
gsea/%.rnk: deseq/$$(subst _,-,%).tsv ~/generic/scripts/enrichment/deseq2gsea.R
	mkdir -p gsea 
	Rscript ~/generic/scripts/enrichment/deseq2gsea.R \
		--deseq2-input-file $< \
		--rnk-output-file $@.part \
		--min-pvalue 0.5
	mv $@.part $@
	
 