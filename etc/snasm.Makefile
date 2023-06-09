
SHELL := bash

CP := cp -v -f
MV := mv -v -f

MINIMAP2 := minimap2 -x asm20 -t 1 -c --cs

MAKEFLAGS += --no-builtin-rules
MAKEFLAGS += --no-builtin-variables

.SUFFIXES:
.DELETE_ON_ERROR:
.SECONDARY:
.ONESHELL:
.DEFAULT: all
.PHONY: all clean ids

IDS := $(shell cat $(IDFILE))
VCF_GZ_FILES := $(addsuffix /vcf.gz,$(IDS))
CNS_FILES := $(addsuffix /cns,$(IDS))
ALIGNED_BED_FILES := $(addsuffix /aligned.bed,$(IDS))

all : $(OUT).core.aln $(OUT).cov
	csvtk pretty -H -t $(OUT).cov

$(OUT).cov : $(REF).dict $(ALIGNED_BED_FILES)
	parallel -a $(IDFILE) -j $(CPUS) \
	"bedtools genomecov -g $< -i {}/aligned.bed -max 1 \
	 | tail -n 1| cut -f1,5 \
         | sed 's/^genome/{}/' \
        " \
	> $@

$(OUT).core.aln : $(OUT).aln
	snp-sites $< > $@

#	goalign compress \
#	--threads $(CPUS) < $< > $@


$(OUT).nwk : $(OUT).treefile
	$(CP) $< $@
	gotree draw text < $<

$(OUT).treefile : $(OUT).aln
	iqtree2 \
	-T $(CPUS) \
	-s $< \
	--prefix $(OUT) \
	-m GTR \
	--ufboot 1000 \
	-redo

# add golign clean sites here
$(OUT).aln : $(OUT).afa
	$(CP) $< $@

$(OUT).afa : $(CNS_FILES)
	cat $^ > $@

$(OUT).vcf.fofn : $(IDFILE) $(VCF_GZ_FILES)
	parallel -a $< -j 1 \
	"echo {}/vcf.gz" > $@

$(OUT).vcf : $(OUT).vcf.fofn
	bcftools merge \
	--threads $(CPUS) \
	--file-list $< \
	--output $@

#$(OUT).vcf : $(VCF_GZ_FILES)
#	bcftools merge \
#	--threads $(CPUS) \
#	--missing-to-ref \
#	$^ > $@

%/cns : %/cns1
	( echo ">$*" \
	; grep -v '>' $< ) \
	| seqkit seq \
	> $@

# don't include insertions

%/cns1 : $(REF) %/vcf.gz %/missing.bed
	bcftools consensus \
	-f $< \
	-s $* \
	--exclude '%ILEN>0' \
	--mask $(word 3,$^) \
	--mask-with '-' \
	--mark-snv lc \
	--mark-del '-' \
	--output $@ \
	$(word 2,$^)

%/vcf.gz : %/vcf
	bcftools convert -Oz $< > $@
	bcftools index $@

$(REF).dict : $(REF)
	samtools faidx $<
	cut -f1,2 $(<).fai | sort > $@

%/aligned.bed : %/paf
	cut -f 6,8,9 $< \
	| bedtools sort -i stdin \
	> $@

##	| csvtk mutate2 -t -H -e '$$2 + 1' \

%/missing.bed : $(REF).dict %/aligned.bed
	bedtools complement \
	-g $(word 1,$^) \
	-i $(word 2,$^) \
	> $@

%/vcf : $(REF) %/paf
	-paftools.js call \
	-L $(ALEN) -l $(ALEN) -s $* \
	-f $^ > $@

%/paf : $(REF) %/fna
	$(MINIMAP2) $+ \
	| sort -k6,6 -k8,8n \
	> $@

%/fna : %/input
	any2fasta -n -u $< > $@

clean :
	#echo $(RM) $(addsuffix /cms1,$(IDS))
	$(RM) */*.bed */cns* */vcf*

ids : $(IDFILE)
	@echo $(VCF_GZ_FILES)
	
#echo $(ALIGNED_BED_FILES)
#echo $(IDS)
#echo $(CNS_FILES)
