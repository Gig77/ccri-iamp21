use warnings FATAL => qw( all );
use strict;

use lib "/mnt/projects/generic/scripts";
use List::Util qw(min max);
use Vcf;
use Tabix;
use Data::Dumper;
use Getopt::Long;
use Carp;

my ($vcf_out, $header, $sample_identifier, $vcf_in);
my ($rmsk_file, $simplerepeat_file, $blacklist_file, $segdup_file, $g1k_accessible_file, $ucsc_retro_file, $remission_variants_file, $evs_file, $exac_file, $clinvar_file, $cosmic_mutation_file);
GetOptions
(
	"sample=s" => \$sample_identifier, # e.g. 314_rem_dia
	"vcf-in=s" => \$vcf_in,  # VCF input file
	"vcf-out=s" => \$vcf_out,  # filtered VCF output file
	"header" => \$header,  # if set, write header line to output
	"rmsk-file=s" => \$rmsk_file, # TABIX indexed UCSC table rmsk
	"simpleRepeat-file=s" => \$simplerepeat_file, # TABIX indexed UCSC table rmsk
	"blacklist-file=s" => \$blacklist_file, # TABIX indexed UCSC table wgEncodeDacMapabilityConsensusExcludable
	"segdup-file=s" => \$segdup_file, # TABIX indexed UCSC table genomicSuperDups
	"g1k-accessible=s" => \$g1k_accessible_file, # TABIX indexed UCSC table tgpPhase1AccessibilityPilotCriteria
	"ucscRetro=s" => \$ucsc_retro_file, # TABIX indexed UCSC table ucscRetroAli5
	"remission-variants-file=s" => \$remission_variants_file, # TABIX indexed file with variants found in remission samples (GATK)
	"evs-file=s" => \$evs_file, # TABIX indexed file with wariants from Exome Variant Server (http://evs.gs.washington.edu/EVS/)
	"exac-file=s" => \$exac_file, # TABIX indexed VCF file from Exome Aggregation Consortium
	"clinvar-file=s" => \$clinvar_file, # TABIX indexed VCF file for ClinVar variants
	"cosmic-mutation-file=s" => \$cosmic_mutation_file
);

# TABLE: filtered-variants
if ($header)
{
	print "sample\t";
	print "var_type\t";
	print "status\t";
	print "rejected_because\t";
	print "chr\t";
	print "pos\t";
	print "ref\t";
	print "alt\t";
	print "gene\t";
	print "add_genes\t";
	print "dbSNP\t";
	print "AF_1000G\t";
	print "AF_EVS\t";
	print "AF_ExAC\t";	
	print "ClinVar\t";
	print "impact\t";
	print "effect\t";
	print "non_silent\t";
	print "deleterious\t";
	print "exons\t";
	print "dp_tot\t";
	print "dp_ref\t";
	print "dp_var\t";
	print "freq\t";
	print "qual\t";
	print "aa_change\t";
	print "snpeff_effect\t";
	print "Polyphen2\t";
	print "SIFT\t";
	print "GERP++\t";
	print "SiPhy\t";
	print "InterPro\t";
	print "cosmic_hits_nt\t";
	print "cosmic_hits_aa\t";
	print "repeat\t";
	print "segdup\t";
	print "blacklist\t";
	print "g1k-accessible\t";	
	print "retro\t";
	print "rem_samples\n";	
}

my $debug = 1;

croak "ERROR: --sample not specified" if (!$sample_identifier);
croak "ERROR: --vcf-in not specified" if (!$vcf_in);
croak "ERROR: --rmsk-file not specified" if (!$rmsk_file);
croak "ERROR: --simpleRepeat-file not specified" if (!$simplerepeat_file);
croak "ERROR: --blacklist-file not specified" if (!$blacklist_file);
croak "ERROR: --segdup-file not specified" if (!$segdup_file);
croak "ERROR: --g1k-accessible not specified" if (!$g1k_accessible_file);
croak "ERROR: --ucscRetro not specified" if (!$ucsc_retro_file);
croak "ERROR: --evs-file not specified" if (!$evs_file);
croak "ERROR: --exac-file not specified" if (!$exac_file);
croak "ERROR: --clinvar-file not specified" if (!$clinvar_file);
croak "ERROR: --cosmic-mutation-file not specified" if (!$cosmic_mutation_file);

my %clnsig = (
	 0 => "Uncertain significance",
	 1 => "not provided",
	 2 => "Benign",
	 3 => "Likely benign",
	 4 => "Likely pathogenic",
	 5 => "Pathogenic",
	 6 => "drug response",
	 7 => "histocompatibility",
	 255 => "other"
);

my $vcf_sample_id = "Sample1";

#-------------------------
# read in auxiliary files 
#-------------------------

# read cosmic mutations
my (%cosmic, %cosmic_leuk);
my $entries_read = 0;
open(C, "$cosmic_mutation_file") or croak "ERROR: Could not open file $cosmic_mutation_file";
<C>; # skip header
while(<C>)
{
	chomp;
	
	my ($gene_name, $accession_number, $gene_cds_length, $hgnc_id, $sample_name, $id_sample, $id_tumour, $primary_site, $site_subtype, $primary_histology,	
		$histology_subtype, $genome_wide_screen, $mutation_id, $mutation_cds, $mutation_aa, $mutation_description, $mutation_zygosity,	
		$GRCh, $mutation_GRCh37_genome_position, $mutation_strand, $snp, $fathmm_prediction, $fathmm_score, $mutation_somatic_status, $pubmed_id, $id_study, $sample_source,
		$tumour_origin, $age, $comments) = split /\t/;
	
	next if ($mutation_somatic_status ne "Confirmed somatic variant");
	$gene_name =~ s/_ENST.*//;
	
	$cosmic{$mutation_GRCh37_genome_position} = defined $cosmic{$mutation_GRCh37_genome_position} ? $cosmic{$mutation_GRCh37_genome_position} + 1 : 1;
	$cosmic_leuk{$mutation_GRCh37_genome_position} = defined $cosmic_leuk{$mutation_GRCh37_genome_position} ? $cosmic_leuk{$mutation_GRCh37_genome_position} + 1 : 1
		if ($histology_subtype =~ /leukaemia/);
		
	if ($mutation_aa =~ /p\.(.)(\d+)(.+)/)
	{
		my ($prev_aa, $aa_pos, $after_aa) = ($1, $2, $3);
		$cosmic{"$gene_name:$prev_aa:$aa_pos"} = defined $cosmic{"$gene_name:$prev_aa:$aa_pos"} ? $cosmic{"$gene_name:$prev_aa:$aa_pos"} + 1 : 1;
		$cosmic_leuk{"$gene_name:$prev_aa:$aa_pos"} = defined $cosmic_leuk{"$gene_name:$prev_aa:$aa_pos"} ? $cosmic_leuk{"$gene_name:$prev_aa:$aa_pos"} + 1 : 1
			if ($histology_subtype =~ /leukaemi/);
	}
	$entries_read ++;
}
close(C);
print STDERR "$entries_read mutations read from file $cosmic_mutation_file\n";

my $rmsk = Tabix->new(-data => $rmsk_file);
my $simpleRepeat = Tabix->new(-data => $simplerepeat_file);
my $blacklistdb = Tabix->new(-data => $blacklist_file);
my $segdup = Tabix->new(-data => $segdup_file);
my $g1kAcc = Tabix->new(-data => $g1k_accessible_file);
my $ucscRetro = Tabix->new(-data => $ucsc_retro_file);
my $remission = Tabix->new(-data => $remission_variants_file) if ($remission_variants_file);
my $evs = Tabix->new(-data => $evs_file);
my $exac = Tabix->new(-data => $exac_file);
my $clinvar = Tabix->new(-data => $clinvar_file);

#-----------
# parse VCF 
#-----------

$| = 1; # turn on autoflush

my %variant_stati = 
(
	0 => 'wildtype',
	1 => 'germline',
	2 => 'somatic',
	3 => 'LOH',
	4 => 'post-transcriptional modification',
	5 => 'unknown'
);

print STDERR "Processing file $vcf_in...\n";

my $vcf = Vcf->new(file => "$vcf_in");
$vcf->parse_header();

my (@samples) = $vcf->get_samples();

# sanity checks
die "ERROR: Sample name $vcf_sample_id not found!\n" if ($vcf_sample_id ne $samples[0]);

if ($vcf_out) 
{
	my $cmd = "grep -P '^#' $vcf_in > $vcf_out";
	system($cmd) == 0 or die "ERROR: grep vcf header failed: $cmd\n";
	open(VCFOUT,">>$vcf_out") or die "ERROR: could not write to file $vcf_out\n";
}

my ($tot_var, $filtered_qual, $filtered_alt) = (0, 0, 0);
my ($numrep, $num_blacklist, $numsegdup, $num_not_accessible, $num_retro, $num_remission, $num_evs, $num_exac) = (0, 0, 0, 0, 0, 0, 0, 0);
my %qual_num;

while (my $vcfline = $vcf->next_line())
{
	my $x = $vcf->next_data_hash($vcfline);

	$tot_var ++;
	$qual_num{$x->{FILTER}->[0]} = $qual_num{$x->{FILTER}->[0]} ? $qual_num{$x->{FILTER}->[0]} + 1 : 1;

	if (@{$x->{ALT}} != 1) # more than one alternative allele?
	{
		$filtered_alt ++;
		next;
	}		
	
	my $ref_allele = $x->{REF};
	my $alt_allele = $x->{ALT}->[0];
	my $var_type;
	if (length($ref_allele) == length($alt_allele))
	{
		$var_type = 'snp';
	}
	elsif (length($ref_allele) < length($alt_allele))
	{
		$var_type = 'ins';
	}
	else
	{
		$var_type = 'del';
	}
	
	my $status = $x->{FILTER}->[0];
		
	if ($status ne "PASS") # rejected by caller
	{
		$filtered_qual ++;
		next;
	}

	##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">
	if (!defined $x->{gtypes}{$vcf_sample_id}{RD}) {
		print STDERR "ERROR: Could not read RD:\n$vcfline";
		exit(1);
	}
	my $ad_ref = $x->{gtypes}{$vcf_sample_id}{RD};
	
	##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">
	if (!defined $x->{gtypes}{$vcf_sample_id}{AD}) {
		print STDERR "ERROR: Could not read AD:\n$vcfline";
		exit(1);
	}
	my $ad_alt = $x->{gtypes}{$vcf_sample_id}{AD};

	if ($ad_alt == 0) {
		print STDERR "ERROR: AD of alternative allele is zero!\n$vcfline";
		exit(1);
	}
	
	##FORMAT=<ID=ABQ,Number=1,Type=Integer,Description="Average quality of variant-supporting bases (qual2)">
	my $avg_qual = $x->{gtypes}{$vcf_sample_id}{ABQ};
		
	my $freq = $ad_alt / ($ad_alt + $ad_ref);
		

	my (@repeats, @dups, @blacklist, @retro, @rem_samples, %evss, @clinvars);
	my ($chr, $pos) = ("chr".$x->{CHROM}, $x->{POS});

	# ----- annotate overlapping repeat regions
	{
		my $iter = $rmsk->query($chr, $pos-1, $pos+1);
		if ($iter and $iter->{_})
		{
			while (my $line = $rmsk->read($iter)) 
			{
				my @s = split("\t", $line);
				push(@repeats, "$s[10]:$s[11]:$s[12]");
			}
		}		
	}

	{
		my $iter = $simpleRepeat->query($chr, $pos-1, $pos+1);
		if ($iter and $iter->{_})
		{
			while (my $line = $simpleRepeat->read($iter)) 
			{
				my @s = split("\t", $line);
				push(@repeats, "$s[16]($s[6])");
			}
		}		
	}
	$numrep ++ if (@repeats > 0);

	# ----- annotate segmental duplications
	{
		my $iter = $segdup->query($chr, $pos-1, $pos+1);
		if ($iter and $iter->{_})
		{
			while (my $line = $segdup->read($iter)) 
			{
				my @s = split("\t", $line);
				push(@dups, "$s[4]:".sprintf("%.2f", $s[26]));
			}		
		}		
		$numsegdup ++ if (@dups > 0);
	}	
	
	# ----- annotate overlapping DAC blacklisted regions
	{
		my $iter = $blacklistdb->query($chr, $pos-1, $pos+1);
		if ($iter and $iter->{_})
		{
			while (my $line = $blacklistdb->read($iter)) 
			{
				my @s = split("\t", $line);
				push(@blacklist, "$s[4]");
			}		
		}
		$num_blacklist ++ if (@blacklist > 0);		
	}

	# ----- annotate overlapping g1k accessible regions
	my $accessible = "no";
	{
		my $iter = $g1kAcc->query($chr, $pos-1, $pos+1);
		if ($iter and $iter->{_})
		{
			while (my $line = $g1kAcc->read($iter)) 
			{
				$accessible = "";
				last;
			}		
		}
		$num_not_accessible ++ if ($accessible eq "no");		
	}

	# ----- annotate overlapping retrotransposed (pseudo) genes
	{
		my $iter = $ucscRetro->query($chr, $pos-1, $pos+1);
		if ($iter and $iter->{_})
		{
			while (my $line = $ucscRetro->read($iter)) 
			{
				my @s = split("\t", $line);
				push(@retro, $s[10]);
			}		
		}
		$num_retro ++ if (@retro > 0);
	}

	# ----- annotate variants found in normal samples
	if ($remission_variants_file)
	{
		my $iter = $remission->query($x->{CHROM}, $pos-1, $pos);
		if ($iter and $iter->{_})
		{
			while (my $line = $remission->read($iter)) 
			{
				my ($sample, $rchr, $rpos, $ref_allele, $alt_allele, $dp, $ad) = split("\t", $line);
				if ($pos eq $rpos and $x->{REF} eq $ref_allele and $x->{ALT}->[0] eq $alt_allele and ($ad >= 4 or ($ad/$dp >= 0.2 and $ad >= 2)))
				{
					push(@rem_samples, "$sample($ad)");
				}
			}		
		}
		$num_remission ++ if (@rem_samples > 0);
	}

	# ----- annotate variants found in Exome Variant Server
	{
		my $iter = $evs->query($chr, $pos-1, $pos);
		if ($iter and $iter->{_})
		{
			while (my $line = $evs->read($iter)) 
			{
				my ($echr, $epos, $rsID, $dbSNPVersion, $alleles, $europeanAmericanAlleleCount, $africanAmericanAlleleCount, $allAlleleCount, $MAFinPercent_EA_AA_All, $europeanAmericanGenotypeCount, 
					$africanAmericanGenotypeCount, $allGenotypeCount, $avgSampleReadDepth, $genes, $geneAccession, $functionGVS, $hgvsProteinVariant, $hgvsCdnaVariant, $codingDnaSize, 
					$conservationScorePhastCons, $conservationScoreGERP, $granthamScore, $polyphen2_score, $refBaseNCBI37, $chimpAllele, $clinicalInfo, $filterStatus, $onIlluminaHumanExomeChip,
					$gwasPubMedInfo, $EA_EstimatedAge_kyrs, $AA_EstimatedAge_kyrs) = split(/\s/, $line);
					
				next if ($echr ne $chr or $epos ne $pos);
				foreach my $allele (split(";", $alleles))
				{
					my ($ref, $alt) = $allele =~ /(.+)\>(.+)/;
					if ($ref eq $x->{REF} and $alt eq $x->{ALT}->[0])
					{
						my ($alt_count, $ref_count) = $allAlleleCount =~ /(\d+).+?(\d+)/;
						my $alt_percent = sprintf("%.3f", $alt_count/($alt_count+$ref_count)*100);
						$evss{$alt_percent} = 1;					
					}
				}
			}		
		}
		$num_evs ++ if (keys(%evss) > 0);
	}
	
	# ----- G1K frequency
	my $g1k_max_freq = 0;
	foreach my $f ($x->{INFO}{'dbNSFP_1000Gp1_AF'}, $x->{INFO}{'dbNSFP_1000Gp1_AFR_AF'}, $x->{INFO}{'dbNSFP_1000Gp1_EUR_AF'}, $x->{INFO}{'dbNSFP_1000Gp1_AMR_AF'}, $x->{INFO}{'dbNSFP_1000Gp1_ASN_AF'})
	{
		$g1k_max_freq = $f if (defined $f and $g1k_max_freq < $f);
	}
	
	# ----- ESP6500 frequency
	my $evs_freq = (defined $x->{INFO}{dbNSFP_ESP6500_AA_AF} and defined $x->{INFO}{dbNSFP_ESP6500_EA_AF}) 
					? max($x->{INFO}{dbNSFP_ESP6500_AA_AF}, $x->{INFO}{dbNSFP_ESP6500_EA_AF}) 
					: defined $x->{INFO}{dbNSFP_ESP6500_AA_AF} 
						? $x->{INFO}{dbNSFP_ESP6500_AA_AF} 
						: defined $x->{INFO}{dbNSFP_ESP6500_EA_AF} 
							? $x->{INFO}{dbNSFP_ESP6500_EA_AF} 
							: 0;
	my $max_evs_freq = max($evs_freq, keys(%evss) > 0 ? (keys(%evss))[0] : 0);
	$num_evs ++ if ($max_evs_freq > 0);

	# ----- annotate ExAC variant allele frequencies
	my $max_exac_freq = 0;
	{
		my $iter = $exac->query($x->{CHROM}, $pos-1, $pos);
		if ($iter and $iter->{_})
		{
			while (my $line = $exac->read($iter)) 
			{
				my ($rchr, $rpos, $rid, $ref_allele, $alt_allele, $rqual, $rfilter, $rinfo) = split("\t", $line);
				next if ($pos ne $rpos or $x->{REF} ne $ref_allele);
				
				my @alt_alleles = split(",", $alt_allele);
				for (my $i = 0; $i < @alt_alleles; $i ++)
				{
					if ($x->{ALT}->[0] eq $alt_alleles[$i])
					{
						my ($afs) = $rinfo =~ /[;\t]AF=([^;\n]+)/;
						$max_exac_freq = (split(",", $afs))[$i];
						$num_exac ++;
						last;
					}	
				}
				last if ($max_exac_freq);
			}		
		}
	}

	# ----- annotate ClinVar variants
	my $pathogenic = 0;
	{
		my $iter = $clinvar->query($x->{CHROM}, $pos-1, $pos);
		if ($iter and $iter->{_})
		{
			while (my $line = $clinvar->read($iter)) 
			{
				my ($rchr, $rpos, $rid, $ref_allele, $alt_allele, $rqual, $rfilter, $rinfo) = split("\t", $line);
				next if ($pos ne $rpos or $x->{REF} ne $ref_allele);
				
				my @alt_alleles = split(",", $alt_allele);
				for (my $i = 0; $i < @alt_alleles; $i ++)
				{
					if ($x->{ALT}->[0] eq $alt_alleles[$i])
					{
						my ($sigcode) = $rinfo =~ /CLNSIG=(\d+)/;
						my @sigcodes = split(",", $sigcode);
						$sigcode = $i < @sigcodes ? $sigcodes[$i] : $sigcodes[0];
						$pathogenic = 1 if (defined $sigcode and $sigcode == 5);
						$rid .= " (".$clnsig{$sigcode}.")" if (defined $sigcode);
						push(@clinvars, $rid);
					}	
				}
			}
		}		
	}

	my $loc = $x->{CHROM}.":".$x->{POS}."-".$x->{POS};
	$loc =~ s/^chr//;
	
	my @rejected_because;
	my $reject = 0;

	if ($g1k_max_freq > 0.01) { push(@rejected_because, "G1K"); $reject = 1; }
	if ($max_evs_freq > 0.01) { push(@rejected_because, "EVS AF >= 0.01"); $reject = 1; }
	if ($max_exac_freq >= 0.01) { push(@rejected_because, "ExAC AF >= 0.01"); $reject = 1; }

	if ($reject and $pathogenic and $g1k_max_freq < 0.05 and $max_evs_freq < 0.05 and $max_exac_freq < 0.05) { push(@rejected_because, "rescued rare polymorphism because pathogenic"); $reject = 0; }
	
	if (@repeats > 0) { push(@rejected_because, "repetitive region"); $reject = 1; }
	if (@dups > 0) { push(@rejected_because, "segmental duplication"); $reject = 1; }
	if (@blacklist > 0) { push(@rejected_because, "blacklisted region"); $reject = 1; }
	if (@retro > 0) { push(@rejected_because, "retrotransposon"); }
	if (@rem_samples) { push(@rejected_because, "present remissions"); $reject = 1; }
	if ($x->{ID} and $x->{ID} ne ".")  { push(@rejected_because, "dbSNP"); $reject = 1; }  
	
	if ($x->{CHROM} eq "hs37d5") { push(@rejected_because, "decoy genome"); $reject = 1; }

	my ($gene, $add_genes, $impact, $effect, $aa_change) = get_impact($x->{INFO}{EFF});

	# keep known pathogenic CRLF2 mutations overlapping with segmental duplication (PAR1 on Y chromosome!)
	if ($gene eq "CRLF2" && $aa_change eq "F232C") {
		push(@rejected_because, "rescued CRLF2 PAR region");
		$reject = 0;
	} 
	
	$vcfline =~ s/^([^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t)[^\t]+/$1REJECT/ if ($reject);
	$vcfline =~ s/^([^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t)[^\t]+/$1MISSED/ if ($status eq "MISSED");
	
	my $polyphen = $x->{INFO}{'dbNSFP_Polyphen2_HVAR_pred'};
	my $sift = undef;
	if (defined $x->{INFO}{'dbNSFP_SIFT_score'})
	{
		foreach my $s (split(",", $x->{INFO}{'dbNSFP_SIFT_score'}))
		{
			next if (!defined $s or $s eq ".");
			$sift = $s if (!defined $sift or $s < $sift);	
		}
	} 
	my $siphy = $x->{INFO}{'dbNSFP_SiPhy_29way_logOdds'};
	
	my $is_deleterious = "n/d";
	$is_deleterious = "yes" if ($effect eq "NON_SYNONYMOUS_CODING" and $polyphen and $polyphen =~ /D/ and defined $sift and $sift < 0.05); # polyphen & sift
	$is_deleterious = "yes" if ($effect eq "NON_SYNONYMOUS_CODING" and $polyphen and $polyphen =~ /D/ and defined $siphy and $siphy >= 12); # polyphen & siphy
	$is_deleterious = "yes" if ($effect eq "NON_SYNONYMOUS_CODING" and defined $sift and $sift < 0.05 and defined $siphy and $siphy >= 12); # sift and siphy
	$is_deleterious = "yes" if ($effect eq "NON_SYNONYMOUS_CODING" and defined $siphy and $siphy > 20); # siphy only, highly conserved (keeps GNAQ)
	$is_deleterious = "yes" if ($effect eq "FRAME_SHIFT" or $effect eq "SPLICE_SITE_ACCEPTOR" or $effect eq "SPLICE_SITE_DONOR" or $effect eq "STOP_GAINED");
	$is_deleterious = "no" if ($is_deleterious ne "yes" and defined $polyphen and defined $sift);
	$is_deleterious = "no" if ($effect eq "DOWNSTREAM" or $effect eq "UPSTREAM" or $effect eq "INTRON" or $effect eq "INTERGENIC" or $effect eq "SYNONYMOUS_CODING" or $effect eq "SYNONYMOUS_STOP" or $effect eq "SYNONYMOUS_START" or $effect eq "UTR_3_PRIME" or $effect eq "UTR_5_PRIME" or $effect eq "UTR_5_DELETED" or $effect eq "UTR_3_DELETED" or $effect eq "START_GAINED");
	
	my $non_silent = 0;
	$non_silent = 1 if ($effect eq "STOP_GAINED" or $effect eq "STOP_LOST" or $effect eq "SPLICE_SITE_DONOR" or $effect eq "SPLICE_SITE_ACCEPTOR" or $effect eq "FRAME_SHIFT" or $effect eq "CODON_CHANGE_PLUS_CODON_INSERTION" or $effect eq "CODON_DELETION" or $effect eq "NON_SYNONYMOUS_CODING" or $effect eq "CODON_INSERTION" or $effect eq "CODON_CHANGE_PLUS_CODON_DELETION" or $effect eq "NON_SYNONYMOUS_START" or $effect eq "START_LOST");
	
	print VCFOUT "$vcfline" if ($vcf_out);
	
	print "$sample_identifier\t";		
	print "$var_type\t";
	print $reject ? "REJECT" : "$status", "\t";
	print @rejected_because > 0 ? join(";", @rejected_because) : "", "\t";
	print $x->{CHROM},"\t";
	print $x->{POS},"\t";
	print $x->{REF},"\t";
	print $x->{ALT}->[0],"\t";
	print "$gene\t";
	print "$add_genes\t";
	print $x->{ID},"\t";
	print $g1k_max_freq > 0 ? $g1k_max_freq : "", "\t";  # alternative allele frequency in the whole 1000Gp1 data
	print $evs_freq > 0 ? $evs_freq : "", "\t"; 
	print $max_exac_freq > 0 ? sprintf("%.4f", $max_exac_freq) : 0, "\t";
	print @clinvars > 0 ? join(",", @clinvars) : "", "\t";
	print "$impact\t";
	print "$effect\t";
	print "$non_silent\t";
	print "$is_deleterious\t";
	print $x->{INFO}{SNPEFF_EXON_ID} ? $x->{INFO}{SNPEFF_EXON_ID} : "", "\t"; 
	print "".($ad_ref+$ad_alt)."\t";
	print "$ad_ref\t";
	print "$ad_alt\t";
	print "$freq\t";
	print "$avg_qual\t";
	print "$aa_change\t";
	print "EFF=",$x->{INFO}{EFF},"\t";
	print defined $polyphen ? $polyphen : "", "\t"; # Polyphen2 prediction based on HumVar, 'D' ('porobably damaging'), 'P' ('possibly damaging') and 'B' ('benign'). Multiple entries separated by ';' 
	print defined $sift ? $sift : "", "\t"; # SIFT score, If a score is smaller than 0.05 the corresponding NS is predicted as 'D(amaging)'; otherwise it is predicted as 'T(olerated)'
	print defined $x->{INFO}{'dbNSFP_GERP++_RS'} ? $x->{INFO}{'dbNSFP_GERP++_RS'} : "", "\t"; # GERP++ RS score, the larger the score, the more conserved the site 
	print defined $siphy ? $siphy : "", "\t"; # SiPhy score based on 29 mammals genomes. The larger the score, the more conserved the site.
	my $domains = $x->{INFO}{'dbNSFP_Interpro_domain'}; # domain or conserved site on which the variant locates. Domain annotations come from Interpro database. The number in the brackets following a specific domain is the count of times Interpro assigns the variant position to that domain, typically coming from different predicting databases
	if ($domains)
	{
		$domains =~ s/\),/\)\|/g;
		$domains =~ s/\|$//;
		print "$domains\t";
	}
	else
	{
		print "\t";
	}
	print $cosmic{$loc} ? $cosmic{$loc} : "0", "\t";
	print aa_hits([$gene, split(",", $add_genes)], "EFF=".$x->{INFO}{EFF}), "\t";	
	print join(',', @repeats), "\t", join(',', @dups), "\t", join(',', @blacklist), "\t$accessible\t", join(",", @retro), "\t", join(",", @rem_samples);
	print "\n";
			
#	print "\n"; print Dumper($x); exit;
}

$vcf->close();
close(VCFOUT) if ($vcf_out);
	
if ($debug)
{
	print STDERR "  Total number of variants: $tot_var\n";
	print STDERR "  Variants by quality:\n";
	foreach my $k (keys(%qual_num))
	{
		print STDERR "    $k: ", $qual_num{$k}, "\n";
	}
	print STDERR "  Rejected by caller: $filtered_qual\n";
	print STDERR "  Excluded due to missing alternative allele: $filtered_alt\n";
	print STDERR "  $numrep variants annotated with overlapping repeat.\n";
	print STDERR "  $num_blacklist variants annotated with overlapping blacklisted region.\n";
	print STDERR "  $numsegdup variants annotated with overlapping segmental duplication.\n";
	print STDERR "  $num_not_accessible variants fall into G1K non-accessible region.\n";
	print STDERR "  $num_retro variants annotated with overlapping retrotransposed (pseudo)gene.\n";
	print STDERR "  $num_remission variants present in remission sample(s).\n";
	print STDERR "  $num_evs variants present in Exome Variant Server (EVS).\n";
	print STDERR "  $num_exac variants present in Exome Aggregation Consortium (ExAC) dataset.\n";
}

# ------------------------------------------

sub get_impact
{
	my $effs = shift or die "ERROR: effect not specified";

	# determine all genes impacted by variants
	my (%genes_by_impact, %all_genes, $combined_effect, $combined_impact, %aa_changes);
	foreach my $eff (split(",", $effs))
	{
		my ($effect, $rest) = $eff =~ /([^\(]+)\(([^\)]+)\)/
			or die "ERROR: could not parse SNP effect: $effs\n";

		my ($impact, $class, $codon_change, $aa_change, $aa_length, $gene_name, $gene_biotype, 
			$coding, $transcript, $exon, $genotype_num) = split('\|', $rest)
				or die "ERROR: could not parse SNP effect: $eff\n";
		 
		$aa_changes{$aa_change} = 1 if ($aa_change);

		if ($exon and $transcript and $gene_name)
		{
			$transcript =~ s/\.\d+$//; # remove version number from accession
			$transcript =~ s/\.\d+$//; 
		}
			
		# gene impacted by variant?
		if ($gene_name)
		{
			$genes_by_impact{$impact}{$gene_name} = $effect;
			$all_genes{$gene_name} = 1;
		}

		$combined_impact = $impact;		
		$combined_effect = $effect;
	}
	
	# if multiple genes are affected, preferentially chose gene with the predicted higher impact
	if ($genes_by_impact{HIGH})
	{
		$combined_impact = "HIGH";
	}
	elsif ($genes_by_impact{MODERATE})
	{
		$combined_impact = "MODERATE";
	}
	elsif ($genes_by_impact{LOW})
	{
		$combined_impact = "LOW";
	}
	elsif ($genes_by_impact{MODIFIER})
	{
		$combined_impact = "MODIFIER";
	}
	
	my ($gene, $add_genes) = ("", "");
	if (keys(%all_genes) > 0)
	{
		my @sorted_genes = sort keys(%{$genes_by_impact{$combined_impact}});
		$gene = $sorted_genes[0]; # first choice is first in alphabetically sorted list
		if ($gene =~ /^LOC/) # if this is a generic gene name, try to find non-generic one instead
		{
			foreach my $g (@sorted_genes)
			{
				if ($g !~ /^LOC/)
				{
					$gene = $g;
					last;
				}	
			}
		}
		$combined_effect = $genes_by_impact{$combined_impact}{$gene};
		delete $all_genes{$gene};
		$add_genes = join(",", keys(%all_genes));
	}

	return ($gene, $add_genes, $combined_impact, $combined_effect, join(";", keys(%aa_changes)));
}

sub aa_hits
{
	my $genes = shift;
	my $snpeff = shift;
	my $leuk = shift;
	
	return "non-coding" if (@$genes == 0);
	
	my $aa_change_found = 0; 
	foreach my $gene (@$genes) # check each gene
	{
		foreach my $eff (split(",", $snpeff)) # check all isoforms for cosmic match
		{
			my ($effect, $rest) = $eff =~ /([^\(]+)\(([^\)]+)\)/
				or croak "ERROR: could not parse SNP effect: $snpeff";
	
			my ($impact, $class, $codon_change, $aa_change, $aa_length, $gene_name, $gene_biotype, 
				$coding, $transcript, $exon, $genotype_num) = split('\|', $rest)
					or croak "ERROR: could not parse SNP effect: $eff";
					 
			if ($aa_change =~ /(.)(\d+)(.+)/)
			{
				$aa_change_found = 1;			
				my ($prev_aa, $aa_pos, $after_aa) = ($1, $2, $3);
				return $cosmic{"$gene:$prev_aa:$aa_pos"} if (!$leuk and defined $cosmic{"$gene:$prev_aa:$aa_pos"});
				return $cosmic_leuk{"$gene:$prev_aa:$aa_pos"} if ($leuk and defined $cosmic_leuk{"$gene:$prev_aa:$aa_pos"});
			}
		}				
	}
	
	return "non-coding" if (!$aa_change_found);
	return "0";
}
