use strict;
use warnings;

my $outfile = "/mnt/projects/iamp/results/geneset_annotation.tsv";
open(OUT, ">$outfile") or die "ERROR: Could not write to output file $outfile.\n";

print OUT "geneset\tcategory\tdescription\tlinkout\n";

# MSigDB 5.0

open(IN, "/mnt/projects/generic/data/msigdb5.0/h.all.v5.0.symbols.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	print OUT uc($gs),"\tMSIGDB_H_HALLMARK\t\thttp://www.broadinstitute.org/gsea/msigdb/cards/$gs.html\n";
}
close(IN);

open(IN, "/mnt/projects/generic/data/msigdb5.0/c1.all.v5.0.symbols.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	print OUT uc($gs),"\tMSIGDB_C1_POSITIONAL\t\thttp://www.broadinstitute.org/gsea/msigdb/cards/$gs.html\n";
}
close(IN);

open(IN, "/mnt/projects/generic/data/msigdb5.0/c2.cgp.v5.0.symbols.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	my $cat = "MSIGDB_C2_CGP";
	$cat .= ",HEMATOLOGIC" if ($gs =~ /(MULLIGHAN_MLL_SIGNATURE|HEMATOP)/);
	print OUT uc($gs),"\t$cat\t\thttp://www.broadinstitute.org/gsea/msigdb/cards/$gs.html\n";
}
close(IN);

open(IN, "/mnt/projects/generic/data/msigdb5.0/c2.cp.v5.0.symbols.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	print OUT uc($gs),"\tMSIGDB_C2_PATHWAYS\t\thttp://www.broadinstitute.org/gsea/msigdb/cards/$gs.html\n";
}
close(IN);

open(IN, "/mnt/projects/generic/data/msigdb5.0/c3.mir.v5.0.symbols.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	print OUT uc($gs),"\tMSIGDB_C3_MOTIF_MIR\t\thttp://www.broadinstitute.org/gsea/msigdb/cards/$gs.html\n";
}
close(IN);

open(IN, "/mnt/projects/generic/data/msigdb5.0/c3.tft.v5.0.symbols.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	print OUT uc($gs),"\tMSIGDB_C3_MOTIF_TFT\t\thttp://www.broadinstitute.org/gsea/msigdb/cards/$gs.html\n";
}
close(IN);

open(IN, "/mnt/projects/generic/data/msigdb5.0/c4.cgn.v5.0.symbols.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	print OUT uc($gs),"\tMSIGDB_C4_COMPUTATIONAL_CGN\t\thttp://www.broadinstitute.org/gsea/msigdb/cards/$gs.html\n";
}
close(IN);

open(IN, "/mnt/projects/generic/data/msigdb5.0/c4.cm.v5.0.symbols.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	print OUT uc($gs),"\tMSIGDB_C4_COMPUTATIONAL_CM\t\thttp://www.broadinstitute.org/gsea/msigdb/cards/$gs.html\n";
}
close(IN);

open(IN, "/mnt/projects/generic/data/msigdb5.0/c5.bp.v5.0.symbols.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	print OUT uc($gs),"\tMSIGDB_C5_GO_BP\t\thttp://www.broadinstitute.org/gsea/msigdb/cards/$gs.html\n";
}
close(IN);

open(IN, "/mnt/projects/generic/data/msigdb5.0/c5.cc.v5.0.symbols.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	print OUT uc($gs),"\tMSIGDB_C5_GO_CC\t\thttp://www.broadinstitute.org/gsea/msigdb/cards/$gs.html\n";
}
close(IN);

open(IN, "/mnt/projects/generic/data/msigdb5.0/c5.mf.v5.0.symbols.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	print OUT uc($gs),"\tMSIGDB_C5_GO_MF\t\thttp://www.broadinstitute.org/gsea/msigdb/cards/$gs.html\n";
}
close(IN);

open(IN, "/mnt/projects/generic/data/msigdb5.0/c6.all.v5.0.symbols.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	print OUT uc($gs),"\tMSIGDB_C6_ONCOGENIC_SIGNATURES\t\thttp://www.broadinstitute.org/gsea/msigdb/cards/$gs.html\n";
}
close(IN);

open(IN, "/mnt/projects/generic/data/msigdb5.0/c7.all.v5.0.symbols.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	my $cat = "MSIGDB_C7_IMMUNOLOGIC_SIGNATURES";
	$cat .= ",HEMATOLOGIC" if ($gs =~ /(BCELL_VS|VS_BCELL|TCELL_VS|VS_TCELL|NEUTROPHIL_VS|VS_MONOCYTE|_DC_VS)/ and $gs !~ /LUPUS|_FLU_|CD4_TCELL|CD8_TCELL|GONDII|IGG_IGA|DONOVANI|_STIM_|_MALAYI_/);
	print OUT uc($gs),"\t$cat\t\thttp://www.broadinstitute.org/gsea/msigdb/cards/$gs.html\n";
}
close(IN);

# DSigDB 1.0

my %dsigdb;
open(IN, "/mnt/projects/generic/data/DSigDB/DSigDB_All.txt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	if (/^Compound\(([^\)]+)\) : (.*)/) {
		my ($gs, $cat) = ($2, $1);
		$cat =~ s/ /_/g;
		
		if ($cat eq "TTD") {
			$cat = "D4_TTD";
		}
		elsif ($cat eq "CTD") {
			$cat = "D4_CTD";
		}
		elsif ($cat eq "BOSS") {
			$cat = "D4_BOSS";
		}
		elsif ($cat eq "D1") {
			$cat = "D1_FDA";
		}
		elsif ($cat eq "D3") {
			$cat = "D3_CMAP_PERTURBATION";
		}
		$dsigdb{$gs} = $cat;
	}
#	elsif (/^Compound : ([^\(]+)\(([^\)]+)\)/) {
#		my ($gs, $cat) = ($1, $2);
#		$cat =~ s/ /_/g;
#		$dsigdb{$gs} = $cat;
#	}
}
close(IN);

open(IN, "/mnt/projects/generic/data/DSigDB/DSigDB_v1.0_All.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs, $url) = split /\t/;
	next if (!$gs);

	my $gs_short = $gs;
	$gs_short =~ s/_TTD_\d+//;
	$gs_short =~ s/_CTD_\d+//;
		
	if ($gs =~ /(.*)_MRC$/) {
		print OUT uc($gs),"\tDSIGDB_D2_KINASE_INH_MRC\t\t$url\n";
		next;
	}
	elsif ($gs =~ /(.*)_LINCS$/) {
		print OUT uc($gs),"\tDSIGDB_D2_KINASE_INH_LINCS\t\t$url\n";
		next;
	}
	elsif ($gs =~ /(.*)_Roche$/) {
		print OUT uc($gs),"\tDSIGDB_D2_KINASE_INH_ROCHE\t\t$url\n";
		next;
	}
	elsif ($gs =~ /(.*)_GSK$/) {
		print OUT uc($gs),"\tDSIGDB_D2_KINASE_INH_GSK\t\t$url\n";
		next;
	}
	elsif ($gs =~ /(.*)_FDA$/) {
		print OUT uc($gs),"\tDSIGDB_D2_KINASE_INH_FDA\t\t$url\n";
		next;
	}
	elsif ($gs =~ /(.*)_RBC$/) {
		print OUT uc($gs),"\tDSIGDB_D2_KINASE_INH_RBC\t\t$url\n";
		next;
	}
	elsif ($gs =~ /(.*)_Kinome Scan$/) {
		print OUT uc($gs),"\tDSIGDB_D2_KINASE_INH_KINOMESCAN\t\t$url\n";
		next;
	}
	elsif ($gs =~ /(.*)_BOSS$/) {
		print OUT uc($gs),"\tDSIGDB_D4_BOSS\t\t$url\n";
		next;
	}
	elsif ($dsigdb{$gs_short}) {
		print OUT uc($gs),"\tDSIGDB_",$dsigdb{$gs_short}, "\t\t$url\n";
	}
	else {
		print STDERR "WARNING: Category for DSigDB gene set $gs unknown.\n";
	}
}
close(IN);

# GeneSigDB 4.0

my %id2title;
open(IN, "/mnt/projects/generic/data/GeneSigDB/ALL_SIGSv4.nodup.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs) = split /\t/;
	next if (!$gs);

	my ($pmid, $rest) = $gs =~ /([^-_]+)[-_](.*)/;
	my $xmlfile = "/mnt/projects/generic/data/GeneSigDB/details/$pmid/$gs-index.xml";
	
	# fetch description
	my $desc = "";
	open(F, $xmlfile) or die "ERROR: Could not open file $xmlfile\n";
	while (<F>) {
		if (/<sigDescription>(.*)<\/sigDescription>/) {
			$desc = $1;
			last;
		}
	}
	close(F);
	
	# fetch paper title
	my $title = "";
	if (!defined $id2title{$pmid}) {
		print "Fetching pubmed ID $pmid... ";
		open(W, "curl -s \"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=$pmid&rettype=medline\"|") or die "ERROR curl: $!\n";
		while(<W>) {
			chomp;
			my ($key, $value) = /(\S+)\s*?-\s+(.*)/;
			if (!$key and $title) {
				my ($nline) = /^\s*(.*)/;			
				$title .= $nline;
			} elsif ($title) {
				last;
			} elsif ($key and $key eq "TI") {
				$title = $value;
			}
		}
		close(W);		
		print "$title\n";
		$id2title{$pmid} = $title;
	} else {
		$title = $id2title{$pmid};
	}
	
	print OUT uc($gs),"\tGeneSigDB\t$title $desc\thttp://compbio.dfci.harvard.edu/genesigdb/signaturedetail.jsp?signatureId=$gs\n";
}
close(IN);

# Jaspar

open(IN, "/mnt/projects/generic/data/opossum3/jaspar_core.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs, $url) = split /\t/;
	next if (!$gs);

	print OUT uc($gs),"\tJaspar\t\t$url\n";
}
close(IN);

# PAZAR

#open(IN, "/mnt/projects/generic/data/pazar/pazar.gmt") or die "ERROR: Could not read input.\n";
#while (<IN>) {
#	my ($gs, $url) = split /\t/;
#	next if (!$gs);
#
#	print OUT uc($gs),"\tPAZAR\t\t$url\n";
#}
#close(IN);

# Laurenti et al. (2013)
# http://www.nature.com.ez.srv.meduniwien.ac.at/ni/journal/v14/n7/full/ni.2615.html

open(IN, "/mnt/projects/generic/data/laurenti_2013_hematopoietic_lineages.gmt") or die "ERROR: Could not read input.\n";
while (<IN>) {
	my ($gs, $url) = split /\t/;
	next if (!$gs);
	$url = "http://www.ncbi.nlm.nih.gov/pubmed/23708252";
	print OUT uc($gs),"\tHEMATOLOGIC\t\t$url\n";
}
close(IN);

close(OUT);

print "Success. Output written to $outfile\n"; 

