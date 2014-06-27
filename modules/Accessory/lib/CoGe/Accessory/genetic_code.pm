package CoGe::Accessory::genetic_code;

use strict;
use POSIX;
BEGIN {
    use Exporter ();
    use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $code $aa_info);
    $VERSION     = 0.1;
    @ISA         = (@ISA, qw (Exporter));
    #Give a hoot don't pollute, do not export more than needed by default
    @EXPORT      = qw (code all_codes);
    @EXPORT_OK   = qw ();
    %EXPORT_TAGS = ();
    $code = {
	     'transl_table=22' => {
				   'name' => 'The Scenedesmus obliquus mitochondrial Code (transl_table=22)',
				   'code' => {
					      'GCC' => 'A',
					      'AGT' => 'S',
					      'TGA' => '*',
					      'TGT' => 'C',
					      'CGA' => 'R',
					      'ATC' => 'I',
					      'AAC' => 'N',
					      'AGC' => 'S',
					      'TAC' => 'Y',
					      'ACA' => 'T',
					      'TCG' => 'S',
					      'CCG' => 'P',
					      'CTG' => 'L',
					      'GCA' => 'A',
					      'GTG' => 'V',
					      'AAG' => 'K',
					      'GTT' => 'V',
					      'CAC' => 'H',
					      'AGA' => 'R',
					      'ACC' => 'T',
					      'CCA' => 'P',
					      'TGG' => 'W',
					      'CTC' => 'L',
					      'CGC' => 'R',
					      'TTG' => 'L',
					      'TAA' => '*',
					      'CAG' => 'Q',
					      'ACG' => 'T',
					      'AAA' => 'K',
					      'ATG' => 'M',
					      'GTA' => 'V',
					      'TAG' => 'L',
					      'CTT' => 'L',
					      'GGA' => 'G',
					      'GTC' => 'V',
					      'TGC' => 'C',
					      'TCA' => '*',
					      'ATT' => 'I',
					      'TAT' => 'Y',
					      'AAT' => 'N',
					      'ACT' => 'T',
					      'CAA' => 'Q',
					      'GAC' => 'D',
					      'GGT' => 'G',
					      'TCC' => 'S',
					      'TTT' => 'F',
					      'AGG' => 'R',
					      'CGT' => 'R',
					      'CGG' => 'R',
					      'CAT' => 'H',
					      'ATA' => 'I',
					      'CCC' => 'P',
					      'GGG' => 'G',
					      'TTA' => 'L',
					      'GAG' => 'E',
					      'CTA' => 'L',
					      'GAT' => 'D',
					      'TCT' => 'S',
					      'TTC' => 'F',
					      'GCG' => 'A',
					      'GGC' => 'G',
					      'GCT' => 'A',
					      'GAA' => 'E',
					      'CCT' => 'P'
					     }
				  },
	     'transl_table=15' => {
				   'name' => 'The Blepharisma Nuclear Code (transl_table=15)',
				   'code' => {
					      'GCC' => 'A',
					      'AGT' => 'S',
					      'TGA' => '*',
					      'TGT' => 'C',
					      'CGA' => 'R',
					      'ATC' => 'I',
					      'AAC' => 'N',
					      'AGC' => 'S',
					      'TAC' => 'Y',
					      'ACA' => 'T',
					      'TCG' => 'S',
					      'CCG' => 'P',
					      'CTG' => 'L',
					      'GCA' => 'A',
					      'GTG' => 'V',
					      'AAG' => 'K',
					      'GTT' => 'V',
					      'CAC' => 'H',
					      'AGA' => 'R',
					      'ACC' => 'T',
					      'CCA' => 'P',
					      'TGG' => 'W',
					      'CTC' => 'L',
					      'CGC' => 'R',
					      'TTG' => 'L',
					      'TAA' => '*',
					      'CAG' => 'Q',
					      'ACG' => 'T',
					      'AAA' => 'K',
					      'ATG' => 'M',
					      'GTA' => 'V',
					      'TAG' => 'Q',
					      'CTT' => 'L',
					      'GGA' => 'G',
					      'GTC' => 'V',
					      'TGC' => 'C',
					      'TCA' => 'S',
					      'ATT' => 'I',
					      'TAT' => 'Y',
					      'AAT' => 'N',
					      'ACT' => 'T',
					      'CAA' => 'Q',
					      'GAC' => 'D',
					      'GGT' => 'G',
					      'TCC' => 'S',
					      'TTT' => 'F',
					      'AGG' => 'R',
					      'CGT' => 'R',
					      'CGG' => 'R',
					      'CAT' => 'H',
					      'ATA' => 'I',
					      'CCC' => 'P',
					      'GGG' => 'G',
					      'TTA' => 'L',
					      'GAG' => 'E',
					      'CTA' => 'L',
					      'GAT' => 'D',
					      'TCT' => 'S',
					      'TTC' => 'F',
					      'GCG' => 'A',
					      'GGC' => 'G',
					      'GCT' => 'A',
					      'GAA' => 'E',
					      'CCT' => 'P'
					     }
				  },
	     'transl_table=9' => {
				  'name' => 'The Echinoderm and Flatworm Mitochondrial Code (transl_table=9)',
				  'code' => {
					     'GCC' => 'A',
					     'AGT' => 'S',
					     'TGA' => 'W',
					     'TGT' => 'C',
					     'CGA' => 'R',
					     'ATC' => 'I',
					     'AAC' => 'N',
					     'AGC' => 'S',
					     'TAC' => 'Y',
					     'ACA' => 'T',
					     'TCG' => 'S',
					     'CCG' => 'P',
					     'CTG' => 'L',
					     'GCA' => 'A',
					     'GTG' => 'V',
					     'AAG' => 'K',
					     'GTT' => 'V',
					     'CAC' => 'H',
					     'AGA' => 'S',
					     'ACC' => 'T',
					     'CCA' => 'P',
					     'TGG' => 'W',
					     'CTC' => 'L',
					     'CGC' => 'R',
					     'TTG' => 'L',
					     'TAA' => '*',
					     'CAG' => 'Q',
					     'ACG' => 'T',
					     'AAA' => 'N',
					     'ATG' => 'M',
					     'GTA' => 'V',
					     'TAG' => '*',
					     'CTT' => 'L',
					     'GGA' => 'G',
					     'GTC' => 'V',
					     'TGC' => 'C',
					     'TCA' => 'S',
					     'ATT' => 'I',
					     'TAT' => 'Y',
					     'AAT' => 'N',
					     'ACT' => 'T',
					     'CAA' => 'Q',
					     'GAC' => 'D',
					     'GGT' => 'G',
					     'TCC' => 'S',
					     'TTT' => 'F',
					     'AGG' => 'S',
					     'CGT' => 'R',
					     'CGG' => 'R',
					     'CAT' => 'H',
					     'ATA' => 'I',
					     'CCC' => 'P',
					     'GGG' => 'G',
					     'TTA' => 'L',
					     'GAG' => 'E',
					     'CTA' => 'L',
					     'GAT' => 'D',
					     'TCT' => 'S',
					     'TTC' => 'F',
					     'GCG' => 'A',
					     'GGC' => 'G',
					     'GCT' => 'A',
					     'GAA' => 'E',
					     'CCT' => 'P'
					    }
				 },
	     'transl_table=12' => {
				   'name' => 'The Alternative Yeast Nuclear Code (transl_table=12)',
				   'code' => {
					      'GCC' => 'A',
					      'AGT' => 'S',
					      'TGA' => '*',
					      'TGT' => 'C',
					      'CGA' => 'R',
					      'ATC' => 'I',
					      'AAC' => 'N',
					      'AGC' => 'S',
					      'TAC' => 'Y',
					      'ACA' => 'T',
					      'TCG' => 'S',
					      'CCG' => 'P',
					      'CTG' => 'S',
					      'GCA' => 'A',
					      'GTG' => 'V',
					      'AAG' => 'K',
					      'GTT' => 'V',
					      'CAC' => 'H',
					      'AGA' => 'R',
					      'ACC' => 'T',
					      'CCA' => 'P',
					      'TGG' => 'W',
					      'CTC' => 'L',
					      'CGC' => 'R',
					      'TTG' => 'L',
					      'TAA' => '*',
					      'CAG' => 'Q',
					      'ACG' => 'T',
					      'AAA' => 'K',
					      'ATG' => 'M',
					      'GTA' => 'V',
					      'TAG' => '*',
					      'CTT' => 'L',
					      'GGA' => 'G',
					      'GTC' => 'V',
					      'TGC' => 'C',
					      'TCA' => 'S',
					      'ATT' => 'I',
					      'TAT' => 'Y',
					      'AAT' => 'N',
					      'ACT' => 'T',
					      'CAA' => 'Q',
					      'GAC' => 'D',
					      'GGT' => 'G',
					      'TCC' => 'S',
					      'TTT' => 'F',
					      'AGG' => 'R',
					      'CGT' => 'R',
					      'CGG' => 'R',
					      'CAT' => 'H',
					      'ATA' => 'I',
					      'CCC' => 'P',
					      'GGG' => 'G',
					      'TTA' => 'L',
					      'GAG' => 'E',
					      'CTA' => 'L',
					      'GAT' => 'D',
					      'TCT' => 'S',
					      'TTC' => 'F',
					      'GCG' => 'A',
					      'GGC' => 'G',
					      'GCT' => 'A',
					      'GAA' => 'E',
					      'CCT' => 'P'
					     }
				  },
	     'transl_table=2' => {
				  'name' => 'The Vertebrate Mitochondrial Code (transl_table=2)',
				  'code' => {
					     'GCC' => 'A',
					     'AGT' => 'S',
					     'TGA' => 'W',
					     'TGT' => 'C',
					     'CGA' => 'R',
					     'ATC' => 'I',
					     'AAC' => 'N',
					     'AGC' => 'S',
					     'TAC' => 'Y',
					     'ACA' => 'T',
					     'TCG' => 'S',
					     'CCG' => 'P',
					     'CTG' => 'L',
					     'GCA' => 'A',
					     'GTG' => 'V',
					     'AAG' => 'K',
					     'GTT' => 'V',
					     'CAC' => 'H',
					     'AGA' => '*',
					     'ACC' => 'T',
					     'CCA' => 'P',
					     'TGG' => 'W',
					     'CTC' => 'L',
					     'CGC' => 'R',
					     'TTG' => 'L',
					     'TAA' => '*',
					     'CAG' => 'Q',
					     'ACG' => 'T',
					     'AAA' => 'K',
					     'ATG' => 'M',
					     'GTA' => 'V',
					     'TAG' => '*',
					     'CTT' => 'L',
					     'GGA' => 'G',
					     'GTC' => 'V',
					     'TGC' => 'C',
					     'TCA' => 'S',
					     'ATT' => 'I',
					     'TAT' => 'Y',
					     'AAT' => 'N',
					     'ACT' => 'T',
					     'CAA' => 'Q',
					     'GAC' => 'D',
					     'GGT' => 'G',
					     'TCC' => 'S',
					     'TTT' => 'F',
					     'AGG' => '*',
					     'CGT' => 'R',
					     'CGG' => 'R',
					     'CAT' => 'H',
					     'ATA' => 'M',
					     'CCC' => 'P',
					     'GGG' => 'G',
					     'TTA' => 'L',
					     'GAG' => 'E',
					     'CTA' => 'L',
					     'GAT' => 'D',
					     'TCT' => 'S',
					     'TTC' => 'F',
					     'GCG' => 'A',
					     'GGC' => 'G',
					     'GCT' => 'A',
					     'GAA' => 'E',
					     'CCT' => 'P'
					    }
				 },
	     'transl_table=13' => {
				   'name' => 'The Ascidian Mitochondrial Code (transl_table=13)',
				   'code' => {
					      'GCC' => 'A',
					      'AGT' => 'S',
					      'TGA' => 'W',
					      'TGT' => 'C',
					      'CGA' => 'R',
					      'ATC' => 'I',
					      'AAC' => 'N',
					      'AGC' => 'S',
					      'TAC' => 'Y',
					      'ACA' => 'T',
					      'TCG' => 'S',
					      'CCG' => 'P',
					      'CTG' => 'L',
					      'GCA' => 'A',
					      'GTG' => 'V',
					      'AAG' => 'K',
					      'GTT' => 'V',
					      'CAC' => 'H',
					      'AGA' => 'G',
					      'ACC' => 'T',
					      'CCA' => 'P',
					      'TGG' => 'W',
					      'CTC' => 'L',
					      'CGC' => 'R',
					      'TTG' => 'L',
					      'TAA' => '*',
					      'CAG' => 'Q',
					      'ACG' => 'T',
					      'AAA' => 'K',
					      'ATG' => 'M',
					      'GTA' => 'V',
					      'TAG' => '*',
					      'CTT' => 'L',
					      'GGA' => 'G',
					      'GTC' => 'V',
					      'TGC' => 'C',
					      'TCA' => 'S',
					      'ATT' => 'I',
					      'TAT' => 'Y',
					      'AAT' => 'N',
					      'ACT' => 'T',
					      'CAA' => 'Q',
					      'GAC' => 'D',
					      'GGT' => 'G',
					      'TCC' => 'S',
					      'TTT' => 'F',
					      'AGG' => 'G',
					      'CGT' => 'R',
					      'CGG' => 'R',
					      'CAT' => 'H',
					      'ATA' => 'M',
					      'CCC' => 'P',
					      'GGG' => 'G',
					      'TTA' => 'L',
					      'GAG' => 'E',
					      'CTA' => 'L',
					      'GAT' => 'D',
					      'TCT' => 'S',
					      'TTC' => 'F',
					      'GCG' => 'A',
					      'GGC' => 'G',
					      'GCT' => 'A',
					      'GAA' => 'E',
					      'CCT' => 'P'
					     }
				  },
	     'transl_table=23' => {
				   'name' => 'The Thraustochytrium Mitochondrial Code (transl_table=23)',
				   'code' => {
					      'GCC' => 'A',
					      'AGT' => 'S',
					      'TGA' => '*',
					      'TGT' => 'C',
					      'CGA' => 'R',
					      'ATC' => 'I',
					      'AAC' => 'N',
					      'AGC' => 'S',
					      'TAC' => 'Y',
					      'ACA' => 'T',
					      'TCG' => 'S',
					      'CCG' => 'P',
					      'CTG' => 'L',
					      'GCA' => 'A',
					      'GTG' => 'V',
					      'AAG' => 'K',
					      'GTT' => 'V',
					      'CAC' => 'H',
					      'AGA' => 'R',
					      'ACC' => 'T',
					      'CCA' => 'P',
					      'TGG' => 'W',
					      'CTC' => 'L',
					      'CGC' => 'R',
					      'TTG' => 'L',
					      'TAA' => '*',
					      'CAG' => 'Q',
					      'ACG' => 'T',
					      'AAA' => 'K',
					      'ATG' => 'M',
					      'GTA' => 'V',
					      'TAG' => '*',
					      'CTT' => 'L',
					      'GGA' => 'G',
					      'GTC' => 'V',
					      'TGC' => 'C',
					      'TCA' => 'S',
					      'ATT' => 'I',
					      'TAT' => 'Y',
					      'AAT' => 'N',
					      'ACT' => 'T',
					      'CAA' => 'Q',
					      'GAC' => 'D',
					      'GGT' => 'G',
					      'TCC' => 'S',
					      'TTT' => 'F',
					      'AGG' => 'R',
					      'CGT' => 'R',
					      'CGG' => 'R',
					      'CAT' => 'H',
					      'ATA' => 'I',
					      'CCC' => 'P',
					      'GGG' => 'G',
					      'TTA' => '*',
					      'GAG' => 'E',
					      'CTA' => 'L',
					      'GAT' => 'D',
					      'TCT' => 'S',
					      'TTC' => 'F',
					      'GCG' => 'A',
					      'GGC' => 'G',
					      'GCT' => 'A',
					      'GAA' => 'E',
					      'CCT' => 'P'
					     }
				  },
	     'transl_table=6' => {
				  'name' => 'The Ciliate, Dasycladacean and Hexamita Nuclear Code (transl_table=6)',
				  'code' => {
					     'GCC' => 'A',
					     'AGT' => 'S',
					     'TGA' => '*',
					     'TGT' => 'C',
					     'CGA' => 'R',
					     'ATC' => 'I',
					     'AAC' => 'N',
					     'AGC' => 'S',
					     'TAC' => 'Y',
					     'ACA' => 'T',
					     'TCG' => 'S',
					     'CCG' => 'P',
					     'CTG' => 'L',
					     'GCA' => 'A',
					     'GTG' => 'V',
					     'AAG' => 'K',
					     'GTT' => 'V',
					     'CAC' => 'H',
					     'AGA' => 'R',
					     'ACC' => 'T',
					     'CCA' => 'P',
					     'TGG' => 'W',
					     'CTC' => 'L',
					     'CGC' => 'R',
					     'TTG' => 'L',
					     'TAA' => 'Q',
					     'CAG' => 'Q',
					     'ACG' => 'T',
					     'AAA' => 'K',
					     'ATG' => 'M',
					     'GTA' => 'V',
					     'TAG' => 'Q',
					     'CTT' => 'L',
					     'GGA' => 'G',
					     'GTC' => 'V',
					     'TGC' => 'C',
					     'TCA' => 'S',
					     'ATT' => 'I',
					     'TAT' => 'Y',
					     'AAT' => 'N',
					     'ACT' => 'T',
					     'CAA' => 'Q',
					     'GAC' => 'D',
					     'GGT' => 'G',
					     'TCC' => 'S',
					     'TTT' => 'F',
					     'AGG' => 'R',
					     'CGT' => 'R',
					     'CGG' => 'R',
					     'CAT' => 'H',
					     'ATA' => 'I',
					     'CCC' => 'P',
					     'GGG' => 'G',
					     'TTA' => 'L',
					     'GAG' => 'E',
					     'CTA' => 'L',
					     'GAT' => 'D',
					     'TCT' => 'S',
					     'TTC' => 'F',
					     'GCG' => 'A',
					     'GGC' => 'G',
					     'GCT' => 'A',
					     'GAA' => 'E',
					     'CCT' => 'P'
					    }
				 },
	     'transl_table=16' => {
				   'name' => 'The Chlorophycean Mitochondrial Code (transl_table=16)',
				   'code' => {
					      'GCC' => 'A',
					      'AGT' => 'S',
					      'TGA' => '*',
					      'TGT' => 'C',
					      'CGA' => 'R',
					      'ATC' => 'I',
					      'AAC' => 'N',
					      'AGC' => 'S',
					      'TAC' => 'Y',
					      'ACA' => 'T',
					      'TCG' => 'S',
					      'CCG' => 'P',
					      'CTG' => 'L',
					      'GCA' => 'A',
					      'GTG' => 'V',
					      'AAG' => 'K',
					      'GTT' => 'V',
					      'CAC' => 'H',
					      'AGA' => 'R',
					      'ACC' => 'T',
					      'CCA' => 'P',
					      'TGG' => 'W',
					      'CTC' => 'L',
					      'CGC' => 'R',
					      'TTG' => 'L',
					      'TAA' => '*',
					      'CAG' => 'Q',
					      'ACG' => 'T',
					      'AAA' => 'K',
					      'ATG' => 'M',
					      'GTA' => 'V',
					      'TAG' => 'L',
					      'CTT' => 'L',
					      'GGA' => 'G',
					      'GTC' => 'V',
					      'TGC' => 'C',
					      'TCA' => 'S',
					      'ATT' => 'I',
					      'TAT' => 'Y',
					      'AAT' => 'N',
					      'ACT' => 'T',
					      'CAA' => 'Q',
					      'GAC' => 'D',
					      'GGT' => 'G',
					      'TCC' => 'S',
					      'TTT' => 'F',
					      'AGG' => 'R',
					      'CGT' => 'R',
					      'CGG' => 'R',
					      'CAT' => 'H',
					      'ATA' => 'I',
					      'CCC' => 'P',
					      'GGG' => 'G',
					      'TTA' => 'L',
					      'GAG' => 'E',
					      'CTA' => 'L',
					      'GAT' => 'D',
					      'TCT' => 'S',
					      'TTC' => 'F',
					      'GCG' => 'A',
					      'GGC' => 'G',
					      'GCT' => 'A',
					      'GAA' => 'E',
					      'CCT' => 'P'
					     }
				  },
	     'transl_table=5' => {
				  'name' => 'The Invertebrate Mitochondrial Code (transl_table=5) ',
				  'code' => {
					     'GCC' => 'A',
					     'AGT' => 'S',
					     'TGA' => 'W',
					     'TGT' => 'C',
					     'CGA' => 'R',
					     'ATC' => 'I',
					     'AAC' => 'N',
					     'AGC' => 'S',
					     'TAC' => 'Y',
					     'ACA' => 'T',
					     'TCG' => 'S',
					     'CCG' => 'P',
					     'CTG' => 'L',
					     'GCA' => 'A',
					     'GTG' => 'V',
					     'AAG' => 'K',
					     'GTT' => 'V',
					     'CAC' => 'H',
					     'AGA' => 'S',
					     'ACC' => 'T',
					     'CCA' => 'P',
					     'TGG' => 'W',
					     'CTC' => 'L',
					     'CGC' => 'R',
					     'TTG' => 'L',
					     'TAA' => '*',
					     'CAG' => 'Q',
					     'ACG' => 'T',
					     'AAA' => 'K',
					     'ATG' => 'M',
					     'GTA' => 'V',
					     'TAG' => '*',
					     'CTT' => 'L',
					     'GGA' => 'G',
					     'GTC' => 'V',
					     'TGC' => 'C',
					     'TCA' => 'S',
					     'ATT' => 'I',
					     'TAT' => 'Y',
					     'AAT' => 'N',
					     'ACT' => 'T',
					     'CAA' => 'Q',
					     'GAC' => 'D',
					     'GGT' => 'G',
					     'TCC' => 'S',
					     'TTT' => 'F',
					     'AGG' => 'S',
					     'CGT' => 'R',
					     'CGG' => 'R',
					     'CAT' => 'H',
					     'ATA' => 'M',
					     'CCC' => 'P',
					     'GGG' => 'G',
					     'TTA' => 'L',
					     'GAG' => 'E',
					     'CTA' => 'L',
					     'GAT' => 'D',
					     'TCT' => 'S',
					     'TTC' => 'F',
					     'GCG' => 'A',
					     'GGC' => 'G',
					     'GCT' => 'A',
					     'GAA' => 'E',
					     'CCT' => 'P'
					    }
				 },
	     'transl_table=14' => {
				   'name' => 'The Alternative Flatworm Mitochondrial Code (transl_table=14)',
				   'code' => {
					      'GCC' => 'A',
					      'AGT' => 'S',
					      'TGA' => 'W',
					      'TGT' => 'C',
					      'CGA' => 'R',
					      'ATC' => 'I',
					      'AAC' => 'N',
					      'AGC' => 'S',
					      'TAC' => 'Y',
					      'ACA' => 'T',
					      'TCG' => 'S',
					      'CCG' => 'P',
					      'CTG' => 'L',
					      'GCA' => 'A',
					      'GTG' => 'V',
					      'AAG' => 'K',
					      'GTT' => 'V',
					      'CAC' => 'H',
					      'AGA' => 'S',
					      'ACC' => 'T',
					      'CCA' => 'P',
					      'TGG' => 'W',
					      'CTC' => 'L',
					      'CGC' => 'R',
					      'TTG' => 'L',
					      'TAA' => 'Y',
					      'CAG' => 'Q',
					      'ACG' => 'T',
					      'AAA' => 'N',
					      'ATG' => 'M',
					      'GTA' => 'V',
					      'TAG' => '*',
					      'CTT' => 'L',
					      'GGA' => 'G',
					      'GTC' => 'V',
					      'TGC' => 'C',
					      'TCA' => 'S',
					      'ATT' => 'I',
					      'TAT' => 'Y',
					      'AAT' => 'N',
					      'ACT' => 'T',
					      'CAA' => 'Q',
					      'GAC' => 'D',
					      'GGT' => 'G',
					      'TCC' => 'S',
					      'TTT' => 'F',
					      'AGG' => 'S',
					      'CGT' => 'R',
					      'CGG' => 'R',
					      'CAT' => 'H',
					      'ATA' => 'I',
					      'CCC' => 'P',
					      'GGG' => 'G',
					      'TTA' => 'L',
					      'GAG' => 'E',
					      'CTA' => 'L',
					      'GAT' => 'D',
					      'TCT' => 'S',
					      'TTC' => 'F',
					      'GCG' => 'A',
					      'GGC' => 'G',
					      'GCT' => 'A',
					      'GAA' => 'E',
					      'CCT' => 'P'
					     }
				  },
	     'transl_table=21' => {
				   'name' => 'The Trematode Mitochondrial Code (transl_table=21)',
				   'code' => {
					      'GCC' => 'A',
					      'AGT' => 'S',
					      'TGA' => 'W',
					      'TGT' => 'C',
					      'CGA' => 'R',
					      'ATC' => 'I',
					      'AAC' => 'N',
					      'AGC' => 'S',
					      'TAC' => 'Y',
					      'ACA' => 'T',
					      'TCG' => 'S',
					      'CCG' => 'P',
					      'CTG' => 'L',
					      'GCA' => 'A',
					      'GTG' => 'V',
					      'AAG' => 'K',
					      'GTT' => 'V',
					      'CAC' => 'H',
					      'AGA' => 'S',
					      'ACC' => 'T',
					      'CCA' => 'P',
					      'TGG' => 'W',
					      'CTC' => 'L',
					      'CGC' => 'R',
					      'TTG' => 'L',
					      'TAA' => '*',
					      'CAG' => 'Q',
					      'ACG' => 'T',
					      'AAA' => 'N',
					      'ATG' => 'M',
					      'GTA' => 'V',
					      'TAG' => '*',
					      'CTT' => 'L',
					      'GGA' => 'G',
					      'GTC' => 'V',
					      'TGC' => 'C',
					      'TCA' => 'S',
					      'ATT' => 'I',
					      'TAT' => 'Y',
					      'AAT' => 'N',
					      'ACT' => 'T',
					      'CAA' => 'Q',
					      'GAC' => 'D',
					      'GGT' => 'G',
					      'TCC' => 'S',
					      'TTT' => 'F',
					      'AGG' => 'S',
					      'CGT' => 'R',
					      'CGG' => 'R',
					      'CAT' => 'H',
					      'ATA' => 'M',
					      'CCC' => 'P',
					      'GGG' => 'G',
					      'TTA' => 'L',
					      'GAG' => 'E',
					      'CTA' => 'L',
					      'GAT' => 'D',
					      'TCT' => 'S',
					      'TTC' => 'F',
					      'GCG' => 'A',
					      'GGC' => 'G',
					      'GCT' => 'A',
					      'GAA' => 'E',
					      'CCT' => 'P'
					     }
				  },
	     'transl_table=11' => {
				   'name' => 'The Bacterial and Plant Plastid Code (transl_table=11) ',
				   'code' => {
					      'GCC' => 'A',
					      'AGT' => 'S',
					      'TGA' => '*',
					      'TGT' => 'C',
					      'CGA' => 'R',
					      'ATC' => 'I',
					      'AAC' => 'N',
					      'AGC' => 'S',
					      'TAC' => 'Y',
					      'ACA' => 'T',
					      'TCG' => 'S',
					      'CCG' => 'P',
					      'CTG' => 'L',
					      'GCA' => 'A',
					      'GTG' => 'V',
					      'AAG' => 'K',
					      'GTT' => 'V',
					      'CAC' => 'H',
					      'AGA' => 'R',
					      'ACC' => 'T',
					      'CCA' => 'P',
					      'TGG' => 'W',
					      'CTC' => 'L',
					      'CGC' => 'R',
					      'TTG' => 'L',
					      'TAA' => '*',
					      'CAG' => 'Q',
					      'ACG' => 'T',
					      'AAA' => 'K',
					      'ATG' => 'M',
					      'GTA' => 'V',
					      'TAG' => '*',
					      'CTT' => 'L',
					      'GGA' => 'G',
					      'GTC' => 'V',
					      'TGC' => 'C',
					      'TCA' => 'S',
					      'ATT' => 'I',
					      'TAT' => 'Y',
					      'AAT' => 'N',
					      'ACT' => 'T',
					      'CAA' => 'Q',
					      'GAC' => 'D',
					      'GGT' => 'G',
					      'TCC' => 'S',
					      'TTT' => 'F',
					      'AGG' => 'R',
					      'CGT' => 'R',
					      'CGG' => 'R',
					      'CAT' => 'H',
					      'ATA' => 'I',
					      'CCC' => 'P',
					      'GGG' => 'G',
					      'TTA' => 'L',
					      'GAG' => 'E',
					      'CTA' => 'L',
					      'GAT' => 'D',
					      'TCT' => 'S',
					      'TTC' => 'F',
					      'GCG' => 'A',
					      'GGC' => 'G',
					      'GCT' => 'A',
					      'GAA' => 'E',
					      'CCT' => 'P'
					     }
				  },
	     'transl_table=4' => {
				  'name' => 'The Mycoplasma/Spiroplasma Code (transl_table=4)',
				  'code' => {
					     'GCC' => 'A',
					     'AGT' => 'S',
					     'TGA' => 'W',
					     'TGT' => 'C',
					     'CGA' => 'R',
					     'ATC' => 'I',
					     'AAC' => 'N',
					     'AGC' => 'S',
					     'TAC' => 'Y',
					     'ACA' => 'T',
					     'TCG' => 'S',
					     'CCG' => 'P',
					     'CTG' => 'L',
					     'GCA' => 'A',
					     'GTG' => 'V',
					     'AAG' => 'K',
					     'GTT' => 'V',
					     'CAC' => 'H',
					     'AGA' => 'R',
					     'ACC' => 'T',
					     'CCA' => 'P',
					     'TGG' => 'W',
					     'CTC' => 'L',
					     'CGC' => 'R',
					     'TTG' => 'L',
					     'TAA' => '*',
					     'CAG' => 'Q',
					     'ACG' => 'T',
					     'AAA' => 'K',
					     'ATG' => 'M',
					     'GTA' => 'V',
					     'TAG' => '*',
					     'CTT' => 'L',
					     'GGA' => 'G',
					     'GTC' => 'V',
					     'TGC' => 'C',
					     'TCA' => 'S',
					     'ATT' => 'I',
					     'TAT' => 'Y',
					     'AAT' => 'N',
					     'ACT' => 'T',
					     'CAA' => 'Q',
					     'GAC' => 'D',
					     'GGT' => 'G',
					     'TCC' => 'S',
					     'TTT' => 'F',
					     'AGG' => 'R',
					     'CGT' => 'R',
					     'CGG' => 'R',
					     'CAT' => 'H',
					     'ATA' => 'I',
					     'CCC' => 'P',
					     'GGG' => 'G',
					     'TTA' => 'L',
					     'GAG' => 'E',
					     'CTA' => 'L',
					     'GAT' => 'D',
					     'TCT' => 'S',
					     'TTC' => 'F',
					     'GCG' => 'A',
					     'GGC' => 'G',
					     'GCT' => 'A',
					     'GAA' => 'E',
					     'CCT' => 'P'
					    }
				 },
	     'transl_table=1' => {
				  'name' => 'The Standard Code (transl_table=1)',
				  'code' => {
					     'GCC' => 'A',
					     'AGT' => 'S',
					     'TGA' => '*',
					     'TGT' => 'C',
					     'CGA' => 'R',
					     'ATC' => 'I',
					     'AAC' => 'N',
					     'AGC' => 'S',
					     'TAC' => 'Y',
					     'ACA' => 'T',
					     'TCG' => 'S',
					     'CCG' => 'P',
					     'CTG' => 'L',
					     'GCA' => 'A',
					     'GTG' => 'V',
					     'AAG' => 'K',
					     'GTT' => 'V',
					     'CAC' => 'H',
					     'AGA' => 'R',
					     'ACC' => 'T',
					     'CCA' => 'P',
					     'TGG' => 'W',
					     'CTC' => 'L',
					     'CGC' => 'R',
					     'TTG' => 'L',
					     'TAA' => '*',
					     'CAG' => 'Q',
					     'ACG' => 'T',
					     'AAA' => 'K',
					     'ATG' => 'M',
					     'GTA' => 'V',
					     'TAG' => '*',
					     'CTT' => 'L',
					     'GGA' => 'G',
					     'GTC' => 'V',
					     'TGC' => 'C',
					     'TCA' => 'S',
					     'ATT' => 'I',
					     'TAT' => 'Y',
					     'AAT' => 'N',
					     'ACT' => 'T',
					     'CAA' => 'Q',
					     'GAC' => 'D',
					     'GGT' => 'G',
					     'TCC' => 'S',
					     'TTT' => 'F',
					     'AGG' => 'R',
					     'CGT' => 'R',
					     'CGG' => 'R',
					     'CAT' => 'H',
					     'ATA' => 'I',
					     'CCC' => 'P',
					     'GGG' => 'G',
					     'TTA' => 'L',
					     'GAG' => 'E',
					     'CTA' => 'L',
					     'GAT' => 'D',
					     'TCT' => 'S',
					     'TTC' => 'F',
					     'GCG' => 'A',
					     'GGC' => 'G',
					     'GCT' => 'A',
					     'GAA' => 'E',
					     'CCT' => 'P'
					    }
				 },
	     'transl_table=3' => {
				  'name' => 'The Yeast Mitochondrial Code (transl_table=3)',
				  'code' => {
					     'GCC' => 'A',
					     'AGT' => 'S',
					     'TGA' => 'W',
					     'TGT' => 'C',
					     'CGA' => 'R',
					     'ATC' => 'I',
					     'AAC' => 'N',
					     'AGC' => 'S',
					     'TAC' => 'Y',
					     'ACA' => 'T',
					     'TCG' => 'S',
					     'CCG' => 'P',
					     'CTG' => 'T',
					     'GCA' => 'A',
					     'GTG' => 'V',
					     'AAG' => 'K',
					     'GTT' => 'V',
					     'CAC' => 'H',
					     'AGA' => 'R',
					     'ACC' => 'T',
					     'CCA' => 'P',
					     'TGG' => 'W',
					     'CTC' => 'T',
					     'CGC' => 'R',
					     'TTG' => 'L',
					     'TAA' => '*',
					     'CAG' => 'Q',
					     'ACG' => 'T',
					     'AAA' => 'K',
					     'ATG' => 'M',
					     'GTA' => 'V',
					     'TAG' => '*',
					     'CTT' => 'T',
					     'GGA' => 'G',
					     'GTC' => 'V',
					     'TGC' => 'C',
					     'TCA' => 'S',
					     'ATT' => 'I',
					     'TAT' => 'Y',
					     'AAT' => 'N',
					     'ACT' => 'T',
					     'CAA' => 'Q',
					     'GAC' => 'D',
					     'GGT' => 'G',
					     'TCC' => 'S',
					     'TTT' => 'F',
					     'AGG' => 'R',
					     'CGT' => 'R',
					     'CGG' => 'R',
					     'CAT' => 'H',
					     'ATA' => 'M',
					     'CCC' => 'P',
					     'GGG' => 'G',
					     'TTA' => 'L',
					     'GAG' => 'E',
					     'CTA' => 'T',
					     'GAT' => 'D',
					     'TCT' => 'S',
					     'TTC' => 'F',
					     'GCG' => 'A',
					     'GGC' => 'G',
					     'GCT' => 'A',
					     'GAA' => 'E',
					     'CCT' => 'P'
					    }
				 },
	     'transl_table=10' => {
				   'name' => 'The Euplotid Nuclear Code (transl_table=10)',
				   'code' => {
					      'GCC' => 'A',
					      'AGT' => 'S',
					      'TGA' => 'C',
					      'TGT' => 'C',
					      'CGA' => 'R',
					      'ATC' => 'I',
					      'AAC' => 'N',
					      'AGC' => 'S',
					      'TAC' => 'Y',
					      'ACA' => 'T',
					      'TCG' => 'S',
					      'CCG' => 'P',
					      'CTG' => 'L',
					      'GCA' => 'A',
					      'GTG' => 'V',
					      'AAG' => 'K',
					      'GTT' => 'V',
					      'CAC' => 'H',
					      'AGA' => 'R',
					      'ACC' => 'T',
					      'CCA' => 'P',
					      'TGG' => 'W',
					      'CTC' => 'L',
					      'CGC' => 'R',
					      'TTG' => 'L',
					      'TAA' => '*',
					      'CAG' => 'Q',
					      'ACG' => 'T',
					      'AAA' => 'K',
					      'ATG' => 'M',
					      'GTA' => 'V',
					      'TAG' => '*',
					      'CTT' => 'L',
					      'GGA' => 'G',
					      'GTC' => 'V',
					      'TGC' => 'C',
					      'TCA' => 'S',
					      'ATT' => 'I',
					      'TAT' => 'Y',
					      'AAT' => 'N',
					      'ACT' => 'T',
					      'CAA' => 'Q',
					      'GAC' => 'D',
					      'GGT' => 'G',
					      'TCC' => 'S',
					      'TTT' => 'F',
					      'AGG' => 'R',
					      'CGT' => 'R',
					      'CGG' => 'R',
					      'CAT' => 'H',
					      'ATA' => 'I',
					      'CCC' => 'P',
					      'GGG' => 'G',
					      'TTA' => 'L',
					      'GAG' => 'E',
					      'CTA' => 'L',
					      'GAT' => 'D',
					      'TCT' => 'S',
					      'TTC' => 'F',
					      'GCG' => 'A',
					      'GGC' => 'G',
					      'GCT' => 'A',
					      'GAA' => 'E',
					      'CCT' => 'P'
					     }
				  }
	    };
    $aa_info = {
		"*"=>{
		     name=>"Stop",
		     "3-letter"=>"Stop",
		     "1-letter"=>"*",
		     energy=>"*",
		     polarity=>"*",
		     "charge_ph7.4"=>"*",
		     "hydropathy"=>"*"
		    },
		"A"=>{
		     name=>"Alanine",
		     "3-letter"=>"Ala",
		     "1-letter"=>"A",
		     energy=>"11.7",
		     polarity=>"nonpolar",
		     "charge_ph7.4"=>"neutral",
		     "hydropathy"=>"1.8"
		    },
		"R"=>{
		     name=>"Arginine",
		     "3-letter"=>"Arg",
		     "1-letter"=>"R",
		     energy=>"27.3",
		     polarity=>"polar",
		     "charge_ph7.4"=>"positive",
		     "hydropathy"=>"-4.5"
		    },
		"N"=>{
		     name=>"Asparagine",
		     "3-letter"=>"Asn",
		     "1-letter"=>"N",
		     energy=>"14.7",
		     polarity=>"polar",
		     "charge_ph7.4"=>"neutral",
		     "hydropathy"=>"-3.5"
		    },
		"D"=>{
		     name=>"Aspartate",
		     "3-letter"=>"Asp",
		     "1-letter"=>"D",
		     energy=>"12.7",
		     polarity=>"polar",
		     "charge_ph7.4"=>"negative",
		     "hydropathy"=>"-3.5"
		    },
		"C"=>{
		     name=>"Cysteine",
		     "3-letter"=>"Cys",
		     "1-letter"=>"C",
		     energy=>"24.7",
		     polarity=>"nonpolar",
		     "charge_ph7.4"=>"neutral",
		     "hydropathy"=>"2.5"
		    },
		"E"=>{
		     name=>"Glutamate",
		     "3-letter"=>"Glu",
		     "1-letter"=>"E",
		     energy=>"15.3",
		     polarity=>"polar",
		     "charge_ph7.4"=>"negative",
		     "hydropathy"=>"-3.5"
		    },
		"Q"=>{
		     name=>"Glutamine",
		     "3-letter"=>"Gln",
		     "1-letter"=>"Q",
		     energy=>"16.3",
		     polarity=>"polar",
		     "charge_ph7.4"=>"neutral",
		     "hydropathy"=>"-3.5"
		    },
		"G"=>{
		     name=>"Glycine",
		     "3-letter"=>"Gly",
		     "1-letter"=>"G",
		     energy=>"11.7",
		     polarity=>"nonpolar",
		     "charge_ph7.4"=>"neutral",
		     "hydropathy"=>"-0.4"
		    },
		"H"=>{
		     name=>"Histidine",
		     "3-letter"=>"His",
		     "1-letter"=>"H",
		     energy=>"38.3",
		     polarity=>"polar",
		     "charge_ph7.4"=>'10% pos/90% neut',
		     "hydropathy"=>"-3.2"
		    },
		"I"=>{
		     name=>"Isoleucine",
		     "3-letter"=>"Ile",
		     "1-letter"=>"I",
		     energy=>"32.3",
		     polarity=>"nonpolar",
		     "charge_ph7.4"=>"neutral",
		     "hydropathy"=>"4.5"
		    },
		"L"=>{
		     name=>"Leucine",
		     "3-letter"=>"Leu",
		     "1-letter"=>"L",
		     energy=>"27.3",
		     polarity=>"nonpolar",
		     "charge_ph7.4"=>"neutral",
		     "hydropathy"=>"3.8"
		    },
		"K"=>{
		     name=>"Lysine",
		     "3-letter"=>"Lys",
		     "1-letter"=>"K",
		     energy=>"30.3",
		     polarity=>"polar",
		     "charge_ph7.4"=>"positive",
		     "hydropathy"=>"-3.9"
		    },
		"M"=>{
		     name=>"Methionine",
		     "3-letter"=>"Met",
		     "1-letter"=>"M",
		     energy=>"34.3",
		     polarity=>"nonpolar",
		     "charge_ph7.4"=>"neutral",
		     "hydropathy"=>"1.9"
		    },
		"F"=>{
		     name=>"Phenylalanine",
		     "3-letter"=>"Phe",
		     "1-letter"=>"F",
		     energy=>"52.0",
		     polarity=>"nonpolar",
		     "charge_ph7.4"=>"neutral",
		     "hydropathy"=>"2.8"
		    },
		"P"=>{
		     name=>"Proline",
		     "3-letter"=>"Pro",
		     "1-letter"=>"P",
		     energy=>"20.3",
		     polarity=>"nonpolar",
		     "charge_ph7.4"=>"neutral",
		     "hydropathy"=>"-1.6"
		    },
		"S"=>{
		     name=>"Serine",
		     "3-letter"=>"Ser",
		     "1-letter"=>"S",
		     energy=>"11.7",
		     polarity=>"polar",
		     "charge_ph7.4"=>"neutral",
		     "hydropathy"=>"-0.8"
		    },
		"T"=>{
		     name=>"Threonine",
		     "3-letter"=>"Thr",
		     "1-letter"=>"T",
		     energy=>"18.7",
		     polarity=>"polar",
		     "charge_ph7.4"=>"neutral",
		     "hydropathy"=>"-0.7"
		    },
		"W"=>{
		     name=>"Tryptophan",
		     "3-letter"=>"Trp",
		     "1-letter"=>"W",
		     energy=>"74.3",
		     polarity=>"nonpolar",
		     "charge_ph7.4"=>"neutral",
		     "hydropathy"=>"-0.9"
		    },
		"Y"=>{
		     name=>"Tyrosine",
		     "3-letter"=>"Tyr",
		     "1-letter"=>"Y",
		     energy=>"50.0",
		     polarity=>"polar",
		     "charge_ph7.4"=>"neutral",
		     "hydropathy"=>"-1.3"
		    },
		"V"=>{
		     name=>"Valine",
		     "3-letter"=>"Val",
		     "1-letter"=>"V",
		     energy=>"23.3",
		     polarity=>"nonpolar",
		     "charge_ph7.4"=>"neutral",
		     "hydropathy"=>"4.2"
		    },
	       }
  }

sub code
  {
    my $self = shift if ref($_[0]) =~ /genetic_code/ ;
    my $type = shift || "transl_table=1";
    $type = "transl_table=".$type unless $type =~ /transl_table=/;
    $type =~ s/\s//g;
    return $code->{$type} if $code->{$type};
    return $code->{"transl_table=1"};

  }

sub all_codes
  {
    my $self = shift;
    return $code;
  }

sub html_code_table
    {
      my $self = shift;
      my %opts = @_;
      my $data = $opts{data};
      my $trans_table = $opts{trans_table};
      my $code = $opts{code};
      my $counts = $opts{counts};
      my $two_colors = $opts{two_colors};
      my $total_codon_count =0;
      if ($counts)
	{
	  $counts = $data;
	  grep {$total_codon_count+=$_} values %$counts;
	  $data = {map {$_,$counts->{$_}/$total_codon_count} keys %$counts};
	}
      $code = $self->code($trans_table) unless $code;
      $code = $code->{code} if $code->{code};
      my ($max_val) = sort {$b <=> $a} map {$data->{$_}} keys %{$code};
      my ($min_val) = sort {$a <=> $b} map {$data->{$_}} keys %{$code};
      my $range = $max_val-$min_val;
      my $html;
      $html .= "<table class='ui-widget-content ui-corner-all xsmall'><tr><td nowrap>";
      my $count = 0;
      my $codon_set_total =0;
      foreach my $codon (sort { $self->sort_nt1(substr($a, 0, 1)) <=> $self->sort_nt1(substr($b,0, 1)) || $self->sort_nt2(substr($a,1,1)) <=> $self->sort_nt2(substr($b,1,1)) || $self->sort_nt3(substr($a,2,1)) <=> $self->sort_nt3(substr($b,2,1)) } keys %{$code})
	{
	  my $color_str;
	  my $current_val = $data->{$codon} =~ /\d/ ? sprintf("%.2f",100*$data->{$codon}) : $data->{$codon};
	  if ($current_val =~ /^\d/)
	    {
	      if ($two_colors)
		{
		  #fade to white at 100
		  #color green if above 100
		  if ($current_val > 100)
		    {
		      my $color = sprintf("%.0f", 255-255*($current_val-100)/($max_val*100-100));
		      $color_str = "rgb($color, 255, $color);";
		    }
		  #color red if below 100;
		  else
		    {
		      my $color = sprintf("%.0f", 255*$current_val/100);
		      $color_str = "rgb(255, $color, $color);";
		    }
		}
	      else
		{
		  my $color = $self->color_by_usage(100*$max_val, $current_val);
		  $color_str = "rgb($color, 255, $color);";
		}
	    }
#	  my $rel_val = ($data->{$codon}-$min_val)/$range;
#	  my $color = $self->get_color(val=>$rel_val);
#	  my $color_str = join (",",@$color);
	  my $str = "<span style=\"background-color: $color_str\" >".$codon."(".$code->{$codon}.") ".$current_val."%";
	  $str .= " (".$counts->{$codon}.")" if $counts;
	  $str .="</span>";
	  unless ($count % 4)
	    {
	      $html .= "<span class=small>Total: ".100*sprintf("%.4f",$codon_set_total/$total_codon_count)."% ($codon_set_total)</span>" if $counts && $count && $count!=16;
	      $html .= "<td size=1 style=\"background-color: rgb(200,200,200)\"><td nowrap>" if $count && $count != 16;
	      $codon_set_total=0 unless $count == 16;
	    }
	  if ($count == 16)
	    {
	      $count = 0;
	      $html .= "<span class=small>Total: ".100*sprintf("%.4f",$codon_set_total/$total_codon_count)."% ($codon_set_total)</span>" if $counts;
	      $html .= "<tr size=1><td colspan=7 style=\"background-color: rgb(200,200,200)\"><tr><td nowrap>";
	      $codon_set_total=0;
	    }
	  $html .= $str."<br>";
	  $count++;
	  $codon_set_total+=$counts->{$codon} if $counts;
#
	}
      $html .= "<span class=small>Total: ".100*sprintf("%.4f",$codon_set_total/$total_codon_count)."% ($codon_set_total)</span>" if $counts;
      $html .= "</table>";
      return $html;
    }
sub sort_nt1
  {
    my $self = shift;
    my $chr = uc(shift);

    $chr = substr($chr, -1,1) if length($chr)>1;
    my $val = 0;
    if ($chr eq "G")
      {
	$val = 1;
      }
    elsif ($chr eq "A")
      {
	$val = 2;
      }
    elsif ($chr eq "T")
      {
	$val = 3;
      }
    return $val;
  }

sub sort_nt2
  {
    my $self = shift;
    my $chr = uc(shift);

    $chr = substr($chr, -1,1) if length($chr)>1;
    my $val = 0;
    if ($chr eq "G")
      {
	$val = 1;
      }
    elsif ($chr eq "A")
      {
	$val = 2;
      }
    elsif ($chr eq "T")
      {
	$val = 3;
      }
    return $val;
  }

sub sort_nt3
  {
    my $self = shift;
    my $chr = uc(shift);

    $chr = substr($chr, -1,1) if length($chr)>1;
    my $val = 0;
    if ($chr eq "G")
      {
	$val = 1;
      }
    elsif ($chr eq "T")
      {
	$val = 2;
      }
    elsif ($chr eq "C")
      {
	$val = 3;
      }
    return $val;
  }

sub color_by_usage
  {
    my $self = shift;
    my ($max,$value, $opt) = @_;
    $opt = 255 unless $opt;
    return $opt unless $max;
    my $g = $opt*(($max - $value) / $max);
    return int($g + .5);
  }

sub get_color
  {
    my $self = shift;
    my %opts = @_;
    my $val = $opts{val};
    return [0,0,0] unless defined $val;
    return [125,125,125] if $val < 0;
    my @colors = (
                  [255,0,0], #red
                  [255,126,0], #orange
                  [255,255,0], #yellow
                  [0,255,0], # green
                  [0,255,255], # cyan
                  [0,0,255], # blue
#                 [255,0,255], #magenta
                  [126,0,126], #purple
                 );
    @colors = reverse @colors;
    my ($index1, $index2) = ((floor((scalar(@colors)-1)*$val)), ceil((scalar(@colors)-1)*$val));

    my $color=[];
    my $step = 1/(scalar (@colors)-1);
    my $scale = $index1*$step;
    my $dist = ($val-$scale)/($step);
    for (my $i=0; $i<=2; $i++)
      {
        my $diff = ($colors[$index1][$i]-$colors[$index2][$i])*$dist;
        push @$color, sprintf("%.0f", $colors[$index1][$i]-$diff);
      }
    return $color;
  }

sub sort_aa_by_gc
  {
    my $self = shift;
    my %opts = @_;
    my $trans_table = $opts{trans_table};
    my $code = $opts{code};
    $code = $self->code($trans_table) unless $code;
    my %aa_sort = map {$_,{}} values %{$code->{code}};
    foreach my $codon (keys %{$code->{code}})
      {
	my $gc = $codon =~ tr/GCgc/GCgc/;
	my $at = $codon =~ tr/ATat/ATat/;
	$aa_sort{$code->{code}{$codon}}{AT}+=$at;
	$aa_sort{$code->{code}{$codon}}{GC}+=$gc;
      }
    %aa_sort = map {$_,($aa_sort{$_}{GC}/($aa_sort{$_}{AT}+$aa_sort{$_}{GC}))} keys %aa_sort;
    return \%aa_sort;
  }

sub html_aa_new
  {
    my $self = shift;
    my %opts = @_;
    my $trans_table = $opts{trans_table};
    my $code = $opts{code};
    $code = $self->code($trans_table) unless $code;
    my $data = $opts{data};
    my $counts = $opts{counts};
    my $two_colors = $opts{two_colors};
    my $table_name = $opts{table_name};
    $table_name = "aa_table" unless $table_name;
    if ($counts)
      {
	$counts = $data;
	my $total = 0;
	grep {$total+=$_} values %$counts;
	$data = {map {$_,$counts->{$_}/$total} keys %$counts};
      }
    my $aa_sort = $self->sort_aa_by_gc(code=>$code, trans_table=>$trans_table);
    my ($max_val) = sort {$b<=>$a} map{$data->{$_}} keys %$aa_sort;
    my ($min_val) = sort {$a<=>$b} map{$data->{$_}} keys %$aa_sort;
    my $range = $max_val-$min_val;
    my $html;
    $html .= qq{<table id="$table_name" class='ui-widget-content ui-corner-all xsmall'><Thead><tr>
 <th>AA</th>
 <th>Polarity</th>
 <th>Charge pH7.4</th>
 <th>Hydropathy</th>
 <th>\%GC</th>
 <th>ATP Costs</th>
 <TH>\%Usage</th>
 </tr></THEAD>
<tbody align=left valign="top" id="aa_table_body">
};
    my $total = 0;
    foreach (sort {$aa_sort->{$b} <=> $aa_sort->{$a} || $a cmp $b}keys %$aa_sort)
      {
	my $color_str;
	my $current_val = $data->{$_} =~ /\d/ ? sprintf("%.2f",100*$data->{$_}) : $data->{$_};#sprintf("%.2f",100*$data->{$_});
	if ($current_val =~ /^\d/)
	    {
	      if ($two_colors)
		{
		  #fade to white at 100
		  #color green if above 100
		  if ($current_val > 100)
		    {
		      my $color = sprintf("%.0f", 255-255*($current_val-100)/($max_val*100-100));
		      $color_str = "rgb($color, 255, $color);";
		    }
		  #color red if below 100;
		  else
		    {
		      my $color = sprintf("%.0f", 255*$current_val/100);
		      $color_str = "rgb(255, $color, $color);";
		    }
		}
	      else
		{
		  my $color = $self->color_by_usage(100*$max_val, $current_val);
		  $color_str = "rgb($color, 255, $color);";
		}
	    }
	#	my $rel_val = ($data->{$_}-$min_aa)/$range;
#	my $color = $self->get_color(val=>$rel_val);
#	my $color_str = join (",",@$color);
#	my $color = $self->color_by_usage(100*$max_aa,$current_val);
	$html .= qq{<tr style="background-color: $color_str">}.
 "<td nowrap>$_". " (".$aa_info->{$_}{"3-letter"}.")".
  qq{<td nowrap>}.$aa_info->{$_}{"polarity"}.qq{</td>}.
  qq{<td nowrap>}.$aa_info->{$_}{"charge_ph7.4"}.qq{</td>}.
  qq{<td nowrap>}.$aa_info->{$_}{"hydropathy"}.qq{</td>}.
  "<td nowrap>".sprintf("%.0f",100*$aa_sort->{$_})."%".qq{</td>}.
  qq{<td nowrap>}.$aa_info->{$_}{"energy"}.qq{</td>}.
  "<td nowrap>".$current_val."%";
	$html .= " (".$counts->{$_}.")" if $counts;
        $html .= qq{</td>};
	$total++;
        $html .= qq{</tr>};
      }
    $html .= "</tbody></table>";
    return $html;
  }

sub html_aa
  {
    my $self = shift;
    my %opts = @_;
    my $trans_table = $opts{trans_table};
    my $code = $opts{code};
    my $split_table =$opts{split};
    $code = $self->code($trans_table) unless $code;    my $data = $opts{data};
    my $counts = $opts{counts};
    my $two_colors = $opts{two_colors};
    if ($counts)
      {
	$counts = $data;
	my $total = 0;
	grep {$total+=$_} values %$counts;
	$data = {map {$_,$counts->{$_}/$total} keys %$counts};
      }
    my $aa_sort = $self->sort_aa_by_gc(code=>$code, trans_table=>$trans_table);
    my ($max_val) = sort {$b<=>$a} map{$data->{$_}} keys %$aa_sort;
    my ($min_val) = sort {$a<=>$b} map{$data->{$_}} keys %$aa_sort;
    my $range = $max_val-$min_val;
    my $html;
    $html .= "<table class='ui-widget-content ui-corner-all'><tr valign=top><td>" if $split_table;
    $html .= "<table class='ui-widget-content ui-corner-all'>";
    my $total = 0;
    foreach (sort {$aa_sort->{$b} <=> $aa_sort->{$a} || $a cmp $b}keys %$aa_sort)
      {
	my $color_str;
	my $current_val = $data->{$_} =~ /\d/ ? sprintf("%.2f",100*$data->{$_}) : $data->{$_};#sprintf("%.2f",100*$data->{$_});
	if ($current_val =~ /^\d/)
	    {
	      if ($two_colors)
		{
		  #fade to white at 100
		  #color green if above 100
		  if ($current_val > 100)
		    {
		      my $color = sprintf("%.0f", 255-255*($current_val-100)/($max_val*100-100));
		      $color_str = "rgb($color, 255, $color);";
		    }
		  #color red if below 100;
		  else
		    {
		      my $color = sprintf("%.0f", 255*$current_val/100);
		      $color_str = "rgb(255, $color, $color);";
		    }
		}
	      else
		{
		  my $color = $self->color_by_usage(100*$max_val, $current_val);
		  $color_str = "rgb($color, 255, $color);";
		}
	    }
	#	my $rel_val = ($data->{$_}-$min_aa)/$range;
#	my $color = $self->get_color(val=>$rel_val);
#	my $color_str = join (",",@$color);
#	my $color = $self->color_by_usage(100*$max_aa,$current_val);
	$html .= "<tr style=\"background-color: $color_str\"><td nowrap>$_ (GC:".sprintf("%.0f",100*$aa_sort->{$_})."%)<td nowrap>".$current_val."%";
	$html .= " (".$counts->{$_}.")" if $counts;
	$total++;
	$html .= "</table><td><table class='ui-widget-content ui-corner-all'>" if $split_table && $total == 10;
      }
    $html .= "</table>";
    $html .= "</table>" if $split_table;
    return $html;
  }

1;
