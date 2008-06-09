package CoGe::Accessory::genetic_code;

use strict;
BEGIN {
    use Exporter ();
    use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $code);
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
      my $max_val = $opts{max_val};
      my $trans_table = $opts{trans_table};
      my $code = $opts{code};
      my $counts = $opts{counts};
      my $total_codon_count =0;
      if ($counts)
	{
	  $counts = $data;
	  grep {$total_codon_count+=$_} values %$counts;
	  $data = {map {$_,$counts->{$_}/$total_codon_count} keys %$counts};
	}
      $code = $self->code($trans_table) unless $code;
      $code = $code->{code} if $code->{code};
      ($max_val) = sort {$b <=> $a} map {$data->{$_}} keys %{$code};
      my $html;
      $html .= "<table><tr><td nowrap>";
      my $count = 0;
      my $codon_set_total =0;
      foreach my $codon (sort { $self->sort_nt1(substr($a, 0, 1)) <=> $self->sort_nt1(substr($b,0, 1)) || $self->sort_nt2(substr($a,1,1)) <=> $self->sort_nt2(substr($b,1,1)) || $self->sort_nt3(substr($a,2,1)) <=> $self->sort_nt3(substr($b,2,1)) } keys %{$code})
	{
	  my $tmp = substr($codon, 0, 1);
	  my $current_val = sprintf("%.2f",100*$data->{$codon});
	  my $color = $self->color_by_usage(100*$max_val, $current_val);
	  my $str = "<span style=\"background-color: rgb($color, 255, $color);\" >".$codon."(".$code->{$codon}.") ".$current_val."%";
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

sub html_aa
  {
    my $self = shift;
    my %opts = @_;
    my $trans_table = $opts{trans_table};
    my $code = $opts{code};
    my $split_table =$opts{split};
    $code = $self->code($trans_table) unless $code;    my $data = $opts{data};
    my $counts = $opts{counts};
    if ($counts)
      {
	$counts = $data;
	my $total = 0;
	grep {$total+=$_} values %$counts;
	$data = {map {$_,$counts->{$_}/$total} keys %$counts};
      }
    my $aa_sort = $self->sort_aa_by_gc(code=>$code, trans_table=>$trans_table);
    my ($max_aa) = sort {$b<=>$a} map{$data->{$_}} keys %$aa_sort;
    my $html;
    $html .= "<table><tr valign=top><td>" if $split_table;
    $html .= "<table>";
    my $total = 0;
    foreach (sort {$aa_sort->{$b} <=> $aa_sort->{$a} || $a cmp $b}keys %$aa_sort)
      {	
	my $current_val = sprintf("%.2f",100*$data->{$_});
	my $color = $self->color_by_usage(100*$max_aa,$current_val);
	$html .= "<tr style=\"background-color: rgb($color,255,$color)\"><td nowrap>$_ (GC:".sprintf("%.0f",100*$aa_sort->{$_})."%)<td nowrap>".$current_val."%";
	$html .= " (".$counts->{$_}.")" if $counts;
	$total++;
	$html .= "</table><td><table>" if $split_table && $total == 10;
      }
    $html .= "</table>";
    $html .= "</table>" if $split_table;
    return $html;
  }


1;
