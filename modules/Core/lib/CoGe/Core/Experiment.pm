package CoGe::Core::Experiment;
use strict;

use Sort::Versions;
use CoGe::Core::Storage qw($DATA_TYPE_QUANT $DATA_TYPE_ALIGN $DATA_TYPE_POLY $DATA_TYPE_MARKER);

our ( @EXPORT, @EXPORT_OK, @ISA, $VERSION, @QUANT_TYPES, @MARKER_TYPES, 
      @OTHER_TYPES, @SUPPORTED_TYPES );

BEGIN {
    require Exporter;
    $VERSION = 0.1;
    @ISA = qw( Exporter );
    @EXPORT = qw(@QUANT_TYPES @MARKER_TYPES @OTHER_TYPES @SUPPORTED_TYPES);
    @EXPORT_OK = qw(experimentcmp detect_data_type);
    
    # Setup supported experiment file types
    @QUANT_TYPES = qw(csv tsv bed wig);
    @MARKER_TYPES = qw(gff gtf gff3);
    @OTHER_TYPES = qw(bam vcf);
    @SUPPORTED_TYPES = (@QUANT_TYPES, @MARKER_TYPES, @OTHER_TYPES);
}

sub experimentcmp($$) { # for sorting DBI-X objects or DBI hashes
    my ($a, $b) = @_;

    if ( ref($a) eq 'HASH' && ref($b) eq 'HASH' ) { # DBI
        versioncmp( $b->{version}, $a->{version} )
          || $a->{name} cmp $b->{name}
          || $b->{id} cmp $a->{id};
    }
    else { # DBI-X
        versioncmp( $b->version, $a->version )
          || $a->name cmp $b->name
          || $b->id cmp $a->id;
    }
}

sub detect_data_type {
    my $filetype = shift;
    my $filepath = shift;
    #print STDOUT "detect_data_type: $filepath\n";

    if (!$filetype or $filetype eq 'autodetect') {
        # Try to determine type based on file extension
        #print STDOUT "log: Detecting file type\n";
        ($filetype) = lc($filepath) =~ /\.([^\.]+)$/;
    }
    
    $filetype = lc($filetype);

    if ( grep { $_ eq $filetype } @QUANT_TYPES ) {
        print STDOUT "log: Detected a quantitative file ($filetype)\n";
        return ($filetype, $DATA_TYPE_QUANT);
    }
    elsif ( $filetype eq 'bam' ) {
        print STDOUT "log: Detected an alignment file ($filetype)\n";
        return ($filetype, $DATA_TYPE_ALIGN);
    }
    elsif ( $filetype eq 'vcf' ) {
        print STDOUT "log: Detected a polymorphism file ($filetype)\n";
        return ($filetype, $DATA_TYPE_POLY);
    }
    elsif ( grep { $_ eq $filetype } @MARKER_TYPES ) {
        print STDOUT "log: Detected a marker file ($filetype)\n";
        return ($filetype, $DATA_TYPE_MARKER);
    }
    else {
        print STDOUT "detect_data_type: unknown file ext '$filetype'\n";
        return ($filetype);
    }
}

1;
