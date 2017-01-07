package CoGe::Builder::Tools::SynMap3D;

use Moose;
extends 'CoGe::Builder::Buildable';

use CoGe::Accessory::Web qw( get_defaults );
use CoGe::JEX::Workflow;
use CoGe::Accessory::Utils qw(units);
use CoGe::Builder::CommonTasks qw( create_gff_generation_job );
use CoGe::Builder::Tools::SynMap qw( add_jobs defaults gen_org_name );
use CoGe::Core::Storage qw( get_workflow_paths );
use Data::Dumper;
use File::Spec::Functions;
use JSON qw( encode_json );
use POSIX;

sub pre_build { # override superclass method
	my ($self, %params) = @_;

	# Initialize workflow -- NOTE: init => 0 means that a previous identical workflow will be reused when submitted
    $self->workflow( $params{jex}->create_workflow(name => $self->get_name, init => 0 ) );
    return unless $self->workflow;

	# Set site_url attribute
	my %opts = ( %{ defaults() }, %{ $self->params } );
	$self->site_url( $opts{tinylink} || get_query_link( $self->conf, $self->db, %opts ) );

    if ($params{requester}) { # request is from internal web page - external API requests will not have a 'requester' field
        my $page = $params{requester}->{page}; # page name used for logging
        $self->page($page) if $page;
    }

}

sub build {
	my $self = shift;
	my $xid = $self->params->{genome_id1};
	my $yid = $self->params->{genome_id2};
	my $zid = $self->params->{genome_id3};
	my ($dir1, $dir2) = sort($xid, $yid);
	my ($dir3, $dir4) = sort($xid, $zid);
	my ($dir5, $dir6) = sort($yid, $zid);

	my $SYN3DIR = $self->conf->{SYN3DIR};
	my $SCRIPTDIR = catdir( $self->conf->{SCRIPTDIR}, 'synmap' );
	my $PYTHON = $self->conf->{PYTHON} // 'python';

	my $MERGER = 'nice ' . $PYTHON . ' ' . catfile($SCRIPTDIR, 'synmerge_3.py');

	#########################################################################
	# Add SynMap jobs.
	#########################################################################
	my @genome_ids;
	my $i = 1;
	while (1) {
		if ($self->params->{'genome_id' . $i}) {
			push @genome_ids, $self->params->{'genome_id' . $i++};
		}
		else {
			last;
		}
	}
	for (my $j=1; $j<$i-1; $j++) {
		for (my $k=$j+1; $k<$i; $k++) {
			$self->params->{genome_id1} = $genome_ids[$j-1];
			$self->params->{genome_id2} = $genome_ids[$k-1];
			my %opts = ( %{ defaults() }, %{ $self->params } );
			my $resp = add_jobs(
				workflow => $self->workflow,
				db       => $self->db,
				config   => $self->conf,
				user     => $self->user,
				%opts
			);
			if ($resp) { # an error occurred
			   return 0;
			}
		}
	}

	#########################################################################
	# Add merging job.
	#########################################################################
	my $workflow = $self->workflow;

	# Input datafiles. TODO: Make this path specification based on the '$final_dagchainer_file' from SynMap.pm
	# AKB (2016-8-15) changed to reflect new naming conventions.
	#my $dot_xy_path = catfile($self->conf->{DIAGSDIR}, $dir1, $dir2, $dir1 . '_' . $dir2 . "_synteny.json");
    #my $dot_xz_path = catfile($self->conf->{DIAGSDIR}, $dir3, $dir4, $dir3 . '_' . $dir4 . "_synteny.json");
	#my $dot_yz_path = catfile($self->conf->{DIAGSDIR}, $dir5, $dir6, $dir5 . '_' . $dir6 . "_synteny.json");

	my $opts_name = '.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.Dm0.ma1.gcoords.dotplot_dots_synteny.json';  # check that this is right. Maybe there needs to be '.gcoords' after 'ma1'? or no '.Dm0.ma1'? Depends on state of final_dagchainer_file at point of assignment...
	my $dot_xy_path = catfile( $self->conf->{DIAGSDIR}, $dir1, $dir2, $dir1 . '_' . $dir2 . $opts_name );
	my $dot_xz_path = catfile( $self->conf->{DIAGSDIR}, $dir3, $dir4, $dir3 . '_' . $dir4 . $opts_name );
	my $dot_yz_path = catfile( $self->conf->{DIAGSDIR}, $dir5, $dir6, $dir5 . '_' . $dir6 . $opts_name );


	# Options.
	my $sort = $self->params->{sort};
    my $min_length = $self->params->{min_length};
    my $min_synteny = $self->params->{min_synteny};
    my $cluster = $self->params->{cluster};
    my $c_eps;
    my $c_min;
    if ($cluster eq 'true') {
        $c_eps = $self->params->{c_eps};
        $c_min = $self->params->{c_min};
    }
    my $ratio = $self->params->{ratio};
    my $r_by;
    my $r_min;
    my $r_max;
    if ($ratio ne 'false') {
        $r_by = $self->params->{r_by};
        $r_min = $self->params->{r_min};
        $r_max = $self->params->{r_max};
    }
    my $graph_out = $self->params->{graph_out};
    my $log_out = $self->params->{log_out};
    my $data_out = $self->params->{download};

	# Build command line arguments.
	my $merge_ids = ' -xid ' . $xid . ' -yid ' . $yid . ' -zid ' . $zid;
    my $merge_ins = ' -i1 ' . $dot_xy_path . ' -i2 ' . $dot_xz_path . ' -i3 ' . $dot_yz_path;
    my $merge_otp = ' -o ' . $SYN3DIR;
    my $merge_opt = ' -S ' . $sort . ' -ml ' . $min_length . ' -ms ' . $min_synteny;
    if ($cluster eq 'true') {
        $merge_opt .= ' -C -C_eps ' . $c_eps . ' -C_ms ' . $c_min;
    }
    if ($ratio ne 'false') {
        $merge_opt .= ' -R ' . $ratio . ' -Rby ' . $r_by . ' -Rmin ' . $r_min . ' -Rmax ' . $r_max;
    }

	# Add job to workflow.
	$workflow->add_job({
		description => "Identifying common points & building graph object...",
        cmd => $MERGER . $merge_ids . $merge_ins . $merge_opt . $merge_otp,
        inputs => [
			$dot_xy_path,
			$dot_xz_path,
			$dot_yz_path
		],
        outputs => [
			catfile($SYN3DIR, $graph_out),
			catfile($SYN3DIR, $log_out),
			catfile($SYN3DIR, $data_out)
		]
    });

	return 1;
}

sub get_name {
    my $self = shift;
    
    my $description;
    foreach my $key (keys $self->params) {
        next unless $key =~ /^genome_id/;
        
        my ($genome_id) = $key =~ /(\d+)$/;
        my ( $org_name ) = gen_org_name(
            db        => $self->db,
            genome_id => $genome_id
        );
        
        $description .= ($description ? 'v. ' : '') . $org_name . ' ';
    }
    
    $description .= 'Ks' if $self->params->{ks_type};
    
	return $description;
}

__PACKAGE__->meta->make_immutable;

1;
