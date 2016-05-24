#!/usr/bin/perl
use strict;
use warnings FATAL => 'all';

# TODO: Loop out the 3 dotplot dots
sub dotplot_dots {
    my %opts = @_;
    my $genome_idX = $opts{genome_idX};
    my $genome_idY = $opts{genome_idY};
    my $genome_idZ = $opts{genome_idZ};
    my $ksfile_xy = $opts{ksfile_xy};
    my $ksfile_xz = $opts{ksfile_xz};
    my $ksfile_yz = $opts{ksfile_yz};
    my $hide = $opts{hide};
    my $min_len = $opts{min_len};

    my $DIAGSDIR = $CONF->{DIAGSDIR};
    my $SYN3DIR = $CONF->{SYN3DIR};
    my $SCRIPTDIR = $CONF->{SCRIPTDIR};

    my $dotlog = 'dotplot_dots_log.json';
    my $merge_log = $genome_idX . '_' . $genome_idY . '_' . $genome_idZ . '_log.json';

    my ( $dir1, $dir2 ) = sort ( $genome_idX, $genome_idY );
    my ( $dir3, $dir4 ) = sort ( $genome_idX, $genome_idZ );
    my ( $dir5, $dir6 ) = sort ( $genome_idY, $genome_idZ );

    my $workflow = $JEX->create_workflow( name => "Finding Syntenic Points", init => 1 );

    # Build/add dotplot_dots XY job.
    my $cmd_xy = catfile($SCRIPTDIR, 'dotplot_dots.py') . ' ' . $ksfile_xy;
    my $dots_xy = $dir1 . '_' . $dir2 . '_synteny.json';
    my $dots_xy_path = catfile($DIAGSDIR, $dir1, $dir2, $dots_xy);
    my $outputs_xy = [catfile($DIAGSDIR, $dir1, $dir2, $dotlog), $dots_xy_path];
    $workflow->add_job({
        cmd => $cmd_xy,
        outputs => $outputs_xy,
        description => "running XY dotplot_dots...",
    });

    # Build/add dotplot_dots XZ job.
    my $cmd_xz = catfile($SCRIPTDIR, 'dotplot_dots.py') . ' ' . $ksfile_xz;
    my $dots_xz = $dir3 . '_' . $dir4 . '_synteny.json';
    my $dots_xz_path = catfile($DIAGSDIR, $dir3, $dir4, $dots_xz);
    my $outputs_xz = [catfile($DIAGSDIR, $dir3, $dir4, $dotlog), $dots_xz_path];
    $workflow->add_job({
        cmd => $cmd_xz,
        outputs => $outputs_xz,
        description => "running XZ dotplot_dots...",
    });

    # Build/add dotplot_dots YZ job.
    my $cmd_yz = catfile($SCRIPTDIR, 'dotplot_dots.py') . ' ' . $ksfile_yz;
    my $dots_yz = $dir5 . '_' . $dir6 . '_synteny.json';
    my $dots_yz_path = catfile($DIAGSDIR, $dir5, $dir6, $dots_yz);
    my $outputs_yz = [catfile($DIAGSDIR, $dir5, $dir6, $dotlog), $dots_yz_path];
    $workflow->add_job({
        cmd => $cmd_yz,
        outputs => $outputs_yz,
        description => "running YZ dotplot_dots...",
    });

    # Build/add three_dots_merge job.
    my $merge_ins = ' -i1 ' . @$outputs_xy[1] . ' -i2 ' . @$outputs_xz[1] . ' -i3 ' . @$outputs_yz[1];
    my $merge_ots = ' -o ' . $SYN3DIR;
    my $merge_gids = ' -xid ' . $genome_idX . ' -yid ' . $genome_idY . ' -zid ' . $genome_idZ;
    my $merge_opts = '';
    if ($hide eq 'true') { $merge_opts .= ' -P' }
    if ($min_len > 0) { $merge_opts .= ' -M ' . $min_len; }

    my $merge_cmd = catfile($SCRIPTDIR, 'three_dots_merge.py') . $merge_ins . $merge_ots . $merge_gids . $merge_opts;
    my $merge_in = [$dots_xy_path, $dots_xz_path, $dots_yz_path];
    my $dots_xyz = $genome_idX . '_' . $genome_idY . '_' . $genome_idZ . '_dots.json';
    my $hist_xyz = $genome_idX . '_' . $genome_idY . '_' . $genome_idZ . '_histogram.json';
    my $merge_out = [catfile($SYN3DIR, $dots_xyz), catfile($SYN3DIR, $merge_log), catfile($SYN3DIR, $hist_xyz)];
    $workflow->add_job({
        cmd => $merge_cmd,
        inputs => $merge_in,
        outputs => $merge_out,
        description => "merging XYZ dots..."
    });

    my $response = $JEX->submit_workflow($workflow);
    return encode_json(
        {
            id => $response->{id},
            status  => $response->{status},
            success => $JEX->is_successful($response)
            ? JSON::true
            : JSON::false
        }
    );
}