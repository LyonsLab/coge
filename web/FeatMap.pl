#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Ajax;
use CoGe::Accessory::Web;
use HTML::Template;
use Data::Dumper;
use URI::Escape;
use CoGeX;
use CoGeX::Result::Feature;
use POSIX;
use File::Temp;
use File::Basename;
use CoGe::Graphics::GenomeView;
use CoGe::Graphics;
use CoGe::Graphics::Chromosome;
use Digest::MD5 qw(md5_base64);
use File::Path;
no warnings 'redefine';
use vars
  qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $TEMPDIR $TEMPURL $FORM $USER $DATE $coge $cogeweb $COOKIE_NAME);

$P         = CoGe::Accessory::Web::get_defaults();
$ENV{PATH} = $P->{COGEDIR};
$TEMPDIR   = $P->{TEMPDIR} . "FeatMap";
$TEMPURL   = $P->{TEMPURL} . "FeatMap";
mkpath( $TEMPDIR, 0, 0777 ) unless -d $TEMPDIR;

$DATE = sprintf(
    "%04d-%02d-%02d %02d:%02d:%02d",
    sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }
      ->(localtime)
);

$FORM   = new CGI;
$DBNAME = $P->{DBNAME};
$DBHOST = $P->{DBHOST};
$DBPORT = $P->{DBPORT};
$DBUSER = $P->{DBUSER};
$DBPASS = $P->{DBPASS};
$connstr =
  "dbi:mysql:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
$coge = CoGeX->connect( $connstr, $DBUSER, $DBPASS );

$COOKIE_NAME = $P->{COOKIE_NAME};

my ($cas_ticket) = $FORM->param('ticket');
$USER = undef;
($USER) = CoGe::Accessory::Web->login_cas(
    cookie_name => $COOKIE_NAME,
    ticket      => $cas_ticket,
    coge        => $coge,
    this_url    => $FORM->url()
) if ($cas_ticket);
($USER) = CoGe::Accessory::Web->get_user(
    cookie_name => $COOKIE_NAME,
    coge        => $coge
) unless $USER;

my $pj = new CGI::Ajax(
    gen_html           => \&gen_html,
    generate_feat_info => \&generate_feat_info,
);
$pj->js_encode_function('escape');
print $pj->build_html( $FORM, \&gen_html );
#print $FORM->header;print gen_html();

sub gen_html {
    my $html;
    my ($body) = gen_body();
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( TITLE      => 'Feature Map',
                      PAGE_TITLE => 'FeatMap',
                      SUPPORT_EMAIL => $P->{SUPPORT_EMAIL},
                      HOME       => $P->{SERVER},
                      HELP       => 'FeatMap',
                      WIKI_URL   => $P->{WIKI_URL} || '',
                      USER       => $USER->display_name || '' );
    $template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
    $template->param( DATE       => $DATE );
    $template->param( BODY       => $body );
    $template->param( ADMIN_ONLY => $USER->is_admin );
    $template->param( CAS_URL    => $P->{CAS_URL} || '' );
    $html .= $template->output;
}

sub gen_body {
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'FeatMap.tmpl' );
    my $form = $FORM;
    my $no_values;
    my $sort_by_type     = $form->param('sort_type');
    my $sort_by_location = $form->param('sort_loc');
    my $fid_list         = [];
    foreach my $item ( $form->param('fid') ) {
        push @$fid_list, $item;
    }
#    my $dsid = $form->param('dsid') if $form->param('dsid');
#    my $anticodon = $form->param('anticodon') if $form->param('anticodon');
#    my $chr = $form->param('chr') if $form->param('chr');
#    my $clade = $form->param('clade') if $form->param('clade');
#    my $ftid = $form->param('ftid') if $form->param('ftid');
#    push @$feat_list, @{get_fids_from_dataset(dsid=>$dsid, ftid=>$ftid, chr=>$chr)} if $dsid;
    if (@$fid_list) {
        my $chromosome_data =
          generate_chromosome_images( feature_id_list => $fid_list );
        my $table = generate_table_data( feature_id_list => $fid_list );
        if ($chromosome_data) {
            $template->param( CHROMOSOME_LOOP => $chromosome_data );
            $template->param( FEAT_TABLE      => $table );
            return $template->output;
        }
    }
    else {
        return "No feature ids were specified.";
    }
}

sub generate_table_data {
    my %opts      = @_;
    my $feat_list = $opts{feature_id_list};
    #print STDERR Dumper \$feat_list;
    return unless @$feat_list;
    my @table;
    $feat_list = [ map { $coge->resultset("Feature")->find($_) } @$feat_list ];
    $feat_list = [
        sort {
                 $a->organism->name cmp $b->organism->name
              || $a->type->name cmp $b->type->name
              || $a->chromosome cmp $b->chromosome
              || $a->start <=> $b->start
          } @$feat_list
    ];
    foreach my $feat (@$feat_list) {
        unless ($feat) {
        #	  warn "feature id $featid failed to return a valid feature object\n";
            next;
        }
        my ($name) = $feat->names;
        push @table,
          {
            FID   => $feat->id,
            COLOR => 'red',
            FEAT_NAME =>
qq{<span class="link" title="Click for Feature Information" onclick="show_feature_info('}
              . $feat->id
              . qq{')">}
              . $name
              . "</span>",
            LOC => $feat->start,
            CHR => $feat->chr,
            ORG => $feat->organism->name,
          };
    }
    return \@table;
}

sub generate_chromosome_images {
    my %opts      = @_;
    my $feat_list = $opts{feature_id_list};
    #print STDERR Dumper \$feat_list;
    return unless @$feat_list;

    $cogeweb = CoGe::Accessory::Web::initialize_basefile( tempdir => $TEMPDIR );
    my $width          = $opts{width} || 1200;
    my $image_filename = $cogeweb->basefile;
    my $height         = ( $width / 16 ) <= 64 ? ( $width / 16 ) : 64;

    my $scale = $opts{scale} || 'linear';
    my %data;
    my $filename;
    my ( @data, @no_data );
    $feat_list = [ map { $coge->resultset("Feature")->find($_) } @$feat_list ];
    $feat_list = [
        sort {
                 $a->organism->name cmp $b->organism->name
              || $a->type->name cmp $b->type->name
              || $a->chromosome cmp $b->chromosome
              || $a->start <=> $b->start
          } @$feat_list
    ];
    foreach my $feat (@$feat_list) {
        unless ($feat) {
        #	  warn "feature id $featid failed to return a valid feature object\n";
            next;
        }

        my $org = $feat->organism->name;
        push @{ $data{$org}{feats} }, $feat;
        $filename = $image_filename . "_*.png";
    }

    foreach my $org ( sort keys %data ) {
        #first, initialize graphic
        $org =~ s/\s+$//;
        $data{$org}{image} = new CoGe::Graphics::GenomeView(
            {
                color_band_flag   => 1,
                image_width       => $width,
                chromosome_height => $height
            }
        ) unless $data{$org}{image};
        foreach my $feat ( @{ $data{$org}{feats} } ) {
            my ($chr) = $feat->chr;

            my ($dsid) = $feat->dataset->id;
            #add chromosome to graphic
            unless ( $data{$org}{chr}{$chr} ) {
                next unless $dsid;
                my $ds       = $coge->resultset('Dataset')->find($dsid);
                my $last_pos = $ds->last_chromosome_position($chr);
                $data{$org}{image}->add_chromosome(
                    name => "Chr: $chr",
                    end  => $last_pos,
                );
                $data{$org}{chr}{$chr} = 1;
            }
            my $up = $feat->strand =~ /-/ ? 0 : 1;
            my ($name) = $feat->names;
            $data{$org}{image}->add_feature(
                name     => $feat->id,
                start    => $feat->start,
                stop     => $feat->stop,
                chr      => "Chr: $chr",
                imagemap => qq/class="imagemaplink" title="/
                  . $name . qq/ (/
                  . $feat->type->name
                  . qq/)" onclick="show_feature_info('/
                  . $feat->id . qq/')"/,
                up => $up,
                #color=>color_feat_by_type(feat=>$feat),
            );
        }
    }
    my $count = 1;
    foreach my $org ( sort keys %data ) {
        if ( $data{$org}{image} ) {
            my $image_file;
            my $image_map;
            unless ( $data{$org}{skip} ) {
                my $x;
                $image_file = $cogeweb->basefile . "_$count.png";
                ( $x, $image_file ) =
                  CoGe::Accessory::Web::check_taint($image_file);

                $data{$org}{image}->image_width($width);
                $data{$org}{image}->chromosome_height($height);
                $data{$org}{image}->show_count(1);
                $data{$org}{image}->generate_png( filename => $image_file );
                $image_map =
                  $data{$org}{image}->generate_imagemap(
                    mapname => $cogeweb->basefilename . "_" . $count );
                my $map_file = $cogeweb->basefile . "_$count.map";
                ( $x, $map_file ) =
                  CoGe::Accessory::Web::check_taint($map_file);
                open( MAP, ">$map_file" );
                print MAP $image_map;
                close MAP;
            }
            else {
                my $x;
                $image_file = $data{$org}{image} . "_$count.png";
                $image_map  = get_map( $cogeweb->basefile . "_$count.map" );
                #print STDERR $image_map_large,"\n";
                ( $x, $image_file ) =
                  CoGe::Accessory::Web::check_taint($image_file);
            }

            $image_file =~ s/$TEMPDIR/$TEMPURL/;

            push @data,
              {
                ORG_NAME  => $org,
                CHR_IMAGE => "<img src=$image_file ismap usemap='"
                  . $cogeweb->basefilename . "_"
                  . "$count' border=0>$image_map"
              };
            $count++;
        }
        else {
            push @no_data,
              {     ORG_NAME => "No Hits: <a href="
                  . $data{$org}{file}
                  . " target=_new>$org</a>" };
        }
    }

    return return [ @data, @no_data ];
}

sub get_map {
    my $file = shift;
    my $map;
    open( IN, $file ) || die "$!";
    while (<IN>) {
        $map .= $_;
    }
    close IN;
    return $map;
}

sub generate_feat_info {
    my $fid = shift;
    #print STDERR $fid, "\n";
    my $width    = 1000;
    my $feat     = $coge->resultset('Feature')->find($fid);
    my $hpp      = $feat->annotation_pretty_print_html();
    my $cs       = new CoGe::Graphics::Chromosome();
    my $ds       = $coge->resultset('Dataset')->find( $feat->dataset->id );
    my $last_pos = $ds->last_chromosome_position( $feat->chr );
    $cs->chr_length($last_pos);
    $cs->iw($width);
    $cs->ih(200);
    $cs->draw_chromosome(1);
    $cs->draw_ruler(1);
    $cs->draw_chr_end(0);
    $cs->mag(0);
    $cs->mag_off(1);
    $cs->minor_tick_labels(0);
    $cs->major_tick_labels(1);
    $cs->draw_hi_qual(0);
    $cs->padding(2);
    $cs->set_region(
        start    => $feat->start - 1000,
        stop     => $feat->stop + 1000,
        forcefit => 1
    );
    $cs->auto_zoom(0);
    $cs->feature_height(10);
    $cs->overlap_adjustment(0);
    $cs->feature_labels(1);
    my $graphics = new CoGe::Graphics;
    $graphics->process_features(
        c      => $cs,
        layers => {
            features => { gene => 1, cds => 1, mrna => 1, rna => 1, cns => 1 }
        },
        ds   => $ds,
        chr  => $feat->chr,
        coge => $coge
    );
    $cs->overlap_adjustment(1);
    my $sub_file = $cogeweb->basefile . ".s." . $feat->id . ".png";
    $cs->generate_png( file => $sub_file );
    my ($name) = $feat->names;
    my $subject_link = qq{
Feature: $name<br>
<a title='Click for Interactive Genome View' href = 'GenomeView.pl?chr=}
      . $feat->chr
      . qq{&ds=}
      . $feat->dataset->id . qq{&x=}
      . $feat->start
      . qq{&z=5' target=_new border=0><img src=$sub_file border=0></a>
};
    $subject_link =~ s/$TEMPDIR/$TEMPURL/;
    #print STDERR "made it this far!\n";
    return $hpp, $subject_link;
}

sub generate_feat_info {
    my $featid = shift;
    my ($feat) = $coge->resultset("Feature")->find($featid);
    unless ( ref($feat) =~ /Feature/i ) {
        return "Unable to retrieve Feature object for id: $featid";
    }
    my $html =
qq{<a href="#check_$featid" onClick="\$('#feat_info').slideToggle(pageObj.speed);" style="float: right;"><img src='picts/delete.png' width='16' height='16' border='0'></a>};
    $html .= $feat->annotation_pretty_print_html();
    return $html;
}
