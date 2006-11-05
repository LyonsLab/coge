package CoGe::Accessory::Web;

use strict;
use CoGe::Genome;
use Data::Dumper;
use base 'Class::Accessor';

BEGIN {
    use Exporter ();
    use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $coge $Q);
    $VERSION     = 0.1;
    @ISA         = (@ISA, qw (Exporter));
    #Give a hoot don't pollute, do not export more than needed by default
    @EXPORT      = qw ();
    @EXPORT_OK   = qw ();
    %EXPORT_TAGS = ();
    $coge = new CoGe::Genome;
 }


sub feat_name_search
  {
    my ($self, $accn, $num, $min_search) = self_or_default(@_);
    #may want to extend this function to use all the options of feat->power_search.  See FeatView.pl for an example
    $num = 1 unless defined $num;
    $min_search = 2 unless defined $min_search;
    return unless $accn;
    my $blank = qq{<input type="hidden" id="accn_select$num">};
#    print STDERR Dumper @_;
    return $blank unless length($accn) > $min_search;# || $type || $org;
    my $html;
    my %seen;
    my @opts = sort map {"<OPTION>$_</OPTION>"} grep {! $seen{$_}++} map {uc($_)} $coge->get_feature_name_obj->power_search(accn=>$accn."%");
    if (@opts > 5000)
      {
	return $blank."Search results over 1000, please refine your search.\n";
      }
    $html .= "<font class=small>Name count: ".scalar @opts."</font>\n<BR>\n";
    $html .= qq{<SELECT id="accn_select$num" SIZE="5" MULTIPLE onChange="source_search_chain(); " >\n};
    $html .= join ("\n", @opts);
    $html .= "\n</SELECT>\n";
    $html =~ s/<OPTION/<OPTION SELECTED/;
    return $blank."No results found.\n", $num unless $html =~ /OPTION/;
    return $html, $num;
  }

sub dataset_search_for_feat_name
  {
    my ($self, $accn, $num) = self_or_default(@_);
    $num = 1 unless $num;
    return unless $accn;
    my $html;
    my %sources;
    foreach my $feat ($coge->get_feats_by_name($accn))
      {
	my $val = $feat->data_info;
	my $name = $val->name;
	my $ver = $val->version;
	my $desc = $val->description;
	my $sname = $val->data_source->name;
	my $ds_name = $val->name;
	my $org = $val->org->name;
	my $title = "$org: $ds_name ($sname, v$ver)";
#	$sources{$feat->data_info->id} = $feat->data_info;
	$sources{$feat->data_info->id} = $title;
      }
    if (keys %sources)
      {
	$html .= qq{
<SELECT name = "dsid$num" id= "dsid$num" onChange="feat_search(['accn$num','dsid$num', 'args__$num'],['feat$num']);">
};
	foreach my $id (sort {$b <=> $a} keys %sources)
	  {
	    my $val = $sources{$id};
	    $html  .= qq{  <option value="$id">$val\n};
	  }
	$html .= qq{</SELECT>\n};
	my $count = scalar keys %sources;
	$html .= qq{<font class=small>($count)</font>};
      }
    else
      {
	$html .= qq{Accession not found <input type="hidden" id="dsid$num">\n<input type="hidden" id="featid$num">\n};	
      }    
    return ($html,$num);
  }

sub feat_search_for_feat_name
  {
    my ($self, $accn, $dsid, $num) = self_or_default(@_);
    return unless $dsid;
    my @feats;
    foreach my $f ($coge->get_feats_by_name($accn))
      {
	next unless $f->dataset->id == $dsid;
	next if $f->type->name =~ /CDS/i;
	next if $f->type->name =~ /RNA/i;
	push @feats, $f;
      }
    my $html;
    if (@feats)
      {
	$html .= qq{
<SELECT name = "featid$num" id = "featid$num" >
  };
	foreach my $feat (@feats)
	  {
	    my $loc = "Chr:".$feat->chr." ".$feat->genbank_location_string;
	    $loc =~ s/(complement)|(join)//g;
	    $html .= qq {  <option value="$feat->id">$loc \n};
	  }
	$html .= qq{</SELECT>\n};
	my $count = scalar @feats;
	$html .= qq{<font class=small>($count)</font>};
      }
    else
      {
	$html .=  qq{<input type="hidden" id="featid$num">\n}
      }
    return $html;
  }

sub type_search_for_feat_name
  {
    my ($accn, $dsid, $num) = self_or_default(@_);
    $num = 1 unless defined $num;
    my $html;
    my $blank = qq{<input type="hidden" id="type_name$num">};
    my %seen;
    my @opts = sort map {"<OPTION>$_</OPTION>"} grep {! $seen{$_}++} map {$_->type->name} $coge->get_features_by_name_and_dataset_id(name=>$accn, id=>$dsid);
    $html .= "<font class=small>Type count: ".scalar @opts."</font>\n<BR>\n";
    $html .= qq{<SELECT id="Type_name" SIZE="5" MULTIPLE onChange="get_anno(['accn_select$num','type_name$num', 'dsid$num'],['anno$num'])" >\n};
    $html .= join ("\n", @opts);
    $html .= "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/;
    return $blank unless $html =~ /OPTION/;
    return ($html, 1, $num);
  }

sub feat_name_search_box
  {
    my ($self, $num, $default) = self_or_default(@_);
    $num = 1 unless defined $num;
    return qq{
<input type="text" name="accn$num" id="accn$num" tabindex="1"  size="10" onKeyup ="feat_name_search_chain(0);"  onBlur="feat_name_search_chain(1);" value="$default"/>
<DIV class="" id="accn_list$num"><input type="hidden" id="accn_select$num" ></DIV>
<DIV class="" id="ds$num"></DIV>
<DIV class="" id="feat_type$num"></DIV>
<DIV class="" id="anno$num"></DIV>
};
  }

sub feat_search_box
  {
    my ($self, $num, $default) = self_or_default(@_);
    $num = 1 unless defined $num;
    return qq{
<input class="backbox" type="text" name="accn$num" id="accn$num" tabindex="$num"  size="20" onBlur="dataset_search(['accn$num','args__$num'],[feat_search_chain])" value="$default" />
<DIV id="ds$num"><input type="hidden" id="dsid$num"></DIV>
<DIV id="feat$num"><input type="hidden" id="featid$num"></DIV>
};
  }

sub feat_name_search_chain
  {
    my ($self, $lim) = self_or_default(@_);
    $lim = 4 unless defined $lim;
    return qq{
<SCRIPT language="javascript">
function feat_name_search_chain(val, num) {
 if (!num) {num = 1};
 minlen = $lim;
 if ((val == 1) && (kaj.GS('accn'+num).length > minlen)){
   return;
 } 
 kaj.GS('accn_list'+num,'fnsc_accn_list'+num);	
 kaj.GS('ds'+num,'fnsc_ds'+num);	
 kaj.GS('feat_type'+num,'fnsc_feat_type'+num);	
// kaj.GS('anno'+num,'');	
 if (val == 1 || kaj.GS('accn'+num).length > minlen) {
//   kaj.GS('DS_title','<font class="loading">Searching. . .</font>');// for "'+kaj.GS('accn')+'" type '+''+'. . .</font>');
   feat_name_search(['accn'+num, 'args__'+num],[ds_search_chain]); 
 }
}
</SCRIPT>
};
  }

sub ds_search_chain
  {
return qq{
<SCRIPT language="javascript">
function ds_search_chain (val, num) {
 kaj.GS('ds'+num,'<font class="loading">dssc_ds'+num+'Loading. . .</font>');	
 kaj.GS('feat_type'+num,'dssc_feat_type'+num);
 kaj.GS('anno'+num,'dssc_anno'+num);
 if (val) {kaj.GS('accn_list'+num,val);}
 dataset_search(['accn_select'+num, 'args__'+num], [feat_type_search_chain]);
}
</SCRIPT>

};

  }

sub feat_search_chain
  {
    return qq{
<SCRIPT language="javascript">
function feat_search_chain(val, num) {
 if (val) {
  kaj.GS('ds'+num, val);
  feat_search(['accn'+num,'dsid'+num, 'args__'+num],['feat'+num]);
 }
}
</SCRIPT>
};
  }

sub feat_type_search_chain
  {
    return qq{
<SCRIPT language="javascript">
function feat_type_search_chain (val1, num) {
// kaj.GS('FT_title','');
 kaj.GS('feat_type'+num,'<font class="loading">FTC_feat_type Loading. . .</font>');	
 kaj.GS('anno'+num,'FTC_anno');	
 if (val1) {kaj.GS('ds'+num, val1);}
// if (val2) {kaj.GS('DS_title', "Dataset");}	
 type_search(['accn_select'+num, 'dsid'+num, 'args__'+num],[get_anno_chain]);	
}
</SCRIPT>
};
}

sub kaj
  {
    return qq{
<SCRIPT language="JavaScript" type="text/javascript" src="./js/kaj.stable.js"></SCRIPT>
};
  }

sub self_or_default { #from CGI.pm
    return @_ if defined($_[0]) && (!ref($_[0])) &&($_[0] eq 'CoGe::Accessory::Web');
    unless (defined($_[0]) && 
            (ref($_[0]) eq 'CoGe::Accessory::Web' || UNIVERSAL::isa($_[0],'CoGe::Accessory::Web')) # slightly optimized for common case
            ) {
        $Q = CoGe::Accessory::Web->new unless defined($Q);
        unshift(@_,$Q);
    }
    return wantarray ? @_ : $Q;
}

sub ajax_func
  {
    return 
      (
       dataset_search=>\&dataset_search_for_feat_name,
       feat_search=>\&feat_search_for_feat_name,
       feat_name_search=>\&feat_name_search,
       type_search=> \&type_search_for_feat_name,
      );
  }

1;
