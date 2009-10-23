#!/usr/bin/perl -w

use strict;
use GD;
use Getopt::Long;
use CoGeX;
use DBI;
use Data::Dumper;
use DBI;
use POSIX;

use vars qw($dagfile $alignfile $width $link $min_chr_size $dsgid1 $dsgid2 $help $coge $graphics_context $CHR1 $CHR2 $basename $link_type $flip $grid $ks_db $ks_type $log $MAX $MIN $assemble $axis_metric);


GetOptions(
	   "dagfile|d=s"=>\$dagfile, #all dots
	   "alignfile|a=s"=>\$alignfile, #syntenic dots
	   "width|w=i"=>\$width,
	   "link|l=s"=>\$link,
	   "link_type|lt=i"=>\$link_type,
	   "min_chr_size|mcs=i"=>\$min_chr_size,
	   "dsgid1|dsg1=i"=>\$dsgid1,
	   "dsgid2|dsg2=i"=>\$dsgid2,
	   "help|h" =>\$help,
	   "chr1|c1=s"=>\$CHR1,
	   "chr2|c2=s"=>\$CHR2,
	   "basename|b=s"=>\$basename,
	   "flip|f=i"=>\$flip,
	   "grid|g=i"=>\$grid,
	   "ksdb|ks_db=s"=>\$ks_db,
	   "ks_type|kst=s"=>\$ks_type,
	   "max=s"=>\$MAX,
	   "min=s"=>\$MIN,
	   "log=s"=>\$log,
	   "assemble=s"=>\$assemble,
	   "axis_metrix|am=s"=>\$axis_metric,
	   );


usage() if $help;
usage() unless -r $dagfile;

#set a default for this, make sure it is uppercase
$ks_type = "kS" unless $ks_type;
$ks_type = uc($ks_type) if $ks_type;

my $ks_hist = "/opt/apache/CoGe/bin/ks_histogram.pl";

$basename = "test" unless $basename;
$width = 1024 unless $width;

$coge = CoGeX->dbconnect();



my ($org1_order, $org2_order) = parse_syn_blocks(file=>$alignfile) if $assemble;

my $skip_non_ordered = $assemble && $assemble == 2 ? 1 : 0;
my $org1info = get_dsg_info(dsgid=>$dsgid1, chr=>$CHR1, minsize=>$min_chr_size, order=>$org1_order, metric=>$axis_metric, skip_non_ordered=>$skip_non_ordered);
my $org1length =0;
map {$org1length+=$_->{length}} values %$org1info;
my $org2info = get_dsg_info(dsgid=>$dsgid2, chr=>$CHR2, minsize=>$min_chr_size, order=>$org2_order, metric=>$axis_metric, skip_non_ordered=>$skip_non_ordered);
my $org2length =0;
map {$org2length+=$_->{length}} values %$org2info;
($org1info, $org1length, $dsgid1, $org2info, $org2length, $dsgid2) = ($org2info, $org2length, $dsgid2, $org1info, $org1length, $dsgid1) if $flip;
($CHR1, $CHR2) = ($CHR2, $CHR1) if $flip && ($CHR1 || $CHR2);
my $height = sprintf("%.0f", $width*$org2length/$org1length);
$height = $width if ($height > 20*$width) || ($height <  $width/20);
my $x_bp_per_pix = sprintf("%.0f", $org1length/$width);
my $x_pix_per_bp = 1/$x_bp_per_pix;
my $y_bp_per_pix = sprintf("%.0f", $org2length/$height);
my $y_pix_per_bp = 1/$y_bp_per_pix;

#Generate new graphics context and fill the background
my $graphics_context = new GD::Image($width, $height);
my $white = $graphics_context->colorResolve(255,255,255);
$graphics_context->fill(1,1,$white);


if ($ks_db)
  {
    my $cmd = $ks_hist;
    $cmd .= " -db $ks_db";
    $cmd .= " -ks_type $ks_type";
    $cmd .= " -log" if $log;
    $cmd .= " -pf $alignfile";
    $cmd .= " -chr1 $CHR1" if defined $CHR1;
    $cmd .= " -chr2 $CHR2" if defined $CHR2;
    $cmd .= " -min $MIN" if defined $MIN;
    $cmd .= " -max $MAX" if defined $MAX;
    $cmd .= " -o $basename";
    $cmd .= ".$MIN" if defined $MIN;
    $cmd .= ".$MAX" if defined $MAX;
    $cmd .= ".hist.png";
    `$cmd &`;
  }


$basename .= ".$MIN" if defined $MIN;
$basename .= ".$MAX" if defined $MAX;

my $pairs = get_pairs(file=>$alignfile, chr1=>$CHR1, chr2=>$CHR2) if $alignfile && $ks_db && -r $alignfile && -r $ks_db;

#Magic happens here.
#Link type seems to indicate the type of tile; i.e. a 'master' (a large, all chromosome) or a blow up of a two chromosome intersection
#draw_chromosome_grid draws either the black chomosome lines, or the light green tile lines, so its always called in addition to the draw_dots function.
draw_chromosome_grid(gd=>$graphics_context, org1=>$org1info, org2=>$org2info, x_pix_per_bp=>$x_pix_per_bp, y_pix_per_bp=>$y_pix_per_bp, link=>$link, link_type=>$link_type, flip=>$flip, grid=>$grid);

#get syntenic gene pairs for ks_data (if needed)
my $ksdata = get_ksdata(ks_db=>$ks_db, ks_type=>$ks_type, chr1=>$CHR1, chr2=> $CHR2, pairs=>$pairs) if $ks_db && -r $ks_db;

#draw dots for all matches
draw_dots(gd=>$graphics_context, file=>$dagfile, org1=>$org1info, org2=>$org2info, x_pix_per_bp=>$x_pix_per_bp, y_pix_per_bp=>$y_pix_per_bp, link_type => $link_type, dsgid1=>$dsgid1, dsgid2=>$dsgid2, flip=>$flip, metric=>$axis_metric);

#add syntenic gene pairs
my $add = 1 if $dsgid1 eq $dsgid2;
draw_dots(gd=>$graphics_context, file=>$alignfile, org1=>$org1info, org2=>$org2info, x_pix_per_bp=>$x_pix_per_bp, y_pix_per_bp=>$y_pix_per_bp, color=>$graphics_context->colorResolve(0,150,0), rev_color=>$graphics_context->colorResolve(200,0,0), size=>2, add_inverse=>$add, flip=>$flip, ksdata=>$ksdata, ks_type=>$ks_type, log=>$log, metric=>$axis_metric);

#Write out graphics context - the generated dot plot - to a .png file
open (OUT, ">".$basename.".png") || die "$!";
binmode OUT;
print OUT $graphics_context->png;
close OUT;

#generate_historgram of ks values if necessary

#This function appears to parse dagchainer output, generated in SynMap.pl, and draw the results to the GD graphics context.
sub draw_dots
  {
    my %opts = @_;
    my $file = $opts{file};
    my $graphics_context = $opts{gd};
    my $org1 = $opts{org1};
    my $org2 = $opts{org2};
    my $x_pix_per_bp=$opts{x_pix_per_bp};
    my $y_pix_per_bp=$opts{y_pix_per_bp};
    my $size = $opts{size} || 1;
    my $color = $opts{color};
    my $rev_color = $opts{rev_color};
    my $add_inverse = $opts{add_inverse};
    my $link_type = $opts{link_type} || 0;
    my $dsgid1 = $opts{dsgid1};
    my $dsgid2 = $opts{dsgid2};
    my $flip = $opts{flip};
    my $ksdata = $opts{ksdata};
    my $log = $opts{log}; #log normalize ksdata for display?
    my $metric= $opts{metric};
    my $has_ksdata = keys %$ksdata ? 1 : 0;
    #min and max will be log normalized if log flag is set
    my ($max, $min) = get_range(data=>$ksdata, min=>$MIN, max=>$MAX, log=>$log) if $has_ksdata;
    my $range = $max-$min if $has_ksdata;
    $color = $graphics_context->colorResolve(150,150,150) unless $color;
    $rev_color = $color unless $rev_color;
    open (IN, $file) || die "$!";
    my @feats;
    my %points;
    my @points;
    while (<IN>)
      {
	chomp;
	next if /^#/;
	next unless $_;
	my @line = split /\t/;
	my $val;
	my @item1 = split/\|\|/, $line[1];
	my @item2 = split/\|\|/, $line[5];
	my $fid1 = $item1[6];
	my $fid2 = $item2[6];
	my $use_color;
	if ($item1[4] && $item2[4])
	  {
	    $use_color = $item1[4] eq $item2[4]? $color : $rev_color;
	  }
	else
	  {
	    $use_color = $color;
	  }

	if ($has_ksdata)
	  {
	    $val = $ksdata->{$fid1}{$fid2};
	    if (defined $val && $val=~ /\d/) 
	      {
		my $orig_val = $val;
		$val = $min if $log && $val == 0;
		if ($log)
		  {
		    if ($val <= 0)
		      {
			$val = 0; #set to minimum color value
		      }
		    else
		      {
			$val = (log10($val)-$min)/$range;
		      }
		  }
		else
		  {
		    $val = ($val-$min)/$range;
		  }
		$val = sprintf("%.4f", $val);
		$use_color = get_color(val=>$val); #val is 0<=x<=1
		$use_color = $graphics_context->colorResolve(@$use_color);
	      }
	    else
	      {
		#don't have ks data -- skip drawing this dot!
		next;
#		print Dumper $ksdata->{$item1[6]}{$item2[6]};
#		$use_color = $graphics_context->colorResolve(0,0,0);
	      }
	  }
	if ($flip)
	  {
	    my @tmp = @line[0..3];
	    @line[0..3] = @line[4..7];
	    @line[4..7] = @tmp;
	    ($fid1, $fid2) = ($fid2, $fid1);
	  }
	my ($a, $chr1) = split /_/,$line[0],2;
	my ($b, $chr2) = split /_/,$line[4],2;
	my $special = 0; #stupid variable name so that if we are viewing a single chr to single chr comparison within the same organism, this will make collinear matches appear on the inverse section 
	if ($CHR1 && $CHR2)
	  {
	    if ($add_inverse && $CHR1 eq $chr2 && $CHR2 eq $chr1)
	      {
		$special = 1;
		($chr1, $chr2) = ($chr2, $chr1);
	      }
	    else
	      {
		next if $chr1 ne $CHR1;
		next if $chr2 ne $CHR2;
	      }
	  }

	next unless $org1->{$chr1} && $org2->{$chr2}; #sometimes there will be data that is skipped, e.g. where chromosome="random";
	my ($xmin, $ymin);
	my ($midx, $midy);
	if ($metric && $metric =~ /gene/i )
	  {
	    next unless $fid1 && $fid2;
	    next unless $org1->{$chr1} && $org1->{$chr1}{gene} ;
	    next unless $org2->{$chr2} && $org2->{$chr2}{gene} ;
	    $xmin = $org1->{$chr1}{gene}{$fid1};
	    $ymin = $org2->{$chr2}{gene}{$fid2};
	    next unless $xmin && $ymin;
	    $midx = $org1->{$chr1}{rev}? sprintf("%.0f",$org1->{$chr1}{start}+$org1->{$chr1}{length}-($xmin)): sprintf("%.0f",$org1->{$chr1}{start}+$xmin);
	    $midy = $org2->{$chr2}{rev}? sprintf("%.0f",$org2->{$chr2}{start}+$org2->{$chr2}{length}-($ymin)): sprintf("%.0f",$org2->{$chr2}{start}+$ymin);

	  }
	else #using nucleotides
	  {
	    ($xmin) = sort ($line[2], $line[3]);
	    ($ymin) = sort ($line[6], $line[7]);
	    $midx = $org1->{$chr1}{rev}? sprintf("%.0f",$org1->{$chr1}{start}+$org1->{$chr1}{length}-($xmin+abs($line[3]-$line[2])/2)): sprintf("%.0f",$org1->{$chr1}{start}+$xmin+abs($line[3]-$line[2])/2);
	    $midy = $org2->{$chr2}{rev}? sprintf("%.0f",$org2->{$chr2}{start}+$org2->{$chr2}{length}-($ymin+abs($line[7]-$line[6])/2)): sprintf("%.0f",$org2->{$chr2}{start}+$ymin+abs($line[7]-$line[6])/2);

	  }


	my $x = sprintf("%.0f",$midx*$x_pix_per_bp);
	my $y = sprintf("%.0f",$midy*$y_pix_per_bp);
	($x,$y) = ($y, $x) if $special;

#	$graphics_context->arc($x, $graphics_context->height-$y, $size, $size, 0, 360, $use_color);
#	$graphics_context->arc($y, $graphics_context->height-$x, $size, $size, 0, 360, $use_color) if ($add_inverse && !$CHR1 && $x ne $y);
#	$graphics_context->arc($y, $graphics_context->height-$x, $size, $size, 0, 360, $use_color) if ($add_inverse && $chr1 eq $chr2 && $x ne $y);
	$val = 0 unless $val; #give it some value for later sorting
	push @points, [ $x, $graphics_context->height-$y, $size, $size, 0, 360, $use_color, $val];
	push @points, [ $y, $graphics_context->height-$x, $size, $size, 0, 360, $use_color, $val] if ($add_inverse && !$CHR1 && $x ne $y);
	push @points, [ $y, $graphics_context->height-$x, $size, $size, 0, 360, $use_color, $val] if ($add_inverse && $chr1 eq $chr2 && $x ne $y);
	if ($link_type == 1)
	  {
	    #working here.  Need to build a GEvo link using datasets/chr/position if dealing with genomic data.
	    my $link = qq{/CoGe/GEvo.pl?drup1=50000&drdown1=50000&drup2=50000&drdown2=50000};
	    if ($fid1)
	      {
		$link .= qq{;fid1=$fid1};
	      }
	    else
	      {
		#just added the coge->ds object to $org
		$item1[6] = "Chr: ".$item1[0]." ".commify($item1[1])." - ".commify($item1[2]);
		my $chr = $item1[0];
		$link .= qq{;dsgid1=$dsgid1;x1=$midx;chr1=$chr};
	      }
	    if ($fid2)
	      {
		$link .= qq{;fid2=$fid2};
	      }
	    else
	      {
		#just added the coge->ds object to $org
		$item2[6] = "Chr: ".$item2[0]." ".commify($item2[1])." - ".commify($item2[2]);
		my $chr = $item2[0];
		$link .= qq{;dsgid2=$dsgid2;x2=$midy;chr2=$chr};
	      }
	    unless ($points{$x}{$y}) #cuts down on the size of the image map
	      {
		push @feats, [$x, $graphics_context->height-$y, $item1[6],$item2[6], $link];
		$points{$x}{$y}=1;
	      }
	  }
      }
    close IN;
    @points = sort {$b->[-1] <=> $a->[-1]} @points if ($has_ksdata);
    foreach my $point (@points)
      {
	my $val = pop @$point;
	$graphics_context->arc(@$point);
      }
    if ($link_type == 1)
      {
	#Okay, now generate the HTML document that contains the click map for the image
	open (OUT, ">".$basename.".html") || die "$!";
	my ($org1name, $org2name);
	if ($dsgid1)
	  {
	    my $dsg = $coge->resultset('DatasetGroup')->find($dsgid1);
	    $org1name = $dsg->organism->name." (v".$dsg->version.")";
	  }
	if ($dsgid2)
	  {
	    my $dsg = $coge->resultset('DatasetGroup')->find($dsgid2);
	    $org2name = $dsg->organism->name." (v".$dsg->version.")";
	  }
	
#	my $org1name = $coge->resultset('DatasetGroup')->find($dsgid1)->organism->name if $dsgid1;
	my $org1length = commify($org1->{$CHR1}->{length})." bp";
#	my $org2name = $coge->resultset('DatasetGroup')->find($dsgid2)->organism->name if $dsgid2;
	my $org2length = commify($org2->{$CHR2}->{length})." bp";
	
	my ($img) = $basename =~ /([^\/]*$)/;
	
	#Print out the page header, Javascript
	print OUT qq{
<html><head>
<link href="/CoGe/css/styled.css" type="text/css" rel="stylesheet"/>
<SCRIPT type="text/javascript">

var ajax = [];function pjx(args,fname,method) { this.target=args[1]; this.args=args[0]; method=(method)?method:'GET'; if(method=='post'){method='POST';} this.method = method; this.r=ghr(); this.url = this.getURL(fname);}function formDump(){ var all = []; var fL = document.forms.length; for(var f = 0;f<fL;f++){ var els = document.forms[f].elements; for(var e in els){ var tmp = (els[e].id != undefined)? els[e].id : els[e].name; if(typeof tmp != 'string'){continue;} if(tmp){ all[all.length]=tmp} } } return all;}function getVal(id) { if (id.constructor == Function ) { return id(); } if (typeof(id)!= 'string') { return id; } var element = document.getElementById(id); if( !element ) { for( var i=0; i<document.forms.length; i++ ){ element = document.forms[i].elements[id]; if( element ) break; } if( element && !element.type ) element = element[0]; } if(!element){ alert('ERROR: Cant find HTML element with id or name: ' + id+'. Check that an element with name or id='+id+' exists'); return 0; } if(element.type == 'select-one') { if(element.selectedIndex == -1) return; var item = element[element.selectedIndex]; return item.value || item.text; } if(element.type == 'select-multiple') { var ans = []; var k =0; for (var i=0;i<element.length;i++) { if (element[i].selected || element[i].checked ) { ans[k++]= element[i].value || element[i].text; } } return ans; } if(element.type == 'radio' || element.type == 'checkbox'){ var ans =[]; var elms = document.getElementsByTagName('input'); var endk = elms.length ; var i =0; for(var k=0;k<endk;k++){ if(elms[k].type== element.type && elms[k].checked && (elms[k].id==id||elms[k].name==id)){ ans[i++]=elms[k].value; } } return ans; } if( element.value == undefined ){ return element.innerHTML; }else{ return element.value; }}function fnsplit(arg) { var url=""; if(arg=='NO_CACHE'){return '&pjxrand='+Math.random()} if((typeof(arg)).toLowerCase() == 'object'){ for(var k in arg){ url += '&' + k + '=' + arg[k]; } }else if (arg.indexOf('__') != -1) { arga = arg.split(/__/); url += '&' + arga[0] +'='+ escape(arga[1]); } else { var res = getVal(arg) || ''; if(res.constructor != Array){ res = [res] } for(var i=0;i<res.length;i++) { url += '&args=' + escape(res[i]) + '&' + arg + '=' + escape(res[i]); } } return url;}pjx.prototype = { send2perl : function(){ var r = this.r; var dt = this.target; this.pjxInitialized(dt); var url=this.url; var postdata; if(this.method=="POST"){ var idx=url.indexOf('?'); postdata = url.substr(idx+1); url = url.substr(0,idx); } r.open(this.method,url,true); ; if(this.method=="POST"){ r.setRequestHeader("Content-Type", "application/x-www-form-urlencoded"); r.send(postdata); } if(this.method=="GET"){ r.send(null); } r.onreadystatechange = handleReturn; }, pjxInitialized : function(){}, pjxCompleted : function(){}, readyState4 : function(){ var rsp = unescape(this.r.responseText); /* the response from perl */ var splitval = '__pjx__'; /* to split text */ /* fix IE problems with undef values in an Array getting squashed*/ rsp = rsp.replace(splitval+splitval+'g',splitval+" "+splitval); var data = rsp.split(splitval); dt = this.target; if (dt.constructor != Array) { dt=[dt]; } if (data.constructor != Array) { data=[data]; } if (typeof(dt[0])!='function') { for ( var i=0; i<dt.length; i++ ) { var div = document.getElementById(dt[i]); if (div.type =='text' || div.type=='textarea' || div.type=='hidden' ) { div.value=data[i]; } else{ div.innerHTML = data[i]; } } } else if (typeof(dt[0])=='function') { dt[0].apply(this,data); } this.pjxCompleted(dt); }, getURL : function(fname) { var args = this.args; var url= 'fname=' + fname; for (var i=0;i<args.length;i++) { url=url + args[i]; } return url; }};handleReturn = function() { for( var k=0; k<ajax.length; k++ ) { if (ajax[k].r==null) { ajax.splice(k--,1); continue; } if ( ajax[k].r.readyState== 4) { ajax[k].readyState4(); ajax.splice(k--,1); continue; } }};var ghr=getghr();function getghr(){ if(typeof XMLHttpRequest != "undefined") { return function(){return new XMLHttpRequest();} } var msv= ["Msxml2.XMLHTTP.7.0", "Msxml2.XMLHTTP.6.0", "Msxml2.XMLHTTP.5.0", "Msxml2.XMLHTTP.4.0", "MSXML2.XMLHTTP.3.0", "MSXML2.XMLHTTP", "Microsoft.XMLHTTP"]; for(var j=0;j<=msv.length;j++){ try { A = new ActiveXObject(msv[j]); if(A){ return function(){return new ActiveXObject(msv[j]);} } } catch(e) { } } return false;}function jsdebug(){ var tmp = document.getElementById('pjxdebugrequest').innerHTML = "<br><pre>"; for( var i=0; i < ajax.length; i++ ) { tmp += '<a href= '+ ajax[i].url +' target=_blank>' + decodeURI(ajax[i].url) + ' </a><br>'; } document.getElementById('pjxdebugrequest').innerHTML = tmp + "</pre>";}function get_pair_info() { var args = get_pair_info.arguments; for( var i=0; i<args[0].length;i++ ) { args[0][i] = fnsplit(args[0][i]); } var l = ajax.length; ajax[l]= new pjx(args,"get_pair_info",args[2]); ajax[l].url = '/CoGe/SynMap.pl?' + ajax[l].url; ajax[l].send2perl(); ;}
4//]]>



</SCRIPT>
<script src="/CoGe/js/jquery-1.3.2.min.js"></script>
<script src="/CoGe/js/xhairs.js"></script> <!--This needs to back out farther since the HTML file is much deeper.-->
<SCRIPT language="JavaScript" type="text/javascript" src="/CoGe/js/jquery-ui-1.7.2.custom.min.js"></SCRIPT>
<link rel="stylesheet" type="text/css" href="/CoGe/css/jquery-ui-1.7.2.custom.css" />
<link rel="stylesheet" type="text/css" href="/CoGe/css/jquery-ui-coge-supplement.css" />
<script type="text/javascript">
\$(function() {\$("#pair_info").draggable();});

/*This array follows the following format:
type (circle, rect), array of coords, href link, mouseover details
Coords will be in [x,y,radius] for cricles and [x1,y1,x2,y2] (diagonal points) for rectangles*/
};

	#Loop through the features list and print out the click map info
	my $temp_arrayindex_counter = 0;
	foreach my $item (@feats)
	  {
	    my ($x, $y, $f1, $f2, $link) = @$item;
	    #<area shape='circle' coords='$x, $y, 2' href='$link' onMouseOver="get_pair_info(['args__$f1','args__$f2'],['pair_info']);" target='_blank' >
	    print OUT qq{
location_list[$temp_arrayindex_counter] = ['circle', [$x, $y, 2], '$link', 'get_pair_info(["args__$f1","args__$f2"],["pair_info"])'];};
	    $temp_arrayindex_counter++;
	  }
	
	#Print out the canvas tag, onLoad call, and other stuff.
	print OUT qq{
		
/*The body onLoad function in unreliable (sometimes does not run), as is simply running the loadstuff function (DOM may not be fully loaded).
Ergo, we rely on jQuery to detect when the DOM is fully loaded, and then run the loadstuff() function.*/
\$(document).ready(function() {
	loadstuff('$img.png');
});
</script>
</head><body>
<canvas id="myCanvas" border="0" style='position: absolute;left: 0px;top:45px;' width="12" height="12" onmousemove="trackpointer(event);" onmousedown="trackclick(event, 'yes');">
	Your browser does not have support for Canvas.
</canvas>	
<span class=xsmall style='position: absolute;left: 0px;top: 45px;'>y-axis: $org2name: $CHR2 ($org2length)</span>
<DIV id=pair_info class="ui-widget-content" style='position: absolute;left: 0px;top: 0px;'></DIV>

};

	#Print out the nametag
	my $pos = $graphics_context->height+45;
	$pos .='px';
	print OUT qq{</map>};
	if (-r $basename.".hist.png")
	  {
	    print OUT qq{
<span style='position: absolute;left: $width; top: 45px'>
<a class="small" href= "$img.hist.png" target=_new>Histogram of synonymous substitutions</a>
<img src= "$img.hist.png" >
</span>
}; 
	  }
	print OUT qq{
<span class=xsmall style='position: absolute;left: 0px;top: $pos;'>x-axis $org1name: $CHR1 ($org1length)
   <a href = "$img.png" target=_new>Link to image</a>

};
	print OUT qq{
</span></body></html>
};
	close OUT;
      }


  }

#This method draws the grid lines indicating cromosome locations, and write out the click map to the HTML file.
#On 'zoom in' tile, this method draws the light green seperator lines.
sub draw_chromosome_grid
  {
    my %opts = @_;
    my $graphics_context = $opts{gd};
    my $org1 = $opts{org1};
    my $org2 = $opts{org2};
    my $x_pix_per_bp=$opts{x_pix_per_bp};
    my $y_pix_per_bp=$opts{y_pix_per_bp};
    my $link = $opts{link};
    my $link_type = $opts{link_type};
    my $grid = $opts{grid};
    my $height = $graphics_context->height;
    my $width = $graphics_context->width;
    my $black = $graphics_context->colorResolve(0,0,0);
    my $red = $graphics_context->colorResolve(255,0,0);
    my $span_color = $graphics_context->colorResolve(200,255,200);
    $graphics_context->line(0,0, $graphics_context->width, 0, $black);
    $graphics_context->line(0, $graphics_context->height-1, $graphics_context->width, $graphics_context->height-1, $black);
    $graphics_context->line(0,0, 0, $graphics_context->height, $black);
    $graphics_context->line($graphics_context->width-1,0, $graphics_context->width-1, $graphics_context->height, $black);
    #make x-axis chromosome deliniators
    my %data;
    my $pchr;
    my $pv;
    foreach my $chr (sort {$org1->{$a}{start}<=>$org1->{$b}{start} } keys %$org1)
      {
	my $x = sprintf("%.0f",$org1->{$chr}{start}*$x_pix_per_bp);
	next if $pv && $x == $pv-1;
	if ($grid)
	  {
	    my $tmp = $x;
	    my $span = sprintf("%.0f",($org1->{$chr}{length}/10*$x_pix_per_bp));
	    for (1..9)
	      {
		$tmp+=$span;
		$graphics_context->line($tmp, 0, $tmp, $graphics_context->height, $span_color);		#Draw green lines?
	      }
	  }
	$graphics_context->line($x, 0, $x, $graphics_context->height, $black);					#Draw black line
	my $str_color = $org1->{$chr}{rev} ? $red : $black;
	$graphics_context->string(gdSmallFont, $x+2, $graphics_context->height-15, $chr, $str_color);		#Draws name of chromosome
	$data{x}{$pchr}=[$pv,$x] if $pchr;
	$pchr=$chr;
	$pv = $x+1;
      }
    $data{x}{$pchr}=[$pv,$graphics_context->width];
    $pchr = undef;
    $pv=undef;
    foreach my $chr (sort {$org2->{$a}{start}<=>$org2->{$b}{start} } keys %$org2)
      {
	my $y = $graphics_context->height-sprintf("%.0f",$org2->{$chr}{start}*$y_pix_per_bp);
	next if $pv && $y == $pv-1;
	if ($grid)
	  {
	    my $tmp = $y;
	    my $span = sprintf("%.0f",($org2->{$chr}{length}/10*$y_pix_per_bp));
	    for (1..9)
	      {
		$tmp-=$span;
		$graphics_context->line(0, $tmp, $graphics_context->width, $tmp, $span_color);
	      }
	  }

	$graphics_context->line(0, $y, $graphics_context->width, $y, $black);
	my $str_color = $org2->{$chr}{rev} ? $red : $black;
	$graphics_context->string(gdSmallFont, 2, $y-15, $chr, $str_color);
	$data{y}{$pchr}=[$y, $pv] if $pchr;
	$pchr = $chr;
	$pv =$y+1;
      }
    $data{y}{$pchr}=[0,$pv-1];
    
    if ($link_type == 2)
	{
		open (OUT, ">".$basename.".html") || die "$!";

		print OUT qq{
<html><head>
<script src="/CoGe/js/jquery-1.3.2.js"></script>
<script src="/CoGe/js/xhairs.js"></script>
<script type="text/javascript">
/*This array follows the following format:
type (circle, rect), array of coords, href link, mouseover details
Coords will be in [x,y,radius] for cricles and [x1,y1,x2,y2] (diagonal points) for rectangles*/
};

		#Loop through our list of chromosomes and print out the corrosponding click map tag
		my $temp_arrayindex_counter = 0;
		foreach my $xchr (sort keys %{$data{x}})
		{
			my ($x1, $x2) = @{$data{x}{$xchr}};
			next if abs($x1-$x2)<5;
			foreach my $ychr (sort keys %{$data{y}})
			{
				my $tmp = $link;
				$tmp =~ s/XCHR/$xchr/;
				$tmp =~ s/YCHR/$ychr/;
				my ($y1, $y2) = @{$data{y}{$ychr}};
				next if abs($y1-$y2)<3;
				print OUT qq{
location_list[$temp_arrayindex_counter] = ['rect', [$x1, $y1, $x2, $y2], '$tmp', ''];};
				$temp_arrayindex_counter++;
			}
		}

		my ($img) = $basename =~ /([^\/]*$)/;
		print OUT qq{
	//Note: Can not use jQuery docready command here - not entirly sure why, but has something to do with being embbeded as an iframe.
	loadstuff('$img.png');	//This is here instead of in the body tag because when the page is called by SynMap the body tag onLoad is not run for some reason.
</script></head><body>
<canvas id="myCanvas" border="0" width="12" height="12" onmousemove="trackpointer(event);" onmousedown="trackclick(event, 'no');">
	Your browser does not have support for Canvas.
</canvas>	
};


		#Print out the HTML footer.
	print OUT qq{
</map>
</body></html>
};
	close OUT;
      }
  }

sub get_dsg_info
  {
    my %opts = @_;
    my $dsgid = $opts{dsgid};
    my $chr = $opts{chr};
    my $minsize = $opts{minsize};
    my $order = $opts{order};
    my $skip_non_ordered = $opts{skip_non_ordered};
    my $metric=$opts{metric};

    my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);

    unless ($dsg)
      {
	warn "No dataset group found with dbid $dsgid\n";
	return;
      }
    my %data;
    if ($metric && $metric =~ /gene/i)
      {
	my $dbh = DBI->connect($coge->db_connection_string,$coge->db_name,$coge->db_passwd);
	foreach my $gs ($dsg->genomic_sequences)
	  {
	    next if $gs->chromosome =~ /random/i;
	    next if $chr && $chr ne $gs->chromosome;
	    my $tmp_chr = $gs->chromosome;
	    my $query = qq{
SELECT count(distinct(feature_id))
  FROM feature
  JOIN dataset_connector dc using (dataset_id)
 WHERE dataset_group_id = $dsgid
   AND feature_type_id = 3
   AND feature.chromosome = '$tmp_chr'

};
	    my ($res) = $dbh->selectrow_array($query);
	    next unless $res || $res==1;
	    if ($data{$gs->chromosome})
	      {
		warn "Duplicate chromosome:".$gs->chromosome."\n";
	      }
	    $data{$gs->chromosome}{length}=$res;
	    #get gene order
	    $query = qq{
SELECT feature_id
  FROM feature
  JOIN dataset_connector dc using (dataset_id)
 WHERE dataset_group_id = $dsgid
   AND feature_type_id = 3
   AND feature.chromosome = '$tmp_chr'
 ORDER BY feature.start

};
	    my $sth = $dbh->prepare($query);
	    $sth->execute();
	    my $i = 1;
	    while( my $row =$sth->fetchrow_arrayref)
	      {
		$data{$gs->chromosome}{gene}{$row->[0]}=$i;
		$i++;
	      }
	  }
	
      }
    foreach my $gs ($dsg->genomic_sequences)
      {
	next if $gs->chromosome =~ /random/i;
	next if $chr && $chr ne $gs->chromosome;
	my $last = $gs->sequence_length;
	next if $minsize && $minsize > $last;
	if ($data{$gs->chromosome})
	  {
	    warn "Duplicate chromosome:".$gs->chromosome."\n" unless $metric && $metric =~ /gene/i;
	  }
	$data{$gs->chromosome}{length}=$last unless $metric && $metric =~ /gene/i;
	$data{$gs->chromosome}{chr_length} = $last if $data{$gs->chromosome};;
      }
    my %rev;#store chromosomes to be reversed in display
    my @ordered;
    if ($order) #chromsomes have a prespecified order
      {
	my %seen;
	foreach my $item (@$order)
	  {
	    my $chr = $item->{chr};
	    next unless $data{$chr};
	    push @ordered, $chr;
	    $rev{$chr}=1 if $item->{rev};
	    $seen{$chr}=1;
	  }
	foreach my $chr (sort {$data{$b}{chr_length} <=> $data{$a}{chr_length} } keys %data) #get any that were not in @$order
	  {
	    next if $seen{$chr};
	    if ($skip_non_ordered)
	      {
		delete $data{$chr};
	      }
	    else
	      {
		push @ordered, $chr;
	      }
	  }
      }
    else     #chromosomes sorted by length.  Longest first.
      {
	@ordered = sort {$data{$b}{chr_length} <=> $data{$a}{chr_length} } keys %data;
      }
    my $pos = 1;
    foreach my $item (@ordered)
      {
	$data{$item}{start} = $pos;
	$pos += $data{$item}{length};
	$data{$item}{rev}=1 if $rev{$item};
      }
    return \%data;
  }


#this sub is used to try to get a build order of contigs for WGS data against a reference genome
#given a set of scaffolds, it will order them such that make an ordered syntenic path along the reference genome
sub parse_syn_blocks
  {
    my %opts = @_;
    my $file = $opts{file};
    my $blocks1=[];
    my $blocks2=[];
    open (IN, $file) || die "Can't open $file for reading: $!";
    $/ = "\n##";
    while (<IN>) #get blocks
      {
	next unless $_;
	s/##//g;
	my ($block1, $block2) = process_syn_block($_);
	push @$blocks1 , $block1;
	push @$blocks2 , $block2;
      }
    close IN;
    $/="\n";
    #which organism has more pieces (aka chromosomes, contigs, etc.)  we are going to assemble the one with more pieces
#    print Dumper \@blocks;
    my $chrs1={};
    my $chrs2={};
    foreach my $item (@$blocks1)
      {
	$chrs1->{$item->{name}}{$item->{match}}{count}++;
	$chrs1->{$item->{name}}{$item->{match}}{score}+=$item->{score};
      }
    foreach my $item (@$blocks2)
      {
	$chrs2->{$item->{name}}{$item->{match}}{count}++;
	$chrs2->{$item->{name}}{$item->{match}}{score}+=$item->{score};
      }

    #blocks1 will contain fewer chromosomes; blocks2 will be ordered by it.
    my $switched =keys %$chrs1 > keys %$chrs2 ? 1 : 0;
    if ($switched)
      {
	($blocks1, $blocks2) = ($blocks2, $blocks1);
	($chrs1, $chrs2) = ($chrs2, $chrs1);
      }

    my $ordered1 =[]; #storage for ordered chromosomes
    my $ordered2 =[]; #storage for ordered chromosomes
    my %seen;
    #sort blocks for chr1 so that the highest scoring ones are first
    foreach my $chr1 (sort{$chrs1->{$b}{score} <=> $chrs1->{$a}{score}} keys %$chrs1)
      {
	push @$ordered1, {chr=>$chr1};
	my @blocks;
	for (my $i=0; $i< @$blocks1; $i++)
	  {
	    my ($block1, $block2) = ($blocks1->[$i],$blocks2->[$i]);
	    next if $block1->{name} ne $chr1;
	    push @blocks, $block2;
	  }
	#Need to check if a given chromosome in @blocks occurs more than once.  This will happen if there is a segmental duplication, or something along those lines.
	my %block_check;
	foreach my $block (@blocks)
	  {
	    if ($block_check{$block->{name}})
	      {
		$block_check{$block->{name}}=$block if $block->{score} > $block_check{$block->{name}}{score}; #higher score wins
	      }
	    else
	      {
		$block_check{$block->{name}}=$block;
	      }
	  }
	@blocks = values %block_check;#create the non-redundant set
	#print out the blocks in order.  Note ones that are in reverse orientation
	foreach my $block (sort {$a->{match_start} <=> $b->{match_start} }@blocks)
	  {
#	    print $block->{name},"\t", $block->{match_start},"\n";
	    next if $seen{$block->{name}};
	    push @$ordered2, {chr=>$block->{name}, rev=>$block->{rev}};
	    $seen{$block->{name}}=1;
	  }
      }
    ($ordered1, $ordered2) = ($ordered2, $ordered1) if $switched;
    return $ordered1, $ordered2;
  }

sub process_syn_block
  {
    my $block = shift;
    my ($head, @block) = split/\n/, $block;
    my ($seq1, $seq2, $block_num, $score, $num_pairs) = 
      $head =~ /alignment\s+(\S+) vs\. (\S+) .*? #(\d+)\s+score = (\S+).*?pairs: (\d+)/;
    my $rev = $head =~/reverse/ ? 1 : 0;
    my ($seq1_start, $seq1_stop, $seq2_start, $seq2_stop);
    #absolute start and stop can give rise to problems if the ends actually hit something far away from the rest of the sytnenic pairs.  Calculating the "mean" position will circumvent this problem
    my @start1;
    my @stop1;
    my @start2;
    my @stop2;
    foreach my $item (@block)
      {
	chomp $item;
	next unless $item;
	my @item = split /\t/, $item;
	push @start1, $item[2];
	push @stop1, $item[3];
	push @start2, $item[6];
	push @stop2, $item[7];
      }
    #remove the ends;
    @start1 = sort {$a<=>$b} @start1;
    @stop1 = sort {$a<=>$b} @stop1;
    @start2 = sort {$a<=>$b} @start2;
    @stop2 = sort {$a<=>$b} @stop2;
    shift @start1 if scalar(@start1) >3;
    pop @start1 if scalar(@start1) >2;
    shift @stop1 if scalar(@stop1) >3;
    pop @stop1  if scalar(@stop1) >2;
    shift @start2 if scalar(@start2) >3;
    pop @start2 if scalar(@start2) >2;
    shift @stop2 if scalar(@stop2) >3;
    pop @stop2  if scalar(@stop2) >2;

    map {$seq1_start+=$_} @start1;
    map {$seq1_stop+=$_} @stop1;
    map {$seq2_start+=$_} @start2;
    map {$seq2_stop+=$_} @stop2;
    $seq1_start = $seq1_start/scalar(@start1);
    $seq1_stop = $seq1_stop/scalar(@stop1);
    $seq2_start = $seq2_start/scalar(@start2);
    $seq2_stop = $seq2_stop/scalar(@stop2);
    $seq1 =~ s/.*?_//;
    $seq2 =~ s/.*?_//;
    my %seq1 = (
		name=>$seq1,
		start=>$seq1_start,
		stop=>$seq1_stop,
		match_start=>$seq2_start,
		match_stop=>$seq2_stop,
		score=>$score,
		rev=>$rev,
		match=>$seq2,
		);
    my %seq2 = (
		name=>$seq2,
		start=>$seq2_start,
		stop=>$seq2_stop,
		match_start=>$seq1_start,
		match_stop=>$seq1_stop,
		score=>$score,
		rev=>$rev,
		match=>$seq1,
		);
    return \%seq1, \%seq2;
  }



sub get_ksdata
  {
    my %opts = @_;
    my $ks_db = $opts{ks_db};
    my $pairs = $opts{pairs};
    my $type = $opts{ks_type};
    my $min = $opts{min};
    my $max = $opts{max};
    my $log = $opts{log};
    my %data;
    return \%data unless -r $ks_db;
    my $select = "select * from ks_data";
    my $dbh = DBI->connect("dbi:SQLite:dbname=$ks_db","","");
    my $sth = $dbh->prepare($select);
    $sth->execute();
    while (my $data = $sth->fetchrow_arrayref)
      {
	if ($pairs)
	  {
	    next unless $pairs->{$data->[1]}{$data->[2]};
	  }
	my %item = (
		    KS=>$data->[3],
		    KN=>$data->[4],
		    'KN_KS'=>$data->[5],
		   );
	my $val = $item{$type};
	next unless defined $val && $val =~ /\d/;
#	if (defined $min || defined $max)
#	  {
#	  }
	$data{$data->[1]}{$data->[2]}=$val;
	
      }
    $sth->finish();
    undef $sth;
    $dbh->disconnect();
    return \%data;
  }

sub get_pairs
  {
    my %opts = @_;
    my $file = $opts{file};
    my $chr1 = $opts{chr1};
    my $chr2 = $opts{chr2};
    my %data;
    open (IN, $file) || die $!;
    while (<IN>)
      {
	chomp;
	next if /^#/;
	next unless $_;
	my @line = split/\t/;
	my @item1 = split/\|\|/, $line[1];
	my @item2 = split/\|\|/, $line[5];
	next unless $item1[6] && $item2[6];
	if ($chr1)
	  {
	    next unless $item1[0] eq $chr1 || $item2[0] eq $chr1;
	  }
	if ($chr2)
	  {
	    next unless $item1[0] eq $chr2 || $item2[0] eq $chr2;
	  }
	$data{$item1[6]}{$item2[6]}=1;
	$data{$item2[6]}{$item1[6]}=1;
      }
    close IN;
    return \%data;
  }

sub get_range
  {
    my %opts = @_;
    my $data = $opts{data};
    my $min = $opts{min};
    my $max = $opts{max};
    my $log = $opts{log};

    my ($set_max, $set_min, $non_zero_min);
    foreach my $k1 (keys %$data)
      {
	foreach my $k2 (keys %{$data->{$k1}})
	  {
	    next unless defined $data->{$k1}{$k2};
	    $non_zero_min = $data->{$k1}{$k2} unless defined $non_zero_min || $data->{$k1}{$k2} == 0;
	    $set_max = $data->{$k1}{$k2} unless defined $set_max;
	    $set_min = $data->{$k1}{$k2} unless defined $set_min;
	    $set_max = $data->{$k1}{$k2} if $data->{$k1}{$k2} > $set_max;
	    $set_min = $data->{$k1}{$k2} if $data->{$k1}{$k2} < $set_min;
	    $non_zero_min = $data->{$k1}{$k2} if $non_zero_min && $data->{$k1}{$k2} < $non_zero_min && $data->{$k1}{$k2} > 0;
	  }
      }
    #let's minimize the data, if needed
    if (defined $min || defined $max)
      {
	foreach my $k1 (keys %$data)
	  {
	    foreach my $k2 (keys %{$data->{$k1}})
	      {
		my $val = $data->{$k1}{$k2};
		next unless defined $val;
		if ($log)
		  {
		    $val = $non_zero_min if $val == 0;
		    $val = log10($val);
		  }
		if ((defined $min && $val < $min) || (defined $max && $val > $max))
		  {
		    $data->{$k1}{$k2} = undef;
		    delete $data->{$k1}{$k2};
		    delete $data->{$k1} unless keys %{$data->{$k1}};
		  }
	      }
	  }
      }
    if ($log)
      {
	$set_max = log10($set_max) if $set_max > 0;
	$set_min = log10($non_zero_min) if $non_zero_min;
      }
    $set_max = $max if defined $max && $max < $set_max;
    $set_min = $min if defined $min && $min > $set_min;
    return ($set_max, $set_min);
  }

sub get_color
  {
    my %opts = @_;
    my $val = $opts{val};
    unless (defined $val && $val >= 0 && $val <=1)
      {
	print STDERR "in sub get_color val is not [0,1]: $val\n";
	return [0,0,0];
      }
    my @colors = (
		  [255,255,0], #yellow
		  [200,200,0], # green
		  [0,200,0], # green
		  [0,100,100], # green
		  [0,200,200], # cyan
		  [0,0,200], # blue
		  [100,0,100], #purple
		  [200,0,200], #magenta
		  [200,0,0], #red
		  [100,0,0], #red
		  [200,100,0], #orange
		  [255,126,0], #orange

#		  [200,200,0], #yellow
#		  [0,255,0], # green
#		  [0,200,200], # cyan
#		  [0,0,255], # blue
#		  [255,0,255], #magenta
#		  [126,0,126], #purple
#		  [255,0,0], #red
#		  [255,126,0], #orange
		 );
    @colors = reverse @colors;
    my ($index1, $index2) = ((floor((scalar(@colors)-1)*$val)), ceil((scalar(@colors)-1)*$val));

    my $color=[];
    my $step = 1/(scalar (@colors)-1);
    my $scale = $index1*$step;
    my $dist = ($val-$scale)/($step);
#    print join ("\t", $val, $step, $scale, $dist),"\n\t";
    for (my $i=0; $i<=2; $i++)
      {
	my $diff = ($colors[$index1][$i]-$colors[$index2][$i])*$dist;
	push @$color, sprintf("%.0f", $colors[$index1][$i]-$diff);
      }
    return $color;
  }

sub commify
  {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text;
  }


#Print out info on script usage
sub usage
  {
    print qq{
Welcome to $0

dagfile      | d       path to dag file containing all the hits

alignfile    | a       path to .aligncoords file generated by 
                       dagchainer containing just the diags

width        | w       width of image in pixels (1024)

min_chr_size | mcs     minimim size of chromosome to be drawn (0)

dsgid1       | dsg1    database id of dataset group (genome) on x-axis

dsgid2       | dsg2    database id of dataset group (genome) on y-axis

chr1         | c1      only show data for this chromosome on the x axis

chr2         | c2      only show data for this chromosome on the y axis

basename     | b       base path and name for output files (test)

link         | l       link to be used in image map

link_type    | lt      are image map links for chromosome blocks or points:
                       1  ::   blocks  (Use "XCHR","YCHR" which will get the appropriate chr substituted in) 
                       2  ::   points

flip         | f       flip axes (1 flip, 0 don't flip [DEFAULT])

grid         | g       add a positional grid to dotplot

ks_db        | ksdb    specify a sqlite database with synonymous/nonsynonymous data
                       to color syntenic points

ks_type| kst           specify the synonymous data to use (kS, kN, kn_ks) for coloring syntenic points

log                    log10 transform ks datadata (val = 0 set to minimum non-zero val)

max                    max ks val cutoff

min                    min ks val cutoff

assemble               if set to 1, output will try to be assembled based on syntenic thread
                       General assumption is aligning a WGS genome sequence to a reference genome
                       
                       if set to 2, will not add any pieces that are not syntenic

axis_metric  | am      metic (units) for dotplot axes:  nucleotides or genes.  Default: nucleotide distances

help         | h       print this message

};
    exit;
  }
