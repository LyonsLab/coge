#!/usr/bin/perl -w

use strict;
use GD;
use Getopt::Long;
use CoGeX;
use Data::Dumper;

use vars qw($dagfile $alignfile $width $link $min_chr_size $orgid1 $orgid2 $help $coge $gd $CHR1 $CHR2 $basename $link_type);


GetOptions(
	   "dagfile|d=s"=>\$dagfile,
	   "alignfile|a=s"=>\$alignfile,
	   "width|w=i"=>\$width,
	   "link|l=s"=>\$link,
	   "link_type|lt=i"=>\$link_type,
	   "min_chr_size|mcs=i"=>\$min_chr_size,
	   "orgid1|o1=i"=>\$orgid1,
	   "orgid2|o2=i"=>\$orgid2,
	   "help|h" =>\$help,
	   "chr1|c1=s"=>\$CHR1,
	   "chr2|c2=s"=>\$CHR2,
	   "basename|b=s"=>\$basename,
	   
	   );
usage() if $help;
usage() unless -r $dagfile;

$basename = "test" unless $basename;
$width = 1024 unless $width;

$coge = CoGeX->dbconnect();


my $org1info = get_org_info(oid=>$orgid1, chr=>$CHR1, minsize=>$min_chr_size);
my $org1length =0;
map {$org1length+=$_->{length}} values %$org1info;
my $org2info = get_org_info(oid=>$orgid2, chr=>$CHR2, minsize=>$min_chr_size);
my $org2length =0;
map {$org2length+=$_->{length}} values %$org2info;

my $height = sprintf("%.0f", $width*$org2length/$org1length);
$height = $width if ($height > 5*$width) || ($height <  $width/5);
my $x_bp_per_pix = sprintf("%.0f", $org1length/$width);
my $x_pix_per_bp = 1/$x_bp_per_pix;
my $y_bp_per_pix = sprintf("%.0f", $org2length/$height);
my $y_pix_per_bp = 1/$y_bp_per_pix;
#print STDERR join ("\n", $width."x". $height, $org1length."x".$org2length, $x_bp_per_pix, $y_bp_per_pix, $x_pix_per_bp, $y_pix_per_bp);

my $gd = new GD::Image($width, $height);
my $white = $gd->colorResolve(255,255,255);
$gd->fill(1,1,$white);

draw_chromosome_grid(gd=>$gd, org1=>$org1info, org2=>$org2info, x_pix_per_bp=>$x_pix_per_bp, y_pix_per_bp=>$y_pix_per_bp, link=>$link, link_type=>$link_type);
draw_dots(gd=>$gd, file=>$dagfile, org1=>$org1info, org2=>$org2info, x_pix_per_bp=>$x_pix_per_bp, y_pix_per_bp=>$y_pix_per_bp, link_type => $link_type, orgid1=>$orgid1, orgid2=>$orgid2);
my $add = 1 if $orgid1 eq $orgid2;
draw_dots(gd=>$gd, file=>$alignfile, org1=>$org1info, org2=>$org2info, x_pix_per_bp=>$x_pix_per_bp, y_pix_per_bp=>$y_pix_per_bp, color=>$gd->colorResolve(0,150,0), size=>2, add_inverse=>$add);

open (OUT, ">".$basename.".png") || die "$!";
binmode OUT;
print OUT $gd->png;
close OUT;
sub draw_dots
  {
    my %opts = @_;
    my $file = $opts{file};
    my $gd = $opts{gd};
    my $org1 = $opts{org1};
    my $org2 = $opts{org2};
    my $x_pix_per_bp=$opts{x_pix_per_bp};
    my $y_pix_per_bp=$opts{y_pix_per_bp};
    my $size = $opts{size} || 1;
    my $color = $opts{color};
    my $add_inverse = $opts{add_inverse};
    my $link_type = $opts{link_type} || 0;
    my $oid1 = $opts{orgid1};
    my $oid2 = $opts{orgid2};
     
    $color = $gd->colorResolve(150,150,150) unless $color;

    open (IN, $file) || die "$!";
    my @feats;
    my %points;
    while (<IN>)
      {
	chomp;
	next if /^#/;
	next unless $_;
	my @line = split /\t/;
	my ($a, $chr1) = split /_/,$line[0],2;
	my ($b, $chr2) = split /_/,$line[4],2;
	my $special = 0; #stupid variable name so that if we a viewing a single chr to single chr comparison within the same organism, this will make collinear matches appear on the inverse section 
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

	my $org1 = $opts{org1};
	my $org2 = $opts{org2};
	my ($xmin) = sort ($line[2], $line[3]);
	my $x = sprintf("%.0f",($org1->{$chr1}{start}+$xmin+abs($line[3]-$line[2])/2)*$x_pix_per_bp);
	my ($ymin) = sort ($line[6], $line[7]);
	my $y = sprintf("%.0f",($org2->{$chr2}{start}+$ymin+abs($line[7]-$line[6])/2)*$y_pix_per_bp);
	($x,$y) = ($y, $x) if $special;
#	print STDERR $x,"x", $y,"\n";
	$gd->arc($x, $gd->height-$y, $size, $size, 0, 360, $color);
	$gd->arc($y, $gd->height-$x, $size, $size, 0, 360, $color) if ($add_inverse && !$CHR1 && $x ne $y);
	$gd->arc($y, $gd->height-$x, $size, $size, 0, 360, $color) if ($add_inverse && $chr1 eq $chr2 && $x ne $y);
	
	if ($link_type == 1)
	  {
	    my @item1 = split /\|\|/, $line[1];
	    my @item2 = split /\|\|/, $line[5];
	    unless ($points{$x}{$y}) #cuts down on the size of the image map
	      {
		push @feats, [$x, $gd->height-$y, $item1[6],$item2[6]];
		$points{$x}{$y}=1;
	      }
	  }
      }
    close IN;

    if ($link_type == 1)
      {
	open (OUT, ">".$basename.".html") || die "$!";
	my $org1name = $coge->resultset('Organism')->find($oid1)->name if $oid1;
	my $org2name = $coge->resultset('Organism')->find($oid2)->name if $oid2;
	print OUT qq{
<html><head>
<link href="/CoGe/css/styled.css" type="text/css" rel="stylesheet"/>
<SCRIPT type="text/javascript">
var ajax = [];function pjx(args,fname,method) { this.target=args[1]; this.args=args[0]; method=(method)?method:'GET'; if(method=='post'){method='POST';} this.method = method; this.r=ghr(); this.url = this.getURL(fname);}function formDump(){ var all = []; var fL = document.forms.length; for(var f = 0;f<fL;f++){ var els = document.forms[f].elements; for(var e in els){ var tmp = (els[e].id != undefined)? els[e].id : els[e].name; if(typeof tmp != 'string'){continue;} if(tmp){ all[all.length]=tmp} } } return all;}function getVal(id) { if (id.constructor == Function ) { return id(); } if (typeof(id)!= 'string') { return id; } var element = document.getElementById(id); if( !element ) { for( var i=0; i<document.forms.length; i++ ){ element = document.forms[i].elements[id]; if( element ) break; } if( element && !element.type ) element = element[0]; } if(!element){ alert('ERROR: Cant find HTML element with id or name: ' + id+'. Check that an element with name or id='+id+' exists'); return 0; } if(element.type == 'select-one') { if(element.selectedIndex == -1) return; var item = element[element.selectedIndex]; return item.value || item.text; } if(element.type == 'select-multiple') { var ans = []; var k =0; for (var i=0;i<element.length;i++) { if (element[i].selected || element[i].checked ) { ans[k++]= element[i].value || element[i].text; } } return ans; } if(element.type == 'radio' || element.type == 'checkbox'){ var ans =[]; var elms = document.getElementsByTagName('input'); var endk = elms.length ; var i =0; for(var k=0;k<endk;k++){ if(elms[k].type== element.type && elms[k].checked && (elms[k].id==id||elms[k].name==id)){ ans[i++]=elms[k].value; } } return ans; } if( element.value == undefined ){ return element.innerHTML; }else{ return element.value; }}function fnsplit(arg) { var url=""; if(arg=='NO_CACHE'){return '&pjxrand='+Math.random()} if((typeof(arg)).toLowerCase() == 'object'){ for(var k in arg){ url += '&' + k + '=' + arg[k]; } }else if (arg.indexOf('__') != -1) { arga = arg.split(/__/); url += '&' + arga[0] +'='+ escape(arga[1]); } else { var res = getVal(arg) || ''; if(res.constructor != Array){ res = [res] } for(var i=0;i<res.length;i++) { url += '&args=' + escape(res[i]) + '&' + arg + '=' + escape(res[i]); } } return url;}pjx.prototype = { send2perl : function(){ var r = this.r; var dt = this.target; this.pjxInitialized(dt); var url=this.url; var postdata; if(this.method=="POST"){ var idx=url.indexOf('?'); postdata = url.substr(idx+1); url = url.substr(0,idx); } r.open(this.method,url,true); ; if(this.method=="POST"){ r.setRequestHeader("Content-Type", "application/x-www-form-urlencoded"); r.send(postdata); } if(this.method=="GET"){ r.send(null); } r.onreadystatechange = handleReturn; }, pjxInitialized : function(){}, pjxCompleted : function(){}, readyState4 : function(){ var rsp = unescape(this.r.responseText); /* the response from perl */ var splitval = '__pjx__'; /* to split text */ /* fix IE problems with undef values in an Array getting squashed*/ rsp = rsp.replace(splitval+splitval+'g',splitval+" "+splitval); var data = rsp.split(splitval); dt = this.target; if (dt.constructor != Array) { dt=[dt]; } if (data.constructor != Array) { data=[data]; } if (typeof(dt[0])!='function') { for ( var i=0; i<dt.length; i++ ) { var div = document.getElementById(dt[i]); if (div.type =='text' || div.type=='textarea' || div.type=='hidden' ) { div.value=data[i]; } else{ div.innerHTML = data[i]; } } } else if (typeof(dt[0])=='function') { dt[0].apply(this,data); } this.pjxCompleted(dt); }, getURL : function(fname) { var args = this.args; var url= 'fname=' + fname; for (var i=0;i<args.length;i++) { url=url + args[i]; } return url; }};handleReturn = function() { for( var k=0; k<ajax.length; k++ ) { if (ajax[k].r==null) { ajax.splice(k--,1); continue; } if ( ajax[k].r.readyState== 4) { ajax[k].readyState4(); ajax.splice(k--,1); continue; } }};var ghr=getghr();function getghr(){ if(typeof XMLHttpRequest != "undefined") { return function(){return new XMLHttpRequest();} } var msv= ["Msxml2.XMLHTTP.7.0", "Msxml2.XMLHTTP.6.0", "Msxml2.XMLHTTP.5.0", "Msxml2.XMLHTTP.4.0", "MSXML2.XMLHTTP.3.0", "MSXML2.XMLHTTP", "Microsoft.XMLHTTP"]; for(var j=0;j<=msv.length;j++){ try { A = new ActiveXObject(msv[j]); if(A){ return function(){return new ActiveXObject(msv[j]);} } } catch(e) { } } return false;}function jsdebug(){ var tmp = document.getElementById('pjxdebugrequest').innerHTML = "<br><pre>"; for( var i=0; i < ajax.length; i++ ) { tmp += '<a href= '+ ajax[i].url +' target=_blank>' + decodeURI(ajax[i].url) + ' </a><br>'; } document.getElementById('pjxdebugrequest').innerHTML = tmp + "</pre>";}function get_pair_info() { var args = get_pair_info.arguments; for( var i=0; i<args[0].length;i++ ) { args[0][i] = fnsplit(args[0][i]); } var l = ajax.length; ajax[l]= new pjx(args,"get_pair_info",args[2]); ajax[l].url = '/CoGe/SynMap.pl?' + ajax[l].url; ajax[l].send2perl(); ;}
4//]]>
</SCRIPT>
</head><body>
<DIV id=pair_info style='position: absolute;left: 0px;top: 0px;'></DIV>
};
	my ($img) = $basename =~ /([^\/]*$)/;
	print OUT qq{

<IMG SRC="$img.png" usemap="#points" border="0" style='position: absolute;left: 0px;top:45px;'>
<span class=xsmall style='position: absolute;left: 0px;top: 45px;'>$org1name: $CHR1</span>
<map name="points">
};
	foreach my $item (@feats)
	  {
	    my ($x, $y, $f1, $f2) = @$item;
	    print OUT qq{
<area shape='circle' coords='$x, $y, 2' href='/CoGe/GEvo.pl?autogo=1&fid1=$f1;fid2=$f2;&drup1=50000&drdown1=50000&drup2=50000&drdown2=50000' onMouseOver="get_pair_info(['args__$f1','args__$f2'],['pair_info']);" target='_blank' >
};
#
	  }
	my $pos = $gd->height+45;
	$pos .='px';
	print OUT qq{
</map>
<span class=xsmall style='position: absolute;left: 0px;top: $pos;'>$org2name: $CHR2</span>
</body></html>
};
	close OUT;
      }


  }

sub draw_chromosome_grid
  {
    my %opts = @_;
    my $gd = $opts{gd};
    my $org1 = $opts{org1};
    my $org2 = $opts{org2};
    my $x_pix_per_bp=$opts{x_pix_per_bp};
    my $y_pix_per_bp=$opts{y_pix_per_bp};
    my $link = $opts{link};
    my $link_type = $opts{link_type};
    my $black = $gd->colorResolve(0,0,0);
    $gd->line(0,0, $gd->width, 0, $black);
    $gd->line(0, $gd->height-1, $gd->width, $gd->height-1, $black);
    $gd->line(0,0, 0, $gd->height, $black);
    $gd->line($gd->width-1,0, $gd->width-1, $gd->height, $black);
    #make x-axis chromosome deliniators
    my %data;
    my $pchr;
    my $pv;
    foreach my $chr (sort {$org1->{$a}{start}<=>$org1->{$b}{start} } keys %$org1)
      {
	my $x = sprintf("%.0f",$org1->{$chr}{start}*$x_pix_per_bp);
	next if $x == $pv-1;
	$gd->line($x, 0, $x, $gd->height, $black);
	$gd->string(gdSmallFont, $x+2, $gd->height-15, $chr, $black);
	$data{x}{$pchr}=[$pv,$x] if $pchr;
	$pchr=$chr;
	$pv = $x+1;
      }
    $data{x}{$pchr}=[$pv,$gd->width];
    $pchr = undef;
    $pv=undef;
    foreach my $chr (sort {$org2->{$a}{start}<=>$org2->{$b}{start} } keys %$org2)
      {
	my $y = $gd->height-sprintf("%.0f",$org2->{$chr}{start}*$y_pix_per_bp);
	next if $y == $pv-1;
	$gd->line(0, $y, $gd->width, $y, $black);
	$gd->string(gdSmallFont, 2, $y-15, $chr, $black);
	$data{y}{$pchr}=[$y, $pv] if $pchr;
	$pchr = $chr;
	$pv =$y+1;
      }
    $data{y}{$pchr}=[0,$pv-1];
    if ($link_type == 2)
      {
	open (OUT, ">".$basename.".html") || die "$!";
	print OUT "<html><head></head><body>\n";
	my ($img) = $basename =~ /([^\/]*$)/;
	print OUT qq{
<IMG SRC="$img.png" usemap="#points" border="0">
<map name="points">
};
	foreach my $xchr (keys %{$data{x}})
	  {
	    my ($x1, $x2) = @{$data{x}{$xchr}};
	    next if abs($x1-$x2)<5;
	    foreach my $ychr (keys %{$data{y}})
	      {
		my $tmp = $link;
		$tmp =~ s/XCHR/$xchr/;
		$tmp =~ s/YCHR/$ychr/;
		my ($y1, $y2) = @{$data{y}{$ychr}};
		next if abs($y1-$y2)<3;
		print OUT qq{
<area href='$tmp' alt="$xchr vs $ychr" title="$xchr vs $ychr" shap=rect coords="$x1, $y1, $x2, $y2">
};
	      }
	  }
	
	print OUT qq{
</map>
</body></html>
};
	close OUT;
      }
  }

sub get_org_info
  {
    my %opts = @_;
    my $oid = $opts{oid};
    my $chr = $opts{chr};
    my $minsize = $opts{minsize};
    my $org = $coge->resultset('Organism')->find($oid);
    unless ($org)
      {
	warn "No organism found with dbid $oid\n";
	return;
      }
    my %data;
    foreach my $ds ($org->current_datasets)
      {
	foreach my $chrtmp ($ds->chromosomes)
	  {
	    next if $chr && $chr ne $chrtmp;
	    my $last = $ds->last_chromosome_position($chrtmp);
	    next if $minsize && $minsize > $last;
	    if ($data{$chrtmp})
	      {
		warn "Duplicate chromosome: $chrtmp\n";
	      }
	    $data{$chrtmp}{length}=$last;
	  }
      }
    my $pos = 1;
    foreach my $item (sort {$data{$b}{length} <=> $data{$a}{length} } keys %data )
      {
	$data{$item}{start} = $pos;
	$pos += $data{$item}{length};
      }
		      
    return \%data;
  }

sub usage
  {
    print qq{
Welcome to $0

dagfile      | d       path to dag file containing all the hits

alignfile    | a       path to .aligncoords file generated by 
                       dagchainer containing just the diags

width        | w       width of image in pixels (1024)

min_chr_size | mcs     minimim size of chromosome to be drawn (0)

orgid1       | oid1    database id of organism on x-axis

orgid2       | oid2    database id of organism on y-axis

chr1         | c1      only show data for this chromosome on the x axis

chr2         | c2      only show data for this chromosome on the y axis

basename     | b       base path and name for output files (test)

link         | l       link to be used in image map

link_type    | lt      are image map links for chromosome blocks or points:
                       1  ::   blocks  (Use "XCHR","YCHR" which will get the appropriate chr substituted in) 
                       2  ::   points

help         | h       print this message

};
    exit;
  }
