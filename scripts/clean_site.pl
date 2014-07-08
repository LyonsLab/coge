#!/usr/bin/perl -w

use strict;

my $DIR = shift || ".";
my $new_dir = shift || "new";

mkdir $new_dir unless -d $new_dir;
#mkdir "$new_dir/$DIR" unless -d "$new_dir/$DIR";
process_dir();

sub process_dir
  {
    my $dir = shift || ".";
    my @files;
    opendir (DIR, "$DIR/$dir") || die;
    while (my $file = readdir(DIR))
      {
	next if $file =~ /^\.\.?$/;
	push @files, "$dir/$file";
      }
    closedir DIR;
    foreach my $file (@files)
      {
	if (-d "$DIR/$file")
	  {
	    mkdir "$new_dir/$file";
	    process_dir($file);
	  }
	else
	  {
	    process_file($file);
	  }
      }
  }

sub process_file
  {
    my $file = shift;
    if ($file =~ /\.html$/)
      {
	open (IN, "$DIR/$file") || die "can't open $DIR/$file for reading";
	open (OUT, ">$new_dir/$file") || die;
	while (<IN>)
	  {
	    my $js = qq{
<SCRIPT language="JavaScript" type="text/javascript" src="/CoGe/js/jquery-1.2.2.js"></SCRIPT>
<SCRIPT language="JavaScript" type="text/javascript" src="/CoGe/js/jquery.jdMenu.js"></SCRIPT>
<SCRIPT language="JavaScript" type="text/javascript" src="/CoGe/js/jquery.dimensions.js"></SCRIPT>
<SCRIPT language="JavaScript" type="text/javascript" src="/CoGe/js/jquery.cookie.js"></SCRIPT>
};
	    my $menu = menu();
	    my $styled = qq{
<link rel="stylesheet" type="text/css" href="/CoGe/css/jdMenu.css" />
<link rel="stylesheet" type="text/css" href="/CoGe/css/jdMenu.slate.css" />
};
	    my $animate = qq{
<SCRIPT language="JavaScript">

jQuery(document).ready(function(){
  var anim_speed = jQuery.cookie('anim_speed');
  if (! anim_speed){jQuery.cookie('anim_speed',600);}
  pageObj = new Object();
  pageObj.speed = jQuery.cookie('anim_speed') || 600;
  define_animation_speed(pageObj.speed)
  //alert(pageObj.speed);
  if (!jQuery.browser.mozilla)browserAlert();
});

jQuery(function(){
jQuery('ul.jd_menu').jdMenu();
});

function browserAlert()
{
 alert('Unfortunately, CoGe has been designed to work with the FireFox browser and may not function appropriately with your browser.  Please use FireFox, which can be downloaded from http://www.mozilla.com/en-US/firefox/');
}

function animate_section(name, style, speedy, show) {
  if (show) {alert (name+"::"+style)}
  if (speedy == undefined) speedy = pageObj.speed;
  //alert('speedy is '+speedy);
  if (!style) style = "toggle";
  if (!speedy) {
    if ((style == "show") || style == "slideDown" || style == "fadeIn") {jQuery(name).show();}
    else if ((style == "hide") || style == "slideUp" || style == "fadeOut") {jQuery(name).hide();}
    else {jQuery(name).toggle();}
   }
  else if (style == "toggle") {jQuery(name).toggle(speedy);}
  else if (style == "hide") {jQuery(name).hide(speedy);}
  else if (style == "show") {jQuery(name).show(speedy);}
  else if (style == "slideToggle") {jQuery(name).slideToggle(speedy);}
  else if (style == "slideUp") {jQuery(name).slideUp(speedy);}
  else if (style == "slideDown") {jQuery(name).slideDown(speedy);}
  else if (style == "fadeIn") {jQuery(name).fadeIn(speedy);}
  else if (style == "fadeOut") {jQuery(name).fadeOut(speedy);}
  return speedy;
 }

 function define_animation_speed(val){
   //alert(val);
   val *= 1;
   pageObj.speed = val;
   jQuery.cookie('anim_speed',val);
   var menu_html;

   if (val == 0){
     menu_html = '<li><a href="#" target="_self" onClick="define_animation_speed(0)">&radic; Off</a></li><li><a href="#" target="_self" onClick="define_animation_speed(600)">Slow</a></li><li><a href="#" target="_self" onClick="define_animation_speed(400)">Medium</a></li><li><a href="#" target="_self" onClick="define_animation_speed(200)">Fast</a></li>';}
   else{
     if (val == 600){
      menu_html = '<li><a href="#" target="_self" onClick="define_animation_speed(0)">Off</a></li><li><a href="#" target="_self" onClick="define_animation_speed(600)">&radic; Slow</a></li><li><a href="#" target="_self" onClick="define_animation_speed(400)">Medium</a></li><li><a href="#" target="_self" onClick="define_animation_speed(200)">Fast</a></li>';}
     else if (val == 400){
      menu_html = '<li><a href="#" target="_self" onClick="define_animation_speed(0)">Off</a></li><li><a href="#" target="_self" onClick="define_animation_speed(600)">Slow</a></li><li><a href="#" target="_self" onClick="define_animation_speed(400)">&radic; Medium</a></li><li><a href="#" target="_self" onClick="define_animation_speed(200)">Fast</a></li>';}
     else{
      menu_html = '<li><a href="#" target="_self" onClick="define_animation_speed(0)">Off</a></li><li><a href="#" target="_self" onClick="define_animation_speed(600)">Slow</a></li><li><a href="#" target="_self" onClick="define_animation_speed(400)">Medium</a></li><li><a href="#" target="_self" onClick="define_animation_speed(200)">&radic; Fast</a></li>';
     }
     }
     jQuery('.animations').html(menu_html);
}

</SCRIPT>
};
#	    s/(<\/head>)/$js\n$styled\n$1/;
#	    s/(<body.*?>)/$1\n$animate\n$menu/;
	    s/(<\/head>)/$styled\n$1/;
	    s/(<body.*?>)/$1\n$menu/;
	    s/Â//g;
	    s/â€™/'/g;
	    s/â€œ/"/g;
	    s/â€/"/g;
	    s/>\xa0</><br></g;
	    s/\xa0/ /g;
#	    s/å//g;
#	    s/overflow: hidden/overflow: visible/g;
	    s/http:\/\/\/?CoGe/\/CoGe/ig;
#	    s/http:\/\/\/?coge/\/CoGe/ig;
	    if (/<\/body>/i)
	      {
		my $tmp = qq{
<script type="text/javascript">
var gaJsHost = (("https:" == document.location.protocol) ? "https://ssl." : "http://www.");
document.write(unescape("%3Cscript src='" + gaJsHost + "google-analytics.com/ga.js' type='text/javascript'%3E%3C/script%3E"));
</script>
<script type="text/javascript">
var pageTracker = _gat._getTracker("UA-3802985-1");
pageTracker._initData();
pageTracker._trackPageview();
</script>
</BODY>
};
		s/<\/body>/$tmp/i;
	      }
	    print OUT $_;
	  }
	close IN;
	close OUT;
      }
    else
      {
	`cp "$DIR/$file" "$new_dir/$file"`;
      }
  }

sub menu
  {
   return qq{
<TABLE width = 100%>
 <TR>
  <TD align=right>
    <ul class="jd_menu jd_menu_slate" style="position: absolute;right: 0px;top: 0px;">
     <li><a href="/CoGe/" target="_self">Home</a></li>
    </ul>
  </TD>
 </TR>
</TABLE>

};
  }
sub menu_complete
  {
   return qq{
<TABLE width = 100%>
 <TR>
  <TD align=right>
    <ul class="jd_menu jd_menu_slate" style="position: absolute;right: 0px;top: 0px;">
     <li><a href="/CoGe/" target="_self">Home</a></li>
     <li>Applications
      <ul>
      <li><a href="/CoGe/GEvo.pl" target="_self">GEvo</a></li>
      <li><a href="/CoGe/FeatView.pl" target="_self">FeatView</a></li>
      <li><a href="/CoGe/MSAView.pl" target="_self">MSAView</a></li>
      <li><a href="/CoGe/TreeView.pl" target="_self">TreeView</a></li>
      <li><a href="/CoGe/GenomeView.pl" target="_self">GenomeView</a></li>
      <li><a href="/CoGe/CoGeBlast.pl" target="_self">CoGeBlast</a></li>
      </ul>
    </li>
    <li>Downloads
     <ul>
      <li><a href="/CoGe/docs/help/Precomputed_menus_and_lists/CoGe_menus.html" target="_self">Syntenic gene sets</a></li>

     </ul>
    </li>
    <li>Preferences
     <ul>
      <li><a href="#" target="_self">Animation &raquo; </a>
        <ul class="animations">
          <li><a href="#" target="_self" onClick="define_animation_speed(0)">Off</a></li>
          <li><a href="#" target="_self" onClick="define_animation_speed(600)">&radic; Slow</a></li>
          <li><a href="#" target="_self" onClick="define_animation_speed(400)">Medium</a></li>
          <li><a href="#" target="_self" onClick="define_animation_speed(200)">Fast</a></li>
         </ul>
       </li>
      </ul>
     </li>
    <li>Help
     <ul>
          <li><a href="./docs/help/<TMPL_VAR NAME=HELP>" target="_blank">Page Docs</a></li>
          <li><a href="/CoGe/docs/help/CoGe/Overview.html" target="_blank">About CoGe</a></li>
          <li><a href="/CoGe/docs/help/Tutorials/Tutorials.html" target="_blank">Tutorials</a></li>
          <li><a href="/CoGe/docs/help/Contact/Contact.html" target="_blank">Contact or Site Us</a></li>
    </ul>
    </li>
    </ul>
  </TD>
 </TR>
</TABLE>

};
  }
