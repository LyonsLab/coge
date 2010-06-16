#!/usr/bin/perl -w

use strict;
use CGI;


my $FORM = new CGI;

my $text = $FORM->param('text');
$text = qq{<div style="
  float: left;
  width: 300px;
  height: 182px;
  background-image: url(/CoGe/picts/animations/DNA_orbit_animated_small-side.gif);
">
<span style="
  font-style: oblique;
  color: #990000;
  font-size: large;
  text-decoration: blink;
">
Loading. . .
</span></div>} unless $text;

print $FORM->header;
print $text;
