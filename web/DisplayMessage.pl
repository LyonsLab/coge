#!/usr/bin/perl -w

use strict;
use CGI;
use CoGe::Accessory::Web;
no warnings 'redefine';

my $FORM   = new CGI;
my $P      = CoGe::Accessory::Web::get_defaults();
my $server = $P->{SERVER};

my $text = $FORM->param('text');
$text = qq{<div style="
  float: left;
  width: 300px;
  height: 182px;
  background-image: url($server/picts/animations/DNA_orbit_animated_small-side.gif);
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
