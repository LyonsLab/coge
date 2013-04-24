#!/usr/bin/perl

use warnings;
use strict;

use CGI::Application::Dispatch;

CGI::Application::Dispatch->dispatch(
	prefix  => 'CoGe::Services::JBrowse',
	table => [
		'config/refseq'					=> { app => 'Configuration', rm => 'refseq_config' },
		'config/tracks'					=> { app => 'Configuration', rm => 'track_config' },
		'sequence/:genome/stats/global'	=> { app => 'Sequence', 	 rm => 'stats_global' },
		'sequence/:genome/features/:chr'=> { app => 'Sequence', 	 rm => 'features' },
		'annotation/:ds/stats/global'	=> { app => 'Annotation', 	 rm => 'stats_global' },
		'annotation/:ds/features/:chr'	=> { app => 'Annotation', 	 rm => 'features' },
		'experiment/:exp/stats/global'	=> { app => 'Experiment', 	 rm => 'stats_global' },
		'experiment/:exp/features/:chr'	=> { app => 'Experiment', 	 rm => 'features' },
	],
);
