#!/usr/bin/perl

use warnings;
use strict;

use CGI::Application::Dispatch;

CGI::Application::Dispatch->dispatch(
	prefix  => 'CoGe::Services::JBrowse',
	table => [
		'config/refseq'						=> { app => 'Configuration', rm => 'refseq_config' },
		'config/tracks'						=> { app => 'Configuration', rm => 'track_config' },
		'sequence/:gid/stats/global'		=> { app => 'Sequence', 	 rm => 'stats_global' },
		'sequence/:gid/features/:chr'		=> { app => 'Sequence', 	 rm => 'features' },
		'annotation/:dsid/stats/global'		=> { app => 'Annotation', 	 rm => 'stats_global' },
		'annotation/:dsid/features/:chr'	=> { app => 'Annotation', 	 rm => 'features' },
		'experiment/:eid/stats/global'		=> { app => 'Experiment', 	 rm => 'stats_global' },
		'experiment/:eid/features/:chr'		=> { app => 'Experiment', 	 rm => 'features' },
		'experiment/notebook/:nid/stats/global'		=> { app => 'Experiment', 	 rm => 'stats_global' },
		'experiment/notebook/:nid/features/:chr'	=> { app => 'Experiment', 	 rm => 'features' },
		'experiment/genome/:gid/stats/global'		=> { app => 'Experiment', 	 rm => 'stats_global' },
		'experiment/genome/:gid/features/:chr'		=> { app => 'Experiment', 	 rm => 'features' },
	],
);
