#!/usr/bin/perl

use warnings;
use strict;

use CGI::Application::Dispatch;

CGI::Application::Dispatch->dispatch(
    table => [

        # Data Services
        'sequence/:gid/:chr?' => {
            prefix => 'CoGe::Services::Data',
            app    => 'Sequence',
            rm     => 'get'
        },
        'organism/search' => {
            prefix => 'CoGe::Services::Data',
            app    => 'Organism',
            rm     => 'search'
        },
        'genome/search' => {
            prefix => 'CoGe::Services::Data',
            app    => 'Genome',
            rm     => 'search'
        },
        'genome/load' => {
            prefix => 'CoGe::Services::Data',
            app    => 'Genome',
            rm     => 'load'
        },
        'notebook/create' => {
            prefix => 'CoGe::Services::Data',
            app    => 'Notebook',
            rm     => 'create'
        },
        'notebook/delete/:nid' => {
            prefix => 'CoGe::Services::Data',
            app    => 'Notebook',
            rm     => 'delete'
        },
        'feature/experiments/:fid/' => {
            prefix => 'CoGe::Services::Data',
            app    => 'Feature',
            rm     => 'get'
        },

        # JBrowse Services
        'config/refseq' => {
            prefix => 'CoGe::Services::JBrowse',
            app    => 'Configuration',
            rm     => 'refseq_config'
        },
        'config/tracks' => {
            prefix => 'CoGe::Services::JBrowse',
            app    => 'Configuration',
            rm     => 'track_config'
        },
        'sequence/:gid/stats/global' => {
            prefix => 'CoGe::Services::JBrowse',
            app    => 'Sequence',
            rm     => 'stats_global'
        },
        'sequence/:gid/features/:chr' => {
            prefix => 'CoGe::Services::JBrowse',
            app    => 'Sequence',
            rm     => 'features'
        },
        'annotation/:dsid/stats/global' => {
            prefix => 'CoGe::Services::JBrowse',
            app    => 'Annotation',
            rm     => 'stats_global'
        },
        'annotation/:dsid/features/:chr' => {
            prefix => 'CoGe::Services::JBrowse',
            app    => 'Annotation',
            rm     => 'features'
        },
        'experiment/:eid/stats/global' => {
            prefix => 'CoGe::Services::JBrowse',
            app    => 'Experiment',
            rm     => 'stats_global'
        },
        'experiment/:eid/stats/region/:chr' => {
            prefix => 'CoGe::Services::JBrowse',
            app    => 'Experiment',
            rm     => 'stats_region'
        },
        'experiment/:eid/stats/regionFeatureDensities/:chr' => {
            prefix => 'CoGe::Services::JBrowse',
            app    => 'Experiment',
            rm     => 'stats_regionFeatureDensities'
        },
        'experiment/:eid/features/:chr' => {
            prefix => 'CoGe::Services::JBrowse',
            app    => 'Experiment',
            rm     => 'features'
        },
        'experiment/notebook/:nid/stats/global' => {
            prefix => 'CoGe::Services::JBrowse',
            app    => 'Experiment',
            rm     => 'stats_global'
        },
        'experiment/notebook/:nid/stats/regionFeatureDensities/:chr' => {
            prefix => 'CoGe::Services::JBrowse',
            app    => 'Experiment',
            rm     => 'stats_regionFeatureDensities'
        },
        'experiment/notebook/:nid/features/:chr' => {
            prefix => 'CoGe::Services::JBrowse',
            app    => 'Experiment',
            rm     => 'features'
        },
        'experiment/genome/:gid/stats/global' => {
            prefix => 'CoGe::Services::JBrowse',
            app    => 'Experiment',
            rm     => 'stats_global'
        },
        'experiment/genome/:gid/stats/regionFeatureDensities/:chr' => {
            prefix => 'CoGe::Services::JBrowse',
            app    => 'Experiment',
            rm     => 'stats_regionFeatureDensities'
        },
        'experiment/genome/:gid/features/:chr' => {
            prefix => 'CoGe::Services::JBrowse',
            app    => 'Experiment',
            rm     => 'features'
        },
    ],
);
