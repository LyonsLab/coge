define([
		   'dojo/_base/declare',
		   'JBrowse/Store/SeqFeature',
		   'JBrowse/Model/SimpleFeature'
	   ],
	   function( declare, SeqFeatureStore, SimpleFeature ) {

return declare(SeqFeatureStore, {

	constructor: function(args) {
		// perform any steps to initialize your new store.  
	},

	getGlobalStats: function(statsCallback, errorCallback) {
		statsCallback({
			"scoreMin" : -1,
			"scoreMax" : 1
		});
	},

	getFeatures: function(query, featureCallback, finishCallback, errorCallback) {
		var hits = this.config.results.get_hits(query.ref, query.start, query.end);
		if (hits) {
			var coge = this.config.config.coge;
			var from = hits[1];
			var to = hits[2];
			hits = hits[0];
			for (var i=from; i<=to; i++) {
				var hit = hits[i];
				if (!coge.data_type || coge.data_type == 1) {
					var strand = hit[2] == 1 ? 1 : -1;
					featureCallback(new SimpleFeature({ data: {
						id: coge.id,
						start: hit[0],
						end: hit[1] == hit[0] ? hit[1] + 1 : hit[1],
						strand: strand,
						score: hit[3] * strand,
						score2: hit[4]
					}}));
				} else if (coge.data_type == 2)
					featureCallback(new SimpleFeature({ data: {
						id: coge.id,
						start: hit[0],
						end: hit[1],
						type: hit[2].toLowerCase() == 'snp' ? hit[2] + hit[4] + 'to' + hit[5] : hit[2],
						ref: hit[4],
						alt: hit[5],
						score: hit[6],
						info: hit[7],
						name: hit[2] + ' ' + hit[4] + ' > ' + hit[5]
					}}));
				else if (coge.data_type == 3) {
					var strand = hit[2] == '1' ? 1 : -1;
					featureCallback(new SimpleFeature({ data: {
						id: coge.id,
						start: hit[0],
						end: hit[1],
						strand: strand,
						name: hit[3]
					}}));
				} else if (coge.data_type == 4) {
					var strand = hit[2] == '1' ? 1 : -1;
					var name = /ID\=([\w\.\-]+)\;/.exec(hit[5]);
					if (name != null)
						name = name[1];
					else
						name = hit[5];
					featureCallback(new SimpleFeature({ data: {
						id: coge.id,
						start: hit[0],
						end: hit[1],
						strand: strand,
						score: hit[4] * strand,
						info: hit[5],
						name: name
					}}));
				}
			}
		}
		finishCallback();
	}
});
});