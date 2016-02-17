/* 
 * Entrez JS
 * 
 */

function Entrez(options) {
	this.database = options.database;
	this.debug = options.debug;
    this.initialize(options);
}

$.extend(Entrez.prototype, {
    initialize: function(options) {

    },
    
    search: function(term) {
    	var self = this;
    	var url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi";
    	return $.ajax({
    		type: 'GET',
    		url: url,
    		dataType: 'xml',
    		cache: false,
    		data: {
    			db: this.database,
    			term: term
    		}
    	}).then(function(xml) {
    		if (self.debug) console.log(xml);
    		if (xml) {
    			var $xml = $(xml),
			    	$id = $xml.find("Id");
    			
    			if ($id.length == 1) {
    				if (self.debug) console.log($id.text());
    				return $id.text();
    			}
    		}
    	});
    },
    
    fetch: function(id) {
    	var self = this;
    	var url ="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi";
    	return $.ajax({
    		type: 'GET',
    		url: url,
    		dataType: 'xml',
    		cache: false,
    		data: {
    			db: this.database,
    			id: id
    		}
    	}).then(function(xml) {
    		if (self.debug) console.log(xml);
    		if (xml) {
    			$xml = $(xml);
//    			var expXml = $xml.find("Item[Name='ExpXml']").text();
//    			var $expXml = $($.parseXML(expXml));
//    			var title = $expXml.find("Title").text();
    			var runs = $xml.find("Item[Name='Runs']").text();
    			var $runs = $($.parseXML(runs));
    			var $run = $runs.find("Run");
    			var accn = $run.attr("acc");
    			var size = $run.attr("total_bases");
    			return {
    				accn: accn,
    				size: size
    			}
    		}
    	});
    }
    
});