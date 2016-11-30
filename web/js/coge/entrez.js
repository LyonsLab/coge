/* 
 * Entrez.js
 *
 * JS interface to the Entrez EUtils.
 * Only supports the ESearch, ESummary, and EFetch functions at the moment.
 *
 * Requires jQuery (promises and XML) -- TODO use ES6 promises
 *
 */

class Entrez {
    constructor(options) {
        this.database = options.database;
        this.debug    = options.debug;
        this.baseUrl  = options.baseUrl || 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
        this.esearchUrl  = options.esearchUrl  || this.baseUrl + 'esearch.fcgi';
        this.esummaryUrl = options.esummaryUrl || this.baseUrl + 'esummary.fcgi';
        this.efetchUrl   = options.efetchUrl   || this.baseUrl + 'efetch.fcgi';
        this.retmax = 100000; // virtually unlimited number of search results (Entrez max value)
    }

    // Returns a list of integer IDs for the given search terms (database and accession)
    esearch(term) {
    	var self = this;
    	return $.ajax({
    		type: 'GET',
    		url: this.esearchUrl,
    		dataType: 'xml',
    		cache: false,
    		data: {
    			db:     this.database,
    			term:   term,
    			retmax: this.retmax
    		}
    	}).then(function(xml) {
    		if (self.debug) console.log(xml);
    		if (xml) {
    			var $xml = $(xml),
			    	$ids = $xml.find("Id");
    			if (!$ids) return;

    			var result = $ids.map(function() {
    				if (self.debug) console.log(this);
    				return $(this).text();
    			});
    			return result;
    		}
    	});
    }

    // Returns the XML summary for the given database and ID
    esummary(id) {
    	var self = this;
    	return $.ajax({
    		type: 'GET',
    		url: this.esummaryUrl,
    		dataType: 'xml',
    		cache: false,
    		data: {
    			db: this.database,
    			id: id
    		}
    	}).then(function(xml) {
    		if (self.debug) console.log(xml);
    		if (xml)
    			return xml;
    		return;
    	});
    }

    // Returns the XML document for the given database and ID
    efetch(id) {
    	var self = this;
    	return $.ajax({
    		type: 'GET',
    		url: this.efetchUrl,
    		dataType: 'xml',
    		cache: false,
    		data: {
    			db: this.database,
    			id: id
    		}
    	}).then(function(xml) {
    		if (self.debug) console.log(xml);
    		if (xml)
    			return xml;
    		return;
    	});
    }
}

class SRA extends Entrez {
    constructor() {
        super({database: 'sra'});
        this.debug = 1;
        this.SRA_ACCN_TYPE_EXPERIMENT = 1;
        this.SRA_ACCN_TYPE_PROJECT    = 2;
    }

    type(accn) {
        var self = this;
        if (!accn) {
            console.error("SRA.query: Missing required accn");
            return;
        }

        var accn_type = accn.toUpperCase().substring(0, 3);
        switch (accn_type) {
            case 'PRJ' :
            case 'SRP' :
                return self.SRA_ACCN_TYPE_PROJECT;
            case 'SRX' :
            case 'SRR' :
                return self.SRA_ACCN_TYPE_EXPERIMENT;
            default :
                return;
        }
    }

    extract(xml) {
        var $xml = $(xml);

        // Get title
        var exp_xml = $xml.find("Item[Name='ExpXml']").text();
        var $exp = $($.parseXML('<exp>'+exp_xml+'</exp>')).find("Study");
        var title = $exp.attr('name');

        // Get runs
        var runs_xml = $xml.find("Item[Name='Runs']").text();
        var $runs = $($.parseXML('<runs>'+runs_xml+'</runs>')).find("Run"); // mdb added "runs" parent element, COGE-778
        var runs = $runs.map(function() {
            var accn = $(this).attr("acc");
            var size = $(this).attr("total_bases");
            return {
                accn: accn,
                size: size
            }
        });

        return {
            title: title,
            runs: runs
        }
    }
}