# -*- coding: robot -*-
*** Settings ***
Documentation	blah blah blah
...		blah blah blah
Resource	resource.robot

*** Test Cases ***
Organism Search
	Create Session	coge	${API_URL}  verify=true
	${resp}=        Get Request	    coge	 ${ORGANISMS}/search/col-0
	Should Be Equal As Strings	${resp.status_code}	200
	Dictionary Should Contain Key	${resp.json()}		organisms

Organism Fetch
	Create Session  coge    ${API_URL}  verify=true
	${resp}=	Get Request	coge	${ORGANISMS}/1
	Should Be Equal As Strings	${resp.status_code}	200
	Dictionary Should Contain Item	${resp.json()}	id	1
	Dictionary Should Contain Item	${resp.json()}	name	Arabidopsis thaliana Col-0 (thale cress)	
	Dictionary Should Contain Item	${resp.json()}	description	Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliophyta; eudicotyledons; core eudicotyledons; rosids; eurosids II; Brassicales; Brassicaceae; Arabidopsis
	Dictionary Should Contain Key	${resp.json()}	genomes

#Organism Add
#        Create Session  coge    ${API_URL} verify=true
#        ${content}=	Evaluate        json.load(open('organism_add.json', 'r'))       json
#	${headers}=	Create Dictionary	Content-Type=application/json
#	${resp}=	Put Request	coge	${ORGANISMS}	data=${content}	headers=${headers}
#	Should Be Equal As Strings	${resp.status_code}	201
#	Dictionary Should Contain Key	${resp.json()}	id

