# -*- coding: robot -*-
*** Settings ***
Documentation	blah blah blah
...		blah blah blah
Resource	resource.robot

*** Test Cases ***
Genome Search
	Create Session	coge	${API_URL}	
	${resp}=        Get Request	    coge	 ${GENOMES}/search/16911
	Should Be Equal As Strings	${resp.status_code}	200
	Dictionary Should Contain Key	${resp.json()}		genomes

Genome Fetch
	Create Session  coge    ${API_URL}   
	${resp}=	Get Request	coge	${GENOMES}/16911  
	Should Be Equal As Strings	${resp.status_code}	200
	${expected}=	Evaluate	json.load(open('genome_fetch.json', 'r'))	json
	Dictionaries Should Be Equal	${resp.json()}	${expected}

#Genome Fetch Sequence
#        Create Session  coge    ${API_URL}
#        ${resp}=        Get Request     coge    ${GENOMES}/16911/sequence
#        Should Be Equal As Strings      ${resp.status_code}     200
#       Set Test Message        ${resp.json()}
#        ${expected}=    Evaluate        json.load(open('genome_fetch_sequence.txt', 'r'))
#        Dictionaries Should Be Equal    ${resp}  ${expected}

Genome Add
        Create Session  coge    ${API_URL}
        ${content}=	Evaluate        json.load(open('genome_add.json', 'r'))       json
	${headers}=	Create Dictionary	Content-Type=application/json
	${resp}=	Put Request	coge	${GENOMES}	data=${content}	headers=${headers}
	Should Be Equal As Strings	${resp.status_code}	201
	Dictionary Should Contain Key	${resp.json()}	id

