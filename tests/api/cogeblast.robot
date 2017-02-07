# -*- coding: robot -*-
*** Settings ***
Library		HttpLibrary.HTTP
Documentation	blah blah blah
...		blah blah blah
Resource	resource.robot

*** Variables ***

*** Test Cases ***
Submit BLAST Job
	[Tags]  auth-required
	Create Session  coge    ${API_URL}  verify=true
	${document1}=	Catenate
	...	{ 
	...		"type": "blast", 
	...		"parameters": {
	...			"genomes": [ 16911 ], 
	...			"query_seq": "ATAAAACCAGTACAATTTTACTAAATTATGATGTTAAAGTAAGAGACTTTGAAGAAGACG",
	...			"program": "blastn" 
	...		}
	...	}
	Should Be Valid JSON	${document1}
	${content1}=	Parse Json	${document1}	
        ${headers1}=    Create Dictionary       Content-Type=application/json
        ${resp1}=       Put Request     coge	${JOBS}/${AUTH_PARAMS}	data=${content1}	headers=${headers1}
        Should Be Equal As Strings	${resp1.status_code}	201
        Dictionary Should Contain Item	${resp1.json()}	success	True
        Dictionary Should Contain Key	${resp1.json()}	id
        Set Suite Variable	${id1}	${resp1.json()["id"]}
	Log	${id1}
